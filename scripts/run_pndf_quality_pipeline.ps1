param(
    [string]$Preset = "msvc-release",
    [string]$TextureDir = "texture",
    [string]$OutDir = "out_quality",
    [int[]]$Downsamples = @(4),
    [string[]]$Maps = @("brush", "isotropic", "scratched"),
    [int[]]$Targets = @(16384, 8192, 4096),

    # Baseline QEM modes.
    [string[]]$QemModes = @("normal", "4d"),

    # PNM-aniso HE runs. Keeping both modes is useful, but it can be slow.
    [string[]]$PnmModes = @("normal", "4d"),
    [int]$Iterations = 5,
    [int]$FlipPasses = 1,
    [double]$FlipImprovementEps = 1e-9,
    [int]$FlipMaxAcceptedPerPass = 0,
    [int]$RelocatePasses = 1,
    [double]$RelocateStepPixels = 1.0,
    [double]$MinAreaRatio = 0.05,
    [int]$RebuildInterval = 10000,
    [int]$ProgressInterval = 5000,

    # Evaluation settings.
    [int]$NormalBins = 64,
    [string]$SigmaPixels = "8,16,32",
    [int]$RandomQueries = 32,
    [int]$FeatureQueries = 16,
    [int]$SeamQueries = 6,

    [string]$Python = "python",
    [switch]$Build,
    [switch]$Force,
    [switch]$SkipPrepare,
    [switch]$SkipQem,
    [switch]$SkipPnm,
    [switch]$SkipAudit,
    [switch]$SkipEmbeddingAudit,
    [switch]$SkipDenseNormalMetrics,
    [switch]$SkipPndfEval,
    [switch]$MakeVisuals
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$root = Resolve-Path (Join-Path $PSScriptRoot "..")
$rootPath = $root.Path

function Resolve-MaybeRelativePath {
    param([string]$Path)
    if ([System.IO.Path]::IsPathRooted($Path)) {
        return [System.IO.Path]::GetFullPath($Path)
    }
    return [System.IO.Path]::GetFullPath((Join-Path $rootPath $Path))
}

$TextureDir = Resolve-MaybeRelativePath $TextureDir
$OutDir = Resolve-MaybeRelativePath $OutDir

$nxyDir = Join-Path $OutDir "nxy"
$meshDir = Join-Path $OutDir "meshes"
$logDir = Join-Path $OutDir "logs"
$auditDir = Join-Path $OutDir "audit"
$metricDir = Join-Path $OutDir "metrics_dense_normal"
$pndfDir = Join-Path $OutDir "metrics_sampled_pndf"
$visDir = Join-Path $OutDir "visuals"

foreach ($d in @($OutDir, $nxyDir, $meshDir, $logDir, $auditDir, $metricDir, $pndfDir, $visDir)) {
    New-Item -ItemType Directory -Force -Path $d | Out-Null
}

if ($Build) {
    cmake --preset $Preset
    if ($LASTEXITCODE -ne 0) {
        throw "cmake configure failed with exit code $LASTEXITCODE"
    }

    cmake --build --preset $Preset
    if ($LASTEXITCODE -ne 0) {
        throw "cmake build failed with exit code $LASTEXITCODE"
    }
}

$exeQem = Join-Path $rootPath "build\$Preset\pndf_run_qem.exe"
$exePnm = Join-Path $rootPath "build\$Preset\pndf_run_pnm_aniso_he.exe"

if (-not $SkipQem -and -not (Test-Path $exeQem)) {
    throw "Executable not found: $exeQem"
}
if (-not $SkipPnm -and -not (Test-Path $exePnm)) {
    throw "Executable not found: $exePnm"
}

function Invoke-LoggedCommand {
    param(
        [string]$Name,
        [string]$LogPath,
        [string]$Exe,
        [string[]]$ArgsList
    )

    Write-Host "[$Name] $Exe $($ArgsList -join ' ')"
    New-Item -ItemType Directory -Force -Path (Split-Path $LogPath -Parent) | Out-Null

    # Some of our native tools intentionally print progress logs to stderr.
    # With `$ErrorActionPreference = "Stop"`, Windows PowerShell turns each
    # redirected stderr line into a NativeCommandError and stops the script even
    # when the executable is still healthy. Temporarily relax the preference
    # while running the native process, then check `$LASTEXITCODE` explicitly.
    $oldErrorActionPreference = $ErrorActionPreference
    $hadNativePreference = Test-Path Variable:\PSNativeCommandUseErrorActionPreference
    if ($hadNativePreference) {
        $oldNativePreference = $PSNativeCommandUseErrorActionPreference
    }

    try {
        $ErrorActionPreference = "Continue"
        if ($hadNativePreference) {
            $PSNativeCommandUseErrorActionPreference = $false
        }

        # Do not use Tee-Object directly here: on Windows PowerShell, stderr
        # objects can still be surfaced as NativeCommandError records. Convert
        # every stream object to text and write it ourselves.
        & $Exe @ArgsList 2>&1 | ForEach-Object {
            $line = $_.ToString()
            Write-Host $line
            Add-Content -Path $LogPath -Value $line -Encoding UTF8
        }
        $exitCode = $LASTEXITCODE
    }
    finally {
        $ErrorActionPreference = $oldErrorActionPreference
        if ($hadNativePreference) {
            $PSNativeCommandUseErrorActionPreference = $oldNativePreference
        }
    }

    if ($exitCode -ne 0) {
        throw "Command failed with exit code $exitCode. See log: $LogPath"
    }
}

function Invoke-PythonLogged {
    param(
        [string]$Name,
        [string]$LogPath,
        [string[]]$ArgsList
    )
    Invoke-LoggedCommand -Name $Name -LogPath $LogPath -Exe $Python -ArgsList $ArgsList
}

function Get-NxyPath {
    param([string]$Map, [int]$Ds)
    return (Join-Path $nxyDir ("{0}_ds{1}_nxy.bin" -f $Map, $Ds))
}

function Get-MeshPath {
    param([string]$Method, [string]$Mode, [string]$Map, [int]$Ds, [int]$Target)
    return (Join-Path $meshDir ("{0}_{1}_{2}_ds{3}_v{4}.meshbin" -f $Map, $Method, $Mode, $Ds, $Target))
}

$manifestRows = New-Object System.Collections.Generic.List[object]

foreach ($ds in $Downsamples) {
    foreach ($map in $Maps) {
        $nxy = Get-NxyPath -Map $map -Ds $ds

        if (-not $SkipPrepare) {
            if ($Force -or -not (Test-Path $nxy)) {
                $log = Join-Path $logDir ("prepare_{0}_ds{1}.log" -f $map, $ds)
                Invoke-PythonLogged -Name "prepare_nxy" -LogPath $log -ArgsList @(
                    (Join-Path $rootPath "python\prepare_nxy.py"),
                    "--map", $map,
                    "--data-dir", $TextureDir,
                    "--out", $nxy,
                    "--downsample", [string]$ds
                )
            } else {
                Write-Host "[prepare_nxy] skip existing $nxy"
            }
        }

        if (-not (Test-Path $nxy)) {
            Write-Warning "NXY input not found, skipping map=$map ds=$ds : $nxy"
            continue
        }

        $runs = New-Object System.Collections.Generic.List[object]

        if (-not $SkipQem) {
            foreach ($mode in $QemModes) {
                foreach ($target in $Targets) {
                    $mesh = Get-MeshPath -Method "qem" -Mode $mode -Map $map -Ds $ds -Target $target
                    $log = Join-Path $logDir ("run_{0}_qem_{1}_ds{2}_v{3}.log" -f $map, $mode, $ds, $target)
                    if ($Force -or -not (Test-Path $mesh)) {
                        Invoke-LoggedCommand -Name "qem" -LogPath $log -Exe $exeQem -ArgsList @(
                            "--input", $nxy,
                            "--output", $mesh,
                            "--target", [string]$target,
                            "--mode", $mode
                        )
                    } else {
                        Write-Host "[qem] skip existing $mesh"
                    }
                    $runs.Add([pscustomobject]@{ method="qem"; mode=$mode; map=$map; ds=$ds; target=$target; mesh=$mesh }) | Out-Null
                }
            }
        }

        if (-not $SkipPnm) {
            foreach ($mode in $PnmModes) {
                foreach ($target in $Targets) {
                    $mesh = Get-MeshPath -Method "pnm_he" -Mode $mode -Map $map -Ds $ds -Target $target
                    $log = Join-Path $logDir ("run_{0}_pnm_he_{1}_ds{2}_v{3}.log" -f $map, $mode, $ds, $target)
                    if ($Force -or -not (Test-Path $mesh)) {
                        $argsList = @(
                            "--input", $nxy,
                            "--output", $mesh,
                            "--target", [string]$target,
                            "--mode", $mode,
                            "--iterations", [string]$Iterations,
                            "--flip-passes", [string]$FlipPasses,
                            "--flip-improvement-eps", [string]$FlipImprovementEps,
                            "--relocate-passes", [string]$RelocatePasses,
                            "--relocate-step-pixels", [string]$RelocateStepPixels,
                            "--min-area-ratio", [string]$MinAreaRatio,
                            "--rebuild-interval", [string]$RebuildInterval,
                            "--progress-interval", [string]$ProgressInterval
                        )
                        if ($FlipMaxAcceptedPerPass -gt 0) {
                            $argsList += @("--flip-max-accepted-per-pass", [string]$FlipMaxAcceptedPerPass)
                        }
                        Invoke-LoggedCommand -Name "pnm_he" -LogPath $log -Exe $exePnm -ArgsList $argsList
                    } else {
                        Write-Host "[pnm_he] skip existing $mesh"
                    }
                    $runs.Add([pscustomobject]@{ method="pnm_he"; mode=$mode; map=$map; ds=$ds; target=$target; mesh=$mesh }) | Out-Null
                }
            }
        }

        foreach ($r in $runs) {
            if (-not (Test-Path $r.mesh)) {
                Write-Warning "Mesh not found, skipping evaluation: $($r.mesh)"
                continue
            }

            $tag = "{0}_{1}_{2}_ds{3}_v{4}" -f $r.map, $r.method, $r.mode, $r.ds, $r.target
            $auditLog = Join-Path $auditDir ("{0}.txt" -f $tag)
            $denseCsv = Join-Path $metricDir ("{0}.csv" -f $tag)
            $pndfOut = Join-Path $pndfDir $tag
            $visOut = Join-Path $visDir $tag

            if (-not $SkipAudit) {
                if ($Force -or -not (Test-Path $auditLog)) {
                    $auditArgs = @(
                        (Join-Path $rootPath "scripts\audit_mesh_topology.py"),
                        $r.mesh
                    )
                    if (-not $SkipEmbeddingAudit) {
                        $auditArgs += "--embedding"
                    }
                    Invoke-PythonLogged -Name "audit" -LogPath $auditLog -ArgsList $auditArgs
                } else {
                    Write-Host "[audit] skip existing $auditLog"
                }
            }

            if (-not $SkipDenseNormalMetrics) {
                if ($Force -or -not (Test-Path $denseCsv)) {
                    $log = Join-Path $logDir ("metrics_dense_{0}.log" -f $tag)
                    Invoke-PythonLogged -Name "dense_normal_metrics" -LogPath $log -ArgsList @(
                        (Join-Path $rootPath "python\pndf_mesh_metrics.py"),
                        "--nxy-bin", $nxy,
                        "--meshbin", $r.mesh,
                        "--csv", $denseCsv
                    )
                } else {
                    Write-Host "[dense_normal_metrics] skip existing $denseCsv"
                }
            }

            if (-not $SkipPndfEval) {
                $pndfDone = Join-Path $pndfOut ".done"
                if ($Force -or -not (Test-Path $pndfDone)) {
                    New-Item -ItemType Directory -Force -Path $pndfOut | Out-Null
                    $log = Join-Path $logDir ("metrics_pndf_{0}.log" -f $tag)
                    Invoke-PythonLogged -Name "sampled_pndf" -LogPath $log -ArgsList @(
                        (Join-Path $rootPath "python\pndf_eval_pndf.py"),
                        "--nxy-bin", $nxy,
                        "--meshbin", $r.mesh,
                        "--out-dir", $pndfOut,
                        "--normal-bins", [string]$NormalBins,
                        "--sigma-pixels", $SigmaPixels,
                        "--random-queries", [string]$RandomQueries,
                        "--feature-queries", [string]$FeatureQueries,
                        "--seam-queries", [string]$SeamQueries,
                        "--no-images"
                    )
                    "done" | Out-File -Encoding ascii $pndfDone
                } else {
                    Write-Host "[sampled_pndf] skip existing $pndfOut"
                }
            }

            if ($MakeVisuals) {
                $visDone = Join-Path $visOut ".done"
                if ($Force -or -not (Test-Path $visDone)) {
                    New-Item -ItemType Directory -Force -Path $visOut | Out-Null
                    $log = Join-Path $logDir ("visualize_{0}.log" -f $tag)
                    Invoke-PythonLogged -Name "visualize" -LogPath $log -ArgsList @(
                        (Join-Path $rootPath "python\visualize_pndf_mesh.py"),
                        "--nxy-bin", $nxy,
                        "--meshbin", $r.mesh,
                        "--out-dir", $visOut,
                        "--full-max-edges", "0",
                        "--crop-max-edges", "0",
                        "--render-reconstruction",
                        "--error-heatmap"
                    )
                    "done" | Out-File -Encoding ascii $visDone
                } else {
                    Write-Host "[visualize] skip existing $visOut"
                }
            }

            $manifestRows.Add([pscustomobject]@{
                method = $r.method
                mode = $r.mode
                map = $r.map
                downsample = $r.ds
                target = $r.target
                nxy = $nxy
                mesh = $r.mesh
                audit_log = $auditLog
                dense_normal_csv = $denseCsv
                sampled_pndf_dir = $pndfOut
                visual_dir = $visOut
            }) | Out-Null
        }
    }
}

$manifest = Join-Path $OutDir "manifest.csv"
$manifestRows | Export-Csv -NoTypeInformation -Encoding UTF8 $manifest
Write-Host "Done. Manifest: $manifest"
