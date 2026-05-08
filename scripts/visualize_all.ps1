<#
.SYNOPSIS
  Visualize QEM output meshes for all maps / modes / resolutions.

.DESCRIPTION
  This script repeatedly calls python/visualize_pndf_mesh.py.
  It is designed for the current pndf repository output naming convention:

    out256_normal\scratched_ds4_nxy.bin
    out256_normal\scratched_normal_ds4_v16384.meshbin
    out256_4d\scratched_4d_ds4_v16384.meshbin

  It can either discover all existing meshbin files or use an explicit target list.

.EXAMPLES
  # Visualize all N=256 outputs already present in out256_normal/out256_4d.
  .\scripts\visualize_all.ps1 -Downsamples 4

  # Visualize explicit targets for N=256.
  .\scripts\visualize_all.ps1 `
    -Downsamples 4 `
    -Targets 49152,32768,16384,8192,4096

  # Visualize N=256 and N=512 outputs.
  .\scripts\visualize_all.ps1 -Downsamples 4,2

  # List commands without running them.
  .\scripts\visualize_all.ps1 -Downsamples 4 -ListOnly
#>

param(
    # Repository root. If omitted, inferred as the parent directory of this script.
    [string]$RepoRoot = "",

    # Python executable. Use ".venv\Scripts\python.exe" if desired.
    [string]$PythonExe = "python",

    # Visualizer script relative to RepoRoot unless an absolute path is given.
    [string]$Visualizer = "python\visualize_pndf_mesh.py",

    # Output root for visualization images.
    [string]$VisRoot = "vis",

    # Dataset maps.
    [string[]]$Maps = @("isotropic", "brush", "scratched"),

    # QEM modes.
    [string[]]$Modes = @("normal", "4d"),

    # Downsample factors used when preparing nxy.bin.
    # ds=4 means N=256 for original 1024 texture.
    # ds=2 means N=512, ds=1 means N=1024.
    [int[]]$Downsamples = @(4),

    # Optional explicit target vertex counts. If empty, meshbin files are discovered.
    [int[]]$Targets = @(),

    # Crop behavior passed to visualize_pndf_mesh.py: "auto" or "u0,v0,size".
    [string]$Crop = "auto",

    # Crop size in normalized UV when Crop="auto".
    [double]$CropSize = 0.15,

    # Output crop resolution in pixels.
    [int]$CropPixels = 2048,

    # Render reconstructed normal map from simplified mesh.
    [bool]$RenderReconstruction = $true,

    # Render angular-error heatmap. Requires reconstructed normal map.
    [bool]$ErrorHeatmap = $true,

    # Skip visualization if summary.txt already exists in output directory.
    [switch]$SkipExisting,

    # Print commands but do not execute.
    [switch]$ListOnly,

    # Stop at first failure.
    [switch]$FailFast
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

# Resolve repository root robustly.
# Important on Windows: when PowerShell is launched from Explorer/Admin shell,
# the current directory can be C:\WINDOWS\system32.  Do not use the current
# directory as the default repo root.  If this script is placed in scripts/,
# $PSScriptRoot\.. is the repository root.
if ([string]::IsNullOrWhiteSpace($RepoRoot)) {
    if (-not [string]::IsNullOrWhiteSpace($PSScriptRoot)) {
        $RepoRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot ".."))
    } else {
        $RepoRoot = [System.IO.Path]::GetFullPath((Get-Location).Path)
    }
} else {
    $RepoRoot = [System.IO.Path]::GetFullPath($RepoRoot)
}
Write-Host "RepoRoot: $RepoRoot"

function Resolve-RepoPath {
    param([string]$PathLike)
    if ([System.IO.Path]::IsPathRooted($PathLike)) {
        return [System.IO.Path]::GetFullPath($PathLike)
    }
    return [System.IO.Path]::GetFullPath((Join-Path $RepoRoot $PathLike))
}

function Get-TextureSizeFromDownsample {
    param([int]$Downsample)
    if ($Downsample -le 0) {
        throw "Downsample must be positive, got $Downsample"
    }
    $n = [int](1024 / $Downsample)
    if ($n * $Downsample -ne 1024) {
        Write-Warning "Downsample=$Downsample does not divide 1024 exactly; inferred N=$n"
    }
    return $n
}

function Get-RunDir {
    param([int]$Downsample, [string]$Mode)
    $n = Get-TextureSizeFromDownsample -Downsample $Downsample
    return Join-Path $RepoRoot ("out{0}_{1}" -f $n, $Mode)
}

function Find-NxyBin {
    param([string]$Map, [string]$Mode, [int]$Downsample)

    $modeDir = Get-RunDir -Downsample $Downsample -Mode $Mode
    $n = Get-TextureSizeFromDownsample -Downsample $Downsample

    $candidates = New-Object System.Collections.Generic.List[string]
    $candidates.Add((Join-Path $modeDir ("{0}_ds{1}_nxy.bin" -f $Map, $Downsample)))
    $candidates.Add((Join-Path $modeDir ("{0}_{1}_nxy.bin" -f $Map, $n)))
    $candidates.Add((Join-Path $modeDir ("{0}_nxy.bin" -f $Map)))

    # nxy is identical between normal/4d, so also search sibling mode directories.
    foreach ($m in $Modes) {
        $dir = Get-RunDir -Downsample $Downsample -Mode $m
        $candidates.Add((Join-Path $dir ("{0}_ds{1}_nxy.bin" -f $Map, $Downsample)))
        $candidates.Add((Join-Path $dir ("{0}_{1}_nxy.bin" -f $Map, $n)))
        $candidates.Add((Join-Path $dir ("{0}_nxy.bin" -f $Map)))
    }

    foreach ($p in $candidates) {
        if (Test-Path $p) { return [System.IO.Path]::GetFullPath($p) }
    }

    # Last-resort discovery in outN_* dirs.
    foreach ($m in $Modes) {
        $dir = Get-RunDir -Downsample $Downsample -Mode $m
        if (Test-Path $dir) {
            $found = Get-ChildItem -Path $dir -Filter ("{0}*nxy*.bin" -f $Map) -File -ErrorAction SilentlyContinue | Select-Object -First 1
            if ($found) { return $found.FullName }
        }
    }

    return $null
}

function Find-MeshBins {
    param([string]$Map, [string]$Mode, [int]$Downsample)

    $dir = Get-RunDir -Downsample $Downsample -Mode $Mode
    if (-not (Test-Path $dir)) { return @() }

    if (@($Targets).Count -gt 0) {
        $list = @()
        foreach ($t in $Targets) {
            $patterns = @(
                ("{0}_{1}_ds{2}_v{3}.meshbin" -f $Map, $Mode, $Downsample, $t),
                ("{0}_{1}_v{2}.meshbin" -f $Map, $Mode, $t),
                ("{0}*{1}*v{2}.meshbin" -f $Map, $Mode, $t)
            )
            $hit = $null
            foreach ($pat in $patterns) {
                $items = @(Get-ChildItem -Path $dir -Filter $pat -File -ErrorAction SilentlyContinue)
                if ($items.Count -gt 0) { $hit = $items[0]; break }
            }
            if ($hit) {
                $list += $hit
            } else {
                Write-Warning "Missing mesh for map=$Map mode=$Mode ds=$Downsample target=$t in $dir"
            }
        }
        return $list
    }

    # Discover all meshbin files matching map+mode+downsample if possible.
    $patterns = @(
        ("{0}_{1}_ds{2}_v*.meshbin" -f $Map, $Mode, $Downsample),
        ("{0}_{1}_v*.meshbin" -f $Map, $Mode),
        ("{0}*{1}*v*.meshbin" -f $Map, $Mode)
    )
    $files = @()
    foreach ($pat in $patterns) {
        $files = @(Get-ChildItem -Path $dir -Filter $pat -File -ErrorAction SilentlyContinue)
        if ($files.Count -gt 0) { break }
    }

    return $files | Sort-Object Name
}

function Get-TargetNameFromMesh {
    param([System.IO.FileInfo]$MeshFile)
    $base = [System.IO.Path]::GetFileNameWithoutExtension($MeshFile.Name)
    if ($base -match "_v([0-9]+)$") { return "v$($Matches[1])" }
    return $base
}

function Invoke-Visualizer {
    param(
        [string]$NxyBin,
        [string]$MeshBin,
        [string]$OutDir
    )

    $args = @(
        (Resolve-RepoPath $Visualizer),
        "--nxy-bin", $NxyBin,
        "--meshbin", $MeshBin,
        "--out-dir", $OutDir,
        "--crop", $Crop,
        "--crop-size", ([string]$CropSize),
        "--crop-pixels", ([string]$CropPixels)
    )

    if ($RenderReconstruction) { $args += "--render-reconstruction" }
    if ($ErrorHeatmap) { $args += "--error-heatmap" }

    Write-Host "python $($args -join ' ')"
    if (-not $ListOnly) {
        & $PythonExe @args
        if ($LASTEXITCODE -ne 0) {
            throw "visualize_pndf_mesh.py failed with exit code $LASTEXITCODE"
        }
    }
}

$visualizerPath = Resolve-RepoPath $Visualizer
if (-not (Test-Path $visualizerPath)) {
    throw "Visualizer script not found: $visualizerPath"
}

$summaryRows = New-Object System.Collections.Generic.List[object]

foreach ($ds in $Downsamples) {
    $n = Get-TextureSizeFromDownsample -Downsample $ds
    foreach ($mode in $Modes) {
        foreach ($map in $Maps) {
            $nxy = Find-NxyBin -Map $map -Mode $mode -Downsample $ds
            if (-not $nxy) {
                Write-Warning "Missing nxy.bin for map=$map mode=$mode ds=$ds. Skipping."
                continue
            }

            $meshes = @(Find-MeshBins -Map $map -Mode $mode -Downsample $ds)
            if ($meshes.Count -eq 0) {
                Write-Warning "No meshbin files found for map=$map mode=$mode ds=$ds. Skipping."
                continue
            }

            foreach ($mesh in $meshes) {
                $targetName = Get-TargetNameFromMesh -MeshFile $mesh
                $outDir = Join-Path $RepoRoot (Join-Path $VisRoot ("N{0}\{1}\{2}\{3}" -f $n, $mode, $map, $targetName))
                $summaryFile = Join-Path $outDir "summary.txt"

                $row = [PSCustomObject]@{
                    N = $n
                    Downsample = $ds
                    Mode = $mode
                    Map = $map
                    Target = $targetName
                    NxyBin = $nxy
                    MeshBin = $mesh.FullName
                    OutDir = [System.IO.Path]::GetFullPath($outDir)
                    Status = "pending"
                }

                if ($SkipExisting -and (Test-Path $summaryFile)) {
                    Write-Host "SKIP existing: $outDir"
                    $row.Status = "skipped_existing"
                    $summaryRows.Add($row)
                    continue
                }

                New-Item -ItemType Directory -Force -Path $outDir | Out-Null
                try {
                    Invoke-Visualizer -NxyBin $nxy -MeshBin $mesh.FullName -OutDir $outDir
                    $row.Status = "ok"
                } catch {
                    $row.Status = "failed: $($_.Exception.Message)"
                    Write-Error $row.Status
                    if ($FailFast) { throw }
                }
                $summaryRows.Add($row)
            }
        }
    }
}

$visRootAbs = [System.IO.Path]::GetFullPath((Join-Path $RepoRoot $VisRoot))
New-Item -ItemType Directory -Force -Path $visRootAbs | Out-Null
$summaryCsv = Join-Path $visRootAbs "visualize_all_summary.csv"
$summaryRows | Export-Csv -Path $summaryCsv -NoTypeInformation -Encoding UTF8
Write-Host "Wrote summary: $summaryCsv"
