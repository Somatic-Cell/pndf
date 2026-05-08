param(
    [string]$Preset = "msvc-release",
    [string]$DataDir = "data",
    [string]$OutDir = "out_pnm_aniso",
    [int[]]$Downsamples = @(4),
    [string[]]$Maps = @("brush", "scratched", "isotropic"),
    [string[]]$Modes = @("normal"),
    [int[]]$Targets = @(),
    [int]$Iterations = 5,
    [int]$FlipPasses = 1,
    [int]$RelocatePasses = 1,
    [double]$RelocateStepPixels = 1.0,
    [int]$RebuildInterval = 10000,
    [int]$ProgressInterval = 5000,
    [switch]$Build,
    [switch]$ListOnly
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$RepoRoot = Resolve-Path (Join-Path $ScriptDir "..")
$Exe = Join-Path $RepoRoot "build\$Preset\pndf_run_pnm_aniso.exe"
$Prepare = Join-Path $RepoRoot "python\prepare_nxy.py"

if ($Build) {
    cmake --build --preset $Preset
}
if (!(Test-Path $Exe)) {
    throw "Executable not found: $Exe. Build first or pass -Build."
}
if (!(Test-Path $Prepare)) {
    throw "prepare_nxy.py not found: $Prepare"
}

function Get-NFromDownsample([int]$ds) {
    if ($ds -le 0) { throw "Downsample must be positive" }
    return [int](1024 / $ds)
}

function Default-Targets([int]$N) {
    $full = $N * $N
    return @(
        [int][Math]::Round($full * 0.25),
        [int][Math]::Round($full * 0.125),
        [int][Math]::Round($full * 0.0625)
    )
}

New-Item -ItemType Directory -Force -Path $OutDir | Out-Null

foreach ($ds in $Downsamples) {
    $N = Get-NFromDownsample $ds
    $targetList = if (@($Targets).Count -gt 0) { $Targets } else { Default-Targets $N }
    foreach ($map in $Maps) {
        $nxy = Join-Path $OutDir ("{0}_ds{1}_nxy.bin" -f $map, $ds)
        if (!(Test-Path $nxy)) {
            $cmd = @("python", $Prepare, "--map", $map, "--data-dir", $DataDir, "--out", $nxy, "--downsample", "$ds")
            Write-Host ($cmd -join " ")
            if (!$ListOnly) { & python $Prepare --map $map --data-dir $DataDir --out $nxy --downsample $ds }
        }
        foreach ($mode in $Modes) {
            foreach ($target in $targetList) {
                $out = Join-Path $OutDir ("{0}_pnm_aniso_{1}_ds{2}_v{3}.meshbin" -f $map, $mode, $ds, $target)
                $args = @(
                    "--input", $nxy,
                    "--output", $out,
                    "--target", "$target",
                    "--mode", $mode,
                    "--iterations", "$Iterations",
                    "--flip-passes", "$FlipPasses",
                    "--relocate-passes", "$RelocatePasses",
                    "--relocate-step-pixels", "$RelocateStepPixels",
                    "--rebuild-interval", "$RebuildInterval",
                    "--progress-interval", "$ProgressInterval"
                )
                Write-Host ($Exe + " " + ($args -join " "))
                if (!$ListOnly) { & $Exe @args }
            }
        }
    }
}
