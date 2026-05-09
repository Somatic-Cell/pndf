param(
    [string]$Preset = "msvc-release",
    [string]$DataDir = "out",
    [string]$OutDir = "out_pnm_aniso_he",
    [int[]]$Downsamples = @(4),
    [string[]]$Maps = @("scratched", "brush", "isotropic"),
    [string[]]$Modes = @("normal", "4d"),
    [int[]]$Targets = @(16384, 8192, 4096),
    [int]$Iterations = 5,
    [int]$FlipPasses = 1,
    [int]$RelocatePasses = 1,
    [double]$RelocateStepPixels = 1.0,
    [double]$MinAreaRatio = 0.05,
    [int]$RebuildInterval = 10000,
    [int]$ProgressInterval = 5000,
    [switch]$Build
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$root = Resolve-Path (Join-Path $PSScriptRoot "..")
if ($Build) {
    cmake --preset $Preset
    cmake --build --preset $Preset
}

$exe = Join-Path $root "build\$Preset\pndf_run_pnm_aniso_he.exe"
if (-not (Test-Path $exe)) {
    throw "Executable not found: $exe"
}

New-Item -ItemType Directory -Force -Path $OutDir | Out-Null

foreach ($ds in $Downsamples) {
    foreach ($map in $Maps) {
        $input = Join-Path $DataDir ("{0}_ds{1}_nxy.bin" -f $map, $ds)
        if (-not (Test-Path $input)) {
            Write-Warning "Input not found, skipping: $input"
            continue
        }
        foreach ($mode in $Modes) {
            foreach ($target in $Targets) {
                $output = Join-Path $OutDir ("{0}_pnm_aniso_he_{1}_ds{2}_v{3}.meshbin" -f $map, $mode, $ds, $target)
                Write-Host "[pnm_aniso_he] map=$map mode=$mode ds=$ds target=$target"
                & $exe `
                    --input $input `
                    --output $output `
                    --target $target `
                    --mode $mode `
                    --iterations $Iterations `
                    --flip-passes $FlipPasses `
                    --relocate-passes $RelocatePasses `
                    --relocate-step-pixels $RelocateStepPixels `
                    --min-area-ratio $MinAreaRatio `
                    --rebuild-interval $RebuildInterval `
                    --progress-interval $ProgressInterval
            }
        }
    }
}
