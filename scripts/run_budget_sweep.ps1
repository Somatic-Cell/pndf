param(
    [string]$Preset = "msvc-release",
    [string]$DataDir = ".",
    [string]$OutDir = "out",
    [string]$Mode = "normal",
    [int]$Downsample = 1,
    [string]$Targets = "",
    [int]$RebuildInterval = 50000,
    [int]$ProgressInterval = 10000
)

cmake --preset $Preset
cmake --build --preset $Preset
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null

$exe = "build/$Preset/pndf_run_qem.exe"
$eval = "build/$Preset/pndf_eval_mesh.exe"
if (-not (Test-Path $exe)) { $exe = "build/$Preset/pndf_run_qem" }
if (-not (Test-Path $eval)) { $eval = "build/$Preset/pndf_eval_mesh" }

if ($Downsample -lt 1) { throw "Downsample must be >= 1" }

if ([string]::IsNullOrWhiteSpace($Targets)) {
    $baseN = [int](1024 / $Downsample)
    if ($baseN * $Downsample -ne 1024) { throw "Downsample must divide 1024 unless explicit -Targets is provided" }
    $full = $baseN * $baseN
    $targetList = @(
        [int]($full * 3 / 4),
        [int]($full / 2),
        [int]($full / 4),
        [int]($full / 8),
        [int]($full / 16)
    ) | Where-Object { $_ -gt 0 } | Select-Object -Unique
} else {
    $targetList = $Targets.Split(',') | ForEach-Object { [int]($_.Trim()) } | Where-Object { $_ -gt 0 } | Sort-Object -Descending | Select-Object -Unique
}
$targetText = ($targetList -join ",")

$maps = @("isotropic", "brush", "scratched")
foreach ($m in $maps) {
    $bin = Join-Path $OutDir "$m`_ds$Downsample`_nxy.bin"
    python python/prepare_nxy.py --map $m --data-dir $DataDir --out $bin --downsample $Downsample

    $prefix = Join-Path $OutDir "$m`_$Mode`_ds$Downsample"
    $profile = Join-Path $OutDir "$m`_$Mode`_ds$Downsample`_profile.csv"

    & $exe `
        --input $bin `
        --output-prefix $prefix `
        --targets $targetText `
        --mode $Mode `
        --rebuild-interval $RebuildInterval `
        --progress-interval $ProgressInterval `
        --profile-csv $profile

    foreach ($t in $targetList) {
        $mesh = "${prefix}_v${t}.meshbin"
        if (Test-Path $mesh) {
            $csv = "${prefix}_v${t}.vertex_error.csv"
            & $eval --input $bin --mesh $mesh --csv $csv
        } else {
            Write-Warning "checkpoint mesh not found: $mesh"
        }
    }
}
