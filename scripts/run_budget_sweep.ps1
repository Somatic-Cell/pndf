param(
    [string]$Preset = "msvc-release",
    [string]$DataDir = ".",
    [string]$OutDir = "out",
    [string]$Mode = "normal"
)

cmake --preset $Preset
cmake --build --preset $Preset
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null

$exe = "build/$Preset/pndf_run_qem.exe"
$eval = "build/$Preset/pndf_eval_mesh.exe"
if (-not (Test-Path $exe)) { $exe = "build/$Preset/pndf_run_qem" }
if (-not (Test-Path $eval)) { $eval = "build/$Preset/pndf_eval_mesh" }

$maps = @("isotropic", "brush", "scratched")
$targets = @(524288, 262144, 131072, 65536)

foreach ($m in $maps) {
    $bin = Join-Path $OutDir "$m`_nxy.bin"
    python python/prepare_nxy.py --map $m --data-dir $DataDir --out $bin
    foreach ($t in $targets) {
        $mesh = Join-Path $OutDir "$m`_$Mode`_$t.meshbin"
        $csv = Join-Path $OutDir "$m`_$Mode`_$t.vertex_error.csv"
        & $exe --input $bin --output $mesh --target $t --mode $Mode
        & $eval --input $bin --mesh $mesh --csv $csv
    }
}
