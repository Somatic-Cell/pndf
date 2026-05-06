build\msvc-release\pndf_run_qem.exe `
    --input out\scratched_nxy.bin `
    --output out\scratched_normal_262k.meshbin `
    --target 262144 `
    --mode normal

build\msvc-release\pndf_run_qem.exe `
    --input out\isotropic_nxy.bin `
    --output out\isotropic_normal_262k.meshbin `
    --target 262144 `
    --mode normal

build\msvc-release\pndf_run_qem.exe `
    --input out\brush_nxy.bin `
    --output out\brush_normal_262k.meshbin `
    --target 262144 `
    --mode normal