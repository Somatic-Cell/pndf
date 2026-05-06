param(
    [string]$Preset = "msvc-release"
)
cmake --preset $Preset
cmake --build --preset $Preset
ctest --test-dir "build/$Preset" --output-on-failure
