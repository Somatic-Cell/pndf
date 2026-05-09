$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$root = Resolve-Path (Join-Path $PSScriptRoot "..")
$cmakePath = Join-Path $root "CMakeLists.txt"
if (-not (Test-Path $cmakePath)) {
    throw "CMakeLists.txt not found: $cmakePath"
}

$pmpDir = Join-Path $root "extern\pmp-library"
if (-not (Test-Path $pmpDir)) {
    Write-Warning "extern\pmp-library was not found. Add it before configuring CMake:"
    Write-Warning "  git submodule add https://github.com/pmp-library/pmp-library extern/pmp-library"
}

$content = Get-Content -Raw -Path $cmakePath
$markerBegin = "# ---- pndf pnm_aniso_he / PMP begin ----"
$markerEnd   = "# ---- pndf pnm_aniso_he / PMP end ----"

$block = @"
$markerBegin
set(PMP_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(PMP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(PMP_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(PMP_BUILD_VIEWERS OFF CACHE BOOL "" FORCE)
set(PMP_INSTALL OFF CACHE BOOL "" FORCE)
set(PMP_STRICT_COMPILATION OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
if (EXISTS "`${CMAKE_CURRENT_SOURCE_DIR}/extern/pmp-library/CMakeLists.txt")
  add_subdirectory(extern/pmp-library EXCLUDE_FROM_ALL)
else()
  message(FATAL_ERROR "extern/pmp-library is missing. Run: git submodule add https://github.com/pmp-library/pmp-library extern/pmp-library")
endif()
add_executable(pndf_run_pnm_aniso_he apps/pndf_run_pnm_aniso_he.cpp)
target_link_libraries(pndf_run_pnm_aniso_he PRIVATE pndf_core pmp)
$markerEnd
"@

if ($content.Contains($markerBegin)) {
    $pattern = [regex]::Escape($markerBegin) + "(?s).*?" + [regex]::Escape($markerEnd)
    $content = [regex]::Replace($content, $pattern, [System.Text.RegularExpressions.MatchEvaluator]{ param($m) $block })
} else {
    if (-not $content.EndsWith("`n")) { $content += "`n" }
    $content += "`n" + $block + "`n"
}

Set-Content -Path $cmakePath -Value $content -NoNewline
Write-Host "Registered pndf_run_pnm_aniso_he and PMP in CMakeLists.txt"
