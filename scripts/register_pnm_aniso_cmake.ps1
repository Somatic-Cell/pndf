param(
    [string]$RepoRoot = ""
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

if ([string]::IsNullOrWhiteSpace($RepoRoot)) {
    $ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
    $RepoRoot = Resolve-Path (Join-Path $ScriptDir "..")
} else {
    $RepoRoot = Resolve-Path $RepoRoot
}

$cmakePath = Join-Path $RepoRoot "CMakeLists.txt"
if (!(Test-Path $cmakePath)) {
    throw "CMakeLists.txt not found: $cmakePath"
}

$text = Get-Content -Raw -Path $cmakePath

if ($text -match "add_executable\s*\(\s*pndf_run_pnm_aniso\b") {
    Write-Host "pndf_run_pnm_aniso target already exists in CMakeLists.txt"
    exit 0
}

$snippet = " add_executable(pndf_run_pnm_aniso apps/pndf_run_pnm_aniso.cpp) target_link_libraries(pndf_run_pnm_aniso PRIVATE pndf_core)"

# Prefer inserting after pndf_eval_mesh, because that target is present in both
# the original one-line CMakeLists.txt and most reformatted variants.
$anchor = "target_link_libraries(pndf_eval_mesh PRIVATE pndf_core)"
if ($text.Contains($anchor)) {
    $text = $text.Replace($anchor, $anchor + $snippet)
}
elseif ($text.Contains("if (PNDF_BUILD_TESTS)")) {
    $text = $text.Replace("if (PNDF_BUILD_TESTS)", $snippet + " if (PNDF_BUILD_TESTS)")
}
elseif ($text.Contains("if(PNDF_BUILD_TESTS)")) {
    $text = $text.Replace("if(PNDF_BUILD_TESTS)", $snippet + " if(PNDF_BUILD_TESTS)")
}
else {
    # Last-resort fallback: append the target. This is still valid CMake as long
    # as pndf_core is already defined above.
    $text = $text.TrimEnd() + $snippet + [Environment]::NewLine
}

Set-Content -Path $cmakePath -Value $text -NoNewline
Write-Host "Registered pndf_run_pnm_aniso in $cmakePath"
