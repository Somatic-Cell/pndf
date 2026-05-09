$root = "out_quality_degenerate_guard_actual_flip_reloc_unwrap"
$stage = "share_text_only"
$zip = "out_quality_degenerate_guard_actual_flip_reloc_unwrap.zip"

Remove-Item $stage -Recurse -Force -ErrorAction SilentlyContinue
New-Item -ItemType Directory $stage | Out-Null

$keepExt = @(".csv", ".txt", ".log", ".json")

Get-ChildItem $root -Recurse -File | Where-Object {
    $keepExt -contains $_.Extension.ToLower()
} | ForEach-Object {
    $rel = $_.FullName.Substring((Resolve-Path $root).Path.Length).TrimStart('\')
    $dst = Join-Path $stage $rel
    New-Item -ItemType Directory (Split-Path $dst) -Force | Out-Null
    Copy-Item $_.FullName $dst
}

Compress-Archive -Path "$stage\*" -DestinationPath $zip -Force