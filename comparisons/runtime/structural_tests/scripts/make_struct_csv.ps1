param(
  [string]$Tag = 'r4w2'
)

$structRoot = Split-Path -Parent $PSScriptRoot
$rawDir = Join-Path $structRoot 'raw'
$processedDir = Join-Path $structRoot 'processed'

$juliaPath = Join-Path $rawDir ("struct_julia_{0}.txt" -f $Tag)
$rPath = Join-Path $rawDir ("struct_r_{0}.txt" -f $Tag)

$juliaCsvPath = Join-Path $processedDir ("struct_julia_{0}.csv" -f $Tag)
$rCsvPath = Join-Path $processedDir ("struct_r_{0}.csv" -f $Tag)
$combinedCsvPath = Join-Path $processedDir ("struct_combined_{0}.csv" -f $Tag)

function Parse-Lines($path, $impl) {
  $rows = @()
  Get-Content $path | ForEach-Object {
    $line = $_.Trim()
    if (-not $line.StartsWith("impl=$impl family=")) { return }
    $obj = [ordered]@{}
    foreach ($tok in ($line -split ' ')) {
      if ($tok -notmatch '=') { continue }
      $kv = $tok -split '=', 2
      $k = $kv[0]; $v = $kv[1]
      if ($v -match '^(true|false|TRUE|FALSE)$') { $obj[$k] = $v }
      elseif ($v -match '^[0-9]+$') { $obj[$k] = [int]$v }
      elseif ($v -match '^[0-9]+\.[0-9]+$') { $obj[$k] = [double]$v }
      else { $obj[$k] = $v }
    }
    $rows += [pscustomobject]$obj
  }
  return $rows
}

$j = Parse-Lines $juliaPath 'julia_struct'
$r = Parse-Lines $rPath 'r_struct'

$j | Export-Csv -NoTypeInformation -Path $juliaCsvPath
$r | Export-Csv -NoTypeInformation -Path $rCsvPath

$indexR = @{}
foreach ($row in $r) {
  $indexR["$($row.family)|$($row.n)"] = $row
}

$combined = @()
foreach ($jr in $j) {
  $key = "$($jr.family)|$($jr.n)"
  if (-not $indexR.ContainsKey($key)) { continue }
  $rr = $indexR[$key]
  $combined += [pscustomobject][ordered]@{
    family = $jr.family
    n = $jr.n
    julia_total_ms = [double]$jr.warm_total_median_ms
    r_total_ms = [double]$rr.warm_total_median_ms
    speedup_total_r_over_julia = [math]::Round(([double]$rr.warm_total_median_ms / [double]$jr.warm_total_median_ms), 4)
    julia_setup_ms = [double]$jr.warm_setup_median_ms
    r_setup_ms = [double]$rr.warm_setup_median_ms
    speedup_setup_r_over_julia = [math]::Round(([double]$rr.warm_setup_median_ms / [double]$jr.warm_setup_median_ms), 4)
    julia_identify_ms = [double]$jr.warm_identify_median_ms
    r_identify_ms = [double]$rr.warm_identify_median_ms
    speedup_identify_r_over_julia = [math]::Round(([double]$rr.warm_identify_median_ms / [double]$jr.warm_identify_median_ms), 4)
    julia_build_ms = [double]$jr.warm_build_median_ms
    julia_simplify_ms = [double]$jr.warm_simplify_median_ms
    julia_step4_ms = [double]$jr.warm_step4_median_ms
    julia_step5_ms = [double]$jr.warm_step5_median_ms
  }
}

$combined | Sort-Object family, n | Export-Csv -NoTypeInformation -Path $combinedCsvPath

Write-Output ("written={0}" -f $juliaCsvPath)
Write-Output ("written={0}" -f $rCsvPath)
Write-Output ("written={0}" -f $combinedCsvPath)
