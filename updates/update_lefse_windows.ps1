# One-time source update from 5fec64a to the LEfSe formatter fix.
# Windows downloads through its local proxy; WSL imports an authenticated Git bundle.
$ErrorActionPreference = 'Stop'
[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
$targetCommit = '759b947154b67283c2d8a98ec52083c4855cf237'
$expectedHash = 'f1aa0fa62b9e081c54e826e318ec3b9ac80ee3f4aabeb3a53b7b77751b75b151'
$bundleUrl = 'https://raw.githubusercontent.com/kangweihou522-cmyk/16s-/main/updates/lefse-759b947.bundle'
$repoPath = '/root/16S-workstation-v02'
$bundleFile = Join-Path $env:TEMP ('16s-update-' + [guid]::NewGuid().ToString('N') + '.bundle')

function Invoke-CheckedWsl {
    param([string[]]$NativeArguments)
    $ErrorActionPreference = 'Continue'
    & wsl.exe -d Ubuntu-24.04 -u root -- @NativeArguments
    if ($LASTEXITCODE -ne 0) { throw ('WSL command failed; exit code: ' + $LASTEXITCODE) }
}

try {
    $branch = (Invoke-CheckedWsl -NativeArguments @('git', '-C', $repoPath, 'symbolic-ref', '--short', 'HEAD') | Out-String).Trim()
    if ($branch -ne 'main') { throw 'The checkout is not on main; no source was changed.' }
    $changes = (Invoke-CheckedWsl -NativeArguments @('git', '-C', $repoPath, 'status', '--porcelain') | Out-String).Trim()
    if ($changes) { throw 'Local source changes exist; preserved without overwriting. Please report git status.' }
    $head = (Invoke-CheckedWsl -NativeArguments @('git', '-C', $repoPath, 'rev-parse', 'HEAD') | Out-String).Trim()
    if ($head -ne $targetCommit) {
        Write-Host 'Downloading the 8 KB source update using Windows...'
        Invoke-WebRequest -UseBasicParsing -Uri $bundleUrl -OutFile $bundleFile -Proxy 'http://127.0.0.1:10808' -TimeoutSec 30
        if ((Get-FileHash -LiteralPath $bundleFile -Algorithm SHA256).Hash.ToLowerInvariant() -ne $expectedHash) {
            throw 'Update checksum mismatch; no source was changed.'
        }
        $linuxBundle = (Invoke-CheckedWsl -NativeArguments @('wslpath', '-u', $bundleFile) | Out-String).Trim()
        Invoke-CheckedWsl -NativeArguments @('git', '-C', $repoPath, 'bundle', 'verify', '--quiet', $linuxBundle)
        Invoke-CheckedWsl -NativeArguments @('git', '-C', $repoPath, 'fetch', '--quiet', '--no-tags', $linuxBundle, 'refs/heads/main')
        Invoke-CheckedWsl -NativeArguments @('git', '-C', $repoPath, 'merge', '--quiet', '--ff-only', $targetCommit)
        $head = (Invoke-CheckedWsl -NativeArguments @('git', '-C', $repoPath, 'rev-parse', 'HEAD') | Out-String).Trim()
        if ($head -ne $targetCommit) { throw 'Source revision verification failed.' }
    }
    Write-Host 'Source revision verified: 759b947. Starting the workstation...'
} finally {
    if (Test-Path -LiteralPath $bundleFile) { Remove-Item -LiteralPath $bundleFile -Force }
}
Invoke-CheckedWsl -NativeArguments @('python3', ($repoPath + '/launch.py'), '--no-update', '--port', '8766')
