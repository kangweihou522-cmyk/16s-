# Windows downloads the checked source bundle; WSL imports it without network.
$ErrorActionPreference = 'Stop'
[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
$targetCommit = '7febb3ecee31ddffc6412cc91202b8c4add04162'
$expectedHash = '9b3b3c3f4cd2b54bf894864e1b11c5ade31401340c5c870290af9f9009add0cc'
$bundleUrl = 'https://raw.githubusercontent.com/kangweihou522-cmyk/16s-/main/updates/workstation-7febb3e.bundle'
$repoPath = '/root/16S-workstation-v02'
$fixedBlob = '2d9a9dc77f0396de12094f288304846928dc06b4'
$originalBlob = 'd01f90d35d6aa8812d86bb9ca1f44d3469d49c84'
$bundleFile = Join-Path $env:TEMP ('16s-update-' + [guid]::NewGuid().ToString('N') + '.bundle')

function Invoke-CheckedWsl {
    param([string[]]$NativeArguments)
    $ErrorActionPreference = 'Continue'
    & wsl.exe -d Ubuntu-24.04 -u root --exec @NativeArguments
    if ($LASTEXITCODE -ne 0) { throw ('WSL command failed; exit code: ' + $LASTEXITCODE) }
}
function Read-Wsl {
    param([string[]]$NativeArguments)
    return (Invoke-CheckedWsl -NativeArguments $NativeArguments | Out-String).Trim()
}

try {
    if ((Read-Wsl @('git', '-C', $repoPath, 'symbolic-ref', '--short', 'HEAD')) -ne 'main') {
        throw 'The checkout is not on main; no source was changed.'
    }
    $changes = Read-Wsl @('git', '-C', $repoPath, 'status', '--porcelain')
    $knownHotfix = $false
    if ($changes) {
        $workBlob = Read-Wsl @('git', '-C', $repoPath, 'hash-object', '--', ($repoPath + '/lefse_analysis.py'))
        $indexBlob = Read-Wsl @('git', '-C', $repoPath, 'rev-parse', ':lefse_analysis.py')
        if ($changes -notmatch '^M\s+lefse_analysis\.py$' -or $workBlob -ne $fixedBlob -or
            $indexBlob -notin @($originalBlob, $fixedBlob)) {
            throw 'Other local source changes exist; preserved without overwriting. Please report git status.'
        }
        $knownHotfix = $true
    }
    $head = Read-Wsl @('git', '-C', $repoPath, 'rev-parse', 'HEAD')
    if ($head -ne $targetCommit) {
        Write-Host 'Downloading the checked source update using Windows...'
        Invoke-WebRequest -UseBasicParsing -Uri $bundleUrl -OutFile $bundleFile -Proxy 'http://127.0.0.1:10808' -TimeoutSec 30
        if ((Get-FileHash -LiteralPath $bundleFile -Algorithm SHA256).Hash.ToLowerInvariant() -ne $expectedHash) {
            throw 'Update checksum mismatch; no source was changed.'
        }
        $linuxBundle = Read-Wsl @('wslpath', '-u', $bundleFile.Replace('\', '/'))
        Invoke-CheckedWsl @('git', '-C', $repoPath, 'bundle', 'verify', '--quiet', $linuxBundle)
        Invoke-CheckedWsl @('git', '-C', $repoPath, 'fetch', '--quiet', '--no-tags', $linuxBundle, 'refs/heads/main')
        Invoke-CheckedWsl @('git', '-C', $repoPath, 'merge-base', '--is-ancestor', $head, $targetCommit)
        if ((Read-Wsl @('git', '-C', $repoPath, 'status', '--porcelain')) -ne $changes) {
            throw 'Source changed during download; no local edits were replaced.'
        }
        $hotfixBackup = ''
        if ($knownHotfix) {
            if ((Read-Wsl @('git', '-C', $repoPath, 'hash-object', '--', ($repoPath + '/lefse_analysis.py'))) -ne $fixedBlob -or
                (Read-Wsl @('git', '-C', $repoPath, 'rev-parse', ':lefse_analysis.py')) -ne $indexBlob) {
                throw 'The local hotfix changed during download; preserved without overwriting.'
            }
            $hotfixBackup = Read-Wsl @('git', '-C', $repoPath, 'stash', 'create', '16S known LEfSe hotfix before update')
            if ($hotfixBackup -notmatch '^[0-9a-f]{40}$') { throw 'Hotfix backup failed; no source was replaced.' }
            Invoke-CheckedWsl @('git', '-C', $repoPath, 'stash', 'store', '-m', '16S known LEfSe hotfix before update', $hotfixBackup)
            Invoke-CheckedWsl @('git', '-C', $repoPath, 'restore', '--source=HEAD', '--staged', '--worktree', '--', 'lefse_analysis.py')
        }
        try {
            Invoke-CheckedWsl @('git', '-C', $repoPath, 'merge', '--quiet', '--ff-only', $targetCommit)
        } catch {
            if ($hotfixBackup) {
                Invoke-CheckedWsl @('git', '-C', $repoPath, 'restore', ('--source=' + $hotfixBackup), '--worktree', '--', 'lefse_analysis.py')
                Invoke-CheckedWsl @('git', '-C', $repoPath, 'restore', ('--source=' + $hotfixBackup + '^2'), '--staged', '--', 'lefse_analysis.py')
            }
            throw
        }
        if ((Read-Wsl @('git', '-C', $repoPath, 'rev-parse', 'HEAD')) -ne $targetCommit) {
            throw 'Source revision verification failed.'
        }
    }
    Write-Host 'Source revision verified: 7febb3e. Starting the workstation...'
} finally {
    if (Test-Path -LiteralPath $bundleFile) { Remove-Item -LiteralPath $bundleFile -Force }
}
Invoke-CheckedWsl @('python3', ($repoPath + '/launch.py'), '--no-update', '--port', '8766')
