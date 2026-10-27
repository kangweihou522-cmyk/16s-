param([string]$Driver=(Join-Path $PSScriptRoot '../windows/Workstation.ps1'))
# Deterministic downloader orchestration checks; no network or Windows/WSL VM.
$ErrorActionPreference='Stop'
. $Driver -FunctionsOnly
$script:count=0
function Assert-True {
    param($Value,[string]$Message)
    if (!$Value) { throw $Message }
    $script:count++
}
function Assert-Throws {
    param([scriptblock]$Body,[string]$Message)
    $caught=$null
    try { & $Body } catch { $caught=$_ }
    Assert-True ($null -ne $caught) $Message
    return $caught
}

$temporary=Join-Path ([IO.Path]::GetTempPath()) ('16s-download-'+[Guid]::NewGuid().ToString('N'))
$null=New-Item -ItemType Directory -Path $temporary
try {
    $script:expectedProxy=New-Object Net.WebProxy('http://127.0.0.1:10808')
    function Get-DownloadProxy { param([string]$Root); return $script:expectedProxy }
    function Invoke-ArchiveTransfer {
        param([Uri]$Uri,[string]$Destination,$Proxy,[int]$RemainingMs)
        $script:attempts++
        Assert-True ($Uri.Scheme -eq 'https' -and $Uri.Host -eq 'codeload.github.com') 'Direct HTTPS codeload source selected'
        Assert-True ([object]::ReferenceEquals($Proxy,$script:expectedProxy)) 'Configured proxy passed to each transfer'
        Assert-True ($RemainingMs -gt 0 -and $RemainingMs -le 180000) 'Transfer receives a bounded remaining total budget'
        Assert-True (!(Test-Path -LiteralPath $Destination)) 'Each retry starts without prior partial bytes'
        Assert-True ($Destination -ne $script:archive) 'Transfer writes a staging file'
        Assert-True (!(Test-Path -LiteralPath $script:archive)) 'Archive is not published while transfer is incomplete'
        $script:budgets+=@($RemainingMs)
        if ($script:mode -eq 'fail' -or ($script:mode -eq 'retry' -and $script:attempts -eq 1)) {
            [IO.File]::WriteAllText($Destination,'incomplete first attempt')
            throw 'synthetic read timeout'
        }
        [IO.File]::WriteAllText($Destination,'complete archive bytes')
    }

    $root=Join-Path $temporary 'retry'
    $null=New-Item -ItemType Directory -Path $root
    $script:archive=Join-Path $root 'source.zip'
    $script:mode='retry';$script:attempts=0;$script:budgets=@()
    Save-SourceArchive $root $script:archive
    Assert-True ($script:attempts -eq 2) 'Read timeout retries the source download'
    Assert-True ([IO.File]::ReadAllText($script:archive) -ceq 'complete archive bytes') 'Successful retry publishes only its complete bytes'
    Assert-True (@(Get-ChildItem -LiteralPath $root -File).Count -eq 1) 'Successful retry removes staging files'
    Assert-True ($script:budgets[1] -le $script:budgets[0]) 'Retries share one total download budget'

    $root=Join-Path $temporary 'failure'
    $null=New-Item -ItemType Directory -Path $root
    $script:archive=Join-Path $root 'source.zip'
    $script:mode='fail';$script:attempts=0;$script:budgets=@()
    $failure=Assert-Throws { Save-SourceArchive $root $script:archive } 'Repeated read timeouts propagate a terminal error'
    Assert-True ($script:attempts -eq 3) 'Download stops after three attempts'
    Assert-True (!(Test-Path -LiteralPath $script:archive)) 'Failed download never publishes an archive'
    Assert-True (@(Get-ChildItem -LiteralPath $root -File).Count -eq 0) 'Failed download removes partial files'

    [IO.File]::WriteAllText($script:archive,'existing archive bytes')
    $script:attempts=0
    $null=Assert-Throws { Save-SourceArchive $root $script:archive } 'Existing archive cannot be overwritten'
    Assert-True ([IO.File]::ReadAllText($script:archive) -ceq 'existing archive bytes') 'Existing archive survives a failed download'
    Assert-True (@(Get-ChildItem -LiteralPath $root -File).Count -eq 1) 'Existing archive is the only file after failure'

    # Exercise the Update caller with a failed download. Mock host integration so
    # this runs on PowerShell for Linux as well as Windows without service access.
    & {
        $script:events=@()
        $script:Action='Update'
        $script:InstallRoot='C:\16s-download-test'
        $script:Distro='';$script:ExplicitPort=$false
        $resolvedRoot=[IO.Path]::GetFullPath($script:InstallRoot).TrimEnd('\','/')
        $script:testConfig=[pscustomobject]@{
            service_id='16s-workstation';schema=1;root_windows=$resolvedRoot
            distro='Ubuntu-24.04';wsl_user='existing_user';port=8766
        }
        function Get-WslExecutable { return 'synthetic-wsl.exe' }
        function Test-Path { param($LiteralPath,$PathType); return $true }
        function Get-Content { param($LiteralPath,[switch]$Raw,$Encoding); return ($script:testConfig | ConvertTo-Json -Compress) }
        function New-Item { param($ItemType,$Path,[switch]$Force) }
        function Select-Ubuntu { param($Requested); return 'Ubuntu-24.04' }
        function Get-WslUser { param($SelectedDistro,$WslUser); return 'existing_user' }
        function Get-ServiceState {
            param($Config,$Operation='probe')
            $script:events+=@($Operation)
            return @{state='running';health=@{busy=$false}}
        }
        function Save-SourceArchive {
            param($Root,$Archive)
            $script:events+=@('download')
            throw 'synthetic download failure'
        }
        function Install-Station { $script:events+=@('install') }
        $null=Assert-Throws { Invoke-Workstation } 'Update propagates download failure'
        Assert-True (($script:events -join ',') -eq 'probe,download') 'Update checks service then downloads without stopping or installing after download failure'
    }
} finally {
    Remove-Item -LiteralPath $temporary -Recurse -Force
}
"PASS: $script:count download checks (mock transfers and host integration; no network or Windows/WSL VM)"
