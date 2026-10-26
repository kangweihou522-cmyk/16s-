param([string]$Driver=(Join-Path $PSScriptRoot '../windows/Workstation.ps1'))
# PowerShell logic and Linux executable-shim checks; not a Windows/WSL VM test.
$ErrorActionPreference='Stop'
. $Driver -FunctionsOnly
$script:count=0
function Assert-True { param($Value,[string]$Message); if (!$Value) { throw $Message }; $script:count++ }
function Assert-Throws { param([scriptblock]$Body,[string]$Message); $caught=$false; try { & $Body } catch { $caught=$true }; Assert-True $caught $Message }
Assert-True ((ConvertTo-NativeArgument '--list') -ceq '--list') 'Bare management flag'
Assert-True ((ConvertTo-NativeArgument '') -ceq '""') 'Empty argument'
Assert-True ((ConvertTo-NativeArgument 'a"b') -ceq '"a\"b"') 'Embedded quote'
Assert-True ((ConvertTo-NativeArgument 'D:\16S data\') -ceq '"D:\16S data\\"') 'Trailing backslash'
Assert-True ((ConvertTo-NativeArgument '16S 分析') -ceq '"16S 分析"') 'Unicode and space'
Assert-True ((ConvertFrom-WslBytes ([Text.Encoding]::Unicode.GetBytes("Ubuntu-24.04`n"))) -ceq "Ubuntu-24.04`n") 'Legacy UTF16 output'

$temporary=Join-Path ([IO.Path]::GetTempPath()) ('16s-launcher-'+[Guid]::NewGuid().ToString('N'))
$null=New-Item -ItemType Directory -Path $temporary
try {
    if ($env:OS -ne 'Windows_NT') {
        $shim=Join-Path $temporary 'wsl.exe'
        $body=@'
#!/usr/bin/env python3
import json, os, sys
args=sys.argv[1:]
if args[0]=='--list':
    sys.stdout.buffer.write('Ubuntu-24.04\nUbuntu\n'.encode('utf-16-le'))
    sys.stderr.write('wsl: localhost 代理提示\n')
elif args[0]=='echo-argv':
    print(json.dumps(args[1:],ensure_ascii=False))
elif args[0]=='encoding':
    print(os.environ.get('WSL_UTF8',''))
elif args[0]=='large-pipes':
    sys.stdout.write('o'*200000);sys.stdout.flush();sys.stderr.write('e'*200000)
else:
    sys.stderr.write('cannot execute binary\n');sys.exit(126)
'@
        [IO.File]::WriteAllText($shim,$body.Replace("`r`n","`n"),(New-Object Text.UTF8Encoding($false)))
        & chmod +x $shim
        $second=Join-Path $temporary 'second'; $null=New-Item -ItemType Directory $second
        Copy-Item -LiteralPath $shim -Destination (Join-Path $second 'wsl.exe')
        $oldPath=$env:PATH
        try {
            $env:PATH=$temporary+[IO.Path]::PathSeparator+$second+[IO.Path]::PathSeparator+$oldPath
            Assert-True ((Get-WslExecutable) -ceq $shim) 'Exactly one executable selected from multiple PATH matches'
            $previous=$env:WSL_UTF8
            $result=Invoke-Wsl -WslArguments @('encoding') -Capture
            Assert-True ($result.Text.Trim() -eq '1' -and $env:WSL_UTF8 -ceq $previous) 'UTF8 set only for child'
            $result=Invoke-Wsl -WslArguments @('--list','--quiet') -Capture
            Assert-True (@($result).Count -eq 1 -and $result.Text -ceq "Ubuntu-24.04`nUbuntu`n") 'UTF16 stdout remains separate from warning'
            Assert-True ($result.Stderr -match '代理') 'UTF8 stderr decoded'
            $wanted=@('', '16S 分析', 'a"b', 'D:\16S data\', "O'Brien", '$(touch WRONG)', '--')
            $result=Invoke-Wsl -WslArguments (@('echo-argv')+$wanted) -Capture
            $actual=@($result.Text | ConvertFrom-Json)
            Assert-True ($actual.Count -eq $wanted.Count) 'Native argument count'
            for ($i=0;$i -lt $wanted.Count;$i++) { Assert-True ($actual[$i] -ceq $wanted[$i]) "Native argument $i preserved" }
            $result=Invoke-Wsl -WslArguments @('large-pipes') -Capture -AllowFailure
            Assert-True ($result.Text.Length -eq 200000 -and $result.Stderr.Length -eq 200000) 'Both pipes drain concurrently'
            $result=Invoke-Wsl -WslArguments @('fail') -Capture -AllowFailure
            Assert-True ($result.ExitCode -eq 126 -and $result.Stderr -match 'cannot execute') 'Native exit and stderr retained'
            Assert-Throws { Invoke-Wsl -WslArguments @('fail') -Capture } 'Native nonzero rejected'
        } finally { $env:PATH=$oldPath }
    }

    function Invoke-Wsl {
        param($WslArguments,[switch]$Capture,[switch]$AllowFailure,[string]$Step)
        $script:lastArgs=@($WslArguments)
        $text=''
        if ($WslArguments -contains '--list') { $text="Ubuntu-24.04`0`nUbuntu`0" }
        elseif ($WslArguments -contains 'wslpath') { $text='/mnt/d/16S 分析' }
        elseif ($WslArguments -contains 'id') { $text='existing_user' }
        if ($script:badPath) { $text="warning`n/mnt/d/root" }
        return [pscustomobject]@{Text=$text;Stderr='proxy warning';ExitCode=0}
    }
    Assert-True ((Select-Ubuntu '') -eq 'Ubuntu-24.04') 'Ubuntu 24.04 preferred'
    Assert-True ((Select-Ubuntu 'Ubuntu') -eq 'Ubuntu') 'Explicit distro preserved'
    Assert-Throws { Select-Ubuntu 'Missing' } 'Missing distro rejected'
    Assert-True ((Get-WslUser 'Ubuntu-24.04' '') -eq 'existing_user') 'Default Ubuntu user retained'
    Assert-True ($script:lastArgs -notcontains '-u') 'First user discovery does not force root'
    $null=Get-WslUser 'Ubuntu-24.04' 'existing_user'
    Assert-True ($script:lastArgs[3] -eq 'existing_user') 'Existing deployment user retained'
    $path="D:\16S 分析\O'Brien"
    $null=Get-WslPath $path 'Ubuntu-24.04' 'existing_user'
    Assert-True ($script:lastArgs[-1] -ceq $path) 'Windows path remains one argument'
    $script:badPath=$true
    Assert-Throws { Get-WslPath $path 'Ubuntu-24.04' 'existing_user' } 'Contaminated path rejected'
    $script:badPath=$false

    Add-Type -AssemblyName System.IO.Compression.FileSystem
    foreach ($name in @('../escape','/absolute','main/../escape','main\escape','main/C:ads','main/CON.txt','main/trailing.','main/link')) {
        $archive=Join-Path $temporary ([Guid]::NewGuid().ToString('N')+'.zip')
        $zip=[IO.Compression.ZipFile]::Open($archive,[IO.Compression.ZipArchiveMode]::Create)
        try { $entry=$zip.CreateEntry($name); if ($name -eq 'main/link') { $entry.ExternalAttributes=[int](-1577058304) } } finally { $zip.Dispose() }
        Assert-Throws { Assert-SafeArchive $archive } ('Unsafe ZIP rejected: '+$name)
    }
    $archive=Join-Path $temporary 'safe.zip'
    $zip=[IO.Compression.ZipFile]::Open($archive,[IO.Compression.ZipArchiveMode]::Create)
    try { $null=$zip.CreateEntry('16s--main/windows/Workstation.ps1') } finally { $zip.Dispose() }
    Assert-SafeArchive $archive
    Assert-True $true 'Valid archive accepted'

    $config=[pscustomobject]@{distro='Ubuntu-24.04';wsl_user='existing_user';program_wsl='/mnt/d/16S 分析/程序/version';root_wsl='/mnt/d/16S 分析';root_windows='D:\16S 分析';runs_windows='D:\16S 分析\分析结果';port=8766}
    function Get-ServiceState { param($Config,$Operation='probe'); $script:events+=@($Operation); if ($Operation -eq 'wait') { return @{state='running';health=@{program_dir=$Config.program_wsl;instance_id='test'}} }; return $script:state }
    function Show-EnvironmentCheck { param($Config,$PlanFile); $script:events+=@('check') }
    function Start-StationConsole { param($Config); $script:events+=@('console') }
    function Assert-WindowsService { param($Config,$ExpectedHealth); $script:events+=@('windows-health') }
    function Start-Process { param($FilePath); $script:events+=@('browser') }
    $script:state=@{state='occupied'}; $script:events=@()
    Assert-Throws { Start-Station $config } 'Foreign service blocks start'
    Assert-True (($script:events -join ',') -eq 'probe') 'Foreign service receives no stop or browser'
    $script:state=@{state='running';health=@{program_dir=$config.program_wsl;busy=$true}}; $script:events=@()
    Start-Station $config
    Assert-True (($script:events -join ',') -eq 'probe,windows-health,browser') 'Existing busy owned service opens current page only'
    Assert-Throws { Assert-UpdateAllowed $config } 'Busy owned service blocks update'
    $script:state=@{state='offline'}; $script:events=@()
    Start-Station $config
    Assert-True (($script:events -join ',') -eq 'probe,check,console,wait,windows-health,browser') 'Health validated before browser'

    function Get-WslPath { param($WindowsPath,$SelectedDistro,$WslUser); return '/mnt/d/source' }
    function Invoke-Wsl { param($WslArguments,[switch]$Capture,[string]$Step); $script:events+=@('prepare'); if ($script:prepareFail) { throw 'prepare failed' }; return @{Text=(@{config=$config;plan_file='/plan'} | ConvertTo-Json -Compress -Depth 4)} }
    function Show-EnvironmentCheck { param($Config,$PlanFile); $script:events+=@('check'); if ($script:checkFail) { throw 'check failed' } }
    function Invoke-Deploy { param($Config,$Operation,$ExtraArguments,[switch]$Capture); $script:events+=@($Operation); return @{Text=($config | ConvertTo-Json -Compress)} }
    function New-StationShortcut { param($Root); $script:events+=@('shortcut') }
    function Start-Station { param($Config); $script:events+=@('start'); if ($script:startFail) { throw 'startup failed' } }
    $script:events=@();$script:state=@{state='offline'}
    Install-Station 'D:\source' 'D:\16S 分析' 'Ubuntu-24.04' 8766 'existing_user' $config
    Assert-True (($script:events -join ',') -eq 'probe,prepare,check,probe,stop,activate,start,complete,shortcut') 'Prepare and candidate checks precede stopping old service'
    $script:events=@();$script:checkFail=$true
    Assert-Throws { Install-Station 'D:\source' 'D:\16S 分析' 'Ubuntu-24.04' 8766 'existing_user' $config } 'Candidate failure propagated'
    Assert-True ($script:events -notcontains 'stop' -and $script:events -notcontains 'activate') 'Candidate failure leaves service and config intact'
    $script:checkFail=$false;$script:events=@();$script:startFail=$true
    Assert-Throws { Install-Station 'D:\source' 'D:\16S 分析' 'Ubuntu-24.04' 8766 'existing_user' $config } 'Startup failure propagated'
    Assert-True ($script:events -contains 'rollback' -and $script:events -notcontains 'complete' -and $script:events -notcontains 'shortcut') 'Startup failure restores old configuration before shortcut creation'
    $script:startFail=$false;$script:events=@();$script:state=@{state='running';health=@{busy=$true}}
    Assert-Throws { Install-Station 'D:\source' 'D:\16S 分析' 'Ubuntu-24.04' 8766 'existing_user' $config } 'Busy deployment rejected'
    Assert-True (($script:events -join ',') -eq 'probe') 'Busy deployment rejected before prepare'
} finally { Remove-Item -LiteralPath $temporary -Recurse -Force }
"PASS: $script:count launcher checks (PowerShell logic and Linux native shim; no Windows/WSL desktop VM)"
