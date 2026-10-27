[CmdletBinding()]
param(
    [ValidateSet('Install','Start','Check','Update','Stop','Status','InstallAdvanced')][string]$Action='Start',
    [string]$InstallRoot='',
    [string]$Distro='',
    [ValidateRange(1024,65535)][int]$Port=8766,
    [string]$Proxy='',
    [switch]$FunctionsOnly
)

function ConvertTo-NativeArgument {
    param([AllowEmptyString()][string]$Value)
    # WSL parses management options from the raw command line and retains quotes.
    # Keep flags/simple values bare; use CRT quoting only where required.
    if ($Value.Length -gt 0 -and $Value -notmatch '[\s"]') { return $Value }
    $escaped=$Value -replace '(\\*)"','$1$1\"'
    $escaped=$escaped -replace '(\\+)$','$1$1'
    return '"'+$escaped+'"'
}

function Get-WslExecutable {
    # Get-Command -CommandType Application can return multiple PATH matches.
    # Passing their .Source array to FileName joins them into an invalid filename.
    if ($env:OS -eq 'Windows_NT' -and $env:SystemRoot) {
        $folders=@('System32')
        if ([Environment]::Is64BitOperatingSystem -and ![Environment]::Is64BitProcess) { $folders=@('Sysnative','System32') }
        foreach ($folder in $folders) {
            $candidate=Join-Path $env:SystemRoot ($folder+'\wsl.exe')
            if (Test-Path -LiteralPath $candidate -PathType Leaf) { return $candidate }
        }
    }
    foreach ($command in @(Get-Command wsl.exe -CommandType Application -ErrorAction SilentlyContinue)) {
        if ($command.Path -and (Test-Path -LiteralPath $command.Path -PathType Leaf)) { return [string]$command.Path }
    }
    throw '未找到可执行的 wsl.exe。请检查 Windows 的 WSL 安装。'
}

function ConvertFrom-WslBytes {
    param([byte[]]$Bytes)
    # wsl --list emits UTF-16LE on some Windows versions; Linux tools use UTF-8.
    if ($Bytes.Length -ge 2 -and (($Bytes[0] -eq 255 -and $Bytes[1] -eq 254) -or $Bytes[1] -eq 0)) {
        return [Text.Encoding]::Unicode.GetString($Bytes).TrimStart([char]0xfeff)
    }
    return [Text.Encoding]::UTF8.GetString($Bytes).TrimStart([char]0xfeff)
}

function Invoke-Wsl {
    param([string[]]$WslArguments,[switch]$Capture,[switch]$AllowFailure,[string]$Step='WSL 调用')
    # WSL flags remain bare; arguments after --exec follow Windows CRT quoting.
    # Never merge stderr into paths/JSON: WSL startup warnings also use stderr.
    $info=New-Object Diagnostics.ProcessStartInfo
    $info.FileName=Get-WslExecutable
    $info.Arguments=($WslArguments | ForEach-Object { ConvertTo-NativeArgument $_ }) -join ' '
    $info.UseShellExecute=$false
    # Native WSL warnings otherwise use UTF-16 while Linux stdout/stderr use UTF-8.
    $info.EnvironmentVariables['WSL_UTF8']='1'
    $info.RedirectStandardOutput=[bool]$Capture
    $info.RedirectStandardError=[bool]$Capture
    $process=New-Object Diagnostics.Process
    $process.StartInfo=$info
    $stdout=$null; $stderr=$null; $errorText=''; $outputText=''
    try {
        try {
            if (!$process.Start()) { throw '无法启动 wsl.exe。' }
        } catch {
            $detail=$_.Exception.GetBaseException().Message
            if ($script:DiagnosticFile) {
                [pscustomobject]@{time=[DateTime]::Now.ToString('o');step=$Step;executable=$info.FileName;arguments=$WslArguments;start_error=$detail} | ConvertTo-Json -Compress | Add-Content -LiteralPath $script:DiagnosticFile -Encoding UTF8
            }
            throw "$Step 无法启动。可执行文件：$($info.FileName)。原因：$detail。诊断日志：$script:DiagnosticFile"
        }
        if ($Capture) {
            $stdout=New-Object IO.MemoryStream
            $stderr=New-Object IO.MemoryStream
            # Drain both pipes concurrently, including output larger than a pipe buffer.
            $outTask=$process.StandardOutput.BaseStream.CopyToAsync($stdout)
            $errTask=$process.StandardError.BaseStream.CopyToAsync($stderr)
        }
        $process.WaitForExit()
        $code=$process.ExitCode
        if ($Capture) {
            $null=$outTask.GetAwaiter().GetResult(); $null=$errTask.GetAwaiter().GetResult()
            $outputText=ConvertFrom-WslBytes $stdout.ToArray()
            $errorText=ConvertFrom-WslBytes $stderr.ToArray()
        }
    } finally {
        if ($stdout) { $stdout.Dispose() }; if ($stderr) { $stderr.Dispose() }; $process.Dispose()
    }
    if ($script:DiagnosticFile) {
        $record=[pscustomobject]@{time=[DateTime]::Now.ToString('o');step=$Step;executable=$info.FileName;command_line=$info.Arguments;arguments=$WslArguments;exit_code=$code;stderr=$errorText}
        $record | ConvertTo-Json -Compress | Add-Content -LiteralPath $script:DiagnosticFile -Encoding UTF8
    }
    if ($Capture -and $code -eq 0 -and $errorText.Trim() -and !$AllowFailure) {
        if ($errorText -match '(?i)localhost.*(?:代理|proxy)') {
            if (!$script:ProxyNoticeShown) {
                Write-Host 'WSL 提示 Windows localhost 代理未转发；日志已单独记录。更新由 Windows 使用系统代理下载。' -ForegroundColor Yellow
                $script:ProxyNoticeShown=$true
            }
        } else { Write-Host $errorText.Trim() -ForegroundColor Yellow }
    }
    if (!$AllowFailure -and $code -ne 0) {
        if ($Capture -and $errorText.Trim()) { Write-Host $errorText.Trim() -ForegroundColor Red }
        throw "$Step 失败（退出码 $code）。WSL 命令：$($info.Arguments)。诊断日志：$script:DiagnosticFile"
    }
    if ($Capture) { return [pscustomobject]@{ExitCode=$code;Text=$outputText;Stderr=$errorText} }
}

function Select-Ubuntu {
    param([string]$Requested)
    $result=Invoke-Wsl -WslArguments @('--list','--quiet') -Capture
    $names=@(($result.Text -replace "`0",'') -split '\r?\n' | ForEach-Object { $_.Trim() } | Where-Object { $_ })
    if ($Requested) {
        if ($names -notcontains $Requested) { throw "找不到 Ubuntu 发行版 $Requested。已安装：$($names -join ', ')。" }
        return $Requested
    }
    foreach ($preferred in @('Ubuntu-24.04','Ubuntu')) { if ($names -contains $preferred) { return $preferred } }
    $ubuntu=@($names | Where-Object { $_ -match '^Ubuntu(?:-|$)' })
    if ($ubuntu.Count -gt 0) { return $ubuntu[0] }
    throw '没有找到 Ubuntu。请先在 Windows 安装 WSL2 Ubuntu 24.04，然后再次运行部署。'
}

function Get-WslPath {
    param([string]$WindowsPath,[string]$SelectedDistro,[string]$WslUser)
    if ($WindowsPath -match '^\\\\wsl(?:\.localhost|\$)\\([^\\]+)\\(.*)$') {
        if ($Matches[1] -ne $SelectedDistro) { throw '源码位于另一 Ubuntu 发行版，请先将源码解压到普通 Windows 文件夹。' }
        return '/'+($Matches[2] -replace '\\','/')
    }
    $result=Invoke-Wsl -WslArguments @('-d',$SelectedDistro,'-u',$WslUser,'--exec','wslpath','-u',$WindowsPath.TrimEnd('\')) -Capture -Step '转换 Windows 路径'
    $path=$result.Text.Trim()
    if (!$path.StartsWith('/') -or $path -match '[\r\n\x00]') { throw "WSL 返回了无效路径，尚未调用分析脚本。Windows 路径：$WindowsPath。诊断日志：$script:DiagnosticFile" }
    return $path
}

function Select-InstallRoot {
    param([string]$Default='D:\16S分析工作站')
    Add-Type -AssemblyName System.Windows.Forms
    Add-Type -AssemblyName System.Drawing
    $form=New-Object Windows.Forms.Form
    $form.Text='部署 16S 分析工作站'
    $form.ClientSize=New-Object Drawing.Size(580,148)
    $form.StartPosition='CenterScreen'; $form.FormBorderStyle='FixedDialog'
    $form.MaximizeBox=$false; $form.MinimizeBox=$false
    $label=New-Object Windows.Forms.Label
    $label.Text='选择工作站文件夹（程序、原始数据、参考文件、分析结果均保存在这里）'
    $label.SetBounds(14,15,550,22)
    $text=New-Object Windows.Forms.TextBox
    $text.Text=$Default; $text.SetBounds(14,43,452,25)
    $browse=New-Object Windows.Forms.Button
    $browse.Text='浏览…'; $browse.SetBounds(476,42,88,27)
    $browse.Add_Click({
        $picker=New-Object Windows.Forms.FolderBrowserDialog
        try {
            $picker.Description='选择或新建 16S 工作站文件夹'
            $picker.ShowNewFolderButton=$true
            if (Test-Path -LiteralPath $text.Text -PathType Container) { $picker.SelectedPath=$text.Text }
            if ($picker.ShowDialog($form) -eq [Windows.Forms.DialogResult]::OK) { $text.Text=$picker.SelectedPath }
        } finally { $picker.Dispose() }
    })
    $ok=New-Object Windows.Forms.Button
    $ok.Text='部署'; $ok.SetBounds(374,99,90,29); $ok.DialogResult=[Windows.Forms.DialogResult]::OK
    $cancel=New-Object Windows.Forms.Button
    $cancel.Text='取消'; $cancel.SetBounds(474,99,90,29); $cancel.DialogResult=[Windows.Forms.DialogResult]::Cancel
    $form.Controls.AddRange(@($label,$text,$browse,$ok,$cancel))
    $form.AcceptButton=$ok; $form.CancelButton=$cancel
    try {
        if ($form.ShowDialog() -ne [Windows.Forms.DialogResult]::OK) { throw '已取消部署；没有修改现有工作站。' }
        if (!$text.Text.Trim()) { throw '请填写工作站文件夹。' }
        return $text.Text.Trim()
    } finally { $form.Dispose() }
}

function Get-WslUser {
    param([string]$SelectedDistro,[string]$ExistingUser)
    $native=@('-d',$SelectedDistro)
    if ($ExistingUser) { $native+=@('-u',$ExistingUser) }
    $native+=@('--exec','id','-un')
    $result=Invoke-Wsl -WslArguments $native -Capture -Step '识别现有 Ubuntu 用户'
    $user=$result.Text.Trim()
    if (!$user -or $user -match '[\s\x00]' -or $user.StartsWith('-')) { throw '无法识别现有 Ubuntu 用户，请检查 WSL 默认用户设置。' }
    return $user
}

function Invoke-Deploy {
    param($Config,[string]$Operation,[string[]]$ExtraArguments=@(),[switch]$Capture)
    $native=@('-d',$Config.distro,'-u',$Config.wsl_user,'--exec','python3',
        ($Config.program_wsl+'/windows/deploy.py'),$Operation,$Config.root_wsl)+$ExtraArguments
    return Invoke-Wsl -WslArguments $native -Capture:$Capture -Step ('工作站操作：'+$Operation)
}

function Get-ServiceState {
    param($Config,[ValidateSet('probe','wait','stop')][string]$Operation='probe')
    $result=Invoke-Wsl -WslArguments @('-d',$Config.distro,'-u',$Config.wsl_user,'--exec','python3',
        ($Config.program_wsl+'/windows/station_runtime.py'),$Operation,$Config.root_wsl,[string]$Config.port) -Capture -Step ('检查服务：'+$Operation)
    $state=$result.Text | ConvertFrom-Json
    if ($state.state -notin @('running','offline','occupied')) { throw '服务检查返回了未知状态，未启动或停止任何服务。' }
    return $state
}

function Assert-UpdateAllowed {
    param($Config)
    $state=Get-ServiceState $Config
    if ($state.state -eq 'occupied') { throw "端口 $($Config.port) 已被其他程序占用。请先自行关闭原有终端，或部署时使用 -Port 选择其他端口。" }
    if ($state.state -eq 'running' -and $state.health.busy) { throw '工作站正在分析；请等待当前分析完成后再更新或部署。' }
}

function Show-EnvironmentCheck {
    param($Config,[string]$PlanFile='')
    $extra=@()
    if ($PlanFile) { $extra=@('--plan',$PlanFile) }
    $result=Invoke-Deploy $Config 'check' -ExtraArguments $extra -Capture
    $report=$result.Text | ConvertFrom-Json
    Write-Host ('复用现有 Conda：'+$report.conda_wsl)
    foreach ($entry in $report.envs.PSObject.Properties) {
        $status='可用'; if (!$entry.Value.available) { $status='未安装（使用相应分析功能前需准备）' }
        Write-Host ('  '+$entry.Value.name+'：'+$status)
    }
    if ($report.optional_envs.advanced) {
        $advanced=$report.optional_envs.advanced
        if ($advanced.available) { Write-Host ('  '+$advanced.name+'：高级分析可用') }
        else { Write-Host ('  '+$advanced.name+'：可选，首次需要时双击安装高级分析环境.bat') }
    }
    foreach ($warning in $report.warnings) { Write-Warning $warning }
}

function Assert-WindowsService {
    param($Config,$ExpectedHealth)
    for ($i=0;$i -lt 5;$i++) {
        try {
            # Local health checks must bypass any Windows/system HTTP proxy.
            $request=[Net.HttpWebRequest]::Create("http://127.0.0.1:$($Config.port)/api/health")
            $request.Proxy=$null; $request.Timeout=4000; $request.ReadWriteTimeout=4000
            $response=$request.GetResponse()
            try {
                $reader=New-Object IO.StreamReader($response.GetResponseStream(),[Text.Encoding]::UTF8)
                try { $health=$reader.ReadToEnd() | ConvertFrom-Json } finally { $reader.Dispose() }
            } finally { $response.Dispose() }
            if ($health.service_id -ne '16s-workstation' -or $health.install_root -ne $Config.root_wsl -or
                $health.program_dir -ne $Config.program_wsl -or !$health.instance_id -or
                ($ExpectedHealth -and $health.instance_id -ne $ExpectedHealth.instance_id)) {
                throw 'Windows 本地端口返回的服务不属于此工作站当前版本。'
            }
            return
        } catch {
            if ($i -eq 4) { throw ('Ubuntu 服务已启动，但 Windows 检查未通过：'+$_.Exception.Message+' 请检查 WSL localhost 转发及端口设置。') }
        }
        Start-Sleep -Seconds 1
    }
}

function Start-StationConsole {
    param($Config)
    $native=@('-d',$Config.distro,'-u',$Config.wsl_user,'--exec','python3',
        ($Config.program_wsl+'/windows/deploy.py'),'start',$Config.root_wsl)
    $info=New-Object Diagnostics.ProcessStartInfo
    $info.FileName=Get-WslExecutable
    $info.Arguments=($native | ForEach-Object { ConvertTo-NativeArgument $_ }) -join ' '
    $info.UseShellExecute=$true
    $process=[Diagnostics.Process]::Start($info)
    if (!$process) { throw '无法打开工作站 Ubuntu 终端。' }
    $process.Dispose()
}

function Start-Station {
    param($Config)
    $state=Get-ServiceState $Config
    if ($state.state -eq 'occupied') { throw "端口 $($Config.port) 已被其他程序占用；没有停止其他程序。部署时可用 -Port 指定新端口。" }
    if ($state.state -eq 'running') {
        if ($state.health.program_dir -ne $Config.program_wsl) { throw '此安装目录的旧版本仍在运行；请在分析结束后运行更新工作站.bat。' }
        Write-Host '工作站已经运行，打开现有页面。'
    } else {
        Show-EnvironmentCheck $Config
        Start-StationConsole $Config
        $state=Get-ServiceState $Config 'wait'
        if ($state.state -ne 'running') { throw '工作站未能就绪，请查看 Ubuntu 终端及工作站日志。' }
    }
    Assert-WindowsService $Config $state.health
    $url="http://localhost:$($Config.port)"
    Write-Host ('工作站位置：'+$Config.root_windows)
    Write-Host ('分析结果：'+$Config.runs_windows)
    Write-Host "网页：$url。请保持新打开的 Ubuntu 工作站终端运行。"
    Start-Process -FilePath $url
}

function New-StationShortcut {
    param([string]$Root)
    try {
        $desktop=[Environment]::GetFolderPath('Desktop')
        $path=Join-Path $desktop '16S分析工作站.lnk'
        $shell=New-Object -ComObject WScript.Shell
        $shortcut=$shell.CreateShortcut($path)
        $arguments='/c ""'+(Join-Path $Root '启动工作站.bat')+'""'
        if (Test-Path -LiteralPath $path) {
            if ($shortcut.Arguments -ne $arguments -or $shortcut.TargetPath -ne $env:ComSpec) {
                Write-Warning '桌面已有同名的其他快捷方式，已保留。请从所选目录双击启动工作站.bat。'
            } else { Write-Host '桌面已有本工作站快捷方式，已保留。' }
            return
        }
        $shortcut.TargetPath=$env:ComSpec; $shortcut.Arguments=$arguments
        $shortcut.WorkingDirectory=$Root; $shortcut.Description='本地 16S 分析工作站'
        $shortcut.Save()
        Write-Host "桌面快捷方式已创建：$path"
    } catch { Write-Warning '无法创建桌面快捷方式；仍可从所选目录双击启动工作站.bat。' }
}

function Get-DownloadProxy {
    param([string]$Root)
    # Use the Windows proxy API (including PAC) by default. An optional local file
    # 配置/download.json accepts {"proxy":"http://127.0.0.1:10808"} or {"proxy":"direct"}.
    $value=$Proxy
    $file=Join-Path $Root '配置\download.json'
    if (!$value -and (Test-Path -LiteralPath $file -PathType Leaf)) {
        $settings=Get-Content -LiteralPath $file -Raw -Encoding UTF8 | ConvertFrom-Json
        $value=[string]$settings.proxy
    }
    if ($value -eq 'direct') { return $null }
    if ($value) {
        $uri=$null
        if (![Uri]::TryCreate($value,[UriKind]::Absolute,[ref]$uri) -or $uri.Scheme -notin @('http','https') -or $uri.UserInfo) {
            throw '下载代理格式无效，请在 配置/download.json 中填写 http(s)://主机:端口。'
        }
        return New-Object Net.WebProxy($uri)
    }
    $system=[Net.WebRequest]::GetSystemWebProxy()
    if ($system) { $system.Credentials=[Net.CredentialCache]::DefaultNetworkCredentials }
    return $system
}

function Invoke-ArchiveTransfer {
    param([Uri]$Uri,[string]$Destination,$Proxy,[int]$RemainingMs)
    if ($Uri.Scheme -ne 'https' -or $RemainingMs -le 0) { throw '源码下载地址或剩余时间无效。' }
    [Net.ServicePointManager]::SecurityProtocol=[Net.ServicePointManager]::SecurityProtocol -bor [Net.SecurityProtocolType]::Tls12
    $request=[Net.HttpWebRequest]::Create($Uri)
    $request.Timeout=[Math]::Min(30000,$RemainingMs)
    $request.ReadWriteTimeout=[Math]::Min(30000,$RemainingMs)
    $request.UserAgent='16s-workstation/2'; $request.Proxy=$Proxy
    $response=$null; $input=$null; $output=$null
    $watch=[Diagnostics.Stopwatch]::StartNew()
    try {
        $response=$request.GetResponse()
        if ($response.ResponseUri.Scheme -ne 'https' -or [int]$response.StatusCode -ne 200) {
            throw '源码下载未返回完整HTTPS响应。'
        }
        $expected=[long]$response.ContentLength
        if ($expected -gt 104857600) { throw '源码归档超过100MB限制，已拒绝。' }
        $input=$response.GetResponseStream()
        $output=[IO.File]::Open($Destination,[IO.FileMode]::CreateNew,[IO.FileAccess]::Write,[IO.FileShare]::None)
        $buffer=New-Object byte[] 65536; [long]$received=0
        while ($true) {
            $remaining=$RemainingMs-[int]$watch.ElapsedMilliseconds
            if ($remaining -le 0) { throw '源码下载已达到总时间上限。' }
            if ($input.CanTimeout) { $input.ReadTimeout=[Math]::Min(30000,$remaining) }
            $read=$input.Read($buffer,0,$buffer.Length)
            if ($read -le 0) { break }
            $received+=$read
            if ($received -gt 104857600) { throw '源码归档超过100MB限制，已拒绝。' }
            $output.Write($buffer,0,$read)
            $status=('{0:N1} MB 已下载' -f ($received/1MB))
            if ($expected -gt 0) { Write-Progress -Id 17 -Activity '下载16S工作站源码' -Status $status -PercentComplete ([Math]::Min(100,[int](100*$received/$expected))) }
            else { Write-Progress -Id 17 -Activity '下载16S工作站源码' -Status $status }
        }
        if ($received -le 0 -or ($expected -ge 0 -and $received -ne $expected)) { throw '源码下载不完整，未使用此文件。' }
    } finally {
        if ($output) { $output.Dispose() }
        if ($input) { $input.Dispose() }
        if ($response) { $response.Dispose() }
        $request.Abort(); $watch.Stop()
        Write-Progress -Id 17 -Activity '下载16S工作站源码' -Completed
    }
}

function Save-SourceArchive {
    param([string]$Root,[string]$Archive)
    if (Test-Path -LiteralPath $Archive) { throw '目标源码归档已存在，未覆盖。' }
    # Avoid the extra github.com redirect. Keep HTTPS and the configured Windows proxy.
    $uri=[Uri]'https://codeload.github.com/kangweihou522-cmyk/16s-/zip/refs/heads/main'
    $downloadProxy=Get-DownloadProxy $Root
    $watch=[Diagnostics.Stopwatch]::StartNew(); $lastError='网络不可用'; $attempts=0
    try {
        for ($attempt=1;$attempt -le 3;$attempt++) {
            $remaining=180000-[int]$watch.ElapsedMilliseconds
            if ($remaining -le 0) { break }
            $attempts=$attempt
            $partial=$Archive+'.'+[Guid]::NewGuid().ToString('N')+'.part'
            Write-Host "下载尝试 $attempt/3（总等待上限180秒；当前版本保留）…"
            try {
                Invoke-ArchiveTransfer -Uri $uri -Destination $partial -Proxy $downloadProxy -RemainingMs $remaining
                $length=(Get-Item -LiteralPath $partial -ErrorAction Stop).Length
                if ($length -le 0 -or $length -gt 104857600) { throw '源码归档为空或超过100MB限制。' }
                [IO.File]::Move($partial,$Archive)
                return
            } catch {
                $lastError=$_.Exception.GetBaseException().Message
                if ($lastError -match 'timed?\s*out|timeout|超时') { $lastError='连接或读取超时' }
                if ($attempt -lt 3 -and $watch.ElapsedMilliseconds -lt 180000) {
                    Write-Warning ("第 $attempt 次下载未完成："+$lastError+'。将重新下载，不使用残缺文件。')
                }
            } finally {
                if (Test-Path -LiteralPath $partial) { Remove-Item -LiteralPath $partial -Force -ErrorAction SilentlyContinue }
            }
        }
    } finally { $watch.Stop() }
    throw ("GitHub源码下载失败（已尝试 $attempts 次，总等待上限180秒）："+$lastError+'。当前程序、环境和分析结果保持原样；可直接使用启动工作站.bat。也可在浏览器下载完整源码，解压后运行部署工作站.bat并选择原安装目录进行离线源码更新。')
}

function Assert-SafeArchive {
    param([string]$Archive)
    Add-Type -AssemblyName System.IO.Compression.FileSystem
    $zip=[IO.Compression.ZipFile]::OpenRead($Archive)
    try {
        $seen=New-Object 'Collections.Generic.HashSet[string]' ([StringComparer]::OrdinalIgnoreCase)
        foreach ($entry in $zip.Entries) {
            $name=$entry.FullName
            if (!$name -or $name.Contains('\') -or $name.Contains(':') -or $name -match '^/|(^|/)\.\.?(/|$)|[\x00-\x1f<>"|?*]') {
                throw '源码 ZIP 含不安全路径，未解压。'
            }
            foreach ($part in ($name.TrimEnd('/') -split '/')) {
                if (!$part -or $part -match '[. ]$|^(?i:CON|PRN|AUX|NUL|COM[1-9]|LPT[1-9])(?:\.|$)') {
                    throw '源码 ZIP 含 Windows 无效路径，未解压。'
                }
            }
            $mode=($entry.ExternalAttributes -shr 16) -band 0xf000
            if ($mode -eq 0xa000 -or ($entry.ExternalAttributes -band 0x400) -ne 0) { throw '源码 ZIP 含符号链接，未解压。' }
            if (!$seen.Add($name.TrimEnd('/'))) { throw '源码 ZIP 含重复路径，未解压。' }
        }
    } finally { $zip.Dispose() }
}

function Get-UpdateSource {
    param([string]$Root)
    $folder=Join-Path $Root ('日志\updates\'+[Guid]::NewGuid().ToString('N'))
    New-Item -ItemType Directory -Path $folder -Force | Out-Null
    $archive=Join-Path $folder 'source.zip'
    Write-Host '由 Windows 下载 GitHub 最新源码，保留现有程序与分析结果…'
    Save-SourceArchive $Root $archive
    $hash=(Get-FileHash -LiteralPath $archive -Algorithm SHA256).Hash
    $hash | Set-Content -LiteralPath (Join-Path $folder 'source.sha256') -Encoding ASCII
    Write-Host ('源码 SHA256：'+$hash)
    Assert-SafeArchive $archive
    $expanded=Join-Path $folder 'source'
    Expand-Archive -LiteralPath $archive -DestinationPath $expanded
    $sources=@(Get-ChildItem -LiteralPath $expanded -Directory | Where-Object {
        Test-Path -LiteralPath (Join-Path $_.FullName 'windows\deploy.py') -PathType Leaf
    })
    if ($sources.Count -ne 1) { throw '下载的源码不完整，当前程序和数据保持原样。' }
    return $sources[0].FullName
}

function Install-Station {
    param([string]$Source,[string]$Root,[string]$SelectedDistro,[int]$SelectedPort,[string]$WslUser,$Existing)
    if ($Existing) { Assert-UpdateAllowed $Existing }
    $sourceWsl=Get-WslPath $Source $SelectedDistro $WslUser
    $rootWsl=Get-WslPath $Root $SelectedDistro $WslUser
    Write-Host "Ubuntu：$SelectedDistro；现有用户：$WslUser"
    Write-Host "源码：$sourceWsl；工作站：$rootWsl"
    $prepared=Invoke-Wsl -WslArguments @('-d',$SelectedDistro,'-u',$WslUser,'--exec','python3',($sourceWsl+'/windows/deploy.py'),
        'prepare',$sourceWsl,$rootWsl,$Root,$SelectedDistro,[string]$SelectedPort,$WslUser) -Capture -Step '准备候选版本（复用现有 Conda）'
    $proposal=$prepared.Text | ConvertFrom-Json
    if (!$proposal.config -or !$proposal.plan_file) { throw '候选版本准备结果不完整，未切换版本。' }
    Show-EnvironmentCheck $proposal.config $proposal.plan_file
    Assert-UpdateAllowed $proposal.config
    # Preparation and checks finish before the authenticated, busy-aware stop.
    if ($Existing) { $null=Get-ServiceState $Existing 'stop' }
    $activated=Invoke-Deploy $proposal.config 'activate' -ExtraArguments @($proposal.plan_file) -Capture
    $config=$activated.Text | ConvertFrom-Json
    try {
        Start-Station $config
        $null=Invoke-Deploy $config 'complete' -ExtraArguments @($proposal.plan_file) -Capture
    } catch {
        $startError=$_
        try {
            $null=Invoke-Deploy $proposal.config 'rollback' -ExtraArguments @($proposal.plan_file) -Capture
            Write-Warning '新版本启动未通过，已恢复部署前的配置和启动入口；原分析数据保持原样。'
        } catch {
            throw ($startError.Exception.Message+' 自动恢复未完成：'+$_.Exception.Message+' 请保留日志并从原安装目录检查服务状态。')
        }
        throw $startError
    }
    New-StationShortcut $Root
    Write-Host '部署完成。使用现有 Ubuntu 和 Conda，工作站保存在所选文件夹。'
}

function Invoke-Workstation {
    $ErrorActionPreference='Stop'
    [Console]::OutputEncoding=New-Object Text.UTF8Encoding($false)
    $OutputEncoding=[Console]::OutputEncoding
    Write-Host ('WSL：'+(Get-WslExecutable))
    $root=$InstallRoot
    if (!$root) {
        if ($Action -eq 'Install') { $root=Select-InstallRoot }
        else { $root='D:\16S分析工作站' }
    }
    if ($root -match '[\r\n\x00]' -or $root -notmatch '^[A-Za-z]:[\\/]') { throw '请选择本机 Windows 磁盘上的绝对路径，例如 D:\16S分析工作站。' }
    $root=[IO.Path]::GetFullPath($root).TrimEnd('\','/')
    if ($root -eq [IO.Path]::GetPathRoot($root).TrimEnd('\','/')) { throw '请使用独立工作站文件夹，不要直接选择磁盘根目录。' }
    if (!(Test-Path -LiteralPath ([IO.Path]::GetPathRoot($root)))) { throw '没有找到目标磁盘，请用 -InstallRoot 指定其他工作站文件夹。' }
    $file=Join-Path $root '配置\workstation.json'
    $existing=$null
    if (Test-Path -LiteralPath $file -PathType Leaf) {
        $existing=Get-Content -LiteralPath $file -Raw -Encoding UTF8 | ConvertFrom-Json
        if ($existing.service_id -ne '16s-workstation' -or $existing.schema -ne 1 -or $existing.root_windows -ine $root) {
            throw '目标目录包含其他或不匹配的配置，未覆盖，请选择独立安装目录。'
        }
    }
    if ($Action -ne 'Install' -and !$existing) { throw '尚未部署。请先运行部署工作站.bat，然后使用安装目录或桌面上的启动入口。' }
    $requested=$Distro
    if (!$requested -and $existing) { $requested=$existing.distro }
    if (!$requested -and $env:WSL_DISTRO) { $requested=$env:WSL_DISTRO }
    $selected=Select-Ubuntu $requested
    if ($existing -and $selected -ne $existing.distro) { throw '已安装工作站沿用原 Ubuntu 发行版；请不要更改 -Distro。' }
    $priorUser=''; if ($existing) { $priorUser=$existing.wsl_user }
    $wslUser=Get-WslUser $selected $priorUser
    $selectedPort=$Port
    if (!$script:ExplicitPort -and $existing) { $selectedPort=$existing.port }
    if ($existing -and $Action -notin @('Install','Update') -and $selectedPort -ne $existing.port) { throw '修改端口请重新运行部署工作站.bat -InstallRoot 路径 -Port 新端口。' }
    $hasher=[Security.Cryptography.SHA256]::Create()
    try { $key=([BitConverter]::ToString($hasher.ComputeHash([Text.Encoding]::UTF8.GetBytes($root.ToLowerInvariant())))).Replace('-','').Substring(0,16) } finally { $hasher.Dispose() }
    $mutex=New-Object Threading.Mutex($false,('Local\16sWorkstation-'+$key))
    if (!$mutex.WaitOne(0)) { $mutex.Dispose(); throw '工作站正在部署、启动或停止，请等待当前操作完成。' }
    try {
        $logDir=Join-Path $root '日志'
        New-Item -ItemType Directory -Path $logDir -Force | Out-Null
        $script:DiagnosticFile=Join-Path $logDir 'wsl-launcher.log'
        if ($Action -in @('Install','Update')) {
            # Busy checks precede network downloads and candidate preparation.
            if ($existing) { Assert-UpdateAllowed $existing }
            $source=[IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
            if ($Action -eq 'Update') { $source=Get-UpdateSource $root }
            Install-Station $source $root $selected $selectedPort $wslUser $existing
        } elseif ($Action -eq 'InstallAdvanced') {
            Assert-UpdateAllowed $existing
            $before=Get-ServiceState $existing
            if ($before.state -eq 'running') { $null=Get-ServiceState $existing 'stop' }
            try { Invoke-Deploy $existing 'install-advanced' }
            finally { if ($before.state -eq 'running') { Start-Station $existing } }
        } elseif ($Action -eq 'Check') {
            Show-EnvironmentCheck $existing
            $state=Get-ServiceState $existing
            Write-Host ('服务状态：'+$state.state)
            if ($state.state -eq 'running') { Write-Host ('正在分析：'+[bool]$state.health.busy) }
        } elseif ($Action -eq 'Status') {
            $state=Get-ServiceState $existing
            Write-Host ('服务状态：'+$state.state)
            if ($state.state -eq 'running') { Write-Host ('正在分析：'+[bool]$state.health.busy) }
        } elseif ($Action -eq 'Stop') {
            $null=Get-ServiceState $existing 'stop'
            Write-Host '工作站服务已停止（或原本未运行）。'
        } else { Start-Station $existing }
    } finally { $mutex.ReleaseMutex(); $mutex.Dispose() }
}

$script:ExplicitPort=$PSBoundParameters.ContainsKey('Port')
if (!$FunctionsOnly) {
    try { Invoke-Workstation; exit 0 }
    catch { Write-Host $_.Exception.Message -ForegroundColor Red; exit 1 }
}
