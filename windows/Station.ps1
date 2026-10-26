param(
    [ValidateSet('Install','Start','Check','Update','Stop','Status','InstallAdvanced')][string]$Action='Start',
    [Nullable[int]]$Port,
    [string]$Proxy=''
)
$ErrorActionPreference='Stop'
try {
    $root=$PSScriptRoot
    $config=Get-Content -LiteralPath (Join-Path $root '配置\workstation.json') -Raw -Encoding UTF8 | ConvertFrom-Json
    if ($config.service_id -ne '16s-workstation' -or $config.schema -ne 1 -or $config.root_windows -ine $root) {
        throw '此目录不是本工作站的安装目录。'
    }
    $driver=Join-Path $config.program_windows 'windows\Workstation.ps1'
    if (!(Test-Path -LiteralPath $driver -PathType Leaf)) { throw '程序文件缺失，请重新下载源码并运行部署工作站.bat，选择当前工作站目录。' }
    $arguments=@{Action=$Action;InstallRoot=$root;Distro=$config.distro}
    if ($PSBoundParameters.ContainsKey('Port')) { $arguments.Port=[int]$Port }
    if ($PSBoundParameters.ContainsKey('Proxy')) { $arguments.Proxy=$Proxy }
    & $driver @arguments
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
} catch {
    Write-Host $_.Exception.Message -ForegroundColor Red
    exit 1
}
