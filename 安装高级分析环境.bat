@echo off
chcp 65001 >nul
setlocal
if exist "%~dp0Station.ps1" (
    powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File "%~dp0Station.ps1" -Action InstallAdvanced %*
) else (
    powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File "%~dp0windows\Workstation.ps1" -Action InstallAdvanced %*
)
set "station_exit=%errorlevel%"
pause
exit /b %station_exit%
