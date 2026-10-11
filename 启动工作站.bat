@echo off
chcp 65001 >nul
wsl.exe -d Ubuntu-24.04 --cd "%~dp0." bash "启动.sh"
pause
