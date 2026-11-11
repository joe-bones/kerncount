@echo off
REM One-time install. No admin rights needed.
REM Uses uv (a single-file Python manager) to fetch Python 3.12 into your user folder,
REM so a system-wide Python install is not required.
setlocal
cd /d "%~dp0"

REM The environment lives at a short path: TensorFlow's files break Windows' 260-character
REM path limit when nested deep (e.g. under Downloads\...), and lifting that limit needs admin.
set "VENV=%LOCALAPPDATA%\kc"

REM 1. Find uv: next to this script, on PATH, or install it for this user only.
set "UV="
if exist "%~dp0uv.exe" set "UV=%~dp0uv.exe"
if not defined UV for %%U in (uv.exe) do if not "%%~$PATH:U"=="" set "UV=%%~$PATH:U"
if not defined UV if exist "%USERPROFILE%\.local\bin\uv.exe" set "UV=%USERPROFILE%\.local\bin\uv.exe"
if not defined UV (
    echo Installing uv for this user...
    powershell -NoProfile -ExecutionPolicy Bypass -Command "irm https://astral.sh/uv/install.ps1 | iex"
    if exist "%USERPROFILE%\.local\bin\uv.exe" set "UV=%USERPROFILE%\.local\bin\uv.exe"
)
if not defined UV (
    echo.
    echo Could not get uv. Download uv-x86_64-pc-windows-msvc.zip from
    echo    https://github.com/astral-sh/uv/releases/latest
    echo unzip it, put uv.exe in this folder ^(%~dp0^) and run setup.bat again.
    pause
    exit /b 1
)

REM 2. Python 3.12 + environment + kerncount and its dependencies.
"%UV%" venv --python 3.12 --seed --allow-existing "%VENV%" || goto :fail
"%UV%" pip install --python "%VENV%\Scripts\python.exe" -e . || goto :fail

REM 3. Smoke test: TensorFlow is the usual failure point.
"%VENV%\Scripts\python.exe" -c "import tensorflow, stardist, kerncount; print('kerncount OK, TensorFlow', tensorflow.__version__)" || goto :tf_fail

echo.
echo Done. Double-click kerncount-gui.bat to work in your browser, or from your workspace folder (containing data\ and analysis\), run e.g.:
echo    "%~dp0kerncount.bat" init my_experiment
echo    "%~dp0kerncount.bat" run my_experiment
pause
exit /b 0

:tf_fail
echo.
echo Installed, but TensorFlow will not load. If the message above mentions msvcp140 or a
echo missing DLL, this PC lacks the Microsoft Visual C++ Redistributable (needs IT/admin):
echo    https://aka.ms/vs/17/release/vc_redist.x64.exe
pause
exit /b 1

:fail
echo.
echo Install failed - see the message above. Behind a university proxy, try first:
echo    set UV_NATIVE_TLS=1
echo and run setup.bat again from the same window.
pause
exit /b 1
