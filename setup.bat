@echo off
REM One-time install. Needs Python 3.12 from python.org.
cd /d "%~dp0"
py -3.12 -m venv .venv || python -m venv .venv
.venv\Scripts\python -m pip install --upgrade pip
.venv\Scripts\python -m pip install -e .
echo.
echo Done. From your workspace folder (containing data\ and analysis\), run e.g.:
echo    "%~dp0kerncount.bat" init my_experiment
echo    "%~dp0kerncount.bat" run my_experiment
pause
