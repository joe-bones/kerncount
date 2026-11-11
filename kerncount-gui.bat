@echo off
REM Double-click to open kerncount in your browser. Keep this window open while you work.
"%LOCALAPPDATA%\kc\Scripts\kerncount.exe" gui %*
if errorlevel 1 pause
