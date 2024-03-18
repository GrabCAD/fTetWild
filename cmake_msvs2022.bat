mkdir _build
cd _build
rem For some confounded reason, running once doesn't set /MT properly.
cmake .. -G "Visual Studio 17 2022"
cmake .. -G "Visual Studio 17 2022"
pause
