SET ROPTS=--no-save --no-environ --no-init-file --no-restore --no-Rconsole

:: Where I got the idea to check if Rscript is already running --- https://stackoverflow.com/a/1329790
tasklist /fi "ImageName eq Rscript.exe" /fo csv 2>NUL | find /I "Rscript.exe">NUL
IF "%ERRORLEVEL%"=="0" (start http://127.0.0.1:7777 & R-4.5.0\bin\x64\Rscript.exe runApp_NPPS.R 1> log_NPPS.log 2>&1) ELSE (R-4.5.0\bin\x64\Rscript.exe runApp_NPPS.R 1> log_NPPS.log 2>&1)