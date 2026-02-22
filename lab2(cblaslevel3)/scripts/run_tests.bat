@echo off
chcp 65001 >nul
echo === Запуск тестов CBLAS Level 3 ===
echo.

if not exist "build\Release\cblas_tests.exe" (
    echo Ошибка: скомпилируйте проект сначала!
    echo Запустите: scripts\build.bat
    exit /b 1
)

build\Release\cblas_tests.exe
exit /b %ERRORLEVEL%