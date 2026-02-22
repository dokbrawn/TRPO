@echo off
chcp 65001 >nul
echo === Сборка проекта CBLAS Level 3 Tests ===
echo.

set OPENBLAS_PATH=C:/OpenBLAS

if not exist "build" mkdir build
cd build

echo [1/3] Генерация проекта...
cmake .. -G "MinGW Makefiles" -DOPENBLAS_ROOT="%OPENBLAS_PATH%"

if %ERRORLEVEL% neq 0 (
    echo Ошибка генерации CMake!
    exit /b 1
)

echo.
echo [2/3] Сборка Release...
cmake --build . --config Release --parallel

if %ERRORLEVEL% neq 0 (
    echo Ошибка сборки!
    exit /b 1
)

echo.
echo [3/3] Готово!
cd ..
echo Запуск: build\Release\cblas_tests.exe