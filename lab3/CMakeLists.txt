cmake_minimum_required(VERSION 3.15)
project(SYRC_Lab3 C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

if(WIN32)
    if(DEFINED ENV{OPENBLAS_PATH})
        set(OPENBLAS_ROOT "$ENV{OPENBLAS_PATH}")
    else()
        set(OPENBLAS_ROOT "C:/OpenBLAS")
    endif()
    message(STATUS "OpenBLAS Root: ${OPENBLAS_ROOT}")
    include_directories("${OPENBLAS_ROOT}/include")
    link_directories("${OPENBLAS_ROOT}/lib")
    set(OPENBLAS_LIB "libopenblas")
endif()

find_package(OpenMP)
if(OpenMP_C_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    message(STATUS "OpenMP Found: TRUE")
else()
    message(STATUS "OpenMP Found: FALSE")
endif()

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/tests)

set(SYRC_SOURCES
    ${CMAKE_SOURCE_DIR}/src/csyrk_zsyrk.c
)

add_executable(syrc_tests 
    ${CMAKE_SOURCE_DIR}/tests/main.c
    ${CMAKE_SOURCE_DIR}/src/test_syrc_errors.c
    ${SYRC_SOURCES}
)
target_link_libraries(syrc_tests ${OPENBLAS_LIB})
if(OpenMP_C_FOUND)
    target_link_libraries(syrc_tests OpenMP::OpenMP_C)
endif()

add_executable(syrc_benchmark 
    ${CMAKE_SOURCE_DIR}/src/benchmark_syrc.c
    ${SYRC_SOURCES}
)
target_link_libraries(syrc_benchmark ${OPENBLAS_LIB})
if(OpenMP_C_FOUND)
    target_link_libraries(syrc_benchmark OpenMP::OpenMP_C)
endif()

if(WIN32)
    add_custom_command(TARGET syrc_tests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${OPENBLAS_ROOT}/bin/libopenblas.dll"
        $<TARGET_FILE_DIR:syrc_tests>)
    
    add_custom_command(TARGET syrc_benchmark POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${OPENBLAS_ROOT}/bin/libopenblas.dll"
        $<TARGET_FILE_DIR:syrc_benchmark>)
endif()

message(STATUS "Using OpenBLAS: ${OPENBLAS_LIB}")