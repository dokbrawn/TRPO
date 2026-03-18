cmake_minimum_required(VERSION 3.20)
project(trpo_syrk C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS ON)

add_library(syrk_lib
    src/syrk.c
)

target_include_directories(syrk_lib PUBLIC include)

target_compile_options(syrk_lib PRIVATE -O3 -Wall -Wextra -Wpedantic)

add_executable(test_syrk
    src/test_syrk.c
)

target_link_libraries(test_syrk PRIVATE syrk_lib m)

target_compile_options(test_syrk PRIVATE -O2 -Wall -Wextra -Wpedantic)

add_executable(benchmark_syrk
    src/benchmark_syrk.c
)

target_link_libraries(benchmark_syrk PRIVATE syrk_lib dl m)

target_compile_options(benchmark_syrk PRIVATE -O3 -Wall -Wextra -Wpedantic)

include(CTest)
if(BUILD_TESTING)
    add_test(NAME syrk_functional COMMAND test_syrk)
endif()
