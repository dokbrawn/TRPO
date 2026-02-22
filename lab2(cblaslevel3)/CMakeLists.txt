cmake_minimum_required(VERSION 3.15)
project(cblas_level3_tests C)

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


option(USE_MOCK_BLAS "Use mock library with wrong interface" OFF)
option(USE_OPENBLAS_THREADS "Enable OpenBLAS threading" ON)


set(TEST_SOURCES
    ${CMAKE_SOURCE_DIR}/src/main.c
    ${CMAKE_SOURCE_DIR}/src/test_gemm.c
    ${CMAKE_SOURCE_DIR}/src/test_symm_hemm.c
    ${CMAKE_SOURCE_DIR}/src/test_trmm_trsm.c
    ${CMAKE_SOURCE_DIR}/src/test_syrk_herk.c
    ${CMAKE_SOURCE_DIR}/src/test_syr2k_her2k.c
)


add_executable(cblas_tests WIN32 ${TEST_SOURCES})


if(WIN32 AND CMAKE_C_COMPILER_ID MATCHES "GNU")
    set_target_properties(cblas_tests PROPERTIES 
        LINK_FLAGS "-mconsole -Wl,--subsystem,console")
endif()


include_directories(${CMAKE_SOURCE_DIR}/include)


if(USE_OPENBLAS_THREADS)
    target_compile_definitions(cblas_tests PRIVATE USE_THREADS)
endif()


if(USE_MOCK_BLAS)
    target_sources(cblas_tests PRIVATE ${CMAKE_SOURCE_DIR}/lib/mock_openblas.c)
    target_compile_definitions(cblas_tests PRIVATE USE_MOCK_BLAS)
    message(STATUS "Using MOCK BLAS library")
else()
    target_link_libraries(cblas_tests ${OPENBLAS_LIB})
    message(STATUS "Using OpenBLAS library")
endif()


if(WIN32 AND NOT USE_MOCK_BLAS)
    add_custom_command(TARGET cblas_tests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${OPENBLAS_ROOT}/bin/libopenblas.dll"
        $<TARGET_FILE_DIR:cblas_tests>
        COMMENT "Copying OpenBLAS DLL")
endif()


enable_testing()
add_test(NAME cblas_level3_tests COMMAND cblas_tests)