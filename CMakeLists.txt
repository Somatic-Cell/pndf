cmake_minimum_required(VERSION 3.23)
project(pndf_qem_phase1 LANGUAGES CXX)

option(PNDF_BUILD_TESTS "Build tests" ON)

add_library(pndf_core
    src/normal_map.cpp
    src/torus_mesh.cpp
    src/qem.cpp
    src/edge_collapse.cpp
    src/mesh_io.cpp
    src/eval.cpp
)

target_include_directories(pndf_core PUBLIC include)
target_compile_features(pndf_core PUBLIC cxx_std_20)

if (MSVC)
    target_compile_options(pndf_core PRIVATE /W4 /permissive- /EHsc)
else()
    target_compile_options(pndf_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(pndf_run_qem apps/pndf_run_qem.cpp)
target_link_libraries(pndf_run_qem PRIVATE pndf_core)

add_executable(pndf_eval_mesh apps/pndf_eval_mesh.cpp)
target_link_libraries(pndf_eval_mesh PRIVATE pndf_core) 

add_executable(pndf_run_pnm_aniso apps/pndf_run_pnm_aniso.cpp) 
target_link_libraries(pndf_run_pnm_aniso PRIVATE pndf_core)

if (PNDF_BUILD_TESTS)
    enable_testing()
    add_executable(pndf_tests
        tests/test_main.cpp
        tests/test_torus.cpp
        tests/test_qem.cpp
        tests/test_collapse.cpp
    )
    target_link_libraries(pndf_tests PRIVATE pndf_core)
    add_test(NAME pndf_tests COMMAND pndf_tests)
endif()
