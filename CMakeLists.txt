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

# ---- pndf pnm_aniso_he / PMP begin ----
set(PMP_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(PMP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(PMP_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(PMP_BUILD_VIEWERS OFF CACHE BOOL "" FORCE)
set(PMP_INSTALL OFF CACHE BOOL "" FORCE)
set(PMP_STRICT_COMPILATION OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/extern/pmp-library/CMakeLists.txt")
  add_subdirectory(extern/pmp-library EXCLUDE_FROM_ALL)
else()
  message(FATAL_ERROR "extern/pmp-library is missing. Run: git submodule add https://github.com/pmp-library/pmp-library extern/pmp-library")
endif()
add_executable(pndf_run_pnm_aniso_he apps/pndf_run_pnm_aniso_he.cpp)
target_link_libraries(pndf_run_pnm_aniso_he PRIVATE pndf_core pmp)
# ---- pndf pnm_aniso_he / PMP end ----
