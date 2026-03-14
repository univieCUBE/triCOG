cmake_minimum_required(VERSION 3.20)
project(new-COG LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Option to disable Twog Singletons for expensive calculations
option(NO_TWOG_SINGLETON "Disable Twog Singletons for expensive calculations" OFF)

if(NO_TWOG_SINGLETON)
    add_compile_definitions(NO_TWOG_SINGLETON)
    message(STATUS "Optimization enabled: NO_TWOG_SINGLETON is defined.")
endif()

# Find Python and pybind11
set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)
# Find spdlog (conda-installed)
find_package(spdlog CONFIG REQUIRED)

# Fetch Catch2 for testing
include(FetchContent)
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.5.2
)
FetchContent_MakeAvailable(Catch2)

# Create triangle library
add_library(triangles STATIC cogmaker.cpp)
target_include_directories(triangles PUBLIC .)
# Link spdlog
target_link_libraries(triangles PRIVATE spdlog::spdlog)

# Create extension module with bindings
pybind11_add_module(cogmaker bindings.cpp)
# Link triangle library
target_link_libraries(cogmaker PRIVATE triangles)

# Create test executable
add_executable(unit_tests ../tests/COGtriangle/C++/test_cogmaker.cpp)
target_link_libraries(unit_tests PRIVATE triangles Catch2::Catch2WithMain)

# Enable testing
enable_testing()
add_test(NAME AllTests COMMAND unit_tests)

