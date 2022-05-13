cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

# ---- Project ----

project(
    olkmc
    VERSION 0.3.0
    DESCRIPTION "An off-lattice kinetic Monte Carlo simulation framework"
    LANGUAGES CXX
)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ---- Include guards ----

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
    message(
        FATAL_ERROR
            "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there."
    )
endif()

# System packages
find_package(Threads)

# ---- Add dependencies via CPM ----

# see https://github.com/TheLartians/CPM.cmake for more info

include(cmake/CPM.cmake)
include(cmake/tools.cmake)

# Fetch all dependcies using CPM

CPMAddPackage("gh:ConorWilliams/Threadpool#v2.1.1")

CPMAddPackage("gh:p-ranav/indicators#v2.0")

CPMAddPackage(
    NAME cereal
    GIT_TAG 48fda3f0a6b58f9a9278ec15c15e7b7c3dceeceb
    GITHUB_REPOSITORY USCiLab/cereal
    OPTIONS "SKIP_PORTABILITY_TEST ON" "JUST_INSTALL_CEREAL ON"
)

CPMAddPackage(
    NAME structopt
    GITHUB_REPOSITORY p-ranav/structopt
    # Version 0.1.1 commit
    GIT_TAG 78299d63f2d0c316b7509d32530fda8026637712
)

CPMAddPackage(
    NAME Eigen
    GITLAB_REPOSITORY libeigen/eigen
    # Need version 3.4+, not yet released
    GIT_TAG 0784d9f87b6602160e8f3e3d507d69999b92d539
    # Eigen's CMakelists are not intended for library use
    DOWNLOAD_ONLY YES
)
if(Eigen_ADDED)
    add_library(Eigen INTERFACE IMPORTED)
    target_include_directories(Eigen INTERFACE ${Eigen_SOURCE_DIR})
endif()

CPMAddPackage(
    NAME Format.cmake
    VERSION 1.6
    GITHUB_REPOSITORY TheLartians/Format.cmake
    OPTIONS # enable cmake formatting (optional)
            "FORMAT_CHECK_CMAKE ON"
)

CPMAddPackage(
    NAME tomlplusplus
    GITHUB_REPOSITORY marzer/tomlplusplus
    VERSION 2.3.0
)
if(tomlplusplus_ADDED)
    add_library(tomlplusplus INTERFACE IMPORTED)
    target_include_directories(tomlplusplus INTERFACE ${tomlplusplus_SOURCE_DIR}/include)
    target_compile_features(tomlplusplus INTERFACE cxx_std_17)
endif()

CPMAddPackage(
    NAME pcg-cpp
    GITHUB_REPOSITORY imneme/pcg-cpp
    VERSION 0.98.1
)
if(pcg-cpp_ADDED)
    add_library(pcg-cpp INTERFACE IMPORTED)
    target_include_directories(pcg-cpp INTERFACE ${pcg-cpp_SOURCE_DIR}/include)
    target_compile_features(pcg-cpp INTERFACE cxx_std_17)
endif()

# ---- Create executable ----

set(sources
    "src/potentials/spline.cpp"
    "src/potentials/EAM/data.cpp"
    "src/potentials/EAM/potential.cpp"
    "src/potentials/ADP/data.cpp"
    "src/potentials/ADP/potential.cpp"
    "src/potentials/potential_base.cpp"
    "src/sp_search/dimer/dimer.cpp"
    "src/sp_search/dimer/shrinking.cpp"
    "src/sp_search/dimer/l_shrink.cpp"
    "src/sp_search/sp_search_base.cpp"
    "src/minimise/LBFGS/lbfgs.cpp"
    "src/minimise/BB/bb.cpp"
    "src/minimise/HYBRID/hybrid.cpp"
    "src/minimise/minimiser_base.cpp"
    "src/local/classify.cpp"
    "src/local/environment.cpp"
    "src/local/geometry.cpp"
    "src/local/catalogue.cpp"
    "src/supercell.cpp"
    "src/utility.cpp"
    "src/main.cpp"
    "src/sp_search/vineyard.cpp"
    "src/sp_search/find_mech.cpp"
    "src/package/package.cpp"
    "src/kinetics/basin.cpp"
    "src/kinetics/superbasin.cpp"
    "src/kinetics/supercache.cpp"
    "src/discrete/discrete_supercell.cpp"
    "src/discrete/discrete_classify.cpp"
)

add_executable(olkmc ${sources})

target_compile_features(olkmc PUBLIC cxx_std_17)

target_compile_options(olkmc PRIVATE -Wall -Wextra -Wpedantic)

target_link_libraries(
    olkmc
    PUBLIC ${CMAKE_THREAD_LIBS_INIT}
           Eigen
           cereal
           structopt
           tomlplusplus
           pcg-cpp
           RiftenThiefpool::RiftenThiefpool
           indicators::indicators
)

target_include_directories(
    olkmc PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
                 $<INSTALL_INTERFACE:src/${PROJECT_NAME}-${PROJECT_VERSION}>
)

# ///
