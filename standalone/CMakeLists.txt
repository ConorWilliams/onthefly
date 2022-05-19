cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(OnTheFlyStandalone LANGUAGES CXX)

# --- Import tools ----

include(../cmake/tools.cmake)

# ---- Dependencies ----

include(../cmake/CPM.cmake)

CPMAddPackage(NAME OnTheFly SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/..)

# ---- Create executables ----

add_executable(${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17 OUTPUT_NAME "OnTheFly")
target_link_libraries(${PROJECT_NAME} OnTheFly::OnTheFly)

add_executable(lattice ${CMAKE_CURRENT_SOURCE_DIR}/src/lattice.cpp)
set_target_properties(lattice PROPERTIES CXX_STANDARD 17 OUTPUT_NAME "lattice")
target_link_libraries(lattice OnTheFly::OnTheFly)
