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


add_executable(env_analysis ${CMAKE_CURRENT_SOURCE_DIR}/src/env_analysis.cpp)
set_target_properties(env_analysis PROPERTIES CXX_STANDARD 17 OUTPUT_NAME "env_analysis")
target_link_libraries(env_analysis OnTheFly::OnTheFly)

add_executable(env_analysis2 ${CMAKE_CURRENT_SOURCE_DIR}/src/env_analysis2.cpp)
set_target_properties(env_analysis2 PROPERTIES CXX_STANDARD 17 OUTPUT_NAME "env_analysis2")
target_link_libraries(env_analysis2 OnTheFly::OnTheFly)

add_executable(sps_test ${CMAKE_CURRENT_SOURCE_DIR}/src/sps_test.cpp)
set_target_properties(sps_test PROPERTIES CXX_STANDARD 17 OUTPUT_NAME "sps_test")
target_link_libraries(sps_test OnTheFly::OnTheFly)



