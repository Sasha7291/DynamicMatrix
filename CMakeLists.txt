cmake_minimum_required(VERSION 3.16)

set(MATRIX_LIB_NAME dynamic_matrix)

add_library(${MATRIX_LIB_NAME} INTERFACE)
target_include_directories(${MATRIX_LIB_NAME} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
