project(MCF-SSP)

cmake_minimum_required(VERSION 2.8.12)

# C++14
add_compile_options(-std=c++14)

# write locations of test instances to file
set(GTE_BAD_20 "${CMAKE_CURRENT_SOURCE_DIR}/instances/gte/gte_bad.20")
set(GTE_BAD_40 "${CMAKE_CURRENT_SOURCE_DIR}/instances/gte/gte_bad.40")
set(GTE_BAD_60 "${CMAKE_CURRENT_SOURCE_DIR}/instances/gte/gte_bad.60")
set(GTE_BAD_160 "${CMAKE_CURRENT_SOURCE_DIR}/instances/gte/gte_bad.160")
set(GTE_BAD_200 "${CMAKE_CURRENT_SOURCE_DIR}/instances/gte/gte_bad.200")
set(GTE_BAD_460 "${CMAKE_CURRENT_SOURCE_DIR}/instances/gte/gte_bad.460")
set(GTE_BAD_510 "${CMAKE_CURRENT_SOURCE_DIR}/instances/gte/gte_bad.510")
set(GTE_BAD_1160 "${CMAKE_CURRENT_SOURCE_DIR}/instances/gte/gte_bad.1160")
set(GTE_BAD_15100 "${CMAKE_CURRENT_SOURCE_DIR}/instances/gte/gte_bad.15100")
set(GTE_BAD_15710 "${CMAKE_CURRENT_SOURCE_DIR}/instances/gte/gte_bad.15710")
set(GTE_BAD_176280 "${CMAKE_CURRENT_SOURCE_DIR}/instances/gte/gte_bad.176280")
set(GTE_BAD_298300 "${CMAKE_CURRENT_SOURCE_DIR}/instances/gte/gte_bad.298300")
set(GTE_BAD_35620 "${CMAKE_CURRENT_SOURCE_DIR}/instances/gte/gte_bad.35620")
set(GTE_BAD_451760 "${CMAKE_CURRENT_SOURCE_DIR}/instances/gte/gte_bad.451760")
set(GTE_BAD_49320 "${CMAKE_CURRENT_SOURCE_DIR}/instances/gte/gte_bad.49320")
set(GTE_BAD_508829 "${CMAKE_CURRENT_SOURCE_DIR}/instances/gte/gte_bad.508829")
set(GTE_BAD_60090 "${CMAKE_CURRENT_SOURCE_DIR}/instances/gte/gte_bad.60090")
set(GTE_BAD_6410 "${CMAKE_CURRENT_SOURCE_DIR}/instances/gte/gte_bad.6410")
set(GTE_BAD_65330 "${CMAKE_CURRENT_SOURCE_DIR}/instances/gte/gte_bad.65330")
set(GTE_BAD_6830 "${CMAKE_CURRENT_SOURCE_DIR}/instances/gte/gte_bad.6830")

configure_file(test/instances.h.in test/instances.h)

add_subdirectory(test)
