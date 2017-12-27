project(unittest)

include_directories(../src)
add_subdirectory(lib/googletest-release-1.8.0)

add_executable(runTest test.cpp)
target_link_libraries(runTest gtest gtest_main)
