project(unittest)

add_subdirectory(lib/googletest-release-1.8.0)
include_directories(lib/googletest-release-1.8.0/googletest/include)
include_directories(lib/googletest-release-1.8.0/googlemock/include)

add_executable(runTest test.cpp)
target_link_libraries(runTest gtest gtest_main)
