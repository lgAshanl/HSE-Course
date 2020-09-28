cmake_minimum_required(VERSION 3.16)

include(GoogleTest)

project("BigInteger")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)

execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .

  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )

if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories("${gtest_SOURCE_DIR}/include")
endif()

# Now simply link against gtest or gtest_main as needed. Eg
add_executable(biginteger tests.cpp biginteger.h biginteger.cpp)
target_link_libraries(biginteger PRIVATE gtest_main)
target_compile_features(biginteger PRIVATE cxx_std_20)
target_compile_options(biginteger PUBLIC -Wall -Wextra)
add_test(NAME biginteger_test COMMAND biginteger)

option(ENABLE_SANITIZERS OFF)
if(ENABLE_SANITIZERS)
  add_library(sanitizers INTERFACE)
  target_compile_options(sanitizers INTERFACE -fsanitize=address,undefined,leak)
  target_link_options(sanitizers INTERFACE -fsanitize=address,undefined,leak)
  target_link_libraries(biginteger PRIVATE sanitizers)
endif()