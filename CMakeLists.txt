cmake_minimum_required(VERSION 3.1)
project(osptBWT)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) ## on g++ this ensures: -std=c++11 and not -std=gnu++11

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if(NOT MYMOD_DIR)
  set(MYMOD_SOURCE_DIR ${CMAKE_SOURCE_DIR}/OnlineRlbwt/module)
  set(MYMOD_BINARY_DIR ${CMAKE_BINARY_DIR}/OnlineRlbwt/module)
  add_subdirectory(${MYMOD_SOURCE_DIR}/Basics ${MYMOD_BINARY_DIR}/Basics)
  add_subdirectory(${MYMOD_SOURCE_DIR}/BTree ${MYMOD_BINARY_DIR}/BTree)
  add_subdirectory(OnlineRlbwt)
  include_directories(OnlineRlbwt)
endif()

add_executable(osptBWT main.cpp)

# target_link_libraries(osptBWT OnlineRlbwt)
