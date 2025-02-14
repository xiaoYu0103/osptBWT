cmake_minimum_required(VERSION 3.1)
project(osptBWT)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# 强制覆盖MYMOD_DIR指向项目B的内部模块路径
if (NOT MYNOD_DIR)
    set(MYMOD_DIR ${CMAKE_SOURCE_DIR}/OnlineRlbwt/module)
    set(MYMOD_SOURCE_DIR ${CMAKE_SOURCE_DIR}/OnlineRlbwt/module)
    set(MYMOD_BINARY_DIR ${CMAKE_BINARY_DIR}/OnlineRlbwt/module)
    add_subdirectory(${MYMOD_SOURCE_DIR}/Basics ${MYMOD_BINARY_DIR}/Basics)
    add_subdirectory(${MYMOD_SOURCE_DIR}/BTree ${MYMOD_BINARY_DIR}/BTree)
endif()
include_directories(${MYMOD_SOURCE_DIR}/Basics)
include_directories(${MYMOD_SOURCE_DIR}/BTree)

# 包含项目B的目录
add_subdirectory(OnlineRlbwt)
include_directories(${CMAKE_SOURCE_DIR}/OnlineRlbwt)

# 链接依赖库（假设项目B生成库名为OnlineRlbwt）
add_executable(osptBWT main.cpp)
# target_link_libraries(osptBWT PRIVATE OnlineRlbwt)