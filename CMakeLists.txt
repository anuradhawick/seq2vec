cmake_minimum_required(VERSION 3.0.0)
project(seq2vec VERSION 1.0.0)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp -O3 ")

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror ${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS}")
set(CMAKE_CXX_STANDARD 20)
set(Boost_NO_BOOST_CMAKE TRUE)
set(Boost_NO_SYSTEM_PATHS TRUE)

include(CTest)
enable_testing()

list(APPEND EXTRA_INCLUDE_DIRS "")
list(APPEND EXTRA_LINKS "")

find_package( Boost COMPONENTS program_options iostreams REQUIRED)

add_executable(${PROJECT_NAME} main.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

include_directories(${Boost_INCLUDE_DIRS} ./include)

target_include_directories(${PROJECT_NAME} 
    PUBLIC 
    ${EXTRA_INCLUDE_DIRS}
)

target_link_libraries(${PROJECT_NAME} 
    ${EXTRA_LINKS} 
    ${Boost_LIBRARIES}
    pthread
    z
)
