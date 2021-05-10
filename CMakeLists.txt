cmake_minimum_required(VERSION 3.0.0)
project(seq2vec VERSION 0.1.0)

include(CTest)
enable_testing()

list(APPEND EXTRA_LIB_DIRS "")
list(APPEND EXTRA_INCLUDE_DIRS "")
list(APPEND EXTRA_LINKS "")

find_package( Boost COMPONENTS program_options iostreams REQUIRED )

add_executable(${PROJECT_NAME} main.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

include_directories(${Boost_INCLUDE_DIRS} ./include)

target_include_directories(${PROJECT_NAME} 
    PUBLIC 
	${EXTRA_INCLUDE_DIRS}
)

target_link_directories(${PROJECT_NAME} 
	PUBLIC
	${EXTRA_LIB_DIRS}
)
target_link_libraries(${PROJECT_NAME} 
    ${EXTRA_LINKS} 
    ${Boost_LIBRARIES}
    pthread
)