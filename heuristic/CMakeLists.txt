cmake_minimum_required(VERSION 3.10)

# set the project name
project(heuristic)
SET(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -O3 -D NDEBUG -D NUSE_INIT_APPROX")

# add the executable
add_executable(heuristic main.c graph.c WELL512a.c sa_lib.c sa_state.c simulatedAnnealing.c
    reorganize.c splayTree.c queue.c darray.c)
target_include_directories(heuristic PUBLIC "${PROJECT_BINARY_DIR}")

target_link_libraries(heuristic m)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED True)
