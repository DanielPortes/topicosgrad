cmake_minimum_required(VERSION 3.16)
project(MM C)

set(CMAKE_C_STANDARD 11)

include_directories(.)


add_executable(MM
        example_read.c
        mmio.c
        mmio.h)

