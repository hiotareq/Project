cmake_minimum_required(VERSION 3.0.0)
project(MyProject VERSION 0.1.0)

add_executable(MyProject main.cpp Particle.cpp Particle.h Src.txt PartDist.cpp PartDist.h)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)