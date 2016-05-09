if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -qopenmp -std=c++11 -ipo -prof-use")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -fopenmp -fcilkplus")
endif()

add_executable(readmapper readmapper.cpp
                          fastareader.cpp
                          sw.cpp)