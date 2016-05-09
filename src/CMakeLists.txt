if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -ipo")
endif()

add_executable(readmapper readmapper.cpp
                          backtrack.cpp
                          fastareader.cpp
                          sw.cpp)
