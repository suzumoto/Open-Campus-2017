cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

set(PROGS classical_MonteCarlo)
project(Classical_MonteCarlo CXX)
include_directories(${PROJECT_BINARY_DIR} ${PROJECT_SOURCE_DIR})

find_package(ALPS REQUIRED PATHS ${ALPS_ROOT_DIR} $ENV{ALPS_HOME} NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "Found ALPS: ${ALPS_ROOT_DIR} (version: ${ALPS_VERSION})")
include(${ALPS_USE_FILE})

foreach(name ${PROGS})
  set(target_name ${name})
  add_executable(${target_name} ${name}.cpp)
  target_link_libraries(${name}	${ALPS_LIBRARIES})
  add_alps_test(${name})	
endforeach(name)

