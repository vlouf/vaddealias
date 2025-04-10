cmake_minimum_required(VERSION 3.17)

# default to a release build (must be before project())
if (DEFINED CMAKE_BUILD_TYPE)
  set(DEF_BUILD_TYPE ${CMAKE_BUILD_TYPE})
else()
  set(DEF_BUILD_TYPE Release)
endif()
set(CMAKE_BUILD_TYPE ${DEF_BUILD_TYPE} CACHE STRING "Build type: <none>|Debug|Release|RelWithDebInfo|MinSizeRel")

# initialize our project
project(vad-dealias)

# setup standard install locations
include(GNUInstallDirs)

# external dependencies
find_package(bom-util REQUIRED)
list(APPEND DEPENDENCY_LIBRARIES bom-util::bom-util)

find_package(bom-core REQUIRED)
list(APPEND DEPENDENCY_LIBRARIES bom-core::bom-core)

# fails link on el8 without this
find_package(Threads)
list(APPEND DEPENDENCY_LIBRARIES Threads::Threads)

# setup our compilation environment
setup_cplusplus()

# build our executables
add_executable(vad-dealias src/main.cc src/array_operations.cc src/cappi.cc src/corrections.cc src/metadata.cc src/io.cc)
target_link_libraries(vad-dealias ${DEPENDENCY_LIBRARIES} stdc++fs)
install(TARGETS vad-dealias DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT runtime)
