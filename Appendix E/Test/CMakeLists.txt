#----------------------------------------------------------------------------
# Setup the project
set(name MuonE)
cmake_minimum_required(VERSION 3.16...3.21)
project(${name})

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(${name} ${name}.cc ${sources} ${headers})
target_link_libraries(${name} ${Geant4_LIBRARIES} )

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build MscHanson. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS ${name} DESTINATION bin)

#
# add to testing
#
#GEANT4_ADD_TEST(em-testing-suite-${name} BUILD ${name} 
#  COMMAND ${name} ${CMAKE_CURRENT_SOURCE_DIR}/${name}.in
#  DEPENDS ${name} ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
#
#GEANT4_ADD_TEST(${name} BUILD ${name} 
#  COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/run.py
#  DEPENDS ${name} ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
