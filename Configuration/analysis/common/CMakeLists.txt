PROJECT( common )

#################
# modules/macros
#################
SET( CMAKE_MODULE_PATH ${${PROJECT_NAME}_SOURCE_DIR}/../cmake ${CMAKE_MODULE_PATH} CACHE PATH "Path to CMake's Modules")


################
# cmake options
################
CMAKE_MINIMUM_REQUIRED( VERSION 2.4.6 )
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)

#################
# project version
#################


SET( ${PROJECT_NAME}_MAJOR_VERSION 1 )
SET( ${PROJECT_NAME}_MINOR_VERSION 0 )
SET( ${PROJECT_NAME}_PATCH_LEVEL   0 )

# library Unix style versioning
SET( ${PROJECT_NAME}_SOVERSION
  "${${PROJECT_NAME}_MAJOR_VERSION}" )
SET( ${PROJECT_NAME}_VERSION
  "${${PROJECT_NAME}_MAJOR_VERSION}.${${PROJECT_NAME}_MINOR_VERSION}.${${PROJECT_NAME}_PATCH_LEVEL}" )

############
# options
############

# build shared library by default
OPTION( BUILD_SHARED_LIBS "Create shared libraries" ON )

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)


###############################
# find packages (dependencies)
###############################

FIND_PACKAGE( ROOT REQUIRED )
# Next line needed to be replaced by the second line, because changed version of FindROOT.cmake does not work with first line
#SET( ROOT_EXTRA_LIBS -lGenVector -lHistPainter)
#SET( ROOT_LIBRARIES ${ROOT_LIBRARIES} -lGenVector -lHistPainter )

FIND_PACKAGE( ZTopUtils REQUIRED )
# now you can use ${ZTopUtils_INCLUDE_DIRS} and ${ZTopUtils_LIBRARIES}

#######################
# installation options
#######################

# change default install prefix
IF( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )
  SET( CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install" CACHE PATH "Install prefix" FORCE )
ENDIF( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )


# default destination for header files: ${CMAKE_INSTALL_PREFIX}/install
SET( INCLUDE_INSTALL_DIR "include" CACHE PATH "Directory to install the header files" )
MARK_AS_ADVANCED( INCLUDE_INSTALL_DIR )

# default destination for library files: ${CMAKE_INSTALL_PREFIX}/install
SET( LIB_INSTALL_DIR "lib" CACHE PATH "Directory to install the libraries" )
MARK_AS_ADVANCED( LIB_INSTALL_DIR )

# default destination for .cmake files: ${CMAKE_INSTALL_PREFIX}/install
SET( CONFIG_INSTALL_DIR "cmake" CACHE PATH "Directory to install the XXXConfig.cmake files" )
MARK_AS_ADVANCED( CONFIG_INSTALL_DIR )


# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# the RPATH to be used when installing
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# provide nicer directory layout in build directory
SET( EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin" CACHE PATH
  "EXECUTABLE_OUTPUT_PATH" FORCE )
SET( LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib" CACHE PATH
  "LIBRARY_OUTPUT_PATH" FORCE )
MARK_AS_ADVANCED( EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH )



############
# subdirs
############

ADD_SUBDIRECTORY( src )

##########
# add doxygen documentation
##########
INCLUDE (DOCUMENTATION)


##########
# install the scripts into the bin directory
##########
#FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/scripts)
#FILE(GLOB scripts RELATIVE  ${PROJECT_SOURCE_DIR}/scripts scripts/*)
#MESSAGE(STATUS "scripts found: ${scripts}")
#FOREACH( script ${scripts})
#  MESSAGE(STATUS "configuring ${script} to  ${CMAKE_BINARY_DIR}/scripts")
#  CONFIGURE_FILE( ${PROJECT_SOURCE_DIR}/scripts/${script} ${CMAKE_BINARY_DIR}/scripts/${script} @ONLY )
#  INSTALL(PROGRAMS  ${CMAKE_BINARY_DIR}/scripts/${script}
#    DESTINATION bin )
#ENDFOREACH( script ${scripts} )


########################
# print status report
########################
#Todo: write fancy macro to generate this

MESSAGE( STATUS )
MESSAGE( STATUS "-------------------------------------------------------------------------------" )
MESSAGE( STATUS "BUILD_SHARED_LIBS                       = ${BUILD_SHARED_LIBS}" )
MESSAGE( STATUS "CMAKE_INSTALL_PREFIX                    = ${CMAKE_INSTALL_PREFIX}" )
MESSAGE( STATUS "CMAKE_BUILD_TYPE                        = ${CMAKE_BUILD_TYPE}" )
MESSAGE( STATUS "BUILD_WITH_DEBUG                        = ${BUILD_WITH_DEBUG} ")


MESSAGE( STATUS "Change a value with: cmake -D<Variable>=<Value>" )
MESSAGE( STATUS "-------------------------------------------------------------------------------" )
MESSAGE( STATUS )

