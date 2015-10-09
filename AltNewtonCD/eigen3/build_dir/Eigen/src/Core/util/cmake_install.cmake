# Install script for directory: /usr/local/include/eigen3/Eigen/src/Core/util

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "Release")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/include/eigen3/Eigen/src/Core/util/BlasUtil.h;/usr/local/include/eigen3/Eigen/src/Core/util/Constants.h;/usr/local/include/eigen3/Eigen/src/Core/util/DisableStupidWarnings.h;/usr/local/include/eigen3/Eigen/src/Core/util/ForwardDeclarations.h;/usr/local/include/eigen3/Eigen/src/Core/util/Macros.h;/usr/local/include/eigen3/Eigen/src/Core/util/Memory.h;/usr/local/include/eigen3/Eigen/src/Core/util/Meta.h;/usr/local/include/eigen3/Eigen/src/Core/util/MKL_support.h;/usr/local/include/eigen3/Eigen/src/Core/util/NonMPL2.h;/usr/local/include/eigen3/Eigen/src/Core/util/ReenableStupidWarnings.h;/usr/local/include/eigen3/Eigen/src/Core/util/StaticAssert.h;/usr/local/include/eigen3/Eigen/src/Core/util/XprHelper.h")
  IF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
  IF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
FILE(INSTALL DESTINATION "/usr/local/include/eigen3/Eigen/src/Core/util" TYPE FILE FILES
    "/usr/local/include/eigen3/Eigen/src/Core/util/BlasUtil.h"
    "/usr/local/include/eigen3/Eigen/src/Core/util/Constants.h"
    "/usr/local/include/eigen3/Eigen/src/Core/util/DisableStupidWarnings.h"
    "/usr/local/include/eigen3/Eigen/src/Core/util/ForwardDeclarations.h"
    "/usr/local/include/eigen3/Eigen/src/Core/util/Macros.h"
    "/usr/local/include/eigen3/Eigen/src/Core/util/Memory.h"
    "/usr/local/include/eigen3/Eigen/src/Core/util/Meta.h"
    "/usr/local/include/eigen3/Eigen/src/Core/util/MKL_support.h"
    "/usr/local/include/eigen3/Eigen/src/Core/util/NonMPL2.h"
    "/usr/local/include/eigen3/Eigen/src/Core/util/ReenableStupidWarnings.h"
    "/usr/local/include/eigen3/Eigen/src/Core/util/StaticAssert.h"
    "/usr/local/include/eigen3/Eigen/src/Core/util/XprHelper.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")

