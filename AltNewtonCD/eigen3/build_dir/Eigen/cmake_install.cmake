# Install script for directory: /usr/local/include/eigen3/Eigen

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
   "/usr/local/include/eigen3/Eigen/Array;/usr/local/include/eigen3/Eigen/Cholesky;/usr/local/include/eigen3/Eigen/CholmodSupport;/usr/local/include/eigen3/Eigen/Core;/usr/local/include/eigen3/Eigen/Dense;/usr/local/include/eigen3/Eigen/Eigen;/usr/local/include/eigen3/Eigen/Eigen2Support;/usr/local/include/eigen3/Eigen/Eigenvalues;/usr/local/include/eigen3/Eigen/Geometry;/usr/local/include/eigen3/Eigen/Householder;/usr/local/include/eigen3/Eigen/IterativeLinearSolvers;/usr/local/include/eigen3/Eigen/Jacobi;/usr/local/include/eigen3/Eigen/LeastSquares;/usr/local/include/eigen3/Eigen/LU;/usr/local/include/eigen3/Eigen/MetisSupport;/usr/local/include/eigen3/Eigen/OrderingMethods;/usr/local/include/eigen3/Eigen/PardisoSupport;/usr/local/include/eigen3/Eigen/PaStiXSupport;/usr/local/include/eigen3/Eigen/QR;/usr/local/include/eigen3/Eigen/QtAlignedMalloc;/usr/local/include/eigen3/Eigen/Sparse;/usr/local/include/eigen3/Eigen/SparseCholesky;/usr/local/include/eigen3/Eigen/SparseCore;/usr/local/include/eigen3/Eigen/SparseLU;/usr/local/include/eigen3/Eigen/SparseQR;/usr/local/include/eigen3/Eigen/SPQRSupport;/usr/local/include/eigen3/Eigen/StdDeque;/usr/local/include/eigen3/Eigen/StdList;/usr/local/include/eigen3/Eigen/StdVector;/usr/local/include/eigen3/Eigen/SuperLUSupport;/usr/local/include/eigen3/Eigen/SVD;/usr/local/include/eigen3/Eigen/UmfPackSupport")
  IF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
  IF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
FILE(INSTALL DESTINATION "/usr/local/include/eigen3/Eigen" TYPE FILE FILES
    "/usr/local/include/eigen3/Eigen/Array"
    "/usr/local/include/eigen3/Eigen/Cholesky"
    "/usr/local/include/eigen3/Eigen/CholmodSupport"
    "/usr/local/include/eigen3/Eigen/Core"
    "/usr/local/include/eigen3/Eigen/Dense"
    "/usr/local/include/eigen3/Eigen/Eigen"
    "/usr/local/include/eigen3/Eigen/Eigen2Support"
    "/usr/local/include/eigen3/Eigen/Eigenvalues"
    "/usr/local/include/eigen3/Eigen/Geometry"
    "/usr/local/include/eigen3/Eigen/Householder"
    "/usr/local/include/eigen3/Eigen/IterativeLinearSolvers"
    "/usr/local/include/eigen3/Eigen/Jacobi"
    "/usr/local/include/eigen3/Eigen/LeastSquares"
    "/usr/local/include/eigen3/Eigen/LU"
    "/usr/local/include/eigen3/Eigen/MetisSupport"
    "/usr/local/include/eigen3/Eigen/OrderingMethods"
    "/usr/local/include/eigen3/Eigen/PardisoSupport"
    "/usr/local/include/eigen3/Eigen/PaStiXSupport"
    "/usr/local/include/eigen3/Eigen/QR"
    "/usr/local/include/eigen3/Eigen/QtAlignedMalloc"
    "/usr/local/include/eigen3/Eigen/Sparse"
    "/usr/local/include/eigen3/Eigen/SparseCholesky"
    "/usr/local/include/eigen3/Eigen/SparseCore"
    "/usr/local/include/eigen3/Eigen/SparseLU"
    "/usr/local/include/eigen3/Eigen/SparseQR"
    "/usr/local/include/eigen3/Eigen/SPQRSupport"
    "/usr/local/include/eigen3/Eigen/StdDeque"
    "/usr/local/include/eigen3/Eigen/StdList"
    "/usr/local/include/eigen3/Eigen/StdVector"
    "/usr/local/include/eigen3/Eigen/SuperLUSupport"
    "/usr/local/include/eigen3/Eigen/SVD"
    "/usr/local/include/eigen3/Eigen/UmfPackSupport"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/usr/local/include/eigen3/build_dir/Eigen/src/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

