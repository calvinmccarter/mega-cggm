# Install script for directory: /usr/local/include/eigen3/Eigen/src/SparseCore

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
   "/usr/local/include/eigen3/Eigen/src/SparseCore/AmbiVector.h;/usr/local/include/eigen3/Eigen/src/SparseCore/CompressedStorage.h;/usr/local/include/eigen3/Eigen/src/SparseCore/ConservativeSparseSparseProduct.h;/usr/local/include/eigen3/Eigen/src/SparseCore/MappedSparseMatrix.h;/usr/local/include/eigen3/Eigen/src/SparseCore/SparseBlock.h;/usr/local/include/eigen3/Eigen/src/SparseCore/SparseColEtree.h;/usr/local/include/eigen3/Eigen/src/SparseCore/SparseCwiseBinaryOp.h;/usr/local/include/eigen3/Eigen/src/SparseCore/SparseCwiseUnaryOp.h;/usr/local/include/eigen3/Eigen/src/SparseCore/SparseDenseProduct.h;/usr/local/include/eigen3/Eigen/src/SparseCore/SparseDiagonalProduct.h;/usr/local/include/eigen3/Eigen/src/SparseCore/SparseDot.h;/usr/local/include/eigen3/Eigen/src/SparseCore/SparseFuzzy.h;/usr/local/include/eigen3/Eigen/src/SparseCore/SparseMatrix.h;/usr/local/include/eigen3/Eigen/src/SparseCore/SparseMatrixBase.h;/usr/local/include/eigen3/Eigen/src/SparseCore/SparsePermutation.h;/usr/local/include/eigen3/Eigen/src/SparseCore/SparseProduct.h;/usr/local/include/eigen3/Eigen/src/SparseCore/SparseRedux.h;/usr/local/include/eigen3/Eigen/src/SparseCore/SparseSelfAdjointView.h;/usr/local/include/eigen3/Eigen/src/SparseCore/SparseSparseProductWithPruning.h;/usr/local/include/eigen3/Eigen/src/SparseCore/SparseTranspose.h;/usr/local/include/eigen3/Eigen/src/SparseCore/SparseTriangularView.h;/usr/local/include/eigen3/Eigen/src/SparseCore/SparseUtil.h;/usr/local/include/eigen3/Eigen/src/SparseCore/SparseVector.h;/usr/local/include/eigen3/Eigen/src/SparseCore/SparseView.h;/usr/local/include/eigen3/Eigen/src/SparseCore/TriangularSolver.h")
  IF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
  IF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  ENDIF (CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
FILE(INSTALL DESTINATION "/usr/local/include/eigen3/Eigen/src/SparseCore" TYPE FILE FILES
    "/usr/local/include/eigen3/Eigen/src/SparseCore/AmbiVector.h"
    "/usr/local/include/eigen3/Eigen/src/SparseCore/CompressedStorage.h"
    "/usr/local/include/eigen3/Eigen/src/SparseCore/ConservativeSparseSparseProduct.h"
    "/usr/local/include/eigen3/Eigen/src/SparseCore/MappedSparseMatrix.h"
    "/usr/local/include/eigen3/Eigen/src/SparseCore/SparseBlock.h"
    "/usr/local/include/eigen3/Eigen/src/SparseCore/SparseColEtree.h"
    "/usr/local/include/eigen3/Eigen/src/SparseCore/SparseCwiseBinaryOp.h"
    "/usr/local/include/eigen3/Eigen/src/SparseCore/SparseCwiseUnaryOp.h"
    "/usr/local/include/eigen3/Eigen/src/SparseCore/SparseDenseProduct.h"
    "/usr/local/include/eigen3/Eigen/src/SparseCore/SparseDiagonalProduct.h"
    "/usr/local/include/eigen3/Eigen/src/SparseCore/SparseDot.h"
    "/usr/local/include/eigen3/Eigen/src/SparseCore/SparseFuzzy.h"
    "/usr/local/include/eigen3/Eigen/src/SparseCore/SparseMatrix.h"
    "/usr/local/include/eigen3/Eigen/src/SparseCore/SparseMatrixBase.h"
    "/usr/local/include/eigen3/Eigen/src/SparseCore/SparsePermutation.h"
    "/usr/local/include/eigen3/Eigen/src/SparseCore/SparseProduct.h"
    "/usr/local/include/eigen3/Eigen/src/SparseCore/SparseRedux.h"
    "/usr/local/include/eigen3/Eigen/src/SparseCore/SparseSelfAdjointView.h"
    "/usr/local/include/eigen3/Eigen/src/SparseCore/SparseSparseProductWithPruning.h"
    "/usr/local/include/eigen3/Eigen/src/SparseCore/SparseTranspose.h"
    "/usr/local/include/eigen3/Eigen/src/SparseCore/SparseTriangularView.h"
    "/usr/local/include/eigen3/Eigen/src/SparseCore/SparseUtil.h"
    "/usr/local/include/eigen3/Eigen/src/SparseCore/SparseVector.h"
    "/usr/local/include/eigen3/Eigen/src/SparseCore/SparseView.h"
    "/usr/local/include/eigen3/Eigen/src/SparseCore/TriangularSolver.h"
    )
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Devel")

