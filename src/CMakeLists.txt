CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(Kmatrix)

SET(SRC alternatives.h ellipse_operations.h Hmatrix.h jac_eqs.h Kmatrix.h SVD.h ccmath_svd.cpp kmatrix_main.cpp matrix.cpp matrix.h vector.cpp LM.h)

ADD_EXECUTABLE(Kmatrix ${SRC})
