INC= -I../include -I/home/khadjiyiannakou/qudaQKXTM_package/lib/external_libs/lime-1.3.2/include -I/home/khadjiyiannakou/CudaArnoldi/Cuda-Arnoldi-master/lapack/lapack-3.5.0/lapacke/include
CC = gcc
CFLAGS = -Wall -O3 -std=c99 $(INC)
CXX = g++
CXXFLAGS = -Wall -O3 $(INC)

LIBSFLAGS = -L/home/khadjiyiannakou/CudaArnoldi/Cuda-Arnoldi-master/lapack/lapack-3.5.0
LIBS = -llapack -lrefblas -ltmglib -lgfortran
QKXTM = libqkxTM.a
QKXTM_OBJS = lattice_util.o qkxTM_read_gauge.o blas_qkxTM.o solvers.o arlib.o
QKXTM_HDRS = lattice_util.h qkxTM.h blas_qkxTM.h solvers.h
