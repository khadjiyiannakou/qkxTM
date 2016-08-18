QKXTM = ../lib/libqkxTM.a

INC= -I../include -I/fhgfs/prometheus/home/khadjiyiannakou/CudaArnoldi/Cuda-Arnoldi-master/lapack/lapack-3.5.0/lapacke/include
CC = gcc
CFLAGS = -Wall -O3 -std=c99 $(INC)
CXX = g++
CXXFLAGS = -Wall -O3 $(INC) 

LIBSFLAGS = /fhgfs/prometheus/home/khadjiyiannakou/qudaQKXTM_package/lib/external_libs/lime-1.3.2/lib/liblime.a -L/fhgfs/prometheus/home/khadjiyiannakou/CudaArnoldi/Cuda-Arnoldi-master/lapack/lapack-3.5.0
LIBS = -llapack -ltmglib -lrefblas -lgfortran
EXECUTABLES = test test_smear test_applyDslash test_invert test_checkInvert test_plaq test_invert_end_trick test_QR IRAM eigTMO eigTMO_lapack eigTMO_chiral_lapack eigTMO_chiral_phases_lapack eigTMO_lapack_check

all: $(EXECUTABLES)

eigTMO_chiral_phases_lapack: eigTMO_chiral_phases_lapack.o $(QKXTM)
	$(CXX) $(CXXFLAGS) $^ -o $@  $(LIBSFLAGS) $(LIBS)

eigTMO_chiral_lapack: eigTMO_chiral_lapack.o $(QKXTM)
	$(CXX) $(CXXFLAGS) $^ -o $@  $(LIBSFLAGS) $(LIBS)

eigTMO_lapack: eigTMO_lapack.o $(QKXTM)
	$(CXX) $(CXXFLAGS) $^ -o $@  $(LIBSFLAGS) $(LIBS)

eigTMO_lapack_check: eigTMO_lapack_check.o $(QKXTM)
	$(CXX) $(CXXFLAGS) $^ -o $@  $(LIBSFLAGS) $(LIBS)

eigTMO: eigTMO.o $(QKXTM)
	$(CXX) $(CXXFLAGS) $^ -o $@  $(LIBSFLAGS) $(LIBS)

IRAM: IRAM.o $(QKXTM)
	$(CXX) $(CXXFLAGS) $^ -o $@  $(LIBSFLAGS) $(LIBS)

test: test.o $(QKXTM)
	$(CXX) $(CXXFLAGS) $^ -o $@  $(LIBSFLAGS) $(LIBS)

test_QR: test_QR.o $(QKXTM)
	$(CXX) $(CXXFLAGS) $^ -o $@  $(LIBSFLAGS) $(LIBS)

test: test.o $(QKXTM)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBSFLAGS) $(LIBS)

test_gauge: test_gauge.o $(QKXTM)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBSFLAGS) $(LIBS)

test_plaq: test_plaq.o $(QKXTM)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBSFLAGS) $(LIBS)

test_smear: test_smear.o $(QKXTM)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBSFLAGS) $(LIBS)

test_applyDslash: test_applyDslash.o $(QKXTM)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBSFLAGS) $(LIBS)

test_invert: test_invert.o $(QKXTM)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBSFLAGS) $(LIBS)

test_checkInvert: test_checkInvert.o $(QKXTM)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBSFLAGS) $(LIBS)

test_invert_end_trick: test_invert_end_trick.o $(QKXTM)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBSFLAGS) $(LIBS)

clean:
	-rm -f $(EXECUTABLES) *.o