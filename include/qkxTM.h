#ifndef _QKXTM_H
#define _QKXTM_H

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <complex>
#include <string.h>


//#ifndef __COMPLEX__KRIKITOS
//#define __COMPLEX__KRIKITOS
//
//#endif



typedef std::complex<double> qkxTMComplex;

enum MATRIXTYPE{matrixNxN,matrixNxK,matrixKxK,vectorN};

inline qkxTMComplex operator*(const double a , const qkxTMComplex b)
{
  qkxTMComplex c;
  c.real() = a*b.real();
  c.imag() = a*b.imag();
  return c;
}

inline qkxTMComplex operator*(const qkxTMComplex b , const double a)
{
  qkxTMComplex c;
  c.real() = a*b.real();
  c.imag() = a*b.imag();
  return c;
}

inline qkxTMComplex operator+(const double a , const qkxTMComplex b)
{
  qkxTMComplex c;
  c.real() = a + b.real();
  c.imag() = b.imag();
  return c;
}




#define NDIM 4  // the timeSpace has always 4 dimensions
#define NSPINS 4
#define NCOLORS 3
#define NDF 2 // for real and imag part

#define MS(spin,color) ( (spin)*NCOLORS + color )                          // macro for spinor
#define MG(v,color1,color2) ( (v)*NCOLORS*NCOLORS + (color1)*NCOLORS + color2 )                           // macro for gauge field
#define MU(color1,color2) ( (color1)*NCOLORS + color2 )


typedef struct {
  int L[4]; // stores the lattice length in each direction
  double kappa;
  double mu;
  int twistSign; // +1 or -1 only
  double tol;        // need it when try inverter
  double maxIter;
  double reliableDelta;
  int NsmearAPE;
  int NsmearGauss;
  double alphaAPE;
  double alphaGauss;
  int boundaryT;   // +1 for periodic boundary conditions -1 for anti-periodic boundary conditions
} LatticeInfo;

typedef struct{
  int dimensionMatrix; // dimension of the matrix we want to work on                             
  int arnoldiTotalVectors; // total number of arnoldi vectors                                    
  int arnoldiWantedVectors; // number of wanted arnoldi vectors                                  
  int arnoldiUnwantedVectors;
  int maxArnoldiIter; // maximum number of arnoldi iterations                                    
  bool doPolyAccel; // boolean to choose if want polynomial acceleration                         
  double tolerance; // tolerance of the arnoldi solver    
  short unsigned int seed;
}ArnoldiInfo;
//////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////

int qkxTM_getLimeMessage(char *fname);
int qkxTM_getGaugeLime(char *fname, qkxTMComplex **u, LatticeInfo *latInfo);
int qkxTM_getGaugeAndreas(char *fname, qkxTMComplex **u, LatticeInfo *latInfo);


#endif
