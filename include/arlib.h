#ifndef ARLIB_H
#define ARLIB_H

#include <qkxTM.h>
#include <lattice_util.h>

namespace qkxTM{

  class Arnoldi;

  void initArnold(Arnoldi &A,Arnoldi &V, Arnoldi &H);
  void initArnold(Arnoldi &V, Arnoldi &H, GaugeField &gauge);
  void arnold(int kstart,int kstop,Arnoldi &A,Arnoldi &V, Arnoldi &H);
  void arnold(int kstart,int kstop,Arnoldi &V, Arnoldi &H,GaugeField &gauge);
  void QR_hs_v2(qkxTMComplex *Q, qkxTMComplex *R, qkxTMComplex *x);
  void qrShiftsRotations(Arnoldi &V, Arnoldi &H, Arnoldi &W);
  void qrShiftsRotations_v2(Arnoldi &V, Arnoldi &H, Arnoldi &W);
  void sortEigenValues(Arnoldi &Q, Arnoldi &W);


  class Arnoldi : public LatticeField{
  private:
    int dimensionMatrix; // dimension of the matrix we want to work on
    int arnoldiTotalVectors; // total number of arnoldi vectors
    int arnoldiWantedVectors; // number of wanted arnoldi vectors
    int arnoldiUnwantedVectors;
    int maxArnoldiIter; // maximum number of arnoldi iterations
    bool doPolyAccel; // boolean to choose if want polynomial acceleration
    double tolerance; // tolerance of the arnoldi solver
    short unsigned int seed;
    LatticeInfo latInfo;
    ArnoldiInfo arInfo;

    qkxTMComplex *p_matrixNxN;
    qkxTMComplex *p_matrixNxK;
    qkxTMComplex *p_vectorN;
    qkxTMComplex *p_matrixKxK;
  public:
    qkxTMComplex *P_matrixNxN() { return p_matrixNxN;}
    qkxTMComplex *P_matrixNxK() { return p_matrixNxK;}
    qkxTMComplex *P_matrixKxK() { return p_matrixKxK;}
    qkxTMComplex *P_vectorN() { return p_vectorN;}

    int DimensionMatrix() { return dimensionMatrix;}
    int ArnoldiTotalVectors() { return arnoldiTotalVectors;}
    int ArnoldiWantedVectors() { return arnoldiWantedVectors;}
    int ArnoldiUnwantedVectors() { return arnoldiUnwantedVectors;}
    int MaxArnoldiIter() { return maxArnoldiIter;}
    bool DoPolyAccel() { return doPolyAccel;}
    double Tolerance() { return tolerance;}
    
    LatticeInfo LatInfo() { return latInfo;}
    ArnoldiInfo ArInfo() { return arInfo;}

    Arnoldi(ArnoldiInfo ,LatticeInfo , enum MATRIXTYPE matrixtype);
    ~Arnoldi();
    void QR_gs(Arnoldi &Q, Arnoldi &R);
    void QR_gs_m(Arnoldi &Q,Arnoldi &R);
    void QR_hs(Arnoldi &H, Arnoldi &Tmp);
    void QR_gvns(Arnoldi &Q, Arnoldi &R);
    void QRLapack(Arnoldi &Q, Arnoldi &R);
    void eig(Arnoldi &Q, Arnoldi &W);
    void eigLapack(Arnoldi &Q, Arnoldi &W);
    void Iram(Arnoldi &X, Arnoldi &H , Arnoldi &W);
    void Iram(Arnoldi &X, Arnoldi &W, GaugeField &gauge);
  };


}  


#endif
