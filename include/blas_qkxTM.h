#include <qkxTM.h>
#include <lattice_util.h>



class SU3 {

 private:
  int ncolors;
 public:
  qkxTMComplex *M;

 SU3():ncolors(NCOLORS){
    M = (qkxTMComplex*)malloc(ncolors*ncolors*sizeof(qkxTMComplex));
    if( M == NULL){
      fprintf(stderr,"Error allocate memory for SU3 matrix");
      exit(EXIT_FAILURE);
    }

    for(int i = 0 ; i< ncolors*ncolors ; i++){
      M[i].real()=0.;
      M[i].imag()=0.;
    }
  } // inline constructor      

  ~SU3(){;}

inline void zero(){

  for(int i = 0 ; i < ncolors*ncolors ; i++){
    M[i].real()=0.;
    M[i].imag()=0.;
  }
    
}



SU3& operator*(const SU3 &B){
    
  qkxTMComplex C[NCOLORS][NCOLORS];
  for(int c1 = 0 ; c1 < NCOLORS ; c1++)
    for(int c2 = 0 ; c2 < NCOLORS ; c2++)
      for(int c = 0 ; c < NCOLORS ; c++)
	C[c1][c2] = C[c1][c2] + this->M[MU(c1,c)] * B.M[MU(c,c2)];
  
  for(int c1 = 0 ; c1 < NCOLORS ; c1++)
    for(int c2 = 0 ; c2 < NCOLORS ; c2++)
      this->M[MU(c1,c2)] = C[c1][c2];
  
  return *this;
}


SU3& operator*(const double a){

  for(int c1 = 0 ; c1 < NCOLORS ; c1++)
    for(int c2 = 0 ; c2 < NCOLORS ; c2++)
      this->M[MU(c1,c2)] = a *  this->M[MU(c1,c2)];
  return *this;  
}
  
SU3& operator+(const SU3 &B){

  qkxTMComplex C[NCOLORS][NCOLORS];
  for(int c1 = 0 ; c1 < NCOLORS ; c1++)
    for(int c2 = 0 ; c2 < NCOLORS ; c2++)
      C[c1][c2] = this->M[MU(c1,c2)] + B.M[MU(c1,c2)];

  for(int c1 = 0 ; c1 < NCOLORS ; c1++)
    for(int c2 = 0 ; c2 < NCOLORS ; c2++)
      this->M[MU(c1,c2)] = C[c1][c2];
  
  return *this;
}

SU3& operator=(const SU3 &B){

  for(int c1 = 0 ; c1 < NCOLORS ; c1++)
    for(int c2 = 0 ; c2 < NCOLORS ; c2++)
      this->M[MU(c1,c2)] = B.M[MU(c1,c2)];
  return *this;
}
  
void dagger(){
  qkxTMComplex C[ncolors][ncolors];
  for(int c1 = 0 ; c1 < NCOLORS ; c1++)
    for(int c2 = 0 ; c2 < NCOLORS ; c2++)
      C[c2][c1] = this->M[MU(c1,c2)];

  for(int c1 = 0 ; c1 < NCOLORS ; c1++)
    for(int c2 = 0 ; c2 < NCOLORS ; c2++)
      this->M[MU(c1,c2)] = conj(C[c1][c2]);
  //  return this;
}

qkxTMComplex traceColor(){
  return this->M[MU(0,0)] + this->M[MU(1,1)] + this->M[MU(2,2)];
}

qkxTMComplex det(){

  return this->M[MU(0,0)]*this->M[MU(1,1)]*this->M[MU(2,2)] + this->M[MU(0,1)]*this->M[MU(1,2)]*this->M[MU(2,0)] + this->M[MU(0,2)]*this->M[MU(1,0)]*this->M[MU(2,1)] - this->M[MU(2,0)]*this->M[MU(1,1)]*this->M[MU(0,2)] - this->M[MU(2,1)]*this->M[MU(1,2)]*this->M[MU(0,0)] - this->M[MU(2,2)]*this->M[MU(1,0)]*this->M[MU(0,1)];

}


};

void projectSU3(qkxTMComplex **gauge, LatticeInfo *latInfo);



class spinColor {


 private:
  int ncolors;
  int nspins;
 public:
  qkxTMComplex *M;

 spinColor():ncolors(NCOLORS),nspins(NSPINS){
    M = (qkxTMComplex*)malloc(ncolors*nspins*sizeof(qkxTMComplex));
    if( M == NULL){
      fprintf(stderr,"Error allocate memory for spinColor matrix");
      exit(EXIT_FAILURE);
    }

    for(int i = 0 ; i< ncolors*nspins ; i++){
      M[i].real()=0.;
      M[i].imag()=0.;
    }
  } // inline constructor      

  ~spinColor(){;}

  void zero(){

    for(int i = 0 ; i < ncolors*nspins ; i++){
      M[i].real()=0.;
      M[i].imag()=0.;
    }
    
  }

  spinColor& operator*(const SU3 &u){
    qkxTMComplex spincolor[4][3];
    for(int mu = 0 ; mu < nspins ; mu++)
      for(int ic = 0 ; ic < ncolors ; ic++){
	spincolor[mu][ic] = u.M[MU(ic,0)] * this->M[MS(mu,0)] + u.M[MU(ic,1)] * this->M[MS(mu,1)] + u.M[MU(ic,2)] * this->M[MS(mu,2)];
    }
    for(int mu = 0 ; mu < nspins ; mu++)
      for(int ic = 0 ; ic < ncolors ; ic++)
	this->M[MS(mu,ic)] = spincolor[mu][ic];

    return *this;
  }

  spinColor& operator=(const spinColor &B){

    for(int mu = 0 ; mu < nspins ; mu++)
      for(int c = 0 ; c < NCOLORS ; c++)
	this->M[MS(mu,c)] = B.M[MS(mu,c)];
    return *this;
  }

  spinColor& operator*(const double a){

    for(int mu = 0 ; mu < nspins ; mu++)
      for(int c = 0 ; c < NCOLORS ; c++)
	this->M[MS(mu,c)] = a*this->M[MS(mu,c)];
    return *this;
  }

  spinColor& operator+(const spinColor &B){

    for(int mu = 0 ; mu < nspins ; mu++)
      for(int c = 0 ; c < NCOLORS ; c++)
	this->M[MS(mu,c)] = this->M[MS(mu,c)] + B.M[MS(mu,c)];
    return *this;
  }




};


double reCdotXdagY(qkxTM::ColorSpinorField &x , qkxTM::ColorSpinorField &y);
void X_eq_X_p_aY(qkxTM::ColorSpinorField &x , double a , qkxTM::ColorSpinorField &y);
void X_eq_aX_p_Y(qkxTM::ColorSpinorField &x , double a , qkxTM::ColorSpinorField &y);

qkxTMComplex dotProduct( qkxTMComplex *a, qkxTMComplex *b, int n);
void matrixXvec(qkxTMComplex *a, qkxTMComplex *vin, qkxTMComplex *vout,int N);
void matrixXmatrix(qkxTMComplex *a, qkxTMComplex *b, qkxTMComplex *c,int N);
void matrixDagXmatrix(qkxTMComplex *a, qkxTMComplex *b, qkxTMComplex *c,int N);
void vecXmatrix(qkxTMComplex *v1, qkxTMComplex *m, qkxTMComplex *v2,int N);


