#include <lattice_util.h>
#include <arlib.h>
#include <ctime>
#include <blas_qkxTM.h>
//#include <lapacke.h>


extern "C" void zgemm_(char *, char *, int * , int *, int *, qkxTMComplex *,qkxTMComplex *, int *, qkxTMComplex *, int *, qkxTMComplex* , qkxTMComplex *, int*);
//extern "C" void zgeev_(char*,char*,int*,qkxTMComplex *,int*,qkxTMComplex *,qkxTMComplex *,int*,qkxTMComplex *,int*);
extern "C" void zgeev_( char* jobvl, char* jobvr, int* n, qkxTMComplex* a,
		   int* lda, qkxTMComplex* w, qkxTMComplex* vl, int* ldvl, qkxTMComplex* vr, int* ldvr,
		   qkxTMComplex* work, int* lwork, double* rwork, int* info );

#define NVEC 4

int main(){

  using namespace qkxTM;

  LatticeInfo latInfo;
  latInfo.L[0] = 4;
  latInfo.L[1] = 4;
  latInfo.L[2] = 4;
  latInfo.L[3] = 4;
  latInfo.mu = 0.01;
  latInfo.kappa = 1./(2*(4.+latInfo.mu));
  latInfo.twistSign = +1;

  latInfo.tol = 1e-8;
  latInfo.maxIter = 1000;
  latInfo.reliableDelta = 0.01;

  ArnoldiInfo arInfo;
  arInfo.dimensionMatrix = latInfo.L[0] * latInfo.L[1] * latInfo.L[2] * latInfo.L[3] * 12;
  arInfo.arnoldiTotalVectors = 11;
  arInfo.arnoldiWantedVectors = 10;
  arInfo.maxArnoldiIter = 1000000;
  arInfo.doPolyAccel = false;
  arInfo.tolerance=1e-8;
  arInfo.seed = 12345;


  GaugeField *gaugefield = new GaugeField(&latInfo);
  printf("Gauge length is %d\n",gaugefield->TotalBytes());
  gaugefield->create();
  qkxTMComplex **u = gaugefield->P_gauge();


  char filename[] = "conf.4x4x4x4";
  qkxTM_getGaugeAndreas(filename, u, &latInfo);


  double plaq = gaugefield->calculatePlaquette();
  printf("The plaquette is %8.7f\n",plaq);

  int NL = arInfo.dimensionMatrix;  
  qkxTMComplex *A = (qkxTMComplex*)malloc(NL*NL*sizeof(qkxTMComplex));
  qkxTMComplex *V = (qkxTMComplex*)malloc(NL*NL*sizeof(qkxTMComplex));
#define VOL (latInfo.L[0]*latInfo.L[1]*latInfo.L[2]*latInfo.L[3])  
#define DELTA(mu,nu) ((mu)==(nu) ? (qkxTMComplex) {1.,0.} : (qkxTMComplex) {0.,0.}) 
#define IV(x,y,z,t)  ( (x) + latInfo.L[0] * (y) + latInfo.L[0]*latInfo.L[1] * (z) + latInfo.L[0]*latInfo.L[1]*latInfo.L[2] * (t))
#define MA(mu,a,iv,nu,b,ivp)  ((mu)*3*VOL*4*3*VOL + (a)*VOL*4*3*VOL + (iv)*4*3*VOL + (nu)*3*VOL + (b)*VOL + ivp )

#define IUNIT ((qkxTMComplex) {0.,1.})
#define RUNIT ((qkxTMComplex) {1.,0.})
#define ZERO ((qkxTMComplex) {0.,0.})
  qkxTMComplex g1[4][4];
  qkxTMComplex g2[4][4];
  qkxTMComplex g3[4][4];
  qkxTMComplex g0[4][4];
  qkxTMComplex g5[4][4];

  for(int mu = 0 ; mu < 4 ; mu++)
    for(int nu = 0 ; nu < 4 ; nu++){
      g1[mu][nu] = ZERO;
      g2[mu][nu] = ZERO;
      g3[mu][nu] = ZERO;
      g0[mu][nu] = ZERO;
      g5[mu][nu] = ZERO;
    }
  /*      
  g1[0][3] = IUNIT;
  g1[1][2] = IUNIT;
  g1[2][1] = -IUNIT;
  g1[3][0] = -IUNIT;

  g2[0][3] = RUNIT;
  g2[1][2] = -RUNIT;
  g2[2][1] = -RUNIT;
  g2[3][0] = RUNIT;

  g3[0][2] = IUNIT;
  g3[1][3] = -IUNIT;
  g3[2][0] = -IUNIT;
  g3[3][1] = IUNIT;

  g4[0][0] = RUNIT;
  g4[1][1] = RUNIT;
  g4[2][2] = -RUNIT;
  g4[3][3] = -RUNIT;

  g5[0][2] = RUNIT;
  g5[1][3] = RUNIT;
  g5[2][0] = RUNIT;
  g5[3][1] = RUNIT;
  */

  g1[0][3] = -IUNIT;
  g1[1][2] = -IUNIT;
  g1[2][1] = IUNIT;
  g1[3][0] = IUNIT;

  g2[0][3] = -RUNIT;
  g2[1][2] = RUNIT;
  g2[2][1] = RUNIT;
  g2[3][0] = -RUNIT;

  g3[0][2] = -IUNIT;
  g3[1][3] = IUNIT;
  g3[2][0] = IUNIT;
  g3[3][1] = -IUNIT;

  g0[0][2] = -RUNIT;
  g0[1][3] = -RUNIT;
  g0[2][0] = -RUNIT;
  g0[3][1] = -RUNIT;

  g5[0][0] = RUNIT;
  g5[1][1] = RUNIT;
  g5[2][2] = -RUNIT;
  g5[3][3] = -RUNIT;


  for(int mu = 0 ; mu < 4 ; mu++)
  for(int a = 0 ; a < 3 ; a++)
  for(int t = 0 ; t < latInfo.L[3] ; t++)
  for(int z = 0 ; z < latInfo.L[2] ; z++)
  for(int y = 0 ; y < latInfo.L[1]; y++)
  for(int x = 0 ; x < latInfo.L[0]; x++)
  for(int nu = 0 ; nu < 4 ; nu++)
  for(int b = 0 ; b < 3 ; b++)
  for(int tp = 0 ; tp < latInfo.L[3] ; tp++)
  for(int zp = 0 ; zp < latInfo.L[2] ; zp++)
  for(int yp = 0 ; yp < latInfo.L[1]; yp++)
  for(int xp = 0 ; xp < latInfo.L[0]; xp++){
    int iv = IV(x,y,z,t);
    int xm1 = ((x-1) + latInfo.L[0])%latInfo.L[0];
    int ym1 = ((y-1) + latInfo.L[1])%latInfo.L[1];
    int zm1 = ((z-1) + latInfo.L[2])%latInfo.L[2];
    int tm1 = ((t-1) + latInfo.L[3])%latInfo.L[3];
    
    int xp1 = ((x+1) )%latInfo.L[0];
    int yp1 = ((y+1) )%latInfo.L[1];
    int zp1 = ((z+1) )%latInfo.L[2];
    int tp1 = ((t+1) )%latInfo.L[3];
    
    int ivmx = IV(xm1,y,z,t);
    int ivmy = IV(x,ym1,z,t);
    int ivmz = IV(x,y,zm1,t);
    int ivmt = IV(x,y,z,tm1);
    
    int ivp = IV(xp,yp,zp,tp);

    A[MA(mu,a,iv,nu,b,ivp)] = -latInfo.kappa *( (DELTA(mu,nu) - g1[mu][nu])*u[0][MG(iv,a,b)]*DELTA(xp1,xp)*DELTA(y,yp)*DELTA(z,zp)*DELTA(t,tp)  + (DELTA(mu,nu) + g1[mu][nu])*conj(u[0][MG(ivmx,b,a)])*DELTA(xm1,xp)*DELTA(y,yp)*DELTA(z,zp)*DELTA(t,tp)  )
      -latInfo.kappa *( (DELTA(mu,nu) - g2[mu][nu])*u[1][MG(iv,a,b)] * DELTA(yp1,yp)*DELTA(x,xp)*DELTA(z,zp)*DELTA(t,tp)  + (DELTA(mu,nu) + g2[mu][nu])*conj(u[1][MG(ivmy,b,a)])*DELTA(ym1,yp)*DELTA(x,xp)*DELTA(z,zp)*DELTA(t,tp)) 	       
      -latInfo.kappa *( (DELTA(mu,nu) - g3[mu][nu])*u[2][MG(iv,a,b)] * DELTA(zp1,zp)*DELTA(x,xp)*DELTA(y,yp)*DELTA(t,tp)  + (DELTA(mu,nu) + g3[mu][nu])*conj(u[2][MG(ivmz,b,a)])*DELTA(zm1,zp)*DELTA(x,xp)*DELTA(y,yp)*DELTA(t,tp))
      -latInfo.kappa *( (DELTA(mu,nu) - g0[mu][nu])*u[3][MG(iv,a,b)] * DELTA(tp1,tp)*DELTA(x,xp)*DELTA(y,yp)*DELTA(z,zp)  + (DELTA(mu,nu) + g0[mu][nu])*conj(u[3][MG(ivmt,b,a)])*DELTA(tm1,tp)*DELTA(x,xp)*DELTA(y,yp)*DELTA(z,zp))
      + ( DELTA(mu,nu) * DELTA(a,b) * DELTA(iv,ivp) + 2.*latInfo.kappa*latInfo.mu*IUNIT*g5[mu][nu]*DELTA(a,b) * DELTA(iv,ivp) );

  }
			  

  qkxTMComplex *AdagA = (qkxTMComplex*)malloc(NL*NL*sizeof(qkxTMComplex));
  qkxTMComplex *Aconj = (qkxTMComplex*)malloc(NL*NL*sizeof(qkxTMComplex));
  qkxTMComplex *w = (qkxTMComplex*)malloc(NL*sizeof(qkxTMComplex));

  for(int i = 0 ; i < NL ; i++)
    for(int j = 0 ; j < NL ; j++)
      Aconj[i*NL+j] = conj(A[i*NL+j]);

  int nl = NL;
  char transA[] = "N";
  char transB[] = "T";
  qkxTMComplex alpha = (qkxTMComplex) {1.,0.};
  qkxTMComplex beta = (qkxTMComplex) {0.,0.};

  zgemm_(transA,transB,&nl,&nl,&nl,&alpha,Aconj,&nl,A,&nl,&beta,AdagA,&nl);

  //  LAPACKE_zgeev(LAPACK_ROW_MAJOR,'N','N',nl,(__complex__ double*) AdagA,nl,(__complex__ double*)w,NULL,nl,NULL,nl);

  char Cvec1[] = "N";
  char Cvec2[] = "V";

  qkxTMComplex wkopt;
  qkxTMComplex *work;
  double *rwork = (double*)malloc(2*NL*sizeof(double));
  int lwork=-1;
  int info;
  zgeev_(Cvec1, Cvec2, &nl, AdagA, &nl, w, NULL, &nl, V, &nl, &wkopt, &lwork, rwork, &info );
  lwork = (int)real(wkopt);
  work = (qkxTMComplex*)malloc(lwork*sizeof(qkxTMComplex));
  zgeev_(Cvec1, Cvec2, &nl, AdagA, &nl, w, NULL, &nl, V, &nl, work, &lwork, rwork, &info );
  //  zgeev_(Cvec1,Cvec2,&nl,AdagA,&nl,w,NULL,&nl,NULL,&nl);

  FILE *ptr;
  ptr = fopen("EIG_chiral.dat","w");
  for(int i = 0 ; i < NVEC ; i++)
    fprintf(ptr,"%+f\n",w[i].real());

  FILE *ptr2;
  ptr2 = fopen("EIGVECTORS_chiral.dat","w");
  for(int i = 0 ; i < NVEC ; i++)
    for(int iv = 0 ; iv < latInfo.L[0]*latInfo.L[1]*latInfo.L[2]*latInfo.L[3]  ; iv++)
      for(int mu = 0 ; mu <4 ; mu++)
	for(int ic = 0 ; ic < 3 ; ic++)
	  fprintf(ptr2,"%+16.15e %+16.15e\n",V[i*4*3*latInfo.L[0]*latInfo.L[1]*latInfo.L[2]*latInfo.L[3] + mu*3*latInfo.L[0]*latInfo.L[1]*latInfo.L[2]*latInfo.L[3] + ic*latInfo.L[0]*latInfo.L[1]*latInfo.L[2]*latInfo.L[3] + iv].real(),V[ i*4*3*latInfo.L[0]*latInfo.L[1]*latInfo.L[2]*latInfo.L[3] + mu*3*latInfo.L[0]*latInfo.L[1]*latInfo.L[2]*latInfo.L[3] + ic*latInfo.L[0]*latInfo.L[1]*latInfo.L[2]*latInfo.L[3] + iv].imag());

  delete gaugefield;

  free(rwork);
  free(work);
  free(w);
  free(AdagA);
  free(A);
  free(Aconj);
  free(V);
 return 0;
}
