#include <lattice_util.h>
#include <arlib.h>
#include <ctime>
#include <blas_qkxTM.h>

int main(){

  using namespace qkxTM;

  LatticeInfo latInfo;
  latInfo.L[0] = 4;
  latInfo.L[1] = 4;
  latInfo.L[2] = 4;
  latInfo.L[3] = 4;
  latInfo.mu = 0.005;
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


  Arnoldi *X = new Arnoldi(arInfo,latInfo,matrixNxK);
  Arnoldi *W = new Arnoldi(arInfo,latInfo,matrixKxK);


  qkxTMComplex *x = X->P_matrixNxK();
  qkxTMComplex *w = W->P_matrixKxK();

  GaugeField *gaugefield = new GaugeField(&latInfo);
  printf("Gauge length is %d\n",gaugefield->TotalBytes());
  gaugefield->create();
  qkxTMComplex **u = gaugefield->P_gauge();


  char filename[] = "conf.0000.dat";
  qkxTM_getGaugeAndreas(filename, u, &latInfo);


  double plaq = gaugefield->calculatePlaquette();
  printf("The plaquette is %8.7f\n",plaq);
  
  int M = arInfo.arnoldiTotalVectors;
  int MP1 = M+1;
  int NL = arInfo.dimensionMatrix;
  int K = arInfo.arnoldiWantedVectors;
  //  int VOL = latInfo.L[0] * latInfo.L[1] * latInfo.L[2] * latInfo.L[3];

    X->Iram(*X,*W,*gaugefield);

    for(int i = 0 ; i < K ; i++){
      printf("%e %e\n",w[i*M+i].real(),w[i*M+i].imag());
    }


  /*
#define DELTA(mu,nu) ((mu)==(nu) ? (qkxTMComplex) {1.,0.} : (qkxTMComplex) {0.,0.}) 
#define IV(x,y,z,t)  ( (x) + latInfo.L[0] * (y) + latInfo.L[0]*latInfo.L[1] * (z) + latInfo.L[0]*latInfo.L[1]*latInfo.L[2] * (t))
#define MA(mu,a,iv,nu,b,ivp)  ((mu)*3*VOL*4*3*VOL + (a)*VOL*4*3*VOL + (iv)*4*3*VOL + (nu)*3*VOL + (b)*VOL + ivp )

#define IUNIT ((qkxTMComplex) {0.,1.})
#define RUNIT ((qkxTMComplex) {1.,0.})
#define ZERO ((qkxTMComplex) {0.,0.})
  qkxTMComplex g1[4][4];
  qkxTMComplex g2[4][4];
  qkxTMComplex g3[4][4];
  qkxTMComplex g4[4][4];
  qkxTMComplex g5[4][4];

  for(int mu = 0 ; mu < 4 ; mu++)
    for(int nu = 0 ; nu < 4 ; nu++){
      g1[mu][nu] = ZERO;
      g2[mu][nu] = ZERO;
      g3[mu][nu] = ZERO;
      g4[mu][nu] = ZERO;
      g5[mu][nu] = ZERO;
    }
      
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
    
    //    aa[MA(mu,a,iv,nu,b,ivp)] =-latInfo.kappa * ( (DELTA(mu,nu) - g1[mu][nu])*u[0][MG(iv,a,b)] * DELTA(xp1,xp)  + (DELTA(mu,nu) + g1[mu][nu])*conj(u[0][MG(ivmx,b,a)]) * DELTA(xm1,xp))
    //  -latInfo.kappa * ( (DELTA(mu,nu) - g2[mu][nu])*u[1][MG(iv,a,b)] * DELTA(yp1,yp)  + (DELTA(mu,nu) + g2[mu][nu])*conj(u[1][MG(ivmy,b,a)]) * DELTA(ym1,yp)) 	       
    //  -latInfo.kappa * ( (DELTA(mu,nu) - g3[mu][nu])*u[2][MG(iv,a,b)] * DELTA(zp1,zp)  + (DELTA(mu,nu) + g3[mu][nu])*conj(u[2][MG(ivmz,b,a)]) * DELTA(zm1,zp))
    //  -latInfo.kappa * ( (DELTA(mu,nu) - g4[mu][nu])*u[3][MG(iv,a,b)] * DELTA(tp1,tp)  + (DELTA(mu,nu) + g4[mu][nu])*conj(u[3][MG(ivmt,b,a)]) * DELTA(tm1,tp)) 	       ;
      //  + ( DELTA(mu,nu) * DELTA(a,b) * DELTA(iv,ivp) - 2.*latInfo.kappa*IUNIT*g5[mu][nu]*DELTA(a,b) * DELTA(iv,ivp) );

    aa[MA(mu,a,iv,nu,b,ivp)] = -latInfo.kappa *( (DELTA(mu,nu) - g1[mu][nu])*u[0][MG(iv,a,b)]*DELTA(xp1,xp)*DELTA(y,yp)*DELTA(z,zp)*DELTA(t,tp)  + (DELTA(mu,nu) + g1[mu][nu])*conj(u[0][MG(ivmx,b,a)])*DELTA(xm1,xp)*DELTA(y,yp)*DELTA(z,zp)*DELTA(t,tp)  )
      -latInfo.kappa *( (DELTA(mu,nu) - g2[mu][nu])*u[1][MG(iv,a,b)] * DELTA(yp1,yp)*DELTA(x,xp)*DELTA(z,zp)*DELTA(t,tp)  + (DELTA(mu,nu) + g2[mu][nu])*conj(u[1][MG(ivmy,b,a)])*DELTA(ym1,yp)*DELTA(x,xp)*DELTA(z,zp)*DELTA(t,tp)) 	       
      -latInfo.kappa *( (DELTA(mu,nu) - g3[mu][nu])*u[2][MG(iv,a,b)] * DELTA(zp1,zp)*DELTA(x,xp)*DELTA(y,yp)*DELTA(t,tp)  + (DELTA(mu,nu) + g3[mu][nu])*conj(u[2][MG(ivmz,b,a)])*DELTA(zm1,zp)*DELTA(x,xp)*DELTA(y,yp)*DELTA(t,tp))
      -latInfo.kappa *( (DELTA(mu,nu) - g4[mu][nu])*u[3][MG(iv,a,b)] * DELTA(tp1,tp)*DELTA(x,xp)*DELTA(y,yp)*DELTA(z,zp)  + (DELTA(mu,nu) + g4[mu][nu])*conj(u[3][MG(ivmt,b,a)])*DELTA(tm1,tp)*DELTA(x,xp)*DELTA(y,yp)*DELTA(z,zp))
      + ( DELTA(mu,nu) * DELTA(a,b) * DELTA(iv,ivp) + 2.*latInfo.kappa*latInfo.mu*IUNIT*g5[mu][nu]*DELTA(a,b) * DELTA(iv,ivp) );

  }
			  
  qkxTMComplex *AdagA = (qkxTMComplex*)malloc(NL*NL*sizeof(qkxTMComplex));
  //matrixDagXmatrix(aa,aa,AdagA,NL);
    //memcpy(aa,AdagA,NL*NL*sizeof(qkxTMComplex));

    */
  /*    
  FILE *ptr;
  ptr = fopen("kale.dat","w");

    for(int mu = 0 ; mu < 4 ; mu++)
      for(int a = 0 ; a < 3 ; a++)
	for(int t = 0 ; t < latInfo.L[3] ; t++)
	  for(int z = 0 ; z < latInfo.L[2] ; z++)
	    for(int y = 0 ; y < latInfo.L[1]; y++)
	      for(int x = 0 ; x < latInfo.L[0]; x++){
		int iv = IV(x,y,z,t);
		fprintf(ptr,"%d %d %d %d %d %d %+e %+e\n",mu,a,t,z,y,x,aa[MA(mu,a,iv,0,0,0)].real(), aa[MA(mu,a,iv,0,0,0)].imag());
    }
  */;
  /*
    for(int nu = 0 ; nu < 4 ; nu++)
      for(int b = 0 ; b < 3 ; b++)
	for(int tp = 0 ; tp < latInfo.L[3] ; tp++)
	  for(int zp = 0 ; zp < latInfo.L[2] ; zp++)
	    for(int yp = 0 ; yp < latInfo.L[1]; yp++)
	      for(int xp = 0 ; xp < latInfo.L[0]; xp++){
		int ivp = IV(xp,yp,zp,tp);
		fprintf(ptr,"%d %d %d %d %d %d %+e %+e\n",nu,b,tp,zp,yp,xp,aa[MA(0,0,0,nu,b,ivp)].real(), aa[MA(0,0,0,nu,b,ivp)].imag());
	      }
  */
   //    A->Iram(*X,*H,*W);



  /*
  qkxTMComplex *AdagA = (qkxTMComplex*)malloc(NL*NL*sizeof(qkxTMComplex));

  for(int i = 0 ; i < NL ; i++)
    for(int j = 0 ; j < NL ; j++){
      if( (i==j) || (i == j-1) || (i == j+1))
	a[i*NL+j] = (qkxTMComplex) {drand48()-0.5,drand48()-0.5};
      else
	a[i*NL+j] = (qkxTMComplex) {0,0.};
    }

  matrixDagXmatrix(a,a,AdagA,NL);
  memcpy(a,AdagA,NL*NL*sizeof(qkxTMComplex));


  A->Iram(*X,*H,*W);
  */
  /*  
  for(int i = 0 ; i < NL ; i++){
    for(int j = 0 ; j < K ; j++)
      printf("%f %f\t",x[i*MP1+j].real(),x[i*MP1+j].imag());
    printf("\n");
      }

  printf("\n");
  */



  /*
  H->P_matrixKxK()[0*M+0] = (qkxTMComplex) {+2.769230e-01,+7.093648e-01};
  H->P_matrixKxK()[0*M+1] = (qkxTMComplex) {+6.948286e-01,+6.550980e-01};
  H->P_matrixKxK()[0*M+2] = (qkxTMComplex) {+4.387444e-01,+9.597440e-01};
  H->P_matrixKxK()[0*M+3] = (qkxTMComplex) {+1.868726e-01,+7.512671e-01};
  H->P_matrixKxK()[1*M+0] = (qkxTMComplex) {+4.617139e-02,+7.546867e-01};
  H->P_matrixKxK()[1*M+1] = (qkxTMComplex) {+3.170995e-01,+1.626117e-01};
  H->P_matrixKxK()[1*M+2] = (qkxTMComplex) {+3.815585e-01,+3.403857e-01};
  H->P_matrixKxK()[1*M+3] = (qkxTMComplex) {+4.897644e-01,+2.550951e-01};
  H->P_matrixKxK()[2*M+0] = (qkxTMComplex) {+9.713178e-02,+2.760251e-01};
  H->P_matrixKxK()[2*M+1] = (qkxTMComplex) {+9.502220e-01,+1.189977e-01};
  H->P_matrixKxK()[2*M+2] = (qkxTMComplex) {+7.655168e-01,+5.852678e-01};
  H->P_matrixKxK()[2*M+3] = (qkxTMComplex) {+4.455862e-01,+5.059571e-01};
  H->P_matrixKxK()[3*M+0] = (qkxTMComplex) {+8.234578e-01,+6.797027e-01};
  H->P_matrixKxK()[3*M+1] = (qkxTMComplex) {+3.444608e-02,+4.983641e-01};
  H->P_matrixKxK()[3*M+2] = (qkxTMComplex) {+7.951999e-01,+2.238119e-01};
  H->P_matrixKxK()[3*M+3] = (qkxTMComplex) {+6.463130e-01,+6.990767e-01};
    

  
  H->QRLapack(*Q,*W);
  for(int i = 0 ; i < K ; i++){
    for(int j = 0 ; j < K ; j++)
      printf("%f %f\t",w[i*M+j].real(),w[i*M+j].imag());
    printf("\n");
      }

  printf("\n");

  for(int i = 0 ; i < K ; i++){
    for(int j = 0 ; j < K ; j++)
      printf("%f %f\t",q[i*M+j].real(),q[i*M+j].imag());
    printf("\n");
      }

  */
  /*
  for(int i = 0 ; i < NL ; i++)
    for(int j = 0 ; j < NL ; j++){
      if( (i==j) || (i == j-1) || (i == j+1))
	a[i*NL+j] = (qkxTMComplex) {1.,0};
      else
	a[i*NL+j] = (qkxTMComplex) {0,0.};
    }

  A->Iram(*V,*H,*W);
  
  for(int i = 0 ; i < NL ; i++){
    for(int j = 0 ; j < K ; j++)
      printf("%f %f\t",v[i*MP1+j].real(),v[i*MP1+j].imag());
    printf("\n");
      }

  printf("\n");

  for(int i = 0 ; i < K ; i++){
    for(int j = 0 ; j < K ; j++)
      printf("%f %f\t",w[i*M+j].real(),w[i*M+j].imag());
    printf("\n");
      }
  */
  /*
  h[0*M+0] = (qkxTMComplex) {-3.84,  2.25};
  h[0*M+1] = (qkxTMComplex) {-8.94, -4.75};
  h[0*M+2] = (qkxTMComplex) {8.95, -6.53};
  h[0*M+3] = (qkxTMComplex) {-9.87,  4.82};

  h[1*M+0] = (qkxTMComplex) {-0.66,  0.83};
  h[1*M+1] = (qkxTMComplex) {-4.40, -3.82};
  h[1*M+2] = (qkxTMComplex) {-3.50, -4.26};
  h[1*M+3] = (qkxTMComplex) {-3.15,  7.36};

  h[2*M+0] = (qkxTMComplex) {-3.99, -4.73};
  h[2*M+1] = (qkxTMComplex) {-5.88, -6.60};
  h[2*M+2] = (qkxTMComplex) {-3.36, -0.40};
  h[2*M+3] = (qkxTMComplex) {-0.75,  5.23};

  h[3*M+0] = (qkxTMComplex) {7.74,  4.18};
  h[3*M+1] = (qkxTMComplex) {3.66, -7.53};
  h[3*M+2] = (qkxTMComplex) {2.58,  3.60};
  h[3*M+3] = (qkxTMComplex) {4.59,  5.41};
  */
  //  H->eigLapack(*Q,*W);

  /*
  matrixDagXmatrix(a,a,AdagA,NL);
  memcpy(a,AdagA,NL*NL*sizeof(qkxTMComplex));

  for(int i = 0 ; i < NL ; i++){
    for(int j = 0 ; j < NL ; j++)
      printf("(%3.1f) (%3.1f)i\t",a[i*NL+j].real(),a[i*NL+j].imag());
    printf("\n");
    }
  */

  //  A->Iram(*V,*H,*W);

  



  /*
  h[0*M+0] = (qkxTMComplex) {0.6324 , 0.4854};
  h[0*M+1] = (qkxTMComplex) {0.5469 , 0.4218};
  h[0*M+2] = (qkxTMComplex) {0.1576 , 0.9595};
  h[1*M+0] = (qkxTMComplex) {0.0975 , 0.8003};
  h[1*M+1] = (qkxTMComplex) {0.9575 , 0.9157};
  h[1*M+2] = (qkxTMComplex) {0.9706 , 0.6557};
  h[2*M+0] = (qkxTMComplex) {0.,0.};
  h[2*M+1] = (qkxTMComplex) {0.9649 , 0.7922};
  h[2*M+2] = (qkxTMComplex) {0.9572 , 0.0357};

  v[0*MP1+0] = (qkxTMComplex) {0.8491 , 0.7655};
  v[0*MP1+1] = (qkxTMComplex) {0.3922 , 0.6463};
  v[0*MP1+2] = (qkxTMComplex) {0.2769 , 0.6551};
  v[0*MP1+3] = (qkxTMComplex) {0.3171 , 0.3404};

  v[1*MP1+0] = (qkxTMComplex) {0.9340 , 0.7952};
  v[1*MP1+1] = (qkxTMComplex) {0.6555 , 0.7094};
  v[1*MP1+2] = (qkxTMComplex) {0.0462 , 0.1626};
  v[1*MP1+3] = (qkxTMComplex) {0.9502 , 0.5853};

  v[2*MP1+0] = (qkxTMComplex) {0.6787 , 0.1869};
  v[2*MP1+1] = (qkxTMComplex) {0.1712 , 0.7547};
  v[2*MP1+2] = (qkxTMComplex) {0.0971 , 0.1190};
  v[2*MP1+3] = (qkxTMComplex) {0.0344 , 0.2238};

  v[3*MP1+0] = (qkxTMComplex) {0.7577 , 0.4898};
  v[3*MP1+1] = (qkxTMComplex) {0.7060 , 0.2760};
  v[3*MP1+2] = (qkxTMComplex) {0.8235 , 0.4984};
  v[3*MP1+3] = (qkxTMComplex) {0.4387 , 0.7513};

  v[4*MP1+0] = (qkxTMComplex) {0.7431 , 0.4456};
  v[4*MP1+1] = (qkxTMComplex) {0.0318 , 0.6797};
  v[4*MP1+2] = (qkxTMComplex) {0.6948 , 0.9597};
  v[4*MP1+3] = (qkxTMComplex) {0.3816 , 0.2551};

  memset(w,0,M*M*sizeof(qkxTMComplex));

  w[0*M+0] = (qkxTMComplex) {0.8407 , 0.2435};
  w[1*M+1] = (qkxTMComplex) {0.2543 , 0.9293};
  w[2*M+2] = (qkxTMComplex) {0.8143 , 0.3500};


  for(int i = 0 ; i < NL ; i++){
    for(int j = 0 ; j < MP1 ; j++)
      printf("%f %f\t",V->P_matrixNxK()[i*MP1 + j].real(),V->P_matrixNxK()[i*MP1 + j].imag());
    printf("\n");
  }      
  printf("\n");

  for(int i = 0 ; i < M ; i++){
    for(int j = 0 ; j < M ; j++)
      printf("%f %f\t",H->P_matrixKxK()[i*M + j].real(),H->P_matrixKxK()[i*M + j].imag());
    printf("\n");
  }      
  printf("\n");

  for(int i = 0 ; i < M ; i++){
    for(int j = 0 ; j < M ; j++)
      printf("%f %f\t",W->P_matrixKxK()[i*M + j].real(),W->P_matrixKxK()[i*M + j].imag());
    printf("\n");
  }      
  printf("\n");

  */
  //  sortEigenValues(*H,*W);
  //qrShiftsRotations_v2(*V,*H,*W);

  /*
  for(int i = 0; i < arInfo.dimensionMatrix ; i++)
    for(int j = 0; j < arInfo.dimensionMatrix ; j++){
      A->P_matrixKxK()[i*arInfo.dimensionMatrix+j].real() = rand()/(double) RAND_MAX;
      A->P_matrixKxK()[i*arInfo.dimensionMatrix+j].imag() = rand()/(double) RAND_MAX;
    }

  A->eig(*X,*W);


  //  initArnold(*A,*V,*H);
  //  arnold(1,arInfo.arnoldiWantedVectors,*A,*V,*H);
  // arnold(arInfo.arnoldiWantedVectors,arInfo.arnoldiTotalVectors,*A,*V,*H);

    printf("A = \n");
    for(int i = 0; i < arInfo.dimensionMatrix ; i++){
      for(int j = 0; j < arInfo.dimensionMatrix ; j++){
	printf("%f %f\t",A->P_matrixKxK()[i*arInfo.dimensionMatrix+j].real(), A->P_matrixKxK()[i*arInfo.dimensionMatrix+j].imag());
      }
      printf("\n");
    }


    printf("X = \n");
    for(int i = 0; i < arInfo.dimensionMatrix ; i++){
      for(int j = 0; j < arInfo.dimensionMatrix ; j++){
	printf("%f %f\t",X->P_matrixKxK()[i*(arInfo.dimensionMatrix)+j].real(), X->P_matrixKxK()[i*(arInfo.dimensionMatrix)+j].imag());
      }
      printf("\n");
    }

    printf("W = \n");
    for(int i = 0; i < arInfo.dimensionMatrix ; i++){
      for(int j = 0; j < arInfo.dimensionMatrix ; j++){
	printf("%f %f\t",W->P_matrixKxK()[i*(arInfo.dimensionMatrix)+j].real(), W->P_matrixKxK()[i*(arInfo.dimensionMatrix)+j].imag());
      }
      printf("\n");
    }

    qkxTMComplex R[2];
    qkxTMComplex Q[4];
    qkxTMComplex x[2];

    x[0] = (qkxTMComplex) {0.8147,0.1270};
    x[1] = (qkxTMComplex) {0.9058,0.9134};

    QR_hs_v2(Q,R,x);

    for(int i = 0 ; i < 4 ; i++)
      printf("%f %f\n",Q[i].real(),Q[i].imag());

    printf("\n");
    for(int i = 0 ; i < 2 ; i++)
      printf("%f %f\n",R[i].real(),R[i].imag());


    printf("Q = \n");
    for(int i = 0; i < arInfo.arnoldiTotalVectors ; i++){
      for(int j = 0; j < arInfo.arnoldiTotalVectors ; j++){
	printf("%f %f\t",Q->P_matrixKxK()[i*arInfo.arnoldiTotalVectors+j].real(), Q->P_matrixKxK()[i*arInfo.arnoldiTotalVectors+j].imag());
      }
      printf("\n");
    }

    FILE *ptr_out;
    ptr_out = fopen("kale.dat","w");
    for(int i = 0; i < arInfo.arnoldiTotalVectors ; i++){
      for(int j = 0; j < arInfo.arnoldiTotalVectors ; j++){
	fprintf(ptr_out,"%f %f ",Q->P_matrixKxK()[i*arInfo.arnoldiTotalVectors+j].real(), Q->P_matrixKxK()[i*arInfo.arnoldiTotalVectors+j].imag());
      }
      fprintf(ptr_out,"\n");
    }

  */

  delete gaugefield;
  delete X;
  delete W;

  //  delete A;
  // delete H;

 return 0;
}
