#include <lattice_util.h>
#include <arlib.h>
#include <ctime>
#include <blas_qkxTM.h>

int main(int argc,char *argv[]){

  using namespace qkxTM;

  if(argc != 3){
    fprintf(stderr,"Error wrong number of inputs\n");
    exit(-1);
  }

  LatticeInfo latInfo;
  latInfo.L[0] = 4;
  latInfo.L[1] = 4;
  latInfo.L[2] = 4;
  latInfo.L[3] = 4;
  latInfo.mu = 0.005;
  latInfo.kappa = 0.124843945;
  latInfo.twistSign = +1;

  latInfo.tol = 1e-8;
  latInfo.maxIter = 1000;
  latInfo.reliableDelta = 0.01;

  ArnoldiInfo arInfo;
  arInfo.dimensionMatrix = latInfo.L[0] * latInfo.L[1] * latInfo.L[2] * latInfo.L[3] * 12;
  arInfo.arnoldiTotalVectors = atoi(argv[1]) + atoi(argv[2]);
  arInfo.arnoldiWantedVectors = atoi(argv[1]);
  arInfo.maxArnoldiIter = 1000;
  arInfo.doPolyAccel = false;
  arInfo.tolerance=1e-8;
  arInfo.seed = 12345;


  Arnoldi *X = new Arnoldi(arInfo,latInfo,matrixNxK);
  Arnoldi *W = new Arnoldi(arInfo,latInfo,matrixKxK);

  qkxTMComplex *w = W->P_matrixKxK();

  GaugeField *gaugefield = new GaugeField(&latInfo);
  //  printf("Gauge length is %d\n",gaugefield->TotalBytes());
  gaugefield->create();
  qkxTMComplex **u = gaugefield->P_gauge();
  char filename[] = "conf.4x4x4x4";
  qkxTM_getGaugeAndreas(filename, u, &latInfo);

  FILE *ptr;
  ptr = fopen("gauge.dat","w");
    for(int iv = 0 ; iv < 4*4*4*4 ; iv++)
      for(int mu = 0 ; mu < 4 ; mu++)
      for(int c1 = 0 ; c1 < 3 ; c1++)
	for(int c2 = 0 ; c2 < 3 ; c2++)
	  fprintf(ptr,"%+e %+e\n",u[mu][MG(iv,c1,c2)].real(),u[mu][MG(iv,c1,c2)].imag());

  double plaq = gaugefield->calculatePlaquette();
  printf("The plaquette is %8.7f\n",plaq);
  
  int M = arInfo.arnoldiTotalVectors;
  int K = arInfo.arnoldiWantedVectors;
  //  X->Iram(*X,*W,*gaugefield);

    for(int i = 0 ; i < K ; i++){
      printf("%e %e\n",w[i*M+i].real(),w[i*M+i].imag());
    }



  delete gaugefield;
  delete X;
  delete W;

  //  delete A;
  // delete H;

 return 0;
}
