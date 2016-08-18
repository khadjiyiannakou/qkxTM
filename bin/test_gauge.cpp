#include <lattice_util.h>

int main(){

  using namespace qkxTM;

  LatticeInfo latInfo;
  latInfo.L[0] = 24;
  latInfo.L[1] = 24;
  latInfo.L[2] = 24;
  latInfo.L[3] = 48;
  latInfo.kappa = 0.161231;
  latInfo.mu = 0.008500;
  latInfo.twistSign = +1;
  latInfo.tol = 1e-8;
  latInfo.maxIter = 1000;
  latInfo.reliableDelta = 0.01;

  

  GaugeField gaugefield(&latInfo);

  printf("Gauge length is %d\n",gaugefield.TotalBytes());

  gaugefield.create(); // allocates memory for gauge field
  qkxTMComplex **gauge = gaugefield.P_gauge();
  char filename[] = "/home/krikitos3/Desktop/QCD/qkxTM/confs/conf.1000";

  qkxTM_getGaugeLime(filename, gauge, &latInfo);
  int  volume = gaugefield.Volume();

  for(int mu = 0 ; mu < 4 ; mu++){
    printf("\n");
    for(int c1 =0 ; c1 < 3 ; c1++)
      for(int c2 = 0 ; c2 < 3 ; c2++)
	printf("%e %e \n",gauge[mu][MG(volume-1,c1,c2)].real() , gauge[mu][MG(volume-1,c1,c2)].imag());
    }
 return 0;
}
