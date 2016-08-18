#include <lattice_util.h>
#include <blas_qkxTM.h>

int main(){

  using namespace qkxTM;

  LatticeInfo latInfo;
  latInfo.L[0] = 16;
  latInfo.L[1] = 16;
  latInfo.L[2] = 16;
  latInfo.L[3] = 32;
  latInfo.kappa = 0.160856;
  latInfo.mu = 0.004;

  latInfo.twistSign = +1;
  latInfo.tol = 1e-8;
  latInfo.maxIter = 1000;
  latInfo.reliableDelta = 0.01;

  latInfo.NsmearAPE = 20;
  latInfo.alphaAPE = 0.5;
  latInfo.NsmearGauss = 50;
  latInfo.alphaGauss = 4.0;

  FILE *ptr_test;
  ptr_test = fopen("kale.dat","w");

  GaugeField *gaugefield = new GaugeField(&latInfo);
  GaugeField *tmpGauge = new GaugeField(&latInfo);
  GaugeField *gaugefieldAPE = new GaugeField(&latInfo);

  ColorSpinorField *psi = new ColorSpinorField(&latInfo);
  ColorSpinorField *psi_sm = new ColorSpinorField(&latInfo);

  printf("Gauge length is %d\n",gaugefield->TotalBytes());

  gaugefield->create(); // allocates memory for gauge field
  tmpGauge->create();
  gaugefieldAPE->create();

  psi->create();
  psi_sm->create();

  qkxTMComplex **gauge = gaugefield->P_gauge();
  qkxTMComplex **gaugeTmp = tmpGauge->P_gauge();
  qkxTMComplex **gaugeAPE = gaugefieldAPE->P_gauge();

  qkxTMComplex **ptr_psi = psi->P_colorSpinor();
  qkxTMComplex **ptr_psi_sm = psi_sm->P_colorSpinor();
  
  char filename[] = "/home/khadjiyiannakou/confs/L16T32/conf.1505";
  qkxTM_getGaugeLime(filename, gauge, &latInfo);

  double plaq = gaugefield->calculatePlaquette();  
  fprintf(stdout,"The plaquette is %e\n",plaq);

  tmpGauge->copy(*gaugefield);

  fprintf(stdout,"Start APE smearing \n");

  gaugefieldAPE->APEsmearing(gaugeTmp,gaugeAPE,&latInfo); // in the object gaugefieldAPE we have the smeared gauge field

  double plaq_ape = gaugefieldAPE->calculatePlaquette();
  
  fprintf(stdout,"The APE plaquette is %e\n",plaq_ape);
  ptr_psi[0][0].real() = 1.;
  ptr_psi[0][0].imag() = 0.;

   psi_sm->gaussianSmearing(ptr_psi,ptr_psi_sm,gaugeAPE,&latInfo);


  
    for(int z = 0 ; z < 16 ; z++)                                                                                                                                          
      for(int y = 0 ; y < 16 ; y++)
	for(int x =0 ; x < 16 ; x++){
	  int iv = x + 16*y + 16*16*z;
	  for(int mu = 0 ; mu < 4 ; mu++)
	    for(int ic = 0 ; ic < 3 ; ic++)
	      fprintf(ptr_test,"%d %d %d %e %e\n",x,y,z,ptr_psi_sm[MS(mu,ic)][iv].real() , ptr_psi_sm[MS(mu,ic)][iv].imag());
	}
  
 return 0;
}
