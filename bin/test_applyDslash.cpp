#include <lattice_util.h>
#include <blas_qkxTM.h>

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

  latInfo.NsmearAPE = 1;
  latInfo.alphaAPE = 0.5;
  latInfo.NsmearGauss = 1;
  latInfo.alphaGauss = 4.0;

  FILE *ptr_test;
  ptr_test = fopen("kale.dat","w");


  GaugeField *gaugefield = new GaugeField(&latInfo);
  GaugeField *tmpGauge = new GaugeField(&latInfo);
  GaugeField *gaugefieldAPE = new GaugeField(&latInfo);

  ColorSpinorField *psi = new ColorSpinorField(&latInfo);
  ColorSpinorField *psi_Dslash = new ColorSpinorField(&latInfo);
  ColorSpinorField *psi_TM = new ColorSpinorField(&latInfo);

  printf("Gauge length is %d\n",gaugefield->TotalBytes());

  gaugefield->create(); // allocates memory for gauge field
  tmpGauge->create();
  gaugefieldAPE->create();

  psi->create();
  psi_Dslash->create();
  psi_TM->create();

  qkxTMComplex **gauge = gaugefield->P_gauge();
  qkxTMComplex **gaugeTmp = tmpGauge->P_gauge();
  qkxTMComplex **gaugeAPE = gaugefieldAPE->P_gauge();

  qkxTMComplex **ptr_psi = psi->P_colorSpinor();
  qkxTMComplex **ptr_psi_Dslash = psi_Dslash->P_colorSpinor();
  qkxTMComplex **ptr_TM = psi_TM->P_colorSpinor();


  char filename[] = "/home/krikitos3/Desktop/QCD/qkxTM/confs/conf.1000";
  qkxTM_getGaugeLime(filename, gauge, &latInfo);
  //  gaugefield->unitField();
  double plaq = gaugefield->calculatePlaquette();  
  fprintf(stdout,"The plaquette is %e\n",plaq);

  tmpGauge->copy(*gaugefield);

  /*
  for(int mu = 0 ; mu < 4 ;mu++)
    for(int it = 0 ; it < 48 ; it++)
      for(int iz = 0 ; iz < 24 ; iz++)
	for(int iy = 0 ; iy < 24 ; iy++)
	  for(int ix = 0 ; ix < 24 ; ix++)
	    for(int c1 = 0 ; c1 < 3 ; c1++)
	      for(int c2 = 0 ; c2 < 3 ; c2++){
		int iv = ix + 24*iy + 24*24*iz + 24*24*24*it;
		fprintf(ptr_test," %d %d %d %d %d %e %e\n",mu,ix,iy,iz,it, gauge[mu][MG(iv,c1,c2)].real() , gauge[mu][MG(iv,c1,c2)].imag() );
	      }
  */

  // fprintf(stdout,"Start APE smearing \n");

  //  gaugefieldAPE->APEsmearing(gaugeTmp,gaugeAPE,&latInfo); // in the object gaugefieldAPE we have the smeared gauge field

  // double plaq_ape = gaugefieldAPE->calculatePlaquette();
  
  //  fprintf(stdout,"The APE plaquette is %e\n",plaq_ape);

  
  for(int mu = 0 ; mu < 4 ; mu++)
    for(int ic = 0 ; ic < 3 ; ic++){
      ptr_psi[MS(mu,ic)][0].real() = 1.;
      ptr_psi[MS(mu,ic)][0].imag() = 0.;
    }
  //   psi_sm->gaussianSmearing(ptr_psi,ptr_psi_sm,gaugeAPE,&latInfo);

  psi_Dslash->applyDslash(ptr_psi, gauge);
  psi_TM->applyTwistAddDslash(ptr_psi_Dslash, ptr_psi, &latInfo);
  psi_Dslash->applyDslashDag(ptr_TM, gauge);
  psi->applyTwistDagAddDslashDag(ptr_psi_Dslash,ptr_TM,&latInfo);
  
  //  psi_Dslash->applyDslashDag(ptr_psi, gauge);
  // psi_TM->applyTwistDagAddDslashDag(ptr_psi_Dslash, ptr_psi, &latInfo); 

  for(int t = 0 ; t < 48 ; t++)
    for(int z = 0 ; z < 24 ; z++)                                                                       
      for(int y = 0 ; y < 24 ; y++)
	for(int x =0 ; x < 24 ; x++){
	  int iv = x + 24*y + 24*24*z + 24*24*24*t;
	  for(int mu = 0 ; mu < 4 ; mu++)
	    for(int ic = 0 ; ic < 3 ; ic++)
	      fprintf(ptr_test,"%+d %+d %+d %+d %+e %+e\n",x,y,z,t,ptr_psi[MS(mu,ic)][iv].real() , ptr_psi[MS(mu,ic)][iv].imag());
	}

  

 return 0;
}
