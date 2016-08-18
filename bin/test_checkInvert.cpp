#include <lattice_util.h>
#include <blas_qkxTM.h>
#include <solvers.h>

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
  latInfo.tol = 1e-6;
  latInfo.maxIter = 2000;
  latInfo.reliableDelta = 0.01;

  latInfo.NsmearAPE = 1;
  latInfo.alphaAPE = 0.5;
  latInfo.NsmearGauss = 1;
  latInfo.alphaGauss = 4.0;

  FILE *ptr_file_in;
  ptr_file_in = fopen("inversion.dat","r");
  FILE *ptr_file_out;
  ptr_file_out = fopen("out.dat","w");


  GaugeField *gaugefield = new GaugeField(&latInfo);
  GaugeField *tmpGauge = new GaugeField(&latInfo);
  GaugeField *gaugefieldAPE = new GaugeField(&latInfo);

  ColorSpinorField *psi = new ColorSpinorField(&latInfo);
  ColorSpinorField *psi_in = new ColorSpinorField(&latInfo);
  ColorSpinorField *psi_out = new ColorSpinorField(&latInfo);

  printf("Gauge length is %d\n",gaugefield->TotalBytes());

  gaugefield->create(); // allocates memory for gauge field
  tmpGauge->create();
  gaugefieldAPE->create();

  psi->create();
  psi_in->create();
  psi_out->create();

  psi->zero();
  psi_out->zero();

  qkxTMComplex **gauge = gaugefield->P_gauge();
  qkxTMComplex **gaugeTmp = tmpGauge->P_gauge();
  qkxTMComplex **gaugeAPE = gaugefieldAPE->P_gauge();

  qkxTMComplex **ptr_psi = psi->P_colorSpinor();
  qkxTMComplex **ptr_psi_in = psi_in->P_colorSpinor();
  qkxTMComplex **ptr_psi_out = psi_out->P_colorSpinor();


  char filename[] = "/home/krikitos3/Desktop/QCD/qkxTM/confs/conf.1000";
  qkxTM_getGaugeLime(filename, gauge, &latInfo);
  //  gaugefield->unitField();
  double plaq = gaugefield->calculatePlaquette();  
  fprintf(stdout,"The plaquette is %e\n",plaq);

  tmpGauge->copy(*gaugefield);


  // fprintf(stdout,"Start APE smearing \n");

  //  gaugefieldAPE->APEsmearing(gaugeTmp,gaugeAPE,&latInfo); // in the object gaugefieldAPE we have the smeared gauge field

  // double plaq_ape = gaugefieldAPE->calculatePlaquette();
  
  //  fprintf(stdout,"The APE plaquette is %e\n",plaq_ape);

  
  //  for(int mu = 0 ; mu < 4 ; mu++)
  //  for(int ic = 0 ; ic < 3 ; ic++){
  //    ptr_psi[MS(mu,ic)][0].real() = 1.;
  //    ptr_psi[MS(mu,ic)][0].imag() = 0.;
  //  }

  //  ptr_psi[0][0].real() = 1.;

  //   psi_sm->gaussianSmearing(ptr_psi,ptr_psi_sm,gaugeAPE,&latInfo);

  int dummy;

  for(int t = 0 ; t < 48 ; t++)
    for(int z = 0 ; z < 24 ; z++)                                                                       
      for(int y = 0 ; y < 24 ; y++)
	for(int x =0 ; x < 24 ; x++){
	  int iv = x + 24*y + 24*24*z + 24*24*24*t;
	  for(int mu = 0 ; mu < 4 ; mu++)
	    for(int ic = 0 ; ic < 3 ; ic++)
	      fscanf(ptr_file_in,"%d %d %d %d %lf %lf\n",&dummy,&dummy,&dummy,&dummy,&(ptr_psi[MS(mu,ic)][iv].real()) , &(ptr_psi[MS(mu,ic)][iv].imag()) );
	}


  // apply D^+ for cg inversion

  //  psi_Dslash->applyDslash(ptr_psi, gauge);
  //  psi_TM->applyTwistAddDslash(ptr_psi_Dslash, ptr_psi, &latInfo);
  psi_out->applyDslash(ptr_psi, gauge);
  psi_out->applyTwistAddDslash(ptr_psi_out,ptr_psi,&latInfo);
  
  //  psi_Dslash->applyDslashDag(ptr_psi, gauge);
  // psi_TM->applyTwistDagAddDslashDag(ptr_psi_Dslash, ptr_psi, &latInfo); 

  // prepare cg object

  //  CG *cg_P = new CG(&latInfo, gauge); // initialize cg object
  //  CG &cg = *cg_P;
  // printf("start cg\n");
  //cg(*psi_out, *psi_in, &latInfo); // preform inversion

  
  for(int t = 0 ; t < 48 ; t++)
    for(int z = 0 ; z < 24 ; z++)                                                                       
      for(int y = 0 ; y < 24 ; y++)
	for(int x =0 ; x < 24 ; x++){
	  int iv = x + 24*y + 24*24*z + 24*24*24*t;
	  for(int mu = 0 ; mu < 4 ; mu++)
	    for(int ic = 0 ; ic < 3 ; ic++)
	      fprintf(ptr_file_out,"%+d %+d %+d %+d %+e %+e\n",x,y,z,t,ptr_psi_out[MS(mu,ic)][iv].real() , ptr_psi_out[MS(mu,ic)][iv].imag());
	}
  
  

 return 0;
}
