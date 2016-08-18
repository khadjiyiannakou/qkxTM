#include <lattice_util.h>

int main(){

  using namespace qkxTM;

  LatticeInfo latInfo;
  latInfo.L[0] = 5;
  latInfo.L[1] = 5;
  latInfo.L[2] = 5;
  latInfo.L[3] = 10;
  latInfo.kappa = 0.161231;
  latInfo.mu = 0.008500;
  latInfo.twistSign = +1;
  latInfo.tol = 1e-8;
  latInfo.maxIter = 1000;
  latInfo.reliableDelta = 0.01;

  

  GaugeField *gaugefield = new GaugeField(&latInfo);

  
  printf("Gauge length is %d\n",gaugefield->TotalBytes());

  gaugefield->create(); // allocates memory for gauge field
  qkxTMComplex **gauge = gaugefield->P_gauge();
  char filename[] = "conf.0000.dat";
  qkxTM_getGaugeAndreas(filename, gauge, &latInfo);


  double plaq = gaugefield->calculatePlaquette();

  /*  
  FILE *kale;
  kale = fopen("kale.dat","w");
  for(int mu = 0 ; mu < 4 ; mu++)
    for(int iv = 0 ; iv < 24*24*24*48 ; iv++)
      for(int c1 = 0 ; c1 < 3 ; c1++)
	for(int c2 = 0 ; c2 < 3 ; c2++)
	  fprintf(kale,"%d %d %d %d \t %e %e\n",mu,iv,c1,c2,gauge[mu][iv*3*3 + c1*3 + c2].real(),gauge[mu][iv*3*3 + c1*3 + c2].imag());
  */
  printf("The plaquette is %8.7f\n",plaq);
  

  /*
  SU3 *A = new SU3();
  SU3 *B = new SU3();
  SU3 *C = new SU3();
  SU3 *D = new SU3();

  // for(int mu = 0 ; mu < 4 ; mu++){
  //  printf("\n");
  //    for(int c1 =0 ; c1 < 3 ; c1++)
  //     for(int c2 = 0 ; c2 < 3 ; c2++)
  //	printf("%e %e \n",gauge[0][MG(0,c1,c2)].real() , gauge[0][MG(0,c1,c2)].imag());
    //}
  

  memcpy(A->M,&(gauge[0][MG(0,0,0)]),9*sizeof(qkxTMComplex));
  memcpy(B->M,&(gauge[1][MG(0,0,0)]),9*sizeof(qkxTMComplex));
  memcpy(C->M,&(gauge[2][MG(0,0,0)]),9*sizeof(qkxTMComplex));

  //  A->PM = &(gauge[0][MG(0,0,0)]);
  // B->PM = &(gauge[1][MG(0,0,0)]);

  printf("\n");
  for(int c1 =0 ; c1 < 3 ; c1++)
    for(int c2 = 0 ; c2 < 3 ; c2++)
      printf("%+e %+e \n",A->M[MU(c1,c2)].real() , A->M[MU(c1,c2)].imag());

  printf("\n");
  for(int c1 =0 ; c1 < 3 ; c1++)
    for(int c2 = 0 ; c2 < 3 ; c2++)
      printf("%+e %+e \n",B->M[MU(c1,c2)].real() , B->M[MU(c1,c2)].imag());

  printf("\n");
  for(int c1 =0 ; c1 < 3 ; c1++)
    for(int c2 = 0 ; c2 < 3 ; c2++)
      printf("%+e %+e \n",C->M[MU(c1,c2)].real() , C->M[MU(c1,c2)].imag());
  

  *D = (*A)*(*B)*(*(C->dagger()));

  qkxTMComplex plaq;
 
  printf("\n");
  for(int c1 =0 ; c1 < 3 ; c1++)
    for(int c2 = 0 ; c2 < 3 ; c2++)
      printf("%+e %+e \n",D->M[MU(c1,c2)].real() , D->M[MU(c1,c2)].imag());

  plaq = D->traceColor();

  std::cout << std::endl;
  std::cout << plaq;

    delete gaugefield;
    delete A;
    delete B;
    delete C;
  */
 return 0;
}
