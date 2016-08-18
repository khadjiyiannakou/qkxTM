#include <blas_qkxTM.h>
#include <lattice_util.h>






void projectSU3(qkxTMComplex **gauge, LatticeInfo *latInfo){


  SU3 *M = new SU3();
  SU3 *M_dag = new SU3();
  SU3 *H = new SU3();
  SU3 *U = new SU3();
  SU3 *v = new SU3();
  //  SU3 *v_dag = new SU3();
  SU3 *vr = new SU3();

  // constants to be used

  const qkxTMComplex thirdRootOneFirst (1.,sqrt(3.)); 
  const qkxTMComplex thirdRootOneSecond (1.,-sqrt(3.)); 
  const double thirdRoot_12 = pow(12.,1./3.);
  const double thirdRoot_18 = pow(18.,1./3.);
  const double thirdRoot_2_3 = pow((2./3.) , 1./3.);
  const int volume = latInfo->L[0] * latInfo->L[1] * latInfo->L[2] * latInfo->L[3];
  
  qkxTMComplex detM;

  double phase, check_sum;
  double eigenValue[3];
  double de;

  double trace;
  double norm_V;

  double a;
  qkxTMComplex b,w,D;
  qkxTMComplex temp_cmplx , temp_cmplx2 , temp_cmplx3, temp_cmplx4;
  
  for(int iv = 0; iv < volume ; iv++)
    for(int mu = 0 ; mu < 3 ; mu++){                  // only x,y,z not t
      
      memcpy(M->M,&(gauge[mu][MG(iv,0,0)]), 9*sizeof(qkxTMComplex) );
      detM = M->det();

      phase = ( atan2(detM.imag() , detM.real()) )/3.;


      *M_dag = *M;
      M_dag->dagger();
      *H = (*M_dag)*(*M);

      // Assure Hermiticity

      H->M[MU(0,1)].real() = (H->M[MU(0,1)].real() + H->M[MU(1,0)].real())/2.;
      H->M[MU(0,1)].imag() = (H->M[MU(0,1)].imag() - H->M[MU(1,0)].imag())/2.;

      H->M[MU(1,0)] = conj(H->M[MU(0,1)]);

      H->M[MU(0,2)].real() = (H->M[MU(0,2)].real() + H->M[MU(2,0)].real())/2.;
      H->M[MU(0,2)].imag() = (H->M[MU(0,2)].imag() - H->M[MU(2,0)].imag())/2.;

      H->M[MU(2,0)] = conj(H->M[MU(0,2)]);

      H->M[MU(1,2)].real() = (H->M[MU(1,2)].real() + H->M[MU(2,1)].real())/2.;
      H->M[MU(1,2)].imag() = (H->M[MU(1,2)].imag() - H->M[MU(2,1)].imag())/2.;

      H->M[MU(2,1)] = conj(H->M[MU(1,2)]);


      ///////////////////////////////////////////////////////////

      check_sum = norm(H->M[MU(0,1)]) + norm(H->M[MU(0,2)]) + norm(H->M[MU(1,2)]);

      // check if its already diagonal

      if(check_sum <= 1e-08){

	eigenValue[0] = sqrt(1./H->M[MU(0,0)].real());
	eigenValue[1] = sqrt(1./H->M[MU(1,1)].real());
        eigenValue[2] = sqrt(1./H->M[MU(2,2)].real());
	
	U->M[MU(0,0)].real() =  H->M[MU(0,0)].real() *eigenValue[0] ;
	U->M[MU(0,0)].imag() =  H->M[MU(0,0)].imag() *eigenValue[0] ;
	U->M[MU(0,1)].real() =  H->M[MU(0,1)].real() *eigenValue[0] ;
	U->M[MU(0,1)].imag() =  H->M[MU(0,1)].imag() *eigenValue[0] ;
	U->M[MU(0,2)].real() =  H->M[MU(0,2)].real() *eigenValue[0] ;
	U->M[MU(0,2)].imag() =  H->M[MU(0,2)].imag() *eigenValue[0] ;

	U->M[MU(1,0)].real() =  H->M[MU(1,0)].real() *eigenValue[1] ;
	U->M[MU(1,0)].imag() =  H->M[MU(1,0)].imag() *eigenValue[1] ;
	U->M[MU(1,1)].real() =  H->M[MU(1,1)].real() *eigenValue[1] ;
	U->M[MU(1,1)].imag() =  H->M[MU(1,1)].imag() *eigenValue[1] ;
	U->M[MU(1,2)].real() =  H->M[MU(1,2)].real() *eigenValue[1] ;
	U->M[MU(1,2)].imag() =  H->M[MU(1,2)].imag() *eigenValue[1] ;

	U->M[MU(2,0)].real() =  H->M[MU(2,0)].real() *eigenValue[2] ;
	U->M[MU(2,0)].imag() =  H->M[MU(2,0)].imag() *eigenValue[2] ;
	U->M[MU(2,1)].real() =  H->M[MU(2,1)].real() *eigenValue[2] ;
	U->M[MU(2,1)].imag() =  H->M[MU(2,1)].imag() *eigenValue[2] ;
	U->M[MU(2,2)].real() =  H->M[MU(2,2)].real() *eigenValue[2] ;
	U->M[MU(2,2)].imag() =  H->M[MU(2,2)].imag() *eigenValue[2] ;
	
	memcpy(&(gauge[mu][MG(iv,0,0)]),U->M, 9*sizeof(qkxTMComplex) );

      }
      else{

	// do traceless
	trace = ( H->M[MU(0,0)].real() + H->M[MU(1,1)].real() + H->M[MU(2,2)].real() )/3.;
	
	for(int i = 0 ; i<3 ; i++) H->M[MU(i,i)].real() -= trace;

	a = -( H->M[MU(2,2)].real()*H->M[MU(2,2)].real() - H->M[MU(0,0)].real()*H->M[MU(1,1)].real() + ( H->M[MU(0,1)] * conj(H->M[MU(0,1)]) ).real() + ( H->M[MU(0,2)] * conj(H->M[MU(0,2)]) ).real() + ( H->M[MU(1,2)] * conj(H->M[MU(1,2)]) ).real() );

	b.real()  = - H->M[MU(0,0)].real() * H->M[MU(1,1)].real() * H->M[MU(2,2)].real() + H->M[MU(2,2)].real() * (H->M[MU(0,1)] * conj(H->M[MU(0,1)])).real() - (H->M[MU(0,1)] * H->M[MU(1,2)] * conj(H->M[MU(0,2)]) ).real() + H->M[MU(1,1)].real() * ( H->M[MU(0,2)]*conj(H->M[MU(0,2)]) ).real() ;

	b.imag()  =   H->M[MU(2,2)].real()*(H->M[MU(0,1)]*conj(H->M[MU(0,1)])).imag() - (H->M[MU(0,1)]*(H->M[MU(1,2)]*conj(H->M[MU(0,2)]))).imag() + H->M[MU(1,1)].real()*(H->M[MU(0,2)]*conj(H->M[MU(0,2)])).imag();

	b.real() +=   H->M[MU(0,0)].real() * ( H->M[MU(1,2)] * conj(H->M[MU(1,2)]) ).real() - ( H->M[MU(0,2)] * conj( H->M[MU(0,1)] ) * conj( H->M[MU(1,2)] ) ).real() ;
	b.imag() +=   H->M[MU(0,0)].real() * ( H->M[MU(1,2)] * conj(H->M[MU(1,2)]) ).imag() - ( H->M[MU(0,2)] * conj( H->M[MU(0,1)] ) * conj( H->M[MU(1,2)] ) ).imag() ;

	temp_cmplx.real() = 12.0 * a * a * a + 81.0 * (b*b).real();
	temp_cmplx.imag() = 81.0 * (b*b).imag()  ;

	w.real() = ( pow( temp_cmplx , 0.5 ) ).real();
	w.imag() = ( pow( temp_cmplx , 0.5 ) ).imag();

	temp_cmplx2.real() = -9.*b.real() + w.real();
	temp_cmplx2.imag() = -9. * b.imag() + w.imag();
	D = pow( temp_cmplx2 , 1./3. );

	temp_cmplx3.real() = a*thirdRoot_2_3;
	temp_cmplx3.imag() = 0.;

	eigenValue[0] = D.real() / (thirdRoot_18) - ( temp_cmplx3 / D ).real();
 
	temp_cmplx4.real() = D.real() * thirdRoot_12;
	temp_cmplx4.imag() = D.imag() * thirdRoot_12;

	eigenValue[1] = a * ( thirdRootOneFirst / temp_cmplx4 ).real() - (thirdRootOneSecond * D).real() / (thirdRoot_18*2.) ;
	eigenValue[2] = -eigenValue[0]-eigenValue[1];

	eigenValue[0] += trace;
	eigenValue[1] += trace;
	eigenValue[2] += trace;

	H->M[MU(0,0)].real() += trace;
	H->M[MU(1,1)].real() += trace;
	H->M[MU(2,2)].real() += trace;

	// eigen Vectors

	v->M[MU(0,0)].real() = -(eigenValue[0] * H->M[MU(2,0)].real() - H->M[MU(2,0)].real() * H->M[MU(1,1)].real() + (H->M[MU(1,0)] * H->M[MU(2,1)]).real() );
	v->M[MU(0,0)].imag() = -(eigenValue[0] * H->M[MU(2,0)].imag() - H->M[MU(2,0)].imag() * H->M[MU(1,1)].real() + (H->M[MU(1,0)] * H->M[MU(2,1)]).imag() );

	v->M[MU(0,1)].real() = -( (H->M[MU(2,0)] * H->M[MU(0,1)]).real() + eigenValue[0] * H->M[MU(2,1)].real() - H->M[MU(0,0)].real() * H->M[MU(2,1)].real());
	v->M[MU(0,1)].imag() = -( (H->M[MU(2,0)] * H->M[MU(0,1)]).imag() + eigenValue[0] * H->M[MU(2,1)].imag() - H->M[MU(0,0)].real() * H->M[MU(2,1)].imag());

	v->M[MU(0,2)].real() = -eigenValue[0] * eigenValue[0] + eigenValue[0] * H->M[MU(0,0)].real() + ( H->M[MU(0,1)] * conj(H->M[MU(0,1)]) ).real() + eigenValue[0] * H->M[MU(1,1)].real() - H->M[MU(0,0)].real() * H->M[MU(1,1)].real();
	v->M[MU(0,2)].imag() = 0.;    

	v->M[MU(1,0)].real() = -(eigenValue[1] * H->M[MU(2,0)].real() - H->M[MU(2,0)].real() * H->M[MU(1,1)].real() + (H->M[MU(1,0)] * H->M[MU(2,1)]).real() );
	v->M[MU(1,0)].imag() = -(eigenValue[1] * H->M[MU(2,0)].imag() - H->M[MU(2,0)].imag() * H->M[MU(1,1)].real() + (H->M[MU(1,0)] * H->M[MU(2,1)]).imag() );

	v->M[MU(1,1)].real() = -( (H->M[MU(2,0)] * H->M[MU(0,1)]).real() + eigenValue[1] * H->M[MU(2,1)].real() - H->M[MU(0,0)].real() * H->M[MU(2,1)].real());
	v->M[MU(1,1)].imag() = -( (H->M[MU(2,0)] * H->M[MU(0,1)]).imag() + eigenValue[1] * H->M[MU(2,1)].imag() - H->M[MU(0,0)].real() * H->M[MU(2,1)].imag());

	v->M[MU(1,2)].real() =-eigenValue[1] * eigenValue[1] + eigenValue[1] * H->M[MU(0,0)].real() + ( H->M[MU(0,1)] * conj(H->M[MU(0,1)]) ).real() + eigenValue[1] * H->M[MU(1,1)].real() - H->M[MU(0,0)].real() * H->M[MU(1,1)].real();
	v->M[MU(1,2)].imag() = 0.;

	// assure eigen vectors orthogonality

	norm_V  = ( v->M[MU(0,0)] * conj(v->M[MU(0,0)]) ).real() + ( v->M[MU(0,1)] * conj(v->M[MU(0,1)]) ).real()  + ( v->M[MU(0,2)] * conj(v->M[MU(0,2)]) ).real()  ;
	w.real()  = ( v->M[MU(0,0)] * conj(v->M[MU(1,0)]) ).real() + ( v->M[MU(0,1)] * conj(v->M[MU(1,1)]) ).real() + ( v->M[MU(0,2)] * conj(v->M[MU(1,2)] ) ).real()  ;
	w.imag()  = ( v->M[MU(0,0)] * conj(v->M[MU(1,0)]) ).imag() + ( v->M[MU(0,1)] * conj(v->M[MU(1,1)]) ).imag() + ( v->M[MU(0,2)] * conj(v->M[MU(1,2)] ) ).imag()  ;

	w.real() /= norm_V;
	w.imag() /= norm_V;

	v->M[MU(1,0)].real() -= (w * v->M[MU(0,0)]).real();
	v->M[MU(1,0)].imag() -= (w * v->M[MU(0,0)]).imag();

	v->M[MU(1,1)].real() -= (w * v->M[MU(0,1)]).real() ;
	v->M[MU(1,1)].imag() -= (w * v->M[MU(0,1)]).imag();

	v->M[MU(1,2)].real() -= (w * v->M[MU(0,2)]).real();
	v->M[MU(1,2)].imag() -= (w * v->M[MU(0,2)]).imag();

	norm_V=1./sqrt(norm_V);

	// norm_Valize eigenVectors

	v->M[MU(0,0)].real() *= norm_V;
	v->M[MU(0,0)].imag() *= norm_V;

	v->M[MU(0,1)].real() *= norm_V;
	v->M[MU(0,1)].imag() *= norm_V;

	v->M[MU(0,2)].real() *= norm_V;
	v->M[MU(0,2)].imag() *= norm_V;


	norm_V = ( v->M[MU(1,0)] * conj(v->M[MU(1,0)]) ).real() + ( v->M[MU(1,1)] * conj(v->M[MU(1,1)]) ).real() + ( v->M[MU(1,2)] * conj(v->M[MU(1,2)]) ).real();

	norm_V=1./sqrt(norm_V);

	v->M[MU(1,0)].real() *= norm_V;
	v->M[MU(1,0)].imag() *= norm_V;

	v->M[MU(1,1)].real() *= norm_V;
	v->M[MU(1,1)].imag() *= norm_V;

	v->M[MU(1,2)].real() *= norm_V;
	v->M[MU(1,2)].imag() *= norm_V;


	v->M[MU(2,0)].real() = +(v->M[MU(0,1)] * v->M[MU(1,2)]).real() - (v->M[MU(0,2)] * v->M[MU(1,1)]).real(); 
	v->M[MU(2,0)].imag() = -(v->M[MU(0,1)] * v->M[MU(1,2)]).imag() + (v->M[MU(0,2)] * v->M[MU(1,1)]).imag();

	v->M[MU(2,1)].real() = -(v->M[MU(0,0)] * v->M[MU(1,2)]).real() + (v->M[MU(0,2)] * v->M[MU(1,0)]).real();
	v->M[MU(2,1)].imag() = +(v->M[MU(0,0)] * v->M[MU(1,2)]).imag() - (v->M[MU(0,2)] * v->M[MU(1,0)]).imag();

	v->M[MU(2,2)].real() = +(v->M[MU(0,0)] * v->M[MU(1,1)]).real() - (v->M[MU(0,1)] * v->M[MU(1,0)]).real();
	v->M[MU(2,2)].imag() = -(v->M[MU(0,0)] * v->M[MU(1,1)]).imag() + (v->M[MU(0,1)] * v->M[MU(1,0)]).imag();

	//////////////////////////////

	for(int i = 0 ; i < 3 ; i++){
	  de = 1./sqrt(eigenValue[i]);                                                                                                                  
	  b.real() = +de * cos(phase);
	  b.imag() = -de * sin(phase);

	  vr->M[MU(i,0)] = b * v->M[MU(i,0)] ;
	  vr->M[MU(i,1)] = b * v->M[MU(i,1)] ;
	  vr->M[MU(i,2)] = b * v->M[MU(i,2)] ;
	}


	v->dagger();
	*H = (*M)*(*v);
	*U = (*H)*(*vr);

	////////////////////////////////

	norm_V  = ( U->M[MU(0,0)] * conj(U->M[MU(0,0)]) ).real() + ( U->M[MU(1,0)] * conj(U->M[MU(1,0)]) ).real() + ( U->M[MU(2,0)] * conj(U->M[MU(2,0)]) ).real();
	w.real()  = ( U->M[MU(0,0)] * conj(U->M[MU(0,1)]) ).real() + ( U->M[MU(1,0)] * conj(U->M[MU(1,1)]) ).real() + ( U->M[MU(2,0)] * conj(U->M[MU(2,1)]) ).real();
	w.imag()  = ( U->M[MU(0,0)] * conj(U->M[MU(0,1)]) ).imag() + ( U->M[MU(1,0)] * conj(U->M[MU(1,1)]) ).imag() + ( U->M[MU(2,0)] * conj(U->M[MU(2,1)]) ).imag();
	
	w.real() /= norm_V;
	w.imag() /= norm_V;

	for(int i = 0 ; i < 3 ; i++){
	  U->M[MU(i,1)].real() -= (w * U->M[MU(i,0)]).real();
	  U->M[MU(i,1)].imag() -= (w * U->M[MU(i,0)]).imag();
	}
	
	norm_V = 1./sqrt(norm_V);

	U->M[MU(0,0)].real() *= norm_V;
	U->M[MU(0,0)].imag() *= norm_V;
	U->M[MU(1,0)].real() *= norm_V;
	U->M[MU(1,0)].imag() *= norm_V;
	U->M[MU(2,0)].real() *= norm_V;
	U->M[MU(2,0)].imag() *= norm_V;

	norm_V = (U->M[MU(0,1)] * conj(U->M[MU(0,1)])).real() + (U->M[MU(1,1)] * conj(U->M[MU(1,1)])).real() + (U->M[MU(2,1)] * conj(U->M[MU(2,1)])).real();
	norm_V = 1./sqrt(norm_V);

	U->M[MU(0,1)].real() *= norm_V;
	U->M[MU(0,1)].imag() *= norm_V;
	U->M[MU(1,1)].real() *= norm_V;
	U->M[MU(1,1)].imag() *= norm_V;
	U->M[MU(2,1)].real() *= norm_V;
	U->M[MU(2,1)].imag() *= norm_V;


	U->M[MU(0,2)].real() = +(U->M[MU(1,0)] * U->M[MU(2,1)]).real() - (U->M[MU(2,0)] * U->M[MU(1,1)]).real();
	U->M[MU(0,2)].imag() = -(U->M[MU(1,0)] * U->M[MU(2,1)]).imag() + (U->M[MU(2,0)] * U->M[MU(1,1)]).imag();

	U->M[MU(1,2)].real() = -(U->M[MU(0,0)] * U->M[MU(2,1)]).real() + (U->M[MU(2,0)] * U->M[MU(0,1)]).real();
	U->M[MU(1,2)].imag() = +(U->M[MU(0,0)] * U->M[MU(2,1)]).imag() - (U->M[MU(2,0)] * U->M[MU(0,1)]).imag();

	U->M[MU(2,2)].real() = +(U->M[MU(0,0)] * U->M[MU(1,1)]).real() - (U->M[MU(1,0)] * U->M[MU(0,1)]).real();
	U->M[MU(2,2)].imag() = -(U->M[MU(0,0)] * U->M[MU(1,1)]).imag() + (U->M[MU(1,0)] * U->M[MU(0,1)]).imag();


	memcpy(&(gauge[mu][MG(iv,0,0)]),U->M , 9*sizeof(qkxTMComplex) );

	//	if(iv==0)printf("%18.16e\n",gauge[mu][MG(iv,0,0)].real());
	//	exit(-1);

      } // close if statement


    } // close for loop

  delete M;
  delete M_dag;
  delete H;
  delete U;
  delete v;
  //  delete v_dag;
  delete vr;
  
  return;
}



double reCdotXdagY(qkxTM::ColorSpinorField &x , qkxTM::ColorSpinorField &y){

  qkxTMComplex res;

  res.real() = 0.;
  res.imag() = 0.;

  qkxTMComplex **X = x.P_colorSpinor();
  qkxTMComplex **Y = y.P_colorSpinor();

  int volume;
  volume = x.Volume();

  for(int iv = 0 ; iv < volume ; iv++)
    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int ic = 0 ; ic < NCOLORS ; ic++){
	res = res + conj(X[MS(mu,ic)][iv]) * Y[MS(mu,ic)][iv]; 
      }
  
  

  return res.real();
}

void X_eq_X_p_aY(qkxTM::ColorSpinorField &x , double a , qkxTM::ColorSpinorField &y){

  int volume;
  volume = x.Volume();

  qkxTMComplex **X = x.P_colorSpinor();
  qkxTMComplex **Y = y.P_colorSpinor();

  for(int iv = 0 ; iv < volume ; iv++)
    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int ic = 0 ; ic < NCOLORS ; ic++){
	X[MS(mu,ic)][iv] = X[MS(mu,ic)][iv] + a * Y[MS(mu,ic)][iv];
      }

}


void X_eq_aX_p_Y(qkxTM::ColorSpinorField &x , double a , qkxTM::ColorSpinorField &y){

  int volume;
  volume = x.Volume();

  qkxTMComplex **X = x.P_colorSpinor();
  qkxTMComplex **Y = y.P_colorSpinor();

  for(int iv = 0 ; iv < volume ; iv++)
    for(int mu = 0 ; mu < NSPINS ; mu++)
      for(int ic = 0 ; ic < NCOLORS ; ic++){
	X[MS(mu,ic)][iv] = a * X[MS(mu,ic)][iv] + Y[MS(mu,ic)][iv];
      }

}

qkxTMComplex dotProduct( qkxTMComplex *a, qkxTMComplex *b, int n){
  qkxTMComplex res= (qkxTMComplex) {0.,0.};
  for(int i = 0 ; i < n ; i++)
    res = res + conj(a[i]) * b[i];
  return res;
}

void matrixXvec(qkxTMComplex *a, qkxTMComplex *vin, qkxTMComplex *vout,int N){
  memset(vout,0,N*sizeof(qkxTMComplex));
  for(int i = 0 ; i < N ; i++)
    for(int j = 0 ; j < N ; j++)
      vout[i] = vout[i] + a[i*N+j]*vin[j];
}

void matrixXmatrix(qkxTMComplex *a, qkxTMComplex *b, qkxTMComplex *c,int N){
  memset(c,0,N*N*sizeof(qkxTMComplex));
  for(int i = 0 ; i < N ; i++)
    for(int j = 0 ; j < N ; j++)
      for(int k = 0 ; k < N ; k++)
	c[i*N+j] = c[i*N+j] + a[i*N+k]*b[k*N+j];
}

void matrixDagXmatrix(qkxTMComplex *a, qkxTMComplex *b, qkxTMComplex *c,int N){
  memset(c,0,N*N*sizeof(qkxTMComplex));
  for(int i = 0 ; i < N ; i++)
    for(int j = 0 ; j < N ; j++)
      for(int k = 0 ; k < N ; k++)
	c[i*N+j] = c[i*N+j] + conj(a[k*N+i])*b[k*N+j];
}

void vecXmatrix(qkxTMComplex *v1, qkxTMComplex *m, qkxTMComplex *v2,int N){
  memset(v2,0,N*sizeof(qkxTMComplex));
  for(int i = 0 ; i < N ; i++)
    for(int j = 0 ; j < N ; j++)
      v2[i] = v2[i] + v1[j]*m[j*N+i];
}
