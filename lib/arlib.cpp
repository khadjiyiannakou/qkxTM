#include <qkxTM.h>
#include <lattice_util.h>
#include <arlib.h>
#include <blas_qkxTM.h>
#include <lapacke.h>
#include <ctime>
#include <limits>



using namespace qkxTM;

Arnoldi::Arnoldi(ArnoldiInfo in_arInfo,LatticeInfo in_latInfo, enum MATRIXTYPE matrixType):
  LatticeField(&in_latInfo), p_matrixNxN(NULL),p_matrixNxK(NULL),p_vectorN(NULL),p_matrixKxK(NULL)
{
  latInfo = in_latInfo;
  arInfo = in_arInfo;

  dimensionMatrix = arInfo.dimensionMatrix;
  arnoldiTotalVectors = arInfo.arnoldiTotalVectors;
  arnoldiWantedVectors = arInfo.arnoldiWantedVectors;
  arnoldiUnwantedVectors = arnoldiTotalVectors - arnoldiWantedVectors;
  maxArnoldiIter = arInfo.maxArnoldiIter;
  doPolyAccel = arInfo.doPolyAccel;
  tolerance = arInfo.tolerance;
  seed = arInfo.seed;

  if(matrixType == matrixNxN){
    p_matrixNxN = (qkxTMComplex*)malloc(dimensionMatrix*dimensionMatrix*sizeof(qkxTMComplex));
    if(p_matrixNxN == NULL){
      fprintf(stderr,"Error allocate matrixNxN\n");
      exit(-1);
    }
  }


  if(matrixType == matrixNxK){
    p_matrixNxK = (qkxTMComplex*)malloc(dimensionMatrix*(arnoldiTotalVectors+1)*sizeof(qkxTMComplex)); // last column holds residual
    if(p_matrixNxK == NULL){
      fprintf(stderr,"Error allocate matrixNxK\n");
      exit(-1);
    }    
  }

  if(matrixType == matrixKxK){
    p_matrixKxK = (qkxTMComplex*)malloc(arnoldiTotalVectors*arnoldiTotalVectors*sizeof(qkxTMComplex));
    if(p_matrixKxK == NULL){
      fprintf(stderr,"Error allocate matrixKxK\n");
      exit(-1);
    }    
  }

  if(matrixType == vectorN){
    p_vectorN = (qkxTMComplex*)malloc(dimensionMatrix*sizeof(qkxTMComplex));
    if(p_vectorN == NULL){
      fprintf(stderr,"Error allocate vectorN\n");
      exit(-1);
    }      
  }

}

Arnoldi::~Arnoldi(){
  if(p_matrixNxN != NULL) free(p_matrixNxN);
  if(p_matrixNxK != NULL) free(p_matrixNxK);
  if(p_matrixKxK != NULL) free(p_matrixKxK);
  if(p_vectorN != NULL) free(p_vectorN);
}

static void GsProjection(qkxTMComplex *a, qkxTMComplex *Q,qkxTMComplex *p,int k,int n){
  //  p <a,e>/<e,e> * e , e must be normalize so <e,e> = 1
  qkxTMComplex *e = (qkxTMComplex*)calloc(n,sizeof(qkxTMComplex));
  memset(p,0,n*sizeof(qkxTMComplex));
  for(int j = 0 ; j < k ; j++){
    for(int i = 0 ; i < n ; i++) e[i] = Q[i*n+j];
    qkxTMComplex dot = dotProduct(a,e,n);
    for(int i = 0 ; i < n ; i++)
      p[i] = p[i] + dot * e[i];
    
  }
  free(e);
}

void Arnoldi::QR_gs(Arnoldi &Q,Arnoldi &R){
  int N = arnoldiTotalVectors;
  qkxTMComplex *p = (qkxTMComplex*)calloc(N,sizeof(qkxTMComplex));
  qkxTMComplex *e = (qkxTMComplex*)calloc(N,sizeof(qkxTMComplex));
  qkxTMComplex *a = (qkxTMComplex*)calloc(N,sizeof(qkxTMComplex));
  qkxTMComplex *u = (qkxTMComplex*)calloc(N,sizeof(qkxTMComplex));

  qkxTMComplex *A = this->P_matrixKxK();
  qkxTMComplex *q = Q.P_matrixKxK();
  qkxTMComplex *r = R.P_matrixKxK();

  for(int k =0 ; k < N ; k++){ // project using Household all columns
    if( k == 0){
      for(int i = 0 ; i < N ; i++) a[i] = A[i*N+k];
      for(int i = 0 ; i < N ; i++) u[i] = a[i];
      double mangitute = sqrt(dotProduct(u,u,N).real());
      for(int i = 0 ; i < N ; i++) e[i] = u[i]*(1./mangitute);
      for(int i = 0 ; i < N ; i++) q[i*N+k] = e[i]; 
    }
    else{
  
      for(int i = 0 ; i < N ; i++) a[i] = A[i*N+k];
      GsProjection(a,q,p,k,N);
      for(int i = 0 ; i < N ; i++) u[i] = a[i] - p[i];
      double mangitute = sqrt(dotProduct(u,u,N).real());
      for(int i = 0 ; i < N ; i++) e[i] = u[i]*(1./mangitute);
      for(int i = 0 ; i < N ; i++) q[i*N+k] = e[i];       
  
    }
    
  }

  memset(r,0,N*N*sizeof(qkxTMComplex));
  for(int i = 0 ; i < N ; i++)
    for(int j = 0 ; j < N ; j++)
      for(int l = 0 ; l < N ; l++)
	r[i*N+j] = r[i*N+j] + conj(q[l*N+i])*this->P_matrixKxK()[l*N+j];
  

  free(p);
  free(e);
  free(a);
  free(u);
}


void Arnoldi::QR_gs_m(Arnoldi &Q,Arnoldi &R){
  
  int N = arnoldiTotalVectors;
  qkxTMComplex *tmp = (qkxTMComplex*)malloc(N*N*sizeof(qkxTMComplex));
  memcpy(tmp, this->P_matrixKxK(),N*N*sizeof(qkxTMComplex));
  qkxTMComplex *a = (qkxTMComplex*)calloc(N,sizeof(qkxTMComplex));
  qkxTMComplex *qr = (qkxTMComplex*)calloc(N,sizeof(qkxTMComplex));
  qkxTMComplex *q = Q.P_matrixKxK();
  qkxTMComplex *r = R.P_matrixKxK();
  
  for(int k = 0 ; k < N ; k++){
    for(int i=0 ; i < N ; i++) a[i] = tmp[i*N+k];
    double alpha = sqrt(dotProduct(a,a,N).real());
    for(int i=0 ; i < N ; i++) q[i*N+k] = tmp[i*N+k]*(1./alpha);
    for(int j = k+1 ; j <N ; j++){
      for(int i=0 ; i < N ; i++) a[i] = tmp[i*N+j];
      for(int i=0 ; i < N ; i++) qr[i] = q[i*N+k];
      qkxTMComplex beta = dotProduct(qr,a,N);
      for(int i=0 ; i < N ; i++) tmp[i*N+j] = tmp[i*N+j] - beta*q[i*N+k];
    } 
  }
  
  
  memset(r,0,N*N*sizeof(qkxTMComplex));

  for(int i = 0 ; i < N ; i++)
    for(int j = 0 ; j < N ; j++)
      for(int l = 0 ; l < N ; l++)
	r[i*N+j] = r[i*N+j] + conj(q[l*N+i])*this->P_matrixKxK()[l*N+j];

  free(a);
  free(qr);
  free(tmp);
}


void qkxTM::QR_hs_v2(qkxTMComplex *Q, qkxTMComplex *R, qkxTMComplex *x){
  // this QR works only in a special case where we want to transform a vector with 2 elements
  qkxTMComplex v[2];
  qkxTMComplex u[2];
  qkxTMComplex w;

    qkxTMComplex alpha = - exp( (qkxTMComplex) {0.,atan2(imag(x[0]),real(x[0]))} ) * sqrt(dotProduct(x,x,2).real());

    u[0] = x[0] + alpha;
    u[1] = x[1];
    qkxTMComplex beta = sqrt(dotProduct(u,u,2).real());
    v[0] = u[0]/beta;
    v[1] = u[1]/beta;    


    w = dotProduct(x,v,2) / dotProduct(v,x,2);

    for(int i = 0 ; i < 2  ; i++)
      for(int j = 0 ; j < 2  ; j++)
	Q[i*2+j] = (i==j ? (qkxTMComplex) {1.,0.} : (qkxTMComplex) {0.,0.}) - (1.+w)*v[i]*conj(v[j]);
      

    memset(R,0,2*sizeof(qkxTMComplex));
    
    for(int i = 0 ; i < 2 ; i++)
      for(int j = 0 ; j < 2 ; j++)
	R[i] = R[i] + conj(Q[j*2+i])*x[j];

}

void Arnoldi::QR_hs(Arnoldi &H, Arnoldi &Tmp){
  int N = arnoldiTotalVectors;
  memcpy(Tmp.P_matrixKxK(), H.P_matrixKxK(),N*N*sizeof(qkxTMComplex));
  qkxTMComplex *x = (qkxTMComplex*)calloc(N,sizeof(qkxTMComplex));
  qkxTMComplex *u = (qkxTMComplex*)calloc(N,sizeof(qkxTMComplex));
  qkxTMComplex *v = (qkxTMComplex*)calloc(N,sizeof(qkxTMComplex));

  qkxTMComplex *Qi = (qkxTMComplex*)calloc(N*N,sizeof(qkxTMComplex));
  qkxTMComplex *tmpOld = (qkxTMComplex*)calloc(N*N,sizeof(qkxTMComplex));

  qkxTMComplex *Q = this->P_matrixKxK();
  qkxTMComplex *tmp = Tmp.P_matrixKxK();
  
  for(int i = 0 ; i < N ; i++)
    for(int j = 0 ; j < N ; j++){
      Q[i*N+j] = i==j ? (qkxTMComplex) {1.,0.} : (qkxTMComplex) {0.,0.}; 
    }

  for(int k = 0 ; k < N-1 ; k++){
    for(int i = k ; i < N ; i++) x[i-k] = tmp[i*N+k];
    qkxTMComplex alpha = - exp( (qkxTMComplex) {0.,atan2(imag(x[0]),real(x[0]))} ) * sqrt(dotProduct(x,x,N-k).real());
    for(int i = k ; i < N ; i++) u[i-k] = x[i-k] + alpha*( (i-k) == 0 ? 1 : 0);

    for(int i = k ; i < N ; i++) v[i-k] = u[i-k] *(1./sqrt(dotProduct(u,u,N-k).real()));
    qkxTMComplex w = dotProduct(x,v,N-k) / dotProduct(v,x,N-k);

    for(int i = 0 ; i < N  ; i++)
      for(int j = 0 ; j < N  ; j++){
	if( i < k || j < k) {Qi[i*N+j] = i==j ? (qkxTMComplex) {1.,0.} : (qkxTMComplex) {0.,0.};}
	else{
	  Qi[i*N+j] = (i==j ? (qkxTMComplex) {1.,0.} : (qkxTMComplex) {0.,0.}) - (1.+w)*v[i-k]*conj(v[j-k]);
	}
      }


    memcpy(tmpOld,tmp,N*N*sizeof(qkxTMComplex));
    memset(tmp,0,N*N*sizeof(qkxTMComplex));
    
    for(int i = 0 ; i < N ; i++)
      for(int j = 0 ; j < N ; j++)
	for(int l = 0 ; l < N ; l++)
	  tmp[i*N+j] = tmp[i*N+j] + Qi[i*N+l]*tmpOld[l*N+j];
	  
    memset(tmpOld,0,N*N*sizeof(qkxTMComplex));
    for(int i = 0 ; i < N ; i++)
      for(int j = 0 ; j < N ; j++){
	for(int l = 0 ; l < N ; l++)
	  tmpOld[i*N+j] = tmpOld[i*N+j] + Q[i*N+l]*Qi[l*N+j];
      }
    
    memcpy(Q,tmpOld,N*N*sizeof(qkxTMComplex));
  }
  free(x);
  free(u);
  free(v);
  free(Qi);
  free(tmpOld);
}


static void givensRotation(qkxTMComplex x,qkxTMComplex y, qkxTMComplex *vec){
  qkxTMComplex *c = &vec[0];
  qkxTMComplex *s = &vec[1];
  qkxTMComplex *t = &vec[2];
  qkxTMComplex *r = &vec[3];
  qkxTMComplex u;

  if(x.real() != 0 || x.imag() != 0 ){
    *t = conj(y/x);
    u = sqrt(1 + conj(*t)*(*t));
    *r = u*x;
    *c = ((qkxTMComplex) {1.,0.} )/u ;
    *s = (*c)*(*t);
  }
  else{
    *t = 1/0.;
    *r = y;
    *c = (qkxTMComplex) {0.,0.};
    *s = (qkxTMComplex) {1.,0.};
  }
}


void Arnoldi::QR_gvns(Arnoldi &Q, Arnoldi &R){

  int N = arnoldiTotalVectors;
  qkxTMComplex *q = Q.P_matrixKxK();
  qkxTMComplex *r = R.P_matrixKxK();
  qkxTMComplex vec[4];
  qkxTMComplex *tmp = (qkxTMComplex*)malloc(N*N*sizeof(qkxTMComplex));
  memcpy(tmp, this->P_matrixKxK(),N*N*sizeof(qkxTMComplex));
  qkxTMComplex *tmpOld = (qkxTMComplex*)calloc(N*N,sizeof(qkxTMComplex));
  qkxTMComplex matrix[2][2];


  for(int j = 0 ; j < N-1 ; j++)
    for(int i = N-2 ; i >= j ; i--){

      givensRotation(tmp[i*N+j],tmp[(i+1)*N+j],vec);
      tmp[(i+1)*N+j] = vec[2];
      tmp[i*N+j] = vec[3];
      matrix[0][0] = vec[0] ; matrix[0][1] = vec[1] ; matrix[1][0] = -conj(vec[1]); ; matrix[1][1] = vec[0];
    
      if( j < N-1){
      memset(tmpOld,0,N*N*sizeof(qkxTMComplex));
      for(int ir = i ; ir <= i+1 ; ir++)
	for(int ic = j+1 ; ic < N ; ic++)
	  for(int ll = i ; ll <= i+1 ; ll++)
	    tmpOld[ir*N+ic] = tmpOld[ir*N+ic] + matrix[ir-i][ll-i] * tmp[ll*N+ic];
      
      
	  for(int ir = i ; ir <= i+1 ; ir++)
	    for(int ic = j+1 ; ic < N ; ic++)
	      tmp[ir*N+ic] = tmpOld[ir*N+ic];
      }
    }
  // second round
  qkxTMComplex t,c,s;
 
  for(int j = N-1 ; j>= 0 ; j--){
    for(int l = 0 ; l <= j ; l++) tmp[l*N+j] = (qkxTMComplex) {0.,0.};
    tmp[j*N+j] = (qkxTMComplex){1.,0.};

    for(int i = j ; i < N-1 ; i++){

      t = tmp[(i+1)*N+j];
      tmp[(i+1)*N+j] = (qkxTMComplex) {0.,0.,};
      c = ( (qkxTMComplex) {1.,0.} ) / sqrt(1. + conj(t)*t);
      if( c.real() != 0 || c.imag() != 0)
	s = c*t;
      else 
	s = (qkxTMComplex) {1.,0.};
      matrix[0][0] = c; matrix[0][1] = -s; matrix[1][0] = conj(s); matrix[1][1] = c;


      memset(tmpOld,0,N*N*sizeof(qkxTMComplex));
      for(int ir = i ; ir <= i+1 ; ir++)
	for(int ic = j ; ic < N ; ic++)
	  for(int ll = i ; ll <= i+1 ; ll++)
	    tmpOld[ir*N+ic] = tmpOld[ir*N+ic] + matrix[ir-i][ll-i] * tmp[ll*N+ic];
	  

	  for(int ir = i ; ir <= i+1 ; ir++)
	    for(int ic = j ; ic < N ; ic++)
	      tmp[ir*N+ic] = tmpOld[ir*N+ic];
	  
    }

  }


  memcpy(q,tmp,N*N*sizeof(qkxTMComplex));

  memset(r,0,N*N*sizeof(qkxTMComplex));
  for(int i = 0 ; i < N ; i++)
    for(int j = 0 ; j < N ; j++)
      for(int l = 0 ; l < N ; l++)
	r[i*N+j] = r[i*N+j] + conj(q[l*N+i])*this->P_matrixKxK()[l*N+j];
 

  free(tmpOld);
  free(tmp);
}

/*
void Arnoldi::QRLapack(Arnoldi &Q, Arnoldi &R){

  int M = arnoldiTotalVectors;
  qkxTMComplex *q = Q.P_matrixKxK();
  qkxTMComplex *r = R.P_matrixKxK();
 
  qkxTMComplex *mh_f = (qkxTMComplex*)malloc(M*M*sizeof(qkxTMComplex));
  for(int i = 0 ; i < M ; i++)
    for(int j = 0 ; j < M ; j++)
      mh_f[i*M+j] = this->P_matrixKxK()[j*M+i];  // change to fortran conventions

  qkxTMComplex wkopt;
  qkxTMComplex *work;
  int lwork=-1;
  int info;
  qkxTMComplex *tau = (qkxTMComplex*)malloc(M*sizeof(qkxTMComplex));

  zgeqrf_(&M,&M,(__complex__ double*)mh_f,&M,(__complex__ double*)tau,(__complex__ double*)&wkopt, &lwork,&info);
  lwork = (int)wkopt.real();
  work = (qkxTMComplex*)malloc( lwork*sizeof(qkxTMComplex) );
  zgeqrf_(&M,&M,(__complex__ double*)mh_f,&M,(__complex__ double*)tau,(__complex__ double*)work, &lwork,&info);


  qkxTMComplex *mh_c = (qkxTMComplex*)malloc(M*M*sizeof(qkxTMComplex));

  for(int i = 0 ; i < M ; i++)
    for(int j = 0 ; j < M ; j++)
      mh_c[j*M+i] = mh_f[i*M+j];
  free(mh_f);  

  qkxTMComplex *vv = (qkxTMComplex*)malloc(M*sizeof(qkxTMComplex));
  qkxTMComplex *HH = (qkxTMComplex*)malloc(M*M*sizeof(qkxTMComplex));
  qkxTMComplex *HHtmp = (qkxTMComplex*)calloc(M*M,sizeof(qkxTMComplex));

  
  for(int i = 0 ; i < M ; i++) HHtmp[i*M+i] = (qkxTMComplex) {1.,0.};

  for(int i = 0 ; i < M ; i++){
    memset(vv,0,M*sizeof(qkxTMComplex));
    vv[i] = (qkxTMComplex) {1.,0.};
    for(int j = i+1 ; j < M ; j++)vv[j] = mh_c[j*M+i];

    // construct H
    for(int ii = 0 ; ii < M ; ii++)
      for(int jj = 0 ; jj < M ; jj++)
	HH[ii*M+jj] = (qkxTMComplex) {ii==jj ? 1. : 0., 0.} - tau[i] * vv[ii] * conj(vv[jj]);

    matrixXmatrix(HHtmp,HH,q,M);
    memcpy(HHtmp,q,M*M*sizeof(qkxTMComplex));
  }
  

  memset(r,0,M*M*sizeof(qkxTMComplex));

  for(int j = 0 ; j < M ; j++)
    for(int i = 0 ; i <= j ; i++)
      r[i*M+j] = mh_c[i*M+j];

  free(HHtmp);
  free(HH);
  free(vv);
  free(mh_c);
  free(work);
  free(tau);
}
*/

void Arnoldi::QRLapack(Arnoldi &Q, Arnoldi &R){

  int M = arnoldiTotalVectors;
  qkxTMComplex *q = Q.P_matrixKxK();
  qkxTMComplex *r = R.P_matrixKxK();
 
  qkxTMComplex *mh_f = (qkxTMComplex*)malloc(M*M*sizeof(qkxTMComplex));
  for(int i = 0 ; i < M ; i++)
    for(int j = 0 ; j < M ; j++)
      mh_f[i*M+j] = this->P_matrixKxK()[j*M+i];  // change to fortran conventions

  qkxTMComplex wkopt;
  qkxTMComplex *work;
  int lwork=-1;
  int info;
  qkxTMComplex *tau = (qkxTMComplex*)malloc(M*sizeof(qkxTMComplex));

  zgeqrf_(&M,&M,(__complex__ double*)mh_f,&M,(__complex__ double*)tau,(__complex__ double*)&wkopt, &lwork,&info);
  lwork = (int)wkopt.real();
  work = (qkxTMComplex*)malloc( lwork*sizeof(qkxTMComplex) );
  zgeqrf_(&M,&M,(__complex__ double*)mh_f,&M,(__complex__ double*)tau,(__complex__ double*)work, &lwork,&info);

  if( info != 0) {
    fprintf(stderr,"Error in qr function\n");
    exit(-1);
  }

  qkxTMComplex *mh_c = (qkxTMComplex*)malloc(M*M*sizeof(qkxTMComplex));

  for(int i = 0 ; i < M ; i++)
    for(int j = 0 ; j < M ; j++)
      mh_c[j*M+i] = mh_f[i*M+j];

  memset(r,0,M*M*sizeof(qkxTMComplex));
  for(int j = 0 ; j < M ; j++)
    for(int i = 0 ; i <= j ; i++)
      r[i*M+j] = mh_c[i*M+j];


  zungqr_(&M,&M,&M,(__complex__ double*)mh_f,&M,(__complex__ double*)tau,(__complex__ double*)work, &lwork,&info);

  if( info != 0) {
    fprintf(stderr,"Error in qr function\n");
    exit(-1);
  }

  for(int i = 0 ; i < M ; i++)
    for(int j = 0 ; j < M ; j++)
      mh_c[j*M+i] = mh_f[i*M+j];


  for(int i = 0 ; i < M ; i++)
    for(int j = 0 ; j < M ; j++)
      q[i*M+j] = mh_c[i*M+j];

  free(mh_f);    
  free(mh_c);
  free(work);
  free(tau);
}


void Arnoldi::eig(Arnoldi &X,Arnoldi &W){

  int N = arnoldiTotalVectors;
  double precision = N*2.22e-16;
  double residual = 1;
  Arnoldi *q = new Arnoldi(arInfo,latInfo,matrixKxK);
  Arnoldi *r = new Arnoldi(arInfo,latInfo,matrixKxK);
  Arnoldi *Qtmp = new Arnoldi(arInfo,latInfo,matrixKxK);
  Arnoldi *Q = new Arnoldi(arInfo,latInfo,matrixKxK);
  Arnoldi *htmp = new Arnoldi(arInfo,latInfo,matrixKxK);
  
  qkxTMComplex *C1 = (qkxTMComplex*)malloc(N*N*sizeof(qkxTMComplex)); 
  qkxTMComplex *C2 = (qkxTMComplex*)malloc(N*N*sizeof(qkxTMComplex));
  qkxTMComplex *D = (qkxTMComplex*)malloc(N*N*sizeof(qkxTMComplex));
  qkxTMComplex *w = (qkxTMComplex*)calloc(N,sizeof(qkxTMComplex));

  memset(W.P_matrixKxK(),0,N*N*sizeof(qkxTMComplex));
  memcpy(htmp->P_matrixKxK(),this->P_matrixKxK(),N*N*sizeof(qkxTMComplex));
  for(int i = 0 ; i < N ; i++)
    for(int j = 0 ; j < N ; j++)
      Qtmp->P_matrixKxK()[i*N+j] = (qkxTMComplex) { i==j?1.:0.,0.};

  while( residual > precision ){
    htmp->QR_gs_m(*q,*r);
    matrixXmatrix(r->P_matrixKxK(),q->P_matrixKxK(),htmp->P_matrixKxK(),N);
    matrixXmatrix(Qtmp->P_matrixKxK(),q->P_matrixKxK(),Q->P_matrixKxK(),N);
    memcpy(Qtmp->P_matrixKxK(),Q->P_matrixKxK(),N*N*sizeof(qkxTMComplex));    
    for(int i=0;i<N;i++) W.P_matrixKxK()[i*N+i] = htmp->P_matrixKxK()[i*N+i];
    residual=0.;
    for(int i = 0 ; i< N ;i++)residual+=sqrt(norm(W.P_matrixKxK()[i*N+i] - w[i]));
    for(int i=0;i<N;i++) w[i] = W.P_matrixKxK()[i*N+i];
    printf("%e\n",residual);
  }

  // now htmp has the triangular matrix we need for the eigenvectors

  memset(X.P_matrixKxK(),0,N*N*sizeof(qkxTMComplex));

  for(int i = 0 ; i < N ; i++) X.P_matrixKxK()[i*N+i] = (qkxTMComplex) {1.,0.};
  for(int i =N-1; i >=0 ; i--)
    for(int j = i-1; j>= 0 ; j--){
      qkxTMComplex sum = (qkxTMComplex) {0.,0.};
      for(int k = j+1 ; k <= i ; k++) sum = sum +htmp->P_matrixKxK()[j*N+k]*X.P_matrixKxK()[k*N+i];
      X.P_matrixKxK()[j*N+i] = -sum/(W.P_matrixKxK()[j*N+j] - W.P_matrixKxK()[i*N+i]);
  }

     
  for(int i = 0 ; i < N ; i++){
    for(int j = 0 ; j < N ; j++) w[j] = X.P_matrixKxK()[j*N+i];
    double norma = sqrt(dotProduct(w,w,N).real());
    for(int j = 0 ; j < N ; j++) X.P_matrixKxK()[j*N+i] = X.P_matrixKxK()[j*N+i]*(1./norma);
  }

  
  memcpy(htmp->P_matrixKxK(),X.P_matrixKxK(),N*N*sizeof(qkxTMComplex));
  matrixXmatrix(Q->P_matrixKxK(),htmp->P_matrixKxK(),X.P_matrixKxK(),N);
  
  free(C1);
  free(C2);
  free(D);
  free(w);
  delete Q;
  delete q;
  delete r;
  delete Qtmp;
  delete htmp;

}

static void noiseCleaner(qkxTMComplex* A, int M, int N){
  double machine_epsilon = std::numeric_limits<double>::epsilon();
  double matrixOrder = 0.;

  for(int i = 0 ; i < M ; i++)
    for(int j = 0 ; j < N ; j++)
      matrixOrder += sqrt(norm(A[i*M+j]));
  for(int i = 0 ; i < M ; i++)
    for(int j = 0 ; j < N ; j++){
      if( sqrt(norm(A[i*M+j])) < sqrt(machine_epsilon)/matrixOrder )
	A[i*M+j] = (qkxTMComplex) {0.,0.};
    }


  for(int i = 0 ; i < M ; i++)
    for(int j = 0 ; j < N ; j++){
      if( fabs(A[i*M+j].real()) < sqrt(machine_epsilon)/matrixOrder) A[i*M+j].real() = 0;
      if( fabs(A[i*M+j].imag()) < sqrt(machine_epsilon)/matrixOrder) A[i*M+j].imag() = 0;
    }


}


void Arnoldi::eigLapack(Arnoldi &Q,Arnoldi &W){

  int M = arnoldiTotalVectors;
  qkxTMComplex *w = (qkxTMComplex*)calloc(M,sizeof(qkxTMComplex));
  qkxTMComplex *lq = (qkxTMComplex*)calloc(M*M,sizeof(qkxTMComplex));
  qkxTMComplex *rq = (qkxTMComplex*)calloc(M*M,sizeof(qkxTMComplex));
  qkxTMComplex *h = (qkxTMComplex*)calloc(M*M,sizeof(qkxTMComplex));
  
  for(int i = 0 ; i < M ; i++)
    for(int j = 0 ; j < M ; j++)
      h[i*M+j] = this->P_matrixKxK()[j*M+i];  // change to fortran conventions

  noiseCleaner(h,M,M);

  qkxTMComplex wkopt;
  qkxTMComplex *work;
  double *rwork = (double*)calloc(2*M,sizeof(double));
  int lwork = -1;
  int info;
  int n = M, lda = M, ldvl = M, ldvr = M ;


  zgeev_((char*) "Vectors",(char*)"Vectors", &n,(double __complex__*) h, &lda,(double __complex__*) w,(double __complex__*) lq, &ldvl,(double __complex__*) rq, &ldvr,(double __complex__*) &wkopt, &lwork, rwork, &info );  

  lwork = (int)wkopt.real();
  work = (qkxTMComplex*)malloc( lwork*sizeof(qkxTMComplex) );

  zgeev_((char*)"Vectors",(char*)"Vectors", &n,(double __complex__*) h, &lda,(double __complex__*) w,(double __complex__*) lq, &ldvl,(double __complex__*) rq, &ldvr,(double __complex__*) work, &lwork, rwork, &info );  


  if( info > 0 ){
    fprintf(stderr,"Error eigLapack failed to converge\n");
    exit(-1);
  }
  
  memset(W.P_matrixKxK(),0,M*M*sizeof(qkxTMComplex));
  for(int i = 0 ; i < M ; i++) W.P_matrixKxK()[i*M+i] = w[i];

  for(int i = 0 ; i < M ; i++)
    for(int j = 0 ; j < M ; j++)
      Q.P_matrixKxK()[j*M+i] = rq[i*M+j];


  free(work);
  free(rwork);
  free(w);
  free(lq);
  free(rq);
  free(h);
}


void qkxTM::initArnold(Arnoldi &A, Arnoldi &V, Arnoldi &H){
  int NL = A.DimensionMatrix();
  int M = A.ArnoldiTotalVectors();
  int MP1 = A.ArnoldiTotalVectors() +1;



  qkxTMComplex *v = V.P_matrixNxK();
  qkxTMComplex *h = H.P_matrixKxK();
  qkxTMComplex *a = A.P_matrixNxN();
  qkxTMComplex *vectmp = (qkxTMComplex*)malloc(NL*sizeof(qkxTMComplex));
  qkxTMComplex *vin = (qkxTMComplex*)malloc(NL*sizeof(qkxTMComplex));

  for(int i = 0 ; i < NL ; i++) v[i*MP1+0] =(qkxTMComplex) {drand48()-0.5,0.};
  for(int i = 0 ; i < NL ; i++) vectmp[i] = v[i*MP1+0];
  qkxTMComplex beta = sqrt(dotProduct(vectmp,vectmp,NL).real());
  for(int i = 0 ; i < NL ; i++) v[i*MP1+0] = v[i*MP1+0]/beta;

  for(int i = 0 ; i < NL ; i++) vin[i] = v[i*MP1+0];
  matrixXvec(a,vin,vectmp,NL);

  for(int i = 0 ; i < NL ; i++) v[i*MP1+1] = vectmp[i];
  qkxTMComplex alpha;
  alpha = dotProduct(vin,vectmp,NL);
  h[0*M+0] = alpha;

  for(int i = 0 ; i < NL ; i++) v[i*MP1+1] = v[i*MP1+1] - alpha*v[i*MP1+0];

  // check orthogonality
  for(int i = 0 ; i < NL ; i++) vectmp[i] = v[i*MP1+1];
  alpha = dotProduct(vin,vectmp,NL);
  h[0*M+0] = h[0*M+0] + alpha;
  for(int i = 0 ; i < NL ; i++) v[i*MP1+1] = v[i*MP1+1] - alpha*v[i*MP1+0];
  

  // initialize vector
  free(vectmp);
  free(vin);
}


void qkxTM::initArnold(Arnoldi &V, Arnoldi &H,GaugeField &gauge){
  int NL = V.DimensionMatrix();
  int M = V.ArnoldiTotalVectors();
  int MP1 = V.ArnoldiTotalVectors() +1;
  int vol = V.Volume();

  LatticeInfo latInfo = V.LatInfo();
  qkxTMComplex *v = V.P_matrixNxK();
  qkxTMComplex *h = H.P_matrixKxK();
 
  qkxTMComplex *vectmp = (qkxTMComplex*)malloc(NL*sizeof(qkxTMComplex));
  qkxTMComplex *vin = (qkxTMComplex*)malloc(NL*sizeof(qkxTMComplex));

  ColorSpinorField *vecN1 = new ColorSpinorField(&latInfo);
  ColorSpinorField *vecN2 = new ColorSpinorField(&latInfo);

  vecN1->create();
  vecN2->create();


  //  for(int i = 0 ; i < NL ; i++) v[i*MP1+0] =(qkxTMComplex) {drand48()-0.5,0.};

  for(int i = 0 ; i < NL ; i++) v[i*MP1+0] = (qkxTMComplex) {1.,0.};

  for(int i = 0 ; i < NL ; i++) vectmp[i] = v[i*MP1+0];
 
 qkxTMComplex beta = sqrt(dotProduct(vectmp,vectmp,NL).real());

  for(int i = 0 ; i < NL ; i++) v[i*MP1+0] = v[i*MP1+0]/beta;


 
  for(int i = 0 ; i < NL ; i++) vin[i] = v[i*MP1+0];

  for(int mu = 0 ; mu < 4 ; mu++)
    for(int c = 0 ; c < 3 ; c++)
      for(int iv = 0 ; iv < vol ; iv++)
	vecN1->P_colorSpinor()[MS(mu,c)][iv] = v[(mu*3*vol + c*vol + iv)*MP1+0];

  vecN2->MdagM(vecN1->P_colorSpinor(),gauge.P_gauge(),&latInfo);
    
  for(int mu = 0 ; mu < 4 ; mu++)
    for(int c = 0 ; c < 3 ; c++)
      for(int iv = 0 ; iv < vol ; iv++)
	vectmp[mu*3*vol + c*vol + iv] = vecN2->P_colorSpinor()[MS(mu,c)][iv]; 

  for(int i = 0 ; i < NL ; i++) v[i*MP1+1] = vectmp[i];
  qkxTMComplex alpha;
  alpha = dotProduct(vin,vectmp,NL);

  
 h[0*M+0] = alpha;

  for(int i = 0 ; i < NL ; i++) v[i*MP1+1] = v[i*MP1+1] - alpha*v[i*MP1+0];

  // check orthogonality
  for(int i = 0 ; i < NL ; i++) vectmp[i] = v[i*MP1+1];
  alpha = dotProduct(vin,vectmp,NL);
  h[0*M+0] = h[0*M+0] + alpha;
  for(int i = 0 ; i < NL ; i++) v[i*MP1+1] = v[i*MP1+1] - alpha*v[i*MP1+0];
   

  delete vecN1;
  delete vecN2;

  // initialize vector
  free(vectmp);
  free(vin);
}



void qkxTM::arnold(int kstart,int kstop,Arnoldi &V, Arnoldi &H,GaugeField &gauge){
  // V is nx(m+1), H is mxm
  int NL = V.DimensionMatrix();
  int M = V.ArnoldiTotalVectors();
  int MP1 = V.ArnoldiTotalVectors() +1;

  int vol = V.Volume();
  LatticeInfo latInfo = V.LatInfo();


  qkxTMComplex *v = V.P_matrixNxK();
  qkxTMComplex *h = H.P_matrixKxK();

  ColorSpinorField *vecN1 = new ColorSpinorField(&latInfo);
  ColorSpinorField *vecN2 = new ColorSpinorField(&latInfo);

  vecN1->create();
  vecN2->create();


  qkxTMComplex *vectmp = (qkxTMComplex*)malloc(NL*sizeof(qkxTMComplex));
  qkxTMComplex *vin = (qkxTMComplex*)malloc(NL*sizeof(qkxTMComplex));
  qkxTMComplex *s = (qkxTMComplex*)malloc(M*M*sizeof(qkxTMComplex));

  for(int j = kstart; j < kstop ; j++){
    int jm1 = j-1;
    int jp1 = j+1;
    for(int i = 0 ; i < NL ; i++) vectmp[i] = v[i*MP1+j];
    qkxTMComplex beta = sqrt(dotProduct(vectmp,vectmp,NL).real());
    h[j*M+jm1] = beta;

    for(int i = 0 ; i < NL ; i++) v[i*MP1+j] = v[i*MP1+j]/beta;


  for(int mu = 0 ; mu < 4 ; mu++)
    for(int c = 0 ; c < 3 ; c++)
      for(int iv = 0 ; iv < vol ; iv++)
	vecN1->P_colorSpinor()[MS(mu,c)][iv] = v[(mu*3*vol + c*vol + iv)*MP1+j];
  
  vecN2->MdagM(vecN1->P_colorSpinor(),gauge.P_gauge(),&latInfo);

  for(int mu = 0 ; mu < 4 ; mu++)
    for(int c = 0 ; c < 3 ; c++)
      for(int iv = 0 ; iv < vol ; iv++)
	v[(mu*3*vol + c*vol + iv)*MP1+jp1] = vecN2->P_colorSpinor()[MS(mu,c)][iv];



    for(int l = 0 ; l < NL ; l++) vectmp[l] = v[l*MP1+jp1];
    for(int i = 0 ; i <= j ; i++){
      for(int l = 0 ; l < NL ; l++) vin[l] = v[l*MP1+i];
      h[i*M+j] = dotProduct(vin,vectmp,NL);
    }

    //    printf("%+e %+e\n",v[0*MP1+2].real(),v[0*MP1+2].imag());

    for(int i = 0 ; i < NL ; i++){
      qkxTMComplex sum = (qkxTMComplex) {0.,0.};
      for(int l = 0 ; l <= j ; l++)
	sum = sum + v[i*MP1+l]*h[l*M+j];
      v[i*MP1+jp1] = v[i*MP1+jp1] - sum;
    }  

    
    for(int l = 0 ; l < NL ; l++) vectmp[l] = v[l*MP1+jp1];
    for(int i = 0 ; i <= j ; i++){
      for(int l = 0 ; l < NL ; l++) vin[l] = v[l*MP1+i];
      s[i*M+j] = dotProduct(vin,vectmp,NL);
    }

    for(int i = 0 ; i < NL ; i++){
      qkxTMComplex sum = (qkxTMComplex) {0.,0.};
      for(int l = 0 ; l <= j ; l++)
	sum = sum + v[i*MP1+l]*s[l*M+j];
      v[i*MP1+jp1] = v[i*MP1+jp1] - sum;
    }
    
    for(int i = 0 ; i <= j ; i++)
      h[i*M+j] = h[i*M+j] + s[i*M+j];
  }

  
  
  //  exit(-1);
  delete vecN1;
  delete vecN2;

  free(vectmp);
  free(vin);
  free(s);
}

void qkxTM::arnold(int kstart,int kstop,Arnoldi &A,Arnoldi &V, Arnoldi &H){
  // V is nx(m+1), H is mxm
  int NL = A.DimensionMatrix();
  int M = A.ArnoldiTotalVectors();
  int MP1 = A.ArnoldiTotalVectors() +1;

  qkxTMComplex *v = V.P_matrixNxK();
  qkxTMComplex *h = H.P_matrixKxK();
  qkxTMComplex *a = A.P_matrixNxN();

  qkxTMComplex *vectmp = (qkxTMComplex*)malloc(NL*sizeof(qkxTMComplex));
  qkxTMComplex *vin = (qkxTMComplex*)malloc(NL*sizeof(qkxTMComplex));
  qkxTMComplex *s = (qkxTMComplex*)malloc(M*M*sizeof(qkxTMComplex));

  for(int j = kstart; j < kstop ; j++){
    int jm1 = j-1;
    int jp1 = j+1;
    for(int i = 0 ; i < NL ; i++) vectmp[i] = v[i*MP1+j];
    qkxTMComplex beta = sqrt(dotProduct(vectmp,vectmp,NL).real());
    h[j*M+jm1] = beta;

    for(int i = 0 ; i < NL ; i++) v[i*MP1+j] = v[i*MP1+j]/beta;
    for(int i = 0 ; i < NL ; i++) vin[i] = v[i*MP1+j];
    matrixXvec(a,vin,vectmp,NL);
    for(int i = 0 ; i < NL ; i++) v[i*MP1+jp1] = vectmp[i];

    for(int l = 0 ; l < NL ; l++) vectmp[l] = v[l*MP1+jp1];
    for(int i = 0 ; i <= j ; i++){
      for(int l = 0 ; l < NL ; l++) vin[l] = v[l*MP1+i];
      h[i*M+j] = dotProduct(vin,vectmp,NL);
    }

    for(int i = 0 ; i < NL ; i++){
      qkxTMComplex sum = (qkxTMComplex) {0.,0.};
      for(int l = 0 ; l <= j ; l++)
	sum = sum + v[i*MP1+l]*h[l*M+j];
      v[i*MP1+jp1] = v[i*MP1+jp1] - sum;
    }

    
    for(int l = 0 ; l < NL ; l++) vectmp[l] = v[l*MP1+jp1];
    for(int i = 0 ; i <= j ; i++){
      for(int l = 0 ; l < NL ; l++) vin[l] = v[l*MP1+i];
      s[i*M+j] = dotProduct(vin,vectmp,NL);
    }

    for(int i = 0 ; i < NL ; i++){
      qkxTMComplex sum = (qkxTMComplex) {0.,0.};
      for(int l = 0 ; l <= j ; l++)
	sum = sum + v[i*MP1+l]*s[l*M+j];
      v[i*MP1+jp1] = v[i*MP1+jp1] - sum;
    }
    
    for(int i = 0 ; i <= j ; i++)
      h[i*M+j] = h[i*M+j] + s[i*M+j];
  }


    
  free(vectmp);
  free(vin);
  free(s);
}


void qkxTM::qrShiftsRotations(Arnoldi &V, Arnoldi &H, Arnoldi &W){
  int M = V.ArnoldiTotalVectors();
  int MP1 = M+1;
  int K = V.ArnoldiWantedVectors();
  int P = V.ArnoldiUnwantedVectors();


  qkxTMComplex x[2];
  qkxTMComplex G[4];
  qkxTMComplex R[2];
  qkxTMComplex *tmp = (qkxTMComplex*)calloc(2*M,sizeof(qkxTMComplex));
  qkxTMComplex tmp2[2];
  qkxTMComplex *e = (qkxTMComplex*)calloc(M,sizeof(qkxTMComplex));
  e[M-1] = (qkxTMComplex) {1.,0.};

  int NL = V.DimensionMatrix();


  qkxTMComplex *h = H.P_matrixKxK();
  qkxTMComplex *v = V.P_matrixNxK();
  qkxTMComplex *w = W.P_matrixKxK();


  for(int jj = 0 ; jj < P ; jj++){
    x[0] = h[0*M+0] - w[jj*M+jj];
    x[1] = h[1*M+0];
    QR_hs_v2(G,R,x);    
    for(int ii = 0 ; ii < M-1 ; ii++){
      if(ii > 0){
	x[0] = h[ii*M+ii-1];
	x[1] = h[(ii+1)*M+ii-1];
	QR_hs_v2(G,R,x);
	h[ii*M+ii-1] = R[0];
	h[(ii+1)*M+ii-1] = R[1];
      }

      // transform H
      memset(tmp,0,2*M*sizeof(qkxTMComplex));
      for(int j = ii ; j < M ; j++)
	for(int i = 0 ; i < 2 ; i++)
	  for(int l = 0 ; l < 2 ; l++)
	    tmp[i*M+j] = tmp[i*M+j] + conj(G[l*2+i])*h[(ii+l)*M+j];
      for(int i = 0 ; i < 2 ; i++)
	for(int j = ii ; j < M ; j++)
	  h[(ii+i)*M+j] = tmp[i*M+j];


      memset(tmp,0,2*M*sizeof(qkxTMComplex));
      for(int i = 0 ; i < M ; i++)
	for(int j = 0 ; j < 2 ; j++)
	  for(int l = 0 ; l < 2 ; l++)
	    tmp[i*2+j] = tmp[i*2+j] + h[i*M+ii+l]*G[l*2+j];
      for(int i = 0 ; i < M ; i++)
	for(int j = 0 ; j < 2 ; j++)
	  h[i*M+ii+j] = tmp[i*2+j]; 
      ////////////////////////////////
      //transform V
      for(int i = 0 ; i < NL ; i++){
	memset(tmp2,0,2*sizeof(qkxTMComplex));
	for(int j = 0 ; j < 2 ; j++)
	  for(int l = 0 ; l < 2 ; l++)
	    tmp2[j] = tmp2[j] + v[i*MP1+ii+l]*G[l*2+j];
	for(int j = 0 ; j < 2 ; j++)
	  v[i*MP1+ii+j] = tmp2[j];
      }
      /////////////////
      // transform e
      memset(tmp2,0,2*sizeof(qkxTMComplex));
      for(int i = 0 ; i < 2 ; i++)
	for(int j = 0 ; j < 2 ; j++)
	  tmp2[i] = tmp2[i] + e[ii+j]*G[j*2+i];
      for(int i = 0 ; i < 2 ; i++) e[ii+i] = tmp2[i];

    }
  }

  // reculculate residual
  for(int i = 0 ;i < NL ; i++)
    v[i*MP1+M] = v[i*MP1+M]*e[K-1];
  for(int i = 0 ; i < NL ; i++)
    v[i*MP1+K] = v[i*MP1+M] + v[i*MP1+K]*h[K*M+K-1];

  free(tmp);
  free(e);
}



void qkxTM::qrShiftsRotations_v2(Arnoldi &V, Arnoldi &H, Arnoldi &W){
  int M = V.ArnoldiTotalVectors();
  int MP1 = M+1;
  int K = V.ArnoldiWantedVectors();
  int P = V.ArnoldiUnwantedVectors();
  int NL = V.DimensionMatrix();
  qkxTMComplex *e = (qkxTMComplex*)calloc(M,sizeof(qkxTMComplex));
  e[M-1] = (qkxTMComplex) {1.,0.};
  Arnoldi *QR = new Arnoldi(V.ArInfo(),V.LatInfo(),matrixKxK);
  Arnoldi *RR = new Arnoldi(V.ArInfo(),V.LatInfo(),matrixKxK);
  Arnoldi *Hshifted = new Arnoldi(V.ArInfo(),V.LatInfo(),matrixKxK);
  
  qkxTMComplex *C1 = (qkxTMComplex*)calloc(M*M,sizeof(qkxTMComplex));
  qkxTMComplex *vec1 = (qkxTMComplex*)calloc(M,sizeof(qkxTMComplex));
  qkxTMComplex *vec2 = (qkxTMComplex*)calloc(M,sizeof(qkxTMComplex));



  for(int jj = 0 ; jj < P ; jj++){


    
    memcpy(Hshifted->P_matrixKxK(),H.P_matrixKxK(),M*M*sizeof(qkxTMComplex));
    for(int i = 0 ; i < M ; i++) Hshifted->P_matrixKxK()[i*M+i] = Hshifted->P_matrixKxK()[i*M+i] - W.P_matrixKxK()[jj*M+jj];

    noiseCleaner(Hshifted->P_matrixKxK(),M,M);


    Hshifted->QRLapack(*QR,*RR);


    noiseCleaner(H.P_matrixKxK(),M,M);
    matrixDagXmatrix(QR->P_matrixKxK(),H.P_matrixKxK(),C1,M);
    matrixXmatrix(C1,QR->P_matrixKxK(),H.P_matrixKxK(),M);
    noiseCleaner(H.P_matrixKxK(),M,M);



    
    for(int i = 0 ; i < NL ; i++){
      for(int ii = 0 ; ii < M ; ii++) vec1[ii] = V.P_matrixNxK()[i*MP1+ii];
      vecXmatrix(vec1,QR->P_matrixKxK(),vec2,M);
      for(int ii = 0 ; ii < M ; ii++) V.P_matrixNxK()[i*MP1+ii] = vec2[ii];
    }
    

    for(int ii = 0 ; ii < M ; ii++) vec1[ii] = e[ii];
    vecXmatrix(vec1,QR->P_matrixKxK(),vec2,M);
    for(int ii = 0 ; ii < M ; ii++) e[ii] = vec2[ii];

  }

  
  for(int i = 0 ; i < NL ; i++) V.P_matrixNxK()[i*MP1 + M] = V.P_matrixNxK()[i*MP1 + M]*e[K-1];
  for(int i = 0 ; i < NL ; i++) V.P_matrixNxK()[i*MP1 + K] = V.P_matrixNxK()[i*MP1 + M] + V.P_matrixNxK()[i*MP1 + K] * H.P_matrixKxK()[K*M+K-1];

  delete QR;
  delete RR;
  delete Hshifted;
  free(C1);
  free(e);
  free(vec1);
  free(vec2);
}



void qkxTM::sortEigenValues(Arnoldi &Q, Arnoldi &W){
  int M = W.ArnoldiTotalVectors();

  qkxTMComplex *w = (qkxTMComplex*)malloc(M*sizeof(qkxTMComplex));
  double *w_mag = (double*)malloc(M*sizeof(double));
  double *w_mag_tmp = (double*)malloc(M*sizeof(double));

  int *pos = (int*)malloc(M*sizeof(int));
  qkxTMComplex *q = (qkxTMComplex*)malloc(M*M*sizeof(qkxTMComplex));
  memcpy(q,Q.P_matrixKxK(),M*M*sizeof(qkxTMComplex));

  for(int i = 0 ; i < M ; i++) w[i] = W.P_matrixKxK()[i*M+i];
  for(int i = 0 ; i < M ; i++) w_mag[i] = sqrt(norm(w[i]));
  memcpy(w_mag_tmp,w_mag,M*sizeof(double));

  for(int i = 0 ; i < M-1 ; i++)
    for(int j = i + 1 ; j < M ; j++){
      if(w_mag[i] < w_mag[j]){
	double a = w_mag[i];
	w_mag[i] = w_mag[j];
	w_mag[j] = a;
      }
    }



  for(int i = 0 ; i < M ; i++)
    for(int j = 0 ; j < M ; j++)
      if(w_mag[i] == w_mag_tmp[j])
	pos[i] = j;
  


  for(int i = 0 ; i < M ; i++) W.P_matrixKxK()[i*M+i] = w[pos[i]];

  for(int i = 0 ; i < M ; i++)
    for(int j = 0 ; j < M ; j++)
      Q.P_matrixKxK()[i*M+j] = q[i*M+pos[j]];

  free(w);
  free(w_mag);
  free(w_mag_tmp);
  free(pos);
  free(q);
}

void Arnoldi::Iram(Arnoldi &X, Arnoldi &H , Arnoldi &W){

  // this function works for Hermitian matrices and it calculates the smallest eigenvalues
  int NL = DimensionMatrix();
  int M = ArnoldiTotalVectors();
  int MP1 = M+1;
  int K = ArnoldiWantedVectors();


  int iterMax = MaxArnoldiIter();
  double tol = Tolerance();
  int iter;
  double res;
  qkxTMComplex *vectmp = (qkxTMComplex*)malloc(NL*sizeof(qkxTMComplex));
  qkxTMComplex *vi1 = (qkxTMComplex*)malloc(NL*sizeof(qkxTMComplex));
  qkxTMComplex *vi2 = (qkxTMComplex*)malloc(NL*sizeof(qkxTMComplex));

  qkxTMComplex *ww = (qkxTMComplex*)malloc(M*sizeof(qkxTMComplex));
  Arnoldi *Qeig = new Arnoldi(arInfo,latInfo,matrixKxK);
  Arnoldi *V = new Arnoldi(arInfo,latInfo,matrixNxK);
  seed48(&seed); // initialize the random number generator

  // initialize arnoldi
  initArnold(*this,*V,H);
  arnold(1,K,*this,*V,H);
  iter = 0; res = 1.;
  while ( (res > tol) && (iter<iterMax) ){
    arnold(K,M,*this,*V,H);
    H.eigLapack(*Qeig,W);
    sortEigenValues(*Qeig,W);
    qrShiftsRotations_v2(*V, H, W);
    H.eigLapack(*Qeig,W);
    sortEigenValues(*Qeig,W);
    memset(X.P_matrixNxK(),0,NL*MP1*sizeof(qkxTMComplex));
    for(int iN = 0 ; iN < NL ; iN++)
      for(int iK = 0 ; iK < K ; iK++)
	for(int iM = 0 ; iM < M ; iM++)
	  X.P_matrixNxK()[iN*MP1+iK] = X.P_matrixNxK()[iN*MP1+iK] + V->P_matrixNxK()[iN*MP1+iM] * Qeig->P_matrixKxK()[iM*M + M-iK-1];
    
    for(int k = 0 ; k < K ; k++){
      for(int l = 0 ; l < NL ; l++) vectmp[l] = X.P_matrixNxK()[l*MP1+k];
      double norma = sqrt(dotProduct(vectmp,vectmp,NL).real());
      for(int l = 0 ; l < NL ; l++)
	X.P_matrixNxK()[l*MP1+k] = X.P_matrixNxK()[l*MP1+k] *(1./norma);
    }
    
    
    for(int i = 0 ; i < M ; i++) ww[i] = W.P_matrixKxK()[i*M+i];
    for(int k = 0 ; k < K ; k++) W.P_matrixKxK()[k*M+k] = ww[M-k-1];

    for(int l = 0 ; l < NL ; l++) vectmp[l] = V->P_matrixNxK()[l*MP1+M];
    res = sqrt(dotProduct(vectmp,vectmp,NL).real()); 
    printf("iter = %d, res = %e\n",iter,res);
    iter++;
  }

  free(ww);
  free(vectmp);
  free(vi1);
  free(vi2);
  
  delete Qeig;
  delete V;

}

/*
static void test(Arnoldi &V, Arnoldi &H,GaugeField &gauge){
  int NL = V.DimensionMatrix();
  int M = V.ArnoldiTotalVectors();
  int MP1 = V.ArnoldiTotalVectors() +1;
  int vol = V.Volume();

  LatticeInfo latInfo = V.LatInfo();
  qkxTMComplex *v = V.P_matrixNxK();
  qkxTMComplex *h = H.P_matrixKxK();
 
  qkxTMComplex *vectmp = (qkxTMComplex*)malloc(NL*sizeof(qkxTMComplex));
  qkxTMComplex *vin = (qkxTMComplex*)malloc(NL*sizeof(qkxTMComplex));

  ColorSpinorField *vecN1 = new ColorSpinorField(&latInfo);
  ColorSpinorField *vecN2 = new ColorSpinorField(&latInfo);

  vecN1->create();
  vecN2->create();


  //  for(int i = 0 ; i < NL ; i++) v[i*MP1+0] =(qkxTMComplex) {drand48()-0.5,0.};

  for(int i = 0 ; i < NL ; i++) v[i*MP1+0] = (qkxTMComplex) {1.,0.} ;

  for(int i = 0 ; i < 100 ; i++){
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int c = 0 ; c < 3 ; c++)
	for(int iv = 0 ; iv < vol ; iv++)
	  vecN1->P_colorSpinor()[MS(mu,c)][iv] = v[(mu*3*vol + c*vol + iv)*MP1+0];

    vecN2->MdagM(vecN1->P_colorSpinor(),gauge.P_gauge(),&latInfo);
    
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int c = 0 ; c < 3 ; c++)
	for(int iv = 0 ; iv < vol ; iv++)
	  v[(mu*3*vol + c*vol + iv)*MP1+0] = vecN2->P_colorSpinor()[MS(mu,c)][iv]; 
  }

  FILE *ptr;
  ptr = fopen("kale.dat","w");
  for(int i = 0 ; i < NL ; i++)
    fprintf(ptr,"%+16.15e %+16.15e\n",v[i*MP1+0].real(),v[i*MP1+0].imag());
   


  delete vecN1;
  delete vecN2;

  // initialize vector
  free(vectmp);
  free(vin);
}
*/

void Arnoldi::Iram(Arnoldi &X, Arnoldi &W, GaugeField &gauge){
  std::clock_t    start;
  std::clock_t    stop;


  start = std::clock();
  // this function works for Hermitian matrices and it calculates the smallest eigenvalues
  int NL = DimensionMatrix();
  int M = ArnoldiTotalVectors();
  int MP1 = M+1;
  int K = ArnoldiWantedVectors();
  //  int P = ArnoldiUnwantedVectors();

  int iterMax = MaxArnoldiIter();
  double tol = Tolerance();
  int iter;
  double res;
  qkxTMComplex *vectmp = (qkxTMComplex*)malloc(NL*sizeof(qkxTMComplex));
  qkxTMComplex *vi1 = (qkxTMComplex*)malloc(NL*sizeof(qkxTMComplex));
  qkxTMComplex *vi2 = (qkxTMComplex*)malloc(NL*sizeof(qkxTMComplex));

  qkxTMComplex *ww = (qkxTMComplex*)malloc(M*sizeof(qkxTMComplex));
  Arnoldi *Qeig = new Arnoldi(arInfo,latInfo,matrixKxK);
  Arnoldi *V = new Arnoldi(arInfo,latInfo,matrixNxK);
  Arnoldi *H = new Arnoldi(arInfo,latInfo,matrixKxK);

  //test(*V,*H,gauge);  

  seed48(&seed); // initialize the random number generator
  /*
  FILE *ptr;
  char filename[257];
  sprintf(filename,"arnoldiRun_K%d_P%d.dat",K,M-K);
  ptr = fopen(filename,"w");
  */
  // initialize arnoldi
  /*
  initArnold(*V,*H,gauge);
  arnold(1,K,*V,*H,gauge);
  arnold(K,M,*V,*H,gauge);
  H->eigLapack(*Qeig,W);
  sortEigenValues(*Qeig,W);
  qrShiftsRotations_v2(*V, *H, W);    

  
  for(int i = 0 ; i < M ; i++){
    for(int j = 0 ; j < M ; j++)
      printf("%+e %+e\t",H->P_matrixKxK()[i*M+j].real(),H->P_matrixKxK()[i*M+j].imag());
    printf("\n");
  }
  printf("\n");
  for(int i = 0 ; i < M ; i++){
    for(int j = 0 ; j < M ; j++)
      printf("%+e %+e\t",W.P_matrixKxK()[i*M+j].real(),W.P_matrixKxK()[i*M+j].imag());
    printf("\n");
  }
  printf("\n");
  for(int i = 0 ; i < M ; i++){
    for(int j = 0 ; j < M ; j++)
      printf("%+e %+e\t",Qeig->P_matrixKxK()[i*M+j].real(),Qeig->P_matrixKxK()[i*M+j].imag());
    printf("\n");
  }
  


   FILE *ptr;
   ptr = fopen("kale.dat","w");
   for(int i = 0 ; i < M+1 ; i++){
     for(int j = 0 ; j < NL ; j++)
       fprintf(ptr,"%e %e\n",V->P_matrixNxK()[j*MP1+i].real(),V->P_matrixNxK()[j*MP1+i].imag());OA
   }
  */


  initArnold(*V,*H,gauge);
  arnold(1,K,*V,*H,gauge);

  
  iter = 0; res = 1.;
  while ( (res > tol) && (iter<iterMax) ){
    arnold(K,M,*V,*H,gauge);
    H->eigLapack(*Qeig,W);
    sortEigenValues(*Qeig,W);

    /*
    printf("iter = %d ",iter);
    for(int i = 0 ; i < M ; i++)printf("%+e %+e\t",W.P_matrixKxK()[i*M+i].real(),W.P_matrixKxK()[i*M+i].imag());
    printf("\n");
    */
    qrShiftsRotations_v2(*V, *H, W);    


    /*
  for(int i = 0 ; i < M ; i++){
    for(int j = 0 ; j < M ; j++)
      printf("%+16.15e %+16.15e\t",H->P_matrixKxK()[i*M+j].real(),H->P_matrixKxK()[i*M+j].imag());
    printf("\n");
  }
    */
    //  H->eigLapack(*Qeig,W);
    //  sortEigenValues(*Qeig,W);
  
    /*
    memset(X.P_matrixNxK(),0,NL*MP1*sizeof(qkxTMComplex));
    for(int iN = 0 ; iN < NL ; iN++)
      for(int iK = 0 ; iK < K ; iK++)
	for(int iM = 0 ; iM < M ; iM++)
	  X.P_matrixNxK()[iN*MP1+iK] = X.P_matrixNxK()[iN*MP1+iK] + V->P_matrixNxK()[iN*MP1+iM] * Qeig->P_matrixKxK()[iM*M + M-iK-1];
    
    for(int k = 0 ; k < K ; k++){
      for(int l = 0 ; l < NL ; l++) vectmp[l] = X.P_matrixNxK()[l*MP1+k];
      double norma = sqrt(dotProduct(vectmp,vectmp,NL).real());
      for(int l = 0 ; l < NL ; l++)
	X.P_matrixNxK()[l*MP1+k] = X.P_matrixNxK()[l*MP1+k] *(1./norma);
    }
    

    for(int i = 0 ; i < M ; i++) ww[i] = W.P_matrixKxK()[i*M+i];
    for(int k = 0 ; k < K ; k++) W.P_matrixKxK()[k*M+k] = ww[M-k-1];
    */
  
    for(int l = 0 ; l < NL ; l++) vectmp[l] = V->P_matrixNxK()[l*MP1+M];
    res = sqrt(dotProduct(vectmp,vectmp,NL).real()); 
    printf("iter = %d, res = %16.15e\n",iter+1,res);
    // fprintf(ptr,"%d %e\n",iter+1,res);
    //  fflush(ptr);
    iter++;
  }
  

  stop = std::clock();
  double time = (stop-start)/ (double)CLOCKS_PER_SEC ;

  /*
  FILE *ptr;
  char filename[257];
  sprintf(filename,"time_K%d_P%d.dat",K,M-K);
  ptr = fopen(filename,"w");

  fprintf(ptr,"%e\n",time);
  */
  free(ww);
  free(vectmp);
  free(vi1);
  free(vi2);
  //  fclose(ptr);
  delete Qeig;
  delete V;
  delete H;
}
