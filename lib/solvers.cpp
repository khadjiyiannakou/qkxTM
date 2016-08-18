#include <qkxTM.h> 
#include <lattice_util.h>
#include <blas_qkxTM.h>
#include <solvers.h>

using namespace qkxTM;

CG::CG(LatticeInfo *param , qkxTMComplex **gaugePointer){

  maxIter = param->maxIter;
  tol = param->tol;
  reliableDelta = param->reliableDelta;
  gauge = gaugePointer;
}

void CG::operator()(qkxTM::ColorSpinorField &out , qkxTM::ColorSpinorField &in , LatticeInfo *param){

  ColorSpinorField *r_P = new ColorSpinorField(param);
  ColorSpinorField *p_P = new ColorSpinorField(param);
  ColorSpinorField *Ap_P = new ColorSpinorField(param);
  ColorSpinorField *tmp_P = new ColorSpinorField(param);

  // take reference  
  ColorSpinorField &r = *r_P;
  ColorSpinorField &p = *p_P;
  ColorSpinorField &Ap = *Ap_P;
  ColorSpinorField &tmp = *tmp_P;

  // create the fields
  r.create();
  p.create();
  Ap.create();
  tmp.create();

  r.copy(in);  // r_ 0 = b - A x_0  initial guess zero
  p.copy(r);   // p_0 = r_0

  double r2;
  double r2_old;
  double stop;
  double src_norm;
  double pAp;
  double alpha;
  double beta;
  int k = 0;           // number of iterations

  r2 = reCdotXdagY(r,r);
  src_norm = reCdotXdagY(in,in);

  stop = src_norm*tol*tol;

  while( (r2 > stop) && (k < maxIter) ){

    fprintf(stdout,"CG : iteration = %d : r2 = %e\n",k , r2);
    // apply D^+ D                                             // A p 
    Ap.applyDslash( p.P_colorSpinor() , gauge);
    tmp.applyTwistAddDslash( Ap.P_colorSpinor() , p.P_colorSpinor() , param);
    Ap.applyDslashDag(tmp.P_colorSpinor(), gauge);
    Ap.applyTwistDagAddDslashDag(Ap.P_colorSpinor(), tmp.P_colorSpinor() , param);
    /////////////////////////////

    pAp = reCdotXdagY(p,Ap); // p A p
    
    alpha = r2 / pAp;         // a = r^T r / ( p^T A p)
    r2_old = r2;

    X_eq_X_p_aY(out, alpha, p);  // x = x + a p
    
    X_eq_X_p_aY(r, -alpha, Ap);   // r = r - a Ap

    r2 = reCdotXdagY(r,r);

    beta = r2/r2_old;       // b = r_(k+1)^T r_(k+1) / ( r_k^T r_k)

    X_eq_aX_p_Y(p, beta, r);  // p = r + beta p

    k++;

  } // close while loop

  printf("CG converged after iteration = %d\n",k);

  delete r_P;
  delete p_P;
  delete Ap_P;
  delete tmp_P;

}
