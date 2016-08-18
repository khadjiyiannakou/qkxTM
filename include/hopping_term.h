#define ZERO(R)              \
  R->M[MS(0,0)].real() = 0.; \
  R->M[MS(0,0)].imag() = 0.; \
  R->M[MS(0,1)].real() = 0.; \
  R->M[MS(0,1)].imag() = 0.; \
  R->M[MS(0,2)].real() = 0.; \
  R->M[MS(0,2)].imag() = 0.; \
                             \
  R->M[MS(1,0)].real() = 0.; \
  R->M[MS(1,0)].imag() = 0.; \
  R->M[MS(1,1)].real() = 0.; \
  R->M[MS(1,1)].imag() = 0.; \
  R->M[MS(1,2)].real() = 0.; \
  R->M[MS(1,2)].imag() = 0.; \
                             \
  R->M[MS(2,0)].real() = 0.; \
  R->M[MS(2,0)].imag() = 0.; \
  R->M[MS(2,1)].real() = 0.; \
  R->M[MS(2,1)].imag() = 0.; \
  R->M[MS(2,2)].real() = 0.; \
  R->M[MS(2,2)].imag() = 0.; \
                             \
  R->M[MS(3,0)].real() = 0.; \
  R->M[MS(3,0)].imag() = 0.; \
  R->M[MS(3,1)].real() = 0.; \
  R->M[MS(3,1)].imag() = 0.; \
  R->M[MS(3,2)].real() = 0.; \
  R->M[MS(3,2)].imag() = 0.; 

#define PROJ_MINUS_X(phi,psi)						\
  phi->M[MS(0,0)].real() = psi[MS(0,0)][fwd[0]].real() + psi[MS(3,0)][fwd[0]].imag(); \
  phi->M[MS(0,0)].imag() = psi[MS(0,0)][fwd[0]].imag() - psi[MS(3,0)][fwd[0]].real(); \
  phi->M[MS(0,1)].real() = psi[MS(0,1)][fwd[0]].real() + psi[MS(3,1)][fwd[0]].imag(); \
  phi->M[MS(0,1)].imag() = psi[MS(0,1)][fwd[0]].imag() - psi[MS(3,1)][fwd[0]].real(); \
  phi->M[MS(0,2)].real() = psi[MS(0,2)][fwd[0]].real() + psi[MS(3,2)][fwd[0]].imag(); \
  phi->M[MS(0,2)].imag() = psi[MS(0,2)][fwd[0]].imag() - psi[MS(3,2)][fwd[0]].real(); \
									\
  phi->M[MS(1,0)].real() = psi[MS(1,0)][fwd[0]].real() + psi[MS(2,0)][fwd[0]].imag(); \
  phi->M[MS(1,0)].imag() = psi[MS(1,0)][fwd[0]].imag() - psi[MS(2,0)][fwd[0]].real(); \
  phi->M[MS(1,1)].real() = psi[MS(1,1)][fwd[0]].real() + psi[MS(2,1)][fwd[0]].imag(); \
  phi->M[MS(1,1)].imag() = psi[MS(1,1)][fwd[0]].imag() - psi[MS(2,1)][fwd[0]].real(); \
  phi->M[MS(1,2)].real() = psi[MS(1,2)][fwd[0]].real() + psi[MS(2,2)][fwd[0]].imag(); \
  phi->M[MS(1,2)].imag() = psi[MS(1,2)][fwd[0]].imag() - psi[MS(2,2)][fwd[0]].real(); 


#define PROJ_PLUS_X(phi,psi)						\
  phi->M[MS(0,0)].real() = psi[MS(0,0)][bwd[0]].real() - psi[MS(3,0)][bwd[0]].imag(); \
  phi->M[MS(0,0)].imag() = psi[MS(0,0)][bwd[0]].imag() + psi[MS(3,0)][bwd[0]].real(); \
  phi->M[MS(0,1)].real() = psi[MS(0,1)][bwd[0]].real() - psi[MS(3,1)][bwd[0]].imag(); \
  phi->M[MS(0,1)].imag() = psi[MS(0,1)][bwd[0]].imag() + psi[MS(3,1)][bwd[0]].real(); \
  phi->M[MS(0,2)].real() = psi[MS(0,2)][bwd[0]].real() - psi[MS(3,2)][bwd[0]].imag(); \
  phi->M[MS(0,2)].imag() = psi[MS(0,2)][bwd[0]].imag() + psi[MS(3,2)][bwd[0]].real(); \
									\
  phi->M[MS(1,0)].real() = psi[MS(1,0)][bwd[0]].real() - psi[MS(2,0)][bwd[0]].imag(); \
  phi->M[MS(1,0)].imag() = psi[MS(1,0)][bwd[0]].imag() + psi[MS(2,0)][bwd[0]].real(); \
  phi->M[MS(1,1)].real() = psi[MS(1,1)][bwd[0]].real() - psi[MS(2,1)][bwd[0]].imag(); \
  phi->M[MS(1,1)].imag() = psi[MS(1,1)][bwd[0]].imag() + psi[MS(2,1)][bwd[0]].real(); \
  phi->M[MS(1,2)].real() = psi[MS(1,2)][bwd[0]].real() - psi[MS(2,2)][bwd[0]].imag(); \
  phi->M[MS(1,2)].imag() = psi[MS(1,2)][bwd[0]].imag() + psi[MS(2,2)][bwd[0]].real(); 



#define PROJ_MINUS_Y(phi,psi)						\
  phi->M[MS(0,0)].real() = psi[MS(0,0)][fwd[1]].real() - psi[MS(3,0)][fwd[1]].real(); \
  phi->M[MS(0,0)].imag() = psi[MS(0,0)][fwd[1]].imag() - psi[MS(3,0)][fwd[1]].imag(); \
  phi->M[MS(0,1)].real() = psi[MS(0,1)][fwd[1]].real() - psi[MS(3,1)][fwd[1]].real(); \
  phi->M[MS(0,1)].imag() = psi[MS(0,1)][fwd[1]].imag() - psi[MS(3,1)][fwd[1]].imag(); \
  phi->M[MS(0,2)].real() = psi[MS(0,2)][fwd[1]].real() - psi[MS(3,2)][fwd[1]].real(); \
  phi->M[MS(0,2)].imag() = psi[MS(0,2)][fwd[1]].imag() - psi[MS(3,2)][fwd[1]].imag(); \
									\
  phi->M[MS(1,0)].real() = psi[MS(1,0)][fwd[1]].real() + psi[MS(2,0)][fwd[1]].real(); \
  phi->M[MS(1,0)].imag() = psi[MS(1,0)][fwd[1]].imag() + psi[MS(2,0)][fwd[1]].imag(); \
  phi->M[MS(1,1)].real() = psi[MS(1,1)][fwd[1]].real() + psi[MS(2,1)][fwd[1]].real(); \
  phi->M[MS(1,1)].imag() = psi[MS(1,1)][fwd[1]].imag() + psi[MS(2,1)][fwd[1]].imag(); \
  phi->M[MS(1,2)].real() = psi[MS(1,2)][fwd[1]].real() + psi[MS(2,2)][fwd[1]].real(); \
  phi->M[MS(1,2)].imag() = psi[MS(1,2)][fwd[1]].imag() + psi[MS(2,2)][fwd[1]].imag(); 


#define PROJ_PLUS_Y(phi,psi)						\
  phi->M[MS(0,0)].real() = psi[MS(0,0)][bwd[1]].real() + psi[MS(3,0)][bwd[1]].real(); \
  phi->M[MS(0,0)].imag() = psi[MS(0,0)][bwd[1]].imag() + psi[MS(3,0)][bwd[1]].imag(); \
  phi->M[MS(0,1)].real() = psi[MS(0,1)][bwd[1]].real() + psi[MS(3,1)][bwd[1]].real(); \
  phi->M[MS(0,1)].imag() = psi[MS(0,1)][bwd[1]].imag() + psi[MS(3,1)][bwd[1]].imag(); \
  phi->M[MS(0,2)].real() = psi[MS(0,2)][bwd[1]].real() + psi[MS(3,2)][bwd[1]].real(); \
  phi->M[MS(0,2)].imag() = psi[MS(0,2)][bwd[1]].imag() + psi[MS(3,2)][bwd[1]].imag(); \
									\
  phi->M[MS(1,0)].real() = psi[MS(1,0)][bwd[1]].real() - psi[MS(2,0)][bwd[1]].real(); \
  phi->M[MS(1,0)].imag() = psi[MS(1,0)][bwd[1]].imag() - psi[MS(2,0)][bwd[1]].imag(); \
  phi->M[MS(1,1)].real() = psi[MS(1,1)][bwd[1]].real() - psi[MS(2,1)][bwd[1]].real(); \
  phi->M[MS(1,1)].imag() = psi[MS(1,1)][bwd[1]].imag() - psi[MS(2,1)][bwd[1]].imag(); \
  phi->M[MS(1,2)].real() = psi[MS(1,2)][bwd[1]].real() - psi[MS(2,2)][bwd[1]].real(); \
  phi->M[MS(1,2)].imag() = psi[MS(1,2)][bwd[1]].imag() - psi[MS(2,2)][bwd[1]].imag(); 


#define PROJ_MINUS_Z(phi,psi)		      \
  phi->M[MS(0,0)].real() = psi[MS(0,0)][fwd[2]].real() + psi[MS(2,0)][fwd[2]].imag(); \
  phi->M[MS(0,0)].imag() = psi[MS(0,0)][fwd[2]].imag() - psi[MS(2,0)][fwd[2]].real(); \
  phi->M[MS(0,1)].real() = psi[MS(0,1)][fwd[2]].real() + psi[MS(2,1)][fwd[2]].imag(); \
  phi->M[MS(0,1)].imag() = psi[MS(0,1)][fwd[2]].imag() - psi[MS(2,1)][fwd[2]].real(); \
  phi->M[MS(0,2)].real() = psi[MS(0,2)][fwd[2]].real() + psi[MS(2,2)][fwd[2]].imag(); \
  phi->M[MS(0,2)].imag() = psi[MS(0,2)][fwd[2]].imag() - psi[MS(2,2)][fwd[2]].real(); \
									\
  phi->M[MS(1,0)].real() = psi[MS(1,0)][fwd[2]].real() - psi[MS(3,0)][fwd[2]].imag(); \
  phi->M[MS(1,0)].imag() = psi[MS(1,0)][fwd[2]].imag() + psi[MS(3,0)][fwd[2]].real(); \
  phi->M[MS(1,1)].real() = psi[MS(1,1)][fwd[2]].real() - psi[MS(3,1)][fwd[2]].imag(); \
  phi->M[MS(1,1)].imag() = psi[MS(1,1)][fwd[2]].imag() + psi[MS(3,1)][fwd[2]].real(); \
  phi->M[MS(1,2)].real() = psi[MS(1,2)][fwd[2]].real() - psi[MS(3,2)][fwd[2]].imag(); \
  phi->M[MS(1,2)].imag() = psi[MS(1,2)][fwd[2]].imag() + psi[MS(3,2)][fwd[2]].real(); 


#define PROJ_PLUS_Z(phi,psi)		      \
  phi->M[MS(0,0)].real() = psi[MS(0,0)][bwd[2]].real() - psi[MS(2,0)][bwd[2]].imag(); \
  phi->M[MS(0,0)].imag() = psi[MS(0,0)][bwd[2]].imag() + psi[MS(2,0)][bwd[2]].real(); \
  phi->M[MS(0,1)].real() = psi[MS(0,1)][bwd[2]].real() - psi[MS(2,1)][bwd[2]].imag(); \
  phi->M[MS(0,1)].imag() = psi[MS(0,1)][bwd[2]].imag() + psi[MS(2,1)][bwd[2]].real(); \
  phi->M[MS(0,2)].real() = psi[MS(0,2)][bwd[2]].real() - psi[MS(2,2)][bwd[2]].imag(); \
  phi->M[MS(0,2)].imag() = psi[MS(0,2)][bwd[2]].imag() + psi[MS(2,2)][bwd[2]].real(); \
									\
  phi->M[MS(1,0)].real() = psi[MS(1,0)][bwd[2]].real() + psi[MS(3,0)][bwd[2]].imag(); \
  phi->M[MS(1,0)].imag() = psi[MS(1,0)][bwd[2]].imag() - psi[MS(3,0)][bwd[2]].real(); \
  phi->M[MS(1,1)].real() = psi[MS(1,1)][bwd[2]].real() + psi[MS(3,1)][bwd[2]].imag(); \
  phi->M[MS(1,1)].imag() = psi[MS(1,1)][bwd[2]].imag() - psi[MS(3,1)][bwd[2]].real(); \
  phi->M[MS(1,2)].real() = psi[MS(1,2)][bwd[2]].real() + psi[MS(3,2)][bwd[2]].imag(); \
  phi->M[MS(1,2)].imag() = psi[MS(1,2)][bwd[2]].imag() - psi[MS(3,2)][bwd[2]].real(); 


#define PROJ_MINUS_T(phi,psi)			\
  phi->M[MS(0,0)].real() = 2*psi[MS(2,0)][fwd[3]].real();	\
  phi->M[MS(0,0)].imag() = 2*psi[MS(2,0)][fwd[3]].imag();	\
  phi->M[MS(0,1)].real() = 2*psi[MS(2,1)][fwd[3]].real();	\
  phi->M[MS(0,1)].imag() = 2*psi[MS(2,1)][fwd[3]].imag();	\
  phi->M[MS(0,2)].real() = 2*psi[MS(2,2)][fwd[3]].real();	\
  phi->M[MS(0,2)].imag() = 2*psi[MS(2,2)][fwd[3]].imag();	\
								\
  phi->M[MS(1,0)].real() = 2*psi[MS(3,0)][fwd[3]].real();	\
  phi->M[MS(1,0)].imag() = 2*psi[MS(3,0)][fwd[3]].imag();	\
  phi->M[MS(1,1)].real() = 2*psi[MS(3,1)][fwd[3]].real();	\
  phi->M[MS(1,1)].imag() = 2*psi[MS(3,1)][fwd[3]].imag();	\
  phi->M[MS(1,2)].real() = 2*psi[MS(3,2)][fwd[3]].real();	\
  phi->M[MS(1,2)].imag() = 2*psi[MS(3,2)][fwd[3]].imag();	


#define PROJ_PLUS_T(phi,psi)			\
  phi->M[MS(0,0)].real() = 2*psi[MS(0,0)][bwd[3]].real();	\
  phi->M[MS(0,0)].imag() = 2*psi[MS(0,0)][bwd[3]].imag();	\
  phi->M[MS(0,1)].real() = 2*psi[MS(0,1)][bwd[3]].real();	\
  phi->M[MS(0,1)].imag() = 2*psi[MS(0,1)][bwd[3]].imag();	\
  phi->M[MS(0,2)].real() = 2*psi[MS(0,2)][bwd[3]].real();	\
  phi->M[MS(0,2)].imag() = 2*psi[MS(0,2)][bwd[3]].imag();	\
								\
  phi->M[MS(1,0)].real() = 2*psi[MS(1,0)][bwd[3]].real();	\
  phi->M[MS(1,0)].imag() = 2*psi[MS(1,0)][bwd[3]].imag();	\
  phi->M[MS(1,1)].real() = 2*psi[MS(1,1)][bwd[3]].real();	\
  phi->M[MS(1,1)].imag() = 2*psi[MS(1,1)][bwd[3]].imag();	\
  phi->M[MS(1,2)].real() = 2*psi[MS(1,2)][bwd[3]].real();	\
  phi->M[MS(1,2)].imag() = 2*psi[MS(1,2)][bwd[3]].imag();	



//////////////////////////////////////////////////////////////////////

#define APPLY_LINK(xi,u,phi,mu)			\
  xi->M[MS(0,0)] = u[mu][MG(iv,0,0)] * phi->M[MS(0,0)] + u[mu][MG(iv,0,1)] * phi->M[MS(0,1)] + u[mu][MG(iv,0,2)] * phi->M[MS(0,2)]; \
  xi->M[MS(0,1)] = u[mu][MG(iv,1,0)] * phi->M[MS(0,0)] + u[mu][MG(iv,1,1)] * phi->M[MS(0,1)] + u[mu][MG(iv,1,2)] * phi->M[MS(0,2)]; \
  xi->M[MS(0,2)] = u[mu][MG(iv,2,0)] * phi->M[MS(0,0)] + u[mu][MG(iv,2,1)] * phi->M[MS(0,1)] + u[mu][MG(iv,2,2)] * phi->M[MS(0,2)]; \
									\
  xi->M[MS(1,0)] = u[mu][MG(iv,0,0)] * phi->M[MS(1,0)] + u[mu][MG(iv,0,1)] * phi->M[MS(1,1)] + u[mu][MG(iv,0,2)] * phi->M[MS(1,2)]; \
  xi->M[MS(1,1)] = u[mu][MG(iv,1,0)] * phi->M[MS(1,0)] + u[mu][MG(iv,1,1)] * phi->M[MS(1,1)] + u[mu][MG(iv,1,2)] * phi->M[MS(1,2)]; \
  xi->M[MS(1,2)] = u[mu][MG(iv,2,0)] * phi->M[MS(1,0)] + u[mu][MG(iv,2,1)] * phi->M[MS(1,1)] + u[mu][MG(iv,2,2)] * phi->M[MS(1,2)]; 

#define APPLY_LINK_DAG(xi,u,phi,mu)			\
  xi->M[MS(0,0)] = conj(u[mu][MG(bwd[mu],0,0)]) * phi->M[MS(0,0)] + conj(u[mu][MG(bwd[mu],1,0)]) * phi->M[MS(0,1)] + conj(u[mu][MG(bwd[mu],2,0)]) * phi->M[MS(0,2)]; \
  xi->M[MS(0,1)] = conj(u[mu][MG(bwd[mu],0,1)]) * phi->M[MS(0,0)] + conj(u[mu][MG(bwd[mu],1,1)]) * phi->M[MS(0,1)] + conj(u[mu][MG(bwd[mu],2,1)]) * phi->M[MS(0,2)]; \
  xi->M[MS(0,2)] = conj(u[mu][MG(bwd[mu],0,2)]) * phi->M[MS(0,0)] + conj(u[mu][MG(bwd[mu],1,2)]) * phi->M[MS(0,1)] + conj(u[mu][MG(bwd[mu],2,2)]) * phi->M[MS(0,2)]; \
									\
  xi->M[MS(1,0)] = conj(u[mu][MG(bwd[mu],0,0)]) * phi->M[MS(1,0)] + conj(u[mu][MG(bwd[mu],1,0)]) * phi->M[MS(1,1)] + conj(u[mu][MG(bwd[mu],2,0)]) * phi->M[MS(1,2)]; \
  xi->M[MS(1,1)] = conj(u[mu][MG(bwd[mu],0,1)]) * phi->M[MS(1,0)] + conj(u[mu][MG(bwd[mu],1,1)]) * phi->M[MS(1,1)] + conj(u[mu][MG(bwd[mu],2,1)]) * phi->M[MS(1,2)]; \
  xi->M[MS(1,2)] = conj(u[mu][MG(bwd[mu],0,2)]) * phi->M[MS(1,0)] + conj(u[mu][MG(bwd[mu],1,2)]) * phi->M[MS(1,1)] + conj(u[mu][MG(bwd[mu],2,2)]) * phi->M[MS(1,2)]; 



////////////////////////////////////////////////////////////////////////////////

#define COLLECT_MINUS_X(R,xi)			\
  R->M[MS(0,0)] = R->M[MS(0,0)] + xi->M[MS(0,0)];	\
  R->M[MS(0,1)] = R->M[MS(0,1)] + xi->M[MS(0,1)];	\
  R->M[MS(0,2)] = R->M[MS(0,2)] + xi->M[MS(0,2)];	\
							\
  R->M[MS(1,0)] = R->M[MS(1,0)] + xi->M[MS(1,0)];	\
  R->M[MS(1,1)] = R->M[MS(1,1)] + xi->M[MS(1,1)];	\
  R->M[MS(1,2)] = R->M[MS(1,2)] + xi->M[MS(1,2)];	\
							\
  R->M[MS(2,0)].real() = R->M[MS(2,0)].real() - xi->M[MS(1,0)].imag();	\
  R->M[MS(2,0)].imag() = R->M[MS(2,0)].imag() + xi->M[MS(1,0)].real();	\
  R->M[MS(2,1)].real() = R->M[MS(2,1)].real() - xi->M[MS(1,1)].imag();	\
  R->M[MS(2,1)].imag() = R->M[MS(2,1)].imag() + xi->M[MS(1,1)].real();	\
  R->M[MS(2,2)].real() = R->M[MS(2,2)].real() - xi->M[MS(1,2)].imag();	\
  R->M[MS(2,2)].imag() = R->M[MS(2,2)].imag() + xi->M[MS(1,2)].real();	\
									\
  R->M[MS(3,0)].real() = R->M[MS(3,0)].real() - xi->M[MS(0,0)].imag();  \
  R->M[MS(3,0)].imag() = R->M[MS(3,0)].imag() + xi->M[MS(0,0)].real();  \
  R->M[MS(3,1)].real() = R->M[MS(3,1)].real() - xi->M[MS(0,1)].imag();  \
  R->M[MS(3,1)].imag() = R->M[MS(3,1)].imag() + xi->M[MS(0,1)].real();  \
  R->M[MS(3,2)].real() = R->M[MS(3,2)].real() - xi->M[MS(0,2)].imag();  \
  R->M[MS(3,2)].imag() = R->M[MS(3,2)].imag() + xi->M[MS(0,2)].real();  

#define COLLECT_PLUS_X(R,xi)			\
  R->M[MS(0,0)] = R->M[MS(0,0)] + xi->M[MS(0,0)];	\
  R->M[MS(0,1)] = R->M[MS(0,1)] + xi->M[MS(0,1)];	\
  R->M[MS(0,2)] = R->M[MS(0,2)] + xi->M[MS(0,2)];	\
							\
  R->M[MS(1,0)] = R->M[MS(1,0)] + xi->M[MS(1,0)];	\
  R->M[MS(1,1)] = R->M[MS(1,1)] + xi->M[MS(1,1)];	\
  R->M[MS(1,2)] = R->M[MS(1,2)] + xi->M[MS(1,2)];	\
							\
  R->M[MS(2,0)].real() = R->M[MS(2,0)].real() + xi->M[MS(1,0)].imag();	\
  R->M[MS(2,0)].imag() = R->M[MS(2,0)].imag() - xi->M[MS(1,0)].real();	\
  R->M[MS(2,1)].real() = R->M[MS(2,1)].real() + xi->M[MS(1,1)].imag();	\
  R->M[MS(2,1)].imag() = R->M[MS(2,1)].imag() - xi->M[MS(1,1)].real();	\
  R->M[MS(2,2)].real() = R->M[MS(2,2)].real() + xi->M[MS(1,2)].imag();	\
  R->M[MS(2,2)].imag() = R->M[MS(2,2)].imag() - xi->M[MS(1,2)].real();	\
									\
  R->M[MS(3,0)].real() = R->M[MS(3,0)].real() + xi->M[MS(0,0)].imag();  \
  R->M[MS(3,0)].imag() = R->M[MS(3,0)].imag() - xi->M[MS(0,0)].real();  \
  R->M[MS(3,1)].real() = R->M[MS(3,1)].real() + xi->M[MS(0,1)].imag();  \
  R->M[MS(3,1)].imag() = R->M[MS(3,1)].imag() - xi->M[MS(0,1)].real();  \
  R->M[MS(3,2)].real() = R->M[MS(3,2)].real() + xi->M[MS(0,2)].imag();  \
  R->M[MS(3,2)].imag() = R->M[MS(3,2)].imag() - xi->M[MS(0,2)].real();  


#define COLLECT_MINUS_Y(R,xi)			\
  R->M[MS(0,0)] = R->M[MS(0,0)] + xi->M[MS(0,0)];	\
  R->M[MS(0,1)] = R->M[MS(0,1)] + xi->M[MS(0,1)];	\
  R->M[MS(0,2)] = R->M[MS(0,2)] + xi->M[MS(0,2)];	\
							\
  R->M[MS(1,0)] = R->M[MS(1,0)] + xi->M[MS(1,0)];	\
  R->M[MS(1,1)] = R->M[MS(1,1)] + xi->M[MS(1,1)];	\
  R->M[MS(1,2)] = R->M[MS(1,2)] + xi->M[MS(1,2)];	\
							\
  R->M[MS(2,0)].real() = R->M[MS(2,0)].real() + xi->M[MS(1,0)].real();	\
  R->M[MS(2,0)].imag() = R->M[MS(2,0)].imag() + xi->M[MS(1,0)].imag();	\
  R->M[MS(2,1)].real() = R->M[MS(2,1)].real() + xi->M[MS(1,1)].real();	\
  R->M[MS(2,1)].imag() = R->M[MS(2,1)].imag() + xi->M[MS(1,1)].imag();	\
  R->M[MS(2,2)].real() = R->M[MS(2,2)].real() + xi->M[MS(1,2)].real();	\
  R->M[MS(2,2)].imag() = R->M[MS(2,2)].imag() + xi->M[MS(1,2)].imag();	\
									\
  R->M[MS(3,0)].real() = R->M[MS(3,0)].real() - xi->M[MS(0,0)].real();	\
  R->M[MS(3,0)].imag() = R->M[MS(3,0)].imag() - xi->M[MS(0,0)].imag();	\
  R->M[MS(3,1)].real() = R->M[MS(3,1)].real() - xi->M[MS(0,1)].real();	\
  R->M[MS(3,1)].imag() = R->M[MS(3,1)].imag() - xi->M[MS(0,1)].imag();	\
  R->M[MS(3,2)].real() = R->M[MS(3,2)].real() - xi->M[MS(0,2)].real();	\
  R->M[MS(3,2)].imag() = R->M[MS(3,2)].imag() - xi->M[MS(0,2)].imag();	


#define COLLECT_PLUS_Y(R,xi)			\
  R->M[MS(0,0)] = R->M[MS(0,0)] + xi->M[MS(0,0)];	\
  R->M[MS(0,1)] = R->M[MS(0,1)] + xi->M[MS(0,1)];	\
  R->M[MS(0,2)] = R->M[MS(0,2)] + xi->M[MS(0,2)];	\
							\
  R->M[MS(1,0)] = R->M[MS(1,0)] + xi->M[MS(1,0)];	\
  R->M[MS(1,1)] = R->M[MS(1,1)] + xi->M[MS(1,1)];	\
  R->M[MS(1,2)] = R->M[MS(1,2)] + xi->M[MS(1,2)];	\
							\
  R->M[MS(2,0)].real() = R->M[MS(2,0)].real() - xi->M[MS(1,0)].real();	\
  R->M[MS(2,0)].imag() = R->M[MS(2,0)].imag() - xi->M[MS(1,0)].imag();	\
  R->M[MS(2,1)].real() = R->M[MS(2,1)].real() - xi->M[MS(1,1)].real();	\
  R->M[MS(2,1)].imag() = R->M[MS(2,1)].imag() - xi->M[MS(1,1)].imag();	\
  R->M[MS(2,2)].real() = R->M[MS(2,2)].real() - xi->M[MS(1,2)].real();	\
  R->M[MS(2,2)].imag() = R->M[MS(2,2)].imag() - xi->M[MS(1,2)].imag();	\
									\
  R->M[MS(3,0)].real() = R->M[MS(3,0)].real() + xi->M[MS(0,0)].real();	\
  R->M[MS(3,0)].imag() = R->M[MS(3,0)].imag() + xi->M[MS(0,0)].imag();	\
  R->M[MS(3,1)].real() = R->M[MS(3,1)].real() + xi->M[MS(0,1)].real();	\
  R->M[MS(3,1)].imag() = R->M[MS(3,1)].imag() + xi->M[MS(0,1)].imag();	\
  R->M[MS(3,2)].real() = R->M[MS(3,2)].real() + xi->M[MS(0,2)].real();	\
  R->M[MS(3,2)].imag() = R->M[MS(3,2)].imag() + xi->M[MS(0,2)].imag();	

#define COLLECT_MINUS_Z(R,xi)			\
  R->M[MS(0,0)] = R->M[MS(0,0)] + xi->M[MS(0,0)];	\
  R->M[MS(0,1)] = R->M[MS(0,1)] + xi->M[MS(0,1)];	\
  R->M[MS(0,2)] = R->M[MS(0,2)] + xi->M[MS(0,2)];	\
							\
  R->M[MS(1,0)] = R->M[MS(1,0)] + xi->M[MS(1,0)];	\
  R->M[MS(1,1)] = R->M[MS(1,1)] + xi->M[MS(1,1)];	\
  R->M[MS(1,2)] = R->M[MS(1,2)] + xi->M[MS(1,2)];	\
									\
  R->M[MS(2,0)].real() = R->M[MS(2,0)].real() - xi->M[MS(0,0)].imag();	\
  R->M[MS(2,0)].imag() = R->M[MS(2,0)].imag() + xi->M[MS(0,0)].real();	\
  R->M[MS(2,1)].real() = R->M[MS(2,1)].real() - xi->M[MS(0,1)].imag();	\
  R->M[MS(2,1)].imag() = R->M[MS(2,1)].imag() + xi->M[MS(0,1)].real();	\
  R->M[MS(2,2)].real() = R->M[MS(2,2)].real() - xi->M[MS(0,2)].imag();	\
  R->M[MS(2,2)].imag() = R->M[MS(2,2)].imag() + xi->M[MS(0,2)].real();	\
									\
  R->M[MS(3,0)].real() = R->M[MS(3,0)].real() + xi->M[MS(1,0)].imag();  \
  R->M[MS(3,0)].imag() = R->M[MS(3,0)].imag() - xi->M[MS(1,0)].real();  \
  R->M[MS(3,1)].real() = R->M[MS(3,1)].real() + xi->M[MS(1,1)].imag();  \
  R->M[MS(3,1)].imag() = R->M[MS(3,1)].imag() - xi->M[MS(1,1)].real();  \
  R->M[MS(3,2)].real() = R->M[MS(3,2)].real() + xi->M[MS(1,2)].imag();  \
  R->M[MS(3,2)].imag() = R->M[MS(3,2)].imag() - xi->M[MS(1,2)].real();  


#define COLLECT_PLUS_Z(R,xi)			\
  R->M[MS(0,0)] = R->M[MS(0,0)] + xi->M[MS(0,0)];	\
  R->M[MS(0,1)] = R->M[MS(0,1)] + xi->M[MS(0,1)];	\
  R->M[MS(0,2)] = R->M[MS(0,2)] + xi->M[MS(0,2)];	\
							\
  R->M[MS(1,0)] = R->M[MS(1,0)] + xi->M[MS(1,0)];	\
  R->M[MS(1,1)] = R->M[MS(1,1)] + xi->M[MS(1,1)];	\
  R->M[MS(1,2)] = R->M[MS(1,2)] + xi->M[MS(1,2)];	\
									\
  R->M[MS(2,0)].real() = R->M[MS(2,0)].real() + xi->M[MS(0,0)].imag();	\
  R->M[MS(2,0)].imag() = R->M[MS(2,0)].imag() - xi->M[MS(0,0)].real();	\
  R->M[MS(2,1)].real() = R->M[MS(2,1)].real() + xi->M[MS(0,1)].imag();	\
  R->M[MS(2,1)].imag() = R->M[MS(2,1)].imag() - xi->M[MS(0,1)].real();	\
  R->M[MS(2,2)].real() = R->M[MS(2,2)].real() + xi->M[MS(0,2)].imag();	\
  R->M[MS(2,2)].imag() = R->M[MS(2,2)].imag() - xi->M[MS(0,2)].real();	\
									\
  R->M[MS(3,0)].real() = R->M[MS(3,0)].real() - xi->M[MS(1,0)].imag();  \
  R->M[MS(3,0)].imag() = R->M[MS(3,0)].imag() + xi->M[MS(1,0)].real();  \
  R->M[MS(3,1)].real() = R->M[MS(3,1)].real() - xi->M[MS(1,1)].imag();  \
  R->M[MS(3,1)].imag() = R->M[MS(3,1)].imag() + xi->M[MS(1,1)].real();  \
  R->M[MS(3,2)].real() = R->M[MS(3,2)].real() - xi->M[MS(1,2)].imag();  \
  R->M[MS(3,2)].imag() = R->M[MS(3,2)].imag() + xi->M[MS(1,2)].real();  
  

#define COLLECT_MINUS_T(R,xi)			\
  R->M[MS(2,0)] = R->M[MS(2,0)] + xi->M[MS(0,0)];	\
  R->M[MS(2,1)] = R->M[MS(2,1)] + xi->M[MS(0,1)];	\
  R->M[MS(2,2)] = R->M[MS(2,2)] + xi->M[MS(0,2)];	\
							\
  R->M[MS(3,0)] = R->M[MS(3,0)] + xi->M[MS(1,0)];	\
  R->M[MS(3,1)] = R->M[MS(3,1)] + xi->M[MS(1,1)];	\
  R->M[MS(3,2)] = R->M[MS(3,2)] + xi->M[MS(1,2)];	

#define COLLECT_PLUS_T(R,xi)			\
  R->M[MS(0,0)] = R->M[MS(0,0)] + xi->M[MS(0,0)];	\
  R->M[MS(0,1)] = R->M[MS(0,1)] + xi->M[MS(0,1)];	\
  R->M[MS(0,2)] = R->M[MS(0,2)] + xi->M[MS(0,2)];	\
							\
  R->M[MS(1,0)] = R->M[MS(1,0)] + xi->M[MS(1,0)];	\
  R->M[MS(1,1)] = R->M[MS(1,1)] + xi->M[MS(1,1)];	\
  R->M[MS(1,2)] = R->M[MS(1,2)] + xi->M[MS(1,2)];	



//////////////////////////////////////////// for dagger /////////////////////////////////////////

#define DAG_PROJ_MINUS_X(phi,psi)						\
  phi->M[MS(0,0)].real() = psi[MS(0,0)][bwd[0]].real() + psi[MS(3,0)][bwd[0]].imag(); \
  phi->M[MS(0,0)].imag() = psi[MS(0,0)][bwd[0]].imag() - psi[MS(3,0)][bwd[0]].real(); \
  phi->M[MS(0,1)].real() = psi[MS(0,1)][bwd[0]].real() + psi[MS(3,1)][bwd[0]].imag(); \
  phi->M[MS(0,1)].imag() = psi[MS(0,1)][bwd[0]].imag() - psi[MS(3,1)][bwd[0]].real(); \
  phi->M[MS(0,2)].real() = psi[MS(0,2)][bwd[0]].real() + psi[MS(3,2)][bwd[0]].imag(); \
  phi->M[MS(0,2)].imag() = psi[MS(0,2)][bwd[0]].imag() - psi[MS(3,2)][bwd[0]].real(); \
									\
  phi->M[MS(1,0)].real() = psi[MS(1,0)][bwd[0]].real() + psi[MS(2,0)][bwd[0]].imag(); \
  phi->M[MS(1,0)].imag() = psi[MS(1,0)][bwd[0]].imag() - psi[MS(2,0)][bwd[0]].real(); \
  phi->M[MS(1,1)].real() = psi[MS(1,1)][bwd[0]].real() + psi[MS(2,1)][bwd[0]].imag(); \
  phi->M[MS(1,1)].imag() = psi[MS(1,1)][bwd[0]].imag() - psi[MS(2,1)][bwd[0]].real(); \
  phi->M[MS(1,2)].real() = psi[MS(1,2)][bwd[0]].real() + psi[MS(2,2)][bwd[0]].imag(); \
  phi->M[MS(1,2)].imag() = psi[MS(1,2)][bwd[0]].imag() - psi[MS(2,2)][bwd[0]].real(); 


#define DAG_PROJ_PLUS_X(phi,psi)						\
  phi->M[MS(0,0)].real() = psi[MS(0,0)][fwd[0]].real() - psi[MS(3,0)][fwd[0]].imag(); \
  phi->M[MS(0,0)].imag() = psi[MS(0,0)][fwd[0]].imag() + psi[MS(3,0)][fwd[0]].real(); \
  phi->M[MS(0,1)].real() = psi[MS(0,1)][fwd[0]].real() - psi[MS(3,1)][fwd[0]].imag(); \
  phi->M[MS(0,1)].imag() = psi[MS(0,1)][fwd[0]].imag() + psi[MS(3,1)][fwd[0]].real(); \
  phi->M[MS(0,2)].real() = psi[MS(0,2)][fwd[0]].real() - psi[MS(3,2)][fwd[0]].imag(); \
  phi->M[MS(0,2)].imag() = psi[MS(0,2)][fwd[0]].imag() + psi[MS(3,2)][fwd[0]].real(); \
									\
  phi->M[MS(1,0)].real() = psi[MS(1,0)][fwd[0]].real() - psi[MS(2,0)][fwd[0]].imag(); \
  phi->M[MS(1,0)].imag() = psi[MS(1,0)][fwd[0]].imag() + psi[MS(2,0)][fwd[0]].real(); \
  phi->M[MS(1,1)].real() = psi[MS(1,1)][fwd[0]].real() - psi[MS(2,1)][fwd[0]].imag(); \
  phi->M[MS(1,1)].imag() = psi[MS(1,1)][fwd[0]].imag() + psi[MS(2,1)][fwd[0]].real(); \
  phi->M[MS(1,2)].real() = psi[MS(1,2)][fwd[0]].real() - psi[MS(2,2)][fwd[0]].imag(); \
  phi->M[MS(1,2)].imag() = psi[MS(1,2)][fwd[0]].imag() + psi[MS(2,2)][fwd[0]].real(); 



#define DAG_PROJ_MINUS_Y(phi,psi)						\
  phi->M[MS(0,0)].real() = psi[MS(0,0)][bwd[1]].real() - psi[MS(3,0)][bwd[1]].real(); \
  phi->M[MS(0,0)].imag() = psi[MS(0,0)][bwd[1]].imag() - psi[MS(3,0)][bwd[1]].imag(); \
  phi->M[MS(0,1)].real() = psi[MS(0,1)][bwd[1]].real() - psi[MS(3,1)][bwd[1]].real(); \
  phi->M[MS(0,1)].imag() = psi[MS(0,1)][bwd[1]].imag() - psi[MS(3,1)][bwd[1]].imag(); \
  phi->M[MS(0,2)].real() = psi[MS(0,2)][bwd[1]].real() - psi[MS(3,2)][bwd[1]].real(); \
  phi->M[MS(0,2)].imag() = psi[MS(0,2)][bwd[1]].imag() - psi[MS(3,2)][bwd[1]].imag(); \
									\
  phi->M[MS(1,0)].real() = psi[MS(1,0)][bwd[1]].real() + psi[MS(2,0)][bwd[1]].real(); \
  phi->M[MS(1,0)].imag() = psi[MS(1,0)][bwd[1]].imag() + psi[MS(2,0)][bwd[1]].imag(); \
  phi->M[MS(1,1)].real() = psi[MS(1,1)][bwd[1]].real() + psi[MS(2,1)][bwd[1]].real(); \
  phi->M[MS(1,1)].imag() = psi[MS(1,1)][bwd[1]].imag() + psi[MS(2,1)][bwd[1]].imag(); \
  phi->M[MS(1,2)].real() = psi[MS(1,2)][bwd[1]].real() + psi[MS(2,2)][bwd[1]].real(); \
  phi->M[MS(1,2)].imag() = psi[MS(1,2)][bwd[1]].imag() + psi[MS(2,2)][bwd[1]].imag(); 


#define DAG_PROJ_PLUS_Y(phi,psi)						\
  phi->M[MS(0,0)].real() = psi[MS(0,0)][fwd[1]].real() + psi[MS(3,0)][fwd[1]].real(); \
  phi->M[MS(0,0)].imag() = psi[MS(0,0)][fwd[1]].imag() + psi[MS(3,0)][fwd[1]].imag(); \
  phi->M[MS(0,1)].real() = psi[MS(0,1)][fwd[1]].real() + psi[MS(3,1)][fwd[1]].real(); \
  phi->M[MS(0,1)].imag() = psi[MS(0,1)][fwd[1]].imag() + psi[MS(3,1)][fwd[1]].imag(); \
  phi->M[MS(0,2)].real() = psi[MS(0,2)][fwd[1]].real() + psi[MS(3,2)][fwd[1]].real(); \
  phi->M[MS(0,2)].imag() = psi[MS(0,2)][fwd[1]].imag() + psi[MS(3,2)][fwd[1]].imag(); \
									\
  phi->M[MS(1,0)].real() = psi[MS(1,0)][fwd[1]].real() - psi[MS(2,0)][fwd[1]].real(); \
  phi->M[MS(1,0)].imag() = psi[MS(1,0)][fwd[1]].imag() - psi[MS(2,0)][fwd[1]].imag(); \
  phi->M[MS(1,1)].real() = psi[MS(1,1)][fwd[1]].real() - psi[MS(2,1)][fwd[1]].real(); \
  phi->M[MS(1,1)].imag() = psi[MS(1,1)][fwd[1]].imag() - psi[MS(2,1)][fwd[1]].imag(); \
  phi->M[MS(1,2)].real() = psi[MS(1,2)][fwd[1]].real() - psi[MS(2,2)][fwd[1]].real(); \
  phi->M[MS(1,2)].imag() = psi[MS(1,2)][fwd[1]].imag() - psi[MS(2,2)][fwd[1]].imag(); 


#define DAG_PROJ_MINUS_Z(phi,psi)		      \
  phi->M[MS(0,0)].real() = psi[MS(0,0)][bwd[2]].real() + psi[MS(2,0)][bwd[2]].imag(); \
  phi->M[MS(0,0)].imag() = psi[MS(0,0)][bwd[2]].imag() - psi[MS(2,0)][bwd[2]].real(); \
  phi->M[MS(0,1)].real() = psi[MS(0,1)][bwd[2]].real() + psi[MS(2,1)][bwd[2]].imag(); \
  phi->M[MS(0,1)].imag() = psi[MS(0,1)][bwd[2]].imag() - psi[MS(2,1)][bwd[2]].real(); \
  phi->M[MS(0,2)].real() = psi[MS(0,2)][bwd[2]].real() + psi[MS(2,2)][bwd[2]].imag(); \
  phi->M[MS(0,2)].imag() = psi[MS(0,2)][bwd[2]].imag() - psi[MS(2,2)][bwd[2]].real(); \
									\
  phi->M[MS(1,0)].real() = psi[MS(1,0)][bwd[2]].real() - psi[MS(3,0)][bwd[2]].imag(); \
  phi->M[MS(1,0)].imag() = psi[MS(1,0)][bwd[2]].imag() + psi[MS(3,0)][bwd[2]].real(); \
  phi->M[MS(1,1)].real() = psi[MS(1,1)][bwd[2]].real() - psi[MS(3,1)][bwd[2]].imag(); \
  phi->M[MS(1,1)].imag() = psi[MS(1,1)][bwd[2]].imag() + psi[MS(3,1)][bwd[2]].real(); \
  phi->M[MS(1,2)].real() = psi[MS(1,2)][bwd[2]].real() - psi[MS(3,2)][bwd[2]].imag(); \
  phi->M[MS(1,2)].imag() = psi[MS(1,2)][bwd[2]].imag() + psi[MS(3,2)][bwd[2]].real(); 


#define DAG_PROJ_PLUS_Z(phi,psi)		      \
  phi->M[MS(0,0)].real() = psi[MS(0,0)][fwd[2]].real() - psi[MS(2,0)][fwd[2]].imag(); \
  phi->M[MS(0,0)].imag() = psi[MS(0,0)][fwd[2]].imag() + psi[MS(2,0)][fwd[2]].real(); \
  phi->M[MS(0,1)].real() = psi[MS(0,1)][fwd[2]].real() - psi[MS(2,1)][fwd[2]].imag(); \
  phi->M[MS(0,1)].imag() = psi[MS(0,1)][fwd[2]].imag() + psi[MS(2,1)][fwd[2]].real(); \
  phi->M[MS(0,2)].real() = psi[MS(0,2)][fwd[2]].real() - psi[MS(2,2)][fwd[2]].imag(); \
  phi->M[MS(0,2)].imag() = psi[MS(0,2)][fwd[2]].imag() + psi[MS(2,2)][fwd[2]].real(); \
									\
  phi->M[MS(1,0)].real() = psi[MS(1,0)][fwd[2]].real() + psi[MS(3,0)][fwd[2]].imag(); \
  phi->M[MS(1,0)].imag() = psi[MS(1,0)][fwd[2]].imag() - psi[MS(3,0)][fwd[2]].real(); \
  phi->M[MS(1,1)].real() = psi[MS(1,1)][fwd[2]].real() + psi[MS(3,1)][fwd[2]].imag(); \
  phi->M[MS(1,1)].imag() = psi[MS(1,1)][fwd[2]].imag() - psi[MS(3,1)][fwd[2]].real(); \
  phi->M[MS(1,2)].real() = psi[MS(1,2)][fwd[2]].real() + psi[MS(3,2)][fwd[2]].imag(); \
  phi->M[MS(1,2)].imag() = psi[MS(1,2)][fwd[2]].imag() - psi[MS(3,2)][fwd[2]].real(); 


#define DAG_PROJ_MINUS_T(phi,psi)			\
  phi->M[MS(0,0)].real() = 2*psi[MS(2,0)][bwd[3]].real();	\
  phi->M[MS(0,0)].imag() = 2*psi[MS(2,0)][bwd[3]].imag();	\
  phi->M[MS(0,1)].real() = 2*psi[MS(2,1)][bwd[3]].real();	\
  phi->M[MS(0,1)].imag() = 2*psi[MS(2,1)][bwd[3]].imag();	\
  phi->M[MS(0,2)].real() = 2*psi[MS(2,2)][bwd[3]].real();	\
  phi->M[MS(0,2)].imag() = 2*psi[MS(2,2)][bwd[3]].imag();	\
								\
  phi->M[MS(1,0)].real() = 2*psi[MS(3,0)][bwd[3]].real();	\
  phi->M[MS(1,0)].imag() = 2*psi[MS(3,0)][bwd[3]].imag();	\
  phi->M[MS(1,1)].real() = 2*psi[MS(3,1)][bwd[3]].real();	\
  phi->M[MS(1,1)].imag() = 2*psi[MS(3,1)][bwd[3]].imag();	\
  phi->M[MS(1,2)].real() = 2*psi[MS(3,2)][bwd[3]].real();	\
  phi->M[MS(1,2)].imag() = 2*psi[MS(3,2)][bwd[3]].imag();	


#define DAG_PROJ_PLUS_T(phi,psi)			\
  phi->M[MS(0,0)].real() = 2*psi[MS(0,0)][fwd[3]].real();	\
  phi->M[MS(0,0)].imag() = 2*psi[MS(0,0)][fwd[3]].imag();	\
  phi->M[MS(0,1)].real() = 2*psi[MS(0,1)][fwd[3]].real();	\
  phi->M[MS(0,1)].imag() = 2*psi[MS(0,1)][fwd[3]].imag();	\
  phi->M[MS(0,2)].real() = 2*psi[MS(0,2)][fwd[3]].real();	\
  phi->M[MS(0,2)].imag() = 2*psi[MS(0,2)][fwd[3]].imag();	\
								\
  phi->M[MS(1,0)].real() = 2*psi[MS(1,0)][fwd[3]].real();	\
  phi->M[MS(1,0)].imag() = 2*psi[MS(1,0)][fwd[3]].imag();	\
  phi->M[MS(1,1)].real() = 2*psi[MS(1,1)][fwd[3]].real();	\
  phi->M[MS(1,1)].imag() = 2*psi[MS(1,1)][fwd[3]].imag();	\
  phi->M[MS(1,2)].real() = 2*psi[MS(1,2)][fwd[3]].real();	\
  phi->M[MS(1,2)].imag() = 2*psi[MS(1,2)][fwd[3]].imag();	


