#ifndef _LATTICE_UTIL_H
#define _LATTICE_UTIL_H

#include <qkxTM.h>

// need to create a class which include all information about the lattice


namespace qkxTM {

  class LatticeField {

  protected:
    int volume;  // lattice volume
    int volumeCB; // check board volume is the half volume
    int surface[NDIM];
    int surfaceCB[NDIM];
    size_t totalBytes; // memory size of the spinor
    int latticeLength[NDIM];
    int nspins;
    int ncolors;
    int ndf;
  public:
    LatticeField(LatticeInfo *param);  // LatticeInfo will be a structure wich encuplsulate all the info 
    ~LatticeField() {;}

    // flag to check normal order or even odd order
    bool evenOdd_order;

    // Functors
    int Volume() const { return volume; } 
    int VolumeCB() const { return volumeCB; }
    int Surface(const int i) const { return surface[i]; } 
    int SurfaceCB(const int i) const { return surfaceCB[i]; }
    int TotalBytes() const { return totalBytes; }
    int LatticeLength(const int i) const { return latticeLength[i]; }

    // Methods
    void getNormalCoords(const int latticePoint, int *x) const;
    void getFwdCoords(const int *x, int *fwd) const;
    void getBwdCoords(const int *x, int *bwd) const;
    void getFwdBwdCoords(const int *x, int fwdBwd[NDIM][NDIM]) const;
    void getBwdFwdCoords(const int *x, int bwdFwd[NDIM][NDIM]) const;


  };

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  class GaugeField : public LatticeField {
    
  private:
    qkxTMComplex** p_gauge;

  public:
    void create(); // this function will return an allocated pointer to gaugefield 
    GaugeField(LatticeInfo *param);
    ~GaugeField();
    
    qkxTMComplex** P_gauge(){ return p_gauge;}
    void copy(GaugeField &src);
    void zero();
    void applyBoundary(LatticeInfo *latInfo);
    double calculatePlaquette();
    void APEsmearing(qkxTMComplex **U_tmp, qkxTMComplex **U_ape, LatticeInfo *latInfo);
    void unitField();

  };
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  class ColorSpinorField : public LatticeField {

  private:
    qkxTMComplex** p_colorSpinor;

  public:
    void create();
    ColorSpinorField(LatticeInfo *param);
    ~ColorSpinorField();

    qkxTMComplex** P_colorSpinor() { return p_colorSpinor; }
    void copy(ColorSpinorField &src);
    void zero();
    void gaussianSmearing(qkxTMComplex **Psi_tmp, qkxTMComplex **Psi_sm, qkxTMComplex **gauge , LatticeInfo *latInfo);
    void applyDslash(qkxTMComplex **Psi, qkxTMComplex **gauge);
    void applyTwistAddDslash(qkxTMComplex **Dslash, qkxTMComplex **Psi, LatticeInfo *latInfo);
    void applyDslashDag(qkxTMComplex **Psi, qkxTMComplex **gauge);
    void applyTwistDagAddDslashDag(qkxTMComplex **DslashDag, qkxTMComplex **Psi, LatticeInfo *latInfo);
    void MdagM(qkxTMComplex **vec,qkxTMComplex **gauge, LatticeInfo *param);

  };

} // close namespace


#endif
