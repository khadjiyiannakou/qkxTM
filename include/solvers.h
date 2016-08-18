
class CG {

 private:
  int maxIter;
  double tol;
  double reliableDelta;
  qkxTMComplex **gauge;

 public:
  CG(LatticeInfo *param , qkxTMComplex **gaugePointer);
  ~CG() { ; }

  /////////// functor
  int MaxIter() { return maxIter; }
  double Tol() { return tol; }
  double ReliableDelta() { return reliableDelta; }
  ///////////////

  void operator()(qkxTM::ColorSpinorField &out , qkxTM::ColorSpinorField &in, LatticeInfo *param);
};
