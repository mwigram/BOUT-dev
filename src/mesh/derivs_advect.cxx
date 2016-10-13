#include "stencil_loops.hxx"
#include <utils.hxx>

/*
 * Define advection operators
 */

/// Upwinding, first order
struct ADVECT_U1 {
  inline BoutReal operator()(BoutReal v, const Stencil &f) const {
    return (v >= 0.0) ? v*(f.c - f.m): v*(f.p - f.c);
  }
};

/// Upwinding, second order
struct ADVECT_U2 {
  inline BoutReal operator()(BoutReal v, const Stencil &f) const {
    return v>=0.0 ? v*(1.5*f.c - 2.0*f.m + 0.5*f.mm): v*(-0.5*f.pp + 2.0*f.p - 1.5*f.c);
  }
};

/// Central, second order
struct ADVECT_C2 {
  inline BoutReal operator()(BoutReal v, const Stencil &f) const {
    return v*0.5*(f.p - f.m);
  }
};

/// Central, fourth order
struct ADVECT_C4 {
  inline BoutReal operator()(BoutReal v, const Stencil &f) const {
    return v*(8.*f.p - 8.*f.m + f.mm - f.pp)/12.;
  }
};

/// WENO, 3rd order
struct ADVECT_W3 {
  inline BoutReal operator()(BoutReal v, const Stencil &f) const {
    BoutReal deriv, w, r;

    if(v > 0.0) {
      // Left-biased stencil
      
      r = (WENO_SMALL + SQ(f.c - 2.0*f.m + f.mm)) / (WENO_SMALL + SQ(f.p - 2.0*f.c + f.m));
      w = 1.0 / (1.0 + 2.0*r*r);
      
      deriv = 0.5*(f.p - f.m) - 0.5*w*(-f.mm + 3.*f.m - 3.*f.c + f.p);
      
    }else {
      // Right-biased
      
      r = (WENO_SMALL + SQ(f.pp - 2.0*f.p + f.c)) / (WENO_SMALL + SQ(f.p - 2.0*f.c + f.m));
      w = 1.0 / (1.0 + 2.0*r*r);
      
      deriv = 0.5*(f.p - f.m) - 0.5*w*( -f.m + 3.*f.c - 3.*f.p + f.pp );
    }
    
    return v*deriv;
  }

  const BoutReal WENO_SMALL = 1.0e-8; // Small number for WENO schemes
};

/*
 * Advection in X
 */

const Field2D indexVDDX_U1(const Field2D &v, const Field2D &f) { loopXRS<ADVECT_U1> op; return op(v, f); }
const Field2D indexVDDX_U2(const Field2D &v, const Field2D &f) { loopXRS<ADVECT_U2> op; return op(v, f); }
const Field2D indexVDDX_C2(const Field2D &v, const Field2D &f) { loopXRS<ADVECT_C2> op; return op(v, f); }
const Field2D indexVDDX_C4(const Field2D &v, const Field2D &f) { loopXRS<ADVECT_C4> op; return op(v, f); }
const Field2D indexVDDX_W3(const Field2D &v, const Field2D &f) { loopXRS<ADVECT_W3> op; return op(v, f); }

const Field3D indexVDDX_U1(const Field3D &v, const Field2D &f) { loopXRS<ADVECT_U1> op; return op(v, f); }
const Field3D indexVDDX_U2(const Field3D &v, const Field2D &f) { loopXRS<ADVECT_U2> op; return op(v, f); }
const Field3D indexVDDX_C2(const Field3D &v, const Field2D &f) { loopXRS<ADVECT_C2> op; return op(v, f); }
const Field3D indexVDDX_C4(const Field3D &v, const Field2D &f) { loopXRS<ADVECT_C4> op; return op(v, f); }
const Field3D indexVDDX_W3(const Field3D &v, const Field2D &f) { loopXRS<ADVECT_W3> op; return op(v, f); }

const Field3D indexVDDX_U1(const Field2D &v, const Field3D &f) { loopXRS<ADVECT_U1> op; return op(v, f); }
const Field3D indexVDDX_U2(const Field2D &v, const Field3D &f) { loopXRS<ADVECT_U2> op; return op(v, f); }
const Field3D indexVDDX_C2(const Field2D &v, const Field3D &f) { loopXRS<ADVECT_C2> op; return op(v, f); }
const Field3D indexVDDX_C4(const Field2D &v, const Field3D &f) { loopXRS<ADVECT_C4> op; return op(v, f); }
const Field3D indexVDDX_W3(const Field2D &v, const Field3D &f) { loopXRS<ADVECT_W3> op; return op(v, f); }

const Field3D indexVDDX_U1(const Field3D &v, const Field3D &f) { loopXRS<ADVECT_U1> op; return op(v, f); }
const Field3D indexVDDX_U2(const Field3D &v, const Field3D &f) { loopXRS<ADVECT_U2> op; return op(v, f); }
const Field3D indexVDDX_C2(const Field3D &v, const Field3D &f) { loopXRS<ADVECT_C2> op; return op(v, f); }
const Field3D indexVDDX_C4(const Field3D &v, const Field3D &f) { loopXRS<ADVECT_C4> op; return op(v, f); }
const Field3D indexVDDX_W3(const Field3D &v, const Field3D &f) { loopXRS<ADVECT_W3> op; return op(v, f); }

/*!
 * Choose one of the indexVDDX functions, based on options
 *
 * Template parameter is the field type (Field2D, Field3D)
 */
template<typename V, typename F>
void defaultIndexVDDX(std::function<const F(const V &, const F &)> *func) {
  // Define the type of the free function needed, for overload resolution
  using func_t = const F(*)(const V &, const F &);

  // Fetch the string specifying the method
  string setting;
  Options::getRoot()->getSection("operators")->get("vddx", setting, "u1");
  setting = lowercase(setting);
  
  if(setting == "u1") {
    // Second order central differencing
    *func = static_cast<func_t>(indexVDDX_U1);
    
  }else if(setting == "u2") {
    // Fourth order central differencing
    *func = static_cast<func_t>(indexVDDX_U2);
    
  }else if(setting == "c2") {
    // Fourth order central differencing
    *func = static_cast<func_t>(indexVDDX_C2);
    
  }else if(setting == "c4") {
    // Fourth order central differencing
    *func = static_cast<func_t>(indexVDDX_C4);
    
  }else if(setting == "w3") {
    // Fourth order central differencing
    *func = static_cast<func_t>(indexVDDX_W3);
  }else {
    throw BoutException("Unrecognised option for DDX: '%s'", setting.c_str());
  }
}

const Field2D indexVDDX(const Field2D &v, const Field2D &f) {
  static std::function<const Field2D(const Field2D &, const Field2D &)> func = nullptr;
  if(func == nullptr) {
    defaultIndexVDDX(&func);
  }
  return func(v,f);
}
