/*
 * Functions which calculate partial derivatives in index space
 */

#include <field3d.hxx>
#include <utils.hxx>

#include <functional>

#include "stencil_loops.hxx"

/*
 * Define first derivative operators
 * 
 *
 * NOTE: Each scheme is a different struct, so that they
 * appear to compilers as distinct types. If they were functions,
 * then all would have the same type, so when used as a template
 * parameter would probably not lead to specialisation.
 */

/// central, 2nd order
struct DERIV_C2 {
  inline BoutReal operator()(const Stencil &f) const {
    return 0.5*(f.p - f.m);
  }
};
  
/// central, 4th order
struct DERIV_C4 {
  inline BoutReal operator()(const Stencil &f) const {
    return (8.*f.p - 8.*f.m + f.mm - f.pp)/12.;
  }
};
  
/*
 * First derivatives in X
 */
  
const Field2D indexDDX_C2(const Field2D &f) { 
  loopX1<DERIV_C2> op; return op(f); }
  
const Field3D indexDDX_C2(const Field3D &f) { 
  loopX1<DERIV_C2> op; return op(f); }

const Field2D indexDDX_C4(const Field2D &f) { 
  loopX1<DERIV_C4> op; return op(f); }

const Field3D indexDDX_C4(const Field3D &f) { 
  loopX1<DERIV_C4> op; return op(f); }

/*
 * First derivatives in Y
 */

const Field2D indexDDY_C2(const Field2D &f) { 
  loopY1<DERIV_C2> op; return op(f); }
  
const Field3D indexDDY_C2(const Field3D &f) { 
  loopY1<DERIV_C2> op; return op(f); }
  
/*
 * First derivatives in Z
 */
  
const Field2D indexDDZ_C2(const Field2D &f) { 
  loopZ1<DERIV_C2> op; return op(f); }
  
const Field3D indexDDZ_C2(const Field3D &f) { 
  loopZ1<DERIV_C2> op; return op(f); }

const Field2D indexDDZ_C4(const Field2D &f) { 
  loopZ1<DERIV_C4> op; return op(f); }

const Field3D indexDDZ_C4(const Field3D &f) { 
  loopZ1<DERIV_C4> op; return op(f); }
  
/*!
 * Index derivative in X, staggered fields
 */ 
const Field3D indexDDX_stag(const Field3D &f, CELL_LOC outloc) {
  throw BoutException("Staggered differencing not implemented yet");
}

/*!
 * Index derivative in X
 * 
 *
 */
typedef const Field3D(*deriv_func3d)(const Field3D &);

const Field3D indexDDX(const Field3D &f, CELL_LOC outloc) {
    
  static std::function<const Field3D(const Field3D &)> func = nullptr;
    
  if(mesh->StaggerGrids && (outloc != f.getLocation())) {
    // Some kind of staggered differencing
    return indexDDX_stag(f, outloc);
  }
  // Non-staggered differencing
    
  if(func == nullptr) {
    // Not defined yet
    string setting;
    Options::getRoot()->getSection("operators")->get("ddx", setting, "c2");
    setting = lowercase(setting);
      
    if(setting == "c2") {
      func = static_cast<deriv_func3d>(indexDDX_C2);
    }else if(setting == "c4") {
      func = static_cast<deriv_func3d>(indexDDX_C4);
    }else {
      throw BoutException("Unrecognised option for DDX: '%s'", setting.c_str());
    }
  }
  // Call the function
  return func(f);
}


