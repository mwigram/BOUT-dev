/*
 * Functions which calculate partial derivatives in index space
 */

#include <field3d.hxx>
#include <utils.hxx>

#include <functional>

#include "stencil_loops.hxx"

#include "derivs_first.hxx"

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
  loopXS<DERIV_C2> op; return op(f); }
  
const Field3D indexDDX_C2(const Field3D &f) { 
  loopXS<DERIV_C2> op; return op(f); }

const Field2D indexDDX_C4(const Field2D &f) { 
  loopXS<DERIV_C4> op; return op(f); }

const Field3D indexDDX_C4(const Field3D &f) { 
  loopXS<DERIV_C4> op; return op(f); }

/*!
 * Choose one of the indexDDX functions, based on options
 *
 * Template parameter is the field type (Field2D, Field3D)
 */
template<typename F>
void defaultIndexDDX(std::function<const F(const F &)> *func) {
  // Define the type of the free function needed, for overload resolution
  using func_t = const F(*)(const F &);

  // Fetch the string specifying the method
  string setting;
  Options::getRoot()->getSection("operators")->get("ddx", setting, "c2");
  setting = lowercase(setting);
  
  if(setting == "c2") {
    // Second order central differencing
    *func = static_cast<func_t>(indexDDX_C2);
    
  }else if(setting == "c4") {
    // Fourth order central differencing
    *func = static_cast<func_t>(indexDDX_C4);
    
  }else {
    throw BoutException("Unrecognised option for DDX: '%s'", setting.c_str());
  }
}

/*
 * First derivatives in Y
 */

const Field2D indexDDY_C2(const Field2D &f) { 
  loopY1<DERIV_C2> op; return op(f); }
  
const Field3D indexDDY_C2(const Field3D &f) { 
  loopY1<DERIV_C2> op; return op(f); }

/*!
 * Choose one of the indexDDY functions, based on options
 *
 * Template parameter is the field type (Field2D, Field3D)
 */
template<typename F>
void defaultIndexDDY(std::function<const F(const F &)> *func) {
  using func_t = const F(*)(const F &);
  
  // Fetch the string specifying the method
  string setting;
  Options::getRoot()->getSection("operators")->get("ddy", setting, "c2");
  setting = lowercase(setting);
  
  if(setting == "c2") {
    // Second order central differencing
    *func = static_cast<func_t>(indexDDY_C2);
    
  }else {
    throw BoutException("Unrecognised option for DDY: '%s'", setting.c_str());
  }
}

/*
 * First derivatives in Z
 *
 * Derivatives of Field2D defined in header (return 0.0)
 */
  
const Field3D indexDDZ_C2(const Field3D &f) { 
  loopZ1<DERIV_C2> op; return op(f); }

const Field3D indexDDZ_C4(const Field3D &f) { 
  loopZ1<DERIV_C4> op; return op(f); }

/*!
 * Choose one of the indexDDZ functions, based on options
 *
 * Template parameter is the field type (Field2D, Field3D)
 */
template<typename F>
void defaultIndexDDZ(std::function<const F(const F &)> *func) {
  // Define the type of the free function needed, for overload resolution
  using func_t = const F(*)(const F &);

  // Fetch the string specifying the method
  string setting;
  Options::getRoot()->getSection("operators")->get("ddz", setting, "c2");
  setting = lowercase(setting);
  
  if(setting == "c2") {
    // Second order central differencing
    *func = static_cast<func_t>(indexDDZ_C2);
    
  }else if(setting == "c4") {
    // Fourth order central differencing
    *func = static_cast<func_t>(indexDDZ_C4);
    
  }else {
    throw BoutException("Unrecognised option for DDX: '%s'", setting.c_str());
  }
}

/*!
 * Index derivative in X, staggered fields
 */ 
const Field3D indexDDX_stag(const Field3D &f, CELL_LOC outloc) {
  throw BoutException("Staggered differencing not implemented yet");
}
const Field2D indexDDX_stag(const Field2D &f, CELL_LOC outloc) {
  throw BoutException("Staggered differencing not implemented yet");
}

const Field3D indexDDY_stag(const Field3D &f, CELL_LOC outloc) {
  throw BoutException("Staggered differencing not implemented yet");
}
const Field2D indexDDY_stag(const Field2D &f, CELL_LOC outloc) {
  throw BoutException("Staggered differencing not implemented yet");
}

const Field3D indexDDZ_stag(const Field3D &f, CELL_LOC outloc) {
  throw BoutException("Staggered differencing not implemented yet");
}
const Field2D indexDDZ_stag(const Field2D &f, CELL_LOC outloc) {
  throw BoutException("Staggered differencing not implemented yet");
}

/*!
 * Index derivatives
 * 
 *
 */

#define F_INDEX_FIRSTDERIV(name, defaultFunc, ftype)                                   \
const ftype name(const ftype &f, CELL_LOC outloc) {                                    \
  static std::function<const ftype(const ftype &)> func = nullptr;                     \
  if(mesh->StaggerGrids && (outloc != CELL_DEFAULT) && (outloc != f.getLocation())) {  \
    /* Some kind of staggered differencing  */                                         \
    return name ## _stag(f, outloc);                                                   \
  }                                                                                    \
  /* Non-staggered differencing */                                                     \
  if(func == nullptr) {                                                                \
    /* Not defined yet, so choose default */                                           \
    defaultFunc(&func);                                                                \
  }                                                                                    \
  return func(f);                                                                      \
}

F_INDEX_FIRSTDERIV(indexDDX, defaultIndexDDX, Field3D); // const Field3D indexDDX(const Field3D &f, CELL_LOC outloc)
F_INDEX_FIRSTDERIV(indexDDX, defaultIndexDDX, Field2D); // const Field2D indexDDX(const Field2D &f, CELL_LOC outloc)

F_INDEX_FIRSTDERIV(indexDDY, defaultIndexDDY, Field3D); // const Field3D indexDDY(const Field3D &f, CELL_LOC outloc)
F_INDEX_FIRSTDERIV(indexDDY, defaultIndexDDY, Field2D); // const Field2D indexDDY(const Field2D &f, CELL_LOC outloc)

F_INDEX_FIRSTDERIV(indexDDZ, defaultIndexDDZ, Field3D); // const Field3D indexDDZ(const Field3D &f, CELL_LOC outloc)
F_INDEX_FIRSTDERIV(indexDDZ, defaultIndexDDZ, Field2D); // const Field2D indexDDZ(const Field2D &f, CELL_LOC outloc)
