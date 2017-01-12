
#include <globals.hxx>

#include <bout/difops/ddy.hxx>

#include <msg_stack.hxx>
#include <boutexception.hxx>

////////////////////////////////////////////////////////////////////


/// Partial derivative in Y
///
/// @param[in] f  The field to be differentiated
const Field2D DDY(const Field2D &f, Difop method) {
  TRACE("DDY(Field2D)");
  
  switch(method) {
  case Difop::C2: // No difference between this and FA, since no yup/down fields
  case Difop::C2_FA: return DDY_C2(f);
  case Difop::C4:
  case Difop::C4_FA: return DDY_C4(f);
  }
  
  throw BoutException("Method not supported");
}

const Field2D DDY_C2(const Field2D &f) {
  TRACE("DDY_C2(Field2D)");
  
  Field2D result;
  result.allocate();
  
  Coordinates *coord = mesh->coordinates();

  for(auto i : f.region(RGN_NOBNDRY)) {
    result[i] = (f[i.yp()] - f[i.ym()]) / (2.*coord->dy[i]);
  }
  return result;
}

const Field2D DDY_C4(const Field2D &f) {
  TRACE("DDY_C4(Field2D)");
  
  Field2D result;
  result.allocate();
  
  Coordinates *coord = mesh->coordinates();

  for(auto i : f.region(RGN_NOBNDRY)) {
    result[i] = (8.*f[i.yp()] - 8.*f[i.ym()] + f[i.offset(0,-2,0)] - f[i.offset(0,2,0)]) / (12.*coord->dy[i]);
  }
  return result;
}

////////////////////////////////////////////////////////////////////
// Field3D operators

const Field3D DDY(const Field3D &f, CELL_LOC outloc, Difop method) {
  
  if(outloc == CELL_DEFAULT) {
    // Output is at the same location as input
    
    switch(method) {
    case Difop::C2:    return DDY_C2(f);
    case Difop::C2_FA: return DDY_C2_FA(f);
    }
  }
  
  switch(f.getLocation()) {
    
  case CELL_YLOW: {
    // Cell ylow -> Cell centre
    
    break;
  }
  }
  throw BoutException("Method not supported");
}

/// Stencil class S
/// 
template<class F, class S, class SD, class Op>
const F applyDeriv(const F &f, const SD &spacing, REGION region) {
  F result;
  result.allocate();
  
  S stencil(f); // Create a stencil

  for(auto i : result.region(RGN_NOBNDRY)) {
    // Update the stencil values (c,p,m)
    stencil.update(i);
    // Update the cell spacing values
    spacing.update(i);
    // Perform operation
    result[i] = Op::op(stencil, spacing);
  }
  
  // Final transformation e.g. from Field Aligned
  return stencil.final(result);
}

const Field3D DDY_C2(const Field3D &f) {
  StencilY2D dy(mesh->coordinates()->dy); // Mesh spacing
  return applyDeriv<Field3D, StencilY3D, StencilY2D, Op_DDX_C2>(f, dy);
}

const Field3D DDY_C2_FA(const Field3D &f) {
  return applyDeriv<Field3D, StencilY3D_FA, Op_DDX_C2>(f);
}



#endif //  __DDY_H__
