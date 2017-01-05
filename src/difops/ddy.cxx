
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

const Field3D DDY(const Field3D &f, Difop method) {
  switch(f.getLocation()) {
  case CELL_CENTRE: {
    // Cell centre -> Cell centre
    switch(method) {
    case Difop::C2:    return DDY_C2(f);
    case Difop::C2_FA: return DDY_C2_FA(f);
    }
    break;
  }
  case CELL_YLOW: {
    // Cell ylow -> Cell centre
    
    break;
  }
  }
  throw BoutException("Method not supported");
}

const Field3D DDY_C2(const Field3D &f) {
  TRACE("DDY_C2(Field3D)");
  
  ASSERT1(f.getLocation() == CELL_CENTRE);

  Field3D result;
  result.allocate();
  result.setLocation(CELL_CENTRE);
  
  Field3D& fup = f.yup();
  Field3D& fdown = f.ydown();
  
  Coordinates *coord = mesh->coordinates();

  for(auto i : f.region(RGN_NOBNDRY)) {
    result[i] = (fup[i.yp()] - fdown[i.ym()]) / (2.*coord->dy[i]);
  }
  return result;
}

const Field3D DDY_C2_FA(const Field3D &f) {
  TRACE("DDY_C2_FA(Field3D)");
  
  ASSERT1(f.getLocation() == CELL_CENTRE);

  Field3D result;
  result.allocate();
  result.setLocation(CELL_CENTRE);
  
  // Transform to field aligned coordinates
  Field3D fa = mesh->toFieldAligned(f);
  
  Coordinates *coord = mesh->coordinates();

  for(auto i : f.region(RGN_NOBNDRY)) {
    result[i] = (fa[i.yp()] - fa[i.ym()]) / (2.*coord->dy[i]);
  }
  // Transform back from field aligned coordinates
  return mesh->fromFieldAligned(result);
}


#endif //  __DDY_H__
