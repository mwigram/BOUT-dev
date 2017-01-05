
#include <globals.hxx>

#include <bout/difops/div_par.hxx>
#include <msg_stack.hxx>
#include <boutexception.hxx>

////////////////////////////////////////////////////////////////////
// Field2D objects

const Field2D Div_par(const Field2D &f, Difop method) {
  TRACE("Div_par(Field2D)");
  
  switch(method) {
  case Difop::C2: // No difference between this and FA, since no yup/down fields
  case Difop::C2_FA: return Div_par_C2(f); // Second order central difference
  case Difop::C4:
  case Difop::C4_FA: return Div_par_C4(f); // Fourth order central difference
  }
  
  throw BoutException("Method not supported");
}

const Field2D Div_par_C2(const Field2D &f) {
  TRACE("Div_par_C2(Field2D)");

  Field2D result;
  result.allocate();

  Coordinates *coord = mesh->coordinates();
  
  for(auto i : f.region(RGN_NOBNDRY)) {
    auto yp = i.yp();
    auto ym = i.ym();
    result[i] = coord->Bxy[i]*( (f[yp]/coord->Bxy[yp]) - (f[ym]/coord->Bxy[ym]) ) / (0.5*(coord->dy[yp]*sqrt(coord->g_22[yp]) + coord->dy[ym]*sqrt(coord->g_22[ym])) + coord->dy[i]*sqrt(coord->g_22[i]));
  }
  return result;
}

const Field2D Div_par_C4(const Field2D &f) {
  TRACE("Div_par_C4(Field2D)");

  Field2D result;
  result.allocate();

  Coordinates *coord = mesh->coordinates();
  
  for(auto i : f.region(RGN_NOBNDRY)) {
    auto yp = i.yp();
    auto ym = i.ym();
    auto ypp = i.offset(0,2,0);
    auto ymm = i.offset(0,-2,0);
    
    result[i] = coord->Bxy[i]*( (8.*f[yp]/coord->Bxy[yp]) - (8.*f[ym]/coord->Bxy[ym]) + (f[ymm]/coord->Bxy[ymm]) - (f[ypp]/coord->Bxy[ypp]) ) / (6.*(0.5*(coord->dy[yp]*sqrt(coord->g_22[yp]) + coord->dy[ym]*sqrt(coord->g_22[ym])) + coord->dy[i]*sqrt(coord->g_22[i])));
  }
  return result;
}

////////////////////////////////////////////////////////////////////
// Field3D objects
// Need to account for different cell locations, and yup/ydown 

const Field3D Div_par(const Field3D &f, Difop method) {
  
  switch(f.getLocation()) {
  case CELL_CENTRE: {
    // Cell centre -> Cell centre
    switch(method) {
    case Difop::C2:    return Div_par_C2(f);
    case Difop::C2_FA: return Div_par_C2_FA(f);
    case Difop::C4_FA: return Div_par_C4_FA(f);
    }
    break;
  }
  case CELL_YLOW: {
    // Cell ylow -> Cell centre
    switch(method) {
    case Difop::C2: return Div_par_LtoC_C2(f);
    case Difop::C2_FA: return Div_par_LtoC_C2_FA(f);
    }
    break;
  }
  }
  throw BoutException("Method not supported");
}

const Field3D Div_par_C2(const Field3D &f) {
  TRACE("Div_par_C2(Field3D)");
  ASSERT1(f.getLocation() == CELL_CENTRE);

  Field3D result;
  result.allocate();
  result.setLocation(CELL_CENTRE);

  Field3D fyup = f.yup();
  Field3D fydown = f.ydown();
  Coordinates *coord = mesh->coordinates();
  
  for(auto i : f.region(RGN_NOBNDRY)) {
    auto yp = i.yp();
    auto ym = i.ym();
    result[i] = coord->Bxy[i]*( (fyup[yp]/coord->Bxy[yp]) - (fydown[ym]/coord->Bxy[ym]) ) / (0.5*(coord->dy[yp]*sqrt(coord->g_22[yp]) + coord->dy[ym]*sqrt(coord->g_22[ym])) + coord->dy[i]*sqrt(coord->g_22[i]));
  }
  return result;
}

const Field3D Div_par_LtoC_C2(const Field3D &f) {
  TRACE("Div_par_LtoC_C2(Field3D)");
  ASSERT1(f.getLocation() == CELL_YLOW);

  Field3D result;
  result.allocate();
  result.setLocation(CELL_CENTRE);

  Field3D fyup = f.yup();
  Coordinates *coord = mesh->coordinates();
  
  for(auto i : f.region(RGN_NOBNDRY)) {
    auto yp = i.yp();
    auto ym = i.ym();
   
    // f[i] is at the lower side of the cell result[i]
    // fyup[yp] is at the upper side of the cell

    // Calculate B at upper and lower side of cell
    BoutReal Bup = 0.5*(coord->Bxy[yp] + coord->Bxy[i]);
    BoutReal Bdown = 0.5*(coord->Bxy[ym] + coord->Bxy[i]);

    result[i] = coord->Bxy[i]*( fyup[yp]/Bup - f[i]/Bdown) / (coord->dy[i]*sqrt(coord->g_22[i]));
  }
  return result;
}

const Field3D Div_par_LtoC_C2_FA(const Field3D &f) {
  TRACE("Div_par_LtoC_C2_FA(Field3D)");
  ASSERT1(f.getLocation() == CELL_YLOW);

  Field3D result;
  result.allocate();
  result.setLocation(CELL_CENTRE);
  
  // Transform to field aligned coordinates
  Field3D fa = mesh->toFieldAligned(f);
  Coordinates *coord = mesh->coordinates();
  
  for(auto i : f.region(RGN_NOBNDRY)) {
    auto yp = i.yp();
    auto ym = i.ym();
   
    // fa[i] is at the lower side of the cell result[i]
    // fa[yp] is at the upper side of the cell

    // Calculate B at upper and lower side of cell
    BoutReal Bup = 0.5*(coord->Bxy[yp] + coord->Bxy[i]);
    BoutReal Bdown = 0.5*(coord->Bxy[ym] + coord->Bxy[i]);

    result[i] = coord->Bxy[i]*( fa[yp]/Bup - fa[i]/Bdown) / (coord->dy[i]*sqrt(coord->g_22[i]));
  }
  // Transform result back to shifted coordinates
  return mesh->fromFieldAligned(result);
}

const Field3D Div_par_C2_FA(const Field3D &f) {
  TRACE("Div_par_C2_FA(Field3D)");
  Field3D result;
  result.allocate();

  // Transform to field aligned coordinates
  Field3D fa = mesh->toFieldAligned(f);
  Coordinates *coord = mesh->coordinates();
  
  for(auto i : f.region(RGN_NOBNDRY)) {
    auto yp = i.yp();
    auto ym = i.ym();
    result[i] = coord->Bxy[i]*( (fa[yp]/coord->Bxy[yp]) - (fa[ym]/coord->Bxy[ym]) ) / (0.5*(coord->dy[yp]*sqrt(coord->g_22[yp]) + coord->dy[ym]*sqrt(coord->g_22[ym])) + coord->dy[i]*sqrt(coord->g_22[i]));
  }
  
  // Transform result back to shifted coordinates
  return mesh->fromFieldAligned(result);
}

const Field3D Div_par_C4_FA(const Field3D &f) {
  TRACE("Div_par_C4_FA(Field3D)");
  Field3D result;
  result.allocate();

  // Transform to field aligned coordinates
  Field3D fa = mesh->toFieldAligned(f);
  Coordinates *coord = mesh->coordinates();
  
  for(auto i : f.region(RGN_NOBNDRY)) {
    auto yp = i.yp();
    auto ym = i.ym();
    auto ypp = i.offset(0,2,0);
    auto ymm = i.offset(0,-2,0);
    
    result[i] = coord->Bxy[i]*( (8.*fa[yp]/coord->Bxy[yp]) - (8.*fa[ym]/coord->Bxy[ym]) + (fa[ymm]/coord->Bxy[ymm]) - (fa[ypp]/coord->Bxy[ypp]) ) / (6.*(0.5*(coord->dy[yp]*sqrt(coord->g_22[yp]) + coord->dy[ym]*sqrt(coord->g_22[ym])) + coord->dy[i]*sqrt(coord->g_22[i])));
  }
  
  // Transform result back to shifted coordinates
  return mesh->fromFieldAligned(result);
}
