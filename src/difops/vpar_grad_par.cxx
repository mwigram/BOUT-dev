
const Field2D Vpar_Grad_par(const Field2D &v, const Field2D &f, Difop method) {
  TRACE("Vpar_Grad_par(Field2D, Field2D)");
  
  switch(method) {
  case Difop::U1:
  case Difop::U1_FA: return Vpar_Grad_par_U1(v,f);
  case Difop::U2:
  case Difop::U2_FA: return Vpar_Grad_par_U2(v,f);
  case Difop::C2:
  case Difop::C2_FA: return Vpar_Grad_par_C2(v,f);
  }
  throw BoutException("method not supported");
}

const Field2D Vpar_Grad_par_U1(const Field2D &v, const Field2D &f) {
  TRACE("Vpar_Grad_par_U1(Field2D, Field2D)");
  
  Field2D result;
  result.allocate();
  
  Coordinates *coord = mesh->coordinates();
  
  for(auto i : f.region(RGN_NOBNDRY)) {
    if(v[i] > 0.0) {
      result[i] = v[i] * ( f[i] - f[i.ym()] ) / (coord->dy[i] * sqrt(coord->g_22[i]));
    }else {
      result[i] = v[i] * ( f[i.yp()] - f[i] ) / (coord->dy[i] * sqrt(coord->g_22[i]));
    }
  }
  return result;
}

const Field2D Vpar_Grad_par_U2(const Field2D &v, const Field2D &f) {
  TRACE("Vpar_Grad_par_U2(Field2D, Field2D)");
  
  Field2D result;
  result.allocate();
  
  Coordinates *coord = mesh->coordinates();
  
  for(auto i : f.region(RGN_NOBNDRY)) {
    if(v[i] > 0.0) {
      result[i] = v[i] * ( 1.5*f[i] - 2.*f[i.ym()] + 0.5*f[i.offset(0,-2,0)] ) / (coord->dy[i] * sqrt(coord->g_22[i]));
    }else {
      result[i] = v[i] * ( -0.5*f[i.offset(0,2,0)] + 2.*f[i.yp()] - 1.5*f[i] ) / (coord->dy[i] * sqrt(coord->g_22[i]));
    }
  }
  return result;
}

const Field2D Vpar_Grad_par_C2(const Field2D &v, const Field2D &f) {
  TRACE("Vpar_Grad_par_C2(Field2D, Field2D)");
  
  Field2D result;
  result.allocate();
  
  Coordinates *coord = mesh->coordinates();
  
  for(auto i : f.region(RGN_NOBNDRY)) {
    result[i] = v[i] * ( f[i.yp()] - f[i.ym()] ) / (2.*coord->dy[i] * sqrt(coord->g_22[i]));
  }
  return result;
}

////////////////////////////////////////////////////////////////////
// Field3D objects
// Need to account for different cell locations, and yup/ydown 

const Field3D Vpar_Grad_par(const Field3D &v, const Field3D &f, Difop method) {
  TRACE("Vpar_Grad_par(Field3D, Field3D)");
  
  switch(v.getLocation()) {
    case CELL_CENTRE: {
      // Assume all fields cell centred. These assumptions
      // should be checked in the individual implementations
      
      switch(method) {
      case Difop::U1: return Vpar_Grad_par_U1(v, f);
      case Difop::U1_FA: return Vpar_Grad_par_U1_FA(v, f);
      case Difop::C2: return Vpar_Grad_par_C2(v, f);
      case Difop::C2_FA: return Vpar_Grad_par_C2_FA(v, f);
      case Difop::W3_FA: return Vpar_Grad_par_W3_FA(v, f);
      }
      break;
    }
  case CELL_YLOW: {
    // Assume velocity at YLOW, f cell centre
    
    break;
  }
  }
  throw BoutException("method not supported");
}

const Field3D Vpar_Grad_par_U1(const Field3D &v, const Field3D &f) {
  TRACE("Vpar_Grad_par_U1(Field3D, Field3D)");
  ASSERT1(v.getLocation() == CELL_CENTRE);
  ASSERT1(f.getLocation() == CELL_CENTRE);
  
  Field3D result;
  result.allocate();
  result.setLocation(CELL_CENTRE);
  
  Field3D &fup = f.yup();
  Field3D &fdown = f.ydown();
  
  Coordinates *coord = mesh->coordinates();
  
  for(auto i : f.region(RGN_NOBNDRY)) {
    if(v[i] > 0.0) {
      result[i] = v[i] * ( f[i] - fdown[i.ym()] ) / (coord->dy[i] * sqrt(coord->g_22[i]));
    }else {
      result[i] = v[i] * ( fup[i.yp()] - f[i] ) / (coord->dy[i] * sqrt(coord->g_22[i]));
    }
  }
  return result;
}

const Field3D Vpar_Grad_par_C2(const Field3D &v, const Field3D &f) {
  TRACE("Vpar_Grad_par_C2(Field3D, Field3D)");
  ASSERT1(v.getLocation() == CELL_CENTRE);
  ASSERT1(f.getLocation() == CELL_CENTRE);
  
  Field3D result;
  result.allocate();
  result.setLocation(CELL_CENTRE);
  
  Field3D &fup = f.yup();
  Field3D &fdown = f.ydown();
  
  Coordinates *coord = mesh->coordinates();
  
  for(auto i : f.region(RGN_NOBNDRY)) {
    result[i] = v[i] * ( fup[i.yp()] - fdown[i.ym()] ) / (2.*coord->dy[i] * sqrt(coord->g_22[i]));
  }
  return result;
}

const Field3D Vpar_Grad_par_U1_FA(const Field3D &v, const Field3D &f) {
  TRACE("Vpar_Grad_par_U1_FA(Field3D, Field3D)");
  ASSERT1(v.getLocation() == CELL_CENTRE);
  ASSERT1(f.getLocation() == CELL_CENTRE);
  
  Field3D result;
  result.allocate();
  result.setLocation(CELL_CENTRE);
  
  // Transform to field aligned coordinates
  Field3D va = mesh->toFieldAligned(v);
  Field3D fa = mesh->toFieldAligned(f);

  Coordinates *coord = mesh->coordinates();
  
  for(auto i : f.region(RGN_NOBNDRY)) {
    if(v[i] > 0.0) {
      result[i] = va[i] * ( fa[i] - fa[i.ym()] ) / (coord->dy[i] * sqrt(coord->g_22[i]));
    }else {
      result[i] = va[i] * ( fa[i.yp()] - fa[i] ) / (coord->dy[i] * sqrt(coord->g_22[i]));
    }
  }
  
  // Transform result back
  return mesh->fromFieldAligned(result);
}

const Field3D Vpar_Grad_par_C2_FA(const Field3D &v, const Field3D &f) {
  TRACE("Vpar_Grad_par_C2_FA(Field3D, Field3D)");
  ASSERT1(v.getLocation() == CELL_CENTRE);
  ASSERT1(f.getLocation() == CELL_CENTRE);
  
  Field3D result;
  result.allocate();
  result.setLocation(CELL_CENTRE);
  
  // Transform to field aligned coordinates
  Field3D va = mesh->toFieldAligned(v);
  Field3D fa = mesh->toFieldAligned(f);

  Coordinates *coord = mesh->coordinates();
  
  for(auto i : f.region(RGN_NOBNDRY)) {
    result[i] = va[i] * ( fa[i.yp()] - fa[i.ym()] ) / (2.*coord->dy[i] * sqrt(coord->g_22[i]));
  }
  
  // Transform result back
  return mesh->fromFieldAligned(result);
}

const Field3D Vpar_Grad_par_W3_FA(const Field3D &v, const Field3D &f) {
  TRACE("Vpar_Grad_par_W3_FA(Field3D, Field3D)");
  ASSERT1(v.getLocation() == CELL_CENTRE);
  ASSERT1(f.getLocation() == CELL_CENTRE);

  const BoutReal WENO_SMALL = 1.0e-8; // Small number for WENO schemes

  Field3D result;
  result.allocate();
  result.setLocation(CELL_CENTRE);
  
  // Transform to field aligned coordinates
  Field3D va = mesh->toFieldAligned(v);
  Field3D fa = mesh->toFieldAligned(f);

  Coordinates *coord = mesh->coordinates();
  
  for(auto i : f.region(RGN_NOBNDRY)) {
    auto ym = i.ym();
    auto yp = i.yp();
    
    if(v[i] > 0.0) {
      // Left-biased stencil
      auto ymm = i.offset(0,-2,0);
      
      BoutReal r = (WENO_SMALL + SQ(fa[i] - 2.0*fa[ym] + fa[ymm])) / (WENO_SMALL + SQ(fa[yp] - 2.0*fa[i] + fa[ym]));
      BoutReal w = 1.0 / (1.0 + 2.0*r*r);
      
      BoutReal deriv = 0.5*(fa[yp] - fa[ym]) - 0.5*w*(-fa[ymm] + 3.*fa[ym] - 3.*fa[i] + fa[yp]);
      
      result[i] = va[i] * deriv / (coord->dy[i] * sqrt(coord->g_22[i]));
    }else {
      // Right-biased stencil
      BoutReal r = (WENO_SMALL + SQ(fa[ypp] - 2.0*fa[yp] + f[i])) / (WENO_SMALL + SQ(fa[yp] - 2.0*fa[i] + fa[ym]));
      BoutReal w = 1.0 / (1.0 + 2.0*r*r);
      
      BoutReal deriv = 0.5*(fa[yp] - fa[ym]) - 0.5*w*( -fa[ym] + 3.*fa[i] - 3.*fa[yp] + fa[ypp] );
      
      result[i] = va[i] * deriv / (coord->dy[i] * sqrt(coord->g_22[i]));
    }
  }
  
  // Transform result back
  return mesh->fromFieldAligned(result);
}


