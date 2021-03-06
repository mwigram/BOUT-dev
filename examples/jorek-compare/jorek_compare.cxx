/**************************************************************************
 * Similar set of equations to JOREK
 * 
 **************************************************************************/

#include <bout.hxx>
#include <boutmain.hxx>

#include <invert_laplace.hxx>
#include <bout/constants.hxx>
#include <math.h>

// Evolving quantities
Field3D rho, Te, Ti, U, Vpar, Apar;
// Derived quantities
Field3D Jpar, phi; // Parallel current, electric potential

// Equilibrium quantities
Field2D rho0, Te0, Ti0; // Equilibrium mass density, electron and ion temperature
Field2D B0, J0, P0;
Vector2D b0xcv; // Curvature term
Vector2D B0vec; // B0 field vector

// Dissipation coefficients
Field2D D_perp; // Particle diffusion coefficient
Field2D chi_eperp, chi_epar; // Electron heat diffusion coefficients
Field2D chi_iperp, chi_ipar; // Ion heat diffusion coefficients

// Collisional terms
BoutReal tau_enorm;
Field3D tau_e; // electron collision time

Field2D eta0;  // Resistivity
Field3D eta;

BoutReal viscos_par, viscos_perp, viscos_coll; // Viscosity coefficients
BoutReal hyperresist; // Hyper-resistivity coefficient

int phi_flags;

// Constants
const BoutReal MU0 = 4.0e-7*PI;
const BoutReal Charge = 1.60217646e-19; // electron charge e (C)
const BoutReal Mi = 2.0*1.67262158e-27; // Ion mass
const BoutReal Me = 9.1093816e-31;  // Electron mass
const BoutReal Me_Mi = Me / Mi; // Electron mass / Ion mass

// Normalisation factors
BoutReal Tnorm, rhonorm; // Partial normalisation to rho and MU0. Temperature normalised

// options

bool nonlinear;
bool full_bfield;   // If true, use divergence-free expression for B
bool flux_method;   // Use flux methods in rho and T equations
int  jpar_bndry_width; // Set jpar = 0 in a boundary region

bool electron_density; // Solve Ne rather than Ni (adds Jpar term to density)

bool vorticity_momentum; // Vorticity is curl of momentum, rather than velocity

bool include_profiles; // Include zero-order equilibrium terms

bool parallel_lc; // Use CtoL and LtoC differencing

int low_pass_z; // Toroidal (Z) filtering of all variables

Vector3D vExB, vD; // Velocities
Field3D divExB;    // Divergence of ExB flow

BoutReal Wei;  // Factor for the electron-ion collision term
bool ohmic_heating;

// Poisson brackets: b0 x Grad(f) dot Grad(g) / B = [f, g]
// Method to use: BRACKET_ARAKAWA, BRACKET_STD or BRACKET_SIMPLE
BRACKET_METHOD bm; // Bracket method for advection terms

// Communication objects
FieldGroup comms;

int physics_init(bool restarting) {
  
  output.write("Solving JOREK-like reduced MHD equations\n");
  output.write("\tFile    : %s\n", __FILE__);
  output.write("\tCompiled: %s at %s\n", __DATE__, __TIME__);

  Options *globalOptions = Options::getRoot();
  Options *options = globalOptions->getSection("jorek");

  //////////////////////////////////////////////////////////////
  // Load data from the grid

  // Load 2D profiles
  mesh->get(J0, "Jpar0");    // A / m^2
  
  if(mesh->get(rho0, "Ni0")) {
    output << "Warning: No density profile available\n";
    BoutReal d0;
    options->get("density", d0, 1.0);
    rho0 = d0;
  }
  rho0 *= 1e20; // Convert to m^[-3]
  
  // Read temperature
  mesh->get(Te0, "Te0");
  mesh->get(Ti0, "Ti0");

  // Try reading pressure profile (in Pascals)
  if(mesh->get(P0, "pressure")) {
    // Just calculate from Temp and density
    P0 = Charge * (Ti0 + Te0) * rho0;
  }else {
    // Make sure that density and temperature are consistent with pressure
    
    Field2D factor = P0 / (Charge * (Ti0 + Te0) * rho0);
    
    output.write("\tPressure factor %e -> %e\n", min(factor,true), max(factor, true));
    
    // Multiply temperatures by this factor
    Te0 *= factor;
    Ti0 *= factor;
  }
  rho0 *= Mi; // Convert density to mass density [kg / m^3]

  // Load dissipation coefficients, override in options file
  if(options->isSet("D_perp")) {
    BoutReal tmp;
    options->get("D_perp", tmp, 0.0);
    D_perp = tmp;
  }else mesh->get(D_perp, "D_perp");

  if(options->isSet("chi_eperp")) {
    BoutReal tmp;
    options->get("chi_eperp", tmp, 0.0);
    chi_eperp = tmp;
  }else mesh->get(chi_eperp, "chi_eperp");

  if(options->isSet("chi_iperp")) {
    BoutReal tmp;
    options->get("chi_iperp", tmp, 0.0);
    chi_iperp = tmp;
  }else mesh->get(chi_iperp, "chi_iperp");

  if(options->isSet("chi_epar")) {
    BoutReal tmp;
    options->get("chi_epar", tmp, 0.0);
    chi_epar = tmp;
  }else mesh->get(chi_epar, "chi_epar");

  if(options->isSet("chi_ipar")) {
    BoutReal tmp;
    options->get("chi_ipar", tmp, 0.0);
    chi_ipar = tmp;
  }else mesh->get(chi_ipar, "chi_ipar");

  if(options->isSet("viscos_perp")) {
    BoutReal tmp;
    options->get("viscos_perp", tmp, -1.0);
    viscos_perp = tmp;
  }else mesh->get(viscos_perp, "viscos_perp");

  if(options->isSet("viscos_par")) {
    BoutReal tmp;
    options->get("viscos_par", tmp, -1.0);
    viscos_par = tmp;
  }else mesh->get(viscos_par, "viscos_par");

  OPTION(options, viscos_coll, -1.0);

  // Load curvature term
  b0xcv.covariant = false; // Read contravariant components
  mesh->get(b0xcv, "bxcv"); // mixed units x: T y: m^-2 z: m^-2
  
  // Metric coefficients
  Field2D Rxy, Bpxy, Btxy, hthe;
  Field2D I; // Shear factor
  
  if(mesh->get(Rxy,  "Rxy")) { // m
    output.write("Error: Cannot read Rxy from grid\n");
    return 1;
  }
  if(mesh->get(Bpxy, "Bpxy")) { // T
    output.write("Error: Cannot read Bpxy from grid\n");
    return 1;
  }
  mesh->get(Btxy, "Btxy"); // T
  mesh->get(B0,   "Bxy");  // T
  mesh->get(hthe, "hthe"); // m
  mesh->get(I,    "sinty");// m^-2 T^-1
  
  OPTION(options, nonlinear,           false);
  OPTION(options, full_bfield,         false);
  OPTION(options, flux_method,         false);
  
  OPTION(options, jpar_bndry_width,    -1);

  OPTION(options, hyperresist,         -1);

  OPTION(options, electron_density,    false);
  OPTION(options, vorticity_momentum,  false);
  OPTION(options, include_profiles,    false);
  OPTION(options, parallel_lc,         true);

  OPTION(options, phi_flags,  0);

  OPTION(options, low_pass_z, -1); // Default is no filtering

  OPTION(options, Wei, 1.0);
  
  OPTION(options, ohmic_heating, true);

  int bracket_method;
  OPTION(options, bracket_method, 0);
  switch(bracket_method) {
  case 0: {
    bm = BRACKET_STD; 
    output << "\tBrackets: default differencing\n";
    break;
  }
  case 1: {
    bm = BRACKET_SIMPLE; 
    output << "\tBrackets: simplified operator\n";
    break;
  }
  case 2: {
    bm = BRACKET_ARAKAWA; 
    output << "\tBrackets: Arakawa scheme\n";
    break;
  }
  case 3: {
    bm = BRACKET_CTU; 
    output << "\tBrackets: Corner Transport Upwind method\n";
    break;
  }
  default:
    output << "ERROR: Invalid choice of bracket method. Must be 0 - 3\n";
    return 1;
  }

  //////////////////////////////////////////////////////////////
  // SHIFTED RADIAL COORDINATES

  if(mesh->ShiftXderivs) {
    if(mesh->IncIntShear) {
      // BOUT-06 style, using d/dx = d/dpsi + I * d/dz
      mesh->IntShiftTorsion = I;
      
    }else {
      // Dimits style, using local coordinate system
      b0xcv.z += I*b0xcv.x;
      I = 0.0;  // I disappears from metric
    }
  }

  //////////////////////////////////////////////////////////////
  // NORMALISE QUANTITIES
  
  rhonorm = max(rho0, true); // Maximum over all grid
  BoutReal Temax = max(Te0,true); // Maximum Te value
  Tnorm = Mi / (MU0 * Charge * rhonorm); // Temperature normalisation

  SAVE_ONCE2(rhonorm, Tnorm); // Save normalisation factors to file

  // Normalise quantities
  
  P0 *= MU0;
  J0 *= MU0;
  rho0 /= rhonorm;
  Te0 /= Tnorm;
  Ti0 /= Tnorm;
  
  viscos_perp *= sqrt(MU0 / rhonorm);
  viscos_par  *= sqrt(MU0 / rhonorm);
  D_perp      *= sqrt(MU0 * rhonorm);
  chi_eperp   *= sqrt(MU0 / rhonorm);
  chi_epar    *= sqrt(MU0 / rhonorm);
  chi_iperp   *= sqrt(MU0 / rhonorm);
  chi_ipar    *= sqrt(MU0 / rhonorm);
  
  // Coulomb logarithm
  BoutReal CoulombLog = 6.6 - 0.5*log(rhonorm/(Mi*1e20)) + 1.5*log(Temax);
  output << "\tCoulomb logarithm = " << CoulombLog << endl;
  
  // Factor in front of tau_e expression
  // tau_e = tau_enorm * Tet^1.5 / rhot
  tau_enorm = 3.44e11*(Mi/rhonorm)*Tnorm*sqrt(Tnorm) / CoulombLog;
  output << "\ttau_enorm = " << tau_enorm ;
  tau_enorm /= sqrt(MU0*rhonorm); // Normalise
  output << "\tNormalised tau_enorm = " << tau_enorm << endl;
  
  // Calculate or read in the resistivity
  if(options->isSet("eta")) {
    BoutReal etafactor;
    options->get("eta", etafactor, 0.0);
    // Calculate in normalised units
    eta0 = etafactor * Me*Mi/(1.96*MU0*rhonorm*Charge*Charge*tau_enorm*rho0);
  }else {
    mesh->get(eta0, "eta0");     // Read in SI units
    eta0 *= sqrt(rhonorm / MU0); // Normalise
  }


  //////////////////////////////////////////////////////////////
  // CALCULATE METRICS
  
  mesh->g11 = (Rxy*Bpxy)^2;
  mesh->g22 = 1.0 / (hthe^2);
  mesh->g33 = (I^2)*mesh->g11 + (B0^2)/mesh->g11;
  mesh->g12 = 0.0;
  mesh->g13 = -I*mesh->g11;
  mesh->g23 = -Btxy/(hthe*Bpxy*Rxy);
  
  mesh->J = hthe / Bpxy;
  mesh->Bxy = B0;
  
  mesh->g_11 = 1.0/mesh->g11 + ((I*Rxy)^2);
  mesh->g_22 = (B0*hthe/Bpxy)^2;
  mesh->g_33 = Rxy*Rxy;
  mesh->g_12 = Btxy*hthe*I*Rxy/Bpxy;
  mesh->g_13 = I*Rxy*Rxy;
  mesh->g_23 = Btxy*hthe*Rxy/Bpxy;
  
  mesh->geometry(); // Calculate quantities from metric tensor

  // Set B field vector
  B0vec.covariant = false;
  B0vec.x = 0.;
  B0vec.y = Bpxy / hthe;
  B0vec.z = 0.;
  
  vExB.setBoundary("v");
  vD.setBoundary("v");
  
  Jpar.setBoundary("Jpar");

  phi.setBoundary("phi");

  // Set starting dissipation terms
  eta = eta0;
  tau_e = tau_enorm * (Te0^1.5)/rho0;

  output.write("\tNormalised tau_e = %e -> %e\n", min(tau_e, true), max(tau_e, true));

  // Set locations for staggered grids
  vD.setLocation(CELL_VSHIFT);

  // SET EVOLVING VARIABLES

  SOLVE_FOR6(rho, Te, Ti, U, Vpar, Apar);
  
  comms.add(rho);
  comms.add(Te);
  comms.add(Ti);
  comms.add(U);
  comms.add(phi);
  comms.add(Vpar);
  comms.add(Apar);
  
  SAVE_ONCE5(P0, J0, rho0, Te0, Ti0); // Save normalised profiles

  if(nonlinear) {
    SAVE_REPEAT(eta);
  }else {
    SAVE_ONCE(eta);
  }

  SAVE_REPEAT2(phi, Jpar); // Save each timestep
  
  SAVE_REPEAT(divExB);

  return 0;
}

// Parallel gradient along perturbed field-line
const Field3D Grad_parP(const Field3D &f, CELL_LOC loc = CELL_DEFAULT) {
  // Derivative along equilibrium field-line
  Field3D result;
  
  if(parallel_lc) {
    if(loc == CELL_YLOW) {
      result = Grad_par_CtoL(f);
    }else
      result = Grad_par_LtoC(f);
  }else
    result = Grad_par(f, loc);
  
  if(nonlinear) {
    if(full_bfield) {
      // Use full expression for perturbed B
      Vector3D Btilde = Curl(B0vec * Apar / B0);
      result += Btilde * Grad(f) / B0;
    }else {
      // Simplified expression
      result -= bracket(Apar, f, BRACKET_ARAKAWA);
    }
  }
  return result;
}

const Field3D Div_parP(const Field3D &f, CELL_LOC loc = CELL_DEFAULT) {
  return B0*Grad_parP(f/B0, loc);
}

int physics_run(BoutReal t) {
  
  int sp = msg_stack.push("Started physics_run(%e)", t);
  
  // Invert laplacian for phi
  if(vorticity_momentum) {
    // Vorticity is b dot curl(rho * v)
    Field2D rprof = rho0;
    if(nonlinear)
      rprof += rho.DC(); // Axisymmetric rho only
    phi = invert_laplace(B0*U/rprof, phi_flags, NULL, &rprof);
  }else {
    // Vorticity is b dot curl(v)
    phi = invert_laplace(B0*U, phi_flags);
  }
  // Apply a boundary condition on phi for target plates
  phi.applyBoundary();
  
  // Communicate variables
  mesh->communicate(comms);
  

  // Get J from Psi
  Jpar = -Delp2(Apar);
  Jpar.applyBoundary();

  if(jpar_bndry_width > 0) {
    // Boundary in jpar
    if(mesh->firstX()) {
      for(int i=jpar_bndry_width;i>=0;i--)
        for(int j=0;j<mesh->LocalNy;j++)
          for(int k=0;k<mesh->LocalNz;k++) {
            Jpar[i][j][k] = 0.5*Jpar[i+1][j][k];
          }
    }
    if(mesh->lastX()) {
      for(int i=mesh->LocalNx-jpar_bndry_width-1;i<mesh->LocalNx;i++)
        for(int j=0;j<mesh->LocalNy;j++)
          for(int k=0;k<mesh->LocalNz;k++) {
            Jpar[i][j][k] = 0.5*Jpar[i-1][j][k];
          }
    }
  }
  
  mesh->communicate(Jpar);
  
  //Jpar = smooth_x(Jpar); // Smooth in x direction

  Field3D rhot = rho0;
  Field3D Tet = Te0;
  Field3D Tit = Ti0;
  Field3D P = rho*(Te0+Ti0) + (Te+Ti)*rho0; // Perturbed pressure
  
  if(nonlinear) {
    rhot += rho;
    Tet += Te;
    Tit += Ti;
    P += rho*(Te+Ti);
    
    eta = eta0*((Tet/Te0)^(-1.5)); // Update resistivity based on Te
    
    tau_e = tau_enorm*(Tet^1.5)/rhot; // Update electron collision rate
  }

  if(flux_method) {
    msg_stack.push("Flux vExB");
    // ExB velocity
    vExB = (B0vec ^ Grad_perp(phi))/(B0*B0);
    vExB.applyBoundary();
    msg_stack.pop();
      
    ////////// Density equation ////////////////
    
    msg_stack.push("Flux Density");
      
    // Diffusive flux (perpendicular)
    vD = -D_perp * Grad_perp(rho);
    vD.applyBoundary();
    
    ddt(rho) = -Div(vExB + vD, rhot);
    
    msg_stack.pop();

    ////////// Temperature equations ////////////
  
    msg_stack.push("Flux Te");

    vD = -chi_eperp * Grad_perp(Te) - Grad_par(Te, CELL_YLOW)*chi_epar*B0vec;
    vD.applyBoundary();
    
    ddt(Te) = 
      - b0xGrad_dot_Grad(phi, Tet)/B0
      - (2./3.)*Tet*Div(vExB)
      - Div(vD, Te)/rhot
      ;
    
    msg_stack.pop();
    
    msg_stack.push("Flux Ti");
    
    vD = -chi_iperp * Grad_perp(Ti) - Grad_par(Ti, CELL_YLOW)*chi_ipar*B0vec;
    vD.applyBoundary();
    
    ddt(Ti) = 
      - b0xGrad_dot_Grad(phi, Tit)/B0
      - (2./3.)*Tit*Div(vExB)
      - Div(vD, Ti)/rhot
      ;

    msg_stack.pop();
  }else {
    // Use analytic expressions, expand terms
    
    // Divergence of ExB velocity (neglecting parallel term)
    msg_stack.push("divExB");
    divExB = b0xcv*Grad(phi)/B0 - b0xGrad_dot_Grad(1./B0, phi);
    msg_stack.pop();
    
    msg_stack.push("density");
    ddt(rho) = 
      - bracket(phi, rhot, bm)      // ExB advection 
      - divExB*rhot                    // Divergence of ExB (compression)
      - Vpar_Grad_par(Vpar, rho)       // Parallel advection
      - rhot*Div_parP(Vpar, CELL_CENTRE) // Parallel compression
      + D_perp * Delp2(rho)            // Perpendicular diffusion
      ;
    
    if(electron_density) {
      // Using electron parallel velocity rather than ion
      ddt(rho) += (Mi/(Charge*sqrt(MU0*rhonorm)))* Div_parP(Jpar, CELL_CENTRE);
    }

    if(low_pass_z > 0)
      ddt(rho) = lowPass(ddt(rho), low_pass_z);   
 
    msg_stack.pop();
    
    msg_stack.push("Te");
    ddt(Te) = 
      -bracket(phi, Tet, bm) - Vpar_Grad_par(Vpar, Tet) // advection
      - (2./3.)*Tet*(divExB  + Div_parP(Vpar, CELL_CENTRE)) // Divergence of flow
      + Div_par_K_Grad_par(chi_epar, Te) / rhot    // Parallel diffusion
      + chi_eperp*Delp2(Te) / rhot   // Perpendicular diffusion
      ;
    
    if(ohmic_heating)
      ddt(Te) += (2./3)*eta*Jpar*Jpar / rhot; // Ohmic heating

    msg_stack.pop();
    
    msg_stack.push("Ti");
    ddt(Ti) = 
      - bracket(phi, Tit, bm) - Vpar_Grad_par(Vpar, Tit)
      - (2./3.)*Tit*(divExB + Div_parP(Vpar, CELL_CENTRE))
      + Div_par_K_Grad_par(chi_ipar, Ti) / rhot
      + chi_iperp*Delp2(Ti) / rhot
      ;
    msg_stack.pop();
    
    
    if(Wei > 0.0) {
      // electron-ion collision term
      // Calculate Wi * (2/3)/rho term. Wei is a scaling factor from options
      Field3D Tei = Wei * 2.*Me_Mi*(Te - Ti) / tau_e;
      
      ddt(Ti) += Tei;
      ddt(Te) -= Tei;
    }

    if(low_pass_z > 0) {
      ddt(Te) = lowPass(ddt(Te), low_pass_z);
      ddt(Ti) = lowPass(ddt(Ti), low_pass_z);
    }
  }
  
  ////////// Vorticity equation ////////////
    
  msg_stack.push("Vorticity");
  if(vorticity_momentum) {
    // Vorticity is b dot curl(rho * v)
    
    ddt(U) = 
      (B0^2)*Grad_parP(Jpar/B0, CELL_CENTRE)  // b0 dot J
      + 2.*b0xcv*Grad(P)
      - rhot*(divExB + Div_parP(Vpar, CELL_CENTRE)) * Delp2(phi)/B0   // drho/dt term
      ;

    // b dot J0
    if(full_bfield) {
      Vector3D Btilde = Curl(B0vec * Apar / B0);
      ddt(U) += B0*Btilde * Grad(J0/B0);
    }else {
      ddt(U) -= B0 * B0*bracket(Apar, J0/B0, BRACKET_ARAKAWA);
    }
    
    if(electron_density) {
      // drho/dt jpar term
      ddt(U) += (Mi/(Charge*sqrt(MU0*rhonorm)))* rhot*Div_parP(Jpar/rhot, CELL_CENTRE) * Delp2(phi)/B0;
    }
    
    if(include_profiles) {
      ddt(U) +=
        (B0^2)*Grad_par(J0/B0) // b0 dot J0
        + 2.*b0xcv*Grad(P0)
        ;
    }
    
    if(nonlinear) {
      ddt(U) -= bracket(phi, U, bm);    // Advection
      ddt(U) -= Vpar_Grad_par(Vpar, U);         // Parallel advection
    }

    // Viscosity terms
    if(viscos_par > 0.0)
      ddt(U) += viscos_par * Grad2_par2(U); // Parallel viscosity
    
    if(viscos_perp > 0.0)
      ddt(U) += viscos_perp * rhot*Delp2(U/rhot); // Perpendicular viscosity
    
  }else {
    // Vorticity is b dot curl(v)
    ddt(U) = (
              (B0^2)*Grad_parP(Jpar/B0, CELL_CENTRE)
              + 2.*b0xcv*Grad(P)  // curvature term
              ) / rhot;
  
    // b dot J0
    if(full_bfield) {
      Vector3D Btilde = Curl(B0vec * Apar / B0);
      ddt(U) += B0*Btilde * Grad(J0/B0) / rhot;
    }else {
      ddt(U) -= B0 * B0* bracket(Apar, J0/B0, BRACKET_ARAKAWA) / rhot;
    }
  
    if(include_profiles) {
      ddt(U) += (
                 (B0^2)*Grad_par(J0/B0) // b0 dot J0
                 + 2.*b0xcv*Grad(P0)
                 ) / rhot;
    }

    if(nonlinear) {
      ddt(U) -= bracket(phi, U);    // Advection
      ddt(U) -= Vpar_Grad_par(Vpar, U);         // Parallel advection
    }
  
    // Viscosity terms
    if(viscos_par > 0.0)
      ddt(U) += viscos_par * Grad2_par2(U) / rhot; // Parallel viscosity
    
    if(viscos_perp > 0.0)
      ddt(U) += viscos_perp * Delp2(U) / rhot;     // Perpendicular viscosity
    
    // Collisional viscosity
    if(viscos_coll > 0.0)
      ddt(U) += viscos_coll / MU0 * eta * Delp2(U) / rhot;
  }
  
  if(low_pass_z > 0)
    ddt(U) = lowPass(ddt(U), low_pass_z);
   
  msg_stack.pop();

  ////////// Parallel velocity equation ////////////
    
  msg_stack.push("Vpar");
  
  ddt(Vpar) = -Grad_parP(P + P0, CELL_YLOW);
  if(nonlinear) {
    ddt(Vpar) -= bracket(phi, Vpar); // Advection
    ddt(Vpar) -= Vpar_Grad_par(Vpar, Vpar); // Parallel advection
  }
  
  if(low_pass_z > 0)
    ddt(Vpar) = lowPass(ddt(Vpar), low_pass_z);
   
  msg_stack.pop();

  ////////// Magnetic potential equation ////////////
  
  msg_stack.push("Apar");
  ddt(Apar) = 
    -Grad_parP(phi, CELL_YLOW)
    - eta*Jpar;
  
  if(hyperresist > 0.0) {
    ddt(Apar) += eta*hyperresist * Delp2(Jpar);
  }

  if(low_pass_z > 0)
    ddt(Apar) = lowPass(ddt(Apar), low_pass_z);

  msg_stack.pop(sp);
  return 0;
}
