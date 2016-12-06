/*******************************************************************
 * Compressible gas-dynamics
 *
 * Changelog
 * ---------
 *
 * November 2016: Ben Dudson <benjamin.dudson@york.ac.uk>
 *   o Updating to PhysicsModel
 *   o Merging R-T and sod-shock examples
 * 
 * December 2007: Ben Dudson <benjamin.dudson@york.ac.uk>
 *   o Initial version
 * 
 *******************************************************************/

#include <bout/physicsmodel.hxx>

/*!
 * Compressible gas dynamics
 */
class GasCompress : public PhysicsModel {
private:

  // Evolving variables 
  Field3D N, P; ///< Density, Pressure
  Vector3D V;   ///< velocity

  // parameters
  BoutReal gamma_ratio;    ///< Ratio of specific heats
  BoutReal nu;             ///< Viscosity
  bool include_viscosity;  ///< Should viscosity be included?
  bool sub_initial;        ///< Subtract forces due to initial profiles?
  
  /// 2D initial profiles
  Field2D N0, P0;
  Vector2D V0;
  
  /// Acceleration due to gravity
  Vector2D g; 
  
protected:
  /*!
   * Initialise the model
   *
   * @param[in] restarting  Is this a restart?
   */
  int init(bool restarting) {
    BoutReal v0_multiply;

    // Read initial conditions
    
    mesh->get(N0, "density");
    mesh->get(P0, "pressure");
    V0.covariant = false; // Read contravariant components
    V.covariant = false; // Evolve contravariant components
    mesh->get(V0, "v");
    g.covariant = false;
    mesh->get(g, "g");
    
    // read options
    
    Options *options = Options::getRoot();
    options = options->getSection("gas");
    
    options->get("gamma",  gamma_ratio, 5.0/3.0);
    options->get("viscosity", nu, 0.1);
    options->get("include_viscosity", include_viscosity, false);
    options->get("v0_multiply", v0_multiply, 1.0);
    options->get("sub_initial", sub_initial, false);
  
    V0 *= v0_multiply;
    
    V.y.setLocation(CELL_YLOW); // Stagger
    
    // Set evolving variables
    SOLVE_FOR3(N, P, V);
    
    if(!restarting) {
      // Set variables to these values (+ the initial perturbation)
      // NOTE: This must be after the calls to bout_solve
      N += N0;
      P += P0;
      V += V0;
    }
    
    return 0;
  }

  /*!
   * Calculate time-derivatives ddt(N), ddt(P), ddt(V) given
   * evolving variables N, P and V
   */
  int rhs(BoutReal t) {
    // Communicate variables
    mesh->communicate(N,P,V);

    // Density
    
    ddt(N) = -V_dot_Grad(V, N) - N*Div(V);
    //ddt(N) = -Div(V, N);
    
    // Velocity 
    
    ddt(V) = -V_dot_Grad(V, V) - Grad(P)/N + g;
    //ddt(V) = -V_dot_Grad(V, V) - Grad(P, CELL_DEFAULT, CELL_YLOW)/N + g;
    
    if(sub_initial) {
      ddt(V) += Grad(P0)/N0 - g;
    }
    
    if(include_viscosity) {
      // Add viscosity
      
      ddt(V).y += nu*Laplace(V.y);
      ddt(V).z += nu*Laplace(V.z);
    }
    
    // Pressure
    
    ddt(P) = -V_dot_Grad(V, P) - gamma_ratio*P*Div(V);
    //ddt(P) = -Div(V, P) - (gamma_ratio-1.)*P*Div(V);
    
    return 0;
  }
};
  
/*!
 * Create a default main() function
 */
BOUTMAIN(GasCompress);
