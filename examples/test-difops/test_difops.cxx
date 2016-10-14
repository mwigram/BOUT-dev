
#include <bout.hxx>

#include <field_factory.hxx>

// New advection operators
#include <bout/index_advect.hxx>

// Timing
#include <chrono>
typedef std::chrono::time_point<std::chrono::steady_clock> SteadyClock;
typedef std::chrono::duration<double> Duration;
using namespace std::chrono;

int main(int argc, char**argv) {
  BoutInitialise(argc, argv);
  
  /*
   * Create fields from input options
   */

  FieldFactory f(mesh);

  Field3D a = f.create3D("a");
  Field3D b = f.create3D("b");

  /*
   * Old operator
   */
  Field3D result = mesh->indexVDDX(a, b, CELL_DEFAULT, DIFF_DEFAULT);  // Run once without timing
  
  SteadyClock start1 = steady_clock::now();
  result = mesh->indexVDDX(a, b, CELL_DEFAULT, DIFF_DEFAULT);
  Duration elapsed1 = steady_clock::now() - start1;
  
  /*
   * New operator in namespace DIFOPS
   */
  
  Field3D result_new = DIFOPS::indexVDDX(a,b);

  SteadyClock start2 = steady_clock::now();
  result_new = DIFOPS::indexVDDX(a,b);  
  Duration elapsed2 = steady_clock::now() - start2;
  
  /*
   * Compare old and new results
   */
  bool pass = true;
  for(auto i : result.region(RGN_NOBNDRY)) {
    BoutReal diff = fabs(result_new[i] - result[i]);
    if(diff > 1e-8) {
      pass = false;
      output.write("(%d,%d,%d) : %e -> %e\n", i.x, i.y, i.z, result[i], result_new[i]);
    }
  }
  
  if(pass) {
    output.write("=> Test passed\n");

    output << "TIMING\n======\n";
    output << "Original : " << elapsed1.count() << std::endl;
    output << "New      : " << elapsed2.count() << std::endl;
  }else {
    output.write("=> Test failed\n");
  }

  /*
   * Clean up and exit
   */  
  BoutFinalise();
  
  if(!pass)
    return 1;
  return 0;
}

