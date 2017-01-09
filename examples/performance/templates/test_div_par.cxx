/*!
 * Testing code in src/difops/
 *
 *
 *  
 * TIMING
 * ======
 *
 * With checks disabled
 * Optimisation flag -O3
 *
 * Hand written      : 0.0102679
 * Templates         : 0.00942375
 * Original          : 0.0133223
 */

#include <chrono>
#include <iostream>

#include <bout.hxx>

// Hand-written version
const Field3D Div_par_C2(const Field3D &f);

// Template version
const Field3D Div_par_C2_TEST(const Field3D &f);

// Original
const Field3D Div_par_orig(const Field3D &f, CELL_LOC outloc=CELL_DEFAULT, DIFF_METHOD method=DIFF_DEFAULT) {
  Coordinates *coord = mesh->coordinates();
  // Need to modify yup and ydown fields
  Field3D f_B = f/coord->Bxy;
  if(&f.yup() == &f) {
    // Identity, yup and ydown point to same field
    f_B.mergeYupYdown();
  }else {
    // Distinct fields
    f_B.splitYupYdown();
    f_B.yup() = f.yup() / coord->Bxy;
    f_B.ydown() = f.ydown() / coord->Bxy;
  }
  Field3D result = coord->Bxy*Grad_par(f_B, outloc, method);
  
  return result;
}

int main(int argc, char **argv) {
  BoutInitialise(argc, argv);

  typedef std::chrono::time_point<std::chrono::steady_clock> SteadyClock;
  typedef std::chrono::duration<double> Duration;
  using namespace std::chrono;

  // Hand written version
  Field3D n;
  GRID_LOAD(n);
  mesh->communicate(n);
  SteadyClock start1 = steady_clock::now();
  Field3D result;
  for(int i=0;i<100;i++) {
     result = Div_par_C2(n);
  }
  Duration elapsed1 = steady_clock::now() - start1;
  
  // Template version
  
  SteadyClock start2 = steady_clock::now();
  Field3D result2;
  for(int i=0;i<100;i++) {
    result2= Div_par_C2_TEST(n);
  }
  Duration elapsed2 = steady_clock::now() - start2;

  // Old version
  
  SteadyClock start3 = steady_clock::now();
  Field3D result3;
  for(int i=0;i<100;i++) {
    result3 = Div_par_orig(n);
  }
  Duration elapsed3 = steady_clock::now() - start3;

  // Check differences
  for(auto i: n.region(RGN_NOBNDRY)) {
    if(fabs(result[i] != result2[i]) > 1e-10 || fabs(result2[i] != result3[i]) > 1e-10) {
      output << i.x << ", " << i.y << ", " << i.z << " : " << result[i] << ", " << result2[i] << ", " << result3[i] << std::endl;
    }
  }
  
  output << "TIMING\n======\n";
  output << "Hand written      : " << elapsed1.count() << std::endl;
  output << "Templates         : " << elapsed2.count() << std::endl;
  output << "Original          : " << elapsed3.count() << std::endl;

  BoutFinalise();
  return 0;
}
