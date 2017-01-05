
#include <bout/difops/delp2.hxx>

const Field3D Delp2(const Field3D &f, Difops method) {
  TRACE("Delp2( Field3D )");
  
  switch(method) {
  case Difops::FFT: return Delp2_FFT(f);
  }
  
  throw BoutException("method not supported");
}


/// Perpendicular Laplacian operator in X-Z
/// 
/// 
///
const Field3D Delp2_FFT(const Field3D &f) {
  TRACE("Delp2_FFT( Field3D )");

  //return mesh->G1*DDX(f) + mesh->G3*DDZ(f) + mesh->g11*D2DX2(f) + mesh->g33*D2DZ2(f); //+ 2.0*mesh->g13*D2DXDZ(f)

  ASSERT2(mesh->xstart > 0); // Need at least one guard cell
  
  Field3D result;
  result.allocate();
  
  int ncz = mesh->LocalNz;
  
  static dcomplex **ft = (dcomplex**) NULL, **delft;
  if(ft == (dcomplex**) NULL) {
    // Allocate memory
    ft = matrix<dcomplex>(mesh->LocalNx, ncz/2 + 1);
    delft = matrix<dcomplex>(mesh->LocalNx, ncz/2 + 1);
  }
  
  // Loop over all y indices
  for(int jy=0;jy<mesh->LocalNy;jy++) {

    // Take forward FFT
    
    for(int jx=0;jx<mesh->LocalNx;jx++)
      rfft(&f(jx,jy,0), ncz, ft[jx]);

    // Loop over kz
    for(int jz=0;jz<=ncz/2;jz++) {
      dcomplex a, b, c;

      // No smoothing in the x direction
      for(int jx=mesh->xstart;jx<=mesh->xend;jx++) {
	// Perform x derivative
	
	laplace_tridag_coefs(jx, jy, jz, a, b, c);

	delft[jx][jz] = a*ft[jx-1][jz] + b*ft[jx][jz] + c*ft[jx+1][jz];
      }
    }
  
    // Reverse FFT
    for(int jx=mesh->xstart;jx<=mesh->xend;jx++) {

      irfft(delft[jx], ncz, &result(jx,jy,0));
    }

    // Boundaries
    for(int jz=0;jz<ncz;jz++) {
      result(0,jy,jz) = 0.0;
      result(mesh->LocalNx-1,jy,jz) = 0.0;
    }
  }
  
  // Set the output location
  result.setLocation(f.getLocation());

  return result;
}
