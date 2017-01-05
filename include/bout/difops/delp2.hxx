
#pragma once

#ifndef __DELP2_H__
#define __DELP2_H__

#include "difop_type.hxx"

#include <field3d.hxx>

/// Perpendicular Laplacian operator
///
/// This version only includes terms in X and Z, dropping
/// derivatives in Y. This is the inverse operation to 
/// the Laplacian inversion class. 
///
/// For the full perpendicular Laplacian, use Laplace_perp
const Field3D Delp2(const Field3D &f, Difops method = Difops::FFT);


const Field3D Delp2_FFT(const Field3D &f);

#endif // __DELP2_H__

