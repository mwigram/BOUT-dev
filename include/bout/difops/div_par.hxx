
#pragma once

#ifndef __DIV_PAR_H__
#define __DIV_PAR_H__

#include "difop_type.hxx"

#include <field2d.hxx>
#include <field3d.hxx>

////////////////////////////////////////////////////////////////////

/// Parallel divergence operator
/// 
/// \f[
///  B \partial_{||}(f/B) = B \nabla\cdot (\mathbf{b}f/B )
/// \f]
/// 
/// @param[in] f  The component of a vector along the magnetic field 
/// @param[in] method  The numerical method to use
///
/// This selects the function to call, based on the
/// method specified at run-time
const Field2D Div_par(const Field2D &f, Difop method = Difop::C2);

/// Parallel divergence
/// 2nd order central differencing
/// Assumes input is cell centred
const Field2D Div_par_C2(const Field2D &f);

/// Parallel divergence
/// 4th order central differencing
/// Assumes input is cell centred
const Field2D Div_par_C4(const Field2D &f);

////////////////////////////////////////////////////////////////////

/// Parallel divergence operator
/// 
/// \f[
///  B \partial_{||}(f/B) = B \nabla\cdot (\mathbf{b}f/B )
/// \f]
/// 
/// @param[in] f  The component of a vector along the magnetic field 
/// @param[in] method  The numerical method to use
///
/// This selects the function to call, based on the
/// cell location, and the method specified at run-time
const Field3D Div_par(const Field3D &f, Difop method = Difop::C2);

/// Parallel divergence
/// 2nd order central differencing
/// Cell centred input and output
/// Uses yup and ydown fields
const Field3D Div_par_C2(const Field3D &f);

/// Parallel divergence
/// 2nd order central differencing
/// CELL_YLOW input, CELL_CENTRE output
/// Uses yup and ydown fields
const Field3D Div_par_LtoC_C2(const Field3D &f);

/// Parallel divergence
/// 2nd order central differencing
/// CELL_YLOW input, CELL_CENTRE output
/// Transforms to and from field aligned coordinates,
/// so doesn't use yup and ydown fields
const Field3D Div_par_LtoC_C2_FA(const Field3D &f);

/// Parallel divergence
/// 2nd order central differencing
/// Cell centred input and output
/// Transforms to and from field aligned coordinates,
/// so doesn't use yup and ydown fields
const Field3D Div_par_C2_FA(const Field3D &f);

/// Parallel divergence
/// 4th order central differencing
/// Cell centred input and output
/// Transforms to and from field aligned coordinates,
/// so doesn't use yup and ydown fields
const Field3D Div_par_C4_FA(const Field3D &f);


#endif // __DIV_PAR_H__
