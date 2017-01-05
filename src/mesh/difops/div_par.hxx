
#pragma once

#ifndef __DIV_PAR_H__
#define __DIV_PAR_H__

/// Choice of numerical method to use
enum class Difop {
  DEFAULT   ///< The default method
    ,C2     ///< Second order central difference. For Y derivatives this uses yup/ydown
    ,C4     ///< Fourth order central difference. For Y derivatives this uses yup/ydown
    ,C2_FA   ///< Second order central difference, Field Aligned
    ,C4_FA   ///< Fourth order central difference, Field Aligned
    };

////////////////////////////////////////////////////////////////////

/// Parallel divergence
/// This selects the function to call, based on the
/// method specified at run-time
const Field2D Div_par(const Field2D &f, Difop method = Difop::DEFAULT);

/// Parallel divergence
/// 2nd order central differencing
/// Assumes input is cell centred
const Field2D Div_par_C2(const Field2D &f);

/// Parallel divergence
/// 4th order central differencing
/// Assumes input is cell centred
const Field2D Div_par_C4(const Field2D &f);

////////////////////////////////////////////////////////////////////

const Field3D Div_par(const Field3D &f, Difop method = Difop::DEFAULT);

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
