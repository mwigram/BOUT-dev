#pragma once

#ifndef __GRAD_PAR_H__
#define __GRAD_PAR_H__

#include "difop_type.hxx"

#include <field2d.hxx>
#include <field3d.hxx>

/// Parallel derivative (central differencing) in Y
/// along unperturbed field 
///
/// @param[in] f  The field to be differentiated
/// @param[in] outloc  The cell location where the output is needed (if staggered grids is enabled)
/// @param[in] method  The method to use
///
const Field2D Grad_par(const Field2D &f, Difop method = Difop::C2);

/// Parallel derivative (central differencing) in Y
/// along unperturbed field 
///
/// @param[in] f  The field to be differentiated
/// @param[in] outloc  The cell location where the output is needed (if staggered grids is enabled)
/// @param[in] method  The method to use
///
const Field3D Grad_par(const Field3D &f, Difop method = Difop::C2);

#endif  // __GRAD_PAR_H__
