
#pragma once

#ifndef __DDY_H__
#define __DDY_H__

#include "difop_type.hxx"

#include <field2d.hxx>
#include <field3d.hxx>

////////////////////////////////////////////////////////////////////
// Field2D operators

/// Calculate first partial derivative in Y
///
///  $\partial / \partial y$
///
/// @param[in] f       The field to be differentiated
/// @param[in] method  Differencing method to use. 
const Field2D DDY(const Field2D &f, Difop method = Difop::C2);

/// Partial derivative in Y
///
/// Second order central differencing
///
/// @param[in] f  The field to be differentiated
const Field2D DDY_C2(const Field2D &f);

/// Partial derivative in Y
///
/// Fourth order central differencing
///
/// @param[in] f  The field to be differentiated
const Field2D DDY_C4(const Field2D &f);

////////////////////////////////////////////////////////////////////
// Field3D operators

/// Calculate first partial derivative in Y
///
///  $\partial / \partial y$
///
/// @param[in] f       The field to be differentiated
/// @param[in] method  Differencing method to use. 
const Field3D DDY(const Field3D &f, Difop method = Difop::C2);

/// Partial derivative in Y
///
/// Second order central difference
/// Input and output cell centred
/// Uses yup and ydown fields
const Field3D DDY_C2(const Field3D &f);

/// Partial derivative in Y
///
/// Second order central difference
/// Input and output cell centred
/// Transforms to/from field aligned coordinates
const Field3D DDY_C2_FA(const Field3D &f);

#endif //  __DDY_H__
