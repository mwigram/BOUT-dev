#pragma once

#ifndef __VPAR_GRAD_PAR_H__
#define __VPAR_GRAD_PAR_H__

#include "difop_type.hxx"

#include <field2d.hxx>
#include <field3d.hxx>

/// Advection of field f by velocity v along magnetic field 
///
/// v * b dot Grad(f)
///
/// All fields cell centred
const Field2D Vpar_Grad_par(const Field2D &v, const Field2D &f, Difop method = Difop::U1);

/// Advection of field f by velocity v along magnetic field 
///
/// v * b dot Grad(f)
///
/// First order upwinding
/// All fields cell centred
const Field2D Vpar_Grad_par_U1(const Field2D &v, const Field2D &f);

/// Advection of field f by velocity v along magnetic field 
///
/// v * b dot Grad(f)
///
/// Second order upwinding
/// All fields cell centred
const Field2D Vpar_Grad_par_U2(const Field2D &v, const Field2D &f);

/// Advection of field f by velocity v along magnetic field 
///
/// v * b dot Grad(f)
///
/// Second order central differencing
/// All fields cell centred
const Field2D Vpar_Grad_par_C2(const Field2D &v, const Field2D &f);

/// Advection of field f by velocity v along magnetic field
/// 
/// v * b dot Grad(f)
const Field3D Vpar_Grad_par(const Field3D &v, const Field3D &f, Difop method = Difop::U1);

/// Advection of field f by velocity v along magnetic field
/// 
/// v * b dot Grad(f)
///
/// First order upwinding
/// All fields cell centred
/// Uses yup() and ydown() fields
const Field3D Vpar_Grad_par_U1(const Field3D &v, const Field3D &f);

/// Advection of field f by velocity v along magnetic field
/// 
/// v * b dot Grad(f)
///
/// Second order central differencing
/// All fields cell centred
/// Uses yup() and ydown() fields
const Field3D Vpar_Grad_par_C2(const Field3D &v, const Field3D &f);

/// Advection of field f by velocity v along magnetic field
/// 
/// v * b dot Grad(f)
///
/// First order upwinding
/// All fields cell centred
/// Transforms to/from Field Aligned coordinates
const Field3D Vpar_Grad_par_U1_FA(const Field3D &v, const Field3D &f);

/// Advection of field f by velocity v along magnetic field
/// 
/// v * b dot Grad(f)
///
/// Second order central differencing
/// All fields cell centred
/// Transforms to/from Field Aligned coordinates
const Field3D Vpar_Grad_par_C2_FA(const Field3D &v, const Field3D &f);

/// Advection of field f by velocity v along magnetic field
/// 
/// v * b dot Grad(f)
///
/// Third order WENO reconstruction
/// All fields cell centred
/// Transforms to/from Field Aligned coordinates
const Field3D Vpar_Grad_par_W3_FA(const Field3D &v, const Field3D &f);

#endif // __VPAR_GRAD_PAR_H__
