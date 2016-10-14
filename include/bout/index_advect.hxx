#pragma once

namespace DIFOPS {

/*
 * Methods chosen in options
 */

const Field2D indexVDDX(const Field2D &v, const Field2D &f);
const Field3D indexVDDX(const Field3D &v, const Field2D &f);
const Field3D indexVDDX(const Field2D &v, const Field3D &f);
const Field3D indexVDDX(const Field3D &v, const Field3D &f);

/*
 * Specific methods
 */

const Field2D indexVDDX_U1(const Field2D &v, const Field2D &f);
const Field2D indexVDDX_U2(const Field2D &v, const Field2D &f);
const Field2D indexVDDX_C2(const Field2D &v, const Field2D &f);
const Field2D indexVDDX_C4(const Field2D &v, const Field2D &f);
const Field2D indexVDDX_W3(const Field2D &v, const Field2D &f);

const Field3D indexVDDX_U1(const Field3D &v, const Field2D &f);
const Field3D indexVDDX_U2(const Field3D &v, const Field2D &f);
const Field3D indexVDDX_C2(const Field3D &v, const Field2D &f);
const Field3D indexVDDX_C4(const Field3D &v, const Field2D &f);
const Field3D indexVDDX_W3(const Field3D &v, const Field2D &f);

const Field3D indexVDDX_U1(const Field2D &v, const Field3D &f);
const Field3D indexVDDX_U2(const Field2D &v, const Field3D &f);
const Field3D indexVDDX_C2(const Field2D &v, const Field3D &f);
const Field3D indexVDDX_C4(const Field2D &v, const Field3D &f);
const Field3D indexVDDX_W3(const Field2D &v, const Field3D &f);

const Field3D indexVDDX_U1(const Field3D &v, const Field3D &f);
const Field3D indexVDDX_U2(const Field3D &v, const Field3D &f);
const Field3D indexVDDX_C2(const Field3D &v, const Field3D &f);
const Field3D indexVDDX_C4(const Field3D &v, const Field3D &f);
const Field3D indexVDDX_W3(const Field3D &v, const Field3D &f);



}; // DIFOPS namespace
