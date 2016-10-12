/*
 * Non-staggered first derivatives in index space
 */
#pragma once

#include <field3d.hxx>
#include <field2d.hxx>

/*
 * These functions use method specified in options
 * so have a small overhead for the dispatch
 */

const Field3D indexDDX(const Field3D &f, CELL_LOC outloc = CELL_DEFAULT);
const Field2D indexDDX(const Field2D &f, CELL_LOC outloc = CELL_DEFAULT);

const Field3D indexDDY(const Field3D &f, CELL_LOC outloc = CELL_DEFAULT);
const Field2D indexDDY(const Field2D &f, CELL_LOC outloc = CELL_DEFAULT);

const Field3D indexDDZ(const Field3D &f, CELL_LOC outloc = CELL_DEFAULT);
const Field2D indexDDZ(const Field2D &f, CELL_LOC outloc = CELL_DEFAULT);

/*
 * These functions use a specific method
 */

const Field2D indexDDX_C2(const Field2D &f);
const Field3D indexDDX_C2(const Field3D &f);
const Field2D indexDDX_C4(const Field2D &f);
const Field3D indexDDX_C4(const Field3D &f);

const Field2D indexDDY_C2(const Field2D &f);
const Field3D indexDDY_C2(const Field3D &f);

const Field2D indexDDZ_C2(const Field2D &f) { return 0.0; }
const Field3D indexDDZ_C2(const Field3D &f);
const Field2D indexDDZ_C4(const Field2D &f) { return 0.0; }
const Field3D indexDDZ_C4(const Field3D &f);

