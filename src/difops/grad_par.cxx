
#include <globals.hxx>

#include <bout/difops/grad_par.hxx>
#include <bout/difops/ddy.cxx>

#include <msg_stack.hxx>

const Field2D Grad_par(const Field2D &f, Difop method) {
  TRACE("Grad_par(Field2D)");
  return DDY(f, method) / sqrt(mesh->coordinates()->g_22);
}

const Field3D Grad_par(const Field3D &f, Difop method) {
  TRACE("Grad_par(Field3D)");
  return DDY(f, method) / sqrt(mesh->coordinates()->g_22);
}
