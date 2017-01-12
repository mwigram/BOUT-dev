#include <globals.hxx>
#include <field3d.hxx>

/// Evaluate a function over a region
///
/// @tparam F   The Field class to return. Field2D or Field3D usually
/// @tparam Op  The operator class. Must have an eval(DataIterator) method
///
/// @param op      Instance of Op
/// @param region  The region to iterate over
template<class F, class Op> 
inline const F evalOperator(const Op &op, REGION region) {
  F result;
  result.allocate();
  
  for(auto i : result.region(region)) {
    // Calculate operation at point i
    result[i] = op.eval(i);
  }
  
  return result;
}

/////////////////////////////////////////////////////////////////////////
// Index stencils

/// Indices along X
struct IndexStencilX {
  IndexStencilX(const DataIterator &i) : c({i.x,i.y,i.z}), p(i.xp()), m(i.xm()) {}
  
  Indices c, p, m; // Centre, plus and minus
};

/// Indices along Y
struct IndexStencilY {
  IndexStencilY(const DataIterator &i) : c({i.x,i.y,i.z}), p(i.yp()), m(i.ym()) {}
  
  Indices c, p, m; // Centre, plus and minus
};

/// Indices along Z
struct IndexStencilZ {
  IndexStencilZ(const DataIterator &i) : c({i.x,i.y,i.z}), p(i.zp()), m(i.zm()) {}
  
  Indices c, p, m; // Centre, plus and minus
};

/////////////////////////////////////////////////////////////////////////
// Stencils from fields and indices

struct Stencil1D {
  BoutReal c,p,m;
};

template<class F, class IS>
class ValueStencil {
public:
  ValueStencil(const F &f) : val(f) {}
  
  inline const Stencil1D eval(const IS &is) const {
    return { val[is.c], val[is.p], val[is.m] };
  }
private:
  const F &val; // Field
};

template<class F, class IS>
class ValueStencilYud {
public:
  ValueStencilYud(const F &f) : val(f), yup(f.yup()), ydown(f.ydown()) {}
  
  inline const Stencil1D eval(const IS &is) const {
    return { val[is.c], yup[is.p], ydown[is.m] };
  }
private:
  const F &val;
  const F &yup;
  const F &ydown;
};


/////////////////////////////////////////////////////////////////////////
// Operators which can be passed into evalOperator
//
// Need to have a method:
//
//   BoutReal eval(const DataIterator &i) const
//


struct derivC2 {
  inline BoutReal eval(Stencil1D s) const {
    return (s.p - s.m)/2.;
  }
};


template<class F, typename Op, class IS, template<typename,typename> class VS = ValueStencil>
class OperatorYderiv {
public:
  /// Initialise with a Field3D or Field2D
  OperatorYderiv(const F &f) : vs(f) {
    coord = mesh->coordinates();
  }
  
  /// Evaluate at point i
  inline BoutReal eval(const DataIterator &i) const {
    IS is(i); // Get the index stencil
    return op.eval(vs.eval(is)) / coord->dy[i];
  }
private:
  Coordinates *coord;
  VS<F,IS> vs; // Value Stencil
  Op op; // Operator
};

const Field3D DDY_C2_TEST(const Field3D &f) {
  OperatorYderiv<Field3D, // Input is a Field3D
                 derivC2, // First derivative, central second order
                 IndexStencilY // Stencil in Y direction
                 > op(f);
  
  // Evaluate operator over a Field3D region
  return evalOperator<Field3D>(op, RGN_NOBNDRY);
}


/////////////////////////////////////////////////////////////////////////


/// Div_par operator, second order central difference
///
/// @tparam F    Field type (Field2D, Field3D). Unfortunately this can't be inferred from constructor
/// @tparam VS   Value Stencil. Calculates values from a set of indices.
template<class F, template<typename,typename> class VS> 
class OpDivParC2 {
public:
  /// Initialise with a Field3D or Field2D
  OpDivParC2(const F &f) : vs(f) {
    coord = mesh->coordinates();
  }
  
  /// Evaluate at point i
  inline BoutReal eval(const DataIterator &i) const {
    // Calculate the index stencil (IS) in the Y direction
    IndexStencilY is(i);
    // Get the stencil of values from field f
    Stencil1D s = vs.eval(is); 
    // Calculate the divergence operator
    return coord->Bxy[i] * ( s.p/coord->Bxy[is.p] - s.m/coord->Bxy[is.m] ) / (2.*coord->dy[i]*sqrt(coord->g_22[i]));
  }
private:
  Coordinates *coord;
  VS<F,IndexStencilY> vs; // Value Stencil
};


const Field2D Div_par_C2_TEST(const Field2D &f) {
  // Create the operator
  OpDivParC2<Field2D,       // Operating on an input Field2D
             ValueStencil   // Simple calculation of stencil values (no yup/ydown) 
             > op(f);
  // Evaluate operator over a Field2D region
  return evalOperator<Field2D>(op, RGN_NOBNDRY);
}

const Field3D Div_par_C2_TEST(const Field3D &f) {
  OpDivParC2<Field3D,       // Operates on Field3D input
             ValueStencilYud  // Use Yup and Ydown to get values from stencils
             > op(f);
  return evalOperator<Field3D>(op, RGN_NOBNDRY);
}

const Field3D Div_par_C2_FA_TEST(const Field3D &f) {
  Field3D fa = mesh->toFieldAligned(f); // Transform to field aligned coordinates
  OpDivParC2<Field3D,      // Operate on Field3D input
             ValueStencil  // No yup/ydown fields
             > op(fa);
  Field3D result = evalOperator<Field3D>(op, RGN_NOBNDRY);
  return mesh->fromFieldAligned(result); // Transform back
}
