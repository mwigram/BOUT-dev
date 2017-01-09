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
const F evalOperator(const Op &op, REGION region) {
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

/*
template<class F, class VS, class Op>
struct OpSingleArg {
  
  
  BoutReal eval(const DataIterator &i) {
    
  }
};
*/

/// 
template<class F, class IS, template<typename,typename> class VS > // 
class OpDivParC2 {
public:
  /// Initialise with a Field3D or Field2D
  OpDivParC2(const F &f) : vs(f) {
    coord = mesh->coordinates();
  }
  
  /// Evaluate at point i
  inline BoutReal eval(const DataIterator &i) const {
    // Calculate the index stencil (IS)
    IS is(i);
    // Get the stencil of values from field f
    Stencil1D s = vs.eval(is); 
    // Calculate the divergence operator
    return coord->Bxy[i] * ( s.p/coord->Bxy[is.p] - s.m/coord->Bxy[is.m] ) / (2.*coord->dy[i]*sqrt(coord->g_22[i]));
  }
private:
  Coordinates *coord;
  VS<F,IS> vs; // Value Stencil
};


const Field2D Div_par_C2_TEST(const Field2D &f) {
  OpDivParC2<Field2D, IndexStencilY, ValueStencil> op(f); // Operator
  return evalOperator<Field2D, OpDivParC2<Field2D, IndexStencilY, ValueStencil> >(op, RGN_NOBNDRY);
}

const Field3D Div_par_C2_TEST(const Field3D &f) {
  OpDivParC2<Field3D, IndexStencilY, ValueStencilYud> op(f); // Operator, using yup and ydown
  return evalOperator<Field3D, OpDivParC2<Field3D, IndexStencilY, ValueStencilYud> >(op, RGN_NOBNDRY);
}

const Field3D Div_par_C2_FA_TEST(const Field3D &f) {
  Field3D fa = mesh->toFieldAligned(f);
  OpDivParC2<Field3D, IndexStencilY, ValueStencil> op(fa);
  Field3D result = evalOperator<Field3D, OpDivParC2<Field3D, IndexStencilY, ValueStencil> >(op, RGN_NOBNDRY);
  return mesh->fromFieldAligned(result);
}
