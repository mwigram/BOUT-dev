/*
 * Loops over the domain, providing stencils to operators
 *
 * The naming convention used is   loop<Direction><Args>
 *
 * e.g. loopXS  -- Loop with stencils in X direction, taking a single Stencil input
 *      
 *      loopYRS -- Stencils in the Y direction. Two inputs, first Real, second Stencil
 * 
 */

#pragma once

#include <globals.hxx>
#include <bout/dataiterator.hxx>

/*!
 * Stencil class, used to shorten the input arguments to derivative operators
 */
struct Stencil {
  BoutReal mm;
  BoutReal m;
  BoutReal c;
  BoutReal p;
  BoutReal pp;
};

/*!
 * Template class FieldOpResult
 *
 * Calculates the return type of an operation involving Field2D and Field3D types
 *
 * FieldOpResult<Field2D, Field2D>::type   -> Field2D
 * FieldOpResult<Field2D, Field3D>::type   -> Field3D
 * 
 */
template<typename F1, typename F2>
struct FieldOpResult {
  typedef Field3D type;  // Default result is Field3D
};

// If both inputs are Field2D, result is Field2D
template<>
struct FieldOpResult<Field2D,Field2D> {
  typedef Field2D type;
};

/*!
 * Template class to loop over the domain RGN_NOBNDRY
 * and call a derivative operator in X
 * 
 */
template<typename Op>
struct loopXS {
  Op op;
  
  template<typename F>
  const F operator()(const F &f) const {
    F result;
    result.allocate();
    
    for(auto i : result.region(RGN_NOBNDRY)) {
      
      result[i] = op({
          f[i.offset(-2,0,0)],
            f[i.xm()],
            f[i],
            f[i.xp()],
            f[i.offset(2,0,0)]
            });
    }
    return result;
  }
};

/*!
 * Template class to loop over the domain RGN_NOBNDRY
 * and call a derivative operator in X which takes two inputs (BoutReal, Stencil)
 * 
 */
template<typename Op>
struct loopXRS {
  Op op;
  
  template<typename V, typename F>
  const typename FieldOpResult<V, F>::type operator()(const V &v, const F &f) const {
    typename FieldOpResult<V, F>::type result;
    result.allocate();
    
    for(auto i : result.region(RGN_NOBNDRY)) {
      
      result[i] = op(v[i], {
          f[i.offset(-2,0,0)],
            f[i.xm()],
            f[i],
            f[i.xp()],
            f[i.offset(2,0,0)]
            });
    }
    return result;
  }
};

/*!
 * Template class to loop over the domain RGN_NOBNDRY
 * and call a derivative operator in Y
 * 
 */
template<typename Op>
struct loopYS {
  Op op;
  
  template<typename F>
  const F operator()(const F &f) const {
    F result;
    result.allocate();
    
    F yup = f.yup();
    F ydown = f.ydown();
    
    for(auto i : result.region(RGN_NOBNDRY)) {
      
      result[i] = op({
          nan(""),       // Note: NaNs used since +/- 2 not yet available
            ydown[i.ym()],
            f[i],
            yup[i.yp()],
            nan("")
            });
    }
    return result;
  }
};

/*!
 * Template class to loop over the domain RGN_NOBNDRY
 * and call a derivative operator in Y which takes two inputs (BoutReal, Stencil)
 * 
 */
template<typename Op>
struct loopYRS {
  Op op;
  
  template<typename V, typename F>
  const typename FieldOpResult<V, F>::type operator()(const V &v, const F &f) const {
    typename FieldOpResult<V, F>::type result;
    result.allocate();
    
    F yup = f.yup();
    F ydown = f.ydown();
    
    for(auto i : result.region(RGN_NOBNDRY)) {
      
      result[i] = op(v[i], {
          nan(""),
            ydown[i.ym()],
            f[i],
            yup[i.yp()],
            nan("")
            });
    }
    return result;
  }
};

/*!
 * Template class to loop over the domain RGN_NOBNDRY
 * and call a derivative operator in Z
 * 
 */
template<typename Op>
struct loopZS {
  Op op;
  
  template<typename F>
  const F operator()(const F &f) const {
    F result;
    result.allocate();
    
    for(auto i : result.region(RGN_NOBNDRY)) {
      
      result[i] = op({
          f[i.offset(0,0,-2)],
            f[i.zm()],
            f[i],
            f[i.zp()],
            f[i.offset(0,0,2)]
            });
    }
    return result;
  }
};

/*!
 * Template class to loop over the domain RGN_NOBNDRY
 * and call a derivative operator in Z which takes two inputs (BoutReal, Stencil)
 * 
 */
template<typename Op>
struct loopZRS {
  Op op;
  
  template<typename V, typename F>
  const typename FieldOpResult<V, F>::type operator()(const V &v, const F &f) const {
    typename FieldOpResult<V, F>::type result;
    result.allocate();
    
    for(auto i : result.region(RGN_NOBNDRY)) {
      
      result[i] = op(v[i], {
          f[i.offset(0,0,-2)],
            f[i.zm()],
            f[i],
            f[i.zp()],
            f[i.offset(0,0,2)]
            });
    }
    return result;
  }
};
