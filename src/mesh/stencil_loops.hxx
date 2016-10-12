/*
 * Loops over the domain, providing stencils to operators
 *
 * Operators are expected to have an operator() function which takes a single Stencil as input
 * 
 */

#pragma once

#include <globals.hxx>
#include <bout/dataiterator.hxx>

struct Stencil {
  BoutReal mm;
  BoutReal m;
  BoutReal c;
  BoutReal p;
  BoutReal pp;
};

/*!
 * Template class to loop over the domain RGN_NOBNDRY
 * and call a derivative operator in X
 * 
 */
template<typename Op>
struct loopX1 {
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
 * and call a derivative operator in Y
 * 
 */
template<typename Op>
struct loopY1 {
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
 * and call a derivative operator in Z
 * 
 */
template<typename Op>
struct loopZ1 {
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

