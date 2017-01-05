#pragma once

#ifndef __DIFOP_TYPE_H__
#define __DIFOP_TYPE_H__

/// Choice of numerical method to use
enum class Difop {
  U1     ///< First order upwinding
    ,U2     ///< Second order upwinding
    ,C2     ///< Second order central difference. For Y derivatives this uses yup/ydown
    ,C4     ///< Fourth order central difference. For Y derivatives this uses yup/ydown
    ,U1_FA   ///< First order upwinding, Field Aligned
    ,U2_FA   ///< First order upwinding, Field Aligned
    ,C2_FA   ///< Second order central difference, Field Aligned
    ,W3_FA   ///< Third order WENO reconstruction, Field Aligned
    ,C4_FA   ///< Fourth order central difference, Field Aligned
    ,FFT     ///< Fourier Transform
    };

#endif // __DIFOP_TYPE_H__
