
#ifndef _MATLABFUNCTIONS_H_
#define _MATLABFUNCTIONS_H_

#include "common.h"

BEGIN_C_DECLS

//-----------------------------------------------------------------------------
// histc() counts the number of values in vector x that fall between the
// elements in the edges vector (which must contain monotonically
// nondecreasing values). n is a length(edges) vector containing these counts.
// No elements of x can be complex.
// http://www.mathworks.co.jp/help/techdoc/ref/histc.html
//
// Input:
//   x              : Input vector
//   x_length       : Length of x
//   edges          : Input matrix (1-dimension)
//   edges_length   : Length of edges
//
// Output:
//   index          : Result counted in vector x
// Caution:
//   Lengths of index and edges must be the same.
//-----------------------------------------------------------------------------
void histc(const float* x, int x_length, const float* edges,
    int edges_length, int* index);

//-----------------------------------------------------------------------------
// interp1() interpolates to find yi, the values of the underlying function Y
// at the points in the vector or array xi. x must be a vector.
// http://www.mathworks.co.jp/help/techdoc/ref/interp1.html
//
// Input:
//   x          : Input vector (Time axis)
//   y          : Values at x[n]
//   x_length   : Length of x (Length of y must be the same)
//   xi         : Required vector
//   xi_length  : Length of xi (Length of yi must be the same)
//
// Output:
//   yi         : Interpolated vector
//-----------------------------------------------------------------------------
void interp1(const float* x, const float* y, int x_length, const float* xi,
    int xi_length, float* yi);

//-----------------------------------------------------------------------------
// decimate() carries out down sampling by both IIR and FIR filters.
// Filter coeffiencts are based on FilterForDecimate().
//
// Input:
//   x          : Input signal
//   x_length   : Length of x
//   r          : Coefficient used for down sampling
//                (fs after down sampling is fs/r)
// Output:
//   y          : Output signal
//-----------------------------------------------------------------------------
void decimate(const float* x, int x_length, int r, float* y);

//-----------------------------------------------------------------------------
// matlab_round() calculates rounding.
//
// Input:
//   x    : Input value
//
// Output:
//   y    : Rounded value
//-----------------------------------------------------------------------------
int matlab_round(float x);

END_C_DECLS

#endif // _MATLABFUNCTIONS_H_
