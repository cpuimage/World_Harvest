
#ifndef _COMMON_H_
#define _COMMON_H_

#include "stb_fft.h"

#ifndef MACRODEFINITIONS_H_
#define MACRODEFINITIONS_H_

#undef BEGIN_C_DECLS
#undef END_C_DECLS
#ifdef __cplusplus
#define BEGIN_C_DECLS extern "C" {
#define END_C_DECLS }
#else // !__cplusplus
#define BEGIN_C_DECLS
#define END_C_DECLS
#endif // __cplusplus

#endif

BEGIN_C_DECLS

#ifndef MAX
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#endif

//-----------------------------------------------------------------------------
// Structs on FFT
//-----------------------------------------------------------------------------
// Forward FFT in the real sequence
typedef struct {
    int fft_size;
    float* waveform;
    cmplx* spectrum;
    stb_fft_real_plan* forward_fft;
} ForwardRealFFT;

//-----------------------------------------------------------------------------
// GetSuitableFFTSize() calculates the suitable FFT size.
// The size is defined as the minimum length whose length is longer than
// the input sample.
//
// Input:
//   sample : Length of the input signal
//
// Output:
//   Suitable FFT size
//-----------------------------------------------------------------------------
int GetSuitableFFTSize(int sample);

const float kFloorF0 = 71.0f;
const float kCeilF0 = 800.0f;
const float kDefaultF0 = 500.0f;
const float kPi = 3.1415926535897932384f;
const float kMySafeGuardMinimum = 0.000000000001f;
const float kLog2 = 0.69314718055994529f;

//-----------------------------------------------------------------------------
// NuttallWindow() calculates the coefficients of Nuttall window whose length
// is y_length and is used in  Harvest().
//-----------------------------------------------------------------------------
void NuttallWindow(int y_length, float* y);

//-----------------------------------------------------------------------------
// These functions are used to speed up the processing.
// Forward FFT
void InitializeForwardRealFFT(int fft_size, ForwardRealFFT* forward_real_fft);

void DestroyForwardRealFFT(ForwardRealFFT* forward_real_fft);

END_C_DECLS

#endif // _COMMON_H_
