
#ifndef _AUDIOIO_H_
#define _AUDIOIO_H_

#ifdef __cplusplus
extern "C" {
#endif

//-----------------------------------------------------------------------------
// wavwrite() write a .wav file.
// Input:
//   x          : Input signal
//   x_ength : Signal length of x [sample]
//   fs         : Sampling frequency [Hz]
//   nbit       : Quantization bit [bit]
//   filename   : Name of the output signal.
// Caution:
//   The variable nbit is not used in this function.
//   This function only supports the 16 bit.
//-----------------------------------------------------------------------------
void wavwrite(const float* x, int x_length, int fs, int nbit,
    const char* filename);

//-----------------------------------------------------------------------------
// wavread() read a .wav file.
// The memory of output x must be allocated in advance.
// Input:
//   filename     : Filename of the input file.
// Output:
//   fs           : Sampling frequency [Hz]
//-----------------------------------------------------------------------------
float* wavread(const char* filename, int* fs, int* sampleCount);
#ifdef __cplusplus
}
#endif

#endif // _AUDIOIO_H_
