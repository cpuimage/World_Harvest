#include "common.h"

#include <math.h>

#define STB_FFT_IMPLEMENTAION

#include "stb_fft.h"

#include "matlabfunctions.h"

//-----------------------------------------------------------------------------
int GetSuitableFFTSize(int sample)
{
    return static_cast<int>(powf(2.0f,
        static_cast<int>(logf(static_cast<float>(sample)) / kLog2) + 1.0f));
}

void NuttallWindow(int y_length, float* y)
{
    float tmp;
    for (int i = 0; i < y_length; ++i) {
        tmp = i / (y_length - 1.0f);
        y[i] = 0.355768f - 0.487396f * cosf(2.0f * kPi * tmp) + 0.144232f * cosf(4.0f * kPi * tmp) - 0.012604f * cosf(6.0f * kPi * tmp);
    }
}

//-----------------------------------------------------------------------------
// FFT, IFFT and minimum phase analysis
void InitializeForwardRealFFT(int fft_size, ForwardRealFFT* forward_real_fft)
{
    forward_real_fft->fft_size = fft_size;
    forward_real_fft->waveform = new float[fft_size];
    forward_real_fft->spectrum = new cmplx[fft_size];
    forward_real_fft->forward_fft = stb_fft_real_plan_dft_1d(fft_size);
}

void DestroyForwardRealFFT(ForwardRealFFT* forward_real_fft)
{
    if (forward_real_fft->forward_fft != NULL)
        free(forward_real_fft->forward_fft);
    delete[] forward_real_fft->spectrum;
    delete[] forward_real_fft->waveform;
}