#include "matlabfunctions.h"

#include <math.h>
#include <stdint.h>

namespace {
//-----------------------------------------------------------------------------
// FilterForDecimate() calculates the coefficients of low-pass filter and
// carries out the filtering. This function is only used for decimate().
//-----------------------------------------------------------------------------
static void FilterForDecimate(const float* x, int x_length, int r, float* y)
{
    float a[3], b[2]; // filter Coefficients
    switch (r) {
    case 11: // fs : 44100 (default)
        a[0] = 2.450743295230728f;
        a[1] = -2.06794904601978f;
        a[2] = 0.59574774438332101f;
        b[0] = 0.0026822508007163792f;
        b[1] = 0.0080467524021491377f;
        break;
    case 12: // fs : 48000
        a[0] = 2.4981398605924205f;
        a[1] = -2.1368928194784025f;
        a[2] = 0.62187513816221485f;
        b[0] = 0.0021097275904709001f;
        b[1] = 0.0063291827714127002f;
        break;
    case 10:
        a[0] = 2.3936475118069387f;
        a[1] = -1.9873904075111861f;
        a[2] = 0.5658879979027055f;
        b[0] = 0.0034818622251927556f;
        b[1] = 0.010445586675578267f;
        break;
    case 9:
        a[0] = 2.3236003491759578f;
        a[1] = -1.8921545617463598f;
        a[2] = 0.53148928133729068f;
        b[0] = 0.0046331164041389372f;
        b[1] = 0.013899349212416812f;
        break;
    case 8: // fs : 32000
        a[0] = 2.2357462340187593f;
        a[1] = -1.7780899984041358f;
        a[2] = 0.49152555365968692f;
        b[0] = 0.0063522763407111993f;
        b[1] = 0.019056829022133598f;
        break;
    case 7:
        a[0] = 2.1225239019534703f;
        a[1] = -1.6395144861046302f;
        a[2] = 0.44469707800587366f;
        b[0] = 0.0090366882681608418f;
        b[1] = 0.027110064804482525f;
        break;
    case 6: // fs : 24000 and 22050
        a[0] = 1.9715352749512141f;
        a[1] = -1.4686795689225347f;
        a[2] = 0.3893908434965701f;
        b[0] = 0.013469181309343825f;
        b[1] = 0.040407543928031475f;
        break;
    case 5:
        a[0] = 1.7610939654280557f;
        a[1] = -1.2554914843859768f;
        a[2] = 0.3237186507788215f;
        b[0] = 0.021334858522387423f;
        b[1] = 0.06400457556716227f;
        break;
    case 4: // fs : 16000
        a[0] = 1.4499664446880227f;
        a[1] = -0.98943497080950582f;
        a[2] = 0.24578252340690215f;
        b[0] = 0.036710750339322612f;
        b[1] = 0.11013225101796784f;
        break;
    case 3:
        a[0] = 0.95039378983237421f;
        a[1] = -0.67429146741526791f;
        a[2] = 0.15412211621346475f;
        b[0] = 0.071221945171178636f;
        b[1] = 0.21366583551353591f;
        break;
    case 2: // fs : 8000
        a[0] = 0.041156734567757189f;
        a[1] = -0.42599112459189636f;
        a[2] = 0.041037215479961225f;
        b[0] = 0.16797464681802227f;
        b[1] = 0.50392394045406674f;
        break;
    default:
        a[0] = 0.0f;
        a[1] = 0.0f;
        a[2] = 0.0f;
        b[0] = 0.0f;
        b[1] = 0.0f;
    }

    // Filtering on time domain.
    float w[3] = { 0.0, 0.0, 0.0 };
    float wt;
    for (int i = 0; i < x_length; ++i) {
        wt = x[i] + a[0] * w[0] + a[1] * w[1] + a[2] * w[2];
        y[i] = b[0] * wt + b[1] * w[0] + b[1] * w[1] + b[0] * w[2];
        w[2] = w[1];
        w[1] = w[0];
        w[0] = wt;
    }
}
} // namespace

void histc(const float* x, int x_length, const float* edges,
    int edges_length, int* index)
{
    int count = 1;

    int i = 0;
    for (; i < edges_length; ++i) {
        index[i] = 1;
        if (edges[i] >= x[0])
            break;
    }
    for (; i < edges_length; ++i) {
        if (edges[i] < x[count]) {
            index[i] = count;
        } else {
            index[i--] = count++;
        }
        if (count == x_length)
            break;
    }
    count--;
    for (i++; i < edges_length; ++i)
        index[i] = count;
}

void interp1(const float* x, const float* y, int x_length, const float* xi,
    int xi_length, float* yi)
{
    float* h = new float[x_length - 1];
    float* p = new float[xi_length];
    float* s = new float[xi_length];
    int* k = new int[xi_length];

    for (int i = 0; i < x_length - 1; ++i)
        h[i] = x[i + 1] - x[i];
    for (int i = 0; i < xi_length; ++i) {
        p[i] = (float)i;
        k[i] = 0;
    }

    histc(x, x_length, xi, xi_length, k);

    for (int i = 0; i < xi_length; ++i)
        s[i] = (xi[i] - x[k[i] - 1]) / h[k[i] - 1];

    for (int i = 0; i < xi_length; ++i)
        yi[i] = y[k[i] - 1] + s[i] * (y[k[i]] - y[k[i] - 1]);

    delete[] k;
    delete[] s;
    delete[] p;
    delete[] h;
}

void decimate(const float* x, int x_length, int r, float* y)
{
    const int kNFact = 9;
    float* tmp1 = new float[x_length + kNFact * 2];
    float* tmp2 = new float[x_length + kNFact * 2];

    for (int i = 0; i < kNFact; ++i)
        tmp1[i] = 2 * x[0] - x[kNFact - i];
    for (int i = kNFact; i < kNFact + x_length; ++i)
        tmp1[i] = x[i - kNFact];
    for (int i = kNFact + x_length; i < 2 * kNFact + x_length; ++i)
        tmp1[i] = 2 * x[x_length - 1] - x[x_length - 2 - (i - (kNFact + x_length))];

    FilterForDecimate(tmp1, 2 * kNFact + x_length, r, tmp2);
    for (int i = 0; i < 2 * kNFact + x_length; ++i)
        tmp1[i] = tmp2[2 * kNFact + x_length - i - 1];
    FilterForDecimate(tmp1, 2 * kNFact + x_length, r, tmp2);
    for (int i = 0; i < 2 * kNFact + x_length; ++i)
        tmp1[i] = tmp2[2 * kNFact + x_length - i - 1];

    int nout = (x_length - 1) / r + 1;
    int nbeg = r - r * nout + x_length;

    int count = 0;
    for (int i = nbeg; i < x_length + kNFact; i += r)
        y[count++] = tmp1[i + kNFact - 1];

    delete[] tmp1;
    delete[] tmp2;
}

int matlab_round(float x)
{
    return x > 0 ? static_cast<int>(x + 0.5) : static_cast<int>(x - 0.5);
}