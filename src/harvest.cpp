//-----------------------------------------------------------------------------
// F0 estimation based on Harvest.
//-----------------------------------------------------------------------------
#include "harvest.h"

#include <math.h>

#include "common.h"

#include "matlabfunctions.h"

//-----------------------------------------------------------------------------
// struct for RawEventByHarvest()
// "negative" means "zero-crossing point going from positive to negative"
// "positive" means "zero-crossing point going from negative to positive"
//-----------------------------------------------------------------------------
typedef struct {
    float* negative_interval_locations;
    float* negative_intervals;
    int number_of_negatives;
    float* positive_interval_locations;
    float* positive_intervals;
    int number_of_positives;
    float* peak_interval_locations;
    float* peak_intervals;
    int number_of_peaks;
    float* dip_interval_locations;
    float* dip_intervals;
    int number_of_dips;
} ZeroCrossings;

#include <unordered_map>

std::unordered_map<int, ForwardRealFFT> g_ForwardRealFFT;

namespace {
//-----------------------------------------------------------------------------
// Since the waveform of beginning and ending after decimate include noise,
// the input waveform is extended. This is the processing for the
// compatibility with MATLAB version.
//-----------------------------------------------------------------------------
static void GetWaveformAndSpectrumSub(const float* x, int x_length,
    int y_length, float actual_fs, int decimation_ratio, float* y)
{
    if (decimation_ratio == 1) {
        for (int i = 0; i < x_length; ++i)
            y[i] = x[i];
        return;
    }

    int lag = static_cast<int>(ceil(140.0 / decimation_ratio) * decimation_ratio);
    int new_x_length = x_length + lag * 2;
    float* new_y = new float[new_x_length];
    memset(new_y, 0, sizeof(*new_y) * new_x_length);
    float* new_x = new float[new_x_length];
    for (int i = 0; i < lag; ++i)
        new_x[i] = x[0];
    for (int i = lag; i < lag + x_length; ++i)
        new_x[i] = x[i - lag];
    for (int i = lag + x_length; i < new_x_length; ++i)
        new_x[i] = x[x_length - 1];

    decimate(new_x, new_x_length, decimation_ratio, new_y);
    for (int i = 0; i < y_length; ++i)
        y[i] = new_y[lag / decimation_ratio + i];

    delete[] new_x;
    delete[] new_y;
}

//-----------------------------------------------------------------------------
// GetWaveformAndSpectrum() calculates the downsampled signal and its spectrum
//-----------------------------------------------------------------------------
static void GetWaveformAndSpectrum(const float* x, int x_length,
    int y_length, float actual_fs, stb_fft_real_plan* fft_plan, int fft_size,
    int decimation_ratio,
    float* y, cmplx* y_spectrum)
{
    // Initialization
    memset(y, 0, fft_size * sizeof(*y));
    // Processing for the compatibility with MATLAB version
    GetWaveformAndSpectrumSub(x, x_length, y_length, actual_fs,
        decimation_ratio, y);

    // Removal of the DC component (y = y - mean value of y)
    float mean_y = 0.0;
    for (int i = 0; i < y_length; ++i)
        mean_y += y[i];
    mean_y /= y_length;
    for (int i = 0; i < y_length; ++i)
        y[i] -= mean_y;
    memset(y + y_length, 0, (fft_size - y_length) * sizeof(*y));
    stb_fft_r2c_exec(fft_plan, y, y_spectrum);
}

//-----------------------------------------------------------------------------
// GetFilteredSignal() calculates the signal that is the convolution of the
// input signal and band-pass filter.
//-----------------------------------------------------------------------------
static void GetFilteredSignal(float boundary_f0, stb_fft_real_plan* fft_plan, int fft_size, float fs,
    const cmplx* y_spectrum, int y_length, float* filtered_signal)
{
    int filter_length_half = matlab_round(fs / boundary_f0 * 2.f);
    float* band_pass_filter = new float[fft_size];
    NuttallWindow(filter_length_half * 2 + 1, band_pass_filter);
    for (int i = -filter_length_half; i <= filter_length_half; ++i)
        band_pass_filter[i + filter_length_half] *= cosf(2 * kPi * boundary_f0 * i / fs);
    memset(band_pass_filter + (filter_length_half * 2 + 1), 0,
        fft_size - (filter_length_half * 2 + 1) * sizeof(*band_pass_filter));
    cmplx* band_pass_filter_spectrum = new cmplx[fft_size];
    stb_fft_r2c_exec(fft_plan, band_pass_filter, band_pass_filter_spectrum);

    // Convolution
    float tmp = y_spectrum[0].real * band_pass_filter_spectrum[0].real - y_spectrum[0].imag * band_pass_filter_spectrum[0].imag;
    band_pass_filter_spectrum[0].imag = y_spectrum[0].real * band_pass_filter_spectrum[0].imag + y_spectrum[0].imag * band_pass_filter_spectrum[0].real;
    band_pass_filter_spectrum[0].real = tmp;
    for (int i = 1; i <= fft_size / 2; ++i) {
        tmp = y_spectrum[i].real * band_pass_filter_spectrum[i].real - y_spectrum[i].imag * band_pass_filter_spectrum[i].imag;
        band_pass_filter_spectrum[i].imag = y_spectrum[i].real * band_pass_filter_spectrum[i].imag + y_spectrum[i].imag * band_pass_filter_spectrum[i].real;
        band_pass_filter_spectrum[i].real = tmp;
        band_pass_filter_spectrum[fft_size - i - 1].real = band_pass_filter_spectrum[i].real;
        band_pass_filter_spectrum[fft_size - i - 1].imag = band_pass_filter_spectrum[i].imag;
    }

    stb_fft_c2r_exec(fft_plan, band_pass_filter_spectrum, filtered_signal);

    // Compensation of the delay.
    int index_bias = filter_length_half + 1;
    for (int i = 0; i < y_length; ++i)
        filtered_signal[i] = filtered_signal[i + index_bias];

    delete[] band_pass_filter_spectrum;
    delete[] band_pass_filter;
}

//-----------------------------------------------------------------------------
// CheckEvent() returns 1, provided that the input value is over 1.
// This function is for RawEventByDio().
//-----------------------------------------------------------------------------
static inline int CheckEvent(int x)
{
    return x > 0 ? 1 : 0;
}

//-----------------------------------------------------------------------------
// ZeroCrossingEngine() calculates the zero crossing points from positive to
// negative.
//-----------------------------------------------------------------------------
static int ZeroCrossingEngine(const float* filtered_signal, int y_length,
    float fs, float* interval_locations, float* intervals)
{
    int* negative_going_points = new int[y_length];

    for (int i = 0; i < y_length - 1; ++i)
        negative_going_points[i] = 0.0 < filtered_signal[i] && filtered_signal[i + 1] <= 0.0 ? i + 1 : 0;
    negative_going_points[y_length - 1] = 0;

    int* edges = new int[y_length];
    int count = 0;
    for (int i = 0; i < y_length; ++i)
        if (negative_going_points[i] > 0)
            edges[count++] = negative_going_points[i];

    if (count < 2) {
        delete[] edges;
        delete[] negative_going_points;
        return 0;
    }

    float* fine_edges = new float[count];
    for (int i = 0; i < count; ++i)
        fine_edges[i] = edges[i] - filtered_signal[edges[i] - 1] / (filtered_signal[edges[i]] - filtered_signal[edges[i] - 1]);

    for (int i = 0; i < count - 1; ++i) {
        intervals[i] = fs / (fine_edges[i + 1] - fine_edges[i]);
        interval_locations[i] = (fine_edges[i] + fine_edges[i + 1]) / 2.0f / fs;
    }

    delete[] fine_edges;
    delete[] edges;
    delete[] negative_going_points;
    return count - 1;
}

//-----------------------------------------------------------------------------
// GetFourZeroCrossingIntervals() calculates four zero-crossing intervals.
// (1) Zero-crossing going from negative to positive.
// (2) Zero-crossing going from positive to negative.
// (3) Peak, and (4) dip. (3) and (4) are calculated from the zero-crossings of
// the differential of waveform.
//-----------------------------------------------------------------------------
static void GetFourZeroCrossingIntervals(float* filtered_signal, int y_length,
    float actual_fs, ZeroCrossings* zero_crossings)
{
    int maximum_number = y_length;
    zero_crossings->negative_interval_locations = new float[maximum_number];
    zero_crossings->positive_interval_locations = new float[maximum_number];
    zero_crossings->peak_interval_locations = new float[maximum_number];
    zero_crossings->dip_interval_locations = new float[maximum_number];
    zero_crossings->negative_intervals = new float[maximum_number];
    zero_crossings->positive_intervals = new float[maximum_number];
    zero_crossings->peak_intervals = new float[maximum_number];
    zero_crossings->dip_intervals = new float[maximum_number];

    zero_crossings->number_of_negatives = ZeroCrossingEngine(filtered_signal,
        y_length, actual_fs,
        zero_crossings->negative_interval_locations,
        zero_crossings->negative_intervals);

    for (int i = 0; i < y_length; ++i)
        filtered_signal[i] = -filtered_signal[i];
    zero_crossings->number_of_positives = ZeroCrossingEngine(filtered_signal,
        y_length, actual_fs,
        zero_crossings->positive_interval_locations,
        zero_crossings->positive_intervals);

    for (int i = 0; i < y_length - 1; ++i)
        filtered_signal[i] = filtered_signal[i] - filtered_signal[i + 1];
    zero_crossings->number_of_peaks = ZeroCrossingEngine(filtered_signal,
        y_length - 1, actual_fs,
        zero_crossings->peak_interval_locations,
        zero_crossings->peak_intervals);

    for (int i = 0; i < y_length - 1; ++i)
        filtered_signal[i] = -filtered_signal[i];
    zero_crossings->number_of_dips = ZeroCrossingEngine(filtered_signal,
        y_length - 1, actual_fs,
        zero_crossings->dip_interval_locations,
        zero_crossings->dip_intervals);
}

static void GetF0CandidateContourSub(const float* const* interpolated_f0_set,
    int f0_length, float f0_floor, float f0_ceil, float boundary_f0,
    float* f0_candidate)
{
    float upper = boundary_f0 * 1.1f;
    float lower = boundary_f0 * 0.9f;
    for (int i = 0; i < f0_length; ++i) {
        f0_candidate[i] = (interpolated_f0_set[0][i] + interpolated_f0_set[1][i] + interpolated_f0_set[2][i] + interpolated_f0_set[3][i]) / 4.0f;

        if (f0_candidate[i] > upper || f0_candidate[i] < lower || f0_candidate[i] > f0_ceil || f0_candidate[i] < f0_floor)
            f0_candidate[i] = 0.0;
    }
}

//-----------------------------------------------------------------------------
// GetF0CandidateContour() calculates the F0 candidate contour in 1-ch signal.
// Calculation of F0 candidates is carried out in GetF0CandidatesSub().
//-----------------------------------------------------------------------------
static void GetF0CandidateContour(const ZeroCrossings* zero_crossings,
    float boundary_f0, float f0_floor, float f0_ceil,
    const float* temporal_positions, int f0_length, float* f0_candidate)
{
    if (0 == CheckEvent(zero_crossings->number_of_negatives - 2) * CheckEvent(zero_crossings->number_of_positives - 2) * CheckEvent(zero_crossings->number_of_peaks - 2) * CheckEvent(zero_crossings->number_of_dips - 2)) {
        memset(f0_candidate, 0, f0_length * sizeof(*f0_candidate));
        return;
    }

    float* interpolated_f0_set[4];
    for (int i = 0; i < 4; ++i)
        interpolated_f0_set[i] = new float[f0_length];

    interp1(zero_crossings->negative_interval_locations,
        zero_crossings->negative_intervals,
        zero_crossings->number_of_negatives,
        temporal_positions, f0_length, interpolated_f0_set[0]);
    interp1(zero_crossings->positive_interval_locations,
        zero_crossings->positive_intervals,
        zero_crossings->number_of_positives,
        temporal_positions, f0_length, interpolated_f0_set[1]);
    interp1(zero_crossings->peak_interval_locations,
        zero_crossings->peak_intervals, zero_crossings->number_of_peaks,
        temporal_positions, f0_length, interpolated_f0_set[2]);
    interp1(zero_crossings->dip_interval_locations,
        zero_crossings->dip_intervals, zero_crossings->number_of_dips,
        temporal_positions, f0_length, interpolated_f0_set[3]);

    GetF0CandidateContourSub(interpolated_f0_set, f0_length, f0_floor,
        f0_ceil, boundary_f0, f0_candidate);
    for (int i = 0; i < 4; ++i)
        delete[] interpolated_f0_set[i];
}

//-----------------------------------------------------------------------------
// DestroyZeroCrossings() frees the memory of array in the struct
//-----------------------------------------------------------------------------
static void DestroyZeroCrossings(ZeroCrossings* zero_crossings)
{
    delete[] zero_crossings->negative_interval_locations;
    delete[] zero_crossings->positive_interval_locations;
    delete[] zero_crossings->peak_interval_locations;
    delete[] zero_crossings->dip_interval_locations;
    delete[] zero_crossings->negative_intervals;
    delete[] zero_crossings->positive_intervals;
    delete[] zero_crossings->peak_intervals;
    delete[] zero_crossings->dip_intervals;
}

//-----------------------------------------------------------------------------
// GetF0CandidateFromRawEvent() f0 candidate contour in 1-ch signal
//-----------------------------------------------------------------------------
static void GetF0CandidateFromRawEvent(float boundary_f0, float fs,
    const cmplx* y_spectrum, int y_length, stb_fft_real_plan* fft_plan,
    int fft_size, float f0_floor,
    float f0_ceil, const float* temporal_positions, int f0_length,
    float* f0_candidate)
{
    float* filtered_signal = new float[fft_size];
    GetFilteredSignal(boundary_f0, fft_plan, fft_size, fs, y_spectrum,
        y_length, filtered_signal);

    ZeroCrossings zero_crossings = { 0 };
    GetFourZeroCrossingIntervals(filtered_signal, y_length, fs,
        &zero_crossings);

    GetF0CandidateContour(&zero_crossings, boundary_f0, f0_floor, f0_ceil,
        temporal_positions, f0_length, f0_candidate);

    DestroyZeroCrossings(&zero_crossings);
    delete[] filtered_signal;
}

//-----------------------------------------------------------------------------
// GetRawF0Candidates() calculates f0 candidates in all channels.
//-----------------------------------------------------------------------------
static void GetRawF0Candidates(const float* boundary_f0_list,
    int number_of_bands, float actual_fs, int y_length,
    const float* temporal_positions, int f0_length,
    const cmplx* y_spectrum, stb_fft_real_plan* fft_plan, int fft_size, float f0_floor,
    float f0_ceil, float** raw_f0_candidates)
{
    for (int i = 0; i < number_of_bands; ++i)
        GetF0CandidateFromRawEvent(boundary_f0_list[i], actual_fs, y_spectrum,
            y_length, fft_plan, fft_size, f0_floor, f0_ceil, temporal_positions, f0_length,
            raw_f0_candidates[i]);
}

//-----------------------------------------------------------------------------
// DetectF0CandidatesSub1() calculates VUV areas.
//-----------------------------------------------------------------------------
static int DetectOfficialF0CandidatesSub1(const int* vuv,
    int number_of_channels, int* st, int* ed)
{
    int number_of_voiced_sections = 0;
    int tmp;
    for (int i = 1; i < number_of_channels; ++i) {
        tmp = vuv[i] - vuv[i - 1];
        if (tmp == 1)
            st[number_of_voiced_sections] = i;
        if (tmp == -1)
            ed[number_of_voiced_sections++] = i;
    }

    return number_of_voiced_sections;
}

//-----------------------------------------------------------------------------
// DetectOfficialF0CandidatesSub2() calculates F0 candidates in a frame
//-----------------------------------------------------------------------------
static int DetectOfficialF0CandidatesSub2(const int* vuv,
    const float* const* raw_f0_candidates, int index,
    int number_of_voiced_sections, const int* st, const int* ed,
    int max_candidates, float* f0_list)
{
    int number_of_candidates = 0;
    float tmp_f0;
    for (int i = 0; i < number_of_voiced_sections; ++i) {
        if (ed[i] - st[i] < 10)
            continue;

        tmp_f0 = 0.0;
        for (int j = st[i]; j < ed[i]; ++j)
            tmp_f0 += raw_f0_candidates[j][index];
        tmp_f0 /= (ed[i] - st[i]);
        f0_list[number_of_candidates++] = tmp_f0;
    }

    memset(f0_list + number_of_candidates, 0, (max_candidates - number_of_candidates) * sizeof(*f0_list));
    return number_of_candidates;
}

//-----------------------------------------------------------------------------
// DetectOfficialF0Candidates() detectes F0 candidates from multi-channel
// candidates.
//-----------------------------------------------------------------------------
static int DetectOfficialF0Candidates(const float* const* raw_f0_candidates,
    int number_of_channels, int f0_length, int max_candidates,
    float** f0_candidates)
{
    int number_of_candidates = 0;

    int* vuv = new int[number_of_channels];
    int* st = new int[number_of_channels];
    int* ed = new int[number_of_channels];
    int number_of_voiced_sections;
    for (int i = 0; i < f0_length; ++i) {
        for (int j = 0; j < number_of_channels; ++j)
            vuv[j] = raw_f0_candidates[j][i] > 0 ? 1 : 0;
        vuv[0] = vuv[number_of_channels - 1] = 0;
        number_of_voiced_sections = DetectOfficialF0CandidatesSub1(vuv,
            number_of_channels, st, ed);
        number_of_candidates = MAX(number_of_candidates,
            DetectOfficialF0CandidatesSub2(vuv, raw_f0_candidates, i,
                number_of_voiced_sections, st, ed,
                max_candidates, f0_candidates[i]));
    }

    delete[] vuv;
    delete[] st;
    delete[] ed;
    return number_of_candidates;
}

//-----------------------------------------------------------------------------
// OverlapF0Candidates() spreads the candidates to anteroposterior frames.
//-----------------------------------------------------------------------------
static void OverlapF0Candidates(int f0_length, int number_of_candidates,
    float** f0_candidates)
{
    int n = 3;
    for (int i = 1; i <= n; ++i)
        for (int j = 0; j < number_of_candidates; ++j) {
            for (int k = i; k < f0_length; ++k)
                f0_candidates[k][j + (number_of_candidates * i)] = f0_candidates[k - i][j];
            for (int k = 0; k < f0_length - i; ++k)
                f0_candidates[k][j + (number_of_candidates * (i + n))] = f0_candidates[k + i][j];
        }
}

//-----------------------------------------------------------------------------
// GetBaseIndex() calculates the temporal positions for windowing.
//-----------------------------------------------------------------------------
static void GetBaseIndex(float current_position, const float* base_time,
    int base_time_length, float fs, int* base_index)
{
    // First-aid treatment
    int basic_index = matlab_round((current_position + base_time[0]) * fs + 0.001f);

    for (int i = 0; i < base_time_length; ++i)
        base_index[i] = basic_index + i;
}

//-----------------------------------------------------------------------------
// GetMainWindow() generates the window function.
//-----------------------------------------------------------------------------
static void GetMainWindow(float current_position, const int* base_index,
    int base_time_length, float fs, float window_length_in_time,
    float* main_window)
{
    float tmp = 0.0;
    for (int i = 0; i < base_time_length; ++i) {
        tmp = (base_index[i] - 1.0f) / fs - current_position;
        main_window[i] = 0.42f + 0.5f * cosf(2.0f * kPi * tmp / window_length_in_time) + 0.08f * cosf(4.0f * kPi * tmp / window_length_in_time);
    }
}

//-----------------------------------------------------------------------------
// GetDiffWindow() generates the differentiated window.
// Diff means differential.
//-----------------------------------------------------------------------------
static void GetDiffWindow(const float* main_window, int base_time_length,
    float* diff_window)
{
    diff_window[0] = -main_window[1] / 2.0f;
    for (int i = 1; i < base_time_length - 1; ++i)
        diff_window[i] = -(main_window[i + 1] - main_window[i - 1]) / 2.0f;
    diff_window[base_time_length - 1] = main_window[base_time_length - 2] / 2.0f;
}

//-----------------------------------------------------------------------------
// GetSpectra() calculates two spectra of the waveform windowed by windows
// (main window and diff window).
//-----------------------------------------------------------------------------
static void GetSpectra(const float* x, int x_length, int fft_size,
    const int* base_index, const float* main_window,
    const float* diff_window, int base_time_length,
    const ForwardRealFFT* forward_real_fft, cmplx* main_spectrum,
    cmplx* diff_spectrum)
{
    int* safe_index = new int[base_time_length];

    for (int i = 0; i < base_time_length; ++i)
        safe_index[i] = MAX(0, MIN(x_length - 1, base_index[i] - 1));
    for (int i = 0; i < base_time_length; ++i)
        forward_real_fft->waveform[i] = x[safe_index[i]] * main_window[i];

    memset(forward_real_fft->waveform + base_time_length, 0,
        (fft_size - base_time_length) * sizeof(*forward_real_fft->waveform));

    stb_fft_r2c_exec(forward_real_fft->forward_fft, forward_real_fft->waveform, forward_real_fft->spectrum);
    for (int i = 0; i <= fft_size / 2; ++i) {
        main_spectrum[i].real = forward_real_fft->spectrum[i].real;
        main_spectrum[i].imag = forward_real_fft->spectrum[i].imag;
    }

    for (int i = 0; i < base_time_length; ++i)
        forward_real_fft->waveform[i] = x[safe_index[i]] * diff_window[i];
    for (int i = base_time_length; i < fft_size; ++i)
        forward_real_fft->waveform[i] = 0.0;
    stb_fft_r2c_exec(forward_real_fft->forward_fft, forward_real_fft->waveform, forward_real_fft->spectrum);
    for (int i = 0; i <= fft_size / 2; ++i) {
        diff_spectrum[i].real = forward_real_fft->spectrum[i].real;
        diff_spectrum[i].imag = forward_real_fft->spectrum[i].imag;
    }

    delete[] safe_index;
}

static void FixF0(const float* power_spectrum, const float* numerator_i,
    int fft_size, float fs, float current_f0, int number_of_harmonics,
    float* refined_f0, float* score)
{
    float* amplitude_list = new float[number_of_harmonics];
    float* instantaneous_frequency_list = new float[number_of_harmonics];

    int index;
    for (int i = 0; i < number_of_harmonics; ++i) {
        index = matlab_round(current_f0 * fft_size / fs * (i + 1));
        instantaneous_frequency_list[i] = power_spectrum[index] == 0.0f ? 0.0f : static_cast<float>(index) * fs / fft_size + numerator_i[index] / power_spectrum[index] * fs / 2.0f / kPi;
        amplitude_list[i] = sqrtf(power_spectrum[index]);
    }
    float denominator = 0.0f;
    float numerator = 0.0f;
    *score = 0.0f;
    for (int i = 0; i < number_of_harmonics; ++i) {
        numerator += amplitude_list[i] * instantaneous_frequency_list[i];
        denominator += amplitude_list[i] * (i + 1.0f);
        *score += fabsf((instantaneous_frequency_list[i] / (i + 1.0f) - current_f0) / current_f0);
    }

    *refined_f0 = numerator / (denominator + kMySafeGuardMinimum);
    *score = 1.0f / (*score / number_of_harmonics + kMySafeGuardMinimum);

    delete[] amplitude_list;
    delete[] instantaneous_frequency_list;
}

//-----------------------------------------------------------------------------
// GetMeanF0() calculates the instantaneous frequency.
//-----------------------------------------------------------------------------

static void GetMeanF0(const float* x, int x_length, float fs,
    float current_position, float current_f0, int fft_size,
    float window_length_in_time, const float* base_time,
    int base_time_length, float* refined_f0, float* refined_score)
{
    if (g_ForwardRealFFT.find(fft_size) == g_ForwardRealFFT.end()) {
        g_ForwardRealFFT[fft_size] = { 0 };
        InitializeForwardRealFFT(fft_size, &g_ForwardRealFFT[fft_size]);
    }

    cmplx* main_spectrum = new cmplx[fft_size];
    cmplx* diff_spectrum = new cmplx[fft_size];

    int* base_index = new int[base_time_length];
    float* main_window = new float[base_time_length];
    float* diff_window = new float[base_time_length];

    GetBaseIndex(current_position, base_time, base_time_length, fs, base_index);
    GetMainWindow(current_position, base_index, base_time_length, fs,
        window_length_in_time, main_window);
    GetDiffWindow(main_window, base_time_length, diff_window);

    GetSpectra(x, x_length, fft_size, base_index, main_window, diff_window,
        base_time_length, &g_ForwardRealFFT[fft_size], main_spectrum, diff_spectrum);

    float* power_spectrum = new float[fft_size / 2 + 1];
    float* numerator_i = new float[fft_size / 2 + 1];
    for (int j = 0; j <= fft_size / 2; ++j) {
        numerator_i[j] = main_spectrum[j].real * diff_spectrum[j].imag - main_spectrum[j].imag * diff_spectrum[j].real;
        power_spectrum[j] = main_spectrum[j].real * main_spectrum[j].real + main_spectrum[j].imag * main_spectrum[j].imag;
    }

    int number_of_harmonics = MIN(static_cast<int>(fs / 2.0 / current_f0), 6);
    FixF0(power_spectrum, numerator_i, fft_size, fs, current_f0,
        number_of_harmonics, refined_f0, refined_score);

    delete[] diff_spectrum;
    delete[] diff_window;
    delete[] main_window;
    delete[] base_index;
    delete[] numerator_i;
    delete[] power_spectrum;
    delete[] main_spectrum;
}

//-----------------------------------------------------------------------------
// GetRefinedF0() calculates F0 and its score based on instantaneous frequency.
//-----------------------------------------------------------------------------
static void GetRefinedF0(const float* x, int x_length, float fs,
    float current_position, float current_f0, float f0_floor, float f0_ceil,
    float* refined_f0, float* refined_score)
{
    if (current_f0 <= 0.0f) {
        *refined_f0 = 0.0f;
        *refined_score = 0.0f;
        return;
    }

    int half_window_length = static_cast<int>(1.5 * fs / current_f0 + 1.0f);
    float window_length_in_time = (2.0f * half_window_length + 1.0f) / fs;
    float* base_time = new float[half_window_length * 2 + 1];
    for (int i = 0; i < half_window_length * 2 + 1; i++)
        base_time[i] = (-half_window_length + i) / fs;
    int fft_size = static_cast<int>(powf(2.0f,
        2.0f + static_cast<int>(logf(half_window_length * 2.0f + 1.0f) / kLog2)));

    GetMeanF0(x, x_length, fs, current_position, current_f0, fft_size,
        window_length_in_time, base_time, half_window_length * 2 + 1,
        refined_f0, refined_score);

    if (*refined_f0 < f0_floor || *refined_f0 > f0_ceil || *refined_score < 2.5f) {
        *refined_f0 = 0.0f;
        *refined_score = 0.0f;
    }

    delete[] base_time;
}

//-----------------------------------------------------------------------------
// RefineF0() modifies the F0 by instantaneous frequency.
//-----------------------------------------------------------------------------
static void RefineF0Candidates(const float* x, int x_length, float fs,
    const float* temporal_positions, int f0_length, int max_candidates,
    float f0_floor, float f0_ceil,
    float** refined_f0_candidates, float** f0_scores)
{
    for (int i = 0; i < f0_length; i++)
        for (int j = 0; j < max_candidates; ++j)
            GetRefinedF0(x, x_length, fs, temporal_positions[i],
                refined_f0_candidates[i][j], f0_floor, f0_ceil,
                &refined_f0_candidates[i][j], &f0_scores[i][j]);
    for (int k = 0; k < g_ForwardRealFFT.size(); ++k) {
        DestroyForwardRealFFT(&g_ForwardRealFFT[k]);
    }
}

//-----------------------------------------------------------------------------
// SelectBestF0() obtains the nearlest F0 in reference_f0.
//-----------------------------------------------------------------------------
static float SelectBestF0(float reference_f0, const float* f0_candidates,
    int number_of_candidates, float allowed_range, float* best_error)
{
    float best_f0 = 0.0;
    *best_error = allowed_range;

    float tmp;
    for (int i = 0; i < number_of_candidates; ++i) {
        tmp = fabsf(reference_f0 - f0_candidates[i]) / reference_f0;
        if (tmp > *best_error)
            continue;
        best_f0 = f0_candidates[i];
        *best_error = tmp;
    }

    return best_f0;
}

static void RemoveUnreliableCandidatesSub(int i, int j,
    const float* const* tmp_f0_candidates, int number_of_candidates,
    float** f0_candidates, float** f0_scores)
{
    float reference_f0 = f0_candidates[i][j];
    float error1, error2, min_error;
    float threshold = 0.05f;
    if (reference_f0 == 0)
        return;
    SelectBestF0(reference_f0, tmp_f0_candidates[i + 1],
        number_of_candidates, 1.0f, &error1);
    SelectBestF0(reference_f0, tmp_f0_candidates[i - 1],
        number_of_candidates, 1.0f, &error2);
    min_error = MIN(error1, error2);
    if (min_error <= threshold)
        return;
    f0_candidates[i][j] = 0;
    f0_scores[i][j] = 0;
}

//-----------------------------------------------------------------------------
// RemoveUnreliableCandidates().
//-----------------------------------------------------------------------------
static void RemoveUnreliableCandidates(int f0_length, int number_of_candidates,
    float** f0_candidates, float** f0_scores)
{
    float** tmp_f0_candidates = new float*[f0_length];
    for (int i = 0; i < f0_length; ++i)
        tmp_f0_candidates[i] = new float[number_of_candidates];
    for (int i = 1; i < f0_length - 1; ++i)
        for (int j = 0; j < number_of_candidates; ++j)
            tmp_f0_candidates[i][j] = f0_candidates[i][j];

    for (int i = 1; i < f0_length - 1; ++i)
        for (int j = 0; j < number_of_candidates; ++j)
            RemoveUnreliableCandidatesSub(i, j, tmp_f0_candidates,
                number_of_candidates, f0_candidates, f0_scores);

    for (int i = 0; i < f0_length; ++i)
        delete[] tmp_f0_candidates[i];
    delete[] tmp_f0_candidates;
}

//-----------------------------------------------------------------------------
// SearchF0Base() gets the F0 with the highest score.
//-----------------------------------------------------------------------------
static void SearchF0Base(const float* const* f0_candidates,
    const float* const* f0_scores, int f0_length, int number_of_candidates,
    float* base_f0_contour)
{
    float tmp_best_score;
    for (int i = 0; i < f0_length; ++i) {
        base_f0_contour[i] = tmp_best_score = 0.0;
        for (int j = 0; j < number_of_candidates; ++j)
            if (f0_scores[i][j] > tmp_best_score) {
                base_f0_contour[i] = f0_candidates[i][j];
                tmp_best_score = f0_scores[i][j];
            }
    }
}

//-----------------------------------------------------------------------------
// Step 1: Rapid change of F0 contour is replaced by 0.
//-----------------------------------------------------------------------------
static void FixStep1(const float* f0_base, int f0_length,
    float allowed_range, float* f0_step1)
{
    f0_step1[0] = f0_step1[1] = 0.0f;
    float reference_f0;
    for (int i = 2; i < f0_length; ++i) {
        if (f0_base[i] == 0.0)
            continue;
        reference_f0 = f0_base[i - 1] * 2 - f0_base[i - 2];
        f0_step1[i] = fabsf((f0_base[i] - reference_f0) / reference_f0) > allowed_range && fabsf((f0_base[i] - f0_base[i - 1])) / f0_base[i - 1] > allowed_range ? 0.0f : f0_base[i];
    }
}

//-----------------------------------------------------------------------------
// GetBoundaryList() detects boundaries between voiced and unvoiced sections.
//-----------------------------------------------------------------------------
static int GetBoundaryList(const float* f0, int f0_length,
    int* boundary_list)
{
    int number_of_boundaries = 0;
    int* vuv = new int[f0_length];
    for (int i = 0; i < f0_length; ++i)
        vuv[i] = f0[i] > 0 ? 1 : 0;
    vuv[0] = vuv[f0_length - 1] = 0;

    for (int i = 1; i < f0_length; ++i)
        if (vuv[i] - vuv[i - 1] != 0) {
            boundary_list[number_of_boundaries] = i - number_of_boundaries % 2;
            number_of_boundaries++;
        }

    delete[] vuv;
    return number_of_boundaries;
}

//-----------------------------------------------------------------------------
// Step 2: Voiced sections with a short period are removed.
//-----------------------------------------------------------------------------
static void FixStep2(const float* f0_step1, int f0_length,
    int voice_range_minimum, float* f0_step2)
{
    for (int i = 0; i < f0_length; ++i)
        f0_step2[i] = f0_step1[i];
    int* boundary_list = new int[f0_length];
    int number_of_boundaries = GetBoundaryList(f0_step1, f0_length, boundary_list);

    for (int i = 0; i < number_of_boundaries / 2; ++i) {
        if (boundary_list[i * 2 + 1] - boundary_list[i * 2] >= voice_range_minimum)
            continue;
        for (int j = boundary_list[i * 2]; j <= boundary_list[(i * 2) + 1]; ++j)
            f0_step2[j] = 0.0;
    }
    delete[] boundary_list;
}

//-----------------------------------------------------------------------------
// GetMultiChannelF0() separates each voiced section into independent channel.
//-----------------------------------------------------------------------------
static void GetMultiChannelF0(const float* f0, int f0_length,
    const int* boundary_list, int number_of_boundaries,
    float** multi_channel_f0)
{
    for (int i = 0; i < number_of_boundaries / 2; ++i) {
        for (int j = 0; j < boundary_list[i * 2]; ++j)
            multi_channel_f0[i][j] = 0.0;
        for (int j = boundary_list[i * 2]; j <= boundary_list[i * 2 + 1]; ++j)
            multi_channel_f0[i][j] = f0[j];
        for (int j = boundary_list[i * 2 + 1] + 1; j < f0_length; ++j)
            multi_channel_f0[i][j] = 0.0;
    }
}

//-----------------------------------------------------------------------------
// abs() often causes bugs, an original function is used.
//-----------------------------------------------------------------------------

static inline int Abs(int x)
{
    return (x ^ (x >> 31)) - (x >> 31);
}

//-----------------------------------------------------------------------------
// ExtendF0() : The Hand erasing the Space.
// The subfunction of Extend().
//-----------------------------------------------------------------------------
static int ExtendF0(const float* f0, int f0_length, int origin,
    int last_point, int shift, const float* const* f0_candidates,
    int number_of_candidates, float allowed_range, float* extended_f0)
{
    int threshold = 4;
    float tmp_f0 = extended_f0[origin];
    int shifted_origin = origin;

    int distance = Abs(last_point - origin);
    int* index_list = new int[distance + 1];
    for (int i = 0; i <= distance; ++i)
        index_list[i] = origin + shift * i;

    int count = 0;
    float dammy;
    for (int i = 0; i <= distance; ++i) {
        extended_f0[index_list[i] + shift] = SelectBestF0(tmp_f0, f0_candidates[index_list[i] + shift],
            number_of_candidates, allowed_range, &dammy);
        if (extended_f0[index_list[i] + shift] == 0.0) {
            count++;
        } else {
            tmp_f0 = extended_f0[index_list[i] + shift];
            count = 0;
            shifted_origin = index_list[i] + shift;
        }
        if (count == threshold)
            break;
    }

    delete[] index_list;
    return shifted_origin;
}

//-----------------------------------------------------------------------------
// Swap the f0 contour and boundary.
// It is used in ExtendSub() and MergeF0();
//-----------------------------------------------------------------------------
static void Swap(int index1, int index2, float** f0, int* boundary)
{
    float* tmp_pointer;
    int tmp_index;
    tmp_pointer = f0[index1];
    f0[index1] = f0[index2];
    f0[index2] = tmp_pointer;
    tmp_index = boundary[index1 * 2];
    boundary[index1 * 2] = boundary[index2 * 2];
    boundary[index2 * 2] = tmp_index;
    tmp_index = boundary[index1 * 2 + 1];
    boundary[index1 * 2 + 1] = boundary[index2 * 2 + 1];
    boundary[index2 * 2 + 1] = tmp_index;
}

static int ExtendSub(const float* const* extended_f0,
    const int* boundary_list, int number_of_sections,
    float** selected_extended_f0, int* selected_boundary_list)
{
    float threshold = 2200.0;
    int count = 0;
    float mean_f0 = 0.0;
    int st, ed;
    for (int i = 0; i < number_of_sections; ++i) {
        st = boundary_list[i * 2];
        ed = boundary_list[i * 2 + 1];
        for (int j = st; j < ed; ++j)
            mean_f0 += extended_f0[i][j];
        mean_f0 /= ed - st;
        if (threshold / mean_f0 < ed - st)
            Swap(count++, i, selected_extended_f0, selected_boundary_list);
    }
    return count;
}

//-----------------------------------------------------------------------------
// Extend() : The Hand erasing the Space.
//-----------------------------------------------------------------------------
static int Extend(const float* const* multi_channel_f0,
    int number_of_sections, int f0_length, const int* boundary_list,
    const float* const* f0_candidates, int number_of_candidates,
    float allowed_range, float** extended_f0, int* shifted_boundary_list)
{
    int threshold = 100;
    for (int i = 0; i < number_of_sections; ++i) {
        shifted_boundary_list[i * 2 + 1] = ExtendF0(multi_channel_f0[i],
            f0_length, boundary_list[i * 2 + 1],
            MIN(f0_length - 2, boundary_list[i * 2 + 1] + threshold),
            1,
            f0_candidates, number_of_candidates, allowed_range,
            extended_f0[i]);
        shifted_boundary_list[i * 2] = ExtendF0(multi_channel_f0[i], f0_length,
            boundary_list[i * 2], MAX(1, boundary_list[i * 2] - threshold),
            -1,
            f0_candidates, number_of_candidates, allowed_range, extended_f0[i]);
    }

    return ExtendSub(multi_channel_f0, shifted_boundary_list,
        number_of_sections, extended_f0, shifted_boundary_list);
}

//-----------------------------------------------------------------------------
// Indices are sorted.
//-----------------------------------------------------------------------------
static void MakeSortedOrder(const int* boundary_list, int number_of_sections,
    int* order)
{
    for (int i = 0; i < number_of_sections; ++i)
        order[i] = i;
    int tmp;
    for (int i = 1; i < number_of_sections; ++i)
        for (int j = i - 1; j >= 0; --j)
            if (boundary_list[order[j] * 2] > boundary_list[order[i] * 2]) {
                tmp = order[i];
                order[i] = order[j];
                order[j] = tmp;
            } else {
                break;
            }
}

//-----------------------------------------------------------------------------
// Serach the highest score with the candidate F0.
//-----------------------------------------------------------------------------
static float SearchScore(float f0, const float* f0_candidates,
    const float* f0_scores, int number_of_candidates)
{
    float score = 0.0;
    for (int i = 0; i < number_of_candidates; ++i)
        if (f0 == f0_candidates[i] && score < f0_scores[i])
            score = f0_scores[i];
    return score;
}

//-----------------------------------------------------------------------------
// Subfunction of MergeF0()
//-----------------------------------------------------------------------------
static int MergeF0Sub(const float* f0_1, int f0_length, int st1, int ed1,
    const float* f0_2, int st2, int ed2, const float* const* f0_candidates,
    const float* const* f0_scores, int number_of_candidates,
    float* merged_f0)
{
    if (st1 <= st2 && ed1 >= ed2)
        return ed1;

    float score1 = 0.0;
    float score2 = 0.0;
    for (int i = st2; i <= ed1; ++i) {
        score1 += SearchScore(f0_1[i], f0_candidates[i], f0_scores[i],
            number_of_candidates);
        score2 += SearchScore(f0_2[i], f0_candidates[i], f0_scores[i],
            number_of_candidates);
    }
    if (score1 > score2)
        for (int i = ed1; i <= ed2; ++i)
            merged_f0[i] = f0_2[i];
    else
        for (int i = st2; i <= ed2; ++i)
            merged_f0[i] = f0_2[i];

    return ed2;
}

//-----------------------------------------------------------------------------
// Overlapped F0 contours are merged by the likability score.
//-----------------------------------------------------------------------------
static void MergeF0(const float* const* multi_channel_f0, int* boundary_list,
    int number_of_channels, int f0_length, const float* const* f0_candidates,
    const float* const* f0_scores, int number_of_candidates,
    float* merged_f0)
{
    int* order = new int[number_of_channels];
    MakeSortedOrder(boundary_list, number_of_channels, order);

    for (int i = 0; i < f0_length; ++i)
        merged_f0[i] = multi_channel_f0[0][i];

    for (int i = 1; i < number_of_channels; ++i)
        if (boundary_list[order[i] * 2] - boundary_list[1] > 0) {
            for (int j = boundary_list[order[i] * 2];
                 j <= boundary_list[order[i] * 2 + 1]; ++j)
                merged_f0[j] = multi_channel_f0[order[i]][j];
            boundary_list[0] = boundary_list[order[i] * 2];
            boundary_list[1] = boundary_list[order[i] * 2 + 1];
        } else {
            boundary_list[1] = MergeF0Sub(merged_f0, f0_length, boundary_list[0], boundary_list[1],
                multi_channel_f0[order[i]], boundary_list[order[i] * 2],
                boundary_list[order[i] * 2 + 1], f0_candidates, f0_scores,
                number_of_candidates, merged_f0);
        }

    delete[] order;
}

//-----------------------------------------------------------------------------
// Step 3: Voiced sections are extended based on the continuity of F0 contour
//-----------------------------------------------------------------------------
static void FixStep3(const float* f0_step2, int f0_length,
    int number_of_candidates, const float* const* f0_candidates,
    float allowed_range, const float* const* f0_scores, float* f0_step3)
{
    for (int i = 0; i < f0_length; ++i)
        f0_step3[i] = f0_step2[i];
    int* boundary_list = new int[f0_length];
    int number_of_boundaries = GetBoundaryList(f0_step2, f0_length, boundary_list);

    float** multi_channel_f0 = new float*[number_of_boundaries / 2];
    for (int i = 0; i < number_of_boundaries / 2; ++i)
        multi_channel_f0[i] = new float[f0_length];
    GetMultiChannelF0(f0_step2, f0_length, boundary_list, number_of_boundaries,
        multi_channel_f0);

    int number_of_channels = Extend(multi_channel_f0, number_of_boundaries / 2, f0_length,
        boundary_list, f0_candidates, number_of_candidates, allowed_range,
        multi_channel_f0, boundary_list);

    if (number_of_channels != 0)
        MergeF0(multi_channel_f0, boundary_list, number_of_channels, f0_length,
            f0_candidates, f0_scores, number_of_candidates, f0_step3);

    for (int i = 0; i < number_of_boundaries / 2; ++i)
        delete[] multi_channel_f0[i];
    delete[] multi_channel_f0;
    delete[] boundary_list;
}

//-----------------------------------------------------------------------------
// Step 4: F0s in short unvoiced section are faked
//-----------------------------------------------------------------------------
static void FixStep4(const float* f0_step3, int f0_length, int threshold,
    float* f0_step4)
{
    for (int i = 0; i < f0_length; ++i)
        f0_step4[i] = f0_step3[i];
    int* boundary_list = new int[f0_length];
    int number_of_boundaries = GetBoundaryList(f0_step3, f0_length, boundary_list);

    int distance;
    float tmp0, tmp1, coefficient;
    int count;
    for (int i = 0; i < number_of_boundaries / 2 - 1; ++i) {
        distance = boundary_list[(i + 1) * 2] - boundary_list[i * 2 + 1] - 1;
        if (distance >= threshold)
            continue;
        tmp0 = f0_step3[boundary_list[i * 2 + 1]] + 1;
        tmp1 = f0_step3[boundary_list[(i + 1) * 2]] - 1;
        coefficient = (tmp1 - tmp0) / (distance + 1.0f);
        count = 1;
        for (int j = boundary_list[i * 2 + 1] + 1;
             j <= boundary_list[(i + 1) * 2] - 1; ++j)
            f0_step4[j] = tmp0 + coefficient * count++;
    }
    delete[] boundary_list;
}

//-----------------------------------------------------------------------------
// FixF0Contour() obtains the likely F0 contour.
//-----------------------------------------------------------------------------
static void FixF0Contour(const float* const* f0_candidates,
    const float* const* f0_scores, int f0_length, int number_of_candidates,
    float* best_f0_contour)
{
    float* tmp_f0_contour1 = new float[f0_length];
    float* tmp_f0_contour2 = new float[f0_length];

    // These parameters are optimized by speech databases.
    SearchF0Base(f0_candidates, f0_scores, f0_length,
        number_of_candidates, tmp_f0_contour1);
    FixStep1(tmp_f0_contour1, f0_length, 0.008f, tmp_f0_contour2);
    FixStep2(tmp_f0_contour2, f0_length, 6, tmp_f0_contour1);
    FixStep3(tmp_f0_contour1, f0_length, number_of_candidates, f0_candidates,
        0.18f, f0_scores, tmp_f0_contour2);
    FixStep4(tmp_f0_contour2, f0_length, 9, best_f0_contour);

    delete[] tmp_f0_contour1;
    delete[] tmp_f0_contour2;
}

//-----------------------------------------------------------------------------
// This function uses zero-lag Butterworth filter.
//-----------------------------------------------------------------------------
static void FilteringF0(const float* a, const float* b, float* x,
    int x_length, int st, int ed, float* y)
{
    float w[2] = { 0.0, 0.0 };
    float wt;
    float* tmp_x = new float[x_length];

    for (int i = 0; i < st; ++i)
        x[i] = x[st];
    for (int i = ed + 1; i < x_length; ++i)
        x[i] = x[ed];

    for (int i = 0; i < x_length; ++i) {
        wt = x[i] + a[0] * w[0] + a[1] * w[1];
        tmp_x[x_length - i - 1] = b[0] * wt + b[1] * w[0] + b[0] * w[1];
        w[1] = w[0];
        w[0] = wt;
    }

    w[0] = w[1] = 0.0;
    for (int i = 0; i < x_length; ++i) {
        wt = tmp_x[i] + a[0] * w[0] + a[1] * w[1];
        y[x_length - i - 1] = b[0] * wt + b[1] * w[0] + b[0] * w[1];
        w[1] = w[0];
        w[0] = wt;
    }

    delete[] tmp_x;
}

//-----------------------------------------------------------------------------
// SmoothF0Contour() uses the zero-lag Butterworth filter for smoothing.
//-----------------------------------------------------------------------------
static void SmoothF0Contour(const float* f0, int f0_length,
    float* smoothed_f0)
{
    const float b[2] = { 0.0078202080334971724f, 0.015640416066994345f };
    const float a[2] = { 1.7347257688092754f, -0.76600660094326412f };
    int lag = 300;
    int new_f0_length = f0_length + lag * 2;
    float* f0_contour = new float[new_f0_length];
    for (int i = 0; i < lag; ++i)
        f0_contour[i] = 0.0f;
    for (int i = lag; i < lag + f0_length; ++i)
        f0_contour[i] = f0[i - lag];
    for (int i = lag + f0_length; i < new_f0_length; ++i)
        f0_contour[i] = 0.0f;

    int* boundary_list = new int[new_f0_length];
    int number_of_boundaries = GetBoundaryList(f0_contour, new_f0_length, boundary_list);
    float** multi_channel_f0 = new float*[number_of_boundaries / 2];
    for (int i = 0; i < number_of_boundaries / 2; ++i)
        multi_channel_f0[i] = new float[new_f0_length];
    GetMultiChannelF0(f0_contour, new_f0_length, boundary_list,
        number_of_boundaries, multi_channel_f0);

    for (int i = 0; i < number_of_boundaries / 2; ++i) {
        FilteringF0(a, b, multi_channel_f0[i], new_f0_length,
            boundary_list[i * 2], boundary_list[i * 2 + 1], f0_contour);
        for (int j = boundary_list[i * 2]; j <= boundary_list[i * 2 + 1]; ++j)
            smoothed_f0[j - lag] = f0_contour[j];
    }

    for (int i = 0; i < number_of_boundaries / 2; ++i)
        delete[] multi_channel_f0[i];
    delete[] multi_channel_f0;
    delete[] f0_contour;
    delete[] boundary_list;
}

//-----------------------------------------------------------------------------
// HarvestGeneralBodySub() is the subfunction of HarvestGeneralBody()
//-----------------------------------------------------------------------------
static int HarvestGeneralBodySub(const float* boundary_f0_list,
    int number_of_channels, int f0_length, float actual_fs, int y_length,
    const float* temporal_positions, const cmplx* y_spectrum,
    stb_fft_real_plan* fft_plan,
    int fft_size, float f0_floor, float f0_ceil, int max_candidates,
    float** f0_candidates)
{
    float** raw_f0_candidates = new float*[number_of_channels];
    for (int i = 0; i < number_of_channels; ++i)
        raw_f0_candidates[i] = new float[f0_length];

    GetRawF0Candidates(boundary_f0_list, number_of_channels,
        actual_fs, y_length, temporal_positions, f0_length, y_spectrum, fft_plan,
        fft_size, f0_floor, f0_ceil, raw_f0_candidates);

    int number_of_candidates = DetectOfficialF0Candidates(raw_f0_candidates,
        number_of_channels, f0_length, max_candidates,
        f0_candidates);

    OverlapF0Candidates(f0_length, number_of_candidates, f0_candidates);

    for (int i = 0; i < number_of_channels; ++i)
        delete[] raw_f0_candidates[i];
    delete[] raw_f0_candidates;
    return number_of_candidates;
}

//-----------------------------------------------------------------------------
// HarvestGeneralBody() estimates the F0 contour based on Harvest.
//-----------------------------------------------------------------------------
static void HarvestGeneralBody(const float* x, int x_length, int fs,
    int frame_period, float f0_floor, float f0_ceil,
    float channels_in_octave, int speed, float* temporal_positions,
    float* f0)
{
    float adjusted_f0_floor = f0_floor * 0.9f;
    float adjusted_f0_ceil = f0_ceil * 1.1f;
    int number_of_channels = 1 + static_cast<int>(logf(adjusted_f0_ceil / adjusted_f0_floor) / kLog2 * channels_in_octave);
    float* boundary_f0_list = new float[number_of_channels];
    for (int i = 0; i < number_of_channels; ++i)
        boundary_f0_list[i] = adjusted_f0_floor * powf(2.0f, (i + 1) / channels_in_octave);

    // normalization
    int decimation_ratio = MAX(MIN(speed, 12), 1);
    int y_length = static_cast<int>(ceil(static_cast<float>(x_length) / decimation_ratio));
    float actual_fs = static_cast<float>(fs) / decimation_ratio;
    int fft_size = GetSuitableFFTSize(y_length + 5 + 2 * static_cast<int>(2.0f * actual_fs / boundary_f0_list[0]));

    stb_fft_real_plan* fft_plan = stb_fft_real_plan_dft_1d(fft_size);
    if (fft_plan == NULL)
        return;

    // Calculation of the spectrum used for the f0 estimation
    float* y = new float[fft_size];
    cmplx* y_spectrum = new cmplx[fft_size];
    GetWaveformAndSpectrum(x, x_length, y_length, actual_fs, fft_plan, fft_size,
        decimation_ratio, y, y_spectrum);

    int f0_length = GetSamplesForHarvest(fs, x_length, (float)frame_period);
    for (int i = 0; i < f0_length; ++i) {
        temporal_positions[i] = i * frame_period / 1000.0f;
    }
    memset(f0, 0, f0_length * sizeof(*f0));
    int overlap_parameter = 7;
    int max_candidates = matlab_round(number_of_channels / 10.0f) * overlap_parameter;
    float** f0_candidates = new float*[f0_length];
    float** f0_candidates_score = new float*[f0_length];
    for (int i = 0; i < f0_length; ++i) {
        f0_candidates[i] = new float[max_candidates];
        f0_candidates_score[i] = new float[max_candidates];
    }

    int number_of_candidates = HarvestGeneralBodySub(boundary_f0_list,
                                   number_of_channels, f0_length, actual_fs, y_length,
                                   temporal_positions,
                                   y_spectrum, fft_plan, fft_size, f0_floor, f0_ceil,
                                   max_candidates,
                                   f0_candidates)
        * overlap_parameter;

    RefineF0Candidates(y, y_length, actual_fs, temporal_positions, f0_length,
        number_of_candidates, f0_floor, f0_ceil, f0_candidates,
        f0_candidates_score);
    RemoveUnreliableCandidates(f0_length, number_of_candidates,
        f0_candidates, f0_candidates_score);

    float* best_f0_contour = new float[f0_length];
    FixF0Contour(f0_candidates, f0_candidates_score, f0_length,
        number_of_candidates, best_f0_contour);
    SmoothF0Contour(best_f0_contour, f0_length, f0);

    delete[] y;
    delete[] best_f0_contour;
    delete[] y_spectrum;
    for (int i = 0; i < f0_length; ++i) {
        delete[] f0_candidates[i];
        delete[] f0_candidates_score[i];
    }
    delete[] f0_candidates;
    delete[] f0_candidates_score;
    delete[] boundary_f0_list;

    free(fft_plan);
}
} // namespace

int GetSamplesForHarvest(int fs, int x_length, float frame_period)
{
    return static_cast<int>(1000.0 * x_length / fs / frame_period) + 1;
}

void Harvest(const float* x, int x_length, int fs,
    const HarvestOption* option, float* temporal_positions, float* f0)
{
    // Several parameters will be controllable for debug.
    float target_fs = 8000.0;
    int dimension_ratio = matlab_round(fs / target_fs);
    float channels_in_octave = 40;

    if (option->frame_period == 1.0) {
        HarvestGeneralBody(x, x_length, fs, 1, option->f0_floor,
            option->f0_ceil, channels_in_octave, dimension_ratio,
            temporal_positions, f0);
        return;
    }

    int basic_frame_period = 1;
    int basic_f0_length = GetSamplesForHarvest(fs, x_length, (float)basic_frame_period);
    float* basic_f0 = new float[basic_f0_length];
    float* basic_temporal_positions = new float[basic_f0_length];
    HarvestGeneralBody(x, x_length, fs, basic_frame_period, option->f0_floor,
        option->f0_ceil, channels_in_octave, dimension_ratio,
        basic_temporal_positions, basic_f0);

    int f0_length = GetSamplesForHarvest(fs, x_length, option->frame_period);
    for (int i = 0; i < f0_length; ++i) {
        temporal_positions[i] = i * option->frame_period / 1000.0f;
        f0[i] = basic_f0[MIN(basic_f0_length - 1,
            matlab_round(temporal_positions[i] * 1000.0f))];
    }

    delete[] basic_f0;
    delete[] basic_temporal_positions;
}

void InitializeHarvestOption(HarvestOption* option)
{
    // You can change default parameters.
    option->f0_ceil = kCeilF0;
    option->f0_floor = kFloorF0;
    option->frame_period = 5;
}