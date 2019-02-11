#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if (defined(__WIN32__) || defined(_WIN32)) && !defined(__MINGW32__)
#include <conio.h>
#include <windows.h>
#pragma comment(lib, "winmm.lib")
#pragma warning(disable : 4996)
#endif
#if (defined(__linux__) || defined(__CYGWIN__) || defined(__APPLE__))

#include <stdint.h>
#include <sys/time.h>

#endif

// For .wav input/output functions.
#include "audioio.h"

// WORLD core functions.
// Note: win.sln uses an option in Additional Include Directories.
// To compile the program, the option "-I $(SolutionDir)..\src" was set.
#include "harvest.h"
#include "matlabfunctions.h"

#if (defined(__linux__) || defined(__CYGWIN__) || defined(__APPLE__))
// Linux porting section: implement timeGetTime() by gettimeofday(),
#ifndef DWORD
#define DWORD uint32_t
#endif

DWORD timeGetTime() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    DWORD ret = static_cast<DWORD>(tv.tv_usec / 1000 + tv.tv_sec * 1000);
    return ret;
}

#endif

//-----------------------------------------------------------------------------
// struct for WORLD
// This struct is an option.
// Users are NOT forced to use this struct.
//-----------------------------------------------------------------------------
typedef struct {
    float frame_period;
    int fs;
    float *f0;
    float *time_axis;
    int f0_length;
    int fft_size;
} HarvestParameters;

namespace {
    void DisplayInformation(int fs, int input_length) {
        printf("File information\n");
        printf("Sampling : %d Hz\n", fs);
        printf("Length %d [sample]\n", input_length);
        printf("Length %f [sec]\n", static_cast<float>(input_length) / fs);
    }

#define DEBUG 1

    void F0EstimationHarvest(float *input, int input_length,
                             HarvestParameters *harvestParameters) {
        HarvestOption option = {0};
        InitializeHarvestOption(&option);

        // You can change the frame period.
        // But the estimation is carried out with 1-ms frame shift.
        option.frame_period = harvestParameters->frame_period;

        // You can set the f0_floor below include::kFloorF0.
        option.f0_floor = 40.0;

        // Parameters setting and memory allocation.
        harvestParameters->f0_length = GetSamplesForHarvest(harvestParameters->fs,
                                                            input_length, harvestParameters->frame_period);
        harvestParameters->f0 = new float[harvestParameters->f0_length];
        harvestParameters->time_axis = new float[harvestParameters->f0_length];

        printf("\nAnalysis\n");
        DWORD elapsed_time = timeGetTime();
        Harvest(input, input_length, harvestParameters->fs, &option,
                harvestParameters->time_axis, harvestParameters->f0);
#if DEBUG
        printf("f0 length %d: \n", harvestParameters->f0_length);
        for (int i = 0; i < harvestParameters->f0_length; ++i) {
            printf("f0 %f: time_axis %f \n", harvestParameters->f0[i], harvestParameters->time_axis[i]);
        }
#endif
        printf("Harvest: %d [msec]\n", timeGetTime() - elapsed_time);
    }

    void DestroyMemory(HarvestParameters *world_parameters) {
        delete[] world_parameters->time_axis;
        delete[] world_parameters->f0;
    }
} // namespace

//-----------------------------------------------------------------------------
int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("usage : %s input.wav \n", argv[0]);
        return -2;
    }

    const char *inFile = argv[1];
    int input_length = 0;
    int fs;
    float *input = wavread(inFile, &fs, &input_length);
    if (input == NULL) {
        printf("error: The file is not .wav format.\n");
        return -1;
    }
    DisplayInformation(fs, input_length);

    //---------------------------------------------------------------------------
    // Analysis part
    //---------------------------------------------------------------------------
    HarvestParameters parameters = {0};
    // You must set fs and frame_period before analysis/synthesis.
    parameters.fs = fs;
    // 5.0 ms is the default value.
    parameters.frame_period = 5.0;

    // F0 estimation
    F0EstimationHarvest(input, input_length, &parameters);

    free(input);
    DestroyMemory(&parameters);

    printf("complete.\n");
    return 0;
}