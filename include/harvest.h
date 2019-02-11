
#ifndef _HARVEST_H_
#define _HARVEST_H_

#include "common.h"

BEGIN_C_DECLS

//-----------------------------------------------------------------------------
// Struct for Harvest
//-----------------------------------------------------------------------------
typedef struct {
    float f0_floor;
    float f0_ceil;
    float frame_period;
} HarvestOption;

//-----------------------------------------------------------------------------
// Harvest
//
// Input:
//   x                    : Input signal
//   x_length             : Length of x
//   fs                   : Sampling frequency
//   option               : Struct to order the parameter for Harvest
//
// Output:
//   temporal_positions   : Temporal positions.
//   f0                   : F0 contour.
//-----------------------------------------------------------------------------
void Harvest(const float* x, int x_length, int fs,
    const HarvestOption* option, float* temporal_positions, float* f0);

//-----------------------------------------------------------------------------
// InitializeHarvestOption allocates the memory to the struct and sets the
// default parameters.
//
// Output:
//   option   : Struct for the optional parameter.
//-----------------------------------------------------------------------------
void InitializeHarvestOption(HarvestOption* option);

//-----------------------------------------------------------------------------
// GetSamplesForHarvest() calculates the number of samples required for
// Harvest().
//
// Input:
//   fs             : Sampling frequency [Hz]
//   x_length       : Length of the input signal [Sample]
//   frame_period   : Frame shift [msec]
//
// Output:
//   The number of samples required to store the results of Harvest().
//-----------------------------------------------------------------------------
int GetSamplesForHarvest(int fs, int x_length, float frame_period);

END_C_DECLS

#endif // _HARVEST_H_
