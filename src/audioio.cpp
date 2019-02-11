#ifndef MAX
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#endif

#include "./audioio.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define DR_MP3_IMPLEMENTATION

#include "dr_mp3.h"

#define DR_WAV_IMPLEMENTATION

#include "dr_wav.h"

#if (defined(__WIN32__) || defined(_WIN32)) && !defined(__MINGW32__)
#pragma warning(disable : 4996)
#endif

void wavwrite(const float* input, int totalSampleCount, int sampleRate,
    const char* filename)
{
    float* output = (float*)malloc(totalSampleCount * sizeof(float));
    if (output == NULL) {
        fprintf(stderr, "ERROR\n");
        exit(1);
    }
    for (int i = 0; i < totalSampleCount; ++i) {
        output[i] = input[i];
    }
    drwav_data_format format;
    format.container = drwav_container_riff; // <-- drwav_container_riff = normal WAV files, drwav_container_w64 = Sony Wave64.
    format.channels = 1;
    format.sampleRate = (drwav_uint32)sampleRate;
    format.bitsPerSample = sizeof(*output) * 8;
    format.format = DR_WAVE_FORMAT_IEEE_FLOAT;
    drwav* pWav = drwav_open_file_write(filename, &format);
    if (pWav) {
        drwav_uint64 samplesWritten = drwav_write(pWav, totalSampleCount, output);
        drwav_uninit(pWav);
        if (samplesWritten != totalSampleCount) {
            fprintf(stderr, "ERROR\n");
            exit(1);
        }
    }
    free(output);
}

float* wavread(const char* filename, int* fs, int* sampleCount)
{
    unsigned int channels = 1;
    unsigned int sampleRate = 0;
    drwav_uint64 totalSampleCount = 0;
    float* input = drwav_open_file_and_read_pcm_frames_f32(filename, &channels, &sampleRate, &totalSampleCount);
    if (input == NULL) {
        drmp3_config pConfig;
        input = drmp3_open_file_and_read_f32(filename, &pConfig, &totalSampleCount);
        if (input != NULL) {
            channels = pConfig.outputChannels;
            sampleRate = pConfig.outputSampleRate;
        }
    }
    if (input == NULL) {
        fprintf(stderr, "read file [%s] error.\n", filename);
        exit(1);
    }
    *fs = sampleRate;
    *sampleCount = (int)totalSampleCount;
    if (channels != 1) {
        float norm = 1.0f / channels;
        float* in = input;
        for (int i = 0; i < totalSampleCount; ++i) {
            float out = 0;
            for (unsigned int c = 0; c < channels; ++c) {
                out += in[c];
            }
            input[i] = out * norm;
            in += channels;
        }
    }
    return input;
}