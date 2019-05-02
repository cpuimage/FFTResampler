
#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define DR_WAV_IMPLEMENTATION

#include "dr_wav.h"

#define DR_MP3_IMPLEMENTATION


#include "dr_mp3.h"

#include "timing.h"


#define STB_FFT_IMPLEMENTAION

#include "stb_fft.h"

#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif


void wavWrite_f32(char *filename, float *buffer, int sampleRate, uint32_t totalSampleCount, uint32_t channels)
{
    drwav_data_format format;
    format.container = drwav_container_riff;
    format.format = DR_WAVE_FORMAT_IEEE_FLOAT;
    format.channels = channels;
    format.sampleRate = (drwav_uint32) sampleRate;
    format.bitsPerSample = 32;
    drwav *pWav = drwav_open_file_write(filename, &format);
    if (pWav) {
        drwav_uint64 samplesWritten = drwav_write(pWav, totalSampleCount, buffer);
        drwav_uninit(pWav);
        if (samplesWritten != totalSampleCount) {
            fprintf(stderr, "write file [%s] error.\n", filename);
            exit(1);
        }
    }
}

float *wavRead_f32(const char *filename, uint32_t *sampleRate, uint64_t *sampleCount, uint32_t *channels)
{
    drwav_uint64 totalSampleCount = 0;
    float *input = drwav_open_file_and_read_pcm_frames_f32(filename, channels, sampleRate, &totalSampleCount);
    if (input == NULL) {
        drmp3_config pConfig;
        input = drmp3_open_file_and_read_f32(filename, &pConfig, &totalSampleCount);
        if (input != NULL) {
            *channels = pConfig.outputChannels;
            *sampleRate = pConfig.outputSampleRate;
        }
    }
    if (input == NULL) {
        fprintf(stderr, "read file [%s] error.\n", filename);
        exit(1);
    }
    *sampleCount = totalSampleCount * (*channels);
    return input;
}

void splitpath(const char *path, char *drv, char *dir, char *name, char *ext)
{
    const char *end;
    const char *p;
    const char *s;
    if (path[0] && path[1] == ':') {
        if (drv) {
            *drv++ = *path++;
            *drv++ = *path++;
            *drv = '\0';
        }
    }
    else if (drv)
        *drv = '\0';
    for (end = path; *end && *end != ':';)
        end++;
    for (p = end; p > path && *--p != '\\' && *p != '/';)
        if (*p == '.') {
            end = p;
            break;
        }
    if (ext)
        for (s = end; (*ext = *s++);)
            ext++;
    for (p = end; p > path;)
        if (*--p == '\\' || *p == '/') {
            p++;
            break;
        }
    if (name) {
        for (s = p; s < end;)
            *name++ = *s++;
        *name = '\0';
    }
    if (dir) {
        for (s = path; s < p;)
            *dir++ = *s++;
        *dir = '\0';
    }
}

typedef struct
{
    int inFrameSize;
    int inWindowSize;
    int inSampleRate;
    float *inWindowing;
    stb_fft_real_plan *inPlan;
    int outFrameSize;
    int outWindowSize;
    int outSampleRate;
    float *outWindowing;
    stb_fft_real_plan *outPlan;
    float *inFifo;
    float *synthesisMem;
    cmplx *samples;
    int pos;
} FFT_Resampler_Handle;

void FFT_Resampler_Free(FFT_Resampler_Handle *handle)
{
    if (handle) {
        if (handle->inFifo) {
            free(handle->inFifo);
            handle->inFifo = NULL;
        }

        if (handle->inPlan) {
            free(handle->inPlan);
            handle->inPlan = NULL;
        }

        if (handle->outPlan) {
            free(handle->outPlan);
            handle->outPlan = NULL;
        }

        if (handle->samples) {
            free(handle->samples);
            handle->samples = NULL;
        }

        if (handle->synthesisMem) {
            free(handle->synthesisMem);
            handle->synthesisMem = NULL;
        }

        if (handle->inWindowing) {
            free(handle->inWindowing);
            handle->inWindowing = NULL;
        }

        if (handle->outWindowing) {
            free(handle->outWindowing);
            handle->outWindowing = NULL;
        }
    }
}

void FFT_Resampler_reset(FFT_Resampler_Handle *handle)
{
    if (handle)
        handle->pos = 0;
}

int FFT_Resampler_Init(FFT_Resampler_Handle *handle, size_t inSampleRate, size_t outSampleRate, size_t nFFT)
{
    if (handle) {
        handle->pos = 0;
        if (outSampleRate < inSampleRate) {
            nFFT = nFFT * inSampleRate * 128 / outSampleRate;
        }
        else {
            nFFT = nFFT * outSampleRate * 128 / inSampleRate;
        }
        nFFT += (nFFT % 2);
        handle->inFrameSize = nFFT;
        handle->inWindowSize = nFFT * 2;
        handle->inSampleRate = inSampleRate;
        handle->inWindowing = (float *) calloc(handle->inFrameSize, sizeof(float));
        handle->inPlan = stb_fft_real_plan_dft_1d(handle->inWindowSize);
        handle->outSampleRate = outSampleRate;
        handle->outFrameSize = handle->inFrameSize * outSampleRate / inSampleRate;
        handle->outWindowSize = handle->outFrameSize * 2;
        handle->outWindowing = (float *) calloc(handle->outFrameSize, sizeof(float));
        handle->outPlan = stb_fft_real_plan_dft_1d(handle->outWindowSize);
        handle->inFifo = (float *) calloc(max(handle->inWindowSize, handle->outWindowSize), sizeof(float));
        handle->samples = (cmplx *) calloc(max(handle->inWindowSize, handle->outWindowSize), sizeof(cmplx));
        handle->synthesisMem = (float *) calloc(handle->outFrameSize, sizeof(float));
        if ((handle->inFifo == NULL) || (handle->inPlan == NULL) || (handle->outPlan == NULL)
            || (handle->samples == NULL)
            || (handle->synthesisMem == NULL) || (handle->inWindowing == NULL) || (handle->outWindowing == NULL)
            ) {
            FFT_Resampler_Free(handle);
            return 0;
        }
        float norm = 1.0f / handle->inWindowSize;
        for (size_t i = 0; i < handle->inFrameSize; i++) {
            double t = sin(.5 * M_PI * (i + .5) / handle->inFrameSize);
            handle->inWindowing[i] = (float) sin(.5 * M_PI * t * t) * norm;
        }
        for (size_t i = 0; i < handle->outFrameSize; i++) {
            double t = sin(.5 * M_PI * (i + .5) / handle->outFrameSize);
            handle->outWindowing[i] = (float) sin(.5 * M_PI * t * t);
        }
        return 1;
    }
    return 0;
}

int FFT_Resampler_Proc(FFT_Resampler_Handle *handle, const float *input, float *output)
{
    if ((input == NULL) || (handle == NULL) || (output == NULL)) {
        return -1;
    }
    float *inFifo = handle->inFifo;
    float *synthesis_mem = handle->synthesisMem;
    for (size_t i = 0; i < handle->inFrameSize; i++) {
        inFifo[i] *= handle->inWindowing[i];
        inFifo[handle->inWindowSize - 1 - i] = input[handle->inFrameSize - 1 - i] * handle->inWindowing[i];
    }
    stb_fft_r2c_exec(handle->inPlan, inFifo, handle->samples);
    if (handle->outWindowSize < handle->inWindowSize) {
        int half_output = (handle->outWindowSize / 2);
        int diff_size = handle->inWindowSize - handle->outWindowSize;
        memset(handle->samples + half_output, 0, diff_size * sizeof(cmplx));
    }
    else if (handle->outWindowSize > handle->inWindowSize) {
        int half_input = handle->inWindowSize / 2;
        int diff_size = handle->outWindowSize - handle->inWindowSize;
        memmove(handle->samples + half_input + diff_size, handle->samples + half_input,
                half_input * sizeof(cmplx));
        memset(handle->samples + half_input, 0, diff_size * sizeof(cmplx));
    }
    stb_fft_c2r_exec(handle->outPlan, handle->samples, inFifo);
    for (size_t i = 0; i < handle->outFrameSize; i++) {
        output[i] = inFifo[i] * handle->outWindowing[i] + synthesis_mem[i];
        inFifo[handle->outWindowSize - 1 - i] *= handle->outWindowing[i];
    }
    memcpy(synthesis_mem, inFifo + handle->outFrameSize, handle->outFrameSize * sizeof(float));
    memcpy(inFifo, input, handle->inFrameSize * sizeof(float));
    if (handle->pos == 0) {
        handle->pos++;
        return 0;
    }
    handle->pos++;
    return 1;
}

void printUsage()
{
    printf("usage:\n");
    printf("./FFTResampler input.wav 48000\n");
    printf("./FFTResampler input.mp3 16000\n");
    printf("or\n");
    printf("./FFTResampler input.wav output.wav 8000\n");
    printf("./FFTResampler input.mp3 output.wav 44100\n");
    printf("press any key to exit.\n");
    getchar();
}

void FFT_Resampler(char *in_file, char *out_file, uint32_t targetSampleRate)
{
    if (targetSampleRate == 0) {
        printUsage();
        return;
    }
    uint32_t sampleRate = 0;
    uint64_t sampleCount = 0;
    uint32_t channels = 0;
    float *input = wavRead_f32(in_file, &sampleRate, &sampleCount, &channels);
    if (input) {
        size_t nFFT = 2;
        FFT_Resampler_Handle *handle = (FFT_Resampler_Handle *) malloc(sizeof(FFT_Resampler_Handle));
        if (handle) {
            double startTime = now();
            if (FFT_Resampler_Init(handle, sampleRate, targetSampleRate, nFFT) == 1) {
                uint64_t frames = (sampleCount / handle->inFrameSize);
                int remainingSample = (sampleCount % handle->inFrameSize);
                int outRemainingSize = remainingSample * targetSampleRate / sampleRate;
                uint64_t targetSampleCount = frames * handle->outFrameSize + outRemainingSize;
                float *output = (float *) calloc(targetSampleCount, sizeof(float));
                if (output) {
                    FFT_Resampler_reset(handle);
                    float *inBuffer = input;
                    float *outBuffer = output;
                    for (int n = 0; n < frames; ++n) {
                        if (FFT_Resampler_Proc(handle, inBuffer, outBuffer) == 1)
                            outBuffer += handle->outFrameSize;
                        inBuffer += handle->inFrameSize;
                    }
                    if (remainingSample == 0) {
                        inBuffer -= handle->inFrameSize;
                        FFT_Resampler_Proc(handle, inBuffer, outBuffer);
                        memcpy(outBuffer, handle->synthesisMem, sizeof(float) * handle->outFrameSize);
                    }
                    else {
                        float *buffer = (float *) calloc(handle->inFrameSize, sizeof(float));
                        if (buffer) {
                            memcpy(buffer, inBuffer, sizeof(float) * remainingSample);
                            if (FFT_Resampler_Proc(handle, buffer, outBuffer) == 1)
                                outBuffer += handle->outFrameSize;
                            free(buffer);

                            memcpy(outBuffer, handle->synthesisMem, sizeof(float) * outRemainingSize);
                        }
                    }
                }
                double time_interval = calcElapsed(startTime, now());
                printf("time interval: %f ms\n ", (time_interval * 1000));
                wavWrite_f32(out_file, output, targetSampleRate, (uint32_t) targetSampleCount, channels);
                free(output);
            }
            FFT_Resampler_Free(handle);
            free(handle);
        }
        free(input);
    }
}

int main(int argc, char *argv[])
{
    printf("Audio Processing\n");
    printf("blog:http://cpuimage.cnblogs.com/\n");
    printf("Audio FFT Resampler\n");
    if (argc < 3) {
        printUsage();
        return -1;
    }
    char *in_file = argv[1];
    if (argc > 3) {
        char *out_file = argv[2];
        uint32_t targetSampleRate = (uint32_t) atoi(argv[3]);
        FFT_Resampler(in_file, out_file, targetSampleRate);
    }
    else {
        int32_t targetSampleRate = (uint32_t) atoi(argv[2]);
        char drive[3];
        char dir[256];
        char fname[256];
        char ext[256];
        char out_file[1024];
        splitpath(in_file, drive, dir, fname, ext);
        sprintf(out_file, "%s%s%s_out.wav", drive, dir, fname);
        FFT_Resampler(in_file, out_file, targetSampleRate);
    }

    return 0;
}

#ifdef __cplusplus
}
#endif
