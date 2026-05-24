#pragma once

/* Generated with cbindgen:0.27.0 */

#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#define MSUTILS_ABI_VERSION 1

typedef struct CPeakOptions {
  double min_integral;
  double min_intensity;
  int min_peak_width_points;
  double noise;
  int auto_noise;
  int auto_baseline;
  int lambda;
  int max_iterations;
  int allow_overlap;
  double min_snr;
} CPeakOptions;

typedef struct Buf {
  uint8_t *ptr;
  uintptr_t len;
} Buf;

extern void js_log(const uint8_t *ptr, uintptr_t len);
