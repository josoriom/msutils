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
  int shape;
  double noise;
  int auto_noise;
  int auto_baseline;
  int lambda;
  int max_iterations;
  int allow_overlap;
  double min_snr;
  double min_r2;
  int kernel_size;
} CPeakOptions;

typedef struct Buf {
  uint8_t *ptr;
  uintptr_t len;
} Buf;

extern void js_log(const uint8_t *ptr, uintptr_t len);

extern int32_t range_read(uint32_t source_id,
                          uint32_t offset_lo,
                          uint32_t offset_hi,
                          uint32_t len,
                          uint8_t *dest_ptr);

extern int32_t parse_mzml(const uint8_t *data, uintptr_t len, void **out);
extern int32_t parse_bin(const uint8_t *data, uintptr_t len, uintptr_t cache_bytes, void **out);
extern int32_t parse_ion_path(const char *path, uintptr_t cache_bytes, void **out);
extern int32_t parse_ion_url(const char *url, uintptr_t cache_bytes, void **out);
extern int32_t plan_open(const uint8_t *header_ptr, uintptr_t header_len, struct Buf *out);
extern int32_t plan_eic(void *h, double target, double from, double to, double ppm, double mz_tol, struct Buf *out);
extern int32_t calculate_eic(void *h, double target, double from, double to, double ppm, double mz_tol, struct Buf *out_x, struct Buf *out_y);
extern int32_t bin_to_json(void *h, struct Buf *out);
extern int32_t bin_to_mzml(void *h, struct Buf *out);
extern int32_t mzml_to_bin(void *h, struct Buf *out, uint8_t level, uint8_t compress);
extern void dispose_mzml(void *h);
extern void free_(uint8_t *ptr, uintptr_t len);
