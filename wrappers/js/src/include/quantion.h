#pragma once

/* Generated with cbindgen:0.29.4 */

#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#define QUANTION_ABI_VERSION 2

typedef struct ImageSession ImageSession;

typedef struct ParsedFile ParsedFile;

typedef struct Buf {
  uint8_t *ptr;
  uintptr_t len;
} Buf;

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

uint32_t quantion_abi_version(void);

uintptr_t quantion_sizeof_peak_options(void);

#if !(defined(QUANTION_TARGET_WASM32) && !defined(QUANTION_TARGET_WASI))
/**
 * Open an ion file whose bytes the caller supplies on demand.
 *
 * The library never reads a file or a socket itself: every time it needs a slice
 * it calls `read`. Use this for a file that lives somewhere the library cannot
 * reach, such as a URL.
 *
 * # Safety
 * `read` must stay callable, and `context` must stay valid, until the returned
 * handle is freed with `free_mzml`.
 */
int parse_ion_source(int32_t (*read)(void *context, uint64_t offset, uint64_t length, uint8_t *dest),
                     void *context,
                     uintptr_t max_cache_size,
                     struct ParsedFile **out);
#endif

#if (defined(QUANTION_TARGET_WASM32) && !defined(QUANTION_TARGET_WASI))
extern int32_t range_read(uint32_t source_id,
                          uint32_t offset_lo,
                          uint32_t offset_hi,
                          uint32_t len,
                          uint8_t *dest_ptr);
#endif

uint8_t *alloc(uintptr_t size);

/**
 * Free memory returned by `alloc`.
 *
 * # Safety
 * `ptr_raw` must be a live pointer from `alloc(size)`.
 * `size` must match the `alloc` call.
 * After this call, `ptr_raw` is invalid.
 */
void free_(uint8_t *ptr_raw, uintptr_t size);

/**
 * Free a `ParsedFile` handle.
 *
 * # Safety
 * `handle` must be a valid pointer returned by this library.
 * Do not free the same handle twice.
 * After this call, `handle` is invalid.
 */
void free_mzml(struct ParsedFile *handle);

/**
 * Parse mzML bytes and store the result in `dest`.
 *
 * # Safety
 * `data_ptr` must point to `data_len` readable bytes.
 * `dest` must be a valid writable pointer to `*mut ParsedFile`.
 * On success, `*dest` must be freed with `free_mzml`.
 */
int parse_mzml(const uint8_t *data_ptr, uintptr_t data_len, struct ParsedFile **dest);

#if (defined(QUANTION_TARGET_WASM32) && !defined(QUANTION_TARGET_WASI))
int parse_ion_url(uint32_t source_id, uintptr_t cache_bytes, struct ParsedFile **dest);
#endif

/**
 * Parse binary data and store the result in `dest`.
 *
 * # Safety
 * `data_ptr` must point to `data_len` readable bytes.
 * `dest` must be a valid writable pointer to `*mut ParsedFile`.
 * On success, `*dest` must be freed with `free_mzml`.
 */
int parse_bin(const uint8_t *data_ptr,
              uintptr_t data_len,
              uintptr_t max_cache_size,
              struct ParsedFile **dest);

#if !(defined(QUANTION_TARGET_WASM32) && !defined(QUANTION_TARGET_WASI))
int parse_ion_path(const char *path_ptr, uintptr_t max_cache_size, struct ParsedFile **dest);
#endif

/**
 * Plan the byte ranges an open needs.
 *
 * # Safety
 * `header_ptr` must point to at least `header_len` readable bytes, and `header_len` must be at
 * least 408 so the total file size can be read.
 */
int plan_open(const uint8_t *header_ptr, uintptr_t header_len, struct Buf *out);

/**
 * Serialize a parsed file to JSON and write it to `out`.
 *
 * # Safety
 * `h` must be a valid `ParsedFile` pointer from this library.
 * `out` must be a valid writable `Buf` pointer.
 */
int bin_to_json(struct ParsedFile *h, struct Buf *out);

/**
 * Convert a parsed file to mzML and write the result to `out`.
 *
 * # Safety
 * `h` must be a valid `ParsedFile` pointer from this library.
 * `out` must be a valid writable `Buf` pointer.
 */
int bin_to_mzml(struct ParsedFile *h, struct Buf *out);

/**
 * Convert a parsed file to binary and write the result to `out`.
 *
 * # Safety
 * `h` must be a valid `ParsedFile` pointer from this library.
 * `out` must be a valid writable `Buf` pointer.
 */
int mzml_to_bin(struct ParsedFile *h,
                struct Buf *out,
                uint8_t compression_level,
                uint8_t f32_compress);

#if !(defined(QUANTION_TARGET_WASM32) && !defined(QUANTION_TARGET_WASI))
int convert_mzml_file_to_ion_file(const char *input_path,
                                  const char *output_path,
                                  uint8_t compression_level,
                                  uint8_t force_f32,
                                  uint8_t section_on_disk);
#endif

/**
 * Find one peak and write the result to `out`.
 *
 * # Safety
 * `x_ptr` and `y_ptr` must point to `len` readable `f64` values.
 * `options` must be null or point to a valid `CPeakOptions`.
 * `out` must be a valid writable `Buf` pointer.
 */
int get_peak(const double *x_ptr,
             const double *y_ptr,
             uintptr_t len,
             double rt,
             double range,
             const struct CPeakOptions *options,
             struct Buf *out);

/**
 * Fit a Gaussian or EMG model to a peak and write the parameters JSON to `out`.
 *
 * `shape` is 0 for Gaussian, 1 for EMG.
 *
 * # Safety
 * `x_ptr` and `y_ptr` must point to `len` readable `f64` values.
 * `out` must be a valid writable `Buf` pointer.
 */
int fit_peak(const double *x_ptr,
             const double *y_ptr,
             uintptr_t len,
             double rt,
             double intensity,
             int shape,
             struct Buf *out);

/**
 * Render a fitted peak model over `x` and write the y curve to `out_y`.
 *
 * `shape` is 0 for Gaussian, 1 for EMG. The output has the same length as `x`.
 *
 * # Safety
 * `x_ptr` must point to `len` readable `f64` values.
 * `out_y` must be a valid writable `Buf` pointer.
 */
int draw_peak(const double *x_ptr,
              uintptr_t len,
              int shape,
              double height,
              double center,
              double fwhm,
              double tail,
              struct Buf *out_y);

/**
 * Find peaks from EIC data and write the result to `out`.
 *
 * # Safety
 * `h` must be a valid `ParsedFile` pointer from this library.
 * `rts_ptr`, `mzs_ptr`, and `ranges_ptr` must point to `n` readable `f64` values.
 * If IDs are provided, `ids_off` and `ids_len` must point to `n` readable `u32` values,
 * and `ids_buf` must point to `ids_buf_len` readable bytes.
 * `opts` must be null or point to a valid `CPeakOptions`.
 * `out` must be a valid writable `Buf` pointer.
 */
int32_t get_peaks_from_eic(const struct ParsedFile *h,
                           const double *rts_ptr,
                           const double *mzs_ptr,
                           const double *ranges_ptr,
                           const uint32_t *ids_off,
                           const uint32_t *ids_len,
                           const uint8_t *ids_buf,
                           uintptr_t ids_buf_len,
                           uintptr_t n,
                           double from,
                           double to,
                           const struct CPeakOptions *opts,
                           uintptr_t cores,
                           struct Buf *out);

/**
 * Find peaks from chromatogram data and write the result to `out`.
 *
 * # Safety
 * `h` must be a valid `ParsedFile` pointer from this library.
 * `sample_indices` must point to `n` readable `u32` values.
 * `target_rts` and `half_widths` must point to `n` readable `f64` values.
 * `opts` must be null or point to a valid `CPeakOptions`.
 * `out` must be a valid writable `Buf` pointer.
 */
int32_t get_peaks_from_chrom(struct ParsedFile *h,
                             const uint32_t *sample_indices,
                             const double *target_rts,
                             const double *ranges,
                             uintptr_t n,
                             const struct CPeakOptions *opts,
                             uintptr_t cores,
                             struct Buf *out);

/**
 * Find peaks and write the result to `out`.
 *
 * # Safety
 * `x_ptr` and `y_ptr` must point to `len` readable `f64` values.
 * `opts` must be null or point to a valid `CPeakOptions`.
 * `out` must be a valid writable `Buf` pointer.
 */
int find_peaks(const double *x_ptr,
               const double *y_ptr,
               uintptr_t len,
               const struct CPeakOptions *opts,
               struct Buf *out);

/**
 * Write the noise window width and the noise intensity for an intensity array.
 *
 * # Safety
 * `y_ptr` must point to `len` readable `f32` values.
 * `out_width` and `out_intensity` must be valid writable pointers.
 */
int find_noise_level(const float *y_ptr,
                     uintptr_t len,
                     uintptr_t *out_width,
                     double *out_intensity);

/**
 * Calculate an EIC and write `x` and `y` to `out_x` and `out_y`.
 *
 * # Safety
 * `h` must be a valid `ParsedFile` pointer from this library.
 * `out_x` and `out_y` must be valid writable `Buf` pointers.
 */
int calculate_eic(struct ParsedFile *h,
                  double target,
                  double from,
                  double to,
                  double ppm,
                  double mz_tol,
                  struct Buf *out_x,
                  struct Buf *out_y);

int plan_eic(struct ParsedFile *h,
             double target,
             double from,
             double to,
             double ppm,
             double mz_tol,
             struct Buf *out);

/**
 * Get scans by query and write the result to `out`.
 *
 * `query_type`: 0=RtRange, 1=ClosestRt, 2=MzRange, 3=ClosestMz.
 * `from_value` and `to_value` are the range bounds for range queries (types 0 and 2).
 * For point queries (types 1 and 3), only `from_value` is used; `to_value` is ignored.
 *
 * # Safety
 * `h` must be a valid `ParsedFile` pointer from this library.
 * `out` must be a valid writable `Buf` pointer.
 */
int get_scans(struct ParsedFile *h,
              uint8_t query_type,
              double from_value,
              double to_value,
              uint8_t level,
              struct Buf *out);

/**
 * Build a 2D ion image for a target m/z by summing intensity in `[target - tolerance, target + tolerance]`
 * per spectrum and scattering the mean into a position_x/position_y grid. Writes a JSON object to `out`.
 *
 * # Safety
 * `h` must be a valid `ParsedFile` pointer from this library.
 * `out` must be a valid writable `Buf` pointer.
 */
int get_ion_image(struct ParsedFile *h,
                  double target,
                  double tolerance,
                  uint8_t level,
                  struct Buf *out);

int image_begin(struct ParsedFile *h,
                double target,
                double tolerance,
                uint8_t level,
                uintptr_t *out_session);

int image_scan_count(struct ImageSession *session, uintptr_t *out_count);

int image_ranges(struct ParsedFile *h,
                 struct ImageSession *session,
                 uintptr_t from,
                 uintptr_t count,
                 struct Buf *out);

int image_fold(struct ParsedFile *h, struct ImageSession *session, uintptr_t from, uintptr_t count);

int image_finish(struct ImageSession *session, struct Buf *out);

void image_free(struct ImageSession *session);

#if !(defined(QUANTION_TARGET_WASM32) && !defined(QUANTION_TARGET_WASI))
/**
 * Find features across all samples in `dir` and write the result to `out`.
 *
 * # Safety
 * `dir` must point to a valid NUL-terminated C string.
 * `peak_opts` must be null or point to a valid `CPeakOptions`.
 * `out` must be a valid writable `Buf` pointer.
 */
int get_features(const char *dir,
                 double from,
                 double to,
                 double eic_ppm,
                 double eic_mz,
                 double grid_min_mz,
                 double grid_max_mz,
                 double grid_step,
                 double align_ppm,
                 double align_mz_absolute,
                 double align_rt_tolerance,
                 int min_samples,
                 const struct CPeakOptions *peak_opts,
                 int cores,
                 struct Buf *out);
#endif

/**
 * Find features in a single sample and write the result to `out`.
 *
 * # Safety
 * `h` must be a valid `ParsedFile` pointer from this library.
 * `peak_opts` must be null or point to a valid `CPeakOptions`.
 * `out` must be a valid writable `Buf` pointer.
 */
int find_features(struct ParsedFile *h,
                  double from,
                  double to,
                  double eic_ppm,
                  double eic_mz,
                  double grid_min_mz,
                  double grid_max_mz,
                  double grid_step,
                  const struct CPeakOptions *peak_opts,
                  int cores,
                  struct Buf *out);

/**
 * Calculate a baseline and write the result to `out`.
 *
 * # Safety
 * `y` must point to `len` readable `f64` values.
 * `out` must be a valid writable `Buf` pointer.
 */
int calculate_baseline(const double *y,
                       uintptr_t len,
                       int lambda,
                       int max_iterations,
                       struct Buf *out);

/**
 * Find a targeted feature for each ROI and write the result to `out`.
 *
 * # Safety
 * `h` must be a valid `ParsedFile` pointer from this library.
 * `target_rts`, `target_mzs`, and `half_widths` must point to `n` readable `f64` values.
 * If IDs are provided, `ids_off` and `ids_len` must point to `n` readable `u32` values,
 * and `ids_buf` must point to `ids_buf_len` readable bytes.
 * `peak_opts` must be null or point to a valid `CPeakOptions`.
 * `out` must be a valid writable `Buf` pointer.
 */
int find_feature(const struct ParsedFile *h,
                 const double *target_rts,
                 const double *target_mzs,
                 const double *half_widths,
                 const uint32_t *ids_off,
                 const uint32_t *ids_len,
                 const uint8_t *ids_buf,
                 uintptr_t ids_buf_len,
                 uintptr_t n,
                 uintptr_t cores,
                 double seed_ppm,
                 double seed_mz,
                 double final_ppm,
                 double final_mz,
                 const struct CPeakOptions *peak_opts,
                 struct Buf *out);
