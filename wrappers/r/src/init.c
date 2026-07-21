#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP C_bind_rust(SEXP path);
SEXP C_unbind_rust(SEXP unused);
SEXP C_parse_mzml(SEXP data);
SEXP C_ion_to_json(SEXP bin);
SEXP C_ion_to_mzml(SEXP bin);
SEXP C_get_peak(SEXP x, SEXP y, SEXP rt, SEXP range, SEXP options);
SEXP C_fit_peak(SEXP x, SEXP y, SEXP rt, SEXP intensity, SEXP shape);
SEXP C_draw_peak(SEXP x, SEXP shape, SEXP height, SEXP center, SEXP fwhm, SEXP tail);
SEXP C_get_peaks_from_eic(SEXP bin, SEXP rts, SEXP mzs, SEXP ranges, SEXP ids, SEXP from_left, SEXP to_right, SEXP options, SEXP cores);
SEXP C_get_peaks_from_chrom(SEXP bin, SEXP idxs, SEXP rts, SEXP ranges, SEXP options, SEXP cores);
SEXP C_calculate_eic(SEXP bin, SEXP targets, SEXP from, SEXP to, SEXP ppm_tol, SEXP mz_tol);
SEXP C_find_peaks(SEXP x, SEXP y, SEXP options);
SEXP C_calculate_baseline(SEXP y, SEXP lambda, SEXP max_iterations);
SEXP C_find_features(SEXP data, SEXP from_time, SEXP to_time, SEXP eic_ppm_tol, SEXP eic_mz_tol, SEXP grid_start, SEXP grid_end, SEXP grid_step_ppm, SEXP options, SEXP cores);
SEXP C_find_feature(SEXP bin, SEXP rts, SEXP mzs, SEXP wins, SEXP ids, SEXP scan_ppm, SEXP scan_mz, SEXP eic_ppm, SEXP eic_mz, SEXP options, SEXP cores);
SEXP C_mzml_to_ion(SEXP bin, SEXP level, SEXP f32_compress);
SEXP C_mzml_to_ion_file(SEXP input_path, SEXP output_path, SEXP level, SEXP f32_compress);
SEXP C_parse_ion(SEXP bin, SEXP max_cache_size);
SEXP C_plan_open(SEXP header);
SEXP C_plan_eic(SEXP ptr, SEXP target, SEXP from, SEXP to, SEXP ppm, SEXP mz_tol);
SEXP C_store_new(void);
SEXP C_store_add(SEXP store_ptr, SEXP offset, SEXP bytes);
SEXP C_parse_ion_source(SEXP store_ptr, SEXP max_cache_size);
SEXP C_parse_ion_path(SEXP path, SEXP max_cache_size);
SEXP C_get_scans(SEXP bin, SEXP query_type, SEXP a, SEXP b, SEXP level);
SEXP C_find_noise_level(SEXP y);
SEXP C_get_features(SEXP dir_path, SEXP from_time, SEXP to_time, SEXP eic_ppm_tol, SEXP eic_mz_tol,
                    SEXP grid_start, SEXP grid_end, SEXP grid_step,
                    SEXP group_ppm_tol, SEXP group_mz_tol, SEXP group_rt_tol,
                    SEXP frequency, SEXP options, SEXP cores);
SEXP C_dispose_mzml(SEXP ptr);

static const R_CallMethodDef CallEntries[] = {
    {"C_bind_rust", (DL_FUNC)&C_bind_rust, 1},
    {"C_unbind_rust", (DL_FUNC)&C_unbind_rust, 1},
    {"C_parse_mzml", (DL_FUNC)&C_parse_mzml, 1},
    {"C_ion_to_json", (DL_FUNC)&C_ion_to_json, 1},
    {"C_ion_to_mzml", (DL_FUNC)&C_ion_to_mzml, 1},
    {"C_get_peak", (DL_FUNC)&C_get_peak, 5},
    {"C_fit_peak", (DL_FUNC)&C_fit_peak, 5},
    {"C_draw_peak", (DL_FUNC)&C_draw_peak, 6},
    {"C_get_peaks_from_eic", (DL_FUNC)&C_get_peaks_from_eic, 9},
    {"C_get_peaks_from_chrom", (DL_FUNC)&C_get_peaks_from_chrom, 6},
    {"C_calculate_eic", (DL_FUNC)&C_calculate_eic, 6},
    {"C_find_peaks", (DL_FUNC)&C_find_peaks, 3},
    {"C_calculate_baseline", (DL_FUNC)&C_calculate_baseline, 3},
    {"C_find_features", (DL_FUNC)&C_find_features, 10},
    {"C_find_feature", (DL_FUNC)&C_find_feature, 11},
    {"C_mzml_to_ion", (DL_FUNC)&C_mzml_to_ion, 3},
    {"C_mzml_to_ion_file", (DL_FUNC)&C_mzml_to_ion_file, 4},
    {"C_parse_ion", (DL_FUNC)&C_parse_ion, 2},
    {"C_plan_open", (DL_FUNC)&C_plan_open, 1},
    {"C_plan_eic", (DL_FUNC)&C_plan_eic, 6},
    {"C_store_new", (DL_FUNC)&C_store_new, 0},
    {"C_store_add", (DL_FUNC)&C_store_add, 3},
    {"C_parse_ion_source", (DL_FUNC)&C_parse_ion_source, 2},
    {"C_parse_ion_path", (DL_FUNC)&C_parse_ion_path, 2},
    {"C_get_scans", (DL_FUNC)&C_get_scans, 5},
    {"C_find_noise_level", (DL_FUNC)&C_find_noise_level, 1},
    {"C_get_features", (DL_FUNC)&C_get_features, 14},
    {"C_dispose_mzml", (DL_FUNC)&C_dispose_mzml, 1},
    {NULL, NULL, 0}};

void R_init_quantion(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
