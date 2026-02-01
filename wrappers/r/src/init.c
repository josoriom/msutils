#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP C_bind_rust(SEXP path);
SEXP C_parse_mzml(SEXP data);
SEXP C_parse_mzmlb(SEXP path);
SEXP C_bin_to_json(SEXP bin);
SEXP C_bin_to_mzml(SEXP bin);
SEXP C_get_peak(SEXP x, SEXP y, SEXP rt, SEXP range, SEXP options);
SEXP C_get_peaks_from_eic(SEXP bin, SEXP rts, SEXP mzs, SEXP ranges, SEXP ids, SEXP from_left, SEXP to_right, SEXP options, SEXP cores);
SEXP C_get_peaks_from_chrom(SEXP bin, SEXP idxs, SEXP rts, SEXP ranges, SEXP options, SEXP cores);
SEXP C_calculate_eic(SEXP bin, SEXP targets, SEXP from, SEXP to, SEXP ppm_tol, SEXP mz_tol);
SEXP C_find_peaks(SEXP x, SEXP y, SEXP options);
SEXP C_calculate_baseline(SEXP y, SEXP baseline_window, SEXP baseline_window_factor);
SEXP C_find_features(SEXP data, SEXP from_time, SEXP to_time, SEXP eic_ppm_tol, SEXP eic_mz_tol, SEXP grid_start, SEXP grid_end, SEXP grid_step_ppm, SEXP options);
SEXP C_find_feature(SEXP bin, SEXP rts, SEXP mzs, SEXP wins, SEXP ids, SEXP scan_ppm, SEXP scan_mz, SEXP eic_ppm, SEXP eic_mz, SEXP options, SEXP cores);
SEXP C_convert_mzml_to_bin(SEXP xml, SEXP level);
SEXP C_parse_bin(SEXP bin);

static const R_CallMethodDef CallEntries[] = {
    {"C_bind_rust", (DL_FUNC)&C_bind_rust, 1},
    {"C_parse_mzml", (DL_FUNC)&C_parse_mzml, 1},
    {"C_parse_mzmlb", (DL_FUNC)&C_parse_mzmlb, 1},
    {"C_bin_to_json", (DL_FUNC)&C_bin_to_json, 1},
    {"C_bin_to_mzml", (DL_FUNC)&C_bin_to_mzml, 1},
    {"C_get_peak", (DL_FUNC)&C_get_peak, 5},
    {"C_get_peaks_from_eic", (DL_FUNC)&C_get_peaks_from_eic, 9},
    {"C_get_peaks_from_chrom", (DL_FUNC)&C_get_peaks_from_chrom, 6},
    {"C_calculate_eic", (DL_FUNC)&C_calculate_eic, 6},
    {"C_find_peaks", (DL_FUNC)&C_find_peaks, 3},
    {"C_calculate_baseline", (DL_FUNC)&C_calculate_baseline, 3},
    {"C_find_features", (DL_FUNC)&C_find_features, 9},
    {"C_find_feature", (DL_FUNC)&C_find_feature, 11},
    {"C_convert_mzml_to_bin", (DL_FUNC)&C_convert_mzml_to_bin, 2},
    {"C_parse_bin", (DL_FUNC)&C_parse_bin, 1},
    {NULL, NULL, 0}};

void R_init_msutils(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
