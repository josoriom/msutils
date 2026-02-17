#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#if defined(_WIN32)
#include <windows.h>
#define DLIB HMODULE
#define DLOPEN(p) LoadLibraryA(p)
#define DLSYM(h, s) GetProcAddress(h, s)
#define DLCLOSE(h) FreeLibrary(h)
static const char *last_err = "LoadLibrary/GetProcAddress failed";
#else
#include <dlfcn.h>
#define DLIB void *
#define DLOPEN(p) dlopen(p, RTLD_NOW | RTLD_GLOBAL)
#define DLSYM(h, s) dlsym(h, s)
#define DLCLOSE(h) dlclose(h)
static const char *last_err = NULL;
#endif

typedef struct MzML MzML;

typedef struct
{
  unsigned char *ptr;
  size_t len;
} Buf;

typedef struct
{
  double integral_threshold;
  double intensity_threshold;
  int32_t width_threshold;
  double noise;
  int32_t auto_noise;
  int32_t auto_baseline;
  int32_t baseline_window;
  int32_t baseline_window_factor;
  int32_t allow_overlap;
  int32_t window_size;
  double sn_ratio;
} CPeakPOptions;

typedef int32_t (*fn_parse_mzml)(const unsigned char *, size_t, MzML **);
typedef int32_t (*fn_bin_to_json)(const unsigned char *, size_t, Buf *);
typedef int32_t (*fn_bin_to_mzml)(const unsigned char *, size_t, Buf *);
typedef int32_t (*fn_get_peak)(const double *, const double *, size_t, double, double, const CPeakPOptions *, Buf *);
typedef int32_t (*fn_calculate_eic)(const unsigned char *, size_t, double, double, double, double, double, Buf *, Buf *);
typedef float (*fn_find_noise_level)(const float *, size_t);
typedef int32_t (*fn_get_peaks_from_eic)(const MzML *, const double *, const double *, const double *, const uint32_t *, const uint32_t *, const unsigned char *, size_t, size_t, double, double, const CPeakPOptions *, size_t, Buf *);
typedef int32_t (*fn_get_peaks_from_chrom)(const MzML *, const uint32_t *, const double *, const double *, size_t, const CPeakPOptions *, size_t, Buf *);
typedef int32_t (*fn_find_peaks)(const double *, const double *, size_t, const CPeakPOptions *, Buf *);
typedef int32_t (*fn_calculate_baseline)(const double *, size_t, int32_t, int32_t, Buf *);
typedef int32_t (*fn_find_features)(const MzML *, double, double, double, double, double, double, double, const CPeakPOptions *, int32_t, Buf *);
typedef int32_t (*fn_find_feature)(const MzML *, const double *, const double *, const double *, const uint32_t *, const uint32_t *, const unsigned char *, size_t, size_t, size_t, double, double, double, double, const CPeakPOptions *, Buf *);
typedef int32_t (*fn_convert_mzml_to_bin)(const unsigned char *, size_t, Buf *, uint8_t, uint8_t);
typedef int32_t (*fn_parse_bin)(const unsigned char *, size_t, MzML **);
typedef void (*fn_free_)(unsigned char *, size_t);
typedef void (*fn_free_mzml)(MzML *);

typedef struct
{
  fn_parse_mzml parse_mzml;
  fn_bin_to_json bin_to_json;
  fn_bin_to_mzml bin_to_mzml;
  fn_get_peak get_peak;
  fn_calculate_eic calculate_eic;
  fn_find_noise_level find_noise_level;
  fn_get_peaks_from_eic get_peaks_from_eic;
  fn_get_peaks_from_chrom get_peaks_from_chrom;
  fn_find_peaks find_peaks;
  fn_calculate_baseline calculate_baseline;
  fn_find_features find_features;
  fn_free_ free_;
  fn_find_feature find_feature;
  fn_convert_mzml_to_bin convert_mzml_to_bin;
  fn_parse_bin parse_bin;
  fn_free_mzml free_mzml;
} abi_type;

static DLIB abi_handle = NULL;
abi_type ABI = (abi_type){0};

static int resolve_required(void **fn, const char *name)
{
  *fn = DLSYM(abi_handle, name);
  return *fn ? 0 : -1;
}

int abi_load(const char *path, const char **err)
{
  if (abi_handle)
  {
    DLCLOSE(abi_handle);
    abi_handle = NULL;
  }
  memset(&ABI, 0, sizeof(ABI));
  abi_handle = DLOPEN(path);
#if !defined(_WIN32)
  if (!abi_handle)
    last_err = dlerror();
#endif
  if (!abi_handle)
  {
    if (err)
      *err = last_err;
    return -1;
  }
  if (resolve_required((void **)&ABI.parse_mzml, "parse_mzml"))
    goto fail;
  if (resolve_required((void **)&ABI.bin_to_json, "bin_to_json"))
    goto fail;
  if (resolve_required((void **)&ABI.bin_to_mzml, "bin_to_mzml"))
    goto fail;
  if (resolve_required((void **)&ABI.get_peak, "get_peak"))
    goto fail;
  if (resolve_required((void **)&ABI.calculate_eic, "calculate_eic"))
    goto fail;
  if (resolve_required((void **)&ABI.find_noise_level, "find_noise_level"))
    goto fail;
  if (resolve_required((void **)&ABI.get_peaks_from_eic, "get_peaks_from_eic"))
    goto fail;
  if (resolve_required((void **)&ABI.get_peaks_from_chrom, "get_peaks_from_chrom"))
    goto fail;
  if (resolve_required((void **)&ABI.find_peaks, "find_peaks"))
    goto fail;
  if (resolve_required((void **)&ABI.calculate_baseline, "calculate_baseline"))
    goto fail;
  if (resolve_required((void **)&ABI.find_features, "find_features"))
    goto fail;
  if (resolve_required((void **)&ABI.find_feature, "find_feature"))
    goto fail;
  if (resolve_required((void **)&ABI.convert_mzml_to_bin, "convert_mzml_to_bin"))
    goto fail;
  if (resolve_required((void **)&ABI.parse_bin, "parse_bin"))
    goto fail;
  if (resolve_required((void **)&ABI.free_mzml, "free_mzml"))
    goto fail;
  ABI.free_ = (fn_free_)DLSYM(abi_handle, "free_");
  if (!ABI.free_)
    goto fail;
  return 0;
fail:
  if (abi_handle)
  {
    DLCLOSE(abi_handle);
    abi_handle = NULL;
  }
  memset(&ABI, 0, sizeof(ABI));
  if (err)
    *err = last_err;
  return -1;
}

void abi_unload(void)
{
  if (abi_handle)
  {
    DLCLOSE(abi_handle);
    abi_handle = NULL;
  }
  memset(&ABI, 0, sizeof(ABI));
}

static void die_code(const char *fname, int code)
{
  const char *msg = "unknown error";
  if (code == 0)
    return;
  if (code == 1)
    msg = "invalid arguments";
  else if (code == 2)
    msg = "panic inside Rust";
  else if (code == 4)
    msg = "parse error";
  error("msutils/%s failed: %s (code=%d)", fname, msg, code);
}

static SEXP mk_string_len(const unsigned char *ptr, size_t len)
{
  SEXP s = PROTECT(Rf_ScalarString(Rf_mkCharLenCE((const char *)ptr, (int)len, CE_UTF8)));
  UNPROTECT(1);
  return s;
}

#define REQUIRE_BOUND(ptr, name)                                            \
  do                                                                        \
  {                                                                         \
    if ((ptr) == NULL)                                                      \
      error("msutils: symbol %s is not bound; did .onLoad() run?", (name)); \
  } while (0)

static SEXP list_get(SEXP lst, const char *name)
{
  if (TYPEOF(lst) != VECSXP)
    return R_NilValue;
  SEXP names = Rf_getAttrib(lst, R_NamesSymbol);
  if (TYPEOF(names) != STRSXP)
    return R_NilValue;
  R_xlen_t n = XLENGTH(lst);
  for (R_xlen_t i = 0; i < n; i++)
  {
    SEXP nm = STRING_ELT(names, i);
    if (nm == R_NilValue)
      continue;
    if (strcmp(CHAR(nm), name) == 0)
      return VECTOR_ELT(lst, i);
  }
  return R_NilValue;
}

static int fill_options(SEXP opts, CPeakPOptions *out)
{
  if (opts == R_NilValue || TYPEOF(opts) != VECSXP || XLENGTH(opts) == 0)
    return 0;
  out->integral_threshold = NAN;
  out->intensity_threshold = NAN;
  out->width_threshold = 0;
  out->noise = NAN;
  out->auto_noise = 0;
  out->auto_baseline = 0;
  out->baseline_window = 0;
  out->baseline_window_factor = 0;
  out->allow_overlap = 0;
  out->window_size = 0;
  out->sn_ratio = NAN;
  SEXP v = R_NilValue;
  v = list_get(opts, "integral_threshold");
  if (v != R_NilValue)
    out->integral_threshold = asReal(v);
  v = list_get(opts, "intensity_threshold");
  if (v != R_NilValue)
    out->intensity_threshold = asReal(v);
  v = list_get(opts, "width_threshold");
  if (v != R_NilValue)
    out->width_threshold = (int32_t)asInteger(v);
  v = list_get(opts, "noise");
  if (v != R_NilValue)
    out->noise = asReal(v);
  v = list_get(opts, "auto_noise");
  if (v != R_NilValue)
    out->auto_noise = (int32_t)asLogical(v);
  v = list_get(opts, "auto_baseline");
  if (v != R_NilValue)
    out->auto_baseline = (int32_t)asLogical(v);
  v = list_get(opts, "baseline_window");
  if (v != R_NilValue)
    out->baseline_window = (int32_t)asInteger(v);
  v = list_get(opts, "baseline_window_factor");
  if (v != R_NilValue)
    out->baseline_window_factor = (int32_t)asInteger(v);
  v = list_get(opts, "allow_overlap");
  if (v != R_NilValue)
    out->allow_overlap = (int32_t)asLogical(v);
  v = list_get(opts, "window_size");
  if (v != R_NilValue)
    out->window_size = (int32_t)asInteger(v);
  v = list_get(opts, "sn_ratio");
  if (v != R_NilValue)
    out->sn_ratio = asReal(v);
  return 1;
}

static int as_opts_ptr(SEXP options, CPeakPOptions *copy, const CPeakPOptions **out_ptr)
{
  *out_ptr = NULL;
  if (options == R_NilValue)
    return 0;
  if (TYPEOF(options) == RAWSXP)
  {
    if ((size_t)XLENGTH(options) != sizeof(CPeakPOptions))
      error("msutils: options raw blob must be length %zu", sizeof(CPeakPOptions));
    memcpy((void *)copy, (const void *)RAW(options), sizeof(CPeakPOptions));
    *out_ptr = copy;
    return 1;
  }
  if (TYPEOF(options) == VECSXP && XLENGTH(options) > 0)
  {
    if (fill_options(options, copy))
    {
      *out_ptr = copy;
      return 1;
    }
  }
  return 0;
}

static MzML *GetHandle(SEXP ptr)
{
  if (TYPEOF(ptr) != EXTPTRSXP)
    error("msutils: expected ExternalPtr");
  MzML *h = (MzML *)R_ExternalPtrAddr(ptr);
  if (!h)
    error("msutils: use of disposed or null pointer");
  return h;
}

SEXP C_bind_rust(SEXP path_)
{
  if (TYPEOF(path_) != STRSXP || LENGTH(path_) != 1)
    error("path");
  const char *path = CHAR(STRING_ELT(path_, 0));
  const char *err = NULL;
  if (abi_load(path, &err) != 0)
    error("dlopen failed: %s", err ? err : "unknown");
  return R_NilValue;
}

static void finalize_mzml(SEXP ptr)
{
  if (R_ExternalPtrAddr(ptr))
  {
    MzML *handle = (MzML *)R_ExternalPtrAddr(ptr);
    if (ABI.free_mzml)
      ABI.free_mzml(handle);
    R_ClearExternalPtr(ptr);
  }
}

SEXP C_dispose_mzml(SEXP ptr)
{
  finalize_mzml(ptr);
  return R_NilValue;
}

SEXP C_parse_mzml(SEXP data)
{
  if (TYPEOF(data) != RAWSXP)
    error("data must be raw");
  REQUIRE_BOUND(ABI.parse_mzml, "parse_mzml");

  MzML *handle = NULL;
  int code = ABI.parse_mzml((const unsigned char *)RAW(data), (size_t)XLENGTH(data), &handle);
  die_code("parse_mzml", code);

  SEXP ptr = PROTECT(R_MakeExternalPtr(handle, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(ptr, finalize_mzml, TRUE);
  UNPROTECT(1);
  return ptr;
}

SEXP C_bin_to_json(SEXP bin)
{
  if (TYPEOF(bin) != RAWSXP)
    error("bin");
  REQUIRE_BOUND(ABI.bin_to_json, "bin_to_json");
  REQUIRE_BOUND(ABI.free_, "free_");
  Buf out = (Buf){0};
  int code = ABI.bin_to_json((const unsigned char *)RAW(bin), (size_t)XLENGTH(bin), &out);
  die_code("bin_to_json", code);
  SEXP res = mk_string_len(out.ptr, out.len);
  ABI.free_(out.ptr, out.len);
  return res;
}

SEXP C_bin_to_mzml(SEXP bin)
{
  if (TYPEOF(bin) != RAWSXP)
    error("bin");
  REQUIRE_BOUND(ABI.bin_to_mzml, "bin_to_mzml");
  REQUIRE_BOUND(ABI.free_, "free_");
  Buf out = (Buf){0};
  int code = ABI.bin_to_mzml((const unsigned char *)RAW(bin), (size_t)XLENGTH(bin), &out);
  die_code("bin_to_mzml", code);
  SEXP res = mk_string_len(out.ptr, out.len);
  ABI.free_(out.ptr, out.len);
  return res;
}

SEXP C_get_peak(SEXP x, SEXP y, SEXP rt, SEXP range, SEXP options)
{
  if (TYPEOF(x) != REALSXP || TYPEOF(y) != REALSXP)
    error("numeric");
  if (XLENGTH(x) != XLENGTH(y) || XLENGTH(x) < 3)
    error("length");
  REQUIRE_BOUND(ABI.get_peak, "get_peak");
  REQUIRE_BOUND(ABI.free_, "free_");
  R_xlen_t n = XLENGTH(y);
  CPeakPOptions opts;
  const CPeakPOptions *opt_ptr = NULL;
  (void)as_opts_ptr(options, &opts, &opt_ptr);
  Buf out = (Buf){0};
  int code = ABI.get_peak(REAL(x), REAL(y), (size_t)n, asReal(rt), asReal(range), opt_ptr, &out);
  die_code("get_peak", code);
  SEXP res = mk_string_len(out.ptr, out.len);
  ABI.free_(out.ptr, out.len);
  return res;
}

SEXP C_get_peaks_from_eic(SEXP bin, SEXP rts, SEXP mzs, SEXP ranges, SEXP ids, SEXP from_left, SEXP to_right, SEXP options, SEXP cores)
{
  MzML *handle = GetHandle(bin);
  REQUIRE_BOUND(ABI.get_peaks_from_eic, "get_peaks_from_eic");
  REQUIRE_BOUND(ABI.free_, "free_");

  if (TYPEOF(rts) != REALSXP || TYPEOF(mzs) != REALSXP || TYPEOF(ranges) != REALSXP)
    error("msutils: bad numeric arguments");
  R_xlen_t n = XLENGTH(rts);
  if (XLENGTH(mzs) != n || XLENGTH(ranges) != n)
    error("msutils: length mismatch");
  uint32_t *offs = (uint32_t *)R_alloc((size_t)n, sizeof(uint32_t));
  uint32_t *lens = (uint32_t *)R_alloc((size_t)n, sizeof(uint32_t));
  unsigned char *ids_buf = NULL;
  size_t ids_len = 0;
  if (ids != R_NilValue)
  {
    if (TYPEOF(ids) != STRSXP)
      error("msutils: ids must be character");
    size_t total = 0;
    for (R_xlen_t i = 0; i < n; i++)
    {
      SEXP s = STRING_ELT(ids, i);
      if (s != R_NilValue)
        total += (size_t)LENGTH(s);
    }
    ids_buf = (unsigned char *)R_alloc(total, 1);
    ids_len = total;
    size_t cur = 0;
    for (R_xlen_t i = 0; i < n; i++)
    {
      SEXP s = STRING_ELT(ids, i);
      if (s == R_NilValue)
      {
        offs[i] = 0;
        lens[i] = 0;
      }
      else
      {
        size_t L = (size_t)LENGTH(s);
        offs[i] = (uint32_t)cur;
        lens[i] = (uint32_t)L;
        memcpy(ids_buf + cur, CHAR(s), L);
        cur += L;
      }
    }
  }

  size_t ncores = (cores == R_NilValue) ? 1 : (size_t)asInteger(cores);
  CPeakPOptions opts;
  const CPeakPOptions *opt_ptr = NULL;
  as_opts_ptr(options, &opts, &opt_ptr);

  Buf out = (Buf){0};
  int code = ABI.get_peaks_from_eic(
      handle, REAL(rts), REAL(mzs), REAL(ranges),
      offs, lens, ids_buf, ids_len,
      (size_t)n, asReal(from_left), asReal(to_right), opt_ptr, ncores, &out);

  die_code("get_peaks_from_eic", code);
  SEXP res = mk_string_len(out.ptr, out.len);
  ABI.free_(out.ptr, out.len);
  return res;
}

SEXP C_get_peaks_from_chrom(SEXP bin, SEXP idxs, SEXP rts, SEXP ranges, SEXP options, SEXP cores)
{
  MzML *handle = GetHandle(bin);

  if (TYPEOF(rts) != REALSXP || TYPEOF(ranges) != REALSXP)
    error("msutils: numeric (double) required for rts/ranges");

  R_xlen_t n = XLENGTH(rts);
  if (XLENGTH(ranges) != n || XLENGTH(idxs) != n)
    error("msutils: length mismatch for chrom parameters");

  REQUIRE_BOUND(ABI.get_peaks_from_chrom, "get_peaks_from_chrom");
  REQUIRE_BOUND(ABI.free_, "free_");

  uint32_t *uidx = (uint32_t *)R_alloc((size_t)n, sizeof(uint32_t));
  if (TYPEOF(idxs) == INTSXP)
  {
    int *ix = INTEGER(idxs);
    for (R_xlen_t i = 0; i < n; i++)
    {
      int v = ix[i];
      uidx[i] = (v == NA_INTEGER || v < 0) ? UINT32_MAX : (uint32_t)v;
    }
  }
  else if (TYPEOF(idxs) == REALSXP)
  {
    double *dx = REAL(idxs);
    for (R_xlen_t i = 0; i < n; i++)
    {
      double v = dx[i];
      uidx[i] = (!R_finite(v) || v < 0) ? UINT32_MAX : (uint32_t)v;
    }
  }
  else
  {
    error("msutils: idx must be integer or numeric");
  }

  size_t ncores = (cores == R_NilValue) ? 1 : (size_t)asInteger(cores);
  if (ncores < 1)
    ncores = 1;

  CPeakPOptions opts;
  const CPeakPOptions *opt_ptr = NULL;
  as_opts_ptr(options, &opts, &opt_ptr);

  Buf out = (Buf){0};
  int code = ABI.get_peaks_from_chrom(
      handle,
      uidx, REAL(rts), REAL(ranges), (size_t)n, opt_ptr, ncores, &out);

  die_code("get_peaks_from_chrom", code);
  SEXP res = mk_string_len(out.ptr, out.len);
  ABI.free_(out.ptr, out.len);
  return res;
}

SEXP C_calculate_eic(SEXP bin, SEXP targets, SEXP from, SEXP to, SEXP ppm_tol, SEXP mz_tol)
{
  if (TYPEOF(bin) != RAWSXP)
    error("bin");
  if ((TYPEOF(targets) != REALSXP && TYPEOF(targets) != INTSXP) || LENGTH(targets) != 1)
    error("targets");
  REQUIRE_BOUND(ABI.calculate_eic, "calculate_eic");
  REQUIRE_BOUND(ABI.free_, "free_");
  double t = asReal(targets);
  Buf bx = (Buf){0}, by = (Buf){0};
  int code = ABI.calculate_eic(
      (const unsigned char *)RAW(bin), (size_t)XLENGTH(bin),
      t, asReal(from), asReal(to), asReal(ppm_tol), asReal(mz_tol),
      &bx, &by);
  die_code("calculate_eic", code);
  size_t nx = bx.len / 8, ny = by.len / 8;
  SEXP Rx = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)nx));
  SEXP Ry = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)ny));
  memcpy(REAL(Rx), bx.ptr, bx.len);
  memcpy(REAL(Ry), by.ptr, by.len);
  ABI.free_(bx.ptr, bx.len);
  ABI.free_(by.ptr, by.len);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, Rx);
  SET_VECTOR_ELT(out, 1, Ry);
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(nms, 0, Rf_mkChar("x"));
  SET_STRING_ELT(nms, 1, Rf_mkChar("y"));
  Rf_setAttrib(out, R_NamesSymbol, nms);
  UNPROTECT(4);
  return out;
}

SEXP C_find_peaks(SEXP x, SEXP y, SEXP options)
{
  if (TYPEOF(x) != REALSXP || TYPEOF(y) != REALSXP)
    error("numeric");
  if (XLENGTH(x) != XLENGTH(y) || XLENGTH(x) < 3)
    error("length");
  REQUIRE_BOUND(ABI.find_peaks, "find_peaks");
  REQUIRE_BOUND(ABI.free_, "free_");
  R_xlen_t n = XLENGTH(y);
  CPeakPOptions opts;
  const CPeakPOptions *opt_ptr = NULL;
  (void)as_opts_ptr(options, &opts, &opt_ptr);
  Buf out = (Buf){0};
  int code = ABI.find_peaks(REAL(x), REAL(y), (size_t)n, opt_ptr, &out);
  die_code("find_peaks", code);
  SEXP res = mk_string_len(out.ptr, out.len);
  ABI.free_(out.ptr, out.len);
  return res;
}

SEXP C_calculate_baseline(SEXP y, SEXP baseline_window, SEXP baseline_window_factor)
{
  if (TYPEOF(y) != REALSXP)
    error("numeric y required");
  REQUIRE_BOUND(ABI.calculate_baseline, "calculate_baseline");
  REQUIRE_BOUND(ABI.free_, "free_");
  R_xlen_t n = XLENGTH(y);
  int bw = asInteger(baseline_window);
  int bf = asInteger(baseline_window_factor);
  Buf out = (Buf){0};
  int code = ABI.calculate_baseline(REAL(y), (size_t)n, (int32_t)bw, (int32_t)bf, &out);
  die_code("calculate_baseline", code);
  size_t m = out.len / 8;
  SEXP Ry = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)m));
  memcpy(REAL(Ry), out.ptr, out.len);
  ABI.free_(out.ptr, out.len);
  UNPROTECT(1);
  return Ry;
}

SEXP C_find_features(SEXP data, SEXP from_time, SEXP to_time, SEXP eic_ppm_tol, SEXP eic_mz_tol, SEXP grid_start, SEXP grid_end, SEXP grid_step_ppm, SEXP options, SEXP cores)
{
  MzML *handle = GetHandle(data);
  REQUIRE_BOUND(ABI.find_features, "find_features");
  REQUIRE_BOUND(ABI.free_, "free_");

  int ncores = asInteger(cores);
  if (ncores < 1)
    ncores = 1;

  CPeakPOptions opts;
  const CPeakPOptions *opt_ptr = NULL;
  as_opts_ptr(options, &opts, &opt_ptr);

  Buf out = (Buf){0};
  int code = ABI.find_features(
      handle,
      asReal(from_time), asReal(to_time),
      asReal(eic_ppm_tol), asReal(eic_mz_tol),
      asReal(grid_start), asReal(grid_end),
      asReal(grid_step_ppm),
      opt_ptr, (int32_t)ncores, &out);

  die_code("find_features", code);
  SEXP res = mk_string_len(out.ptr, out.len);
  ABI.free_(out.ptr, out.len);
  return res;
}

SEXP C_find_feature(SEXP bin, SEXP rts, SEXP mzs, SEXP wins, SEXP ids, SEXP scan_ppm, SEXP scan_mz, SEXP eic_ppm, SEXP eic_mz, SEXP options, SEXP cores)
{
  MzML *handle = GetHandle(bin);
  REQUIRE_BOUND(ABI.find_feature, "find_feature");
  REQUIRE_BOUND(ABI.free_, "free_");

  if (TYPEOF(rts) != REALSXP || TYPEOF(mzs) != REALSXP || TYPEOF(wins) != REALSXP)
    error("msutils: numeric (double) required for rts/mzs/wins");

  R_xlen_t n = XLENGTH(rts);
  if (XLENGTH(mzs) != n || XLENGTH(wins) != n)
    error("msutils: length mismatch for rts/mzs/wins");

  uint32_t *offs = (uint32_t *)R_alloc((size_t)n, sizeof(uint32_t));
  uint32_t *lens = (uint32_t *)R_alloc((size_t)n, sizeof(uint32_t));
  unsigned char *ids_buf = NULL;
  size_t ids_len = 0;

  if (ids != R_NilValue)
  {
    if (TYPEOF(ids) != STRSXP)
      error("msutils: ids must be character");
    size_t total = 0;
    for (R_xlen_t i = 0; i < n; i++)
    {
      SEXP s = STRING_ELT(ids, i);
      if (s != R_NilValue)
        total += (size_t)LENGTH(s);
    }
    ids_buf = (unsigned char *)R_alloc(total, 1);
    ids_len = total;
    size_t cur = 0;
    for (R_xlen_t i = 0; i < n; i++)
    {
      SEXP s = STRING_ELT(ids, i);
      if (s == R_NilValue)
      {
        offs[i] = 0;
        lens[i] = 0;
      }
      else
      {
        size_t L = (size_t)LENGTH(s);
        offs[i] = (uint32_t)cur;
        lens[i] = (uint32_t)L;
        memcpy(ids_buf + cur, CHAR(s), L);
        cur += L;
      }
    }
  }

  size_t ncores = (cores == R_NilValue) ? 1 : (size_t)asInteger(cores);
  CPeakPOptions opts;
  const CPeakPOptions *opt_ptr = NULL;
  as_opts_ptr(options, &opts, &opt_ptr);

  Buf out = (Buf){0};
  int32_t code = ABI.find_feature(
      handle, REAL(rts), REAL(mzs), REAL(wins),
      offs, lens, ids_buf, ids_len,
      (size_t)n, ncores,
      asReal(scan_ppm), asReal(scan_mz), asReal(eic_ppm), asReal(eic_mz),
      opt_ptr, &out);

  die_code("find_feature", code);
  SEXP res = mk_string_len(out.ptr, out.len);
  ABI.free_(out.ptr, out.len);
  return res;
}

SEXP C_convert_mzml_to_bin(SEXP xml, SEXP level, SEXP f32_compress)
{
  if (TYPEOF(xml) != RAWSXP)
    error("xml must be a raw vector");
  if (!(TYPEOF(level) == INTSXP || TYPEOF(level) == REALSXP) || LENGTH(level) != 1)
    error("level must be a scalar number");
  int lv = asInteger(level);
  if (lv < 0 || lv > 22)
    error("level must be in [0,22]");

  int fc = asLogical(f32_compress);
  if (fc == NA_LOGICAL)
    error("f32_compress must be TRUE/FALSE");

  REQUIRE_BOUND(ABI.convert_mzml_to_bin, "convert_mzml_to_bin");
  REQUIRE_BOUND(ABI.free_, "free_");

  Buf out = (Buf){0};
  int32_t code = ABI.convert_mzml_to_bin(
      (const unsigned char *)RAW(xml),
      (size_t)XLENGTH(xml),
      &out,
      (uint8_t)lv,
      (uint8_t)(fc ? 1 : 0));

  die_code("convert_mzml_to_bin", code);

  SEXP res = PROTECT(Rf_allocVector(RAWSXP, (R_xlen_t)out.len));
  memcpy(RAW(res), out.ptr, out.len);
  ABI.free_(out.ptr, out.len);
  UNPROTECT(1);
  return res;
}

SEXP C_parse_bin(SEXP bin)
{
  if (TYPEOF(bin) != RAWSXP)
    error("msutils: data must be a raw vector");
  REQUIRE_BOUND(ABI.parse_bin, "parse_bin");

  MzML *handle = NULL;
  int code = ABI.parse_bin((const unsigned char *)RAW(bin), (size_t)XLENGTH(bin), &handle);
  die_code("parse_bin", code);

  SEXP ptr = PROTECT(R_MakeExternalPtr(handle, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(ptr, finalize_mzml, TRUE);
  UNPROTECT(1);
  return ptr;
}
