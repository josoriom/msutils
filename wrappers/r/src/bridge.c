#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "include/quantion.h"

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
static const char *last_err = "dlopen or dlsym failed";
#endif

typedef struct MzML MzML;

typedef uint32_t (*fn_quantion_abi_version)(void);
typedef size_t (*fn_quantion_sizeof_peak_options)(void);
typedef int32_t (*fn_parse_mzml)(const unsigned char *, size_t, MzML **);
typedef int32_t (*fn_bin_to_json)(const MzML *, Buf *);
typedef int32_t (*fn_bin_to_mzml)(const MzML *, Buf *);
typedef int32_t (*fn_get_peak)(const double *, const double *, size_t, double, double, const CPeakOptions *, Buf *);
typedef int32_t (*fn_fit_peak)(const double *, const double *, size_t, double, double, int32_t, Buf *);
typedef int32_t (*fn_draw_peak)(const double *, size_t, int32_t, double, double, double, double, Buf *);
typedef int32_t (*fn_calculate_eic)(const MzML *, double, double, double, double, double, Buf *, Buf *);
typedef int32_t (*fn_find_noise_level)(const float *, size_t, size_t *, double *);
typedef int32_t (*fn_get_peaks_from_eic)(const MzML *, const double *, const double *, const double *, const uint32_t *, const uint32_t *, const unsigned char *, size_t, size_t, double, double, const CPeakOptions *, size_t, Buf *);
typedef int32_t (*fn_get_peaks_from_chrom)(const MzML *, const uint32_t *, const double *, const double *, size_t, const CPeakOptions *, size_t, Buf *);
typedef int32_t (*fn_find_peaks)(const double *, const double *, size_t, const CPeakOptions *, Buf *);
typedef int32_t (*fn_calculate_baseline)(const double *, size_t, int32_t, int32_t, Buf *);
typedef int32_t (*fn_find_features)(const MzML *, double, double, double, double, double, double, double, const CPeakOptions *, int32_t, Buf *);
typedef int32_t (*fn_find_feature)(const MzML *, const double *, const double *, const double *, const uint32_t *, const uint32_t *, const unsigned char *, size_t, size_t, size_t, double, double, double, double, const CPeakOptions *, Buf *);
typedef int32_t (*fn_mzml_to_bin)(const MzML *, Buf *, uint8_t, uint8_t);
typedef int32_t (*fn_convert_mzml_file_to_ion_file)(const char *, const char *, uint8_t, uint8_t);
typedef int32_t (*fn_parse_bin)(const unsigned char *, size_t, size_t, MzML **);
typedef int32_t (*fn_read_range)(void *, uint64_t, uint64_t, uint8_t *);
typedef int32_t (*fn_parse_ion_source)(fn_read_range, void *, size_t, MzML **);
typedef int32_t (*fn_plan_open)(const uint8_t *, size_t, Buf *);
typedef int32_t (*fn_plan_eic)(MzML *, double, double, double, double, double, Buf *);
typedef int32_t (*fn_parse_ion_path)(const char *, size_t, MzML **);
typedef int32_t (*fn_get_features)(const char *, double, double, double, double, double, double, double, double, double, double, int32_t, const CPeakOptions *, int32_t, Buf *);
typedef int32_t (*fn_get_scans)(const MzML *, uint8_t, double, double, uint8_t, Buf *);
typedef void (*fn_free_)(unsigned char *, size_t);
typedef void (*fn_free_mzml)(MzML *);

typedef struct
{
  fn_quantion_abi_version quantion_abi_version;
  fn_quantion_sizeof_peak_options quantion_sizeof_peak_options;
  fn_parse_mzml parse_mzml;
  fn_bin_to_json bin_to_json;
  fn_bin_to_mzml bin_to_mzml;
  fn_get_peak get_peak;
  fn_fit_peak fit_peak;
  fn_draw_peak draw_peak;
  fn_calculate_eic calculate_eic;
  fn_find_noise_level find_noise_level;
  fn_get_peaks_from_eic get_peaks_from_eic;
  fn_get_peaks_from_chrom get_peaks_from_chrom;
  fn_find_peaks find_peaks;
  fn_calculate_baseline calculate_baseline;
  fn_find_features find_features;
  fn_free_ free_;
  fn_find_feature find_feature;
  fn_mzml_to_bin mzml_to_bin;
  fn_convert_mzml_file_to_ion_file convert_mzml_file_to_ion_file;
  fn_parse_bin parse_bin;
  fn_parse_ion_source parse_ion_source;
  fn_plan_open plan_open;
  fn_plan_eic plan_eic;
  fn_parse_ion_path parse_ion_path;
  fn_free_mzml free_mzml;
  fn_get_features get_features;
  fn_get_scans get_scans;
} abi_type;

static DLIB abi_handle = NULL;
abi_type ABI = (abi_type){0};

static int resolve_required(void **fn, const char *name)
{
#if !defined(_WIN32)
  dlerror();
#endif
  *fn = DLSYM(abi_handle, name);
  if (!*fn)
  {
#if !defined(_WIN32)
    last_err = dlerror();
#endif
    return -1;
  }
  return 0;
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
  if (resolve_required((void **)&ABI.quantion_abi_version, "quantion_abi_version"))
    goto fail;
  if (ABI.quantion_abi_version() != QUANTION_ABI_VERSION)
  {
    last_err = "ABI version mismatch: recompile the R wrapper";
    goto fail;
  }
  if (resolve_required((void **)&ABI.quantion_sizeof_peak_options, "quantion_sizeof_peak_options"))
    goto fail;
  if (ABI.quantion_sizeof_peak_options() != sizeof(CPeakOptions))
  {
    last_err = "PeakOptions size mismatch: native binary and R wrapper are out of sync — reinstall quantion";
    goto fail;
  }
  if (resolve_required((void **)&ABI.parse_mzml, "parse_mzml"))
    goto fail;
  if (resolve_required((void **)&ABI.bin_to_json, "bin_to_json"))
    goto fail;
  if (resolve_required((void **)&ABI.bin_to_mzml, "bin_to_mzml"))
    goto fail;
  if (resolve_required((void **)&ABI.get_peak, "get_peak"))
    goto fail;
  if (resolve_required((void **)&ABI.fit_peak, "fit_peak"))
    goto fail;
  if (resolve_required((void **)&ABI.draw_peak, "draw_peak"))
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
  if (resolve_required((void **)&ABI.mzml_to_bin, "mzml_to_bin"))
    goto fail;
  if (resolve_required((void **)&ABI.convert_mzml_file_to_ion_file, "convert_mzml_file_to_ion_file"))
    goto fail;
  ABI.parse_ion_source = (fn_parse_ion_source)DLSYM(abi_handle, "parse_ion_source");
  ABI.plan_open = (fn_plan_open)DLSYM(abi_handle, "plan_open");
  ABI.plan_eic = (fn_plan_eic)DLSYM(abi_handle, "plan_eic");
  if (resolve_required((void **)&ABI.parse_bin, "parse_bin"))
    goto fail;
  if (resolve_required((void **)&ABI.parse_ion_path, "parse_ion_path"))
    goto fail;
  if (resolve_required((void **)&ABI.free_mzml, "free_mzml"))
    goto fail;
  if (resolve_required((void **)&ABI.get_features, "get_features"))
    goto fail;
  if (resolve_required((void **)&ABI.get_scans, "get_scans"))
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
  else if (code == 5)
    msg = "encode error";
  else if (code == 6)
    msg = "fast EIC path unavailable: this .ion file has no usable spectrum bounds (A3); re-encode it with the current Ionic to use the fast EIC path";
  error("quantion/%s failed: %s (code=%d)", fname, msg, code);
}

struct mk_string_data
{
  const unsigned char *ptr;
  size_t len;
};

static SEXP mk_string_build(void *data)
{
  struct mk_string_data *d = (struct mk_string_data *)data;
  return Rf_ScalarString(Rf_mkCharLenCE((const char *)d->ptr, (int)d->len, CE_UTF8));
}

static void mk_string_free(void *data, Rboolean jump)
{
  (void)jump;
  struct mk_string_data *d = (struct mk_string_data *)data;
  if (d->ptr && ABI.free_)
    ABI.free_((unsigned char *)d->ptr, d->len);
  d->ptr = NULL;
}

static SEXP mk_string_len(const unsigned char *ptr, size_t len)
{
  if (len > (size_t)INT_MAX)
  {
    if (ABI.free_)
      ABI.free_((unsigned char *)ptr, len);
    error("quantion: result too large for an R string (%zu bytes)", len);
  }
  struct mk_string_data d = {ptr, len};
  SEXP cont = PROTECT(R_MakeUnwindCont());
  SEXP s = R_UnwindProtect(mk_string_build, &d, mk_string_free, &d, cont);
  UNPROTECT(1);
  return s;
}

#define REQUIRE_BOUND(ptr, name)                                            \
  do                                                                        \
  {                                                                         \
    if ((ptr) == NULL)                                                      \
      error("quantion: symbol %s is not bound; did .onLoad() run?", (name)); \
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

static int fill_options(SEXP opts, CPeakOptions *out)
{
  if (opts == R_NilValue || TYPEOF(opts) != VECSXP || XLENGTH(opts) == 0)
    return 0;
  out->min_integral = NAN;
  out->min_intensity = NAN;
  out->min_peak_width_points = 0;
  out->noise = NAN;
  out->auto_noise = 0;
  out->auto_baseline = 0;
  out->lambda = 0;
  out->max_iterations = 0;
  out->allow_overlap = 0;
  out->min_snr = NAN;
  out->min_r2 = NAN;
  out->shape = 1;
  out->kernel_size = 0;
  SEXP v = R_NilValue;
  v = list_get(opts, "min_integral");
  if (v != R_NilValue)
    out->min_integral = asReal(v);
  v = list_get(opts, "min_intensity");
  if (v != R_NilValue)
    out->min_intensity = asReal(v);
  v = list_get(opts, "min_peak_width_points");
  if (v != R_NilValue)
    out->min_peak_width_points = (int32_t)asInteger(v);
  v = list_get(opts, "noise");
  if (v != R_NilValue)
    out->noise = asReal(v);
  v = list_get(opts, "auto_noise");
  if (v != R_NilValue)
    out->auto_noise = (int32_t)asLogical(v);
  v = list_get(opts, "auto_baseline");
  if (v != R_NilValue)
    out->auto_baseline = (int32_t)asLogical(v);
  v = list_get(opts, "lambda");
  if (v != R_NilValue)
    out->lambda = (int32_t)asInteger(v);
  v = list_get(opts, "max_iterations");
  if (v != R_NilValue)
    out->max_iterations = (int32_t)asInteger(v);
  v = list_get(opts, "allow_overlap");
  if (v != R_NilValue)
    out->allow_overlap = (int32_t)asLogical(v);
  v = list_get(opts, "min_snr");
  if (v != R_NilValue)
    out->min_snr = asReal(v);
  v = list_get(opts, "min_r2");
  if (v != R_NilValue)
    out->min_r2 = asReal(v);
  v = list_get(opts, "shape");
  if (v != R_NilValue)
  {
    if (TYPEOF(v) == STRSXP && XLENGTH(v) > 0)
      out->shape = (strcmp(CHAR(STRING_ELT(v, 0)), "gaussian") == 0) ? 0 : 1;
    else
      out->shape = (int32_t)asInteger(v);
  }
  v = list_get(opts, "kernel_size");
  if (v != R_NilValue)
    out->kernel_size = (int32_t)asInteger(v);
  return 1;
}

static int as_opts_ptr(SEXP options, CPeakOptions *copy, const CPeakOptions **out_ptr)
{
  *out_ptr = NULL;
  if (options == R_NilValue)
    return 0;
  if (TYPEOF(options) == RAWSXP)
  {
    if ((size_t)XLENGTH(options) != sizeof(CPeakOptions))
      error("quantion: options raw blob must be length %zu", sizeof(CPeakOptions));
    memcpy((void *)copy, (const void *)RAW(options), sizeof(CPeakOptions));
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
    error("quantion: expected ExternalPtr");
  MzML *h = (MzML *)R_ExternalPtrAddr(ptr);
  if (!h)
    error("quantion: use of disposed or null pointer");
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

SEXP C_unbind_rust(SEXP unused)
{
  abi_unload();
  return R_NilValue;
}

typedef struct
{
  uint64_t offset;
  uint64_t length;
  uint8_t *bytes;
} held_range;

typedef struct
{
  held_range *ranges;
  size_t count;
  size_t room;
} range_store;

static range_store *store_new(void)
{
  range_store *store = (range_store *)calloc(1, sizeof(range_store));
  return store;
}

static void store_free(range_store *store)
{
  if (store == NULL)
    return;
  for (size_t index = 0; index < store->count; index++)
    free(store->ranges[index].bytes);
  free(store->ranges);
  free(store);
}

static int store_add(range_store *store, uint64_t offset, const uint8_t *bytes, uint64_t length)
{
  if (store->count == store->room)
  {
    size_t room = store->room == 0 ? 8 : store->room * 2;
    held_range *grown = (held_range *)realloc(store->ranges, room * sizeof(held_range));
    if (grown == NULL)
      return -1;
    store->ranges = grown;
    store->room = room;
  }
  uint8_t *copy = (uint8_t *)malloc(length > 0 ? (size_t)length : 1);
  if (copy == NULL)
    return -1;
  memcpy(copy, bytes, (size_t)length);
  store->ranges[store->count].offset = offset;
  store->ranges[store->count].length = length;
  store->ranges[store->count].bytes = copy;
  store->count++;
  return 0;
}

static int32_t store_serve(void *context, uint64_t offset, uint64_t length, uint8_t *dest)
{
  range_store *store = (range_store *)context;
  if (store == NULL)
    return -1;
  for (size_t index = 0; index < store->count; index++)
  {
    held_range *range = &store->ranges[index];
    if (offset < range->offset)
      continue;
    uint64_t inside = offset - range->offset;
    if (inside > range->length || length > range->length - inside)
      continue;
    memcpy(dest, range->bytes + inside, (size_t)length);
    return 0;
  }
  return -1;
}

static void finalize_store(SEXP ptr)
{
  if (R_ExternalPtrAddr(ptr))
  {
    store_free((range_store *)R_ExternalPtrAddr(ptr));
    R_ClearExternalPtr(ptr);
  }
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

SEXP C_ion_to_json(SEXP bin)
{
  MzML *handle = GetHandle(bin);
  REQUIRE_BOUND(ABI.bin_to_json, "bin_to_json");
  REQUIRE_BOUND(ABI.free_, "free_");
  Buf out = (Buf){0};
  int code = ABI.bin_to_json(handle, &out);
  die_code("bin_to_json", code);
  SEXP res = mk_string_len(out.ptr, out.len);
  return res;
}

SEXP C_ion_to_mzml(SEXP bin)
{
  MzML *handle = GetHandle(bin);
  REQUIRE_BOUND(ABI.bin_to_mzml, "bin_to_mzml");
  REQUIRE_BOUND(ABI.free_, "free_");
  Buf out = (Buf){0};
  int code = ABI.bin_to_mzml(handle, &out);
  die_code("bin_to_mzml", code);
  SEXP res = mk_string_len(out.ptr, out.len);
  return res;
}

SEXP C_mzml_to_ion(SEXP bin, SEXP level, SEXP f32_compress)
{
  MzML *handle = GetHandle(bin);
  if (!(TYPEOF(level) == INTSXP || TYPEOF(level) == REALSXP) || LENGTH(level) != 1)
    error("level must be a scalar number");
  int lv = asInteger(level);
  if (lv < 0 || lv > 22)
    error("level must be in [0,22]");
  int fc = asLogical(f32_compress);
  if (fc == NA_LOGICAL)
    error("f32_compress must be TRUE/FALSE");
  REQUIRE_BOUND(ABI.mzml_to_bin, "mzml_to_bin");
  REQUIRE_BOUND(ABI.free_, "free_");
  Buf out = (Buf){0};
  int32_t code = ABI.mzml_to_bin(
      handle,
      &out,
      (uint8_t)lv,
      (uint8_t)(fc ? 1 : 0));
  die_code("mzml_to_bin", code);
  SEXP res = PROTECT(Rf_allocVector(RAWSXP, (R_xlen_t)out.len));
  memcpy(RAW(res), out.ptr, out.len);
  ABI.free_(out.ptr, out.len);
  UNPROTECT(1);
  return res;
}

SEXP C_mzml_to_ion_file(SEXP input_path, SEXP output_path, SEXP level, SEXP f32_compress)
{
  if (TYPEOF(input_path) != STRSXP || LENGTH(input_path) != 1)
    error("input_path must be a length-1 character string");
  if (TYPEOF(output_path) != STRSXP || LENGTH(output_path) != 1)
    error("output_path must be a length-1 character string");
  if (!(TYPEOF(level) == INTSXP || TYPEOF(level) == REALSXP) || LENGTH(level) != 1)
    error("level must be a scalar number");
  int lv = asInteger(level);
  if (lv < 0 || lv > 22)
    error("level must be in [0,22]");
  int fc = asLogical(f32_compress);
  if (fc == NA_LOGICAL)
    error("f32_compress must be TRUE/FALSE");
  REQUIRE_BOUND(ABI.convert_mzml_file_to_ion_file, "convert_mzml_file_to_ion_file");
  const char *in_path = CHAR(STRING_ELT(input_path, 0));
  const char *out_path = CHAR(STRING_ELT(output_path, 0));
  int32_t code = ABI.convert_mzml_file_to_ion_file(
      in_path, out_path, (uint8_t)lv, (uint8_t)(fc ? 1 : 0));
  die_code("convert_mzml_file_to_ion_file", code);
  return R_NilValue;
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
  CPeakOptions opts;
  const CPeakOptions *opt_ptr = NULL;
  (void)as_opts_ptr(options, &opts, &opt_ptr);
  Buf out = (Buf){0};
  int code = ABI.get_peak(REAL(x), REAL(y), (size_t)n, asReal(rt), asReal(range), opt_ptr, &out);
  die_code("get_peak", code);
  SEXP res = mk_string_len(out.ptr, out.len);
  return res;
}

SEXP C_fit_peak(SEXP x, SEXP y, SEXP rt, SEXP intensity, SEXP shape)
{
  if (TYPEOF(x) != REALSXP || TYPEOF(y) != REALSXP)
    error("numeric");
  if (XLENGTH(x) != XLENGTH(y) || XLENGTH(x) < 5)
    error("length");
  REQUIRE_BOUND(ABI.fit_peak, "fit_peak");
  REQUIRE_BOUND(ABI.free_, "free_");
  R_xlen_t n = XLENGTH(y);
  Buf out = (Buf){0};
  int code = ABI.fit_peak(REAL(x), REAL(y), (size_t)n,
                          asReal(rt), asReal(intensity), (int32_t)asInteger(shape), &out);
  die_code("fit_peak", code);
  SEXP res = mk_string_len(out.ptr, out.len);
  return res;
}

SEXP C_draw_peak(SEXP x, SEXP shape, SEXP height, SEXP center, SEXP fwhm, SEXP tail)
{
  if (TYPEOF(x) != REALSXP)
    error("numeric");
  REQUIRE_BOUND(ABI.draw_peak, "draw_peak");
  REQUIRE_BOUND(ABI.free_, "free_");
  R_xlen_t n = XLENGTH(x);
  Buf out = (Buf){0};
  int code = ABI.draw_peak(REAL(x), (size_t)n, (int32_t)asInteger(shape),
                           asReal(height), asReal(center), asReal(fwhm), asReal(tail), &out);
  die_code("draw_peak", code);
  size_t ny = out.len / 8;
  SEXP Ry = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)ny));
  memcpy(REAL(Ry), out.ptr, ny * sizeof(double));
  ABI.free_(out.ptr, out.len);
  UNPROTECT(1);
  return Ry;
}

SEXP C_get_peaks_from_eic(SEXP bin, SEXP rts, SEXP mzs, SEXP ranges, SEXP ids, SEXP from_left, SEXP to_right, SEXP options, SEXP cores)
{
  MzML *handle = GetHandle(bin);
  REQUIRE_BOUND(ABI.get_peaks_from_eic, "get_peaks_from_eic");
  REQUIRE_BOUND(ABI.free_, "free_");

  if (TYPEOF(rts) != REALSXP || TYPEOF(mzs) != REALSXP || TYPEOF(ranges) != REALSXP)
    error("quantion: bad numeric arguments");
  R_xlen_t n = XLENGTH(rts);
  if (XLENGTH(mzs) != n || XLENGTH(ranges) != n)
    error("quantion: length mismatch");
  uint32_t *offs = (uint32_t *)R_alloc((size_t)n, sizeof(uint32_t));
  uint32_t *lens = (uint32_t *)R_alloc((size_t)n, sizeof(uint32_t));
  unsigned char *ids_buf = NULL;
  size_t ids_len = 0;
  if (ids != R_NilValue)
  {
    if (TYPEOF(ids) != STRSXP)
      error("quantion: ids must be character");
    if (XLENGTH(ids) != n)
      error("quantion: ids must have the same length as rt");
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
  if (ncores < 1)
    ncores = 1;
  CPeakOptions opts;
  const CPeakOptions *opt_ptr = NULL;
  as_opts_ptr(options, &opts, &opt_ptr);

  Buf out = (Buf){0};
  int code = ABI.get_peaks_from_eic(
      handle, REAL(rts), REAL(mzs), REAL(ranges),
      offs, lens, ids_buf, ids_len,
      (size_t)n, asReal(from_left), asReal(to_right), opt_ptr, ncores, &out);

  die_code("get_peaks_from_eic", code);
  SEXP res = mk_string_len(out.ptr, out.len);
  return res;
}

SEXP C_get_peaks_from_chrom(SEXP bin, SEXP idxs, SEXP rts, SEXP ranges, SEXP options, SEXP cores)
{
  MzML *handle = GetHandle(bin);

  if (TYPEOF(rts) != REALSXP || TYPEOF(ranges) != REALSXP)
    error("quantion: numeric (double) required for rts/ranges");

  R_xlen_t n = XLENGTH(rts);
  if (XLENGTH(ranges) != n || XLENGTH(idxs) != n)
    error("quantion: length mismatch for chrom parameters");

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
    error("quantion: idx must be integer or numeric");
  }

  size_t ncores = (cores == R_NilValue) ? 1 : (size_t)asInteger(cores);
  if (ncores < 1)
    ncores = 1;

  CPeakOptions opts;
  const CPeakOptions *opt_ptr = NULL;
  as_opts_ptr(options, &opts, &opt_ptr);

  Buf out = (Buf){0};
  int code = ABI.get_peaks_from_chrom(
      handle,
      uidx, REAL(rts), REAL(ranges), (size_t)n, opt_ptr, ncores, &out);

  die_code("get_peaks_from_chrom", code);
  SEXP res = mk_string_len(out.ptr, out.len);
  return res;
}

SEXP C_calculate_eic(SEXP bin, SEXP targets, SEXP from, SEXP to, SEXP ppm_tol, SEXP mz_tol)
{
  MzML *handle = GetHandle(bin);
  if ((TYPEOF(targets) != REALSXP && TYPEOF(targets) != INTSXP) || LENGTH(targets) != 1)
    error("targets");
  REQUIRE_BOUND(ABI.calculate_eic, "calculate_eic");
  REQUIRE_BOUND(ABI.free_, "free_");
  double t = asReal(targets);
  Buf bx = (Buf){0}, by = (Buf){0};
  int code = ABI.calculate_eic(
      handle,
      t, asReal(from), asReal(to), asReal(ppm_tol), asReal(mz_tol),
      &bx, &by);
  die_code("calculate_eic", code);
  size_t nx = bx.len / 8, ny = by.len / 8;
  SEXP Rx = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)nx));
  SEXP Ry = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)ny));
  memcpy(REAL(Rx), bx.ptr, nx * sizeof(double));
  memcpy(REAL(Ry), by.ptr, ny * sizeof(double));
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
  CPeakOptions opts;
  const CPeakOptions *opt_ptr = NULL;
  (void)as_opts_ptr(options, &opts, &opt_ptr);
  Buf out = (Buf){0};
  int code = ABI.find_peaks(REAL(x), REAL(y), (size_t)n, opt_ptr, &out);
  die_code("find_peaks", code);
  SEXP res = mk_string_len(out.ptr, out.len);
  return res;
}

SEXP C_calculate_baseline(SEXP y, SEXP lambda, SEXP max_iterations)
{
  if (TYPEOF(y) != REALSXP)
    error("numeric y required");
  REQUIRE_BOUND(ABI.calculate_baseline, "calculate_baseline");
  REQUIRE_BOUND(ABI.free_, "free_");
  R_xlen_t n = XLENGTH(y);
  int lam = asInteger(lambda);
  int maxit = asInteger(max_iterations);
  Buf out = (Buf){0};
  int code = ABI.calculate_baseline(REAL(y), (size_t)n, (int32_t)lam, (int32_t)maxit, &out);
  die_code("calculate_baseline", code);
  size_t m = out.len / 8;
  SEXP Ry = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)m));
  memcpy(REAL(Ry), out.ptr, m * sizeof(double));
  ABI.free_(out.ptr, out.len);
  UNPROTECT(1);
  return Ry;
}

SEXP C_find_noise_level(SEXP y)
{
  REQUIRE_BOUND(ABI.find_noise_level, "find_noise_level");
  R_xlen_t n = XLENGTH(y);
  const double *yd = REAL(y);
  float *yf = (float *)R_alloc((size_t)n, sizeof(float));
  for (R_xlen_t i = 0; i < n; i++)
    yf[i] = (float)yd[i];
  size_t width = 0;
  double intensity = 0.0;
  int code = ABI.find_noise_level(yf, (size_t)n, &width, &intensity);
  die_code("find_noise_level", code);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, Rf_ScalarInteger((int)width));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(intensity));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(nms, 0, Rf_mkChar("width"));
  SET_STRING_ELT(nms, 1, Rf_mkChar("intensity"));
  Rf_setAttrib(out, R_NamesSymbol, nms);
  UNPROTECT(2);
  return out;
}

SEXP C_find_features(SEXP data, SEXP from_time, SEXP to_time, SEXP eic_ppm_tol, SEXP eic_mz_tol, SEXP grid_start, SEXP grid_end, SEXP grid_step_ppm, SEXP options, SEXP cores)
{
  MzML *handle = GetHandle(data);
  REQUIRE_BOUND(ABI.find_features, "find_features");
  REQUIRE_BOUND(ABI.free_, "free_");

  int ncores = asInteger(cores);
  if (ncores < 1)
    ncores = 1;

  CPeakOptions opts;
  const CPeakOptions *opt_ptr = NULL;
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
  return res;
}

SEXP C_find_feature(SEXP bin, SEXP rts, SEXP mzs, SEXP wins, SEXP ids, SEXP scan_ppm, SEXP scan_mz, SEXP eic_ppm, SEXP eic_mz, SEXP options, SEXP cores)
{
  MzML *handle = GetHandle(bin);
  REQUIRE_BOUND(ABI.find_feature, "find_feature");
  REQUIRE_BOUND(ABI.free_, "free_");

  if (TYPEOF(rts) != REALSXP || TYPEOF(mzs) != REALSXP || TYPEOF(wins) != REALSXP)
    error("quantion: numeric (double) required for rts/mzs/wins");

  R_xlen_t n = XLENGTH(rts);
  if (XLENGTH(mzs) != n || XLENGTH(wins) != n)
    error("quantion: length mismatch for rts/mzs/wins");

  uint32_t *offs = (uint32_t *)R_alloc((size_t)n, sizeof(uint32_t));
  uint32_t *lens = (uint32_t *)R_alloc((size_t)n, sizeof(uint32_t));
  unsigned char *ids_buf = NULL;
  size_t ids_len = 0;

  if (ids != R_NilValue)
  {
    if (TYPEOF(ids) != STRSXP)
      error("quantion: ids must be character");
    if (XLENGTH(ids) != n)
      error("quantion: ids must have the same length as rt");
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
  if (ncores < 1)
    ncores = 1;
  CPeakOptions opts;
  const CPeakOptions *opt_ptr = NULL;
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
  return res;
}

static SEXP ranges_to_matrix(Buf *out)
{
  size_t count = out->len / 16;
  SEXP result = PROTECT(allocMatrix(REALSXP, (int)count, 2));
  double *values = REAL(result);
  for (size_t index = 0; index < count; index++)
  {
    uint64_t offset = 0;
    uint64_t length = 0;
    memcpy(&offset, out->ptr + index * 16, 8);
    memcpy(&length, out->ptr + index * 16 + 8, 8);
    values[index] = (double)offset;
    values[count + index] = (double)length;
  }
  if (ABI.free_ && out->ptr)
    ABI.free_(out->ptr, out->len);
  UNPROTECT(1);
  return result;
}

SEXP C_plan_open(SEXP header)
{
  if (TYPEOF(header) != RAWSXP)
    error("quantion: header must be a raw vector");
  REQUIRE_BOUND(ABI.plan_open, "plan_open");

  Buf out;
  memset(&out, 0, sizeof(out));
  int code = ABI.plan_open((const uint8_t *)RAW(header), (size_t)XLENGTH(header), &out);
  die_code("plan_open", code);
  return ranges_to_matrix(&out);
}

SEXP C_plan_eic(SEXP ptr, SEXP target, SEXP from, SEXP to, SEXP ppm, SEXP mz_tol)
{
  REQUIRE_BOUND(ABI.plan_eic, "plan_eic");
  MzML *handle = (MzML *)R_ExternalPtrAddr(ptr);
  if (handle == NULL)
    error("quantion: the file handle is not valid");

  Buf out;
  memset(&out, 0, sizeof(out));
  int code = ABI.plan_eic(handle, asReal(target), asReal(from), asReal(to),
                          asReal(ppm), asReal(mz_tol), &out);
  die_code("plan_eic", code);
  return ranges_to_matrix(&out);
}

SEXP C_store_new(void)
{
  range_store *store = store_new();
  if (store == NULL)
    error("quantion: cannot make a range store");
  SEXP ptr = PROTECT(R_MakeExternalPtr(store, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(ptr, finalize_store, TRUE);
  UNPROTECT(1);
  return ptr;
}

SEXP C_store_add(SEXP store_ptr, SEXP offset, SEXP bytes)
{
  range_store *store = (range_store *)R_ExternalPtrAddr(store_ptr);
  if (store == NULL)
    error("quantion: the range store is not valid");
  if (TYPEOF(bytes) != RAWSXP)
    error("quantion: bytes must be a raw vector");
  if (store_add(store, (uint64_t)asReal(offset), (const uint8_t *)RAW(bytes),
                (uint64_t)XLENGTH(bytes)) != 0)
    error("quantion: cannot hold that range");
  return R_NilValue;
}

SEXP C_parse_ion_source(SEXP store_ptr, SEXP max_cache_size)
{
  REQUIRE_BOUND(ABI.parse_ion_source, "parse_ion_source");
  range_store *store = (range_store *)R_ExternalPtrAddr(store_ptr);
  if (store == NULL)
    error("quantion: the range store is not valid");

  size_t cache = (max_cache_size == R_NilValue) ? 0 : (size_t)asReal(max_cache_size);
  MzML *handle = NULL;
  int code = ABI.parse_ion_source(store_serve, store, cache, &handle);
  die_code("parse_ion_source", code);

  SEXP ptr = PROTECT(R_MakeExternalPtr(handle, store_ptr, R_NilValue));
  R_RegisterCFinalizerEx(ptr, finalize_mzml, TRUE);
  UNPROTECT(1);
  return ptr;
}

SEXP C_parse_ion(SEXP bin, SEXP max_cache_size)
{
  if (TYPEOF(bin) != RAWSXP)
    error("quantion: data must be a raw vector");
  REQUIRE_BOUND(ABI.parse_bin, "parse_bin");

  size_t cache = (max_cache_size == R_NilValue) ? 0 : (size_t)asReal(max_cache_size);

  MzML *handle = NULL;
  int code = ABI.parse_bin(
      (const unsigned char *)RAW(bin),
      (size_t)XLENGTH(bin),
      cache,
      &handle);
  die_code("parse_bin", code);

  SEXP ptr = PROTECT(R_MakeExternalPtr(handle, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(ptr, finalize_mzml, TRUE);
  UNPROTECT(1);
  return ptr;
}


SEXP C_parse_ion_path(SEXP path, SEXP max_cache_size)
{
  if (TYPEOF(path) != STRSXP || LENGTH(path) != 1)
    error("quantion: path must be a single string");

  const char *value = CHAR(STRING_ELT(path, 0));
  size_t cache = (max_cache_size == R_NilValue) ? 0 : (size_t)asReal(max_cache_size);

  REQUIRE_BOUND(ABI.parse_ion_path, "parse_ion_path");

  MzML *handle = NULL;
  int code = ABI.parse_ion_path(value, cache, &handle);
  die_code("parse_ion_path", code);

  SEXP ptr = PROTECT(R_MakeExternalPtr(handle, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(ptr, finalize_mzml, TRUE);
  UNPROTECT(1);
  return ptr;
}

SEXP C_get_scans(SEXP bin, SEXP query_type, SEXP a, SEXP b, SEXP level)
{
  MzML *handle = GetHandle(bin);
  REQUIRE_BOUND(ABI.get_scans, "get_scans");
  REQUIRE_BOUND(ABI.free_, "free_");
  Buf out = (Buf){0};
  int code = ABI.get_scans(handle, (uint8_t)asInteger(query_type), asReal(a), asReal(b), (uint8_t)asInteger(level), &out);
  die_code("get_scans", code);
  SEXP res = mk_string_len(out.ptr, out.len);
  return res;
}

SEXP C_get_features(SEXP dir_path, SEXP from_time, SEXP to_time,
                    SEXP eic_ppm_tol, SEXP eic_mz_tol,
                    SEXP grid_start, SEXP grid_end, SEXP grid_step,
                    SEXP group_ppm_tol, SEXP group_mz_tol, SEXP group_rt_tol,
                    SEXP frequency, SEXP options, SEXP cores)
{
  if (TYPEOF(dir_path) != STRSXP || LENGTH(dir_path) != 1)
    error("dir_path must be a length-1 character string");
  REQUIRE_BOUND(ABI.get_features, "get_features");
  REQUIRE_BOUND(ABI.free_, "free_");

  const char *path = CHAR(STRING_ELT(dir_path, 0));

  int ncores = asInteger(cores);
  if (ncores < 1)
    ncores = 1;

  CPeakOptions opts;
  const CPeakOptions *opt_ptr = NULL;
  as_opts_ptr(options, &opts, &opt_ptr);

  Buf out = (Buf){0};
  int code = ABI.get_features(
      path,
      asReal(from_time), asReal(to_time),
      asReal(eic_ppm_tol), asReal(eic_mz_tol),
      asReal(grid_start), asReal(grid_end), asReal(grid_step),
      asReal(group_ppm_tol), asReal(group_mz_tol), asReal(group_rt_tol),
      asInteger(frequency),
      opt_ptr, (int32_t)ncores,
      &out);

  die_code("get_features", code);
  SEXP res = mk_string_len(out.ptr, out.len);
  return res;
}
