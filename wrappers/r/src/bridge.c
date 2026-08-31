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
typedef int32_t (*fn_plan_scans)(MzML *, uint8_t, double, double, uint8_t, Buf *);
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
  fn_plan_scans plan_scans;
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
  ABI.plan_scans = (fn_plan_scans)DLSYM(abi_handle, "plan_scans");
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


typedef struct
{
  uint32_t element_type;
  uint64_t byte_offset;
  uint64_t byte_length;
} bridge_section;

typedef struct
{
  const unsigned char *raw;
  uint64_t total_bytes;
  uint16_t payload_kind;
  uint64_t record_count;
  uint32_t section_count;
  bridge_section sections[64];
  uint32_t section_ids[64];
} bridge_view;

static uint16_t read_u16(const unsigned char *raw, uint64_t at)
{
  uint16_t value;
  memcpy(&value, raw + at, sizeof(value));
  return value;
}

static uint32_t read_u32(const unsigned char *raw, uint64_t at)
{
  uint32_t value;
  memcpy(&value, raw + at, sizeof(value));
  return value;
}

static uint64_t read_u64(const unsigned char *raw, uint64_t at)
{
  uint64_t value;
  memcpy(&value, raw + at, sizeof(value));
  return value;
}

static double read_f64(const unsigned char *raw, uint64_t at)
{
  double value;
  memcpy(&value, raw + at, sizeof(value));
  return value;
}

static size_t size_of_element(uint32_t element_type)
{
  if (element_type == QUANTION_ELEMENT_F64)
    return 8;
  if (element_type == QUANTION_ELEMENT_U32)
    return 4;
  if (element_type == QUANTION_ELEMENT_U64)
    return 8;
  if (element_type == QUANTION_ELEMENT_U8)
    return 1;
  return 0;
}

static uint32_t expected_element_type(uint32_t section_id)
{
  if (section_id == QUANTION_SECTION_POINT_STARTS)
    return QUANTION_ELEMENT_U64;
  if (section_id == QUANTION_SECTION_IMAGE_SHAPE || section_id == QUANTION_SECTION_IMAGE_COUNTS)
    return QUANTION_ELEMENT_U32;
  if (section_id <= QUANTION_SECTION_IMAGE_DATA)
    return QUANTION_ELEMENT_F64;
  return 0;
}

static void read_bridge_view(bridge_view *bridge, const unsigned char *raw, size_t len)
{
  if (len < QUANTION_BRIDGE_HEADER_BYTES)
    error("quantion bridge: buffer is shorter than the header");

  uint64_t total = (uint64_t)len;
  if (read_u32(raw, 0) != QUANTION_BRIDGE_MAGIC)
    error("quantion bridge: magic does not match");
  if (read_u16(raw, 4) != QUANTION_BRIDGE_LAYOUT_VERSION)
    error("quantion bridge: layout version is not supported");

  uint32_t count = read_u32(raw, 8);
  if (read_u32(raw, 12) != QUANTION_BRIDGE_HEADER_BYTES)
    error("quantion bridge: section table does not start after the header");
  if (read_u64(raw, 16) != total)
    error("quantion bridge: total bytes does not match the buffer");
  if (count > 64)
    error("quantion bridge: too many sections");
  if ((uint64_t)count * QUANTION_SECTION_ENTRY_BYTES > total - QUANTION_BRIDGE_HEADER_BYTES)
    error("quantion bridge: section table does not fit");

  bridge->raw = raw;
  bridge->total_bytes = total;
  bridge->payload_kind = read_u16(raw, 6);
  bridge->record_count = read_u64(raw, 24);
  bridge->section_count = count;

  uint64_t reach = QUANTION_BRIDGE_HEADER_BYTES + (uint64_t)count * QUANTION_SECTION_ENTRY_BYTES;
  for (uint32_t index = 0; index < count; index++)
  {
    uint64_t start = QUANTION_BRIDGE_HEADER_BYTES + (uint64_t)index * QUANTION_SECTION_ENTRY_BYTES;
    uint32_t section_id = read_u32(raw, start);
    uint32_t element_type = read_u32(raw, start + 4);
    uint64_t offset = read_u64(raw, start + 8);
    uint64_t length = read_u64(raw, start + 16);

    for (uint32_t seen = 0; seen < index; seen++)
      if (bridge->section_ids[seen] == section_id)
        error("quantion bridge: section %u appears twice", section_id);
    if (offset > total || length > total - offset)
      error("quantion bridge: section %u runs past the end", section_id);
    if (offset % 8 != 0)
      error("quantion bridge: section %u is not aligned to eight bytes", section_id);
    if (offset < reach)
      error("quantion bridge: section %u overlaps the section before it", section_id);

    size_t element_size = size_of_element(element_type);
    if (element_size == 0)
      error("quantion bridge: section %u has an unknown element type", section_id);
    if (length % (uint64_t)element_size != 0)
      error("quantion bridge: section %u is not a whole number of elements", section_id);

    uint32_t expected = expected_element_type(section_id);
    if (expected != 0 && expected != element_type)
      error("quantion bridge: section %u has the wrong element type", section_id);

    reach = offset + length;
    bridge->section_ids[index] = section_id;
    bridge->sections[index].element_type = element_type;
    bridge->sections[index].byte_offset = offset;
    bridge->sections[index].byte_length = length;
  }
}

static const bridge_section *find_section(const bridge_view *bridge, uint32_t section_id)
{
  for (uint32_t index = 0; index < bridge->section_count; index++)
    if (bridge->section_ids[index] == section_id)
      return &bridge->sections[index];
  error("quantion bridge: section %u is missing", section_id);
  return NULL;
}

static SEXP number_column(const bridge_view *bridge, uint32_t section_id, R_xlen_t count)
{
  const bridge_section *section = find_section(bridge, section_id);
  if (section->byte_length / 8 != (uint64_t)count)
    error("quantion bridge: section %u does not match the record count", section_id);
  SEXP column = PROTECT(Rf_allocVector(REALSXP, count));
  if (count > 0)
    memcpy(REAL(column), bridge->raw + section->byte_offset, (size_t)count * 8);
  UNPROTECT(1);
  return column;
}

static SEXP point_column(const bridge_view *bridge, uint32_t section_id,
                         const unsigned char *starts_at, R_xlen_t count)
{
  const bridge_section *values = find_section(bridge, section_id);
  SEXP column = PROTECT(Rf_allocVector(VECSXP, count));
  for (R_xlen_t index = 0; index < count; index++)
  {
    uint64_t from = read_u64(starts_at, (uint64_t)index * 8);
    uint64_t to = read_u64(starts_at, (uint64_t)(index + 1) * 8);
    if (to < from || to > values->byte_length / 8)
      error("quantion bridge: point starts run past section %u", section_id);
    R_xlen_t points = (R_xlen_t)(to - from);
    SEXP scan = PROTECT(Rf_allocVector(REALSXP, points));
    if (points > 0)
      memcpy(REAL(scan), bridge->raw + values->byte_offset + from * 8, (size_t)points * 8);
    SET_VECTOR_ELT(column, index, scan);
    UNPROTECT(1);
  }
  UNPROTECT(1);
  return column;
}

static const uint32_t METADATA_SECTIONS[10] = {
    QUANTION_SECTION_RT_SECONDS,
    QUANTION_SECTION_BASE_PEAK_MZ,
    QUANTION_SECTION_SELECTED_ION_MZ,
    QUANTION_SECTION_BASE_PEAK_INT,
    QUANTION_SECTION_TOTAL_ION_CURRENT,
    QUANTION_SECTION_MS_LEVEL,
    QUANTION_SECTION_POLARITY,
    QUANTION_SECTION_POSITION_X,
    QUANTION_SECTION_POSITION_Y,
    QUANTION_SECTION_POSITION_Z,
};

static const char *METADATA_NAMES[10] = {
    "rt_seconds", "base_peak_mz", "selected_ion_mz", "base_peak_int",
    "total_ion_current", "ms_level", "polarity",
    "position_x", "position_y", "position_z"};

static SEXP set_data_frame(SEXP columns, SEXP names, R_xlen_t rows)
{
  Rf_setAttrib(columns, R_NamesSymbol, names);
  SEXP row_names = PROTECT(Rf_allocVector(INTSXP, 2));
  INTEGER(row_names)[0] = NA_INTEGER;
  INTEGER(row_names)[1] = -(int)rows;
  Rf_setAttrib(columns, R_RowNamesSymbol, row_names);
  UNPROTECT(1);
  Rf_setAttrib(columns, R_ClassSymbol, Rf_mkString("data.frame"));
  return columns;
}

static SEXP build_metadata(const bridge_view *bridge, R_xlen_t count)
{
  SEXP columns = PROTECT(Rf_allocVector(VECSXP, 10));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 10));
  for (int index = 0; index < 10; index++)
  {
    SET_VECTOR_ELT(columns, index, number_column(bridge, METADATA_SECTIONS[index], count));
    SET_STRING_ELT(names, index, Rf_mkChar(METADATA_NAMES[index]));
  }
  set_data_frame(columns, names, count);
  UNPROTECT(2);
  return columns;
}

struct scans_bridge_data
{
  const unsigned char *ptr;
  size_t len;
};

static SEXP scans_bridge_build(void *data)
{
  struct scans_bridge_data *held = (struct scans_bridge_data *)data;
  bridge_view bridge;
  read_bridge_view(&bridge, held->ptr, held->len);
  if (bridge.payload_kind != QUANTION_PAYLOAD_SCANS)
    error("quantion bridge: expected a scans payload");

  R_xlen_t count = (R_xlen_t)bridge.record_count;
  const bridge_section *starts = find_section(&bridge, QUANTION_SECTION_POINT_STARTS);
  if (starts->byte_length / 8 != (uint64_t)count + 1)
    error("quantion bridge: point starts do not match the scan count");

  const bridge_section *mz = find_section(&bridge, QUANTION_SECTION_MZ);
  const bridge_section *intensity = find_section(&bridge, QUANTION_SECTION_INTENSITY);
  if (intensity->byte_length != mz->byte_length)
    error("quantion bridge: intensity length does not match m/z length");

  const unsigned char *starts_at = bridge.raw + starts->byte_offset;
  uint64_t previous = 0;
  for (R_xlen_t index = 0; index <= count; index++)
  {
    uint64_t value = read_u64(starts_at, (uint64_t)index * 8);
    if (value < previous)
      error("quantion bridge: point starts do not increase");
    previous = value;
  }
  if (previous != mz->byte_length / 8)
    error("quantion bridge: point starts do not span the m/z section");

  SEXP columns = PROTECT(Rf_allocVector(VECSXP, 4));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
  SET_VECTOR_ELT(columns, 0, number_column(&bridge, QUANTION_SECTION_RT, count));
  SET_STRING_ELT(names, 0, Rf_mkChar("rt"));
  SET_VECTOR_ELT(columns, 1, point_column(&bridge, QUANTION_SECTION_MZ, starts_at, count));
  SET_STRING_ELT(names, 1, Rf_mkChar("mz"));
  SET_VECTOR_ELT(columns, 2, point_column(&bridge, QUANTION_SECTION_INTENSITY, starts_at, count));
  SET_STRING_ELT(names, 2, Rf_mkChar("intensity"));
  SET_VECTOR_ELT(columns, 3, build_metadata(&bridge, count));
  SET_STRING_ELT(names, 3, Rf_mkChar("metadata"));
  set_data_frame(columns, names, count);
  UNPROTECT(2);
  return columns;
}

static void scans_bridge_free(void *data, Rboolean jump)
{
  (void)jump;
  struct scans_bridge_data *held = (struct scans_bridge_data *)data;
  if (held->ptr && ABI.free_)
    ABI.free_((unsigned char *)held->ptr, held->len);
  held->ptr = NULL;
}


struct buf_copy_data
{
  const unsigned char *first_ptr;
  size_t first_len;
  const unsigned char *second_ptr;
  size_t second_len;
  int as_raw;
};

static SEXP copy_one_buf(const unsigned char *ptr, size_t len, int as_raw)
{
  if (as_raw)
  {
    SEXP out = Rf_allocVector(RAWSXP, (R_xlen_t)len);
    if (len > 0)
      memcpy(RAW(out), ptr, len);
    return out;
  }
  R_xlen_t count = (R_xlen_t)(len / 8);
  SEXP out = Rf_allocVector(REALSXP, count);
  if (count > 0)
    memcpy(REAL(out), ptr, (size_t)count * 8);
  return out;
}

static SEXP buf_copy_build(void *data)
{
  struct buf_copy_data *held = (struct buf_copy_data *)data;
  if (held->second_ptr == NULL)
    return copy_one_buf(held->first_ptr, held->first_len, held->as_raw);

  SEXP x = PROTECT(copy_one_buf(held->first_ptr, held->first_len, held->as_raw));
  SEXP y = PROTECT(copy_one_buf(held->second_ptr, held->second_len, held->as_raw));
  SEXP pair = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(pair, 0, x);
  SET_VECTOR_ELT(pair, 1, y);
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, Rf_mkChar("x"));
  SET_STRING_ELT(names, 1, Rf_mkChar("y"));
  Rf_setAttrib(pair, R_NamesSymbol, names);
  UNPROTECT(4);
  return pair;
}

static void buf_copy_free(void *data, Rboolean jump)
{
  (void)jump;
  struct buf_copy_data *held = (struct buf_copy_data *)data;
  if (held->first_ptr && ABI.free_)
    ABI.free_((unsigned char *)held->first_ptr, held->first_len);
  if (held->second_ptr && ABI.free_)
    ABI.free_((unsigned char *)held->second_ptr, held->second_len);
  held->first_ptr = NULL;
  held->second_ptr = NULL;
}

static SEXP ranges_build(void *data)
{
  struct buf_copy_data *held = (struct buf_copy_data *)data;
  size_t count = held->first_len / 16;
  SEXP result = PROTECT(allocMatrix(REALSXP, (int)count, 2));
  double *values = REAL(result);
  for (size_t index = 0; index < count; index++)
  {
    uint64_t offset = 0;
    uint64_t length = 0;
    memcpy(&offset, held->first_ptr + index * 16, 8);
    memcpy(&length, held->first_ptr + index * 16 + 8, 8);
    values[index] = (double)offset;
    values[count + index] = (double)length;
  }
  UNPROTECT(1);
  return result;
}


typedef struct
{
  const char *name;
  uint32_t id;
  uint32_t second_id;
  int kind;
} record_column;

#define COLUMN_NUMBER 0
#define COLUMN_COUNT 1
#define COLUMN_TEXT 2

static const record_column PEAK_COLUMNS[] = {
    {"from", 18, 0, COLUMN_NUMBER}, {"to", 19, 0, COLUMN_NUMBER},
    {"rt", 20, 0, COLUMN_NUMBER}, {"integral", 21, 0, COLUMN_NUMBER},
    {"intensity", 22, 0, COLUMN_NUMBER}, {"n_points", 23, 0, COLUMN_COUNT},
    {"noise", 24, 0, COLUMN_NUMBER}, {"r2", 25, 0, COLUMN_NUMBER}};

static const record_column FEATURE_COLUMNS[] = {
    {"mz", 26, 0, COLUMN_NUMBER}, {"rt", 27, 0, COLUMN_NUMBER},
    {"from", 28, 0, COLUMN_NUMBER}, {"to", 29, 0, COLUMN_NUMBER},
    {"intensity", 30, 0, COLUMN_NUMBER}, {"integral", 31, 0, COLUMN_NUMBER},
    {"n_points", 32, 0, COLUMN_COUNT}, {"noise", 33, 0, COLUMN_NUMBER}};

static const record_column CHROM_COLUMNS[] = {
    {"index", 34, 0, COLUMN_COUNT}, {"id", 42, 43, COLUMN_TEXT},
    {"ort", 35, 0, COLUMN_NUMBER}, {"rt", 36, 0, COLUMN_NUMBER},
    {"from", 37, 0, COLUMN_NUMBER}, {"to", 38, 0, COLUMN_NUMBER},
    {"intensity", 39, 0, COLUMN_NUMBER}, {"integral", 40, 0, COLUMN_NUMBER},
    {"total_area", 41, 0, COLUMN_NUMBER}, {"timestamp", 44, 45, COLUMN_TEXT}};

static const record_column FIT_COLUMNS[] = {
    {"shape", 46, 0, COLUMN_NUMBER}, {"height", 47, 0, COLUMN_NUMBER},
    {"center", 48, 0, COLUMN_NUMBER}, {"fwhm", 49, 0, COLUMN_NUMBER},
    {"tail", 50, 0, COLUMN_NUMBER}, {"r2", 51, 0, COLUMN_NUMBER}};

static const record_column EIC_PEAK_COLUMNS[] = {
    {"id", 52, 53, COLUMN_TEXT}, {"mz", 54, 0, COLUMN_NUMBER},
    {"ort", 55, 0, COLUMN_NUMBER}, {"rt", 56, 0, COLUMN_NUMBER},
    {"from", 57, 0, COLUMN_NUMBER}, {"to", 58, 0, COLUMN_NUMBER},
    {"intensity", 59, 0, COLUMN_NUMBER}, {"integral", 60, 0, COLUMN_NUMBER},
    {"noise", 61, 0, COLUMN_NUMBER}};

static const record_column CONSENSUS_COLUMNS[] = {
    {"mz", 62, 0, COLUMN_NUMBER}, {"rt", 63, 0, COLUMN_NUMBER},
    {"from", 64, 0, COLUMN_NUMBER}, {"to", 65, 0, COLUMN_NUMBER},
    {"intensity", 66, 0, COLUMN_NUMBER}, {"integral", 67, 0, COLUMN_NUMBER},
    {"frequency", 68, 0, COLUMN_NUMBER}};

static const record_column FOUND_COLUMNS[] = {
    {"id", 69, 70, COLUMN_TEXT}, {"mz", 71, 0, COLUMN_NUMBER},
    {"rt", 72, 0, COLUMN_NUMBER}, {"from", 73, 0, COLUMN_NUMBER},
    {"to", 74, 0, COLUMN_NUMBER}, {"intensity", 75, 0, COLUMN_NUMBER},
    {"integral", 76, 0, COLUMN_NUMBER}, {"n_points", 77, 0, COLUMN_COUNT},
    {"noise", 78, 0, COLUMN_NUMBER}};

static const record_column *columns_for_kind(uint16_t kind, int *count)
{
  switch (kind)
  {
  case QUANTION_PAYLOAD_PEAKS:
    *count = 8;
    return PEAK_COLUMNS;
  case QUANTION_PAYLOAD_FEATURES:
    *count = 8;
    return FEATURE_COLUMNS;
  case QUANTION_PAYLOAD_CHROM_PEAKS:
    *count = 10;
    return CHROM_COLUMNS;
  case QUANTION_PAYLOAD_FIT_RESULT:
    *count = 6;
    return FIT_COLUMNS;
  case QUANTION_PAYLOAD_EIC_PEAKS:
    *count = 9;
    return EIC_PEAK_COLUMNS;
  case QUANTION_PAYLOAD_CONSENSUS_FEATURES:
    *count = 7;
    return CONSENSUS_COLUMNS;
  case QUANTION_PAYLOAD_FOUND_FEATURES:
    *count = 9;
    return FOUND_COLUMNS;
  default:
    *count = 0;
    return NULL;
  }
}

static SEXP count_column(const bridge_view *bridge, uint32_t section_id, R_xlen_t rows)
{
  const bridge_section *section = find_section(bridge, section_id);
  if (section->byte_length / 4 != (uint64_t)rows)
    error("quantion bridge: section %u does not match the record count", section_id);
  SEXP column = PROTECT(Rf_allocVector(REALSXP, rows));
  for (R_xlen_t index = 0; index < rows; index++)
    REAL(column)[index] = (double)read_u32(bridge->raw, section->byte_offset + (uint64_t)index * 4);
  UNPROTECT(1);
  return column;
}

static SEXP text_column(const bridge_view *bridge, uint32_t starts_id, uint32_t bytes_id,
                        R_xlen_t rows)
{
  const bridge_section *starts = find_section(bridge, starts_id);
  const bridge_section *blob = find_section(bridge, bytes_id);
  if (starts->byte_length / 8 != (uint64_t)rows + 1)
    error("quantion bridge: text starts do not match the record count");

  SEXP column = PROTECT(Rf_allocVector(STRSXP, rows));
  uint64_t previous = 0;
  for (R_xlen_t index = 0; index < rows; index++)
  {
    uint64_t from = read_u64(bridge->raw, starts->byte_offset + (uint64_t)index * 8);
    uint64_t to = read_u64(bridge->raw, starts->byte_offset + (uint64_t)(index + 1) * 8);
    if (from < previous || to < from || to > blob->byte_length)
      error("quantion bridge: text starts run past the blob");
    previous = from;
    SET_STRING_ELT(column, index,
                   Rf_mkCharLenCE((const char *)(bridge->raw + blob->byte_offset + from),
                                  (int)(to - from), CE_UTF8));
  }
  UNPROTECT(1);
  return column;
}

static SEXP records_bridge_build(void *data)
{
  struct scans_bridge_data *held = (struct scans_bridge_data *)data;
  bridge_view bridge;
  read_bridge_view(&bridge, held->ptr, held->len);

  int column_count = 0;
  const record_column *spec = columns_for_kind(bridge.payload_kind, &column_count);
  if (spec == NULL)
    error("quantion bridge: payload kind %u is not a record table", bridge.payload_kind);

  R_xlen_t rows = (R_xlen_t)bridge.record_count;
  SEXP columns = PROTECT(Rf_allocVector(VECSXP, column_count));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, column_count));
  for (int index = 0; index < column_count; index++)
  {
    SEXP column;
    if (spec[index].kind == COLUMN_NUMBER)
      column = number_column(&bridge, spec[index].id, rows);
    else if (spec[index].kind == COLUMN_COUNT)
      column = count_column(&bridge, spec[index].id, rows);
    else
      column = text_column(&bridge, spec[index].id, spec[index].second_id, rows);
    SET_VECTOR_ELT(columns, index, column);
    SET_STRING_ELT(names, index, Rf_mkChar(spec[index].name));
  }
  set_data_frame(columns, names, rows);
  UNPROTECT(2);
  return columns;
}

static SEXP take_records_bridge(const char *label, Buf *out, SEXP cont)
{
  (void)label;
  struct scans_bridge_data held = {out->ptr, out->len};
  return R_UnwindProtect(records_bridge_build, &held, scans_bridge_free, &held, cont);
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
  SEXP cont = PROTECT(R_MakeUnwindCont());
  Buf out = (Buf){0};
  int32_t code = ABI.mzml_to_bin(
      handle,
      &out,
      (uint8_t)lv,
      (uint8_t)(fc ? 1 : 0));
  if (code != 0)
  {
    UNPROTECT(1);
    die_code("mzml_to_bin", code);
  }
  struct buf_copy_data held = {out.ptr, out.len, NULL, 0, 1};
  SEXP res = R_UnwindProtect(buf_copy_build, &held, buf_copy_free, &held, cont);
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
  SEXP cont = PROTECT(R_MakeUnwindCont());
  Buf out = (Buf){0};
  int code = ABI.get_peak(REAL(x), REAL(y), (size_t)n, asReal(rt), asReal(range), opt_ptr, &out);
  if (code != 0)
  {
    UNPROTECT(1);
    die_code("get_peak", code);
  }
  SEXP res = take_records_bridge("get_peak", &out, cont);
  UNPROTECT(1);
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
  SEXP cont = PROTECT(R_MakeUnwindCont());
  Buf out = (Buf){0};
  int code = ABI.fit_peak(REAL(x), REAL(y), (size_t)n,
                          asReal(rt), asReal(intensity), (int32_t)asInteger(shape), &out);
  if (code != 0)
  {
    UNPROTECT(1);
    die_code("fit_peak", code);
  }
  SEXP res = take_records_bridge("fit_peak", &out, cont);
  UNPROTECT(1);
  return res;
}

SEXP C_draw_peak(SEXP x, SEXP shape, SEXP height, SEXP center, SEXP fwhm, SEXP tail)
{
  if (TYPEOF(x) != REALSXP)
    error("numeric");
  REQUIRE_BOUND(ABI.draw_peak, "draw_peak");
  REQUIRE_BOUND(ABI.free_, "free_");
  R_xlen_t n = XLENGTH(x);
  SEXP cont = PROTECT(R_MakeUnwindCont());
  Buf out = (Buf){0};
  int code = ABI.draw_peak(REAL(x), (size_t)n, (int32_t)asInteger(shape),
                           asReal(height), asReal(center), asReal(fwhm), asReal(tail), &out);
  if (code != 0)
  {
    UNPROTECT(1);
    die_code("draw_peak", code);
  }
  struct buf_copy_data held = {out.ptr, out.len, NULL, 0, 0};
  SEXP Ry = R_UnwindProtect(buf_copy_build, &held, buf_copy_free, &held, cont);
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

  SEXP cont = PROTECT(R_MakeUnwindCont());
  Buf out = (Buf){0};
  int code = ABI.get_peaks_from_eic(
      handle, REAL(rts), REAL(mzs), REAL(ranges),
      offs, lens, ids_buf, ids_len,
      (size_t)n, asReal(from_left), asReal(to_right), opt_ptr, ncores, &out);

  if (code != 0)
  {
    UNPROTECT(1);
    die_code("get_peaks_from_eic", code);
  }
  SEXP res = take_records_bridge("get_peaks_from_eic", &out, cont);
  UNPROTECT(1);
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

  SEXP cont = PROTECT(R_MakeUnwindCont());
  Buf out = (Buf){0};
  int code = ABI.get_peaks_from_chrom(
      handle,
      uidx, REAL(rts), REAL(ranges), (size_t)n, opt_ptr, ncores, &out);

  if (code != 0)
  {
    UNPROTECT(1);
    die_code("get_peaks_from_chrom", code);
  }
  SEXP res = take_records_bridge("get_peaks_from_chrom", &out, cont);
  UNPROTECT(1);
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
  SEXP cont = PROTECT(R_MakeUnwindCont());
  Buf bx = (Buf){0}, by = (Buf){0};
  int code = ABI.calculate_eic(
      handle,
      t, asReal(from), asReal(to), asReal(ppm_tol), asReal(mz_tol),
      &bx, &by);
  if (code != 0)
  {
    UNPROTECT(1);
    die_code("calculate_eic", code);
  }
  struct buf_copy_data held = {bx.ptr, bx.len, by.ptr, by.len, 0};
  SEXP out = R_UnwindProtect(buf_copy_build, &held, buf_copy_free, &held, cont);
  UNPROTECT(1);
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
  SEXP cont = PROTECT(R_MakeUnwindCont());
  Buf out = (Buf){0};
  int code = ABI.find_peaks(REAL(x), REAL(y), (size_t)n, opt_ptr, &out);
  if (code != 0)
  {
    UNPROTECT(1);
    die_code("find_peaks", code);
  }
  SEXP res = take_records_bridge("find_peaks", &out, cont);
  UNPROTECT(1);
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
  SEXP cont = PROTECT(R_MakeUnwindCont());
  Buf out = (Buf){0};
  int code = ABI.calculate_baseline(REAL(y), (size_t)n, (int32_t)lam, (int32_t)maxit, &out);
  if (code != 0)
  {
    UNPROTECT(1);
    die_code("calculate_baseline", code);
  }
  struct buf_copy_data held = {out.ptr, out.len, NULL, 0, 0};
  SEXP Ry = R_UnwindProtect(buf_copy_build, &held, buf_copy_free, &held, cont);
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

  SEXP cont = PROTECT(R_MakeUnwindCont());
  Buf out = (Buf){0};
  int code = ABI.find_features(
      handle,
      asReal(from_time), asReal(to_time),
      asReal(eic_ppm_tol), asReal(eic_mz_tol),
      asReal(grid_start), asReal(grid_end),
      asReal(grid_step_ppm),
      opt_ptr, (int32_t)ncores, &out);

  if (code != 0)
  {
    UNPROTECT(1);
    die_code("find_features", code);
  }
  SEXP res = take_records_bridge("find_features", &out, cont);
  UNPROTECT(1);
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

  SEXP cont = PROTECT(R_MakeUnwindCont());
  Buf out = (Buf){0};
  int32_t code = ABI.find_feature(
      handle, REAL(rts), REAL(mzs), REAL(wins),
      offs, lens, ids_buf, ids_len,
      (size_t)n, ncores,
      asReal(scan_ppm), asReal(scan_mz), asReal(eic_ppm), asReal(eic_mz),
      opt_ptr, &out);

  if (code != 0)
  {
    UNPROTECT(1);
    die_code("find_feature", code);
  }
  SEXP res = take_records_bridge("find_feature", &out, cont);
  UNPROTECT(1);
  return res;
}

static SEXP ranges_to_matrix(Buf *out)
{
  SEXP cont = PROTECT(R_MakeUnwindCont());
  struct buf_copy_data held = {out->ptr, out->len, NULL, 0, 0};
  SEXP result = R_UnwindProtect(ranges_build, &held, buf_copy_free, &held, cont);
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

SEXP C_plan_scans(SEXP ptr, SEXP query_type, SEXP from_value, SEXP to_value, SEXP level)
{
  REQUIRE_BOUND(ABI.plan_scans, "plan_scans");
  MzML *handle = (MzML *)R_ExternalPtrAddr(ptr);
  if (handle == NULL)
    error("quantion: the file handle is not valid");

  Buf out;
  memset(&out, 0, sizeof(out));
  int code = ABI.plan_scans(handle, (uint8_t)asInteger(query_type),
                            asReal(from_value), asReal(to_value),
                            (uint8_t)asInteger(level), &out);
  die_code("plan_scans", code);
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

  SEXP cont = PROTECT(R_MakeUnwindCont());
  Buf out = (Buf){0};
  int code = ABI.get_scans(handle, (uint8_t)asInteger(query_type), asReal(a), asReal(b), (uint8_t)asInteger(level), &out);
  if (code != 0)
  {
    UNPROTECT(1);
    die_code("get_scans", code);
  }

  struct scans_bridge_data held = {out.ptr, out.len};
  SEXP result = R_UnwindProtect(scans_bridge_build, &held, scans_bridge_free, &held, cont);
  UNPROTECT(1);
  return result;
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

  SEXP cont = PROTECT(R_MakeUnwindCont());
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

  if (code != 0)
  {
    UNPROTECT(1);
    die_code("get_features", code);
  }
  SEXP res = take_records_bridge("get_features", &out, cont);
  UNPROTECT(1);
  return res;
}
