#if __has_include(<napi.h>)
#include <napi.h>
#else
#include "node-addon-api/napi.h"
#endif

#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <cmath>
#include <string>
#include <vector>
#include <atomic>
#include "include/msutils.h"

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

struct MzMLWrapper
{
  MzML *ptr;
  size_t estimated_size;
};

static_assert(sizeof(CPeakOptions) == 80, "CPeakOptions must be 80 bytes");

typedef uint32_t (*fn_msutils_abi_version)(void);
typedef size_t   (*fn_msutils_sizeof_peak_options)(void);
typedef int32_t (*fn_parse_mzml)(const unsigned char *, size_t, MzML **);
typedef int32_t (*fn_parse_bin)(const unsigned char *, size_t, size_t, MzML **);
typedef int32_t (*fn_parse_ion_path)(const char *, size_t, MzML **);
typedef int32_t (*fn_parse_ion_url)(const char *, size_t, MzML **);
typedef void (*fn_free_mzml)(MzML *);
typedef int32_t (*fn_bin_to_json)(const MzML *, Buf *);
typedef int32_t (*fn_bin_to_mzml)(const MzML *, Buf *);
typedef int32_t (*fn_get_peak)(const double *, const double *, size_t, double, double, const CPeakOptions *, Buf *);
typedef int32_t (*fn_calculate_eic)(const MzML *, double, double, double, double, double, Buf *, Buf *);
typedef double (*fn_find_noise_level)(const double *, size_t);
typedef int32_t (*fn_get_peaks_from_eic)(const MzML *, const double *, const double *, const double *, const uint32_t *, const uint32_t *, const unsigned char *, size_t, size_t, double, double, const CPeakOptions *, size_t, Buf *);
typedef int32_t (*fn_get_peaks_from_chrom)(const MzML *, const uint32_t *, const double *, const double *, size_t, const CPeakOptions *, size_t, Buf *);
typedef int32_t (*fn_find_peaks)(const double *, const double *, size_t, const CPeakOptions *, Buf *);
typedef int32_t (*fn_calculate_baseline)(const double *, size_t, int32_t, int32_t, Buf *);
typedef int32_t (*fn_find_features)(const MzML *, double, double, double, double, double, double, double, const CPeakOptions *, int32_t, Buf *);
typedef int32_t (*fn_find_feature)(const MzML *, const double *, const double *, const double *, const uint32_t *, const uint32_t *, const unsigned char *, size_t, size_t, size_t, double, double, double, double, const CPeakOptions *, Buf *);
typedef int32_t (*fn_mzml_to_bin)(const MzML *, Buf *, uint8_t, uint8_t);
typedef int32_t (*fn_get_features)(const char *, double, double, double, double, double, double, double, double, double, double, int32_t, const CPeakOptions *, int32_t, Buf *);
typedef int32_t (*fn_get_scans)(const MzML *, uint8_t, double, double, uint8_t, Buf *);
typedef int32_t (*fn_get_ion_image)(const MzML *, double, double, uint8_t, Buf *);
typedef void (*fn_free_)(unsigned char *, size_t);

typedef struct
{
  fn_msutils_abi_version msutils_abi_version;
  fn_msutils_sizeof_peak_options msutils_sizeof_peak_options;
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
  fn_find_feature find_feature;
  fn_mzml_to_bin mzml_to_bin;
  fn_parse_bin parse_bin;
  fn_parse_ion_path parse_ion_path;
  fn_parse_ion_url parse_ion_url;
  fn_free_ free_;
  fn_free_mzml free_mzml;
  fn_get_features get_features;
  fn_get_scans get_scans;
  fn_get_ion_image get_ion_image;
} msabi_t;

static msabi_t ABI{};
static DLIB LIB_HANDLE = NULL;

static std::atomic<size_t> OUTSTANDING_EXTERNAL{0};

struct OwnedBuf
{
  Buf buf{nullptr, 0};

  void Reset()
  {
    if (buf.ptr && ABI.free_)
      ABI.free_(buf.ptr, buf.len);
    buf.ptr = nullptr;
    buf.len = 0;
  }

  ~OwnedBuf()
  {
    Reset();
  }

  Buf *Out()
  {
    return &buf;
  }

  Buf Release()
  {
    Buf out = buf;
    buf.ptr = nullptr;
    buf.len = 0;
    return out;
  }
};

struct ExternalBuf
{
  unsigned char *ptr;
  size_t len;
};

static void FinalizeExternalBuffer(Napi::Env, uint8_t *, ExternalBuf *ext)
{
  if (ext)
  {
    if (ext->ptr && ABI.free_)
      ABI.free_(ext->ptr, ext->len);
    delete ext;
    OUTSTANDING_EXTERNAL.fetch_sub(1, std::memory_order_relaxed);
  }
}

static void FinalizeExternalArrayBuffer(Napi::Env, void *, ExternalBuf *ext)
{
  if (ext)
  {
    if (ext->ptr && ABI.free_)
      ABI.free_(ext->ptr, ext->len);
    delete ext;
    OUTSTANDING_EXTERNAL.fetch_sub(1, std::memory_order_relaxed);
  }
}

static void FinalizeMzML(Napi::Env env, MzMLWrapper *wrapper)
{
  if (wrapper)
  {
    if (wrapper->ptr && ABI.free_mzml)
    {
      ABI.free_mzml(wrapper->ptr);
      Napi::MemoryManagement::AdjustExternalMemory(env, -(int64_t)wrapper->estimated_size);
    }
    delete wrapper;
  }
}

static Napi::Value DisposeMzML(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  if (info.Length() < 1 || !info[0].IsExternal())
    return env.Undefined();

  MzMLWrapper *w = info[0].As<Napi::External<MzMLWrapper>>().Data();
  if (w && w->ptr)
  {
    if (ABI.free_mzml)
      ABI.free_mzml(w->ptr);
    Napi::MemoryManagement::AdjustExternalMemory(env, -(int64_t)w->estimated_size);
    w->ptr = nullptr;
  }
  return env.Undefined();
}

static MzML *GetHandle(Napi::Value val)
{
  if (!val.IsExternal())
    return nullptr;
  MzMLWrapper *w = val.As<Napi::External<MzMLWrapper>>().Data();
  if (!w || !w->ptr)
    return nullptr;
  return w->ptr;
}

static int resolve_required(void **fn, const char *name)
{
#if !defined(_WIN32)
  dlerror();
#endif
  *fn = DLSYM(LIB_HANDLE, name);
  if (!*fn)
  {
#if !defined(_WIN32)
    last_err = dlerror();
#endif
    return -1;
  }
  return 0;
}

static int abi_load(const char *path, const char **err)
{
  if (LIB_HANDLE)
  {
    DLCLOSE(LIB_HANDLE);
    LIB_HANDLE = NULL;
  }
  memset(&ABI, 0, sizeof(ABI));
#if !defined(_WIN32)
  dlerror();
#endif
  LIB_HANDLE = DLOPEN(path);
#if !defined(_WIN32)
  if (!LIB_HANDLE)
    last_err = dlerror();
#endif
  if (!LIB_HANDLE)
  {
    if (err)
      *err = last_err;
    return -1;
  }

  if (resolve_required((void **)&ABI.msutils_abi_version, "msutils_abi_version"))
    goto fail;
  if (ABI.msutils_abi_version() != MSUTILS_ABI_VERSION)
  {
    last_err = "ABI version mismatch: rebuild the native addon";
    goto fail;
  }
  if (resolve_required((void **)&ABI.msutils_sizeof_peak_options, "msutils_sizeof_peak_options"))
    goto fail;
  if (ABI.msutils_sizeof_peak_options() != sizeof(CPeakOptions))
  {
    last_err = "PeakOptions size mismatch: native binary and JS addon are out of sync — rebuild";
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
  if (resolve_required((void **)&ABI.calculate_eic, "calculate_eic"))
    goto fail;
  if (resolve_required((void **)&ABI.get_peaks_from_eic, "get_peaks_from_eic"))
    goto fail;
  if (resolve_required((void **)&ABI.get_peaks_from_chrom, "get_peaks_from_chrom"))
    goto fail;
  if (resolve_required((void **)&ABI.find_peaks, "find_peaks"))
    goto fail;

  ABI.calculate_baseline = (fn_calculate_baseline)DLSYM(LIB_HANDLE, "calculate_baseline");
  if (!ABI.calculate_baseline)
    ABI.calculate_baseline = (fn_calculate_baseline)DLSYM(LIB_HANDLE, "calculate_baseline_v2");

  if (resolve_required((void **)&ABI.find_features, "find_features"))
    goto fail;
  if (resolve_required((void **)&ABI.find_feature, "find_feature"))
    goto fail;
  if (resolve_required((void **)&ABI.mzml_to_bin, "mzml_to_bin"))
    goto fail;
  if (resolve_required((void **)&ABI.parse_bin, "parse_bin"))
    goto fail;
  ABI.parse_ion_path = (fn_parse_ion_path)DLSYM(LIB_HANDLE, "parse_ion_path");
  ABI.parse_ion_url = (fn_parse_ion_url)DLSYM(LIB_HANDLE, "parse_ion_url");
  if (resolve_required((void **)&ABI.get_features, "get_features"))
    goto fail;
  if (resolve_required((void **)&ABI.get_scans, "get_scans"))
    goto fail;

  ABI.get_ion_image = (fn_get_ion_image)DLSYM(LIB_HANDLE, "get_ion_image");

  ABI.find_noise_level = (fn_find_noise_level)DLSYM(LIB_HANDLE, "find_noise_level");
  ABI.free_ = (fn_free_)DLSYM(LIB_HANDLE, "free_");
  if (!ABI.free_)
    goto fail;
  if (resolve_required((void **)&ABI.free_mzml, "free_mzml"))
    goto fail;

  return 0;

fail:
  if (LIB_HANDLE)
  {
    DLCLOSE(LIB_HANDLE);
    LIB_HANDLE = NULL;
  }
  memset(&ABI, 0, sizeof(ABI));
  if (err)
    *err = last_err;
  return -1;
}

static const char *CodeMessage(int32_t code)
{
  if (code == 0)
    return "ok";
  if (code == 1)
    return "invalid arguments";
  if (code == 2)
    return "panic inside Rust";
  if (code == 4)
    return "parse error";
  if (code == 5)
    return "encode error";
  return "unknown";
}

static bool ThrowIfMissing(Napi::Env env, void *fn_ptr, const char *name)
{
  if (fn_ptr != nullptr)
    return true;
  std::string msg = "native symbol not exported: ";
  msg += name;
  Napi::Error::New(env, msg).ThrowAsJavaScriptException();
  return false;
}

static Napi::Value ThrowRc(Napi::Env env, const char *api, int32_t rc)
{
  std::string msg = api;
  msg += ": ";
  msg += CodeMessage(rc);
  Napi::Error::New(env, msg).ThrowAsJavaScriptException();
  return env.Undefined();
}

static Napi::Buffer<uint8_t> TakeBuffer(Napi::Env env, Buf *buf)
{
  ExternalBuf *ext = new ExternalBuf{buf->ptr, buf->len};
  OUTSTANDING_EXTERNAL.fetch_add(1, std::memory_order_relaxed);
  Napi::Buffer<uint8_t> out = Napi::Buffer<uint8_t>::New(env, (uint8_t *)ext->ptr, ext->len, FinalizeExternalBuffer, ext);
  buf->ptr = nullptr;
  buf->len = 0;
  return out;
}

static Napi::ArrayBuffer TakeArrayBuffer(Napi::Env env, Buf *buf)
{
  ExternalBuf *ext = new ExternalBuf{buf->ptr, buf->len};
  OUTSTANDING_EXTERNAL.fetch_add(1, std::memory_order_relaxed);
  Napi::ArrayBuffer ab = Napi::ArrayBuffer::New(env, (void *)ext->ptr, ext->len, FinalizeExternalArrayBuffer, ext);
  buf->ptr = nullptr;
  buf->len = 0;
  return ab;
}

static Napi::String TakeUtf8String(Napi::Env env, Buf *buf)
{
  Napi::String s = Napi::String::New(env, (const char *)buf->ptr, buf->len);
  if (ABI.free_)
    ABI.free_(buf->ptr, buf->len);
  buf->ptr = nullptr;
  buf->len = 0;
  return s;
}

static double GetF64(const Napi::Object &o, const char *key, double fallback)
{
  if (!o.Has(key))
    return fallback;
  Napi::Value v = o.Get(key);
  if (v.IsNumber())
    return v.As<Napi::Number>().DoubleValue();
  return fallback;
}

static int32_t GetI32(const Napi::Object &o, const char *key, int32_t fallback)
{
  if (!o.Has(key))
    return fallback;
  Napi::Value v = o.Get(key);
  if (v.IsNumber())
    return v.As<Napi::Number>().Int32Value();
  if (v.IsBoolean())
    return v.As<Napi::Boolean>().Value() ? 1 : 0;
  return fallback;
}

static const int32_t SHAPE_GAUSSIAN = 0;
static const int32_t SHAPE_EMG = 1;

static int32_t GetShapeCode(const Napi::Object &o, const char *key)
{
  if (!o.Has(key))
    return SHAPE_EMG;
  Napi::Value v = o.Get(key);
  if (v.IsNumber())
    return v.As<Napi::Number>().Int32Value();
  if (v.IsString())
  {
    std::string name = v.As<Napi::String>().Utf8Value();
    if (name == "gaussian")
      return SHAPE_GAUSSIAN;
    return SHAPE_EMG;
  }
  return SHAPE_EMG;
}

static const CPeakOptions *ReadPeakOptionsObject(Napi::Value value, CPeakOptions *out)
{
  if (value.IsUndefined() || value.IsNull())
    return nullptr;
  if (!value.IsObject() || value.IsArray() || value.IsBuffer())
    return nullptr;
  Napi::Object o = value.As<Napi::Object>();
  out->min_integral          = GetF64(o, "minIntegral",         NAN);
  out->min_intensity         = GetF64(o, "minIntensity",        NAN);
  out->min_peak_width_points = GetI32(o, "minPeakWidthPoints",  0);
  out->shape                 = GetShapeCode(o, "shape");
  out->noise                 = GetF64(o, "noise",               NAN);
  out->auto_noise            = GetI32(o, "autoNoise",           0);
  out->auto_baseline         = GetI32(o, "autoBaseline",        0);
  out->lambda                = GetI32(o, "lambda",              0);
  out->max_iterations        = GetI32(o, "maxIterations",       0);
  out->allow_overlap         = GetI32(o, "allowOverlap",        0);
  out->min_snr               = GetF64(o, "minSnr",              NAN);
  out->min_r2                = GetF64(o, "minR2",               NAN);
  out->kernel_size           = GetI32(o, "kernelSize",         0);
  return out;
}

static const double *Float64Ptr(const Napi::Float64Array &arr)
{
  return (const double *)((const uint8_t *)arr.ArrayBuffer().Data() + arr.ByteOffset());
}

static const uint32_t *Uint32Ptr(const Napi::Uint32Array &arr)
{
  return (const uint32_t *)((const uint8_t *)arr.ArrayBuffer().Data() + arr.ByteOffset());
}

struct PackedIds
{
  const uint32_t *offs_ptr = nullptr;
  const uint32_t *lens_ptr = nullptr;
  const unsigned char *ids_buf_ptr = nullptr;
  size_t ids_buf_len = 0;

  std::vector<uint32_t> offs;
  std::vector<uint32_t> lens;
  std::vector<unsigned char> ids_buf;
};

static bool BuildPackedIds(Napi::Env env, const Napi::Value &value, size_t count, PackedIds *out)
{
  (void)env;
  if (value.IsUndefined() || value.IsNull())
    return true;
  if (!value.IsArray())
    return false;

  Napi::Array ids = value.As<Napi::Array>();
  if ((size_t)ids.Length() != count)
    return false;

  out->offs.assign(count, 0);
  out->lens.assign(count, 0);

  std::vector<std::string> tmp;
  tmp.resize(count);

  size_t total = 0;
  for (size_t i = 0; i < count; i++)
  {
    Napi::Value v = ids.Get((uint32_t)i);
    if (v.IsString())
    {
      tmp[i] = v.As<Napi::String>().Utf8Value();
      total += tmp[i].size();
    }
    else
    {
      tmp[i].clear();
    }
  }

  out->ids_buf.resize(total);

  size_t cur = 0;
  for (size_t i = 0; i < count; i++)
  {
    const std::string &s = tmp[i];
    out->offs[i] = (uint32_t)cur;
    out->lens[i] = (uint32_t)s.size();
    if (!s.empty())
    {
      memcpy(out->ids_buf.data() + cur, s.data(), s.size());
      cur += s.size();
    }
  }

  out->offs_ptr = out->offs.data();
  out->lens_ptr = out->lens.data();
  out->ids_buf_ptr = out->ids_buf.data();
  out->ids_buf_len = out->ids_buf.size();
  return true;
}

static Napi::Value Bind(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  if (info.Length() < 1 || !info[0].IsString())
  {
    Napi::TypeError::New(env, "expected: path string").ThrowAsJavaScriptException();
    return env.Undefined();
  }

  if (OUTSTANDING_EXTERNAL.load(std::memory_order_relaxed) != 0)
  {
    Napi::Error::New(env, "cannot rebind native library while external buffers are still alive").ThrowAsJavaScriptException();
    return env.Undefined();
  }

  std::string path = info[0].As<Napi::String>();
  const char *err = nullptr;
  if (abi_load(path.c_str(), &err) != 0)
  {
    std::string msg = "dlopen failed: ";
    if (err)
      msg += err;
    Napi::Error::New(env, msg).ThrowAsJavaScriptException();
    return env.Undefined();
  }
  return env.Undefined();
}

static Napi::Value ParseMzML(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  if (info.Length() < 1 || !info[0].IsBuffer())
    return env.Undefined();

  Napi::Buffer<uint8_t> input = info[0].As<Napi::Buffer<uint8_t>>();
  size_t input_len = input.Length();

  MzML *handle = nullptr;
  int32_t rc = ABI.parse_mzml(input.Data(), input_len, &handle);

  if (rc != 0)
    return ThrowRc(env, "parse_mzml", rc);

  MzMLWrapper *w = new MzMLWrapper{handle, input_len};
  Napi::MemoryManagement::AdjustExternalMemory(env, (int64_t)input_len);
  return Napi::External<MzMLWrapper>::New(env, w, FinalizeMzML);
}

static Napi::Value BinToJson(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  MzML *handle = GetHandle(info[0]);
  if (!handle)
    return ThrowRc(env, "UseAfterFree/InvalidHandle", 0);

  OwnedBuf out;
  int32_t rc = ABI.bin_to_json(handle, out.Out());
  if (rc != 0)
    return ThrowRc(env, "bin_to_json", rc);

  return TakeUtf8String(env, out.Out());
}

static Napi::Value BinToMzML(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  MzML *handle = GetHandle(info[0]);
  if (!handle)
    return ThrowRc(env, "UseAfterFree/InvalidHandle", 0);

  OwnedBuf out;
  int32_t rc = ABI.bin_to_mzml(handle, out.Out());
  if (rc != 0)
    return ThrowRc(env, "bin_to_mzml", rc);

  return TakeUtf8String(env, out.Out());
}

static Napi::Value GetPeak(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  if (!ThrowIfMissing(env, (void *)ABI.get_peak, "get_peak"))
    return env.Undefined();
  if (!ThrowIfMissing(env, (void *)ABI.free_, "free_"))
    return env.Undefined();
  if (info.Length() < 4 || !info[0].IsTypedArray() || !info[1].IsTypedArray() || !info[2].IsNumber() || !info[3].IsNumber())
  {
    Napi::TypeError::New(env, "expected: (Float64Array x, Float64Array y, number targetRt, number rtRange, Buffer|null options?)").ThrowAsJavaScriptException();
    return env.Undefined();
  }

  Napi::Float64Array x_arr = info[0].As<Napi::Float64Array>();
  Napi::Float64Array y_arr = info[1].As<Napi::Float64Array>();
  double target_rt = info[2].As<Napi::Number>().DoubleValue();
  double rt_range = info[3].As<Napi::Number>().DoubleValue();

  size_t n = x_arr.ElementLength();
  if (y_arr.ElementLength() != n)
  {
    Napi::TypeError::New(env, "x/y length mismatch").ThrowAsJavaScriptException();
    return env.Undefined();
  }

  CPeakOptions opts;
  const CPeakOptions *p_opts = nullptr;
  if (info.Length() > 4)
    p_opts = ReadPeakOptionsObject(info[4], &opts);

  const double *x_ptr = Float64Ptr(x_arr);
  const double *y_ptr = Float64Ptr(y_arr);

  OwnedBuf out;
  int32_t rc = ABI.get_peak(x_ptr, y_ptr, n, target_rt, rt_range, p_opts, out.Out());
  if (rc != 0)
    return ThrowRc(env, "get_peak", rc);

  return TakeUtf8String(env, out.Out());
}

static Napi::Value CalculateEic(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();

  if (!ThrowIfMissing(env, (void *)ABI.calculate_eic, "calculate_eic") || !ThrowIfMissing(env, (void *)ABI.free_, "free_"))
    return env.Undefined();

  MzML *handle = GetHandle(info[0]);
  if (!handle)
    return ThrowRc(env, "UseAfterFree/InvalidHandle", 0);

  if (info.Length() < 6 || !info[1].IsNumber() || !info[2].IsNumber() ||
      !info[3].IsNumber() || !info[4].IsNumber() || !info[5].IsNumber())
  {
    Napi::TypeError::New(env, "expected: (External handle, number target, number fromRt, number toRt, number ppmTol, number mzTol)")
        .ThrowAsJavaScriptException();
    return env.Undefined();
  }

  double targets = info[1].As<Napi::Number>().DoubleValue();
  double from_rt = info[2].As<Napi::Number>().DoubleValue();
  double to_rt = info[3].As<Napi::Number>().DoubleValue();
  double ppm_tol = info[4].As<Napi::Number>().DoubleValue();
  double mz_tol = info[5].As<Napi::Number>().DoubleValue();

  OwnedBuf x_buf;
  OwnedBuf y_buf;

  int32_t rc = ABI.calculate_eic(handle, targets, from_rt, to_rt, ppm_tol, mz_tol, x_buf.Out(), y_buf.Out());

  if (rc != 0)
    return ThrowRc(env, "calculate_eic", rc);

  size_t nx = x_buf.buf.len / 8;
  size_t ny = y_buf.buf.len / 8;

  Napi::ArrayBuffer abx = TakeArrayBuffer(env, x_buf.Out());
  Napi::ArrayBuffer aby = TakeArrayBuffer(env, y_buf.Out());

  Napi::Float64Array X = Napi::Float64Array::New(env, nx, abx, 0);
  Napi::Float64Array Y = Napi::Float64Array::New(env, ny, aby, 0);

  Napi::Object out = Napi::Object::New(env);
  out.Set("x", X);
  out.Set("y", Y);

  return out;
}

static Napi::Value FindNoiseLevel(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  if (!ThrowIfMissing(env, (void *)ABI.find_noise_level, "find_noise_level"))
    return env.Undefined();
  if (info.Length() < 1 || !info[0].IsTypedArray())
  {
    Napi::TypeError::New(env, "expected: (Float64Array y)").ThrowAsJavaScriptException();
    return env.Undefined();
  }

  Napi::Float64Array y_arr = info[0].As<Napi::Float64Array>();
  const double *y_ptr = Float64Ptr(y_arr);
  double value = ABI.find_noise_level(y_ptr, y_arr.ElementLength());
  return Napi::Number::New(env, value);
}

static Napi::Value GetPeaksFromEic(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  if (info.Length() < 8)
    return env.Undefined();
  MzML *handle = GetHandle(info[0]);
  if (!handle)
    return ThrowRc(env, "UseAfterFree", 0);

  if (!info[1].IsTypedArray() || !info[2].IsTypedArray() || !info[3].IsTypedArray())
  {
    Napi::TypeError::New(env, "expected: (External, Float64Array rts, Float64Array mzs, Float64Array ranges, ...)").ThrowAsJavaScriptException();
    return env.Undefined();
  }

  Napi::Float64Array rts = info[1].As<Napi::Float64Array>();
  Napi::Float64Array mzs = info[2].As<Napi::Float64Array>();
  Napi::Float64Array rng = info[3].As<Napi::Float64Array>();
  size_t count = rts.ElementLength();
  PackedIds packed;
  BuildPackedIds(env, info[4], count, &packed);
  double f_l = info[5].As<Napi::Number>().DoubleValue();
  double t_r = info[6].As<Napi::Number>().DoubleValue();
  CPeakOptions opts;
  const CPeakOptions *p_opts = ReadPeakOptionsObject(info[7], &opts);
  size_t cores = info.Length() > 8 ? info[8].As<Napi::Number>().Uint32Value() : 1;
  OwnedBuf out;
  int32_t rc = ABI.get_peaks_from_eic(handle, Float64Ptr(rts), Float64Ptr(mzs), Float64Ptr(rng),
                                      packed.offs_ptr, packed.lens_ptr, packed.ids_buf_ptr, packed.ids_buf_len,
                                      count, f_l, t_r, p_opts, cores, out.Out());
  if (rc != 0)
    return ThrowRc(env, "get_peaks_from_eic", rc);
  return TakeUtf8String(env, out.Out());
}

static Napi::Value GetPeaksFromChrom(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  if (info.Length() < 5)
    return env.Undefined();

  MzML *handle = GetHandle(info[0]);
  if (!handle)
    return ThrowRc(env, "UseAfterFree", 0);

  if (!info[1].IsTypedArray() || !info[2].IsTypedArray() || !info[3].IsTypedArray())
  {
    Napi::TypeError::New(env, "expected: (External, Uint32Array indices, Float64Array rts, Float64Array windows, ...)").ThrowAsJavaScriptException();
    return env.Undefined();
  }

  Napi::Uint32Array idx = info[1].As<Napi::Uint32Array>();
  Napi::Float64Array rts = info[2].As<Napi::Float64Array>();
  Napi::Float64Array rng = info[3].As<Napi::Float64Array>();
  size_t count = rts.ElementLength();
  CPeakOptions opts;
  const CPeakOptions *p_opts = ReadPeakOptionsObject(info[4], &opts);
  size_t cores = info.Length() > 5 ? info[5].As<Napi::Number>().Uint32Value() : 1;
  OwnedBuf out;
  int32_t rc = ABI.get_peaks_from_chrom(handle, Uint32Ptr(idx), Float64Ptr(rts), Float64Ptr(rng), count, p_opts, cores, out.Out());
  if (rc != 0)
    return ThrowRc(env, "get_peaks_from_chrom", rc);
  return TakeUtf8String(env, out.Out());
}

static Napi::Value FindPeaks(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  if (!ThrowIfMissing(env, (void *)ABI.find_peaks, "find_peaks"))
    return env.Undefined();
  if (!ThrowIfMissing(env, (void *)ABI.free_, "free_"))
    return env.Undefined();
  if (info.Length() < 2 || !info[0].IsTypedArray() || !info[1].IsTypedArray())
  {
    Napi::TypeError::New(env, "expected: (Float64Array x, Float64Array y, Buffer|null options?)").ThrowAsJavaScriptException();
    return env.Undefined();
  }

  Napi::Float64Array x_arr = info[0].As<Napi::Float64Array>();
  Napi::Float64Array y_arr = info[1].As<Napi::Float64Array>();
  size_t n = x_arr.ElementLength();
  if (y_arr.ElementLength() != n)
  {
    Napi::TypeError::New(env, "x/y length mismatch").ThrowAsJavaScriptException();
    return env.Undefined();
  }

  CPeakOptions opts;
  const CPeakOptions *p_opts = nullptr;
  if (info.Length() > 2)
    p_opts = ReadPeakOptionsObject(info[2], &opts);

  const double *x_ptr = Float64Ptr(x_arr);
  const double *y_ptr = Float64Ptr(y_arr);

  OwnedBuf out;
  int32_t rc = ABI.find_peaks(x_ptr, y_ptr, n, p_opts, out.Out());
  if (rc != 0)
    return ThrowRc(env, "find_peaks", rc);

  return TakeUtf8String(env, out.Out());
}

static Napi::Value CalculateBaseline(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  if (!ThrowIfMissing(env, (void *)ABI.calculate_baseline, "calculate_baseline"))
    return env.Undefined();
  if (!ThrowIfMissing(env, (void *)ABI.free_, "free_"))
    return env.Undefined();
  if (info.Length() < 1 || !info[0].IsTypedArray())
  {
    Napi::TypeError::New(env, "expected: Float64Array").ThrowAsJavaScriptException();
    return env.Undefined();
  }

  Napi::Float64Array y_arr = info[0].As<Napi::Float64Array>();
  int32_t lambda = 0, max_iterations = 0;

  if (info.Length() >= 2 && info[1].IsObject() && !info[1].IsBuffer() && !info[1].IsTypedArray())
  {
    Napi::Object o = info[1].As<Napi::Object>();
    if (o.Has("lambda"))
      lambda = o.Get("lambda").ToNumber().Int32Value();
    if (o.Has("maxIterations"))
      max_iterations = o.Get("maxIterations").ToNumber().Int32Value();
  }
  else
  {
    if (info.Length() > 1 && info[1].IsNumber())
      lambda = info[1].As<Napi::Number>().Int32Value();
    if (info.Length() > 2 && info[2].IsNumber())
      max_iterations = info[2].As<Napi::Number>().Int32Value();
  }

  const double *y_ptr = Float64Ptr(y_arr);
  size_t n = y_arr.ElementLength();

  OwnedBuf out;
  int32_t rc = ABI.calculate_baseline(y_ptr, n, lambda, max_iterations, out.Out());
  if (rc != 0)
    return ThrowRc(env, "calculate_baseline", rc);

  size_t ny = out.buf.len / 8;
  Napi::ArrayBuffer aby = TakeArrayBuffer(env, out.Out());
  return Napi::Float64Array::New(env, ny, aby, 0);
}

static Napi::Value FindFeatures(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  if (info.Length() < 10)
    return env.Undefined();

  MzML *handle = GetHandle(info[0]);
  if (!handle)
    return ThrowRc(env, "UseAfterFree", 0);

  double from = info[1].As<Napi::Number>().DoubleValue();
  double to = info[2].As<Napi::Number>().DoubleValue();
  double ppm = info[3].As<Napi::Number>().DoubleValue();
  double mz = info[4].As<Napi::Number>().DoubleValue();
  double gs = info[5].As<Napi::Number>().DoubleValue();
  double ge = info[6].As<Napi::Number>().DoubleValue();
  double gst = info[7].As<Napi::Number>().DoubleValue();
  CPeakOptions opts;
  const CPeakOptions *p_opts = ReadPeakOptionsObject(info[8], &opts);
  int32_t cores = info[9].As<Napi::Number>().Int32Value();
  OwnedBuf out;
  int32_t rc = ABI.find_features(handle, from, to, ppm, mz, gs, ge, gst, p_opts, cores, out.Out());
  if (rc != 0)
    return ThrowRc(env, "find_features", rc);
  return TakeUtf8String(env, out.Out());
}

static Napi::Value FindFeature(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  if (info.Length() < 10)
    return env.Undefined();

  MzML *handle = GetHandle(info[0]);
  if (!handle)
    return ThrowRc(env, "UseAfterFree", 0);

  if (!info[1].IsTypedArray() || !info[2].IsTypedArray() || !info[3].IsTypedArray())
  {
    Napi::TypeError::New(env, "expected: (External, Float64Array rts, Float64Array mzs, Float64Array windows, ...)").ThrowAsJavaScriptException();
    return env.Undefined();
  }

  Napi::Float64Array rts = info[1].As<Napi::Float64Array>();
  Napi::Float64Array mzs = info[2].As<Napi::Float64Array>();
  Napi::Float64Array wins = info[3].As<Napi::Float64Array>();
  size_t count = rts.ElementLength();
  PackedIds packed;
  BuildPackedIds(env, info[4], count, &packed);
  size_t cores = info[5].As<Napi::Number>().Uint32Value();
  double s_ppm = info[6].As<Napi::Number>().DoubleValue();
  double s_mz = info[7].As<Napi::Number>().DoubleValue();
  double e_ppm = info[8].As<Napi::Number>().DoubleValue();
  double e_mz = info[9].As<Napi::Number>().DoubleValue();
  CPeakOptions opts;
  const CPeakOptions *p_opts = ReadPeakOptionsObject(info[10], &opts);
  OwnedBuf out;
  int32_t rc = ABI.find_feature(handle, Float64Ptr(rts), Float64Ptr(mzs), Float64Ptr(wins), packed.offs_ptr, packed.lens_ptr, packed.ids_buf_ptr, packed.ids_buf_len, count, cores, s_ppm, s_mz, e_ppm, e_mz, p_opts, out.Out());
  if (rc != 0)
    return ThrowRc(env, "find_feature", rc);
  return TakeUtf8String(env, out.Out());
}

static Napi::Value MzmlToBin(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  MzML *handle = GetHandle(info[0]);
  if (!handle)
    return ThrowRc(env, "UseAfterFree/InvalidHandle", 0);

  int32_t lv = info[1].As<Napi::Number>().Int32Value();
  int32_t c = info[2].As<Napi::Number>().Int32Value();

  OwnedBuf out;
  int32_t rc = ABI.mzml_to_bin(handle, out.Out(), (uint8_t)lv, (uint8_t)(c != 0 ? 1 : 0));
  if (rc != 0)
    return ThrowRc(env, "mzml_to_bin", rc);

  return TakeBuffer(env, out.Out());
}

static Napi::Value ParseBin(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  if (info.Length() < 1 || !info[0].IsBuffer())
    return env.Undefined();

  Napi::Buffer<uint8_t> input = info[0].As<Napi::Buffer<uint8_t>>();
  size_t input_len = input.Length();
  size_t max_cache = info.Length() > 1 && info[1].IsNumber()
                         ? (size_t)info[1].As<Napi::Number>().Int64Value()
                         : 0;

  MzML *handle = nullptr;
  int32_t rc = ABI.parse_bin(input.Data(), input_len, max_cache, &handle);

  if (rc != 0)
    return ThrowRc(env, "parse_bin", rc);

  MzMLWrapper *w = new MzMLWrapper{handle, input_len};
  Napi::MemoryManagement::AdjustExternalMemory(env, (int64_t)input_len);
  return Napi::External<MzMLWrapper>::New(env, w, FinalizeMzML);
}

static Napi::Value ParseIonPath(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  if (!ThrowIfMissing(env, (void *)ABI.parse_ion_path, "parse_ion_path"))
    return env.Undefined();

  if (info.Length() < 1 || !info[0].IsString())
  {
    Napi::TypeError::New(env, "expected: (string path, number maxCacheSize?)").ThrowAsJavaScriptException();
    return env.Undefined();
  }

  std::string file_path = info[0].As<Napi::String>().Utf8Value();
  size_t cache_size = info.Length() > 1 && info[1].IsNumber()
                          ? (size_t)info[1].As<Napi::Number>().Int64Value()
                          : 0;

  MzML *file_handle = nullptr;
  int32_t code = ABI.parse_ion_path(file_path.c_str(), cache_size, &file_handle);
  if (code != 0)
    return ThrowRc(env, "parse_ion_path", code);

  MzMLWrapper *file_wrapper = new MzMLWrapper{file_handle, 0};
  return Napi::External<MzMLWrapper>::New(env, file_wrapper, FinalizeMzML);
}

static Napi::Value ParseIonUrl(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  if (!ThrowIfMissing(env, (void *)ABI.parse_ion_url, "parse_ion_url"))
    return env.Undefined();

  if (info.Length() < 1 || !info[0].IsString())
  {
    Napi::TypeError::New(env, "expected: (string url, number maxCacheSize?)").ThrowAsJavaScriptException();
    return env.Undefined();
  }

  std::string url = info[0].As<Napi::String>().Utf8Value();
  size_t cache_size = info.Length() > 1 && info[1].IsNumber()
                          ? (size_t)info[1].As<Napi::Number>().Int64Value()
                          : 0;

  MzML *file_handle = nullptr;
  int32_t code = ABI.parse_ion_url(url.c_str(), cache_size, &file_handle);
  if (code != 0)
    return ThrowRc(env, "parse_ion_url", code);

  MzMLWrapper *file_wrapper = new MzMLWrapper{file_handle, 0};
  return Napi::External<MzMLWrapper>::New(env, file_wrapper, FinalizeMzML);
}

static Napi::Value GetFeatures(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  if (info.Length() < 14 || !info[0].IsString())
  {
    Napi::TypeError::New(env, "Invalid arguments for getFeatures").ThrowAsJavaScriptException();
    return env.Undefined();
  }

  std::string dir_path = info[0].As<Napi::String>().Utf8Value();
  double from = info[1].As<Napi::Number>().DoubleValue();
  double to = info[2].As<Napi::Number>().DoubleValue();
  double eic_ppm = info[3].As<Napi::Number>().DoubleValue();
  double eic_mz = info[4].As<Napi::Number>().DoubleValue();
  double g_start = info[5].As<Napi::Number>().DoubleValue();
  double g_end = info[6].As<Napi::Number>().DoubleValue();
  double g_step = info[7].As<Napi::Number>().DoubleValue();
  double group_ppm = info[8].As<Napi::Number>().DoubleValue();
  double group_da = info[9].As<Napi::Number>().DoubleValue();
  double group_rt = info[10].As<Napi::Number>().DoubleValue();
  int32_t frequency = info[11].As<Napi::Number>().Int32Value();
  CPeakOptions opts;
  const CPeakOptions *p_opts = ReadPeakOptionsObject(info[12], &opts);
  int32_t cores = info[13].As<Napi::Number>().Int32Value();

  OwnedBuf out;
  int32_t rc = ABI.get_features(
      dir_path.c_str(),
      from, to, eic_ppm, eic_mz,
      g_start, g_end, g_step,
      group_ppm, group_da, group_rt,
      frequency, p_opts, cores,
      out.Out());

  if (rc != 0)
    return ThrowRc(env, "get_features", rc);
  return TakeUtf8String(env, out.Out());
}

static Napi::Value GetScans(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  if (!ThrowIfMissing(env, (void *)ABI.get_scans, "get_scans"))
    return env.Undefined();

  MzML *handle = GetHandle(info[0]);
  if (!handle)
    return ThrowRc(env, "UseAfterFree/InvalidHandle", 0);

  if (info.Length() < 5 || !info[1].IsNumber() || !info[2].IsNumber() || !info[3].IsNumber() || !info[4].IsNumber())
  {
    Napi::TypeError::New(env, "expected: (External handle, number queryType, number a, number b, number level)")
        .ThrowAsJavaScriptException();
    return env.Undefined();
  }

  uint8_t query_type = (uint8_t)info[1].As<Napi::Number>().Uint32Value();
  double a = info[2].As<Napi::Number>().DoubleValue();
  double b = info[3].As<Napi::Number>().DoubleValue();
  uint8_t level = (uint8_t)info[4].As<Napi::Number>().Uint32Value();

  OwnedBuf out;
  int32_t rc = ABI.get_scans(handle, query_type, a, b, level, out.Out());
  if (rc != 0)
    return ThrowRc(env, "get_scans", rc);

  return TakeUtf8String(env, out.Out());
}

static Napi::Value GetIonImage(const Napi::CallbackInfo &info)
{
  Napi::Env env = info.Env();
  if (!ThrowIfMissing(env, (void *)ABI.get_ion_image, "get_ion_image"))
    return env.Undefined();

  MzML *handle = GetHandle(info[0]);
  if (!handle)
    return ThrowRc(env, "UseAfterFree/InvalidHandle", 0);

  if (info.Length() < 4 || !info[1].IsNumber() || !info[2].IsNumber() || !info[3].IsNumber())
  {
    Napi::TypeError::New(env, "expected: (External handle, number mz, number tolerance, number level)")
        .ThrowAsJavaScriptException();
    return env.Undefined();
  }

  double mz = info[1].As<Napi::Number>().DoubleValue();
  double tolerance = info[2].As<Napi::Number>().DoubleValue();
  uint8_t level = (uint8_t)info[3].As<Napi::Number>().Uint32Value();

  OwnedBuf out;
  int32_t rc = ABI.get_ion_image(handle, mz, tolerance, level, out.Out());
  if (rc != 0)
    return ThrowRc(env, "get_ion_image", rc);

  return TakeUtf8String(env, out.Out());
}

static Napi::Object Init(Napi::Env env, Napi::Object exports)
{
  exports.Set("bind", Napi::Function::New(env, Bind));
  exports.Set("parseMzML", Napi::Function::New(env, ParseMzML));
  exports.Set("binToJson", Napi::Function::New(env, BinToJson));
  exports.Set("binToMzML", Napi::Function::New(env, BinToMzML));
  exports.Set("getPeak", Napi::Function::New(env, GetPeak));
  exports.Set("calculateEic", Napi::Function::New(env, CalculateEic));
  exports.Set("findNoiseLevel", Napi::Function::New(env, FindNoiseLevel));
  exports.Set("getPeaksFromEic", Napi::Function::New(env, GetPeaksFromEic));
  exports.Set("getPeaksFromChrom", Napi::Function::New(env, GetPeaksFromChrom));
  exports.Set("findPeaks", Napi::Function::New(env, FindPeaks));
  exports.Set("calculateBaseline", Napi::Function::New(env, CalculateBaseline));
  exports.Set("findFeatures", Napi::Function::New(env, FindFeatures));
  exports.Set("findFeature", Napi::Function::New(env, FindFeature));
  exports.Set("mzmlToBin", Napi::Function::New(env, MzmlToBin));
  exports.Set("parseBin", Napi::Function::New(env, ParseBin));
  exports.Set("parseIonPath", Napi::Function::New(env, ParseIonPath));
  exports.Set("parseIonUrl", Napi::Function::New(env, ParseIonUrl));
  exports.Set("dispose", Napi::Function::New(env, DisposeMzML));
  exports.Set("getFeatures", Napi::Function::New(env, GetFeatures));
  exports.Set("getScans", Napi::Function::New(env, GetScans));
  exports.Set("getIonImage", Napi::Function::New(env, GetIonImage));
  return exports;
}

NODE_API_MODULE(msutils, Init)
