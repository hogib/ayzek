#pragma once

// The few dense kernels inference spends its time in, with a NEON path for
// aarch64 (Raspberry Pi 4/5) and a portable scalar path everywhere else.
//
// Every model layer bottoms out in `dot`, `gemv` or an elementwise op, so this
// is the only file that knows which instruction set is in use. The scalar path
// is written so GCC and Clang auto-vectorise it at -O3 on x86; the NEON path is
// explicit because the Pi is the deployment target and its speed should not
// depend on what a compiler felt like doing.
//
// Both paths compute the same sums in a different association order, so they
// agree to floating-point rounding, not bit for bit. tests/test_simd.cpp holds
// them to that.

#include <cstddef>

#if defined(__ARM_NEON) || defined(__aarch64__)
#include <arm_neon.h>
#define AYZEK_NEON 1
#else
#define AYZEK_NEON 0
#endif

namespace ayzek::simd {

inline constexpr const char* kBackend = AYZEK_NEON ? "neon" : "scalar";

// sum(a[i] * b[i])
inline float dot(const float* a, const float* b, std::size_t n) noexcept {
    std::size_t i = 0;
#if AYZEK_NEON
    float32x4_t s0 = vdupq_n_f32(0), s1 = vdupq_n_f32(0), s2 = vdupq_n_f32(0), s3 = vdupq_n_f32(0);
    for (; i + 16 <= n; i += 16) {
        s0 = vmlaq_f32(s0, vld1q_f32(a + i), vld1q_f32(b + i));
        s1 = vmlaq_f32(s1, vld1q_f32(a + i + 4), vld1q_f32(b + i + 4));
        s2 = vmlaq_f32(s2, vld1q_f32(a + i + 8), vld1q_f32(b + i + 8));
        s3 = vmlaq_f32(s3, vld1q_f32(a + i + 12), vld1q_f32(b + i + 12));
    }
    for (; i + 4 <= n; i += 4) s0 = vmlaq_f32(s0, vld1q_f32(a + i), vld1q_f32(b + i));
    float s = vaddvq_f32(vaddq_f32(vaddq_f32(s0, s1), vaddq_f32(s2, s3)));
#else
    float s = 0.0f;
#endif
    for (; i < n; ++i) s += a[i] * b[i];
    return s;
}

// y[r] = bias[r] + sum_c W[r, c] * x[c], W row-major (rows x cols). bias may be null.
inline void gemv(const float* W, std::size_t rows, std::size_t cols, const float* x,
                 const float* bias, float* y) noexcept {
    for (std::size_t r = 0; r < rows; ++r) {
        y[r] = dot(W + r * cols, x, cols) + (bias ? bias[r] : 0.0f);
    }
}

// y += x
inline void add(const float* x, float* y, std::size_t n) noexcept {
    std::size_t i = 0;
#if AYZEK_NEON
    for (; i + 4 <= n; i += 4) vst1q_f32(y + i, vaddq_f32(vld1q_f32(y + i), vld1q_f32(x + i)));
#endif
    for (; i < n; ++i) y[i] += x[i];
}

// y += a * x
inline void axpy(float a, const float* x, float* y, std::size_t n) noexcept {
    std::size_t i = 0;
#if AYZEK_NEON
    const float32x4_t va = vdupq_n_f32(a);
    for (; i + 4 <= n; i += 4) vst1q_f32(y + i, vmlaq_f32(vld1q_f32(y + i), va, vld1q_f32(x + i)));
#endif
    for (; i < n; ++i) y[i] += a * x[i];
}

// x = x * s + b
inline void affine(float* x, float s, float b, std::size_t n) noexcept {
    std::size_t i = 0;
#if AYZEK_NEON
    const float32x4_t vs = vdupq_n_f32(s), vb = vdupq_n_f32(b);
    for (; i + 4 <= n; i += 4) vst1q_f32(x + i, vaddq_f32(vmulq_f32(vld1q_f32(x + i), vs), vb));
#endif
    for (; i < n; ++i) x[i] = x[i] * s + b;
}

// x = max(x, 0)
inline void relu(float* x, std::size_t n) noexcept {
    std::size_t i = 0;
#if AYZEK_NEON
    const float32x4_t z = vdupq_n_f32(0);
    for (; i + 4 <= n; i += 4) vst1q_f32(x + i, vmaxq_f32(vld1q_f32(x + i), z));
#endif
    for (; i < n; ++i) x[i] = x[i] > 0.0f ? x[i] : 0.0f;
}

// y = a * x  (elementwise), y may alias x
inline void mul(const float* a, const float* x, float* y, std::size_t n) noexcept {
    std::size_t i = 0;
#if AYZEK_NEON
    for (; i + 4 <= n; i += 4) vst1q_f32(y + i, vmulq_f32(vld1q_f32(a + i), vld1q_f32(x + i)));
#endif
    for (; i < n; ++i) y[i] = a[i] * x[i];
}

// sum(x) over doubles, for window statistics.
inline double sum(const double* x, std::size_t n) noexcept {
    std::size_t i = 0;
#if AYZEK_NEON
    float64x2_t s0 = vdupq_n_f64(0), s1 = vdupq_n_f64(0);
    for (; i + 4 <= n; i += 4) {
        s0 = vaddq_f64(s0, vld1q_f64(x + i));
        s1 = vaddq_f64(s1, vld1q_f64(x + i + 2));
    }
    double s = vaddvq_f64(vaddq_f64(s0, s1));
#else
    double s = 0.0;
#endif
    for (; i < n; ++i) s += x[i];
    return s;
}

// sum(x * y) over doubles.
inline double dot(const double* x, const double* y, std::size_t n) noexcept {
    std::size_t i = 0;
#if AYZEK_NEON
    float64x2_t s0 = vdupq_n_f64(0), s1 = vdupq_n_f64(0);
    for (; i + 4 <= n; i += 4) {
        s0 = vmlaq_f64(s0, vld1q_f64(x + i), vld1q_f64(y + i));
        s1 = vmlaq_f64(s1, vld1q_f64(x + i + 2), vld1q_f64(y + i + 2));
    }
    double s = vaddvq_f64(vaddq_f64(s0, s1));
#else
    double s = 0.0;
#endif
    for (; i < n; ++i) s += x[i] * y[i];
    return s;
}

// x *= a  (elementwise) over doubles
inline void mul_d(const double* a, double* x, std::size_t n) noexcept {
    std::size_t i = 0;
#if AYZEK_NEON
    for (; i + 2 <= n; i += 2) vst1q_f64(x + i, vmulq_f64(vld1q_f64(a + i), vld1q_f64(x + i)));
#endif
    for (; i < n; ++i) x[i] *= a[i];
}

// x = x * s + b over doubles
inline void affine(double* x, double s, double b, std::size_t n) noexcept {
    std::size_t i = 0;
#if AYZEK_NEON
    const float64x2_t vs = vdupq_n_f64(s), vb = vdupq_n_f64(b);
    for (; i + 2 <= n; i += 2) vst1q_f64(x + i, vaddq_f64(vmulq_f64(vld1q_f64(x + i), vs), vb));
#endif
    for (; i < n; ++i) x[i] = x[i] * s + b;
}

}  // namespace ayzek::simd
