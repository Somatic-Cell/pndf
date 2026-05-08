#pragma once
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>

namespace pndf {

constexpr double kEps = 1e-12;

struct Vec2 {
    double x = 0.0;
    double y = 0.0;
    constexpr Vec2() = default;
    constexpr Vec2(double x_, double y_) : x(x_), y(y_) {}
    constexpr Vec2 operator+(const Vec2& b) const { return {x + b.x, y + b.y}; }
    constexpr Vec2 operator-(const Vec2& b) const { return {x - b.x, y - b.y}; }
    constexpr Vec2 operator*(double s) const { return {x * s, y * s}; }
    constexpr Vec2 operator/(double s) const { return {x / s, y / s}; }
    Vec2& operator+=(const Vec2& b) { x += b.x; y += b.y; return *this; }
};

struct Vec4 {
    double x = 0.0, y = 0.0, z = 0.0, w = 0.0;
    constexpr Vec4() = default;
    constexpr Vec4(double a, double b, double c, double d) : x(a), y(b), z(c), w(d) {}
    double& operator[](int i) { return (&x)[i]; }
    double operator[](int i) const { return (&x)[i]; }
    Vec4 operator+(const Vec4& b) const { return {x+b.x,y+b.y,z+b.z,w+b.w}; }
    Vec4 operator-(const Vec4& b) const { return {x-b.x,y-b.y,z-b.z,w-b.w}; }
    Vec4 operator*(double s) const { return {x*s,y*s,z*s,w*s}; }
};

inline double dot(const Vec2& a, const Vec2& b) { return a.x*b.x + a.y*b.y; }
inline double norm2(const Vec2& a) { return dot(a,a); }
inline double norm(const Vec2& a) { return std::sqrt(norm2(a)); }
inline double dot4(const Vec4& a, const Vec4& b) { return a.x*b.x + a.y*b.y + a.z*b.z + a.w*b.w; }

inline double wrap01(double x) {
    x = std::fmod(x, 1.0);
    if (x < 0.0) x += 1.0;
    return x;
}

inline double wrap_delta(double d) {
    if (d > 0.5) d -= 1.0;
    if (d < -0.5) d += 1.0;
    return d;
}

inline Vec2 wrap_delta(Vec2 d) { return {wrap_delta(d.x), wrap_delta(d.y)}; }
inline Vec2 fract(Vec2 v) { return {wrap01(v.x), wrap01(v.y)}; }

inline Vec2 unwrap_near(Vec2 p, Vec2 reference) {
    return reference + wrap_delta(p - reference);
}

inline Vec2 integer_chart_shift(Vec2 source_origin, Vec2 target_origin) {
    return {std::round(source_origin.x - target_origin.x),
            std::round(source_origin.y - target_origin.y)};
}

inline Vec2 periodic_midpoint(Vec2 a, Vec2 b) {
    Vec2 d = wrap_delta(b - a);
    return fract(a + d * 0.5);
}

inline Vec2 periodic_midpoint_unwrapped(Vec2 a, Vec2 b) {
    return (a + unwrap_near(b, a)) * 0.5;
}

inline Vec2 clamp_projected_normal(Vec2 n) {
    const double r2 = n.x*n.x + n.y*n.y;
    if (r2 > 0.998001) {
        const double s = 0.999 / std::sqrt(r2);
        return {n.x*s, n.y*s};
    }
    return n;
}

inline double angular_error_deg(Vec2 a, Vec2 b) {
    a = clamp_projected_normal(a);
    b = clamp_projected_normal(b);
    const double za = std::sqrt(std::max(0.0, 1.0 - a.x*a.x - a.y*a.y));
    const double zb = std::sqrt(std::max(0.0, 1.0 - b.x*b.x - b.y*b.y));
    double d = a.x*b.x + a.y*b.y + za*zb;
    d = std::max(-1.0, std::min(1.0, d));
    return std::acos(d) * 180.0 / 3.14159265358979323846;
}

struct Mat5 {
    std::array<double, 25> a{};
    double& operator()(int r, int c) { return a[5*r + c]; }
    double operator()(int r, int c) const { return a[5*r + c]; }
    void add_outer(const std::array<double,5>& p, double weight = 1.0) {
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                (*this)(i,j) += weight * p[i] * p[j];
            }
        }
    }
    Mat5& operator+=(const Mat5& b) {
        for (int i=0;i<25;++i) a[i] += b.a[i];
        return *this;
    }
};

inline Mat5 operator+(Mat5 lhs, const Mat5& rhs) { lhs += rhs; return lhs; }

inline double eval_quadric(const Mat5& q, const std::array<double,5>& z) {
    double v = 0.0;
    for (int i=0;i<5;++i) for (int j=0;j<5;++j) v += z[i] * q(i,j) * z[j];
    return v;
}

inline Mat5 translate_uv_quadric(const Mat5& q, Vec2 t) {
    // Re-express a homogeneous quadric under x_old = x_new + t,
    // where x = (u,v,nx,ny).  If z_old = M z_new, the transformed
    // quadric is M^T Q M.  Only the u/v coordinates are translated.
    double M[5][5]{};
    for (int i = 0; i < 5; ++i) M[i][i] = 1.0;
    M[0][4] = t.x;
    M[1][4] = t.y;

    double tmp[5][5]{};
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            for (int k = 0; k < 5; ++k) tmp[i][j] += q(i,k) * M[k][j];
        }
    }

    Mat5 out;
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            for (int k = 0; k < 5; ++k) out(i,j) += M[k][i] * tmp[k][j];
        }
    }
    return out;
}

inline Mat5 rebase_uv_quadric(const Mat5& q, Vec2 source_origin, Vec2 target_origin) {
    return translate_uv_quadric(q, integer_chart_shift(source_origin, target_origin));
}

inline bool solve2x2(double a00, double a01, double a11, double b0, double b1, double& x0, double& x1) {
    const double det = a00*a11 - a01*a01;
    if (std::abs(det) < 1e-18) return false;
    x0 = ( a11*b0 - a01*b1) / det;
    x1 = (-a01*b0 + a00*b1) / det;
    return true;
}

} // namespace pndf
