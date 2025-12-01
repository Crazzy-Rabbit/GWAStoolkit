#include "StatFunc.hpp"
#include <cmath>
#include <algorithm>

namespace StatFunc {
// φ(x) = exp(-x^2/2)/sqrt(2π)
double dnorm(double x) {
    return 0.3989422804014327 * std::exp(-0.5 * x * x);
}

// Φ̅(x) = P(Z ≥ x)
double pnorm_upper(double x) {
    // 来源：StatFunc::pnorm（upper-tail）
    double z = (x > 0 ? -x : x);
    double t = 1 / (1 + 0.2316419 * std::fabs(z));
    double poly = (((1.330274429 * t - 1.821255978) * t + 1.781477937)
                    * t - 0.356563782) * t + 0.31938153;
    double nd = 0.3989422804014327 * std::exp(-0.5 * z * z);
    double p = nd * t * poly;
    return x >= 0 ? p : 1 - p;
}

// Newton-Raphson iteration step
double qnorm_sub(double x, double y) {
    return y + 0.5 * x * y * y +
           (2 * x * x + 1) * y * y * y / 6.0 +
           (6 * x * x * x + 7 * x) * y * y * y * y / 12.0;
}

// 反正态分布
double qnorm(double p, bool upper) {
    if (upper) p = 1.0 - p;

    double x = 0; // 初始猜测
    for (int i = 0; i < 4; i++)
        x = x + qnorm_sub(x, (pnorm_upper(x) - p) / dnorm(x));

    return x;
}

// =======================
// p→z 转换（两侧）
// =======================
double p2z_two_tailed(double p) {
    if (p <= 0) return 38.0;
    if (p >= 1) return 0.0;
    double z = qnorm(p / 2.0, false);
    return std::fabs(z);
}

// 左尾 lower tail
double p2z_lower(double p) {
    if (p <= 0) return -38;
    if (p >= 1) return 0.0;
    return qnorm(p, false);
}

// 右尾 upper tail
double p2z_upper(double p) {
    if (p <= 0) return 38;
    if (p >= 1) return 0.0;
    return qnorm(p, true);
}

}
