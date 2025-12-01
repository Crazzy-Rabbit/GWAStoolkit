#ifndef GWASTOOLKIT_statfunc_HPP
#define GWASTOOLKIT_statfunc_HPP

#include <cmath>

namespace StatFunc {
    // 标准正态概率密度 φ(x)
    double dnorm(double x);

    // 标准正态分布 upper-tail 概率：P(Z >= x)
    double pnorm_upper(double x);

    // qnorm_sub：牛顿迭代子函数
    double qnorm_sub(double x, double y);

    // 反正态分布（lower-tail / upper-tail 都支持）
    double qnorm(double p, bool upper = false);

    // 双尾 p → z（两侧）
    double p2z_two_tailed(double p);

    // 单尾 p → z（左尾 lower tail）
    double p2z_lower(double p);

    // 单尾 p → z（右尾 upper tail）
    double p2z_upper(double p);
}

#endif