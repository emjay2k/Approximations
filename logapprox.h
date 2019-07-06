/*
 * logapprox.cpp
 *
 *  Created on: Jun 13, 2019
 *      Author: Matthias Jung
 */

#ifndef SRC_LOGAPPROX_H_
#define SRC_LOGAPPROX_H_

#include <cmath>
#include <limits>
#include <boost/config.hpp>

// inspired by the paper "New close form approximations of ln(1+x) by Khattri
// uses approximation (a*x + b)/(c*x + d)
// parameters were fitted using a linear program with 100000 samples from range [0.5-1.0]
// Final evaluation was done using 1e11 values from [1.0, 2.0]
// speedup over log2_avx implementation ~ 3;
// Max error: ~0.00147 (accurate to 9 bits)
// Theoretical max error from the extreme values of the function: 0.001464579127038346
template <typename T>
inline T fastLog2p1(const T value) {
    const T a = 1.4767235475800453;
    const T b = -1.477808113688585;
    const T c = 0.60987486544988612;
    const T d = 0.43559347328148307;

    if(BOOST_UNLIKELY(!std::isfinite(value))) {
        return value == std::numeric_limits<T>::infinity() ? value : nan("1");
    } else if(BOOST_LIKELY(value > 0)) {
        int iExp;
        T dM = frexp(value, &iExp);

        T x = 1.0/(c*dM + d);

        return iExp + x * (a*dM + b);
    } else {
        return value == 0 ? -std::numeric_limits<T>::infinity() : -std::nan("1");
    }
}

// inspired by the paper "New close form approximations of ln(1+x) by Khattri
// uses approximation (a*x^2 + b*x + c)/(d*x^2 + e*x + f)
// parameters were fitted using a linear program with 100000 samples from range [0.5-1.0]
// Final evaluation was done using 1e11 values from [1.0, 2.0]
// speedup over log2_avx implementation ~ 2.7;
// Max error: ~3.46e-06 (accurate to 18 bits)
// Theoretical max error from the extreme values of the function: 3.458795617805599e-06
template <typename T>
inline T fastLog2p2(const T value) {
    const T a = 1.9127166899499954;
    const T b = -0.68851400593499545;
    const T c = -1.22420645509838;
    const T d = 0.49463685172392841;
    const T e = 1.426594307123505;
    const T f = 0.2533316901691966;

    if(BOOST_UNLIKELY(!std::isfinite(value))) {
        return value == std::numeric_limits<T>::infinity() ? value : nan("1");
    } else if(BOOST_LIKELY(value > 0)) {
        int iExp;
        T dM = frexp(value, &iExp);
        T dM2 = dM * dM;

        T x = 1.0/(d*dM2 + e*dM + f);

        return iExp + x * (a*dM2 + b*dM + c);
    } else {
        return value == 0 ? -std::numeric_limits<T>::infinity() : -std::nan("1");
    }
}

// inspired by the paper "New close form approximations of ln(1+x) by Khattri
// uses approximation (a*x^3 + b*x^2 + c*x + d)/(e*x^3 + f*x^2 + g*x + h)
// parameters were fitted using a linear program with 100000 samples from range [0.5-1.0]
// Final evaluation was done using 1e11 values from [1.0, 2.0]
// speedup over log2_avx implementation ~ 2.2;
// Max error: ~7.79e-09 (accurate to almost 27 bits)
// Theoretical max error from the extreme values of the function: 7.786322031577697e-09
template <typename T>
inline T fastLog2p3(const T value) {
    const T a = 1.1098414161667869;
    const T b = 1.4491119665946153;
    const T c = -2.0697678829202806;
    const T d = -0.48918550780729392;
    const T e = 0.22977948696488379;
    const T f = 1.4961611668393175;
    const T g = 1.071708023446889;
    const T h = 0.084444549259932208;

    if(BOOST_UNLIKELY(!std::isfinite(value))) {
        return value == std::numeric_limits<T>::infinity() ? value : nan("1");
    } else if(BOOST_LIKELY(value > 0)) {
        int iExp;
        T dM = frexp(value, &iExp);
        T dM2 = dM * dM;
        T dM3 = dM * dM2;

        T x = 1.0/(e*dM3 + f*dM2 + g*dM + h);

        return iExp + x * (a*dM3 + b*dM2 + c*dM + d);
    } else {
        return value == 0 ? -std::numeric_limits<T>::infinity() : -std::nan("1");
    }
}

// inspired by the paper "New close form approximations of ln(1+x) by Khattri
// uses approximation (a*x^4 + b*x^3 + c*x^2 + d*x + e)/(f*x^4 + g*x^3 + h*x^2 + i*x + j)
// parameters were fitted using a linear program with 100000 samples from range [0.5-1.0]
// Final evaluation was done using 1e11 values from [1.0, 2.0]
// speedup over log2_avx implementation ~ 2;
// Max error: ~1.77e-11 (accurate to almost 36 bits)
// Theoretical max error from the extreme values of the function: 1.772559876656032e-11
template <typename T>
inline T fastLog2p4(const T value) {
    const T a = 0.59329970349044314;
    const T b = 2.3979646338966889;
    const T c = -0.96358966800238843;
    const T d = -1.8439274267589987;
    const T e = -0.18374724264449727;
    const T f = 0.1068562844523792;
    const T g = 1.2392957064266512;
    const T h = 2.0062979261642901;
    const T i = 0.63680961689938775;
    const T j = 0.028211791264274255;

    if(BOOST_UNLIKELY(!std::isfinite(value))) {
        return value == std::numeric_limits<T>::infinity() ? value : nan("1");
    } else if(BOOST_LIKELY(value > 0)) {
        int iExp;
        T dM = frexp(value, &iExp);
        T dM2 = dM * dM;
        T dM3 = dM * dM2;
        T dM4 = dM2 * dM2;

        T x = 1.0/(f*dM4 + g*dM3 + h*dM2 + i*dM + j);

        return iExp + x * (a*dM4 + b*dM3 + c*dM2 + d*dM + e);
    } else {
        return value == 0 ? -std::numeric_limits<T>::infinity() : -std::nan("1");
    }
}

// inspired by the paper "New close form approximations of ln(1+x) by Khattri
// uses approximation (a*x^5 + b*x^4 + c*x^3 + d*x^2 + e*x + f)/(g*x^5 + h*x^4 + i*x^3 + j*x^2 + k*x + l)
// parameters were fitted using ceres on 38 values from range [0.5, 1.0] and preliminary evaluated with 1e7 values from
// the same range. This way a total set of 52 results with a better fitness of 4e-13 was collected and taken as input
// population of differential evolution algorithm for further optimization regarding max error.
// Final evaluation was done using 1e11 values from [1.0, 2.0]
// speedup over log2_avx implementation ~ 1.7;
// Max error: ~1.85e-14 (accurate to 45 bits)
// Theoretical max error from the extreme values of the function: 1.820765760385257e-14
template <typename T>
inline T fastLog2p5(const T value) {
    const T a = 1.000000000000000000000e+00;
    const T b = 7.724129062840519033273e+00;
    const T c = 3.878748300736706067227e+00;
    const T d = -8.633584487586309919038e+00;
    const T e = -3.763821035138447346213e+00;
    const T f = -2.054718408527802797625e-01;
    const T g = 1.636794214704399841320e-01;
    const T h = 2.927552384822291742239e+00;
    const T i = 8.327974821567529772892e+00;
    const T j = 5.887266723437919679895e+00;
    const T k = 1.034256519609611402188e+00;
    const T l = 2.889807696644808646114e-02;

    if(BOOST_UNLIKELY(!std::isfinite(value))) {
        return value == std::numeric_limits<T>::infinity() ? value : nan("1");
    } else if(BOOST_LIKELY(value > 0)) {
        int iExp;
        T dM = frexp(value, &iExp);
        T dM2 = dM * dM;
        T dM3 = dM * dM2;
        T dM4 = dM2 * dM2;
        T dM5 = dM2 * dM3;

        T x = 1.0/(g*dM5 + h*dM4 + i*dM3 + j*dM2 + k*dM + l);

        return iExp + x * (a*dM5 + b*dM4 + c*dM3 + d*dM2 + e*dM + f);
    } else {
        return value == 0 ? -std::numeric_limits<T>::infinity() : -std::nan("1");
    }
}

// inspired by the paper "New close form approximations of ln(1+x) by Khattri
// uses approximation (a*x^6 + b*x^5 + c*x^4 + d*x^3 + e*x^2 + f*x + g)/(h*x^6 + i*x^5 + j*x^4 + k*x^3 + l*x^2 + m*x + n)
// parameters were fitted using ceres on 38 values from range [0.5,1.0] and preliminary evaluated with 1e7 values from
// the same range. Final evaluation was done using 1e11 values from [1.0,2.0]
// speedup over log2_avx implementation ~ 1.6;
// Max error: 5.73e-16 (accurate to 50 bits)
// Theoretical max error from the extreme values of the function: 1.140580685454751e-16
template <typename T>
inline T fastLog2p6(const T value) {
    const T a = 1.000000000000000000000e+00L;
    const T b = 1.268415661999539878479e+01L;
    const T c = 2.113404189541518007900e+01L;
    const T d = -1.096780633851080466457e+01L;
    const T e = -1.945948524176642990824e+01L;
    const T f = -4.241771817238755737378e+00L;
    const T g = -1.491351178945869715431e-01L;
    const T h = 1.515072667524440808151e-01L;
    const T i = 3.930936965120734871704e+00L;
    const T j = 1.761700191252228364647e+01L;
    const T k = 2.234352453493684720343e+01L;
    const T l = 8.920637882427030262988e+00L;
    const T m = 1.007889789099811306272e+00L;
    const T n = 1.966898367167627115215e-02L;

    if(BOOST_UNLIKELY(!std::isfinite(value))) {
        return value == std::numeric_limits<T>::infinity() ? value : nan("1");
    } else if(BOOST_LIKELY(value > 0)) {
        int iExp;
        T dM = frexp(value, &iExp);
        T dM2 = dM * dM;
        T dM3 = dM * dM2;
        T dM4 = dM2 * dM2;
        T dM5 = dM2 * dM3;
        T dM6 = dM3 * dM3;

        T x = 1.0/(h*dM6 + i*dM5 + j*dM4 + k*dM3 + l*dM2 + m*dM + n);

        return iExp + x * (a*dM6 + b*dM5 + c*dM4 + d*dM3 + e*dM2 + f*dM + g);
    } else {
        return value == 0 ? -std::numeric_limits<T>::infinity() : -std::nan("1");
    }
}

template <typename T>
inline T fastLn(const T value, T (*log2func)(const T))
{
    const T g_rlog2_e = 0.693147180559945309417232121458176568075500134360255254120L;

    return g_rlog2_e * log2func(value);
}

template <typename T>
inline T fastLog10(const T value, T(*log2func)(const T))
{
    const T g_rlog2_10 = 0.301029995663981195213738894724493026768189881462108541310L;

    return g_rlog2_10 * log2func(value);
}

#endif /* SRC_LOGAPPROX_H_ */
