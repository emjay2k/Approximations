/*
 * logtest.h
 *
 *  Created on: Jul 3, 2019
 *      Author: Matthias Jung
 */

#ifndef SRC_LOGTEST_H_
#define SRC_LOGTEST_H_

#include <chrono>
#include <cmath>
#include <future>
#include <iostream>
#include <vector>
#include <thread>
#include "logapprox.h"

struct MaxErrorFuture {
    double dMaxError[7] = { 0 };
};

template <typename T>
static void validateWorkerP5(size_t uNumSamples, size_t uStart, size_t uEnd, std::promise<MaxErrorFuture> *promiseObj) {
    T dStepSize = 1.0/static_cast<T>(uNumSamples);

    MaxErrorFuture retVal;
    double *dMaxError = &retVal.dMaxError[0];
    // generate data for mantissa only (since exponent is perfect)
    // actual representation is from [1,2], so we sample that range
    for(size_t i = uStart; i < uEnd; ++i) {
        T x = 1.0 + static_cast<T>(i) * dStepSize;
        T dPrecise = log2(x);
        T dApprox = fastLog2p5<T>(x);

        double dDiff = static_cast<double>(fabs(dPrecise - dApprox));
        if(dDiff > dMaxError[0])
        {
            dMaxError[0] = dDiff;
        }
    }

    promiseObj->set_value(std::move(retVal));
}

template <typename T>
static void validateWorkerP6(size_t uNumSamples, size_t uStart, size_t uEnd, std::promise<MaxErrorFuture> *promiseObj) {
    T dStepSize = 1.0/static_cast<T>(uNumSamples);

    MaxErrorFuture retVal;
    double *dMaxError = &retVal.dMaxError[0];
    // generate data for mantissa only (since exponent is perfect)
    // actual representation is from [1,2], so we sample that range
    for(size_t i = uStart; i < uEnd; ++i) {
        T x = 1.0 + static_cast<T>(i) * dStepSize;
        T dPrecise = log2(x);
        T dApprox = fastLog2p6<T>(x);

        double dDiff = static_cast<double>(fabs(dPrecise - dApprox));
        if(dDiff > dMaxError[0])
        {
            dMaxError[0] = dDiff;
        }
    }

    promiseObj->set_value(std::move(retVal));
}

template <typename T>
static void validateWorker(size_t uNumSamples, size_t uStart, size_t uEnd, std::promise<MaxErrorFuture> *promiseObj) {
    T dStepSize = 1.0/static_cast<T>(uNumSamples);

    MaxErrorFuture retVal;
    double *dMaxError = &retVal.dMaxError[0];
    // generate data for mantissa only (since exponent is perfect)
    // actual representation is from [1,2], so we sample that range
    for(size_t i = uStart; i < uEnd; ++i) {
        T x = 1.0 + static_cast<T>(i) * dStepSize;
        T dPrecise = log2(x);
        T dApprox = fastLog2p1<T>(x);

        double dDiff = static_cast<double>(fabs(dPrecise - dApprox));
        if(dDiff > dMaxError[0])
        {
            dMaxError[0] = dDiff;
        }
        dApprox = fastLog2p2<T>(x);
        dDiff = fabs(dPrecise - dApprox);
        if(dDiff > dMaxError[1])
        {
            dMaxError[1] = dDiff;
        }
        dApprox = fastLog2p3<T>(x);
        dDiff = fabs(dPrecise - dApprox);
        if(dDiff > dMaxError[2])
        {
            dMaxError[2] = dDiff;
        }
        dApprox = fastLog2p4<T>(x);
        dDiff = fabs(dPrecise - dApprox);
        if(dDiff > dMaxError[3])
        {
            dMaxError[3] = dDiff;
        }
        dApprox = fastLog2p5<T>(x);
        dDiff = fabs(dPrecise - dApprox);
        if(dDiff > dMaxError[4])
        {
            dMaxError[4] = dDiff;
        }
        dApprox = fastLog2p6<T>(x);
        dDiff = fabs(dPrecise - dApprox);
        if(dDiff > dMaxError[5])
        {
            dMaxError[5] = dDiff;
        }
//        dApprox = log2f(x);
//        dDiff = fabs(dPrecise - dApprox);
//        if(dDiff > dMaxError[6])
//        {
//            dMaxError[6] = dDiff;
//        }
    }

    promiseObj->set_value(std::move(retVal));
}

template <typename T>
static inline void validateAccuracy(size_t uNumSamples, const size_t uNumThreads=1)
{
    // threads, promises and futures
    std::vector<std::unique_ptr<std::thread>> threadList;
    threadList.reserve(uNumThreads);
    std::vector<std::promise<MaxErrorFuture>> promises(uNumThreads);
    std::vector<std::future<MaxErrorFuture>> futures;
    futures.reserve(uNumThreads);

    // start the workers
    size_t uNumSamplesPerThread = uNumSamples/uNumThreads;
    size_t uStart = 0;
    size_t uEnd = std::min(uStart + uNumSamplesPerThread, uNumSamples);
    for(size_t i = 0; i < uNumThreads; ++i) {
        futures.emplace_back(promises[i].get_future());
        threadList.emplace_back(std::make_unique<std::thread>(validateWorker<T>, uNumSamples, uStart, uEnd, &promises[i]));

        uStart = uEnd;
        uEnd = (i == (uNumThreads - 2)) ? uNumSamples+1 : (uStart + uNumSamplesPerThread);
    }

    // get the futures when threads are done
    double dMaxError[7] = { 0 };
    for(auto &&f : futures) {
        auto tf = f.get();
        for(int i = 0; i < 7; ++i) {
            dMaxError[i] = std::max(dMaxError[i], tf.dMaxError[i]);
        }
    }
    // join all threads
    for(auto &&t : threadList) {
        t->join();
    }

    std::cout.precision(24);
    std::cout << "Max errors: ";
    for(int i = 0; i < 7; ++i) {
        std::cout << dMaxError[i] << ",";
    }
    std::cout << std::endl;
}

template <typename T>
static inline void validatePerformance(size_t uNumSamples) {
    std::vector<T> values(uNumSamples);

    // generate data
    for(size_t i = 1; i < uNumSamples; ++i) {
        values[i] = static_cast<T>(i);
    }

    // compare speed
    T dSum = 0;
    auto start_time = std::chrono::high_resolution_clock::now();
    for(size_t i = 1; i < uNumSamples; ++i) {
        dSum += log2(values[i]);
    }
    auto end_time1 = std::chrono::high_resolution_clock::now();
    for(size_t i = 1; i < uNumSamples; ++i) {
        dSum += fastLog2p1<T>(values[i]);
    }
    auto end_time2 = std::chrono::high_resolution_clock::now();
    for(size_t i = 1; i < uNumSamples; ++i) {
        dSum += fastLog2p2<T>(values[i]);
    }
    auto end_time3 = std::chrono::high_resolution_clock::now();
    for(size_t i = 1; i < uNumSamples; ++i) {
        dSum += fastLog2p3<T>(values[i]);
    }
    auto end_time4 = std::chrono::high_resolution_clock::now();
    for(size_t i = 1; i < uNumSamples; ++i) {
        dSum += fastLog2p4<T>(values[i]);
    }
    auto end_time5 = std::chrono::high_resolution_clock::now();
    for(size_t i = 1; i < uNumSamples; ++i) {
        dSum += fastLog2p5<T>(values[i]);
    }
    auto end_time6 = std::chrono::high_resolution_clock::now();
    for(size_t i = 1; i < uNumSamples; ++i) {
        dSum += fastLog2p6<T>(values[i]);
    }
    auto end_time7 = std::chrono::high_resolution_clock::now();

    // compare to logf
    // generate new data so the static_cast from T to float is removed
    values.clear();
    std::vector<float> valuesFloat(uNumSamples);
    for(size_t i = 1; i < uNumSamples; ++i) {
        valuesFloat[i] = static_cast<T>(i);
    }
    auto start_time8 = std::chrono::high_resolution_clock::now();
    float dSumFloat = 0;
    for(size_t i = 1; i < uNumSamples; ++i) {
        dSumFloat += log2f(values[i]);
    }
    auto end_time8 = std::chrono::high_resolution_clock::now();

    std::cout << dSum << "," << dSumFloat << std::endl;
    std::cout << "speed:" << std::chrono::duration_cast<std::chrono::microseconds>(end_time1-start_time).count()
                    << "," << std::chrono::duration_cast<std::chrono::microseconds>(end_time2-end_time1).count()
                    << "," << std::chrono::duration_cast<std::chrono::microseconds>(end_time3-end_time2).count()
                    << "," << std::chrono::duration_cast<std::chrono::microseconds>(end_time4-end_time3).count()
                    << "," << std::chrono::duration_cast<std::chrono::microseconds>(end_time5-end_time4).count()
                    << "," << std::chrono::duration_cast<std::chrono::microseconds>(end_time6-end_time5).count()
                    << "," << std::chrono::duration_cast<std::chrono::microseconds>(end_time7-end_time6).count()
                    << "," << std::chrono::duration_cast<std::chrono::microseconds>(end_time8-start_time8).count() << std::endl;
}

#endif /* SRC_LOGTEST_H_ */
