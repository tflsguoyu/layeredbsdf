#include <iostream>
#include <chrono>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/render/sampler.h>
#include "stats.h"
#include "sgmixture.hpp"

#define NLOBES 3


MTS_NAMESPACE_BEGIN


class SphericalGaussianBench : public Utility {
public:
    int run(int argc, char **argv) {
        const Float gVal = 0.8, params[] = { 1.0759e-01, 2.0330e+00, 1.0008e+00, 1.4336e+01, 2.3987e+00, 6.5918e+01 };
        const Float params1[] = { 1.0452, 1.7461 };

        Properties propPhase("hg");
        propPhase.setFloat("g", gVal);
        ref<PhaseFunction> phase = static_cast<PhaseFunction *>(PluginManager::getInstance()->createObject(MTS_CLASS(PhaseFunction), propPhase));
        phase->configure();

        ref<Sampler> sampler = static_cast<Sampler *>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), Properties("independent")));
        PhaseFunctionSamplingRecord pRec(MediumSamplingRecord(), Vector(0.0, 0.0, 0.0));

        Vector e0 = normalize(Vector(-0.5, 0.0, 1.0)), e1 = normalize(Vector(0.5, 0.0, 1.0));
        long long N = 10000000, K = 3;

        Float dZ;
        if ( argc >= 2 )
            dZ = atof(argv[1]);
        else
            dZ = 0.5f;
        Log(EInfo, "dZ = %.2lf", dZ);

        int flag;
        if ( argc >= 3 )
            flag = atoi(argv[2]);
        else
            flag = 0;
        Log(EInfo, "flag = %d", flag);

        std::chrono::time_point<std::chrono::steady_clock> start;
        ::Statistics stat;

        start = std::chrono::steady_clock::now();
        stat.reset();
        for ( long long i = 0; i < N; ++i ) {
            pRec.wi = -e0;
            phase->sample(pRec, sampler);

            Float t = dZ/pRec.wo.z;
            if ( t > Epsilon ) {
                Float val = math::fastexp(-t)/std::abs(pRec.wo.z);
                pRec.wi = -e1;
                val *= phase->eval(pRec);
                stat.push(val);
            }
            else
                stat.push(0.0f);
        }
        printf("%.2le +- %.2le (%.2lf secs)\n", stat.getMean(), stat.getCI(),
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count()/1000.0);

        SphericalGaussianMixture<NLOBES> sgm0(params), sgm1(params);
        sgm0.setMu(e0); sgm1.setMu(e1);

        SphericalGaussianMixture<1> sgm2(params1);
        sgm2.setMu(Vector(0.0f, 0.0f, dZ/std::abs(dZ)));

        SphericalGaussianMixture<NLOBES*NLOBES> sgm3(sgm0, sgm1), sgm;
        if ( flag )
            sgm.init(sgm3, sgm2);
        else
            sgm = sgm3;

        start = std::chrono::steady_clock::now();
        stat.reset();
        for ( long long i = 0; i < N/K; ++i ) {
            Point2 sample = sampler->next2D();
            Float pdf;
            Vector x = sgm.sample(sampler->next2D(), pdf);
            //if ( std::abs(pdf - sgm.pdf(x)) > Epsilon ) fprintf(stderr, "WTF??\n");

            Float t = dZ/x.z;
            if ( t > Epsilon ) {
                Float val = math::fastexp(-t)/std::abs(x.z);
                pRec.wi = -e0; pRec.wo = x;
                val *= phase->eval(pRec);
                pRec.wi = -e1; pRec.wo = x;
                val *= phase->eval(pRec);
                stat.push(val/pdf);
            }
            else
                stat.push(0.0f);
        }
        printf("%.2le +- %.2le (%.2lf secs)\n", stat.getMean(), stat.getCI(),
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count()/1000.0);

        return 0;
    }

    MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(SphericalGaussianBench, "SG Test")
MTS_NAMESPACE_END
