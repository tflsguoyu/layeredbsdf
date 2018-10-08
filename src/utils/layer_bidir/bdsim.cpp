#include <iostream>
#include <chrono>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/sampler.h>
#include "stats.h"

#if defined(MTS_OPENMP)
#   include <omp.h>
#endif

#define ANISOTROPIC_MEDIUM
//#define BIDIR_USE_ANALOG


MTS_NAMESPACE_BEGIN

class BidirSim : public Utility {
protected:
    // Estimates BRDF(wi, wo)
    Spectrum simulate1(const Vector &wi, const Vector &wo, int depth, const Medium *med, Sampler *sampler)
    {
        const PhaseFunction *phase = med->getPhaseFunction();
        Spectrum ret(0.0), throughput(1.0);

        Ray ray(Point(0.0, 0.0, 1.0), -wo, 0.0);
        for ( int i = 1; depth < 0 || i <= depth; ++i ) {
            MediumSamplingRecord mRec;
            if ( med->sampleDistance(ray, mRec, sampler) ) {
                ray.o.z += mRec.t*ray.d.z;
                if ( ray.o.z < 0.0 || ray.o.z > 1.0 ) break;

                // Next-event estimation
                throughput *= mRec.sigmaS*mRec.transmittance/mRec.pdfSuccess;

                Float phaseVal;
                if ( depth < 0 || i == depth ) {
                    phaseVal = phase->eval(PhaseFunctionSamplingRecord(mRec, -ray.d, wi));
                    ret += throughput*med->evalTransmittance(Ray(ray.o, wi, 0.0, (1.0 - ray.o.z)/wi.z, 0.0))*phaseVal;
                }

                PhaseFunctionSamplingRecord pRec(mRec, -ray.d);
                phaseVal = phase->sample(pRec, sampler);
                if ( phaseVal < Epsilon ) break;
                throughput *= phaseVal;
                ray.d = pRec.wo;
            }
            else
                break;
        }

        return ret/std::abs(wi.z);
    }


    // Estimates BRDF(wi, wo)
    Spectrum simulate1_bidir(const Vector &wi, const Vector &wo, const Medium *med, Sampler *sampler)
    {
        const PhaseFunction *phase = med->getPhaseFunction();

        typedef std::pair<Float, Vector> _INFO;
        std::vector<_INFO> path0, path1;
        std::vector<Spectrum> ratio0, ratio1, attn0, attn1, thru0, thru1;
        MediumSamplingRecord mRec;

        // "Camera" sub-path
        {
            Ray ray0(Point(0.0, 0.0, 1.0), -wo, 0.0);
            path0.push_back(_INFO(ray0.o.z, ray0.d));
            thru0.push_back(Spectrum(1.0));
            for ( int i = 0; ; ++i ) {
                if ( med->sampleDistance(ray0, mRec, sampler) ) {
                    ray0.o.z += mRec.t*ray0.d.z;
                    if ( ray0.o.z < 0.0 || ray0.o.z > 1.0 ) break;

                    Spectrum albedo = mRec.sigmaS*mRec.transmittance/mRec.pdfSuccess;
#ifdef BIDIR_USE_ANALOG
                    if ( sampler->next1D() < albedo.max() )
                        albedo /= albedo.max();
                    else
                        break;
#endif
                    thru0.push_back(thru0.back()*albedo);
                    Spectrum sigmaT = mRec.sigmaS + mRec.sigmaA;

                    PhaseFunctionSamplingRecord pRec(mRec, -ray0.d);
                    if ( phase->sample(pRec, sampler) < Epsilon ) break;
                    ray0.d = pRec.wo;

                    path0.push_back(_INFO(ray0.o.z, ray0.d));
                    attn0.push_back((-mRec.t*sigmaT).exp());

                    Spectrum r = Spectrum(1.0)/phase->eval(PhaseFunctionSamplingRecord(mRec, path0[i + 1].second, -path0[i].second));
                    if ( i ) r += Spectrum(1.0)/phase->eval(PhaseFunctionSamplingRecord(mRec, -path0[i - 1].second, path0[i].second));
                    r *= std::abs(path0[i].second.z)/attn0[i];
                    if ( i ) r += ratio0[i - 1];
                    ratio0.push_back(r);
                }
                else
                    break;
            }
        }

        // "Light" sub-path
        {
            Ray ray1(Point(0.0, 0.0, 1.0), -wi, 0.0);
            path1.push_back(_INFO(ray1.o.z, ray1.d));
            thru1.push_back(Spectrum(1.0));
            for ( int i = 0; ; ++i ) {
                if ( med->sampleDistance(ray1, mRec, sampler) ) {
                    ray1.o.z += mRec.t*ray1.d.z;
                    if ( ray1.o.z < 0.0 || ray1.o.z > 1.0 ) break;

                    Spectrum albedo = mRec.sigmaS*mRec.transmittance/mRec.pdfSuccess;
#ifdef BIDIR_USE_ANALOG
                    if ( sampler->next1D() < albedo.max() )
                        albedo /= albedo.max();
                    else
                        break;
#endif
                    thru1.push_back(thru1.back()*albedo);
                    Spectrum sigmaT = mRec.sigmaS + mRec.sigmaA;


                    PhaseFunctionSamplingRecord pRec(mRec, -ray1.d);
                    if ( phase->sample(pRec, sampler) < Epsilon ) break;
                    ray1.d = pRec.wo;

                    path1.push_back(_INFO(ray1.o.z, ray1.d));
                    attn1.push_back((-mRec.t*sigmaT).exp());

                    Spectrum r = Spectrum(1.0)/phase->eval(PhaseFunctionSamplingRecord(mRec, path1[i + 1].second, -path1[i].second));
                    if ( i ) r += Spectrum(1.0)/phase->eval(PhaseFunctionSamplingRecord(mRec, -path1[i - 1].second, path1[i].second));
                    r *= std::abs(path1[i].second.z)/attn1[i];
                    if ( i ) r += ratio1[i - 1];
                    ratio1.push_back(r);
                }
                else
                    break;
            }
        }

        Spectrum ret(0.0);
        for ( size_t i = 0; i < path0.size(); ++i ) {
            for ( size_t j = (i == 0 ? 1 : 0); j < path1.size(); ++j ) {
                Vector dirs[2];
                int ndirs = 0;
                if ( j ) dirs[ndirs++] = path0[i].second;
                if ( i ) dirs[ndirs++] = -path1[j].second;

                for ( int k = 0; k < ndirs; ++k ) {
                    const Vector &dir = dirs[k];

                    Float t = (path1[j].first - path0[i].first)/dir.z;
                    if ( t > Epsilon ) {
                        Spectrum r(0.0), r1;
                        if ( i ) r += Spectrum(1.0)/phase->eval(PhaseFunctionSamplingRecord(mRec, -path0[i - 1].second, dir));
                        if ( j ) r += Spectrum(1.0)/phase->eval(PhaseFunctionSamplingRecord(mRec, -path1[j - 1].second, -dir));
                        Spectrum tmp = med->evalTransmittance(Ray(Point(0, 0, path0[i].first), dir, 0, t, 0));
                        r *= std::abs(dir.z) / tmp;

                        if ( i ) {
                            r1 = Spectrum(1.0)/phase->eval(PhaseFunctionSamplingRecord(mRec, dir, -path0[i - 1].second));
                            if ( i > 1 ) {
                                r1 += Spectrum(1.0)/phase->eval(PhaseFunctionSamplingRecord(mRec, -path0[i - 2].second, path0[i - 1].second));
                                r += ratio0[i - 2];
                            }
                            r1 *= std::abs(path0[i - 1].second.z)/attn0[i - 1];
                            r += r1;
                        }

                        if ( j ) {
                            r1 = Spectrum(1.0)/phase->eval(PhaseFunctionSamplingRecord(mRec, -dir, -path1[j - 1].second));
                            if ( j > 1 ) {
                                r1 += Spectrum(1.0)/phase->eval(PhaseFunctionSamplingRecord(mRec, -path1[j - 2].second, path1[j - 1].second));
                                r += ratio1[j - 2];
                            }
                            r1 *= std::abs(path1[j - 1].second.z)/attn1[j - 1];
                            r += r1;
                        }

                        ret += thru0[i]*thru1[j]/r;
                    }
                }
            }
        }

        return ret;
    }


public:
    int run(int argc, char **argv) {
        int nworker;
#if defined(MTS_OPENMP)
        nworker = omp_get_max_threads();
#else
        nworker = 1;
#endif

        ref_vector<Sampler> samplers(nworker);
        samplers[0] = static_cast<Sampler *>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), Properties("independent")));
        for ( int i = 1; i < nworker; ++i ) samplers[i] = samplers[0]->clone();

        ref<Medium> mediumIso, mediumAniso;

        // Initialize homogeneous isotropic medium
        {
            Properties propMedium("homogeneous");
            propMedium.setSpectrum("sigmaT", Spectrum(2.0));
            propMedium.setSpectrum("albedo", Spectrum(0.9));

            Properties propPhase("hg");
            //propPhase.setFloat("g", 0.95);
            propPhase.setFloat("g", 0.8);
            // propPhase.setFloat("g", 0.0);

            ref<PhaseFunction> phase = static_cast<PhaseFunction *>(PluginManager::getInstance()->createObject(MTS_CLASS(PhaseFunction), propPhase));
            phase->configure();

            mediumIso = static_cast<Medium *>(PluginManager::getInstance()->createObject(MTS_CLASS(Medium), propMedium));
            mediumIso->addChild(phase);
            mediumIso->configure();
        }

        // Initialize homogeneous anisotropic medium
        {
            Properties propMicroflake("microflake");
            // propMicroflake.setFloat("stddev", 0.1);
            propMicroflake.setFloat("stddev", 0.25);
            ref<PhaseFunction> phase = static_cast<PhaseFunction *>(PluginManager::getInstance()->createObject(MTS_CLASS(PhaseFunction), propMicroflake));
            phase->configure();

            Properties propMedium("homogeneous_aniso");
            propMedium.setFloat("density", 2.0);
            propMedium.setSpectrum("albedo", Spectrum(0.9));
            propMedium.setVector("orientation", Vector(1.0, 0.0, 0.0));
            mediumAniso = static_cast<Medium *>(PluginManager::getInstance()->createObject(MTS_CLASS(Medium), propMedium));
            mediumAniso->addChild(phase);
            mediumAniso->configure();
        }

        int flag = (argc >= 2 ? atoi(argv[1]) : 0);
        ref<Medium> med = (flag ? mediumAniso : mediumIso);

        std::cout << med->toString() << std::endl;

        Vector d0(3.0, 2.0, 1.0), d1(1.0, 2.0, 3.0);
        d0 = normalize(d0);
        d1 = normalize(d1);

        long long N = (argc >= 3 ? atoll(argv[2]) : 1000000);
        int depth = (argc >= 4 ? atoi(argv[3]) : -1);
        int k = (argc >= 5 ? atoi(argv[4]) : 5);

        std::cout << '\n' << "N = " << N << ", depth = " << depth << ", k = " << k << '\n' << std::endl;

        std::vector<::Statistics> stats(nworker);
        std::chrono::time_point<std::chrono::steady_clock> start;

        // Uni-directional

        start = std::chrono::steady_clock::now();
        for ( int i = 0; i < nworker; ++i ) stats[i].reset();
#if defined(MTS_OPENMP)
#   pragma omp parallel for
#endif
        for ( long long omp_i = 0; omp_i < N; ++omp_i ) {
            int tid;
#if defined(MTS_OPENMP)
            tid = omp_get_thread_num();
#else
            tid = 0;
#endif

            Spectrum val = simulate1(d0, d1, depth, med, samplers[tid]);
            stats[tid].push(val[0]);
        }
        for ( int i = 1; i < nworker; ++i ) stats[0].push(stats[i]);

        printf("Unidir       : %.3le +- %.3le (%.1lf secs)\n", stats[0].getMean(), stats[0].getCI(),
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count()/1000.0);

        // Uni-directional (adjoint)

        start = std::chrono::steady_clock::now();
        for ( int i = 0; i < nworker; ++i ) stats[i].reset();
#if defined(MTS_OPENMP)
#   pragma omp parallel for
#endif
        for ( long long omp_i = 0; omp_i < N; ++omp_i ) {
            int tid;
#if defined(MTS_OPENMP)
            tid = omp_get_thread_num();
#else
            tid = 0;
#endif

            Spectrum val = simulate1(d1, d0, depth, med, samplers[tid]);
            stats[tid].push(val[0]);
        }
        for ( int i = 1; i < nworker; ++i ) stats[0].push(stats[i]);

        printf("Unidir (adj.): %.3le +- %.3le (%.1lf secs)\n", stats[0].getMean(), stats[0].getCI(),
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count()/1000.0);

        // Bidirectional

        if ( depth < 0 ) {
            start = std::chrono::steady_clock::now();
            for ( int i = 0; i < nworker; ++i ) stats[i].reset();
#if defined(MTS_OPENMP)
#   pragma omp parallel for
#endif
            for ( long long omp_i = 0; omp_i < N/k; ++omp_i ) {
                int tid;
#if defined(MTS_OPENMP)
                tid = omp_get_thread_num();
#else
                tid = 0;
#endif

                Spectrum val = simulate1_bidir(d0, d1, med, samplers[tid]);
                stats[tid].push(val[0]);
            }
            for ( int i = 1; i < nworker; ++i ) stats[0].push(stats[i]);

            printf("Bidir        : %.3le +- %.3le (%.1lf secs)\n", stats[0].getMean(), stats[0].getCI(),
                std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count()/1000.0);
        }

        return 0;
    }

    MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(BidirSim, "Bidirectional Simulation for Layered Config (Volume)")
MTS_NAMESPACE_END
