#include <iostream>
#include <chrono>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/sampler.h>
#include <boost/math/special_functions/fpclassify.hpp>
#include "stats.h"

//#undef MTS_OPENMP

#if defined(MTS_OPENMP)
#   include <omp.h>
#endif


MTS_NAMESPACE_BEGIN

class BidirSimSurf : public Utility {
protected:
    Spectrum simulate2(const Vector &wi, const Vector &wo, const BSDF *top, const BSDF *bottom, Sampler *sampler) {
        Spectrum ret(0.0), throughput(1.0);

        Intersection its;
        its.geoFrame = its.shFrame = Frame(Vector(0.0, 0.0, 1.0));
        its.hasUVPartials = false;
        its.uv = Point2(0.0);
        its.time = 0.0;

        BSDFSamplingRecord bRec(its, sampler, ERadiance);
        bRec.typeMask = BSDF::EAll;
        bRec.component = -1;

        Float z = 1.0;
        Vector d = -wo;
        for ( int i = 0; ; ++i ) {
            bRec.wi = -d;
            const BSDF *bsdf = ( z > 1.0 - Epsilon ? top : bottom );
            if ( bsdf == top ) {
                bRec.wo = wi;
                ret += throughput*bsdf->eval(bRec);
            }

            Spectrum val = bsdf->sample(bRec, sampler->next2D());
            if ( val.isZero() ) break;
            throughput *= val;
            d = bRec.wo;

            if ( d.z > Epsilon )
                z += 1.0;
            else
                z -= 1.0;
            if ( z < -Epsilon || z > 1.0 + Epsilon ) break;
        }

        return ret/std::abs(wi.z);
    }


    Spectrum simulate2_bidir(const Vector &wi, const Vector &wo, const BSDF *top, const BSDF *bottom, Sampler *sampler) {
        typedef std::pair<Float, Vector> _INFO;
        std::vector<_INFO> path0, path1;
        std::vector<Spectrum> thru0, thru1;
        std::vector<Float> ratio0, ratio1;

        Intersection its;
        its.geoFrame = its.shFrame = Frame(Vector(0.0, 0.0, 1.0));
        its.hasUVPartials = false;
        its.uv = Point2(0.0);
        its.time = 0.0;

        // "Camera" sub-path
        {
            BSDFSamplingRecord bRec(its, sampler, ERadiance);
            bRec.typeMask = BSDF::EAll;
            bRec.component = -1;

            Float z = 1.0;
            Vector d = -wo;
            thru0.push_back(Spectrum(1.0));
            for ( int i = 0; ; ++i ) {
                bRec.wi = -d;
                const BSDF *bsdf = ( z > 1.0 - Epsilon ? top : bottom );
                Spectrum val = bsdf->sample(bRec, sampler->next2D());
                if ( val.isZero() ) break;
                d = bRec.wo;

                thru0.push_back(thru0.back()*val);
                path0.push_back(_INFO(z, d));

                if ( d.z > Epsilon )
                    z += 1.0;
                else
                    z -= 1.0;
                if ( z < -Epsilon || z > 1.0 + Epsilon ) break;
            }

            ratio0.resize(path0.size());
            if ( !path0.empty() ) {
                BSDFSamplingRecord bRec0 = bRec, bRec1 = bRec;
                ratio0[0] = 0.0;
                for ( size_t i = 1; i < path0.size(); ++i ) {
                    Vector d0 = ( i > 1 ? path0[i - 2].second : -wo ),
                           d1 = path0[i - 1].second,
                           d2 = path0[i].second;
                    const BSDF *bsdf0 = ( path0[i - 1].first > 1.0 - Epsilon ? top : bottom ),
                               *bsdf1 = ( path0[i].first > 1.0 - Epsilon ? top : bottom );

                    bRec0.wi = -d0; bRec0.wo = d1;
                    bRec1.wi = d2; bRec1.wo = -d1;
                    Float r = 1.0/bsdf0->pdf(bRec0) + 1.0/bsdf1->pdf(bRec1);

                    bRec0.wi = -d1; bRec0.wo = d2;
                    Float r1 = bsdf1->pdf(bRec1)/bsdf1->pdf(bRec0);

                    ratio0[i] = r1*(r + ratio0[i - 1]);
                }
            }
        }

        // "Light" sub-path
        {
            BSDFSamplingRecord bRec(its, sampler, EImportance);
            bRec.typeMask = BSDF::EAll;
            bRec.component = -1;

            double z = 1.0;
            Vector d = -wi;
            thru1.push_back(Spectrum(1.0));
            for ( int i = 0; ; ++i ) {
                bRec.wi = -d;
                const BSDF *bsdf = ( z > 1.0 - Epsilon ? top : bottom );
                Spectrum val = bsdf->sample(bRec, sampler->next2D());
                if ( val.isZero() ) break;
                d = bRec.wo;

                thru1.push_back(thru1.back()*val);
                path1.push_back(_INFO(z, d));

                if ( d.z > Epsilon )
                    z += 1.0;
                else
                    z -= 1.0;
                if ( z < -Epsilon || z > 1.0 + Epsilon ) break;
            }

            ratio1.resize(path1.size());
            if ( !path1.empty() ) {
                BSDFSamplingRecord bRec0(bRec), bRec1(bRec);
                ratio1[0] = 0.0;
                for ( size_t i = 1; i < path1.size(); ++i ) {
                    Vector d0 = ( i > 1 ? path1[i - 2].second : -wi ),
                           d1 = path1[i - 1].second,
                           d2 = path1[i].second;
                    const BSDF *bsdf0 = ( path1[i - 1].first > 1.0 - Epsilon ? top : bottom ),
                               *bsdf1 = ( path1[i].first > 1.0 - Epsilon ? top : bottom );

                    bRec0.wi = -d0; bRec0.wo = d1;
                    bRec1.wi = d2; bRec1.wo = -d1;
                    Float r = 1.0/bsdf0->pdf(bRec0) + 1.0/bsdf1->pdf(bRec1);

                    bRec0.wi = -d1; bRec0.wo = d2;
                    Float r1 = bsdf1->pdf(bRec1)/bsdf1->pdf(bRec0);

                    ratio1[i] = r1*(r + ratio1[i - 1]);
                }
            }
        }

        BSDFSamplingRecord bRec0(its, sampler, ERadiance), bRec1(its, sampler, EImportance);
        bRec0.typeMask = BSDF::EAll;
        bRec0.component = -1;
        bRec1.typeMask = BSDF::EAll;
        bRec1.component = -1;

        bRec0.wi = wo; bRec0.wo = wi;
        Spectrum ret = top->eval(bRec0)/std::abs(wi.z);

        for ( size_t i = 0; i < path0.size(); ++i )
            for ( size_t j = 0; j < path1.size(); ++j ) {
                if ( std::abs(path0[i].first - path1[j].first) < Epsilon ) continue;

                const BSDF *bsdf0, *bsdf1, *bsdf2;
                Spectrum f, bsdfVal;
                Float w, w1;

                bsdf0 = ( path0[i].first > 1.0 - Epsilon ? top : bottom );
                bsdf1 = ( path1[j].first > 1.0 - Epsilon ? top : bottom );

                // "i" path

                bRec0.wi = ( i ? -path0[i - 1].second : wo );
                bRec0.wo = path0[i].second;

                bRec1.wi = ( j ? -path1[j - 1].second : wi );
                bRec1.wo = -path0[i].second;

                bsdfVal = bsdf1->eval(bRec1)/std::abs(bRec1.wo.z);
                if ( !bsdfVal.isZero() ) {
                    f = thru0[i + 1]*thru1[j]*bsdfVal;
                    w = 1.0 + bsdf1->pdf(bRec1)/bsdf0->pdf(bRec0);

                    if ( i ) {
                        bsdf2 = ( path0[i - 1].first > 1.0 - Epsilon ? top : bottom );

                        BSDFSamplingRecord bRec2(bRec0);
                        bRec2.reverse();

                        BSDFSamplingRecord bRec3(bRec0);
                        bRec3.wi = ( i > 1 ? -path0[i - 2].second : wo );
                        bRec3.wo = path0[i - 1].second;

                        w1 = 1.0 + bsdf0->pdf(bRec2)/bsdf2->pdf(bRec3);
                        w += w1*bsdf1->pdf(bRec1)/bsdf0->pdf(bRec0);

                        w += ratio0[i - 1]*bsdf0->pdf(bRec2)*bsdf1->pdf(bRec1)/bsdf0->pdf(bRec0);
                    }
                    if ( j ) {
                        bsdf2 = ( path1[j - 1].first > 1.0 - Epsilon ? top : bottom );

                        BSDFSamplingRecord bRec2(bRec1);
                        bRec2.reverse();

                        BSDFSamplingRecord bRec3(bRec1);
                        bRec3.wi = ( j > 1 ? -path1[j - 2].second : wi );
                        bRec3.wo = path1[j - 1].second;

                        w1 = 1.0 + bsdf1->pdf(bRec2)/bsdf2->pdf(bRec3);
                        w += w1;

                        w += ratio1[j - 1]*bsdf1->pdf(bRec2);
                    }
                    ret += f/w;
                }

                // "j" path

                bRec0.wi = ( i ? -path0[i - 1].second : wo );
                bRec0.wo = -path1[j].second;

                bRec1.wi = ( j ? -path1[j - 1].second : wi );
                bRec1.wo = path1[j].second;

                bsdfVal = bsdf0->eval(bRec0)/std::abs(bRec0.wo.z);
                if ( !bsdfVal.isZero() ) {
                    f = thru0[i]*thru1[j + 1]*bsdfVal;
                    w = 1.0 + bsdf0->pdf(bRec0)/bsdf1->pdf(bRec1);
                    if ( i ) {
                        bsdf2 = ( path0[i - 1].first > 1.0 - Epsilon ? top : bottom );

                        BSDFSamplingRecord bRec2(bRec0);
                        bRec2.reverse();

                        BSDFSamplingRecord bRec3(bRec0);
                        bRec3.wi = ( i > 1 ? -path0[i - 2].second : wo );
                        bRec3.wo = path0[i - 1].second;

                        w1 = 1.0 + bsdf0->pdf(bRec2)/bsdf2->pdf(bRec3);
                        w += w1;

                        w += ratio0[i - 1]*bsdf0->pdf(bRec2);
                    }
                    if ( j ) {
                        bsdf2 = ( path1[j - 1].first > 1.0 - Epsilon ? top : bottom );

                        BSDFSamplingRecord bRec2(bRec1);
                        bRec2.reverse();

                        BSDFSamplingRecord bRec3(bRec1);
                        bRec3.wi = ( j > 1 ? -path1[j - 2].second : wi );
                        bRec3.wo = path1[j - 1].second;

                        w1 = 1.0 + bsdf1->pdf(bRec2)/bsdf2->pdf(bRec3);
                        w += w1*bsdf0->pdf(bRec0)/bsdf1->pdf(bRec1);

                        w += ratio1[j - 1]*bsdf0->pdf(bRec0)*bsdf1->pdf(bRec2)/bsdf1->pdf(bRec1);
                    }
                    ret += f/w;
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

        ref<BSDF> bsdfTop, bsdfBottom;
        {
            Properties propTop("roughdielectric");
            propTop.setFloat("intIOR", 1.5);
            propTop.setFloat("alpha", 0.25);

            bsdfTop = static_cast<BSDF *>(PluginManager::getInstance()->createObject(MTS_CLASS(BSDF), propTop));
            bsdfTop->configure();

            Properties propBottom("roughconductor");
            propBottom.setFloat("alpha", 0.25);
            // Properties propBottom("diffuse");

            bsdfBottom = static_cast<BSDF *>(PluginManager::getInstance()->createObject(MTS_CLASS(BSDF), propBottom));
            bsdfBottom->configure();
        }

        Vector d0(3.0, 2.0, 1.0), d1(1.0, 2.0, 3.0);
        d0 = normalize(d0);
        d1 = normalize(d1);

        long long N = ( argc >=2 ? atoll(argv[1]) : 1000000 );

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

            Spectrum val = simulate2(d0, d1, bsdfTop, bsdfBottom, samplers[tid]);
            stats[tid].push(val[0]);
        }
        for ( int i = 1; i < nworker; ++i ) stats[0].push(stats[i]);

        printf("Unidir       : %.2le +- %.2le (%.1lf secs)\n", stats[0].getMean(), stats[0].getCI(),
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count()/1000.0);

#ifdef RUN_ADJOINT
        // Uni-directional, adjoint

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

            Spectrum val = simulate2(d1, d0, bsdfTop, bsdfBottom, samplers[tid]);
            stats[tid].push(val[0]);
        }
        for ( int i = 1; i < nworker; ++i ) stats[0].push(stats[i]);

        printf("Unidir (adj.): %.2le +- %.2le (%.1lf secs)\n", stats[0].getMean(), stats[0].getCI(),
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count()/1000.0);
#endif

        // Bi-directional

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

            Spectrum val = simulate2_bidir(d0, d1, bsdfTop, bsdfBottom, samplers[tid]);
            stats[tid].push(val[0]);
        }
        for ( int i = 1; i < nworker; ++i ) stats[0].push(stats[i]);

        printf("Bidir        : %.2le +- %.2le (%.1lf secs)\n", stats[0].getMean(), stats[0].getCI(),
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count()/1000.0);

        return 0;
    }

    MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(BidirSimSurf, "Bidirectional Simulation for Layered Config (Surface)")
MTS_NAMESPACE_END
