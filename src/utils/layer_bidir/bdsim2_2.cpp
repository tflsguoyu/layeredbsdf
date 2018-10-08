#include <iostream>
#include <chrono>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/core/frame.h>
#include <mitsuba/core/pmf.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/sampler.h>
#include <boost/math/special_functions/fpclassify.hpp>
#include "stats.h"
#include "sgmixture.hpp"

//#undef MTS_OPENMP

#if defined(MTS_OPENMP)
#   include <omp.h>
#endif

#define BIDIR_USE_ANALOG
#define NLOBES 3


MTS_NAMESPACE_BEGIN


struct PathInfo {
    Float x;
    Vector d;
    const BSDF *surf;
    Spectrum thru0, thru1;
#ifdef BIDIR_USE_ANALOG
    Float pSurvival;
#endif
};


class BidirSimFullEx : public Utility {
protected:
    void simulatePath(const Vector &w0, ETransportMode mode,
                      const BSDF *top, const BSDF *bottom, const Medium *med, Sampler *sampler,
                      std::vector<PathInfo> &path, std::vector<Vector2> &vtxProb, std::vector<Vector2> &edgeProb)
    {
        const PhaseFunction *phase = med->getPhaseFunction();
        Spectrum throughput(1.0);

        Intersection its;
        its.geoFrame = its.shFrame = Frame(Vector(0.0, 0.0, 1.0));
        its.hasUVPartials = false;
        its.uv = Point2(0.0);
        its.time = 0.0;

        BSDFSamplingRecord bRec(its, sampler, mode);
        bRec.typeMask = BSDF::EAll;
        bRec.component = -1;

        Ray ray(Point(0.0, 0.0, 1.0), -w0, 0.0);
        bool surface = true;
        for ( int i = 0; ; ++i ) {
            PathInfo info;
            Vector2 prob;

            if ( surface ) {
                const BSDF *bsdf = ( ray.o.z > 1.0 - Epsilon ? top : bottom );
                bRec.wi = -ray.d;
                Spectrum bsdfVal = bsdf->sample(bRec, sampler->next2D());
                if ( bsdfVal.isZero() ) {
                    edgeProb.pop_back(); break;
                }

                info.surf = bsdf;
                info.x = ray.o.z;
                info.thru0 = throughput;
                info.d = ray.d = bRec.wo;
                info.thru1 = throughput *= bsdfVal;
#ifdef BIDIR_USE_ANALOG
                info.pSurvival = 1.0;
#endif
                path.push_back(info);

                prob[0] = bsdf->pdf(bRec);
                bRec.reverse();
                prob[1] = bsdf->pdf(bRec);
                if ( prob[0] < Epsilon || prob[1] < Epsilon )
                    fprintf(stderr, "Badness: zero prob (surface)");
                vtxProb.push_back(prob);

                if ( ray.o.z > 1.0 - Epsilon && ray.d.z > Epsilon || ray.o.z < Epsilon && ray.d.z < Epsilon ) break;
                surface = false;
            }
            else {
                MediumSamplingRecord mRec;
                Float maxt = ((ray.d.z > 0.0 ? 1.0 : 0.0) - ray.o.z)/ray.d.z;
                if ( med->sampleDistance(Ray(ray, 0, maxt), mRec, sampler) ) {
                    info.surf = NULL;
                    info.x = ray.o.z += mRec.t*ray.d.z;
                    if ( ray.o.z < 0.0 || ray.o.z > 1.0 ) {
                        fprintf(stderr, "Badness 0: %.4lf\n", ray.o.z);
                    }

                    Spectrum albedo = mRec.sigmaS*mRec.transmittance/mRec.pdfSuccess;
#ifdef BIDIR_USE_ANALOG
                    Float pSurvival = albedo.max();
                    if ( sampler->next1D() > pSurvival ) break;
                    info.thru0 = throughput *= albedo/pSurvival;
                    info.pSurvival = pSurvival;
#else
                    info.thru0 = throughput *= albedo;
#endif

                    prob[0] = mRec.pdfSuccess/std::abs(ray.d.z);
                    prob[1] = path.back().surf ? mRec.pdfFailure : mRec.pdfSuccessRev/std::abs(ray.d.z);
                    edgeProb.push_back(prob);

                    PhaseFunctionSamplingRecord pRec(mRec, -ray.d);
                    Float phaseVal = phase->sample(pRec, sampler);
                    if ( phaseVal < Epsilon ) {
                        edgeProb.pop_back(); break;
                    }
                    info.d = ray.d = pRec.wo;
                    info.thru1 = throughput *= phaseVal;
                    path.push_back(info);

                    prob[0] = phase->pdf(pRec);
                    pRec.reverse();
                    prob[1] = phase->pdf(pRec);
#ifdef BIDIR_USE_ANALOG
                    prob *= pSurvival;
#endif
                    if ( prob[0] < Epsilon || prob[1] < Epsilon )
                        fprintf(stderr, "Badness: zero prob (volume)");
                    vtxProb.push_back(prob);
                }
                else {
                    throughput *= mRec.transmittance/mRec.pdfFailure;
                    ray.o.z = ray.d.z > 0.0 ? 1.0 : 0.0;

                    prob[0] = mRec.pdfFailure;
                    prob[1] = path.back().surf ? mRec.pdfFailure : mRec.pdfSuccessRev/std::abs(ray.d.z);
                    edgeProb.push_back(prob);

                    surface = true;
                }
            }
        }

        if ( !path.empty() ) {
            if ( path.size() != vtxProb.size() || path.size() != edgeProb.size() + 1 )
                fprintf(stderr, "Badness 1: %zu %zu %zu\n", path.size(), vtxProb.size(), edgeProb.size());
        }
    }


    Spectrum simulate_bidir(const Vector &wi, const Vector &wo, const BSDF *top, const BSDF *bottom, const Medium *med,
        bool flag, const Float *params, Sampler *sampler)
    {
        const PhaseFunction *phase = med->getPhaseFunction();

        std::vector<PathInfo> path0, path1;
        std::vector<Vector2> edgeProb0, edgeProb1, vtxProb0, vtxProb1;
        std::vector<Float> ratio0, ratio1;

        Intersection its;
        its.geoFrame = its.shFrame = Frame(Vector(0.0, 0.0, 1.0));
        its.hasUVPartials = false;
        its.uv = Point2(0.0);
        its.time = 0.0;

        SphericalGaussianMixture<NLOBES> sgm0(params), sgm1(params);
        SphericalGaussianMixture<NLOBES*NLOBES> sgm;

        // "Camera" sub-path
        {
            simulatePath(wo, ERadiance, top, bottom, med, sampler, path0, vtxProb0, edgeProb0);

            ratio0.resize(path0.size());
            if ( !path0.empty() ) {
                ratio0[0] = 0.0;
                for ( size_t i = 1; i < path0.size(); ++i ) {
                    if ( flag && path0[i - 1].surf == NULL && path0[i].surf == NULL ) {
                        sgm0.setMu(path0[i - 2].d);
                        sgm1.setMu(path0[i].d);
                        sgm.init(sgm0, sgm1);
                        ratio0[i] = sgm.pdf(path0[i - 1].d)/(vtxProb0[i - 1][0]*vtxProb0[i][0]*edgeProb0[i - 1][0]) +
                                    ratio0[i - 1]*vtxProb0[i][1]*edgeProb0[i - 1][1]/(vtxProb0[i][0]*edgeProb0[i - 1][0]);
                    }
                    else {
                        Float r = 1.0/vtxProb0[i - 1][0] + 1.0/vtxProb0[i][1];
                        Float r1 = vtxProb0[i][1]/vtxProb0[i][0];
                        ratio0[i] = r1*(r + edgeProb0[i - 1][1]*ratio0[i - 1])/edgeProb0[i - 1][0];
                    }
                    if ( !boost::math::isfinite(ratio0[i]) ) {
                        fprintf(stderr, "Badness ratio0: %zu %lf\n", i, ratio0[i]);
                    }
                }
            }
        }

        // "Light" sub-path
        {
            simulatePath(wi, EImportance, top, bottom, med, sampler, path1, vtxProb1, edgeProb1);

            ratio1.resize(path1.size());
            if ( !path1.empty() ) {
                ratio1[0] = 0.0;
                for ( size_t i = 1; i < path1.size(); ++i ) {
                    if ( flag && path1[i - 1].surf == NULL && path1[i].surf == NULL ) {
                        sgm0.setMu(path1[i - 2].d);
                        sgm1.setMu(path1[i].d);
                        sgm.init(sgm0, sgm1);
                        ratio1[i] = sgm.pdf(path1[i - 1].d)/(vtxProb1[i - 1][0]*vtxProb1[i][0]*edgeProb1[i - 1][0]) +
                                    ratio1[i - 1]*vtxProb1[i][1]*edgeProb1[i - 1][1]/(vtxProb1[i][0]*edgeProb1[i - 1][0]);
                    }
                    else {
                        Float r = 1.0/vtxProb1[i - 1][0] + 1.0/vtxProb1[i][1];
                        Float r1 = vtxProb1[i][1]/vtxProb1[i][0];
                        ratio1[i] = r1*(r + edgeProb1[i - 1][1]*ratio1[i - 1])/edgeProb1[i - 1][0];
                    }
                    if ( !boost::math::isfinite(ratio1[i]) ) {
                        fprintf(stderr, "Badness ratio1: %zu %lf\n", i, ratio1[i]);
                    }
                }
            }
        }

        BSDFSamplingRecord bRec0(its, sampler, ERadiance), bRec1(its, sampler, EImportance);
        bRec0.typeMask = BSDF::EAll;
        bRec0.component = -1;
        bRec1.typeMask = BSDF::EAll;
        bRec1.component = -1;

        MediumSamplingRecord mRec;
        mRec.medium = med;
        PhaseFunctionSamplingRecord pRec(mRec, Vector(0.0));

        bRec0.wi = wo; bRec0.wo = wi;
        Spectrum ret = top->eval(bRec0)/std::abs(wi.z);
        for ( size_t i = 0; i < path0.size(); ++i )
            for ( size_t j = 0; j < path1.size(); ++j ) {
                Spectrum f, funcVal;
                Float w, w1, t;
                Vector2 pdfI, pdfJ;

                if ( flag && path0[i].surf == NULL && path1[j].surf == NULL ) {
                    sgm0.setMu(path0[i - 1].d);
                    sgm1.setMu(-path1[j - 1].d);
                    sgm.init(sgm0, sgm1);
                    Float pdf;
                    Vector dir = sgm.sample(sampler->next2D(), pdf);
                    t = (path1[j].x - path0[i].x)/dir.z;
                    if ( t > Epsilon ) {
                        pRec.wi = -path1[j - 1].d;
                        pRec.wo = -dir;
                        pdfJ[0] = phase->eval(pRec);
                        pRec.reverse();
                        pdfJ[1] = phase->eval(pRec);

                        pRec.wi = -path0[i - 1].d;
                        pRec.wo = dir;
                        pdfI[0] = phase->eval(pRec);
                        pRec.reverse();
                        pdfI[1] = phase->eval(pRec);

                        funcVal = Spectrum(pdfI[0]*pdfJ[0]/pdf);
#ifdef BIDIR_USE_ANALOG
                        pdfJ *= path1[j].pSurvival;
                        pdfI *= path0[i].pSurvival;
#endif

                        if ( !funcVal.isZero() ) {
                            med->eval(Ray(Point(0.0, 0.0, path0[i].x), dir, 0.0, t, 0.0), mRec);
                            mRec.pdfSuccess /= std::abs(dir.z);
                            mRec.pdfSuccessRev /= std::abs(dir.z);

                            f = path0[i].thru0*path1[j].thru0*funcVal*mRec.transmittance/std::abs(dir.z);
                            w = 1.0;

                            if ( i ) {
                                if ( path0[i - 1].surf == NULL ) {
                                    sgm0.setMu(path0[i - 2].d);
                                    sgm1.setMu(dir);
                                    sgm.init(sgm0, sgm1);
                                    w1 = sgm.pdf(path0[i - 1].d);
                                }
                                else
                                    w1 = pdfI[1] + vtxProb0[i - 1][0];
                                w += pdfJ[0]*mRec.pdfSuccessRev*w1/(vtxProb0[i - 1][0]*edgeProb0[i - 1][0]*pdf) +
                                     ratio0[i - 1]*pdfJ[0]*mRec.pdfSuccessRev*pdfI[1]*edgeProb0[i - 1][1]/(pdf*edgeProb0[i - 1][0]);
                            }

                            if ( j ) {
                                if ( path1[j - 1].surf == NULL ) {
                                    sgm0.setMu(path1[j - 2].d);
                                    sgm1.setMu(-dir);
                                    sgm.init(sgm0, sgm1);
                                    w1 = sgm.pdf(path1[j - 1].d);
                                }
                                else
                                    w1 = pdfJ[1] + vtxProb1[j - 1][0];
                                w += pdfI[0]*mRec.pdfSuccess*w1/(vtxProb1[j - 1][0]*edgeProb1[j - 1][0]*pdf) +
                                     ratio1[j - 1]*pdfI[0]*mRec.pdfSuccess*pdfJ[1]*edgeProb1[j - 1][1]/(pdf*edgeProb1[j - 1][0]);
                            }

                            if ( !boost::math::isnormal(w) ) {
                                fprintf(stderr, "WTF (i-j): %lf\n", w);
                            }
                            ret += f/w;
                        }
                    }
                }
                else {
                    // "i" path

                    t = (path1[j].x - path0[i].x)/path0[i].d.z;
                    if ( t > Epsilon ) {
						med->eval(Ray(Point(0.0, 0.0, path0[i].x), path0[i].d, 0.0, t, 0.0), mRec);
						mRec.pdfSuccess /= std::abs(path0[i].d.z);
						mRec.pdfSuccessRev /= std::abs(path0[i].d.z);

                        if ( path1[j].surf ) {
                            bRec1.wi = ( j ? -path1[j - 1].d : wi );
                            bRec1.wo = -path0[i].d;
                            funcVal = path1[j].surf->eval(bRec1);
                        }
                        else {
                            pRec.wi = -path1[j - 1].d;
                            pRec.wo = -path0[i].d;
                            funcVal = Spectrum(phase->eval(pRec));
                        }

                        if ( !funcVal.isZero() ) {
                            f = path0[i].thru1*path1[j].thru0*funcVal*mRec.transmittance/std::abs(path0[i].d.z);

                            if ( path1[j].surf ) {
                                bRec1.wi = ( j ? -path1[j - 1].d : wi );
                                bRec1.wo = -path0[i].d;
                                pdfJ[0] = path1[j].surf->pdf(bRec1);
                                bRec1.reverse();
                                pdfJ[1] = path1[j].surf->pdf(bRec1);
                            }
                            else {
                                pRec.wi = -path1[j - 1].d;
                                pRec.wo = -path0[i].d;
                                pdfJ[0] = phase->pdf(pRec);
                                pRec.reverse();
                                pdfJ[1] = phase->pdf(pRec);
                            }
                            pdfI[0] = vtxProb0[i][0];
#ifdef BIDIR_USE_ANALOG
                            pdfJ *= path1[j].pSurvival;
                            pdfI *= path0[i].pSurvival;
#endif
                            w = 1.0 + pdfJ[0]/pdfI[0];

                            if ( i ) {
                                w += ratio0[i]*( path0[i].surf ? mRec.pdfFailure : mRec.pdfSuccessRev )*pdfJ[0];
                            }
                            if ( j ) {
                                if ( flag && path1[j - 1].surf == NULL && path1[j].surf == NULL ) {
                                    sgm0.setMu(path1[j - 2].d);
                                    sgm1.setMu(-path0[i].d);
                                    sgm.init(sgm0, sgm1);
                                    w += mRec.pdfSuccess*sgm.pdf(path1[j - 1].d)/(vtxProb1[j - 1][0]*edgeProb1[j - 1][0]) +
                                         ratio1[j - 1]*mRec.pdfSuccess*pdfJ[1]*edgeProb1[j - 1][1]/edgeProb1[j - 1][0];
                                }
                                else {
                                    w1 = 1.0 + pdfJ[1]/vtxProb1[j - 1][0];
                                    w += ( path1[j].surf ? mRec.pdfFailure : mRec.pdfSuccess )/edgeProb1[j - 1][0]*w1;
                                    w += ratio1[j - 1]*( path1[j].surf ? mRec.pdfFailure : mRec.pdfSuccess )*pdfJ[1]*edgeProb1[j - 1][1]/edgeProb1[j - 1][0];
                                }
                            }
                            if ( !boost::math::isnormal(w) ) {
                                fprintf(stderr, "WTF (i): %lf\n", w);
                            }
                            ret += f/w;
                        }
                    }

                    // "j" path
                    t = (path1[j].x - path0[i].x)/-path1[j].d.z;
                    if ( t > Epsilon ) {
						med->eval(Ray(Point(0.0, 0.0, path1[j].x), path1[j].d, 0.0, t, 0.0), mRec);
						mRec.pdfSuccess /= std::abs(path1[j].d.z);
						mRec.pdfSuccessRev /= std::abs(path1[j].d.z);

                        if ( path0[i].surf ) {
                            bRec0.wi = ( i ? -path0[i - 1].d : wo );
                            bRec0.wo = -path1[j].d;
                            funcVal = path0[i].surf->eval(bRec0);
                        }
                        else {
                            pRec.wi = ( i ? -path0[i - 1].d : wo );
                            pRec.wo = -path1[j].d;
                            funcVal = Spectrum(phase->eval(pRec));
                        }

                        if ( !funcVal.isZero() ) {
                            f = path0[i].thru0*path1[j].thru1*funcVal*mRec.transmittance/std::abs(path1[j].d.z);

                            if ( path0[i].surf ) {
                                bRec0.wi = ( i ? -path0[i - 1].d : wo );
                                bRec0.wo = -path1[j].d;
                                pdfI[0] = path0[i].surf->pdf(bRec0);
                                bRec0.reverse();
                                pdfI[1] = path0[i].surf->pdf(bRec0);
                            }
                            else {
                                pRec.wi = -path0[i - 1].d;
                                pRec.wo = -path1[j].d;
                                pdfI[0] = phase->pdf(pRec);
                                pRec.reverse();
                                pdfI[1] = phase->pdf(pRec);
                            }
                            pdfJ[0] = vtxProb1[j][0];
#ifdef BIDIR_USE_ANALOG
                            pdfI *= path0[i].pSurvival;
                            pdfJ *= path1[j].pSurvival;
#endif
                            w = 1.0 + pdfI[0]/pdfJ[0];

                            if ( j ) {
                                w += ratio1[j]*( path1[j].surf ? mRec.pdfFailure : mRec.pdfSuccessRev )*pdfI[0];
                            }
                            if ( i ) {
                                if ( flag && path0[i - 1].surf == NULL && path0[i].surf == NULL ) {
                                    sgm0.setMu(path0[i - 2].d);
                                    sgm1.setMu(-path1[j].d);
                                    sgm.init(sgm0, sgm1);
                                    w += mRec.pdfSuccess*sgm.pdf(path0[i - 1].d)/(vtxProb0[i - 1][0]*edgeProb0[i - 1][0]) +
                                         ratio0[i - 1]*mRec.pdfSuccess*pdfI[1]*edgeProb0[i - 1][1]/edgeProb0[i - 1][0];
                                }
                                else {
                                    w1 = 1.0 + pdfI[1]/vtxProb0[i - 1][0];
                                    w += ( path0[i].surf ? mRec.pdfFailure : mRec.pdfSuccess )/edgeProb0[i - 1][0]*w1;
                                    w += ratio0[i - 1]*( path0[i].surf ? mRec.pdfFailure : mRec.pdfSuccess )*pdfI[1]*edgeProb0[i - 1][1]/edgeProb0[i - 1][0];
                                }
                            }
                            if ( !boost::math::isnormal(w) ) {
                                fprintf(stderr, "WTF (j): %lf\n", w);
                            }
                            ret += f/w;
                        }
                    }
                }
            }

        return ret;
    }


public:
    int run(int argc, char **argv) {
        const Float gVal = 0.5, params[] = { 8.4589e-02, 1.6381e+01, 1.1373e-01, 1.0000e+00, 2.7983e-01, 5.0540e+00 };
        //const Float gVal = 0.8, params[] = { 1.0759e-01, 2.0330e+00, 1.0008e+00, 1.4336e+01, 2.3987e+00, 6.5918e+01 };
        //const Float gVal = 0.95, params[] = { 2.7426e-01, 1.0261e+01, 3.8448e+00, 8.0110e+01, 5.7925e+01, 5.3259e+02 };

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

            Properties propBottom("roughconductor");
            propBottom.setFloat("alpha", 0.25);
            // Properties propBottom("diffuse");

            bsdfTop = static_cast<BSDF *>(PluginManager::getInstance()->createObject(MTS_CLASS(BSDF), propTop));
            bsdfTop->configure();
            bsdfBottom = static_cast<BSDF *>(PluginManager::getInstance()->createObject(MTS_CLASS(BSDF), propBottom));
            bsdfBottom->configure();
        }

        ref<Medium> med;
        {
            Properties propMedium("homogeneous");
            propMedium.setSpectrum("sigmaT", Spectrum(2.0));
            // propMedium.setSpectrum("sigmaT", Spectrum(0.0));
            propMedium.setSpectrum("albedo", Spectrum(0.9));

            Properties propPhase("hg");
            propPhase.setFloat("g", gVal);

            ref<PhaseFunction> phase = static_cast<PhaseFunction *>(PluginManager::getInstance()->createObject(MTS_CLASS(PhaseFunction), propPhase));
            phase->configure();

            med = static_cast<Medium *>(PluginManager::getInstance()->createObject(MTS_CLASS(Medium), propMedium));
            med->addChild(phase);
            med->configure();
        }
        std::cout << med->toString() << std::endl;

        Vector d0(3.0, 2.0, 1.0), d1(1.0, 2.0, 3.0);
        d0 = normalize(d0);
        d1 = normalize(d1);

        int flag = (argc >= 2 ? atoi(argv[1]) : 0);
        long long N = (argc >= 3 ? atoll(argv[2]) : 1000000);
        int K = (argc >= 4 ? atoi(argv[3]) : 10);

        std::cout << '\n' << "N = " << N << ", K = " << K << '\n';
        if ( flag ) std::cout << "Using SG sampling\n";
        std::cout << std::endl;

        std::vector<::Statistics> stats(nworker);
        std::chrono::time_point<std::chrono::steady_clock> start;

        start = std::chrono::steady_clock::now();
        for ( int i = 0; i < nworker; ++i ) stats[i].reset();
#ifdef MTS_OPENMP
#   pragma omp parallel for
#endif
        for ( long long omp_i = 0; omp_i < N/K; ++omp_i ) {
            int tid;
#ifdef MTS_OPENMP
            tid = omp_get_thread_num();
#else
            tid = 0;
#endif

            Spectrum val = simulate_bidir(d0, d1, bsdfTop, bsdfBottom, med, flag > 0, params, samplers[tid]);
            stats[tid].push(val[0]);
        }
        for ( int i = 1; i < nworker; ++i ) stats[0].push(stats[i]);

#ifdef BIDIR_USE_ANALOG
        Log(EInfo, "Bidir (ana.) : %.2le +- %.2le (%.1lf secs)", stats[0].getMean(), stats[0].getCI(),
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count()/1000.0);
#else
        Log(EInfo, "Bidir        : %.2le +- %.2le (%.1lf secs)", stats[0].getMean(), stats[0].getCI(),
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count()/1000.0);
#endif

        return 0;
    }

    MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(BidirSimFullEx, "Bidirectional Simulation for Layered Config (Extended Full)")
MTS_NAMESPACE_END
