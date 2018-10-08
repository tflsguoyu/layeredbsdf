#include <fstream>
#include <vector>
#include <string>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/ray.h>
#include <mitsuba/render/util.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/medium.h>

// #undef MTS_OPENMP

#if defined(MTS_OPENMP)
#   include <omp.h>
#endif

MTS_NAMESPACE_BEGIN
class LayeredUtil : public Utility {
public:

    void saveToFile(const std::vector< std::vector<Spectrum> > &bsdfMatrix, const char *filename, const Vector2i &outputSize) {

        ref<Bitmap> outputBitmap = new Bitmap(Bitmap::ERGB, Bitmap::EFloat32, outputSize);
        float *outputData = outputBitmap->getFloat32Data();
        for (int row = 0; row < outputSize.x; row++) {
            for (int col = 0; col < outputSize.y; col++) {
                Float r, g, b;
                bsdfMatrix[row][col].toLinearRGB(r, g, b);
                *outputData++ = (float) r;
                *outputData++ = (float) g;
                *outputData++ = (float) b;
            }
        }
        ref<FileStream> outputFile = new FileStream(filename, FileStream::ETruncReadWrite);
        outputBitmap->write(Bitmap::EOpenEXR, outputFile);
    }

    Spectrum Sample(RayDifferential &ray, RadianceQueryRecord &rRec, int &depth) {
        
        const Scene *scene = rRec.scene;
        Intersection &its = rRec.its;
        MediumSamplingRecord mRec;

        rRec.rayIntersect(ray);
        // scene->rayIntersect(ray, its);
        Spectrum throughput(1.0f);

        int maxDepth = -1;
        while (rRec.depth <= maxDepth || maxDepth < 0) {    
            if (rRec.medium && rRec.medium->sampleDistance(Ray(ray, 0, its.t), mRec, rRec.sampler) && mRec.p.z < 0 && mRec.p.z > -1) {
                const PhaseFunction *phase = rRec.medium->getPhaseFunction();
                throughput *= mRec.sigmaS * mRec.transmittance / mRec.pdfSuccess;

                PhaseFunctionSamplingRecord pRec(mRec, -ray.d);
                Float phaseVal = phase->sample(pRec, rRec.sampler);
                if (std::abs(phaseVal) < Epsilon)
                    break;
                throughput *= phaseVal;

                // Trace a ray
                ray = Ray(mRec.p, pRec.wo, 0.0f);
                ray.mint = 0;
                if ( !scene->rayIntersect(ray, its) ) {
                    throughput = Spectrum(0.0f);
                    break;
                }
            } else {
                if (rRec.medium)
                    throughput *= mRec.transmittance / mRec.pdfFailure;

                if (!its.isValid())
                    break;

                const BSDF *bsdf = its.getBSDF();
                BSDFSamplingRecord bRec(its, rRec.sampler);
                // bRec.mode = EImportance;
                Spectrum bsdfVal = bsdf->sample(bRec, rRec.nextSample2D());
                throughput *= bsdfVal;
                
                const Vector wo = its.toWorld(bRec.wo);
                if (its.isMediumTransition()) {
                    rRec.medium = its.getTargetMedium(wo);
                }
                ray = Ray(its.p, wo, 0.0f);
                if ( !scene->rayIntersect(ray, its) ) {
                    if ( rRec.medium ) {
                        throughput = Spectrum(0.0f);
                        break;
                    }
                }
            }
            rRec.depth++;
        }
        depth = rRec.depth;

        return throughput;
    }

    Spectrum sampleRefraction(RadianceQueryRecord &rRec, const Vector &wo, Vector &wo_refract, 
        bool flag_Jacobian) {
        
        Spectrum weight = Spectrum(1.0f);
 
        const Scene *scene = rRec.scene;
        Point originPoint = wo.z >= 0 ? Point(0.0f,0.0f,0.0f) : Point(0.0f,0.0f,-1.0f);
        RayDifferential ray(originPoint + wo, -wo, 0.0f);
        Intersection its;
        scene->rayIntersect(ray, its);

        const BSDF *bsdf = its.getBSDF();

        if (bsdf->getType() & BSDF::ENull) {
            wo_refract = wo;
            return weight;
        }
        
        BSDFSamplingRecord bRec(its, rRec.sampler);     
        bRec.typeMask = BSDF::ETransmission;
        bRec.mode = EImportance;
        weight = bsdf->sample(bRec, rRec.nextSample2D());   
        // weight *= Spectrum(std::pow(bRec.eta,4)); // If using EImportance, remove this line
        if ( !weight.isZero() ) {
            Float Jacobian1 = std::abs(Frame::cosTheta(bRec.wi)/Frame::cosTheta(bRec.wo))/(bRec.eta*bRec.eta);
            Float Jacobian2 = std::abs(Frame::cosTheta(bRec.wi)/Frame::cosTheta(bRec.wo));
            weight *= flag_Jacobian ? Jacobian2 : Jacobian1; // Jacobian!
            wo_refract = -its.toWorld(bRec.wo);
        }
        else
            wo_refract = Vector(0.0f);
        
        return weight;
    }

    Spectrum Eval(RayDifferential &ray, RadianceQueryRecord &rRec, const Vector wo_query) {
        if (wo_query.z == 0) 
            return Spectrum(0.0f);

        bool flag_type = wo_query.z > 0; // 1:reflection  0:transmission

        const Scene *scene = rRec.scene;
        Intersection &its = rRec.its;
        MediumSamplingRecord mRec;
        Spectrum Li(0.0f);

        bool flag_surface = true; // 1:top   0:bottom
        int topCounter = 0;
        int bottomCounter = 0;

        rRec.rayIntersect(ray);
        // scene->rayIntersect(ray, its);
        Spectrum throughput(1.0f);

        int maxDepth = -1;
        while (rRec.depth <= maxDepth || maxDepth < 0) {            
            
            if (rRec.medium && rRec.medium->sampleDistance(Ray(ray, 0, its.t), mRec, rRec.sampler) && mRec.p.z < 0 && mRec.p.z > -1) {

                const PhaseFunction *phase = rRec.medium->getPhaseFunction();
                throughput *= mRec.sigmaS * mRec.transmittance / mRec.pdfSuccess;

                // Direct
                DirectSamplingRecord dRec(mRec.p, 0.0f);
                int maxInteractions = maxDepth - rRec.depth;
            
                Vector wo_refract;
                bool flag_Jacobian = false;
                Spectrum weight_refract = sampleRefraction(rRec, wo_query, wo_refract, flag_Jacobian); 
                
                if ( !weight_refract.isZero() ) {
                    dRec.d = wo_refract;
                    Intersection its_nee;
                    Ray ray_nee = Ray(dRec.ref, dRec.d, 0.0f);
                    ray_nee.mint = 0.0f;
                    if (!scene->rayIntersect(ray_nee, its_nee)) {
                        cout << "ray.o:" << ray.o.toString() << " ray.d:" << ray.d.toString() << " ray.mint:" << ray.mint << " ray.maxt:" << its.t << endl;
                        cout << "ray_nee.o:" << ray_nee.o.toString() << " ray_nee.d:" << ray_nee.d.toString() << endl << endl;
                        cout << "[GY]: Warning in layeredUtil::Eval()::medium" << endl;
                        break;
                    }
                    dRec.p = its_nee.p; // Only work for single layer
                    
                    Spectrum value = scene->evalTransmittance_refract(dRec.ref, false, dRec.p, true, 
                        0.0f, rRec.medium, maxInteractions, rRec.sampler);

                    if (!value.isZero()) {
                        Li += throughput * value * phase->eval(
                                PhaseFunctionSamplingRecord(mRec, -ray.d, dRec.d)) * weight_refract;
                    }
                }

                // Indirect
                PhaseFunctionSamplingRecord pRec(mRec, -ray.d);
                Float phaseVal = phase->sample(pRec, rRec.sampler);
                if (std::abs(phaseVal) < Epsilon)
                    break;
                throughput *= phaseVal;

                ray = Ray(mRec.p, pRec.wo, 0.0f);
                ray.mint = 0;
                if ( !scene->rayIntersect(ray, its) )
                    break;

            } else {
    
                if (rRec.medium) {
                    throughput *= mRec.transmittance / mRec.pdfFailure;
                }

                if (!its.isValid())
                    break;
                
                const BSDF *bsdf = its.getBSDF();
                
                if (its.shape->getName() == "topsurface") {
                    ++topCounter;
                    flag_surface = true;
                } else if (its.shape->getName() == "bottomsurface") {
                    ++bottomCounter;
                    flag_surface = false;
                } else
                    Log(EWarn, "Unknown surface: %s", its.shape->getName().c_str());
                
                if (bsdf->getType() & BSDF::ESmooth) {
                    // Direct                   
                    
                    if (flag_type && flag_surface && topCounter ==1) { 
                        BSDFSamplingRecord bRec(its, its.toLocal(wo_query));
                        bRec.sampler = rRec.sampler;
                        Li += throughput * bsdf->eval(bRec);
                    } 

                    if ((flag_type && !flag_surface) || (!flag_type && flag_surface)) {                 
                        int maxInteractions = maxDepth - rRec.depth;

                        Vector wo_refract;
                        bool flag_Jacobian = !flag_type && flag_surface && topCounter == 1;
                        Spectrum weight_refract = sampleRefraction(rRec, wo_query, wo_refract, flag_Jacobian);

                        if ( !weight_refract.isZero() ) {
                            DirectSamplingRecord dRec(its);                     
                            dRec.d = wo_refract;
                            Intersection its_nee; 
                            if (!scene->rayIntersect(Ray(its.p, dRec.d, 0.0), its_nee)) {
                                cout << "[GY]: Warning in layeredUtil::Eval()::surface" << endl;
                                break;
                            }
                            dRec.p = its_nee.p; // Only work for single layer

                            if (its.shape && its.isMediumTransition()) {
                                rRec.medium = its.getTargetMedium(dRec.d);
                            }
                            
                            Spectrum value = Spectrum(1.0f);
                            value = scene->evalTransmittance_refract(its.p, true, dRec.p, true,
                                0.0f, rRec.medium, maxInteractions, rRec.sampler);

                            if (!value.isZero()) {
                                BSDFSamplingRecord bRec(its, its.toLocal(dRec.d));
                                bRec.sampler = rRec.sampler;
                                Li += throughput * bsdf->eval(bRec) * value * weight_refract;
                            }
                        }
                    }
                }

                // Indirect
                BSDFSamplingRecord bRec(its, rRec.sampler);
                bRec.mode = EImportance;
                Spectrum bsdfVal = bsdf->sample(bRec, rRec.nextSample2D());
                
                if (bsdfVal.isZero())
                    break;
                throughput *= bsdfVal;
                
                const Vector wo = its.toWorld(bRec.wo);
                if (its.isMediumTransition())
                    rRec.medium = its.getTargetMedium(wo);
                ray = Ray(its.p, wo, ray.time);
                if ( !scene->rayIntersect(ray, its) ) {
                    if ( rRec.medium ) break;
                }
            }
            rRec.depth++;
        }
        return Li;      
    }


    int run(int argc, char **argv) {
        int nworker;
#if defined(MTS_OPENMP)
        nworker = omp_get_max_threads();
#else
        nworker = 1;
#endif

        ref<Scene> scene = loadScene(argv[7]);
        scene->incRef();
        scene->initialize();

        ref_vector<Sampler> samplers(nworker);
        samplers[0] = static_cast<Sampler *>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), Properties("independent")));
        for ( int i = 1; i < nworker; ++i ) samplers[i] = samplers[0]->clone();

        const Point lightSource = Point(static_cast<Float>(std::atof(argv[4])), static_cast<Float>(std::atof(argv[5])), static_cast<Float>(std::atof(argv[6])));
        const Vector wi = normalize(lightSource - Point(0.0f,0.0f,0.0f));
        const RayDifferential ray_init = RayDifferential(lightSource, -wi, 0.0f);

        const Vector2i outputSize = Vector2i(std::atoi(argv[1]),std::atoi(argv[1]));
        const int numberOfBins = outputSize.x * outputSize.y;
        const long long N = std::atoll(argv[2]);

        // Sample
        {       
            std::cout << "Sample() ... " << std::endl;
            // our sample()
            std::vector< std::vector<Spectrum> > bsdfSampleMatrix(outputSize.y, std::vector<Spectrum>(outputSize.x, Spectrum(0.0f)));
            std::vector< std::vector<Spectrum> > bsdfSampleMatrix_B(outputSize.y, std::vector<Spectrum>(outputSize.x, Spectrum(0.0f)));
            
            std::vector<int> avgPathLength(10,0);

            long long tot = 0;
#if defined(MTS_OPENMP)
            #pragma omp parallel for
#endif
            for (int omp_i = 0; omp_i < N; omp_i++) {
                int tid;
#if defined(MTS_OPENMP)
                tid = omp_get_thread_num();
#else
                tid = 0;
#endif
                RayDifferential ray(ray_init);
                RadianceQueryRecord rRec(scene, samplers[tid]);
                rRec.type = RadianceQueryRecord::ERadiance;
                
                int depth = 0;
                Spectrum throughput = Sample(ray, rRec, depth);

#if defined(MTS_OPENMP)
                #pragma omp critical
                {
#endif
                    if (++tot % 5000000 == 0)
                        // std::cout << "Sample count: " << tot/1000000 << "M of " << (N+1)/1000000 << "M" << std::endl;
                        std::cout << "Sample count: " << tot/1000000 << "M of " << (N+1)/1000000 << "M\r" << std::flush;

                    if (depth >= 9) avgPathLength[9] ++;
                    else avgPathLength[depth] ++;

                    if (!throughput.isZero()) {
                        int col, row;
                        col = ray.d.x == 1.0f ? outputSize.x -1 : static_cast<int>(std::floor((ray.d.x + 1.0)/2.0 * outputSize.x));
                        row = ray.d.y == 1.0f ? outputSize.y -1 : static_cast<int>(std::floor((ray.d.y + 1.0)/2.0 * outputSize.y));
                                            
                        if (ray.d.z >= 0)
                            bsdfSampleMatrix[outputSize.x-row-1][col] += throughput*static_cast<Float>(numberOfBins)/static_cast<Float>(4*N)*ray.d.z;
                        else
                            bsdfSampleMatrix_B[outputSize.x-row-1][col] -= throughput*static_cast<Float>(numberOfBins)/static_cast<Float>(4*N)*ray.d.z;
                    }
#if defined(MTS_OPENMP)
                } //omp critical
#endif      
            }

            std::cout << std::endl;
            std::cout << "Sample() Done! Save to " << argv[8] << std::endl;

            ////// Save top Sample as .exr
            saveToFile(bsdfSampleMatrix, argv[8], outputSize);

            ///// Save bottom Sample as .exr
            if (argc >= 11) {
                saveToFile(bsdfSampleMatrix_B, argv[10], outputSize);
            }

            
            ////// print out how many raypaths in different path length
            for (int i = 0; i < 9; i++) {
                std::cout << "NumOfPath = " << i << ": " << avgPathLength[i] << std::endl;
            }
            std::cout << "NumOfPath >= 9: " << avgPathLength[9] << std::endl;
        }

        // Eval
        {
            std::cout << "Eval() ... " << std::endl;
            ///// Our eval()
            std::vector< std::vector<Spectrum> > bsdfEvalMatrix(outputSize.y, std::vector<Spectrum>(outputSize.x));
            std::vector< std::vector<Spectrum> > bsdfEvalMatrix_B(outputSize.y, std::vector<Spectrum>(outputSize.x));

            long long M = std::atoll(argv[3]);
            if (M == -1) {
                if (N >= numberOfBins)
                    M = static_cast<long long>(N / numberOfBins);
                else
                    M = 1;
            }
            std::cout << "Path per direction: " << M << std::endl;

#if defined(MTS_OPENMP)
            #pragma omp parallel for
#endif
            for (int omp_i = 0; omp_i < outputSize.x*outputSize.y; ++omp_i) {
                int tid;
#if defined(MTS_OPENMP)
                tid = omp_get_thread_num();
#else
                tid = 0;
#endif
                const int row = omp_i / outputSize.x, col = omp_i % outputSize.x;

                Spectrum bsdfEval = Spectrum(0.0f); 
                Spectrum bsdfEval_B = Spectrum(0.0f);

                Float x = (col+0.5f)*2.0f/outputSize.y-1;
                Float y = (row+0.5f)*2.0f/outputSize.x-1;
                Float z2 = x*x + y*y;
                if (z2 <= 1) {
                    // light goes out from top surface
                    Float z = std::sqrt(1-z2);  
                    const Vector wo = Point(x,y,z) - Point(0.0f,0.0f,0.0f);             
                    for (int i = 0; i < M; i++) {
                        RayDifferential ray(ray_init);
                        RadianceQueryRecord rRec(scene, samplers[tid]);
                        rRec.type = RadianceQueryRecord::ERadiance;
                        bsdfEval += Eval(ray, rRec, wo);
                    }
                    bsdfEval /= static_cast<Float>(M);

                    // light goes out from bottom surface
                    const Vector wo_B = Point(x,y,-z) - Point(0.0f,0.0f,0.0f);                  
                    for (int i = 0; i < M; i++) {
                        RayDifferential ray_B(ray_init);
                        RadianceQueryRecord rRec_B(scene, samplers[tid]);
                        rRec_B.type = RadianceQueryRecord::ERadiance;
                        bsdfEval_B += Eval(ray_B, rRec_B, wo_B);
                    }
                    bsdfEval_B /= static_cast<Float>(M);
                }

                bsdfEvalMatrix[outputSize.x-row-1][col] = bsdfEval;
                bsdfEvalMatrix_B[outputSize.x-row-1][col] = bsdfEval_B;
            }

            std::cout << "Eval() Done! Save to " << argv[9] << std::endl;

            ////// Save top Eval as .exr
            saveToFile(bsdfEvalMatrix, argv[9], outputSize);
            
            if (argc >= 12) {
                saveToFile(bsdfEvalMatrix_B, argv[11], outputSize);
            }

        }

        std::cout<<"Done!"<<std::endl;
        return 0;
    }

    MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(LayeredUtil, " ")
MTS_NAMESPACE_END