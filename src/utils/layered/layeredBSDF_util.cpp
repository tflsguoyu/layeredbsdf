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
  
#define MAX_DEPTH -1
#define MAX_TOP_SAMPLE -1
#define MAX_BOTTOM_SAMPLE -1
#define MAX_TOP_EVAL -1
#define MAX_BOTTOM_EVAL -1

MTS_NAMESPACE_BEGIN
class layeredBSDFUtil : public Utility {
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

    bool rayIntersect(const Ray &ray, Intersection &its) {

        its.p.x = std::numeric_limits<Float>::infinity();
        its.p.y = std::numeric_limits<Float>::infinity();
        its.p.z = std::numeric_limits<Float>::infinity();
        its.t = std::numeric_limits<Float>::infinity();
        its.wi = normalize(-ray.d);

        its.geoFrame = Frame(Normal(0.0f, 0.0f, 1.0f));
        its.shFrame = its.geoFrame;

        Float z = ray(ray.mint).z;
        
        if (std::abs(ray.d.z) < Epsilon || (z >= 0 && ray.d.z > 0) || (z <= -1 && ray.d.z < 0))
            return false;

        if ((z >= 0 && ray.d.z < 0) || (z < 0 && z > -1 && ray.d.z > 0))
            its.p.z = 0.0f;
        else if ((z <= -1 && ray.d.z > 0) || (z < 0 && z > -1 && ray.d.z < 0))
            its.p.z = -1.0f;
        else
            cout << "[GY]: Warning in layeredBSDFUtil::rayIntersect()" << endl;

        its.p.x = ray.o.x - ray.d.x * ((ray.o.z - its.p.z) / ray.d.z);
        its.p.y = ray.o.y - ray.d.y * ((ray.o.z - its.p.z) / ray.d.z);
        its.t = (its.p - ray.o).length();
        
        
        return true;

    }

    bool rayIntersectAndLookForEmitter(const Ray &ray, Intersection &its, const bool flag_type, bool &isEmitter) {

        isEmitter = false;

        its.p.x = std::numeric_limits<Float>::infinity();
        its.p.y = std::numeric_limits<Float>::infinity();
        its.p.z = std::numeric_limits<Float>::infinity();
        its.t = std::numeric_limits<Float>::infinity();
        its.wi = normalize(-ray.d);

        its.geoFrame = Frame(Normal(0.0f, 0.0f, 1.0f));
        its.shFrame = its.geoFrame;

        Float z = ray(ray.mint).z;

        if (std::abs(ray.d.z) < Epsilon || (z >= 0 && ray.d.z > 0) || (z <= -1 && ray.d.z < 0))
            return false;

        if ((z >= 0 && ray.d.z < 0) || (z < 0 && z > -1 && ray.d.z > 0))
            its.p.z = 0.0f;
        else if ((z <= -1 && ray.d.z > 0) || (z < 0 && z > -1 && ray.d.z < 0))
            its.p.z = -1.0f;
        else
            cout << "[GY]: Warning in layeredBSDFUtil::rayIntersect()" << endl;

        its.p.x = ray.o.x - ray.d.x * ((ray.o.z - its.p.z) / ray.d.z);
        its.p.y = ray.o.y - ray.d.y * ((ray.o.z - its.p.z) / ray.d.z);
        its.t = (its.p - ray.o).length();

        // 
        if ((std::abs(its.p.z - 0.0f) < Epsilon && flag_type) || (std::abs(its.p.z - (-1.0f)) < Epsilon && !flag_type)) {
            isEmitter = true;
        }

        return true;

    }

    Spectrum Sample(RayDifferential &ray, Sampler* sampler, int &depth, 
        const BSDF* bsdf_t, const BSDF* bsdf_b, const Medium* medium, const PhaseFunction* phase,
        const Normal& normal_t, const Normal& normal_b) {

        bool flag_medium = false;
        Intersection its;       
        MediumSamplingRecord mRec;

        rayIntersect(ray,its);
        Spectrum throughput(1.0f);

		int topCounter = 0;
		int bottomCounter = 0;
		while (depth <= MAX_DEPTH || MAX_DEPTH < 0) { 
            if (flag_medium && medium->sampleDistance(Ray(ray, 0, its.t), mRec, sampler) && mRec.p.z < 0 && mRec.p.z > -1) {
                throughput *= mRec.sigmaS * mRec.transmittance / mRec.pdfSuccess;

                PhaseFunctionSamplingRecord pRec(mRec, -ray.d);
                Float phaseVal = phase->sample(pRec, sampler);
                if (std::abs(phaseVal) < Epsilon) {
                    break;
                }
                throughput *= phaseVal;

                // Trace a ray
                ray = Ray(mRec.p, pRec.wo, 0.0f);
                ray.mint = 0.0f;
                if ( !rayIntersect(ray, its) ) {
                    throughput = Spectrum(0.0f);
                    break;
                }
            } else {
                if (flag_medium)
                    throughput *= mRec.transmittance / mRec.pdfFailure;
      
                if (!its.isValid())
                    break;

                const BSDF *bsdf = std::abs(its.p.z - 0.0f) < Epsilon ? bsdf_t : bsdf_b;
                const Normal normal = std::abs(its.p.z - 0.0f) < Epsilon ? normal_t : normal_b;

				if (std::abs(its.p.z - 0.0f) < Epsilon) {
					++topCounter;
					if (MAX_TOP_SAMPLE > -0.0001 && topCounter > MAX_TOP_SAMPLE) break;
				}
				else if (std::abs(its.p.z + 1.0f) < Epsilon) {
					++bottomCounter;
					if (MAX_BOTTOM_SAMPLE > -0.0001 && bottomCounter > MAX_BOTTOM_SAMPLE) break;
				}
				else
					cout << "[GY]: Warning in layeredBSDFUtil::Sample" << endl;

                // perturbe shading normal for its
                its.shFrame = Frame(normalize(normal));
                if (its.wi.z * its.toLocal(its.wi).z <= 0) {
                    throughput = Spectrum(0.0f);
                    break;
                }
                
                BSDFSamplingRecord bRec(its, sampler);
                bRec.wi = its.toLocal(its.wi);
                
                bRec.mode = EImportance;
                Spectrum bsdfVal = bsdf->sample(bRec, sampler->next2D());
                if (bRec.wo.z * its.toWorld(bRec.wo).z <= 0) {
                    throughput = Spectrum(0.0f);
                    break;
                }
                
                throughput *= bsdfVal;
                
                const Vector wo = its.toWorld(bRec.wo);
                if (its.wi.z * wo.z < 0) {
                    flag_medium = !flag_medium;
                }
                ray = Ray(its.p, wo, 0.0f);
                if ( !rayIntersect(ray, its) ) {
                    if ( flag_medium ) {
                        throughput = Spectrum(0.0f);
                        break;
                    }
                }
            }
            depth++;
        }
       
        return throughput;
    }

    inline Float miWeight(Float pdfA, Float pdfB) const {
        pdfA *= pdfA; pdfB *= pdfB;
		//if ((pdfA + pdfB) < Epsilon) {
		//	cout << "[GY]: Warning in LayeredNormalmapMIS_util::miWeight()" << endl;
		//	cout << "pdfA: " << pdfA << endl;
		//	cout << "pdfB: " << pdfB << endl;
		//	cout << "A/(A+B): " << pdfA / (pdfA + pdfB) << endl;
		//}
        return pdfA / (pdfA + pdfB);
    }

    Spectrum sampleRefraction(Sampler* sampler, const BSDF* bsdf, const Normal& normal, 
        const Vector &wo_query, Vector &wo_refract, Float &pdf_refract) {
        
        Spectrum weight = Spectrum(1.0f);
        wo_refract = Vector(0.0f);

        if (bsdf->getType() & BSDF::ENull) {
            wo_refract = wo_query;
            return weight;
        }

        Intersection its;
        its.wi = wo_query;
        its.shFrame = Frame(normalize(normal));
        if (its.wi.z * its.toLocal(its.wi).z <= 0) {
            return Spectrum(0.0f);
        }
        BSDFSamplingRecord bRec(its, sampler);
        bRec.wi = its.toLocal(its.wi);

        //bRec.typeMask = BSDF::ETransmission;
        bRec.mode = EImportance;
		weight = bsdf->sample(bRec, pdf_refract, sampler->next2D());
		if (bRec.wi.z * bRec.wo.z >= 0) {
			return Spectrum(0.0f);
		}
		if (bRec.wo.z * its.toWorld(bRec.wo).z <= 0) {
            return Spectrum(0.0f);
        }
        const Vector wo = its.toWorld(bRec.wo);
        
        if ( !weight.isZero() ) {
			Float Jacobian = std::abs(bRec.wo.z / bRec.wi.z) * (bRec.eta * bRec.eta);
	        weight /= Jacobian; // Jacobian!
            //pdf_refract *= flag_Jacobian ? Jacobian2 : Jacobian1;
            wo_refract = -wo;
        }
        else {
            return Spectrum(0.0f);
        }

        return weight;
    }

    Spectrum evaluateRefraction(Sampler* sampler, const BSDF* bsdf, const Normal& normal,
        const Vector &wo_query, const Intersection &its, Float &pdf_refract) {

        Spectrum weight = Spectrum(1.0f);

        if (bsdf->getType() & BSDF::ENull) {
            pdf_refract = 0.0f;
            return Spectrum(0.0f);
        }

        if (wo_query.z * its.toLocal(wo_query).z <= 0) {
            pdf_refract = 0.0f;
            return Spectrum(0.0f);
        }
		if (its.wi.z * its.toLocal(its.wi).z <= 0) {
			pdf_refract = 0.0f;
			return Spectrum(0.0f);
		}

        BSDFSamplingRecord bRec(its, its.toLocal(wo_query), its.toLocal(its.wi));
		bRec.sampler = sampler;
		//bRec.typeMask = BSDF::ETransmission;
        bRec.mode = EImportance;
      
		weight = bsdf->eval(bRec);
		pdf_refract = bsdf->pdf(bRec);
        Float eta = bsdf->getEta();
		eta = bRec.wo.z < 0 ? eta : 1.0f / eta;
		
        if (!weight.isZero()) {
            Float Jacobian = std::abs(bRec.wo.z / bRec.wi.z) * (eta*eta);
			weight /= Jacobian; // Jacobian!
            //pdf_refract *= Jacobian; // Jacobian!
        }
        else {
            return Spectrum(0.0f);
        }

        return weight;
    }


    Spectrum Eval(RayDifferential &ray, Sampler *sampler, const Vector wo_query, int &depth,
        const BSDF* bsdf_t, const BSDF* bsdf_b, const Medium* medium, const PhaseFunction* phase,
        const Normal& normal_t, const Normal& normal_b, const bool MISenable = true) {
        
		const Vector wi_init = -ray.d;
		if (wo_query.z == 0.0f) 
            return Spectrum(0.0f);

        bool flag_type = (wi_init.z * wo_query.z) > 0; // 1:reflection  0:transmission
		bool flag_incidentDir = wi_init.z > 0;

        bool flag_medium = false;
        Intersection its;
        MediumSamplingRecord mRec;
        Spectrum Li(0.0f);

        bool flag_surface = true; // 1:top   0:bottom
        int topCounter = 0;
        int bottomCounter = 0;

        rayIntersect(ray,its);
        Spectrum throughput(1.0f);

        while (depth <= MAX_DEPTH || MAX_DEPTH < 0) {                     
            if (flag_medium && medium->sampleDistance(Ray(ray, 0.0f, its.t), mRec, sampler) && mRec.p.z < 0 && mRec.p.z > -1) {
                if (mRec.p.z >= 0 || mRec.p.z <= -1) {
                    cout << "[GY]: Warning in layeredBSDFUtil::Eval()::medium->sampleDistance" << endl; 
                    cout << "ray.o.z:" << ray.o.z << " ray.d.z:" << ray.d.z << " its.p.z:" << its.p.z << " mRec.p.z:" << mRec.p.z << endl << endl;
                }
                throughput *= mRec.sigmaS * mRec.transmittance / mRec.pdfSuccess;

                // Direct
                Vector wo_refract;
                Float refractPdf;
                const BSDF* bsdf_nee = (flag_incidentDir && flag_type) || (!flag_incidentDir && !flag_type) ? bsdf_t : bsdf_b;
                const Normal normal_nee = (flag_incidentDir && flag_type) || (!flag_incidentDir && !flag_type) ? normal_t : normal_b;
				Spectrum throughput_refract = sampleRefraction(sampler, bsdf_nee, normal_nee,
                    wo_query, wo_refract, refractPdf); 
                
                if ( !throughput_refract.isZero() ) {
                    Intersection its_nee;
                    Ray ray_nee = Ray(mRec.p, wo_refract, 0.0f);
                    ray_nee.mint = 0.0f;
                    if (!rayIntersect(ray_nee, its_nee)) {
                        cout << "[GY]: Warning in layeredBSDFUtil::Eval()::medium" << endl;
                        cout << "ray.o:" << ray.o.toString() << " ray.d:" << ray.d.toString() << " ray.mint:" << ray.mint << " ray.maxt:" << its.t << endl;
                        cout << "ray_nee.o:" << ray_nee.o.toString() << " ray_nee.d:" << ray_nee.d.toString() << endl << endl;
                        break;
                    }
                    Spectrum value = medium->evalTransmittance(Ray(ray_nee, 0.0f, its_nee.t), sampler);

                    if (!value.isZero()) {
                        PhaseFunctionSamplingRecord pRec(mRec, -ray.d, wo_refract);
                        Float phaseVal = phase->eval(pRec);
                        Float phasePdf = phase->pdf(pRec);

                        /* Weight using the power heuristic */
                        const Float weight = miWeight(refractPdf, phasePdf);  
						if (MISenable)
							Li += throughput * phaseVal * value * throughput_refract * weight;
						else
							Li += throughput * phaseVal * value * throughput_refract;
                    }
                }

                // Indirect
                Float phasePdf;
                PhaseFunctionSamplingRecord pRec(mRec, -ray.d);
                Float phaseVal = phase->sample(pRec, phasePdf, sampler);
                
                if (std::abs(phaseVal) < Epsilon)
                    break;
                throughput *= phaseVal;

                ray = Ray(mRec.p, pRec.wo, 0.0f);
                ray.mint = 0.0f;
                bool isEmitter = false;
                if (!rayIntersectAndLookForEmitter(ray, its, (flag_incidentDir && flag_type) || (!flag_incidentDir && !flag_type), isEmitter)) {
                    break;
                }
                
				if (MISenable) {
					if (isEmitter) {
						const BSDF *bsdf = std::abs(its.p.z - 0.0f) < Epsilon ? bsdf_t : bsdf_b;
						const Normal normal = std::abs(its.p.z - 0.0f) < Epsilon ? normal_t : normal_b;
						its.shFrame = Frame(normalize(normal));

						Spectrum refractVal = evaluateRefraction(sampler, bsdf, normal, 
							wo_query, its, refractPdf);

						if (!refractVal.isZero()) {
							Spectrum value = medium->evalTransmittance(Ray(ray, 0.0f, its.t), sampler);
							const Float weight = miWeight(phasePdf, refractPdf);

							Li += throughput * value * refractVal * weight;
						}
					}
				}

            } else {
    
                if (flag_medium) {
                    throughput *= mRec.transmittance / mRec.pdfFailure;
                }

                if (!its.isValid())
                    break;
                
                const BSDF *bsdf = std::abs(its.p.z - 0.0f) < Epsilon ? bsdf_t : bsdf_b;
                const Normal normal = std::abs(its.p.z - 0.0f) < Epsilon ? normal_t : normal_b;

                if (std::abs(its.p.z - 0.0f) < Epsilon) {
					++topCounter;
					if (MAX_TOP_EVAL > -0.0001 && topCounter > MAX_TOP_EVAL) break;
					flag_surface = true;
                } else if (std::abs(its.p.z + 1.0f) < Epsilon) {
					++bottomCounter;
					if (MAX_BOTTOM_EVAL > -0.0001 && bottomCounter > MAX_BOTTOM_EVAL) break;
					flag_surface = false;
                } else
                    cout << "[GY]: Warning in layeredBSDFUtil::Eval" << endl;
                
                if (bsdf->getType() & BSDF::ESmooth) {
                    // Direct                                    
                    if ((flag_incidentDir && flag_type && flag_surface && topCounter == 1) ||
						(!flag_incidentDir && flag_type && !flag_surface && bottomCounter == 1)) {
                        its.shFrame = Frame(normalize(normal));
                        if (its.wi.z * its.toLocal(its.wi).z <= 0) {
                            break;
                        }
						if (wo_query.z * its.toLocal(wo_query).z <= 0) {
							break;
						}
                        BSDFSamplingRecord bRec(its, its.toLocal(its.wi), its.toLocal(wo_query));
                        bRec.sampler = sampler;
                        Li += throughput * bsdf->eval(bRec);
                    } 

                    if ((flag_incidentDir && flag_type && !flag_surface) ||
						(flag_incidentDir && !flag_type && flag_surface) ||
						(!flag_incidentDir && flag_type && flag_surface) ||
						(!flag_incidentDir && !flag_type && !flag_surface)) {

                        Vector wo_refract;
                        Float refractPdf;
						const BSDF* bsdf_nee = (flag_incidentDir && flag_type) || (!flag_incidentDir && !flag_type) ? bsdf_t : bsdf_b;
                        const Normal normal_nee = (flag_incidentDir && flag_type) || (!flag_incidentDir && !flag_type) ? normal_t : normal_b;
						Float eta_all = bsdf_t->getEta() * bsdf_b->getEta();
						Spectrum throughput_refract = sampleRefraction(sampler, bsdf_nee, normal_nee, 
							wo_query, wo_refract, refractPdf);

                        if ( !throughput_refract.isZero() ) {
                            Intersection its_nee; 
                            Ray ray_nee = Ray(its.p, wo_refract, 0.0);
                            if (!rayIntersect(ray_nee, its_nee)) {
								cout << "throughput_refract: " << throughput_refract.toString() << endl;
								cout << "wo_refract: " << wo_refract.z << endl;
								
                                cout << "[GY]: Warning in layeredBSDFUtil::Eval()::surface" << endl;
                                break;
                            }
                            Spectrum value = medium->evalTransmittance(Ray(ray_nee, 0, its_nee.t), sampler);

                            if (!value.isZero()) {
                                its.shFrame = Frame(normalize(normal));
                                if (its.wi.z * its.toLocal(its.wi).z <= 0) {
                                    break;
                                }
								if (wo_refract.z * its.toLocal(wo_refract).z <= 0) {
									break;
								}
                                BSDFSamplingRecord bRec(its, its.toLocal(its.wi), its.toLocal(wo_refract));
                                bRec.sampler = sampler;
								bRec.mode = EImportance;
								Spectrum bsdfVal = bsdf->eval(bRec);
                                const Float bsdfPdf = bsdf->pdf(bRec);
                                const Float weight = miWeight(refractPdf, bsdfPdf);
								if (MISenable)
									Li += throughput * bsdfVal * value * throughput_refract * weight;
								else
									Li += throughput * bsdfVal * value * throughput_refract;
                            }
                        }
                    }
                }

                // Indirect
                its.shFrame = Frame(normalize(normal));
                if (its.wi.z * its.toLocal(its.wi).z <= 0) {
                    break;
                }
                BSDFSamplingRecord bRec(its, sampler);
                bRec.wi = its.toLocal(its.wi);
                
                bRec.mode = EImportance;
                Float bsdfPdf;
                Spectrum bsdfWeight = bsdf->sample(bRec, bsdfPdf, sampler->next2D());
                if (bsdfWeight.isZero())
                    break;
                
                const Vector wo = its.toWorld(bRec.wo);
                if (bRec.wo.z * its.toWorld(bRec.wo).z <= 0) {
                    break;
                }

                ray = Ray(its.p, wo, 0.0f);

                throughput *= bsdfWeight;
                if (its.wi.z * wo.z < 0) {
                    flag_medium = !flag_medium;
                }
                
                bool isEmitter = false;
                if (!rayIntersectAndLookForEmitter(ray, its, (flag_incidentDir && flag_type) || (!flag_incidentDir && !flag_type), isEmitter)) {
                    if ( flag_medium ) 
                        break;
                }

				if (MISenable) {
					if (isEmitter) {
						const BSDF *bsdf = std::abs(its.p.z - 0.0f) < Epsilon ? bsdf_t : bsdf_b;
						const Normal normal = std::abs(its.p.z - 0.0f) < Epsilon ? normal_t : normal_b;
						its.shFrame = Frame(normalize(normal));

						Float refractPdf;
						Spectrum refractVal = evaluateRefraction(sampler, bsdf, normal, 
							wo_query, its, refractPdf);

						if (!refractVal.isZero()) {
							Spectrum value = medium->evalTransmittance(Ray(ray, 0.0f, its.t), sampler);
							const Float weight = miWeight(bsdfPdf, refractPdf);

							Li += throughput * value * refractVal * weight;
						}
					}
				}

            }
            depth++;
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

        ref_vector<Sampler> samplers(nworker);
        samplers[0] = static_cast<Sampler *>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), Properties("independent")));
        for ( int i = 1; i < nworker; ++i ) samplers[i] = samplers[0]->clone();
        
        const Normal normal_top = Normal(static_cast<Float>(std::atof(argv[7])), static_cast<Float>(std::atof(argv[8])), static_cast<Float>(std::atof(argv[9])));
        const Normal normal_bottom = Normal(static_cast<Float>(std::atof(argv[10])), static_cast<Float>(std::atof(argv[11])), static_cast<Float>(std::atof(argv[12])));
        
        // Incident ray
        const Point lightSource = Point(static_cast<Float>(std::atof(argv[4])), static_cast<Float>(std::atof(argv[5])), static_cast<Float>(std::atof(argv[6])));
		bool flag_incidentDir = lightSource.z > 0;
		Point originPoint = flag_incidentDir ? Point(0.0f, 0.0f, 0.0f) : Point(0.0f, 0.0f, -1.0f);
		const Vector wi = normalize(lightSource - originPoint);
        const RayDifferential ray_init = RayDifferential(lightSource, -wi, 0.0f);

		ref<Scene> scene = loadScene(argv[1]);
		scene->incRef();
		scene->initialize();

		const RayDifferential ray1 = RayDifferential(Point(0.0f, 0.0f, 1.0f), Vector(0.0f, 0.0f, -1.0f), 0.0f);
		Intersection its1;
		scene->rayIntersect(ray1, its1);
		const BSDF* topSurfaceBSDF = its1.getBSDF();
		const Medium* interiorMedium = its1.getTargetMedium(Vector(0.0f, 0.0f, -1.0f));
		const PhaseFunction* interiorPhaseFunction = interiorMedium->getPhaseFunction();

		const RayDifferential ray2 = RayDifferential(Point(0.0f, 0.0f, -9.0f), Vector(0.0f, 0.0f, -1.0f), 0.0f);
		Intersection its2;
		scene->rayIntersect(ray2, its2);
		const BSDF* bottomSurfaceBSDF = its2.getBSDF();


		//// Top surface BSDF
		//Properties topSurfaceProps = Properties("roughdielectric");
		//topSurfaceProps.setFloat("intIOR", Float(2.5f));
		//topSurfaceProps.setFloat("extIOR", Float(1.0f));
		//topSurfaceProps.setFloat("alpha", Float(0.2f));
		//ref<BSDF> topSurfaceBSDF = static_cast<BSDF *> (PluginManager::getInstance()->
		//	createObject(MTS_CLASS(BSDF), topSurfaceProps));
		//topSurfaceBSDF->configure();

		//// Bottom surface BSDF
		//Properties bottomSurfaceProps = Properties("roughdielectric");
		//bottomSurfaceProps.setFloat("intIOR", Float(1.0f));
		//bottomSurfaceProps.setFloat("extIOR", Float(2.5f));
		//bottomSurfaceProps.setFloat("alpha", Float(0.2f));
		//ref<BSDF> bottomSurfaceBSDF = static_cast<BSDF *> (PluginManager::getInstance()->
		//	createObject(MTS_CLASS(BSDF), bottomSurfaceProps));
		//bottomSurfaceBSDF->configure();

		//// Interior Medium
		//Properties interiorProps = Properties("homogeneous");
		//interiorProps.setSpectrum("sigmaT", Spectrum(0.5f));
		//Spectrum albedo; albedo.fromLinearRGB(0.1f, 0.5f, 0.9f);
		//interiorProps.setSpectrum("albedo", albedo);
		////interiorProps.setSpectrum("sigmaS", Spectrum(0.0f));
		////Spectrum sigmaA; sigmaA.fromLinearRGB(0.9f, 0.5f, 0.1f);
		////interiorProps.setSpectrum("sigmaA", sigmaA);
		//ref<Medium> interiorMedium = static_cast<Medium *> (PluginManager::getInstance()->
		//	createObject(MTS_CLASS(Medium), interiorProps));
		//interiorMedium->configure();

		//// Phase
		//Properties interiorPhaseProps = Properties("isotropic");
		//ref<PhaseFunction> interiorPhaseFunction = static_cast<PhaseFunction *> (PluginManager::getInstance()->
		//	createObject(MTS_CLASS(PhaseFunction), interiorPhaseProps));
		//interiorPhaseFunction->configure();

        // Pre-defined
        const Vector2i outputSize = Vector2i(std::atoi(argv[3]),std::atoi(argv[3]));
		//const Vector2i outputSize = Vector2i(512, 512);
        const int numberOfBins = outputSize.x * outputSize.y;
        const long long N = std::atoll(argv[2]);

        // Sample
        {
            std::cout << "Sample() ... " << std::endl;
            // our sample()
            std::vector< std::vector<Spectrum> > bsdfSampleMatrix(outputSize.y, std::vector<Spectrum>(outputSize.x, Spectrum(0.0f)));
            std::vector< std::vector<Spectrum> > bsdfSampleMatrix_B(outputSize.y, std::vector<Spectrum>(outputSize.x, Spectrum(0.0f)));

            std::vector<int> avgPathLength(10, 0);

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

                int depth = 0;
                Spectrum throughput = Sample(ray, samplers[tid], depth,
                    topSurfaceBSDF, bottomSurfaceBSDF, interiorMedium, interiorPhaseFunction,
                    normalize(normal_top), normalize(normal_bottom));

#if defined(MTS_OPENMP)
#pragma omp critical
                {
#endif
                    if (++tot % 5000000 == 0)
                        // std::cout << "Sample count: " << tot/1000000 << "M of " << (N+1)/1000000 << "M" << std::endl;
                        std::cout << "Sample count: " << tot / 1000000 << "M of " << (N + 1) / 1000000 << "M\r" << std::flush;

                    if (depth >= 9) avgPathLength[9] ++;
                    else avgPathLength[depth] ++;

                    if (!throughput.isZero()) {
                        int col, row;
                        col = ray.d.x == 1.0f ? outputSize.x - 1 : static_cast<int>(std::floor((ray.d.x + 1.0) / 2.0 * outputSize.x));
                        row = ray.d.y == 1.0f ? outputSize.y - 1 : static_cast<int>(std::floor((ray.d.y + 1.0) / 2.0 * outputSize.y));

                        if (ray.d.z >= 0)
                            bsdfSampleMatrix[outputSize.x - row - 1][col] += throughput*static_cast<Float>(numberOfBins) / static_cast<Float>(4 * N)*ray.d.z;
                        else
                            bsdfSampleMatrix_B[outputSize.x - row - 1][col] -= throughput*static_cast<Float>(numberOfBins) / static_cast<Float>(4 * N)*ray.d.z;
                    }
#if defined(MTS_OPENMP)
                } //omp critical
#endif      
            }

            std::cout << std::endl;
            std::cout << "Sample() Done! Save to " << argv[13] << std::endl;

            ////// Save top Sample as .exr
            saveToFile(bsdfSampleMatrix, argv[13], outputSize);

            ///// Save bottom Sample as .exr
            if (argc >= 16) {
                saveToFile(bsdfSampleMatrix_B, argv[15], outputSize);
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

            // long long M = std::atoll(argv[3]);
			long long M = -1;
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
                        int depth = 0;
                        bsdfEval += Eval(ray, samplers[tid], wo, depth, 
                            topSurfaceBSDF, bottomSurfaceBSDF, interiorMedium, interiorPhaseFunction, 
                            normalize(normal_top), normalize(normal_bottom));
                    }
                    bsdfEval /= static_cast<Float>(M);

                    // light goes out from bottom surface
                    const Vector wo_B = Point(x,y,-z) - Point(0.0f,0.0f,0.0f);                  
                    for (int i = 0; i < M; i++) {
                        RayDifferential ray_B(ray_init);
                        int depth = 0;
                        bsdfEval_B += Eval(ray_B, samplers[tid], wo_B, depth,
                            topSurfaceBSDF, bottomSurfaceBSDF, interiorMedium, interiorPhaseFunction,
                            normalize(normal_top), normalize(normal_bottom));
                    }
                    bsdfEval_B /= static_cast<Float>(M);
                }

                bsdfEvalMatrix[outputSize.x-row-1][col] = bsdfEval;
                bsdfEvalMatrix_B[outputSize.x-row-1][col] = bsdfEval_B;
            }

            std::cout << "Eval() Done! Save to " << argv[14] << std::endl;

            ////// Save top Eval as .exr
            saveToFile(bsdfEvalMatrix, argv[14], outputSize);
            
            if (argc >= 17) {
                saveToFile(bsdfEvalMatrix_B, argv[16], outputSize);
            }

        }

        std::cout<<"Done!"<<std::endl;
        return 0;
    }

    MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(layeredBSDFUtil, " ")
MTS_NAMESPACE_END