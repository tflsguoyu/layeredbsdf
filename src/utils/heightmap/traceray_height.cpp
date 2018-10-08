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

//#undef MTS_OPENMP

#if defined(MTS_OPENMP)
#   include <omp.h>
#endif

#define MAX_DEPTH -1
//#define SCENE

MTS_NAMESPACE_BEGIN

class traceRay : public Utility {
public:

	void saveToFile(const std::vector< std::vector<Spectrum> > &bsdfMatrix, const char *filename, const Vector2i &outputSize) {

		ref<Bitmap> outputBitmap = new Bitmap(Bitmap::ERGB, Bitmap::EFloat32, outputSize);
		float *outputData = outputBitmap->getFloat32Data();
		for (int row = 0; row < outputSize.x; row++) {
			for (int col = 0; col < outputSize.y; col++) {
				Float r, g, b;
				bsdfMatrix[row][col].toLinearRGB(r, g, b);
				*outputData++ = (float)r;
				*outputData++ = (float)g;
				*outputData++ = (float)b;
			}
		}
		ref<FileStream> outputFile = new FileStream(filename, FileStream::ETruncReadWrite);
		outputBitmap->write(Bitmap::EOpenEXR, outputFile);
	}

	/// Temporary storage for patch-ray intersections
	struct PatchIntersectionRecord {
		Point p;
		int x, y;
	};

	bool rayIntersect(const Ray &ray, Intersection &its, const Shape* shape) {

		//Float nearT, farT;
		//if (!shape->getAABB().rayIntersect(ray, nearT, farT)) {
		//	return false;
		//}

		//Float t;

		PatchIntersectionRecord tmp;
		if (!shape->rayIntersect(ray, ray.mint, ray.maxt, its.t, &tmp)) {
			return false;
		}
		
		shape->fillIntersectionRecord(ray, &tmp, its);

		computeShadingFrame(its.shFrame.n, its.dpdu, its.shFrame);
		its.wi = its.toLocal(-ray.d);

		return true;
	}


#ifdef SCENE
	Spectrum Sample(const Scene* scene, Ray &ray, Sampler* sampler, int &depth,
		const BSDF* bsdf) {
#else
	 Spectrum Sample(const Shape* heightfield, Ray &ray, Sampler* sampler, int &depth,
		const BSDF* bsdf) {
#endif
	    Intersection its;
    
#ifdef SCENE
		if (!scene->rayIntersect(ray, its)) {
#else
		if (!rayIntersect(ray, its, heightfield)) {
#endif
			cout << "Not intersect.." << endl;
		}
		Spectrum throughput(1.0f);

        while (depth <= MAX_DEPTH || MAX_DEPTH < 0) {
            if (!its.isValid())
                break;

            BSDFSamplingRecord bRec(its, sampler);
            Spectrum bsdfVal = bsdf->sample(bRec, sampler->next2D());
            throughput *= bsdfVal;

            const Vector wo = its.toWorld(bRec.wo);
			ray = Ray(its.p, wo, 0.0);
			
#ifdef SCENE
			if (!scene->rayIntersect(ray, its)) {
#else
			if (!rayIntersect(ray, its, heightfield)) {
#endif
				break;
            }
		    depth++;
        }
        
	    return throughput;
    }

#ifdef SCENE
	Spectrum Eval(const Scene* scene, Ray &ray, Sampler *sampler, const Vector wo_query, int &depth,
		const BSDF* bsdf) {
#else
		Spectrum Eval(const Shape* heightfield, Ray &ray, Sampler *sampler, const Vector wo_query, int &depth,
		const BSDF* bsdf) {
#endif

		Intersection its;
        Spectrum Li(0.0f);
		Spectrum throughput(1.0f);

#ifdef SCENE
        scene->rayIntersect(ray, its);
#else
		rayIntersect(ray, its, heightfield);
#endif

        while (depth <= MAX_DEPTH || MAX_DEPTH < 0) {
        
                if (!its.isValid())
                    break;

                // Direct
				Ray ray_nee = Ray(its.p, wo_query, 0.0);
				Intersection its_nee;
#ifdef SCENE
				if (!scene->rayIntersect(ray_nee, its_nee)) {
#else
				if (!rayIntersect(ray_nee, its_nee, heightfield)) {
#endif
					BSDFSamplingRecord bRec_Direct(its, its.toLocal(wo_query));
					bRec_Direct.sampler = sampler;
					Li += throughput * bsdf->eval(bRec_Direct);
				}
                // Indirect
                BSDFSamplingRecord bRec_inDirect(its, sampler);
				Spectrum bsdfVal = bsdf->sample(bRec_inDirect, sampler->next2D());

                if (bsdfVal.isZero())
                    break;
                throughput *= bsdfVal;

                const Vector wo = its.toWorld(bRec_inDirect.wo);
                ray = Ray(its.p, wo, 0.0);

#ifdef SCENE
				if (!scene->rayIntersect(ray, its)) {
#else
				if (!rayIntersect(ray, its, heightfield)) {
#endif		
					break;
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

        const Point lightSource = Point(-0.70710678f, 0.0f, 0.70710678f)*10;
        const Vector wi = normalize(lightSource - Point(0.0f,0.0f,0.0f));
        const Ray ray_init = Ray(lightSource, -wi, 0.0f);

#ifdef SCENE
		// scene
		ref<Scene> scene = loadScene(argv[1]);
		scene->incRef();
		scene->initialize();
#else
		// shape
		Properties shapeProps = Properties("heightfield");
		shapeProps.setString("filename", argv[1]);
		ref<Shape> heightfieldShape = static_cast<Shape *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(Shape), shapeProps));
		heightfieldShape->configure();
#endif

		// Top surface BSDF
		Properties surfaceProps = Properties("diffuse");
		ref<BSDF> surfaceBSDF = static_cast<BSDF *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(BSDF), surfaceProps));
		surfaceBSDF->configure();

        const Vector2i outputSize = Vector2i(400,400);
        const int numberOfBins = outputSize.x * outputSize.y;
        const long long N = std::atoll(argv[2]);

        // Sample
		{
			// our sample() disk
			std::vector< std::vector<Spectrum> > bsdfSampleMatrix(outputSize.y, std::vector<Spectrum>(outputSize.x, Spectrum(0.0f)));
			std::vector< std::vector<Spectrum> > bsdfSampleMatrix_B(outputSize.y, std::vector<Spectrum>(outputSize.x, Spectrum(0.0f)));

			std::vector<int> avgPathLength(10, 0);
			long long tot = 0, totMissed = 0;
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
				Ray ray(ray_init);

				int depth = 0;
#ifdef SCENE
				Spectrum throughput = Sample(scene, ray, samplers[tid], depth,
					surfaceBSDF);
#else			
				Spectrum throughput = Sample(heightfieldShape, ray, samplers[tid], depth,
					surfaceBSDF);
#endif

#if defined(MTS_OPENMP)
#pragma omp critical
				{
#endif
					if (++tot % 5000000 == 0)
						std::cout << "Sample count: " << tot / 1000000 << "M of " << (N + 1) / 1000000 << "M\r" << std::flush;

					if (depth >= 9) avgPathLength[9] ++;
					else avgPathLength[depth] ++;

					if (throughput != Spectrum(0.0f)) {
						int col, row;
						col = ray.d.x == 1.0f ? outputSize.x - 1 : static_cast<int>(std::floor((ray.d.x + 1.0) / 2.0 * outputSize.x));
						row = ray.d.y == 1.0f ? outputSize.y - 1 : static_cast<int>(std::floor((ray.d.y + 1.0) / 2.0 * outputSize.y));

						if (ray.d.z >= 0)
							bsdfSampleMatrix[outputSize.x - row - 1][col] += throughput*numberOfBins / (4 * N)*ray.d.z;
						else
							bsdfSampleMatrix_B[outputSize.x - row - 1][col] -= throughput*numberOfBins / (4 * N)*ray.d.z;
					}
#if defined(MTS_OPENMP)
				} //omp critical
#endif
			}

			std::cout << std::endl;
			std::cout << "Sample() Done! Save to " << argv[3] << std::endl;

			////// Save top Sample as .exr
			saveToFile(bsdfSampleMatrix, argv[3], outputSize);

			////// print out how many raypaths in different path length
			for (int i = 0; i < 9; i++) {
				std::cout << "NumOfPath = " << i << ": " << avgPathLength[i] << std::endl;
			}
			std::cout << "NumOfPath >= 9: " << avgPathLength[9] << std::endl;
		

		}

        // Eval
		{
			///// Our eval() disk
			std::vector< std::vector<Spectrum> > bsdfEvalMatrix(outputSize.y, std::vector<Spectrum>(outputSize.x));
			std::vector< std::vector<Spectrum> > bsdfEvalMatrix_B(outputSize.y, std::vector<Spectrum>(outputSize.x));

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

				Float x = (col + 0.5f)*2.0f / outputSize.y - 1;
				Float y = (row + 0.5f)*2.0f / outputSize.x - 1;
				Float z2 = x*x + y*y;
				if (z2 <= 1) {
					Float z = std::sqrt(1 - z2); // Only consider top hemisphere
					const Vector wo = Point(x, y, z) - Point(0.0f, 0.0f, 0.0f);

					const long long m = N / (outputSize.x*outputSize.y);
					for (int i = 0; i < m; i++) {
						Ray ray(ray_init);
						int depth;
#ifdef SCENE
						bsdfEval += Eval(scene, ray, samplers[tid], wo, depth,
							surfaceBSDF);
#else
						bsdfEval += Eval(heightfieldShape, ray, samplers[tid], wo, depth,
							surfaceBSDF);
#endif					
					}
					bsdfEval /= static_cast<Float>(m);
				}

				bsdfEvalMatrix[outputSize.x - row - 1][col] = bsdfEval;
			}

			std::cout << "Eval() Done! Save to " << argv[4] << std::endl;

			////// Save top Eval as .exr
			saveToFile(bsdfEvalMatrix, argv[4], outputSize);
		}


        std::cout<<"Done!"<<std::endl;
        return 0;
    }

    MTS_DECLARE_UTILITY()
};


MTS_EXPORT_UTILITY(traceRay, " ")
MTS_NAMESPACE_END
