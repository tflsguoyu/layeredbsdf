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
#include <chrono>

//#undef MTS_OPENMP

#if defined(MTS_OPENMP)
#   include <omp.h>
#endif

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

	int run(int argc, char **argv) {
		int nworker;
#if defined(MTS_OPENMP)
		nworker = omp_get_max_threads();
#else
		nworker = 1;
#endif

		ref_vector<Sampler> samplers(nworker);
		samplers[0] = static_cast<Sampler *>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), Properties("independent")));
		for (int i = 1; i < nworker; ++i) samplers[i] = samplers[0]->clone();

		// Incident ray
		const Point p = Point(0.0f, 0.0f, 0.0f);
		const Vector wi = normalize(Vector(static_cast<Float>(std::atof(argv[3])),
			static_cast<Float>(std::atof(argv[4])),
			static_cast<Float>(std::atof(argv[5]))));
		const RayDifferential ray_init = RayDifferential(p + wi * 10.0f, -wi, 0.0f);

		const Vector wo = normalize(Vector(static_cast<Float>(std::atof(argv[6])),
			static_cast<Float>(std::atof(argv[7])),
			static_cast<Float>(std::atof(argv[8]))));

		ref<Scene> scene = loadScene(argv[1]);
		scene->incRef();
		scene->initialize();

		const RayDifferential ray = RayDifferential(Point(0.0f, 0.0f, 1.0f), Vector(0.0f, 0.0f, -1.0f), 0.0f);
		Intersection its;
		scene->rayIntersect(ray, its);
		const BSDF* surfaceBSDF = its.getBSDF();
		
		long long M = std::atoll(argv[2]);

		Spectrum bsdfEval = Spectrum(0.0f);
	
		std::chrono::time_point<std::chrono::steady_clock> start;
		start = std::chrono::steady_clock::now();

#if defined(MTS_OPENMP)
#pragma omp parallel for
#endif
		for (int omp_i = 0; omp_i < M; ++omp_i) {
			int tid;
#if defined(MTS_OPENMP)
			tid = omp_get_thread_num();
#else
			tid = 0;
#endif
					
			BSDFSamplingRecord bRec(its, wi, wo);
			bRec.sampler = samplers[tid];
			bsdfEval += surfaceBSDF->eval(bRec);
				
		}
		bsdfEval /= static_cast<Float>(M);
		bsdfEval /= wo.z;


		std::cout << "wi:" << wi.toString() << endl;
		std::cout << "wo:" << wo.toString() << endl;
		std::cout << "bsdfEval:" << bsdfEval.toString() << endl;
		printf("Time: %.2lf secs)\n", std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count() / 1000.0);

		const Vector2i outputSize = Vector2i(1, 1);
		std::vector< std::vector<Spectrum> > bsdfEvalMatrix(1, std::vector<Spectrum>(1));
		bsdfEvalMatrix[0][0] = bsdfEval;
		saveToFile(bsdfEvalMatrix, "tmp.exr", outputSize);

		return 0;
	}

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(traceRay, " ")
MTS_NAMESPACE_END
