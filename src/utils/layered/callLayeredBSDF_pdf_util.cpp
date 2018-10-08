#include <fstream>
#include <vector>
#include <string>
#include <chrono>

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

const int flag_sample = 1;

MTS_NAMESPACE_BEGIN
class traceRay : public Utility {
public:

	void saveToFile(const std::vector< std::vector<Float> > &pdfMatrix, const char *filename, const Vector2i &outputSize) {

		ref<Bitmap> outputBitmap = new Bitmap(Bitmap::ELuminance, Bitmap::EFloat32, outputSize);
		float *outputData = outputBitmap->getFloat32Data();
		for (int row = 0; row < outputSize.x; row++) {
			for (int col = 0; col < outputSize.y; col++) {
				*outputData++ = (float)pdfMatrix[row][col];
			}
		}
		ref<FileStream> outputFile = new FileStream(filename, FileStream::ETruncReadWrite);
		outputBitmap->write(Bitmap::EOpenEXR, outputFile);
	}

	Spectrum Sample(RayDifferential &ray, Float &pdf, RadianceQueryRecord &rRec) {
		pdf = 0.0f;
		if (std::abs(ray.d.z) < Epsilon)
			return Spectrum(0.0f);

		const Scene *scene = rRec.scene;
		Intersection &its = rRec.its;

		if (!scene->rayIntersect(ray, its)) {
			return Spectrum(0.0f);
		}

		const BSDF *bsdf = its.getBSDF();
		BSDFSamplingRecord bRec(its, rRec.sampler);
		Spectrum sampleVal;
		sampleVal = bsdf->sample(bRec, rRec.nextSample2D());

		const Vector wo = its.toWorld(bRec.wo);
		ray = Ray(its.p, wo, 0.0f);

		return sampleVal;
	}

	Float Pdf(RayDifferential &ray, RadianceQueryRecord &rRec, const Vector wo_query) {
		if (std::abs(ray.d.z) < Epsilon || std::abs(wo_query.z) < Epsilon)
			return 0.0f;

		const Scene *scene = rRec.scene;
		Intersection &its = rRec.its;
		
		if (!scene->rayIntersect(ray, its)) {
			return 0.0f;
		}

		const BSDF *bsdf = its.getBSDF();

		BSDFSamplingRecord bRec_nee(its, its.toLocal(wo_query));
		bRec_nee.sampler = rRec.sampler;
		
		Float pdf = 0.0f;
		//Spectrum bsdfVal(0.0f);

		//bsdfVal = bsdf->eval(bRec_nee, pdf);
		pdf = bsdf->pdf(bRec_nee);

		//Spectrum evalVal(0.0), sampleVal(0.0);
		//Float evalPdf = 0.0, samplePdf = 0.0;
		//bsdf->evalAndSample(bRec_nee, evalVal, evalPdf, sampleVal, samplePdf, rRec.nextSample2D());
		
		return pdf;

	}


	int run(int argc, char **argv) {
		int nworker;
#if defined(MTS_OPENMP)
		nworker = omp_get_max_threads();
#else
		nworker = 1;
#endif
		ParameterMap param;
		param["IOR"] = argv[2];
		param["alphaT"] = argv[3];
		param["alphaB"] = argv[4];
		param["sigmaT"] = argv[5];
		param["hg"] = argv[6];
		ref<Scene> scene = loadScene(argv[1], param);
		scene->incRef();
		scene->initialize();

		ref_vector<Sampler> samplers(nworker);
		samplers[0] = static_cast<Sampler *>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), Properties("independent")));
		for (int i = 1; i < nworker; ++i) samplers[i] = samplers[0]->clone();

		const Point p = Point(static_cast<Float>(std::atof(argv[10])),
			static_cast<Float>(std::atof(argv[11])),
			static_cast<Float>(std::atof(argv[12])));
		const Vector wi = normalize(Vector(static_cast<Float>(std::atof(argv[13])),
			static_cast<Float>(std::atof(argv[14])),
			static_cast<Float>(std::atof(argv[15]))));
		const RayDifferential ray_init = RayDifferential(p + wi * 10.0f, -wi, 0.0f);

		const Vector2i outputSize = Vector2i(std::atoi(argv[9]), std::atoi(argv[9]));
		const int numberOfBins = outputSize.x * outputSize.y;
		const long long N = std::atoll(argv[7]);


		// Sample 1 
		if (flag_sample == 1) {
			std::cout << "Sample1() ... " << std::endl;
			// our sample()
			std::vector< std::vector<Float> > pdfSampleMatrix(outputSize.y, std::vector<Float>(outputSize.x, 0.0f));
			std::vector< std::vector<Float> > pdfSampleMatrix_B(outputSize.y, std::vector<Float>(outputSize.x, 0.0f));

			std::vector<int> avgPathLength(10, 0);
			long long countLocalAll = 0;
			long long countGlobalAll = 0;

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

				Float pdf;
				Spectrum pdfSample = Sample(ray, pdf, rRec);


#if defined(MTS_OPENMP)
#pragma omp critical
				{
#endif
					if (++tot % 5000000 == 0)
						// std::cout << "Sample count: " << tot/1000000 << "M of " << (N+1)/1000000 << "M" << std::endl;
						std::cout << "Sample count: " << tot / 1000000 << "M of " << (N + 1) / 1000000 << "M\r" << std::flush;

					if (depth >= 9) avgPathLength[9] ++;
					else avgPathLength[depth] ++;

					if (!pdfSample.isZero()) {
						int col, row;
						col = ray.d.x == 1.0f ? outputSize.x - 1 : static_cast<int>(std::floor((ray.d.x + 1.0) / 2.0 * outputSize.x));
						row = ray.d.y == 1.0f ? outputSize.y - 1 : static_cast<int>(std::floor((ray.d.y + 1.0) / 2.0 * outputSize.y));

						if (ray.d.z >= 0)
							pdfSampleMatrix[outputSize.x - row - 1][col] += 1 * static_cast<Float>(numberOfBins) / static_cast<Float>(4 * N)*ray.d.z;
						else
							pdfSampleMatrix_B[outputSize.x - row - 1][col] -= 1 * static_cast<Float>(numberOfBins) / static_cast<Float>(4 * N)*ray.d.z;
					}
#if defined(MTS_OPENMP)
				} //omp critical
#endif      
			}

			std::cout << std::endl;
			std::cout << "Sample() Done! Save to " << argv[16] << std::endl;

			////// Save top Sample as .exr
			saveToFile(pdfSampleMatrix, argv[16], outputSize);

			///// Save bottom Sample as .exr
			if (argc >= 14) {
				saveToFile(pdfSampleMatrix_B, argv[18], outputSize);
			}


			////// print out how many raypaths in different path length
			for (int i = 0; i < 9; i++) {
				std::cout << "NumOfPath = " << i << ": " << avgPathLength[i] << std::endl;
			}
			std::cout << "NumOfPath >= 9: " << avgPathLength[9] << std::endl;
		}

		else {
			std::cout << "Sample2() ... " << std::endl;
			// our sample()
			std::vector< std::vector<Float> > pdfSampleMatrix(outputSize.y, std::vector<Float>(outputSize.x, 0.0f));
			std::vector< std::vector<Float> > pdfSampleMatrix_B(outputSize.y, std::vector<Float>(outputSize.x, 0.0f));
			std::vector< std::vector<int> > pdfSampleMatrixCount(outputSize.y, std::vector<int>(outputSize.x, 0));
			std::vector< std::vector<int> > pdfSampleMatrixCount_B(outputSize.y, std::vector<int>(outputSize.x, 0));

			std::vector<int> avgPathLength(10, 0);
			long long countLocalAll = 0;
			long long countGlobalAll = 0;

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

				Float pdf;
				Spectrum pdfSample = Sample(ray, pdf, rRec);


#if defined(MTS_OPENMP)
#pragma omp critical
				{
#endif
					if (++tot % 5000000 == 0)
						// std::cout << "Sample count: " << tot/1000000 << "M of " << (N+1)/1000000 << "M" << std::endl;
						std::cout << "Sample count: " << tot / 1000000 << "M of " << (N + 1) / 1000000 << "M\r" << std::flush;

					if (depth >= 9) avgPathLength[9] ++;
					else avgPathLength[depth] ++;

					if (!pdfSample.isZero()) {
						int col, row;
						col = ray.d.x == 1.0f ? outputSize.x - 1 : static_cast<int>(std::floor((ray.d.x + 1.0) / 2.0 * outputSize.x));
						row = ray.d.y == 1.0f ? outputSize.y - 1 : static_cast<int>(std::floor((ray.d.y + 1.0) / 2.0 * outputSize.y));

						if (ray.d.z >= 0) {
							pdfSampleMatrix[outputSize.x - row - 1][col] += pdf;
							pdfSampleMatrixCount[outputSize.x - row - 1][col] ++;
						}
						else {
							pdfSampleMatrix_B[outputSize.x - row - 1][col] += pdf;
							pdfSampleMatrixCount_B[outputSize.x - row - 1][col] ++;
						}
					}
#if defined(MTS_OPENMP)
				} //omp critical
#endif      
			}

#if defined(MTS_OPENMP)
#pragma omp parallel for
#endif
			for (int omp_i = 0; omp_i < outputSize.x*outputSize.y; ++omp_i) {
				const int row = omp_i / outputSize.x, col = omp_i % outputSize.x;

				if (pdfSampleMatrixCount[row][col] != 0) {
					pdfSampleMatrix[row][col] /= pdfSampleMatrixCount[row][col];
				}
				if (pdfSampleMatrixCount_B[row][col] != 0) {
					pdfSampleMatrix_B[row][col] /= pdfSampleMatrixCount_B[row][col];
				}
			}

			std::cout << std::endl;
			std::cout << "Sample() Done! Save to " << argv[16] << std::endl;

			////// Save top Sample as .exr
			saveToFile(pdfSampleMatrix, argv[16], outputSize);

			///// Save bottom Sample as .exr
			if (argc >= 14) {
				saveToFile(pdfSampleMatrix_B, argv[18], outputSize);
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
			std::vector< std::vector<Float> > pdfEvalMatrix(outputSize.y, std::vector<Float>(outputSize.x, 0.0f));
			std::vector< std::vector<Float> > pdfEvalMatrix_B(outputSize.y, std::vector<Float>(outputSize.x, 0.0f));

			long long M = std::atoll(argv[8]);
			if (M == -1) {
				if (N >= numberOfBins)
					M = static_cast<long long>(N / numberOfBins);
				else
					M = 1;
			}
			std::cout << "Path per direction: " << M << std::endl;

			std::chrono::time_point<std::chrono::steady_clock> start;
			start = std::chrono::steady_clock::now();

#if defined(MTS_OPENMP)
#pragma omp parallel for
#endif
			for (int omp_i = 0; omp_i < outputSize.x*outputSize.y; ++omp_i) {
			//for (int omp_i = 100 * 200 + 100; omp_i <= 100*200 + 100; ++omp_i) {
					int tid;
#if defined(MTS_OPENMP)
				tid = omp_get_thread_num();
#else
				tid = 0;
#endif
				const int row = omp_i / outputSize.x, col = omp_i % outputSize.x;

				Float pdfEval = 0.0f;
				Float pdfEval_B = 0.0f;

				//Float x = 0.70710678f;
				//Float y = 0.0f;
				Float x = (col + 0.5f)*2.0f / outputSize.y - 1;
				Float y = (row + 0.5f)*2.0f / outputSize.x - 1;
				Float z2 = x * x + y * y;

				if (z2 <= 1) {
					// light goes out from top surface
					Float z = std::sqrt(1 - z2);
					const Vector wo = Point(x, y, z) - Point(0.0f, 0.0f, 0.0f);
					for (int i = 0; i < M; i++) {
						RayDifferential ray(ray_init);
						RadianceQueryRecord rRec(scene, samplers[tid]);
						rRec.type = RadianceQueryRecord::ERadiance;
						int depth = 0;
						Float val = Pdf(ray, rRec, wo);
						pdfEval += val;
						//cout << ">>> " << val << endl;
					}
					pdfEval /= static_cast<Float>(M);
					//cout << "eval: " << pdfEval << endl;

					if (argc >= 14) {
						// light goes out from bottom surface
						const Vector wo_B = Point(x, y, -z) - Point(0.0f, 0.0f, 0.0f);
						for (int i = 0; i < M; i++) {
							RayDifferential ray_B(ray_init);
							RadianceQueryRecord rRec_B(scene, samplers[tid]);
							rRec_B.type = RadianceQueryRecord::ERadiance;
							int depth = 0;
							pdfEval_B += Pdf(ray_B, rRec_B, wo_B);
						}
						pdfEval_B /= static_cast<Float>(M);
					}
				}

				pdfEvalMatrix[outputSize.x - row - 1][col] = pdfEval;

				if (argc >= 15) {
					pdfEvalMatrix_B[outputSize.x - row - 1][col] = pdfEval_B;
				}
			}

			std::cout << "Eval() Done! Save to " << argv[17] << std::endl;
			printf("Time: %.2lf secs)\n", std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count() / 1000.0);

			////// Save top Eval as .exr
			saveToFile(pdfEvalMatrix, argv[17], outputSize);

			if (argc >= 15) {
				saveToFile(pdfEvalMatrix_B, argv[19], outputSize);
			}

		}

		std::cout << "Done!" << std::endl;
		return 0;
	}

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(traceRay, " ")
MTS_NAMESPACE_END


