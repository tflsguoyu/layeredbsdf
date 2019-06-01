#include <vector>
#include <ctime>

#include <mitsuba/core/warp.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/ray.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/phase.h>

#include <mitsuba/render/texture.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>

#include <mitsuba/core/statistics.h>

#include <boost/math/special_functions/fpclassify.hpp>


MTS_NAMESPACE_BEGIN

class MultiLayeredBSDF : public BSDF {
public:
	MultiLayeredBSDF(const Properties &props) : BSDF(props) {
		cout << "##################################" << endl;
		m_maxDepth = props.getInteger("maxDepth", -1);
		m_MISenable = props.getBoolean("MIS", true);
		m_multiLayerSupport = props.getBoolean("multilayer", false);
		m_pdfMode = props.getString("pdf", "bidirStochTRT");
		m_stochPdfDepth = props.getInteger("stochPdfDepth", -1);
		m_pdfRepetitive = props.getInteger("pdfRepetitive", 1);
		m_diffusePdf = props.getFloat("diffusePdf", 0.0);
		m_bidirUseAnalog = props.getBoolean("bidirUseAnalog", false);
		m_bidir = props.getBoolean("bidir", true);
		m_maxSurvivalProb = props.getFloat("maxSurvivalProb", 1.0f);

		m_nbLayers = props.getInteger("nbLayers", 2);

		for (int l = 0; l < m_nbLayers-1; ++l) {
			std::string index = std::to_string(l);
			std::string name;

			name = std::string("aniso_") + index;
			m_flag_aniso.push_back(props.getBoolean(name, false));
			if (!m_flag_aniso[l]) std::cout << "[GY]: layer " << l + 1 << " of " << m_nbLayers << " :: medium :: Homogeneous Isotropic" << std::endl;
			else std::cout << "[GY]: layer " << l + 1 << " of " << m_nbLayers << " :: medium :: Homogeneous Anisotropic" << std::endl;

			name = std::string("sigmaT_") + index;
			m_spectrum_sigmaTs.push_back(props.getSpectrum(name, Spectrum(1.0)));
			if (!m_flag_aniso[l]) std::cout << "[GY]: layer " << l + 1 << " of " << m_nbLayers << " :: medium :: sigmaT(Spectrum)" << std::endl;

			m_flag_sigmaTs.push_back(false);

			name = std::string("density_") + index;
			m_float_densities.push_back(props.getFloat(name, 1.0));
			if (m_flag_aniso[l]) std::cout << "[GY]: layer " << l + 1 << " of " << m_nbLayers << " :: medium :: density(Float)" << std::endl;

			m_flag_densities.push_back(false);

			name = std::string("albedo_") + index;
			m_spectrum_albedos.push_back(props.getSpectrum(name, Spectrum(1.0)));
			std::cout << "[GY]: layer " << l + 1 << " of " << m_nbLayers << " :: medium :: albedo(Spectrum)" << std::endl;
			
			m_flag_albedos.push_back(false);

			name = std::string("orientation_") + index;
			m_vector_orientations.push_back(props.getVector(name, Vector(1.0,0.0,0.0)));
			if (m_flag_aniso[l]) std::cout << "[GY]: layer " << l + 1 << " of " << m_nbLayers << " :: medium :: orientation(Vector)" << std::endl;

			m_flag_orientations.push_back(false);

			name = std::string("normal_") + index;
			m_vector_normals.push_back(props.getVector(name, Vector(0.0, 0.0, 1.0)));

			m_flag_normals.push_back(false);
		}
		m_vector_normals.push_back(props.getVector(std::string("normal_") + std::to_string(m_nbLayers - 1), Vector(0.0, 0.0, 1.0)));
		m_flag_normals.push_back(false);

		m_bsdfs.resize(m_nbLayers);
		m_texture_normals.resize(m_nbLayers);

		m_texture_sigmaTs.resize(m_nbLayers - 1);
		m_texture_densities.resize(m_nbLayers - 1);
		m_texture_albedos.resize(m_nbLayers - 1);
		m_texture_orientations.resize(m_nbLayers - 1);

		m_phaseFunctions.resize(m_nbLayers-1);

		cout << "##################################" << endl;
	}

	MultiLayeredBSDF(Stream *stream, InstanceManager *manager)
		: BSDF(stream, manager) {
		Log(EError, "Not implemented.");
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Log(EError, "Not implemented.");
	}

	void configure() {
		cout << "##################################" << endl;

		if (m_bidir) {
			cout << "[GY]: Using bidir Eval" << endl;
		}
		else {
			cout << "[GY]: Using unidir Eval" << endl;
		}

		cout << "[GY]: Pdf mode: " << m_pdfMode << endl;
		if (m_pdfMode == "stoch" || m_pdfMode == "bidirStoch")
			cout << "[GY]: stochPdfDepth is: " << m_stochPdfDepth << endl;

		for (int l = 0; l < m_nbLayers; ++l) {
			m_bsdfs[l]->configure();
		}

		// Medium
		m_mediums.resize(NUM_COUNTERS, ref_vector<Medium>(m_nbLayers - 1));
		m_pdfMediums.resize(NUM_COUNTERS, ref_vector<Medium>(m_nbLayers - 1));

		Properties interiorIsoHomoProps("homogeneous");
		interiorIsoHomoProps.setSpectrum("sigmaT", Spectrum(0.0));
		interiorIsoHomoProps.setSpectrum("albedo", Spectrum(1.0));

		Properties interiorAniHomoProps("homogeneous_aniso");
		interiorAniHomoProps.setFloat("density", 0.0);
		interiorAniHomoProps.setSpectrum("albedo", Spectrum(1.0));
		interiorAniHomoProps.setVector("orientation", Vector(1.0, 0.0, 0.0));

		for (size_t i = 0; i < NUM_COUNTERS; ++i) {
			for (int l = 0; l < m_nbLayers - 1; ++l) {
				if (m_flag_aniso[l]) {
					m_mediums[i][l] = static_cast<Medium *> (PluginManager::getInstance()->
						createObject(MTS_CLASS(Medium), interiorAniHomoProps));
					m_mediums[i][l]->addChild(m_phaseFunctions[l]);
					m_mediums[i][l]->configure();

					m_pdfMediums[i][l] = static_cast<Medium *> (PluginManager::getInstance()->
						createObject(MTS_CLASS(Medium), interiorAniHomoProps));
					m_pdfMediums[i][l]->addChild(m_phaseFunctions[l]);
					m_pdfMediums[i][l]->configure();
				}
				else {
					m_mediums[i][l] = static_cast<Medium *> (PluginManager::getInstance()->
						createObject(MTS_CLASS(Medium), interiorIsoHomoProps));
					m_mediums[i][l]->addChild(m_phaseFunctions[l]);
					m_mediums[i][l]->configure();

					m_pdfMediums[i][l] = static_cast<Medium *> (PluginManager::getInstance()->
						createObject(MTS_CLASS(Medium), interiorIsoHomoProps));
					m_pdfMediums[i][l]->addChild(m_phaseFunctions[l]);
					m_pdfMediums[i][l]->configure();
				}
			}
		}


		// 
		size_t componentCount = m_bsdfs[m_nbLayers-1]->getComponentCount();
		m_components.reserve(componentCount);
		m_components.clear();
		for (int j = 0; j < m_bsdfs[m_nbLayers - 1]->getComponentCount(); ++j) {
			int componentType = m_bsdfs[m_nbLayers - 1]->getType(j);
			m_components.push_back(componentType);
		}

		BSDF::configure();

		if (BSDF::hasComponent(BSDF::ETransmission)) {
			m_eta = 1.0;
			for (int l = 0; l < m_nbLayers; ++l) m_eta *= m_bsdfs[l]->getEta();
			m_invEta = Float(1.0) / m_eta;
			cout << "[GY]: Tranparent layer" << endl;
			cout << "[GY]: eta = " << m_eta << endl;
		}
		else {
			cout << "[GY]: Non-tranparent layer" << endl;
		}

		cout << "[GY]: Configuration Done!" << endl;
		cout << "##################################" << endl;
	}

	Spectrum getDiffuseReflectance(const Intersection &its) const {
		Log(EError, "Not implemented.");
		return Spectrum(0.0);
	}

	Spectrum getSpecularReflectance(const Intersection &its) const {
		Log(EError, "Not implemented.");
		return Spectrum(0.0);
	}

	inline Float miWeight(Float pdfA, Float pdfB) const {
		if (pdfA == 0 && pdfB == 0)
			return 0.0;
		pdfA *= pdfA; pdfB *= pdfB;
		return pdfA / (pdfA + pdfB);
	}

	Normal getNormalFromTexture(const ref<Texture2D> normal_texture, const Point2 &uv) const {
		Normal normal;
		normal_texture->eval(uv).toLinearRGB(normal.x, normal.y, normal.z);
		for (int i = 0; i < 3; ++i) {
			//if (normal[i] > 1 || normal[i] < 0) {
			//	cout << "[GY]: Warning in MultiLayeredBSDF::getNormalFromTexture()" << endl;
			//	cout << "normal: " << normal.toString() << endl;
			//}
			normal[i] = Float(2.0) * normal[i] - Float(1.0);
			normal[i] = normal[i] > Float(1.0)-Epsilon ? Float(1.0) : normal[i];
			normal[i] = normal[i] < Float(-1.0)+Epsilon ? Float(-1.0) : normal[i];
		}
		return normalize(normal);
	}

	Vector getOrientationFromTexture(const ref<Texture2D> orien_texture, const Point2 &uv) const {
		Vector orien;
		orien_texture->eval(uv).toLinearRGB(orien.x, orien.y, orien.z);
		for (int i = 0; i < 3; ++i) {
			//if (normal[i] > 1 || normal[i] < 0) {
			//	cout << "[GY]: Warning in MultiLayeredBSDF::getNormalFromTexture()" << endl;
			//	cout << "normal: " << normal.toString() << endl;
			//}
			orien[i] = Float(2.0) * orien[i] - Float(1.0);
			orien[i] = orien[i] > Float(1.0) - Epsilon ? Float(1.0) : orien[i];
			orien[i] = orien[i] < Float(-1.0) + Epsilon ? Float(-1.0) : orien[i];
		}
		return normalize(orien);
	}

	void setParameters(const BSDFSamplingRecord &_bRec, std::vector<Frame> &frames,
		ref_vector<Medium> &mediums) const {

		Point2 uv = _bRec.its.uv;
		int thread_id = Thread::getID() & NUM_COUNTERS_MASK;

		// medium and phase function for specific position.     
		mediums = m_mediums[thread_id];

		for (int l = 0; l < m_nbLayers-1; ++l) {
			if (m_flag_aniso[l]) {
				Float density = m_flag_densities[l] ? m_texture_densities[l]->eval(uv)[0] * m_float_densities[l] : m_float_densities[l];
				Spectrum albedo = m_spectrum_albedos[l];
				if (m_flag_albedos[l]) albedo *= m_texture_albedos[l]->eval(uv);
				Vector orientation = m_flag_orientations[l] ? getOrientationFromTexture(m_texture_orientations[l], uv) : m_vector_orientations[l];
				mediums[l]->setMediumProp(density, albedo, orientation);
			}
			else {
				Spectrum sigmaT = m_spectrum_sigmaTs[l];
				if (m_flag_sigmaTs[l]) sigmaT *= m_texture_sigmaTs[l]->eval(uv);
				Spectrum albedo = m_spectrum_albedos[l];
				if (m_flag_albedos[l]) albedo *= m_texture_albedos[l]->eval(uv);
				mediums[l]->setSigmaAST(sigmaT, albedo);
			}
		}

		// Shading normal   
		for (int l = 0; l < m_nbLayers; ++l) {
			const Normal normal = m_flag_normals[l] ? getNormalFromTexture(m_texture_normals[l], uv) : m_vector_normals[l];
			frames.push_back(Frame(normalize(normal)));
		}

	}

	void setParametersPdf(const BSDFSamplingRecord &_bRec,
		ref_vector<Medium> &mediums) const {

		Point2 uv = _bRec.its.uv;
		int thread_id = Thread::getID() & NUM_COUNTERS_MASK;

		// medium and phase function for specific position.     
		mediums = m_pdfMediums[thread_id];

		for (int l = 0; l < m_nbLayers - 1; ++l) {
			if (m_flag_aniso[l]) {
				Float density = m_flag_densities[l] ? m_texture_densities[l]->eval(uv)[0] * m_float_densities[l] : m_float_densities[l];
				Spectrum albedo = Spectrum(0.0);
				Vector orientation = m_flag_orientations[l] ? getOrientationFromTexture(m_texture_orientations[l], uv) : m_vector_orientations[l];
				mediums[l]->setMediumProp(density, albedo, orientation);
			}
			else {
				Spectrum sigmaT = m_spectrum_sigmaTs[l];
				if (m_flag_sigmaTs[l]) sigmaT *= m_texture_sigmaTs[l]->eval(uv);
				Spectrum albedo = Spectrum(0.0);
				mediums[l]->setSigmaAST(sigmaT, albedo);
			}
		}

	}

	bool rayIntersect(const Ray &ray, Intersection &its) const {

		its.p.x = std::numeric_limits<Float>::infinity();
		its.p.y = std::numeric_limits<Float>::infinity();
		its.p.z = std::numeric_limits<Float>::infinity();
		its.t = std::numeric_limits<Float>::infinity();
		its.wi = normalize(-ray.d);

		its.geoFrame = Frame(Normal(0.0, 0.0, 1.0));
		its.shFrame = its.geoFrame;

		if (std::abs(ray.d.z) < Epsilon)
			return false;

		bool isUp = ray.d.z > 0 ? true : false;
		
		Float z = ray(Epsilon).z;

		Float z_bot = std::floor(z);
		Float z_top = z_bot + 1;

		if ((z_top > 0 && isUp) || (z_bot < -(m_nbLayers - 1) && !isUp)) {
			return false;
		}

		its.p.z = isUp ? z_top : z_bot;
		if (its.p.z > Epsilon || its.p.z < -(m_nbLayers - 1) - Epsilon) {
			cout << "[GY]: Warning in rayIntersect()" << endl;
			return false;
		}

		its.p.x = ray.o.x + ray.d.x * ((its.p.z - ray.o.z) / ray.d.z);
		its.p.y = ray.o.y + ray.d.y * ((its.p.z - ray.o.z) / ray.d.z);
		its.t = (its.p - ray.o).length();

		return true;
	}

	bool rayIntersectAndLookForEmitter(const Ray &ray, Intersection &its, const bool flag_type, bool &isEmitter) const {

		isEmitter = false;

		its.p.x = std::numeric_limits<Float>::infinity();
		its.p.y = std::numeric_limits<Float>::infinity();
		its.p.z = std::numeric_limits<Float>::infinity();
		its.t = std::numeric_limits<Float>::infinity();
		its.wi = normalize(-ray.d);

		its.geoFrame = Frame(Normal(0.0, 0.0, 1.0));
		its.shFrame = its.geoFrame;

		if (std::abs(ray.d.z) < Epsilon)
			return false;

		bool isUp = ray.d.z > 0 ? true : false;

		Float z = ray(Epsilon).z;

		Float z_bot = std::floor(z);
		Float z_top = z_bot + 1;

		if ((z_top > 0 && isUp) || (z_bot < -(m_nbLayers - 1) && !isUp))
			return false;

		its.p.z = isUp ? z_top : z_bot;
		if (its.p.z > Epsilon || its.p.z < -(m_nbLayers - 1) - Epsilon) {
			cout << "[GY]: Warning in rayIntersect()" << endl;
			return false;
		}

		its.p.x = ray.o.x + ray.d.x * ((its.p.z - ray.o.z) / ray.d.z);
		its.p.y = ray.o.y + ray.d.y * ((its.p.z - ray.o.z) / ray.d.z);
		its.t = (its.p - ray.o).length();

		// 
		int curLayer = -math::roundToInt(its.p.z);
		if ((flag_type && curLayer == 0) || (!flag_type && curLayer == (m_nbLayers-1))) {
			isEmitter = true;
		}

		return true;

	}

	Spectrum sampleRefraction(const BSDFSamplingRecord &_bRec, const BSDF *bsdf, const Frame &frame,
		const Vector &wo_query, Vector &wo_refract, Float &pdf_refract, Float &pdf_refract_re) const {

		Spectrum weight;
		wo_refract = Vector(0.0);

		if (bsdf->getType() & BSDF::ENull) {
			wo_refract = -wo_query;
			pdf_refract = pdf_refract_re = 1.0;
			return Spectrum(1.0);
		}

		if (wo_query.z * frame.toLocal(wo_query).z <= 0) {
			pdf_refract = pdf_refract_re = 0.0;
			return Spectrum(0.0);
		}

		BSDFSamplingRecord bRec(_bRec);
		bRec.mode = EImportance;
		bRec.wi = frame.toLocal(wo_query);

		if (!m_multiLayerSupport) bRec.typeMask = BSDF::ETransmission;
		weight = bsdf->sample(bRec, pdf_refract, _bRec.sampler->next2D());
		BSDFSamplingRecord bRec_re(bRec);
		bRec_re.reverse();
		pdf_refract_re = bsdf->pdf(bRec_re);

		if (m_multiLayerSupport) {
			if (bRec.wi.z * bRec.wo.z >= 0) {
				pdf_refract = pdf_refract_re = 0.0;
				return Spectrum(0.0);
			}
		}

		wo_refract = frame.toWorld(bRec.wo);
		if (bRec.wo.z * wo_refract.z <= 0) {
			pdf_refract = pdf_refract_re = 0.0;
			return Spectrum(0.0);
		}

		if (!weight.isZero()) {
			Float Jacobian = std::abs(bRec.wo.z / bRec.wi.z) * (bRec.eta * bRec.eta);
			weight /= Jacobian; // Jacobian!
								//pdf_refract /= (bRec.eta * bRec.eta);
		}

		return weight;
	}

	Spectrum evaluateRefraction(const BSDFSamplingRecord &_bRec, const BSDF *bsdf, const Frame &frame,
		const Vector &wo_query, const Vector &wo_refract, Float &pdf_refract, Float &pdf_refract_re) const {

		Spectrum weight;

		if (bsdf->getType() & BSDF::ENull) {
			pdf_refract = pdf_refract_re = 1.0;
			return Spectrum(1.0);
		}

		if (wo_query.z * frame.toLocal(wo_query).z <= 0) {
			pdf_refract = pdf_refract_re = 0.0;
			return Spectrum(0.0);
		}
		if (wo_refract.z * frame.toLocal(wo_refract).z <= 0) {
			pdf_refract = pdf_refract_re = 0.0;
			return Spectrum(0.0);
		}

		BSDFSamplingRecord bRec(_bRec);
		bRec.mode = EImportance;
		bRec.wi = frame.toLocal(wo_query);
		bRec.wo = frame.toLocal(wo_refract);

		if (!m_multiLayerSupport) bRec.typeMask = BSDF::ETransmission;

		weight = bsdf->eval(bRec);
		pdf_refract = bsdf->pdf(bRec);
		BSDFSamplingRecord bRec_re(bRec);
		bRec_re.reverse();
		pdf_refract_re = bsdf->pdf(bRec_re);

		Float eta = bsdf->getEta();
		eta = bRec.wo.z < 0 ? eta : Float(1.0) / eta;

		if (!weight.isZero()) {
			Float Jacobian = std::abs(bRec.wo.z / bRec.wi.z) * (eta*eta);
			weight /= Jacobian;
			//pdf_refract /= (eta*eta); // Jacobian!
		}

		return weight;
	}

	struct PathInfo {
		Point p;
		Vector wi, wo;
		int topCounter, bottomCounter;
		bool surf;
		int layerID;
		Spectrum thru0, thru1;
		Vector2 vpdf, epdf;
		Float pSurvival;
	};

	static void printPath(const std::vector<PathInfo> &paths) {
		for (const auto &path : paths) {
			cout << path.p.toString() << ' ' << path.wi.toString() << ' ' << path.wo.toString() << '\n'
				 << path.surf << ' ' << path.layerID << '\n'
				 << path.thru0.toString() << ' ' << path.thru1.toString() << '\n'
				 << path.vpdf.toString() << ' ' << path.epdf.toString() << "\n==========" << endl;
		}
	}

	static void printVector(const std::vector<Float> &data) {
		for (const auto &value : data)
			cout << value << ' ';
		cout << endl;
	}

	Spectrum generatePath(BSDFSamplingRecord &_bRec,
		const std::vector<Frame> &frames, const ref_vector<Medium> &mediums, const int maxDepth,
		std::vector<PathInfo> &path, std::vector<Float> &ratio, std::vector<Float> &ratioPdf, bool flag_backward, bool flag_bidir) const {

		Assert(_bRec.sampler);
		Sampler *sampler = _bRec.sampler;

		// Path tracing
		bool flag_incidentDir = _bRec.wi.z > 0;

		Point originPoint = flag_incidentDir ? Point(0.0, 0.0, 0.0) : Point(0.0, 0.0, Float(1-m_nbLayers));
		Ray ray(originPoint + _bRec.wi, -_bRec.wi, 0.0);

		bool flag_medium = false;
		Intersection its;
		MediumSamplingRecord mRec;

		rayIntersect(ray, its);
		Spectrum throughput(1.0);

		int topCounter = 0;
		int bottomCounter = 0;
		int depth = 0;
		int curLayer;
		while (depth < maxDepth || maxDepth < 0) {
			PathInfo path_this;

			curLayer = -math::ceilToInt(ray(Epsilon).z);
			if (curLayer > m_nbLayers - 2) --curLayer;
			if (flag_medium && mediums[curLayer]->sampleDistance(Ray(ray, 0, its.t), mRec, sampler)) {
				if (mRec.p.z > Epsilon || mRec.p.z < -(m_nbLayers - 1)-Epsilon)
					cout << "[GY]: Warning in BSDF::multilayeredBSDF::generatePath()" << endl;

				if (curLayer < 0 || curLayer > m_nbLayers - 2)
					cout << "[GY]: Warning in BSDF::multilayeredBSDF::generatePath()" << endl;

				const Medium* medium = mediums[curLayer].get();
				const PhaseFunction* phase = medium->getPhaseFunction();

				path_this.layerID = curLayer;
				path_this.surf = false;
				path_this.p = mRec.p;
				path_this.wi = -ray.d;
				path_this.topCounter = 0;
				path_this.bottomCounter = 0;

				Spectrum albedo = mRec.sigmaS * mRec.transmittance / mRec.pdfSuccess;
				
				Float pSurvival;
				if (m_bidir && m_bidirUseAnalog) {
					pSurvival = std::min(albedo.max(), m_maxSurvivalProb);
					if (sampler->next1D() > pSurvival) {
						throughput = Spectrum(0.0);
						break;
					}
					path_this.thru0 = throughput *= albedo / pSurvival;
					path_this.pSurvival = pSurvival;
				}
				else {
					path_this.thru0 = throughput *= albedo;
					path_this.pSurvival = pSurvival = 1.0;
				}

				path_this.epdf[0] = mRec.pdfSuccess / std::abs(ray.d.z);
				if (path_this.epdf[0] < Epsilon) {
					throughput = Spectrum(0.0);
					break;					
				}

				path_this.epdf[1] = path.back().surf ? mRec.pdfFailure : mRec.pdfSuccessRev / std::abs(ray.d.z);

				PhaseFunctionSamplingRecord pRec(mRec, -ray.d);
				Float phaseVal = phase->sample(pRec, sampler);
				if (std::abs(phaseVal) < Epsilon) {
					throughput = Spectrum(0.0);
					break;
				}
				throughput *= phaseVal;
				path_this.thru1 = throughput;

				path_this.wo = pRec.wo;

				path_this.vpdf[0] = phase->pdf(pRec);
				if (path_this.vpdf[0] < Epsilon) {
					throughput = Spectrum(0.0);
					break;					
				}

				PhaseFunctionSamplingRecord pRec_reverse(pRec);
				pRec_reverse.reverse();
				path_this.vpdf[1] = phase->pdf(pRec_reverse);
				path_this.vpdf *= pSurvival;

				// Trace a ray
				ray = Ray(path_this.p, path_this.wo, 0.0);
				ray.mint = 0.0;

				path.push_back(path_this);

				if (!rayIntersect(ray, its)) {
					throughput = Spectrum(0.0);
					break;
				}
			}
			else {
				if (flag_medium) {
					throughput *= mRec.transmittance / mRec.pdfFailure;

					path_this.epdf[0] = mRec.pdfFailure;
					if (path_this.epdf[0] < Epsilon) {
						throughput = Spectrum(0.0);
						break;					
					}

					path_this.epdf[1] = path.back().surf ? mRec.pdfFailure : mRec.pdfSuccessRev / std::abs(ray.d.z);
				}

				if (!its.isValid()) {
					break;
				}
			
				curLayer = -math::roundToInt(its.p.z);
				
				if (curLayer == 0) {
					++topCounter;
				}
				else if (curLayer == (m_nbLayers - 1)) {
					++bottomCounter;
				}
				else
					;
			
				if (curLayer >= m_nbLayers)
					cout << "[GY]: Warning that current layer exceed max layers" << endl;

				const BSDF *bsdf = m_bsdfs[curLayer].get();
				const Frame frame = frames[curLayer];

				path_this.p = its.p;
				path_this.wi = its.wi;
				path_this.topCounter = topCounter;
				path_this.bottomCounter = bottomCounter;

				if (path_this.wi.z * frame.toLocal(path_this.wi).z <= 0) {
					throughput = Spectrum(0.0);
					break;
				}

				if ( m_bidir && m_bidirUseAnalog ) {
					if ( sampler->next1D() > m_maxSurvivalProb ) {
						throughput = Spectrum(0.0);
						break;
					}
					throughput /= m_maxSurvivalProb;
					path_this.pSurvival = m_maxSurvivalProb;
				}
				else
					path_this.pSurvival = 1.0f;

				BSDFSamplingRecord bRec(_bRec);
				bRec.mode = EImportance;
				bRec.wi = frame.toLocal(path_this.wi);
				Spectrum bsdfVal = bsdf->sample(bRec, sampler->next2D());
				if (bsdfVal.isZero()) {
					throughput = Spectrum(0.0);
					break;
				}

				const Vector wo = frame.toWorld(bRec.wo);

				if (bRec.wo.z * wo.z <= 0) {
					throughput = Spectrum(0.0);
					break;
				}

				path_this.wo = wo;

				if (flag_backward)
					bsdfVal *= std::abs((bRec.wi.z / bRec.wo.z)*(path_this.wo.z / path_this.wi.z));

				path_this.layerID = curLayer;
				path_this.surf = true;
				
				path_this.thru0 = throughput;
				if (depth == 0) {
					path_this.epdf[0] = 0.0;
					path_this.epdf[1] = 0.0;
				}

				throughput *= bsdfVal;
				path_this.thru1 = throughput;

				path_this.vpdf[0] = bsdf->pdf(bRec);
				if (path_this.vpdf[0] < Epsilon) {
					throughput = Spectrum(0.0);
					break;					
				}
				BSDFSamplingRecord bRec_reverse(bRec);
				bRec_reverse.reverse();
				path_this.vpdf[1] = bsdf->pdf(bRec_reverse);

				if ((curLayer == 0 || curLayer == (m_nbLayers-1)) && (path_this.wi.z * path_this.wo.z < 0)) {
					flag_medium = !flag_medium;
				}
				ray = Ray(path_this.p, path_this.wo, 0.0);

				path.push_back(path_this);

				if (!rayIntersect(ray, its)) {
					if (flag_medium) {
						throughput = Spectrum(0.0);
						break;
					}
				}

			}
			depth++;
		}
		_bRec.wo = ray.d;

		if (_bRec.wi.z * _bRec.wo.z >= 0)
			_bRec.eta = 1.0;
		else
			_bRec.eta = _bRec.wo.z < 0 ? m_eta : m_invEta;

		ratio.resize(path.size());
		ratioPdf.resize(path.size());
		if (flag_bidir && !path.empty()) {
			ratio[0] = 0.0;
			ratioPdf[0] = path[0].vpdf[1] / path[0].vpdf[0];
			for (size_t i = 1; i < path.size(); ++i) {
				Float r = Float(1.0) / path[i].vpdf[1] + Float(1.0) / path[i - 1].vpdf[0];
				Float r1 = path[i].vpdf[1] / path[i].vpdf[0];

				ratio[i] = r1 * (r + path[i].epdf[1] * ratio[i - 1]) / path[i].epdf[0];
				ratioPdf[i] = ratioPdf[i - 1] * (path[i].vpdf[1] / path[i].vpdf[0]) *(path[i].epdf[1] / path[i].epdf[0]);
			}
		}

		return throughput;
	}

	void unidirEvaluation(const BSDFSamplingRecord &_bRec,
		const std::vector<Frame> &frames, const ref_vector<Medium> &mediums,
		const std::vector<PathInfo> &path, const int mode, Spectrum &_val, Float &_pdf) const {
		
		Assert(_bRec.sampler);
		Sampler *sampler = _bRec.sampler;

		bool flag_type = _bRec.wi.z * _bRec.wo.z > 0; // 1:reflection  0:transmission
		bool flag_incidentDir = Frame::cosTheta(_bRec.wi) > 0;
		bool flag_neeDir = (flag_incidentDir && flag_type) || (!flag_incidentDir && !flag_type); // 1. go out from top  0: go out from bottom

		Intersection its;
		MediumSamplingRecord mRec;
		Spectrum Li(0.0);
		Float pdf = 0.0;

		if (path.empty()) {
			//cout << "[GY]: Warning in layeredBSDF::evaluatePdf1" << endl;
			_val = Spectrum(0.0);
			_pdf = 0.0;
			return;
		}

		int curLayer;
		size_t depth = 0;
		if (mode == 1) depth = path.size();
		if (mode == 2) depth = m_stochPdfDepth < 0 ? path.size() : std::min(path.size(), size_t(m_stochPdfDepth));

		for (size_t i = 0; i < depth; ++i) {
			//cout << path[i].p.z << "|" << path[i].wi.z << "|" << path[i].wo.z << "|" << path[i].surf << "|" << path[i].layerID << endl;
			if (!path[i].surf) {// medium
				if ((flag_neeDir && path[i].layerID == 0) || (!flag_neeDir && path[i].layerID == (m_nbLayers-2))) {
					curLayer = -math::ceilToInt(path[i].p.z);
					int curLayerTest = path[i].layerID;
					if (curLayer != curLayerTest) {
						cout << "[GY]: Warning curLayer not match" << endl;
						curLayer = curLayerTest;
					}

					const Medium *medium = mediums[curLayer].get();
					const PhaseFunction *phase = medium->getPhaseFunction();
					
					// Direct
					Vector wo_refract;
					Float refractPdf, refractPdf_re;

					const BSDF *bsdf_nee = flag_neeDir ? m_bsdfs[0].get() : m_bsdfs[m_nbLayers - 1].get();
					const Frame frame_nee = flag_neeDir ? frames[0] : frames[m_nbLayers - 1];

					Spectrum throughput_refract = sampleRefraction(_bRec, bsdf_nee, frame_nee,
						_bRec.wo, wo_refract, refractPdf, refractPdf_re);

					if (!throughput_refract.isZero()) {
						Intersection its_nee;
						Ray ray_nee = Ray(path[i].p, -wo_refract, 0.0);
						ray_nee.mint = 0.0;
						if (!rayIntersect(ray_nee, its_nee)) {
							break;
						}
						Spectrum value(0.0);
						if (mode == 1) value = medium->evalTransmittance(Ray(ray_nee, 0.0, its_nee.t), sampler);
						medium->eval(Ray(its_nee.p, wo_refract, 0.0, its_nee.t, 0.0), mRec);

						if (!value.isZero()) {
							PhaseFunctionSamplingRecord pRec(mRec, path[i].wi, -wo_refract);
							Float phaseVal;
							if (mode == 1) phaseVal = phase->eval(pRec);
							const Float phasePdf = phase->pdf(pRec);
							const Float weight = miWeight(refractPdf, phasePdf);
							if (m_MISenable) {
								if (mode == 1) Li += path[i].thru0 * phaseVal * value * throughput_refract * weight;
								if (mode == 2) pdf += phasePdf * mRec.pdfFailure * refractPdf_re / refractPdf * weight;
							}
							else {
								if (mode == 1) Li += path[i].thru0 * phaseVal * value * throughput_refract;
								if (mode == 2) pdf += phasePdf * mRec.pdfFailure * refractPdf_re / refractPdf;
							}
						}
					}

					// Indirect
					if (m_MISenable) {

						Ray ray = Ray(path[i].p, path[i].wo, 0.0);
						ray.mint = 0.0;
						bool isEmitter = false;
						if (!rayIntersectAndLookForEmitter(ray, its, flag_neeDir, isEmitter)) {
							//cout << "[GY]: Warning in layeredBSDF::evaluatePdf4" << endl;
							break;
						}

						if (isEmitter) {
							const BSDF *bsdf = flag_neeDir ? m_bsdfs[0].get() : m_bsdfs[m_nbLayers - 1].get();
							const Frame frame = flag_neeDir ? frames[0] : frames[m_nbLayers-1];

							Spectrum refractVal = evaluateRefraction(_bRec, bsdf, frame,
								_bRec.wo, its.wi, refractPdf, refractPdf_re);

							if (!refractVal.isZero()) {
								Spectrum value;
								if (mode == 1) value = medium->evalTransmittance(Ray(ray, 0.0, its.t), sampler);
								medium->eval(Ray(ray.o, ray.d, 0.0, its.t, 0.0), mRec);

								const Float weight = miWeight(path[i].vpdf[0], refractPdf);

								if (mode == 1) Li += path[i].thru1 * value * refractVal * weight;
								if (mode == 2) pdf += refractPdf_re * mRec.pdfFailure * weight;
							}
						}
					}
				}
			}
			else {// surface	
				curLayer = path[i].layerID;
				const BSDF *bsdf = m_bsdfs[curLayer].get();
				const Frame frame = frames[curLayer];

				// Direct
				if ((flag_incidentDir && flag_type && curLayer == 0 && path[i].topCounter == 1) ||
					(!flag_incidentDir && flag_type && curLayer == (m_nbLayers - 1) && path[i].bottomCounter == 1)) {

					BSDFSamplingRecord bRec(_bRec);
					bRec.wi = frame.toLocal(path[i].wi);
					bRec.wo = frame.toLocal(_bRec.wo);
					if (mode == 1) Li += path[i].thru0 * bsdf->eval(bRec);
					if (mode == 2) pdf += bsdf->pdf(bRec);
				}

				if ((flag_incidentDir && flag_type && curLayer == 1) ||
					(flag_incidentDir && !flag_type && curLayer == (m_nbLayers - 2)) ||
					(!flag_incidentDir && flag_type && curLayer == (m_nbLayers - 2)) ||
					(!flag_incidentDir && !flag_type && curLayer == 1)) {

					Vector wo_refract;
					Float refractPdf, refractPdf_re;

					const BSDF *bsdf_nee = flag_neeDir ? m_bsdfs[0].get() : m_bsdfs[m_nbLayers - 1].get();
					const Frame frame_nee = flag_neeDir ? frames[0] : frames[m_nbLayers - 1];
					const Medium *medium = flag_neeDir ? mediums[0].get() : mediums[m_nbLayers - 2].get();
					
					Spectrum throughput_refract = sampleRefraction(_bRec, bsdf_nee, frame_nee,
						_bRec.wo, wo_refract, refractPdf, refractPdf_re);

					if (!throughput_refract.isZero()) {
						Intersection its_nee;
						Ray ray_nee = Ray(path[i].p, -wo_refract, 0.0);
						if (!rayIntersect(ray_nee, its_nee)) {
							cout << "[GY]: Warning in layeredBSDF::evaluatePdf6" << endl;
							break;
						}
						Spectrum value(0.0);
						if (mode == 1) value = medium->evalTransmittance(Ray(ray_nee, 0.0, its_nee.t), sampler);
						if (mode == 2) medium->eval(Ray(its_nee.p, wo_refract, 0.0, its_nee.t, 0.0), mRec);

						if (!value.isZero()) {
							BSDFSamplingRecord bRec(_bRec);
							bRec.mode = EImportance;
							bRec.wi = frame.toLocal(path[i].wi);
							bRec.wo = frame.toLocal(-wo_refract);
							Spectrum bsdfVal;
							if (mode == 1) bsdfVal = bsdf->eval(bRec);
							const Float bsdfPdf = bsdf->pdf(bRec);
							const Float weight = miWeight(refractPdf, bsdfPdf);
							if (m_MISenable) {
								if (mode == 1) Li += path[i].thru0 * bsdfVal * value * throughput_refract * weight;
								if (mode == 2) pdf += bsdfPdf * mRec.pdfFailure * refractPdf_re / refractPdf * weight;
							}
							else {
								if (mode == 1) Li += path[i].thru0 * bsdfVal * value * throughput_refract;
								if (mode == 2) pdf += bsdfPdf * mRec.pdfFailure * refractPdf_re / refractPdf;
							}
						}
					}

					// Indirect
					if (m_MISenable) {

						Ray ray = Ray(path[i].p, path[i].wo, 0.0);

						bool isEmitter = false;
						if (!rayIntersectAndLookForEmitter(ray, its, flag_neeDir, isEmitter)) {
							// cout << "[GY]: Warning in layeredBSDF::evaluatePdf8" << endl;
							break;
						}

						if (isEmitter) {
							const BSDF *bsdf = flag_neeDir ? m_bsdfs[0].get() : m_bsdfs[m_nbLayers - 1].get();
							const Frame frame = flag_neeDir ? frames[0] : frames[m_nbLayers - 1];

							Spectrum refractVal = evaluateRefraction(_bRec, bsdf, frame,
								_bRec.wo, its.wi, refractPdf, refractPdf_re);

							if (!refractVal.isZero()) {
								Spectrum value;
								if (mode == 1) value = medium->evalTransmittance(Ray(ray, 0.0, its.t), sampler);
								if (mode == 2) medium->eval(Ray(ray.o, ray.d, 0.0, its.t, 0.0), mRec);

								const Float weight = miWeight(path[i].vpdf[0], refractPdf);
								if (mode == 1) Li += path[i].thru1 * value * refractVal * weight;
								if (mode == 2) pdf += refractPdf_re * mRec.pdfFailure * weight;
							}
						}
					}
				}
			}
		}

		if (mode == 1) {
			_val = Li;
			_pdf = 0.0;
		}
		if (mode == 2) {
			_val = Spectrum(0.0);
			_pdf = pdf + m_diffusePdf;
		}
	}

	void bidirEvaluation(const BSDFSamplingRecord &_bRec,
		const std::vector<Frame> &frames, const ref_vector<Medium> &mediums,
		const std::vector<PathInfo> path_L, const std::vector<Float> ratio_L, const std::vector<Float> ratioPdf_L,
		const std::vector<PathInfo> path_R, const std::vector<Float> ratio_R, const std::vector<Float> ratioPdf_R,
		const int mode, Spectrum &_val, Float &_pdf) const {

		Assert(_bRec.sampler);

		// Bi-directional start
		bool flag_type = _bRec.wi.z * _bRec.wo.z > 0; // 1:reflection  0:transmission
		bool flag_incidentDir = _bRec.wi.z > 0; // 1: from top surface  0: from bottom surface

		MediumSamplingRecord mRec;
		PhaseFunctionSamplingRecord pRec(mRec, Vector(0.0));

		BSDFSamplingRecord bRec(_bRec);
		BSDFSamplingRecord bRec_L(bRec);
		bRec_L.mode = EImportance;
		BSDFSamplingRecord bRec_R(bRec);
		bRec_R.mode = EImportance;

		_val = Spectrum(0.0);
		_pdf = 0.0;

		const Frame frame_wo = _bRec.wo.z > 0 ? frames[0] : frames[m_nbLayers - 1];
		if (_bRec.wo.z * frame_wo.toLocal(_bRec.wo).z <= 0) return;

		Spectrum Li0(0.0);
		Spectrum Li(0.0);
		Float pdf0 = 0.0;
		Float pdf = 0.0;

		if (flag_type) {
			bRec.wo = frame_wo.toLocal(_bRec.wo);
			bRec.wi = frame_wo.toLocal(_bRec.wi);
			if (mode == 1) Li0 += flag_incidentDir ? m_bsdfs[0]->eval(bRec) : m_bsdfs[m_nbLayers - 1]->eval(bRec);
			if (mode == 2) pdf0 += flag_incidentDir ? m_bsdfs[0]->pdf(bRec) : m_bsdfs[m_nbLayers - 1]->pdf(bRec);
		}

		//if (!flag_type) {
		//	cout << "PATHLLLLLLLLLLLL:" << endl;
		//	for (size_t i = 0; i < path_L.size(); ++i)
		//		cout << path_L[i].p.z << "|" << path_L[i].wi.z << "|" << path_L[i].wo.z << "|" << path_L[i].surf << "|" << path_L[i].layerID << endl;
		//	cout << "PATHRRRRRRRRRRRR:" << endl;
		//	for (size_t i = 0; i < path_R.size(); ++i)
		//		cout << path_R[i].p.z << "|" << path_R[i].wi.z << "|" << path_R[i].wo.z << "|" << path_R[i].surf << "|" << path_R[i].layerID << endl;
		//}

		for (size_t i = 0; i < (mode == 1 || m_stochPdfDepth < 0 ? path_L.size() : std::min(path_L.size(), size_t(m_stochPdfDepth))); ++i) {
			for (size_t j = 0; j < (mode == 1 || m_stochPdfDepth < 0 ? path_R.size() : std::min(size_t(m_stochPdfDepth) - i, path_R.size())); ++j) {
	
				bool validConnection = false;
				int id_L = path_L[i].layerID;
				int id_R = path_R[j].layerID;
				int id_connectMedium;

				Float z1 = path_L[i].p.z;
				Float z2 = path_R[j].p.z;
				if (path_L[i].surf) {
					if (z2 > z1 - 1 - Epsilon && z2 < z1 + 1 + Epsilon && std::abs(z1 - z2) > Epsilon) {
						validConnection = true;
						if (z1 > z2)
							id_connectMedium = id_L;
						else
							id_connectMedium = id_L - 1;
					}
				}
				else {
					if (z2 > std::floor(z1) - Epsilon && z2 < std::ceil(z1) + Epsilon) {
						validConnection = true;
						id_connectMedium = id_L;
					}
				}

				if (validConnection) {
					Spectrum f(0.0), funcVal(0.0);
					Float w, t, f_pdf;
					Vector2 pdf_LL, pdf_RR;
					Frame frame;

					// sample from left
					frame = frames[id_R];

					if (!path_R[j].surf || (path_R[j].wi.z*frame.toLocal(path_R[j].wi).z > 0 && -path_L[i].wo.z*frame.toLocal(-path_L[i].wo).z > 0)) {
						t = (path_R[j].p.z - path_L[i].p.z) / path_L[i].wo.z;
						if (t > Epsilon) {
							const Medium *medium = mediums[id_connectMedium].get();
							medium->eval(Ray(Point(0.0, 0.0, path_L[i].p.z), path_L[i].wo, 0.0, t, 0.0), mRec);
							mRec.pdfSuccess /= std::abs(path_L[i].wo.z);
							mRec.pdfSuccessRev /= std::abs(path_L[i].wo.z);

							if (path_R[j].surf) {
								const BSDF *bsdf = m_bsdfs[id_R].get();
								bRec_R.wi = frame.toLocal(path_R[j].wi);
								bRec_R.wo = frame.toLocal(-path_L[i].wo);
								if (mode == 1) {
									funcVal = bsdf->eval(bRec_R);
									funcVal *= std::abs((bRec_R.wi.z / bRec_R.wo.z)*(-path_L[i].wo.z / path_R[j].wi.z));
								}
								pdf_RR[0] = bsdf->pdf(bRec_R);
								BSDFSamplingRecord bRec_R_reverse(bRec_R);
								bRec_R_reverse.reverse();
								pdf_RR[1] = bsdf->pdf(bRec_R_reverse);
							}
							else {
								const PhaseFunction *phase = mediums[id_R].get()->getPhaseFunction();
								if (j == 0) fprintf(stderr, "Badness i: 0\n");
								pRec.wi = path_R[j].wi;
								pRec.wo = -path_L[i].wo;
								if (mode == 1) {
									funcVal = Spectrum(phase->eval(pRec));
								}
								pdf_RR[0] = phase->pdf(pRec);
								PhaseFunctionSamplingRecord pRec_reverse(pRec);
								pRec_reverse.reverse();
								pdf_RR[1] = phase->pdf(pRec_reverse);
							}
							pdf_RR *= path_R[j].pSurvival;

							if (mode == 2 || !funcVal.isZero()) {
								if (mode == 1) f = path_L[i].thru1 * path_R[j].thru0 * funcVal * mRec.transmittance / std::abs(path_L[i].wo.z);
								if (mode == 2) f_pdf = (path_R[j].surf ? mRec.pdfFailure : mRec.pdfSuccess) * pdf_RR[1] * (j ? ratioPdf_R[j - 1] * (path_R[j].epdf[1] / path_R[j].epdf[0]) : Float(1.0));

								w = Float(1.0) + pdf_RR[0] / (path_L[i].vpdf[0] * path_L[i].pSurvival);

								if (i) {
									w += ratio_L[i] * (path_L[i].surf ? mRec.pdfFailure : mRec.pdfSuccessRev) * pdf_RR[0];
								}
								if (j) {
									w += (path_R[j].surf ? mRec.pdfFailure : mRec.pdfSuccess) / path_R[j].epdf[0] *
										(Float(1.0) + pdf_RR[1] / path_R[j - 1].vpdf[0]);
									w += ratio_R[j - 1] * (path_R[j].surf ? mRec.pdfFailure : mRec.pdfSuccess) * pdf_RR[1] * path_R[j].epdf[1] / path_R[j].epdf[0];
								}
								//if (!boost::math::isnormal(w)) {
								//	fprintf(stderr, "WTF (i): %lf\n", w);
								//}
								if (mode == 1) {
									Float etas = 1.0;
									if ((flag_incidentDir && flag_type) || (!flag_incidentDir && !flag_type)) {
										for (int e = 0; e < id_connectMedium + 1; ++e)
											etas *= (m_bsdfs[e]->getEta()*m_bsdfs[e]->getEta());
										etas = Float(1.0) / etas;
									}
									else if ((flag_incidentDir && !flag_type) || (!flag_incidentDir && flag_type)) {
										for (int e = id_connectMedium + 1; e < m_nbLayers; ++e)
											etas *= (m_bsdfs[e]->getEta()*m_bsdfs[e]->getEta());
									}
									else
										;
									Li += f / w * etas;
								}
								if (mode == 2) pdf += f_pdf / w;
							}
						}
					}

					// sample from right
					frame = frames[id_L];

					if (!path_L[i].surf || (path_L[i].wi.z*frame.toLocal(path_L[i].wi).z > 0 && -path_R[j].wo.z*frame.toLocal(-path_R[j].wo).z > 0)) {

						t = (path_R[j].p.z - path_L[i].p.z) / -path_R[j].wo.z;
						if (t > Epsilon) {
							const Medium *medium = mediums[id_connectMedium].get();
							medium->eval(Ray(Point(0.0, 0.0, path_L[i].p.z), -path_R[j].wo, 0.0, t, 0.0), mRec);
							mRec.pdfSuccess /= std::abs(path_R[j].wo.z);
							mRec.pdfSuccessRev /= std::abs(path_R[j].wo.z);

							if (path_L[i].surf) {
								const BSDF *bsdf = m_bsdfs[id_L].get();
								bRec_L.wi = frame.toLocal(path_L[i].wi);
								bRec_L.wo = frame.toLocal(-path_R[j].wo);
								if (mode == 1) funcVal = bsdf->eval(bRec_L);
								pdf_LL[0] = bsdf->pdf(bRec_L);
								BSDFSamplingRecord bRec_L_reverse(bRec_L);
								bRec_L_reverse.reverse();
								pdf_LL[1] = bsdf->pdf(bRec_L_reverse);
							}
							else {
								const PhaseFunction *phase = mediums[id_L].get()->getPhaseFunction();
								if (i == 0) fprintf(stderr, "Badness j: 0\n");
								pRec.wi = path_L[i].wi;
								pRec.wo = -path_R[j].wo;
								if (mode == 1) funcVal = Spectrum(phase->eval(pRec));
								pdf_LL[0] = phase->pdf(pRec);
								PhaseFunctionSamplingRecord pRec_reverse(pRec);
								pRec_reverse.reverse();
								pdf_LL[1] = phase->pdf(pRec_reverse);
							}
							pdf_LL *= path_L[i].pSurvival;

							if (mode == 2 || !funcVal.isZero()) {
								if (mode == 1) f = path_L[i].thru0 * path_R[j].thru1 * funcVal * mRec.transmittance / std::abs(path_R[j].wo.z);
								if (mode == 2) f_pdf = pdf_LL[0] * (path_R[j].surf ? mRec.pdfFailure : mRec.pdfSuccessRev) * ratioPdf_R[j];

								w = Float(1.0) + pdf_LL[0] / path_R[j].vpdf[0];

								if (j) {
									w += ratio_R[j] * (path_R[j].surf ? mRec.pdfFailure : mRec.pdfSuccessRev)*pdf_LL[0];
								}
								if (i) {
									w += (path_L[i].surf ? mRec.pdfFailure : mRec.pdfSuccess) / path_L[i].epdf[0] *
										(Float(1.0) + pdf_LL[1] / path_L[i - 1].vpdf[0]);
									w += ratio_L[i - 1] * (path_L[i].surf ? mRec.pdfFailure : mRec.pdfSuccess)*pdf_LL[1] * path_L[i].epdf[1] / path_L[i].epdf[0];
								}
								//if (!boost::math::isnormal(w)) {
								//	fprintf(stderr, "WTF (j): %lf\n", w);
								//}
								if (mode == 1) {
									Float etas = 1.0;
									if ((flag_incidentDir && flag_type) || (!flag_incidentDir && !flag_type)) {
										for (int e = 0; e < id_connectMedium + 1; ++e)
											etas *= (m_bsdfs[e]->getEta()*m_bsdfs[e]->getEta());
										etas = Float(1.0) / etas;
									}
									else if ((flag_incidentDir && !flag_type) || (!flag_incidentDir && flag_type)) {
										for (int e = id_connectMedium + 1; e < m_nbLayers; ++e)
											etas *= (m_bsdfs[e]->getEta()*m_bsdfs[e]->getEta());
									}
									else
										;
									Li += f / w * etas;
								}
								if (mode == 2) pdf += f_pdf / w;
							}
						}
					}
				}
			}
		}

		if (mode == 1) _val = Li0 + Li * std::abs(_bRec.wo.z);
		if (mode == 2) _pdf = pdf0 + pdf + m_diffusePdf;

		//cout << Li.toString() << endl;
	}

	Float pdfTRT(const BSDFSamplingRecord &_bRec, const std::vector<Frame> &frames, const ref_vector<Medium> &mediums) const {
		Assert(_bRec.sampler);
		Sampler *sampler = _bRec.sampler;

		Float pdf = 0.0;

		BSDFSamplingRecord bRec(_bRec);
		bRec.typeMask = BSDF::ETransmission;
		BSDFSamplingRecord bRecPdf(_bRec);
		//bRecPdf.typeMask = BSDF::ETransmission;

		std::vector<int> wi_id, wo_id;
		wi_id.resize(m_nbLayers);
		wo_id.resize(m_nbLayers);

		std::vector<Vector> wi,wo;
		wi.resize(m_nbLayers);
		wo.resize(m_nbLayers);

		std::vector<Float> ratio;
		ratio.resize(m_nbLayers);

		if (_bRec.wi.z > 0)	for (int i = 0; i < m_nbLayers; ++i) wi_id.push_back(i); 
		if (_bRec.wi.z < 0)	for (int i = 0; i < m_nbLayers; ++i) wi_id.push_back(m_nbLayers - i - 1);
		if (_bRec.wo.z > 0) for (int i = 0; i < m_nbLayers; ++i) wo_id.push_back(i);
		if (_bRec.wo.z < 0) for (int i = 0; i < m_nbLayers; ++i) wo_id.push_back(m_nbLayers - i - 1);
		
		wi[0] = _bRec.wi;
		wo[0] = _bRec.wo;
		ratio[0] = 1.0;

		for (int i = 0; i < m_nbLayers-1; ++i) {
			// wi
			bRec.wi = frames[wi_id[i]].toLocal(wi[i]);
			if (bRec.wi.z * wi[i].z <= 0) return 0.0;
			m_bsdfs[wi_id[i]]->sample(bRec, sampler->next2D());
			wi[i + 1] = -frames[wi_id[i]].toWorld(bRec.wo);
			if (-wi[i+1].z * bRec.wo.z <= 0) return 0.0;

			// wo
			bRec.wi = frames[wo_id[i]].toLocal(wo[i]);
			if (bRec.wi.z * wo[i].z <= 0) return 0.0;
			m_bsdfs[wo_id[i]]->sample(bRec, sampler->next2D());
			wo[i + 1] = -frames[wo_id[i]].toWorld(bRec.wo);
			if (-wo[i + 1].z * bRec.wo.z <= 0) return 0.0;

			// ratio from wo
			Float pdf0 = m_bsdfs[wo_id[i]]->pdf(bRec);
			bRec.reverse();
			Float pdf1 = m_bsdfs[wo_id[i]]->pdf(bRec);
			
			const Medium *medium = wo[i].z > 0 ? mediums[i].get() : mediums[m_nbLayers-2-i].get();
			MediumSamplingRecord mRec;
			Float t = Float(1.0) / std::abs(wo[i+1].z);
			medium->eval(Ray(Point(0, 0, 0), wo[i + 1], 0.0, t, 0.0), mRec);
			mRec.pdfFailure = 1.0;

			if(pdf0 > Epsilon) ratio[i + 1] = ratio[i] * mRec.pdfFailure * pdf1 / pdf0;
		}

		for (int i = 1; i < m_nbLayers; ++i) {
			//cout << "wi" << i << ":" << wi[i].toString() << endl;
			//cout << "wo" << i << ":" << wo[i].toString() << endl;
			if (wi[i].z > 0 && wo[i].z > 0) {
				bRecPdf.wi = frames[wi_id[i]].toLocal(wi[i]);
				bRecPdf.wo = frames[wi_id[i]].toLocal(wo[i]);
				pdf += m_bsdfs[wi_id[i]]->pdf(bRecPdf) * ratio[i];
				//cout << "pdf" << i << ":" << m_bsdfs[wi_id[i]]->pdf(bRecPdf) << endl;
			}
			else
			{ 
				//cout << "WWWWWWWWWWWWWWWWWWW" << endl;
			}
		}

		return pdf + m_diffusePdf;

	}

	void pdfEvaluation(const BSDFSamplingRecord &_bRec, const BSDFSamplingRecord &bRec,
		const std::vector<Frame> &frames, const ref_vector<Medium> &mediums,
		const int mode, Float &samplePdf, Float &evalPdf) const {

		if (m_pdfMode == "const") {
			samplePdf = m_diffusePdf;
			evalPdf = m_diffusePdf;
		}
		else if (m_pdfMode == "TRT") {
			samplePdf = 0.0;
			evalPdf = 0.0;
			for (int i = 0; i < m_pdfRepetitive; ++i) {
				Float _samplePdf = 0.0;
				Float _evalPdf = 0.0;
				if (mode == 1 || mode == 3) _samplePdf = pdfTRT(_bRec, frames, mediums);
				if (mode == 2 || mode == 3)	_evalPdf = pdfTRT(bRec, frames, mediums);
				samplePdf += _samplePdf;
				evalPdf += _evalPdf;
			}
			samplePdf /= m_pdfRepetitive;
			evalPdf /= m_pdfRepetitive;
		}
		else if (m_pdfMode == "bidirStochTRT") {
			ref_vector<Medium> mediumsForPdf;
			setParametersPdf(_bRec, mediumsForPdf);
			samplePdf = 0.0;
			evalPdf = 0.0;
			for (int i = 0; i < m_pdfRepetitive; ++i) {
				BSDFSamplingRecord bRec_tmp(_bRec);
				std::vector<PathInfo> path, path_R_sample, path_R_eval;
				std::vector<Float> ratio, ratio_R_sample, ratio_R_eval;
				std::vector<Float> ratioPdf, ratioPdf_R_sample, ratioPdf_R_eval;
				generatePath(bRec_tmp, frames, mediumsForPdf, m_stochPdfDepth, path, ratio, ratioPdf, false, true);

				Float _samplePdf = 0.0;
				Float _evalPdf = 0.0;
				if (mode == 1 || mode == 3) {
					bRec_tmp.wi = _bRec.wo;
					generatePath(bRec_tmp, frames, mediumsForPdf, m_stochPdfDepth, path_R_sample, ratio_R_sample, ratioPdf_R_sample, true, true);
					Spectrum sampleVal(0.0);
					bidirEvaluation(_bRec, frames, mediumsForPdf, path, ratio, ratioPdf, path_R_sample, ratio_R_sample, ratioPdf_R_sample, 2, sampleVal, _samplePdf);
				}
				if (mode == 2 || mode == 3) {
					bRec_tmp.wi = bRec.wo;
					generatePath(bRec_tmp, frames, mediumsForPdf, m_stochPdfDepth, path_R_eval, ratio_R_eval, ratioPdf_R_eval, true, true);
					Spectrum evalVal(0.0);
					bidirEvaluation(bRec, frames, mediumsForPdf, path, ratio, ratioPdf, path_R_eval, ratio_R_eval, ratioPdf_R_eval, 2, evalVal, _evalPdf);
				}
				samplePdf += _samplePdf;
				evalPdf += _evalPdf;
			}
			samplePdf /= m_pdfRepetitive;
			evalPdf /= m_pdfRepetitive;
		}
		else if (m_pdfMode == "bidirStoch") {
			
			samplePdf = 0.0;
			evalPdf = 0.0;
			for (int i = 0; i < m_pdfRepetitive; ++i) {
				BSDFSamplingRecord bRec_tmp(_bRec);
				std::vector<PathInfo> path, path_R_sample, path_R_eval;
				std::vector<Float> ratio, ratio_R_sample, ratio_R_eval;
				std::vector<Float> ratioPdf, ratioPdf_R_sample, ratioPdf_R_eval;
				generatePath(bRec_tmp, frames, mediums, m_stochPdfDepth, path, ratio, ratioPdf, false, true);

				Float _samplePdf = 0.0;
				Float _evalPdf = 0.0;
				if (mode == 1 || mode == 3) {
					bRec_tmp.wi = _bRec.wo;
					generatePath(bRec_tmp, frames, mediums, m_stochPdfDepth, path_R_sample, ratio_R_sample, ratioPdf_R_sample, true, true);
					
					Spectrum sampleVal(0.0);
					bidirEvaluation(_bRec, frames, mediums, path, ratio, ratioPdf, path_R_sample, ratio_R_sample, ratioPdf_R_sample, 2, sampleVal, _samplePdf);
				}
				if (mode == 2 || mode == 3) {
					bRec_tmp.wi = bRec.wo;
					generatePath(bRec_tmp, frames, mediums, m_stochPdfDepth, path_R_eval, ratio_R_eval, ratioPdf_R_eval, true, true);
					
					Spectrum evalVal(0.0);
					bidirEvaluation(bRec, frames, mediums, path, ratio, ratioPdf, path_R_eval, ratio_R_eval, ratioPdf_R_eval, 2, evalVal, _evalPdf);
				}
				samplePdf += _samplePdf;
				evalPdf += _evalPdf;
			}
			samplePdf /= m_pdfRepetitive;
			evalPdf /= m_pdfRepetitive;
		}
		else {
			cout << "[GY]: Pdf Mode _" << m_pdfMode << "_ is not implement" << endl;
		}

	}

	void evalAndSample(BSDFSamplingRecord &_bRec, Spectrum &evalVal, Float &evalPdf, Spectrum &sampleVal, Float &samplePdf,
		const Point2 &nextSample, EMeasure measure) const {
		Assert(_bRec.sampler);
	
		// Eval(pdf) and Sample(pdf)
		const BSDFSamplingRecord bRec(_bRec);
		BSDFSamplingRecord bRec_tmp(_bRec);

		std::vector<Frame> frames;
		ref_vector<Medium> mediums;
		setParameters(bRec, frames, mediums);

		Float evalPdf_tmp = 0.0; 
		
		if (!(BSDF::getType() & BSDF::ETransmission) && (Frame::cosTheta(bRec.wi) <= 0 || Frame::cosTheta(bRec.wo) <= 0)) {
			evalVal = Spectrum(0.0);
			evalPdf = 0.0;

			if (Frame::cosTheta(bRec.wi) <= 0) {
				sampleVal = Spectrum(0.0);
				samplePdf = 0.0;
			}
			else {
				// sample 
				{
					std::vector<PathInfo> path;
					std::vector<Float> ratio, ratioPdf;
					sampleVal = generatePath(_bRec, frames, mediums, -1, path, ratio, ratioPdf, false, false);
				}
				// sample pdf
				{
					if (sampleVal.isZero()) {
						samplePdf = 0.0;
					}
					else {
						pdfEvaluation(_bRec, bRec, frames, mediums, 1, samplePdf, evalPdf_tmp);
					}
				}
			}
		}
		else {
			{
				// eval, sample, wo, (eval pdf)
				if (m_bidir) {
					std::vector<PathInfo> path;
					std::vector<Float> ratio, ratioPdf;
					sampleVal = generatePath(_bRec, frames, mediums, -1, path, ratio, ratioPdf, false, true);

					// backward sample
					std::vector<PathInfo> path_R;
					std::vector<Float> ratio_R, ratioPdf_R;
					bRec_tmp.wi = bRec.wo;
					generatePath(bRec_tmp, frames, mediums, -1, path_R, ratio_R, ratioPdf_R, true, true);

					bidirEvaluation(bRec, frames, mediums, path, ratio, ratioPdf, path_R, ratio_R, ratioPdf_R, 1, evalVal, evalPdf_tmp);
				}
				else {
					std::vector<PathInfo> path;
					std::vector<Float> ratio, ratioPdf;
					sampleVal = generatePath(_bRec, frames, mediums, -1, path, ratio, ratioPdf, false, false);
					unidirEvaluation(bRec, frames, mediums, path, 1, evalVal, evalPdf_tmp);
				}
			}
			{
				// sample pdf, eval pdf
				pdfEvaluation(_bRec, bRec, frames, mediums, 3, samplePdf, evalPdf);
			}
		}

		if (evalPdf < 0) {
			cout << "[GY]: eval pdf < 0" << endl;
			evalPdf = 0.0;
		}
		if (samplePdf < 0) {
			cout << "[GY]: sample pdf < 0" << endl;
			samplePdf = 0.0;
		}
	}

	Float pdf(const BSDFSamplingRecord &_bRec, EMeasure measure) const {
		Assert(_bRec.sampler);

		const BSDFSamplingRecord bRec_tmp(_bRec);

		if (!(BSDF::getType() & BSDF::ETransmission) && (Frame::cosTheta(_bRec.wi) <= 0 || Frame::cosTheta(_bRec.wo) <= 0)) {
			return 0.0;
		}

		std::vector<Frame> frames;
		ref_vector<Medium> mediums;
		setParameters(_bRec, frames, mediums);

		
		Float pdf_return = 0.0, pdf_tmp = 0.0;
		pdfEvaluation(_bRec, bRec_tmp, frames, mediums, 1, pdf_return, pdf_tmp);
	
		return pdf_return;
	}

	Spectrum sample(BSDFSamplingRecord &_bRec, const Point2 &sample) const {
		Assert(_bRec.sampler);

		if (!(BSDF::getType() & BSDF::ETransmission) && (Frame::cosTheta(_bRec.wi) <= 0)) {
			return Spectrum(0.0);
		}

		std::vector<Frame> frames;
		ref_vector<Medium> mediums;
		setParameters(_bRec, frames, mediums);

		Spectrum sampleVal(0.0);

		{
			std::vector<PathInfo> path;
			std::vector<Float> ratio, ratioPdf;
			sampleVal = generatePath(_bRec, frames, mediums, -1, path, ratio, ratioPdf, false, false);
		}

		return sampleVal;
	}

	Spectrum sample(BSDFSamplingRecord &_bRec, Float &_pdf, const Point2 &sample) const {
		Assert(_bRec.sampler);

		if (!(BSDF::getType() & BSDF::ETransmission) && (Frame::cosTheta(_bRec.wi) <= 0)) {
			return Spectrum(0.0);
		}

		std::vector<Frame> frames;
		ref_vector<Medium> mediums;
		setParameters(_bRec, frames, mediums);
		
		Spectrum sampleVal(0.0);

		{
			std::vector<PathInfo> path;
			std::vector<Float> ratio, ratioPdf;
			sampleVal = generatePath(_bRec, frames, mediums, -1, path, ratio, ratioPdf, false, false);
		}

		{
			Float evalPdf = 0.0;
			pdfEvaluation(_bRec, _bRec, frames, mediums, 1, _pdf, evalPdf);
		}

		return sampleVal;
	}

	Spectrum eval(const BSDFSamplingRecord &_bRec, EMeasure measure) const {
		Assert(_bRec.sampler);

		BSDFSamplingRecord bRec(_bRec);
		BSDFSamplingRecord bRec_tmp(_bRec);

		std::vector<Frame> frames;
		ref_vector<Medium> mediums;
		setParameters(bRec, frames, mediums);

		Spectrum evalVal(0.0);

		if (!(BSDF::getType() & BSDF::ETransmission) && (Frame::cosTheta(bRec.wi) <= 0 || Frame::cosTheta(bRec.wo) <= 0)) {
			
			return evalVal;
		
		}
		else {
			Float evalPdf_tmp = 0;
			if (m_bidir) {
				std::vector<PathInfo> path;
				std::vector<Float> ratio, ratioPdf;
				generatePath(bRec, frames, mediums, -1, path, ratio, ratioPdf, false, true);

				// backward sample
				std::vector<PathInfo> path_R;
				std::vector<Float> ratio_R, ratioPdf_R;
				bRec_tmp.wi = _bRec.wo;
				generatePath(bRec_tmp, frames, mediums, -1, path_R, ratio_R, ratioPdf_R, true, true);

				bidirEvaluation(_bRec, frames, mediums, path, ratio, ratioPdf, path_R, ratio_R, ratioPdf_R, 1, evalVal, evalPdf_tmp);
			}
			else {
				std::vector<PathInfo> path;
				std::vector<Float> ratio, ratioPdf;
				generatePath(bRec, frames, mediums, -1, path, ratio, ratioPdf, false, false);
				unidirEvaluation(_bRec, frames, mediums, path, 1, evalVal, evalPdf_tmp);
			}
		}

		return evalVal;
	}

	Float getEta() const {
		return m_eta;
	}

	Float getRoughness(const Intersection &its, int component) const {
		Log(EError, "Not implemented.");
		return 0.0;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "MultiLayeredBSDF[]" << endl;
		return oss.str();
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		std::string prefix = name.substr(0, name.size() - 2);
		int index = atoi(name.substr(name.size() - 1).c_str());
		if (child->getClass()->derivesFrom(MTS_CLASS(BSDF))) {
			if (prefix == "surface") {
				std::cout << "[GY]: layer " << index + 1 << " of " << m_nbLayers << " :: surface :: bsdf" << std::endl;
				m_bsdfs[index] = static_cast<BSDF *>(child);
			}
			else
				BSDF::addChild(name, child);
		}
		else if (child->getClass()->derivesFrom(MTS_CLASS(Texture2D))) {
			if (prefix == "normal_tex") {
				std::cout << "[GY]: layer " << index + 1 << " of " << m_nbLayers << " :: surface :: normal(texture)" << std::endl;
				m_flag_normals[index] = true;
				m_texture_normals[index] = static_cast<Texture2D *>(child);
			}
			else if (prefix == "sigmaT_tex") {
				std::cout << "[GY]: layer " << index + 1 << " of " << m_nbLayers << " :: medium :: sigmaT(texture)" << std::endl;
				m_flag_sigmaTs[index] = true;
				m_texture_sigmaTs[index] = static_cast<Texture2D *>(child);
			}
			else if (prefix == "density_tex") {
				std::cout << "[GY]: layer " << index + 1 << " of " << m_nbLayers << " :: medium :: density(texture)" << std::endl;
				m_flag_densities[index] = true;
				m_texture_densities[index] = static_cast<Texture2D *>(child);
			}
			else if (prefix == "albedo_tex") {
				std::cout << "[GY]: layer " << index + 1 << " of " << m_nbLayers << " :: medium :: albedo(texture)" << std::endl;
				m_flag_albedos[index] = true;
				m_texture_albedos[index] = static_cast<Texture2D *>(child);
			}
			else if (prefix == "orientation_tex") {
				std::cout << "[GY]: layer " << index + 1 << " of " << m_nbLayers << " :: medium :: orientation(texture)" << std::endl;
				m_flag_orientations[index] = true;
				m_texture_orientations[index] = static_cast<Texture2D *>(child);
			}
			else
				BSDF::addChild(name, child);
		}
		else if (child->getClass()->derivesFrom(MTS_CLASS(PhaseFunction))) {
			if (prefix == "phase") {
				std::cout << "[GY]: layer " << index + 1 << " of " << m_nbLayers << " :: surface :: bsdf" << std::endl;
				m_phaseFunctions[index] = static_cast<PhaseFunction *>(child);
			}
			else
				BSDF::addChild(name, child);
		}
		else
			BSDF::addChild(name, child);
	}

	MTS_DECLARE_CLASS()

private:
	int m_maxDepth;
	Float m_maxSurvivalProb;

	bool m_MISenable;
	bool m_multiLayerSupport;
	bool m_bidirUseAnalog;
	bool m_bidir;
	std::string m_pdfMode;
	int m_stochPdfDepth;
	int m_pdfRepetitive;
	Float m_diffusePdf;
	Float m_eta, m_invEta;

	int m_nbLayers;
		
	ref_vector<BSDF> m_bsdfs;
	std::vector<ref_vector<Medium>> m_mediums, m_pdfMediums;
	ref_vector<PhaseFunction> m_phaseFunctions;
	
	std::vector<Spectrum> m_spectrum_sigmaTs;
	ref_vector<Texture2D> m_texture_sigmaTs;
	std::vector<bool> m_flag_sigmaTs;

	std::vector<Float> m_float_densities;
	ref_vector<Texture2D> m_texture_densities;
	std::vector<bool> m_flag_densities;

	std::vector<Spectrum> m_spectrum_albedos;
	ref_vector<Texture2D> m_texture_albedos;
	std::vector<bool> m_flag_albedos;

	std::vector<Vector> m_vector_orientations;
	ref_vector<Texture2D> m_texture_orientations;
	std::vector<bool> m_flag_orientations;

	std::vector<bool> m_flag_aniso;
	
	std::vector<Vector> m_vector_normals;
	ref_vector<Texture2D> m_texture_normals;
	std::vector<bool> m_flag_normals;
};


MTS_IMPLEMENT_CLASS_S(MultiLayeredBSDF, false, BSDF)
MTS_EXPORT_PLUGIN(MultiLayeredBSDF, "Multi-Layered Texture BRDF with Shading Normal");
MTS_NAMESPACE_END