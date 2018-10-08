/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/scene.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/core/statistics.h>
#include <boost/algorithm/string.hpp>
#include <algorithm>

MTS_NAMESPACE_BEGIN


class HomogeneousAnisotropicMedium : public Medium {
public:
    HomogeneousAnisotropicMedium(const Properties &props) : Medium(props) {
        m_density = props.getFloat("density");
        m_albedo = props.getSpectrum("albedo");
        m_orientation = props.getVector("orientation");
    }

    /* Unserialize from a binary data stream */
    HomogeneousAnisotropicMedium(Stream *stream, InstanceManager *manager) : Medium(stream, manager) {
        m_density = stream->readFloat();
        m_albedo = Spectrum(stream);
        m_orientation = Vector(stream);
        configure();
    }

    /* Serialize the volume to a binary data stream */
    void serialize(Stream *stream, InstanceManager *manager) const {
        Medium::serialize(stream, manager);
        stream->writeFloat(m_density);
        m_albedo.serialize(stream);
        m_orientation.serialize(stream);
    }

    void configure() {
        Medium::configure();
        if ( !m_phaseFunction->needsDirectionallyVaryingCoefficients() )
            Log(EError, "Phase function needs to be anisotropic!");
        if ( m_orientation.isZero() )
            Log(EError, "Invalid orientation!");
        else
            m_orientation = normalize(m_orientation);
    }

    Spectrum evalTransmittance(const Ray &ray, Sampler *sampler) const {
        Float t = std::max(ray.maxt - ray.mint, static_cast<Float>(0.0));
        return Spectrum(math::fastexp(-t*lookupDensity(ray.d)));
    }

    bool sampleDistance(const Ray &ray, MediumSamplingRecord &mRec, Sampler *sampler) const {
        Float sigmaT = lookupDensity(ray.d);
        Float sampledDistance = -math::fastlog(sampler->next1D())/sigmaT;
        bool success = true;
        if ( sampledDistance < ray.maxt - ray.mint ) {
            mRec.t = sampledDistance + ray.mint;
            mRec.p = ray(mRec.t);
            mRec.sigmaS = m_albedo*sigmaT;
            mRec.sigmaA = Spectrum(sigmaT) - mRec.sigmaS;
            mRec.orientation = m_orientation;

            /* Fail if there is no forward progress
               (e.g. due to roundoff errors) */
            if ( mRec.p == ray.o ) success = false;
        }
        else {
            sampledDistance = ray.maxt - ray.mint;
            success = false;
        }

        Float expVal = math::fastexp(-sigmaT*sampledDistance);
        mRec.pdfFailure = expVal;
        mRec.transmittance = Spectrum(expVal);
        mRec.pdfSuccess = mRec.pdfSuccessRev = sigmaT*expVal;
        mRec.time = ray.time;
        mRec.medium = this;

        if ( mRec.transmittance.max() < 1e-20 ) mRec.transmittance = Spectrum(0.0f);
        return success;
    }

    void eval(const Ray &ray, MediumSamplingRecord &mRec) const {
        Float sigmaT = lookupDensity(ray.d);
        Float distance = ray.maxt - ray.mint;

        mRec.sigmaS = m_albedo*sigmaT;
        mRec.sigmaA = Spectrum(sigmaT) - mRec.sigmaS;
        mRec.orientation = m_orientation;

        Float expVal = math::fastexp(-sigmaT*distance);
        mRec.pdfFailure = expVal;
        mRec.transmittance = Spectrum(expVal);
        mRec.pdfSuccess = mRec.pdfSuccessRev = sigmaT*expVal;
        mRec.time = ray.time;
        mRec.medium = this;

        if ( mRec.transmittance.max() < 1e-20 ) mRec.transmittance = Spectrum(0.0f);
    }

    bool isHomogeneous() const {
        return true;
    }

	void setMediumProp(const Float &density, const Spectrum &albedo, const Vector &orientation) {
		m_density = density;
		m_albedo = albedo;
		m_orientation = normalize(orientation);
	}


    std::string toString() const {
        std::ostringstream oss;
        oss << "HomogeneousAnisotropicMedium[" << endl
            << "  density = " << m_density << "," << endl
            << "  albedo = " << m_albedo.toString() << "," << endl
            << "  orientation = " << m_orientation.toString() << "," << endl
            << "  phase = " << indent(m_phaseFunction->toString()) << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()


protected:
    inline Float lookupDensity(const Vector &d) const {
        return m_density*m_phaseFunction->sigmaDir(dot(d, m_orientation));
    }

    Float m_density;
    Spectrum m_albedo;
    Vector m_orientation;
};

MTS_IMPLEMENT_CLASS_S(HomogeneousAnisotropicMedium, false, Medium)
MTS_EXPORT_PLUGIN(HomogeneousAnisotropicMedium, "Homogeneous anisotropic medium");
MTS_NAMESPACE_END
