/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

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

#include <mitsuba/render/phase.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/core/frame.h>
#include <boost/math/special_functions.hpp>
#include <mitsuba/core/warp.h>
#include <mitsuba/core/vmf.h>

MTS_NAMESPACE_BEGIN

/**
 * von Mises-Fisher Phase Function
 */
class vMFPhaseFunction : public PhaseFunction {
public:
	vMFPhaseFunction(const Properties &props) 
		: PhaseFunction(props) {
		m_kappa = props.getFloat("kappa");
	}

	vMFPhaseFunction(Stream *stream, InstanceManager *manager) 
		: PhaseFunction(stream, manager) {
		m_kappa = stream->readFloat();
        configure();
	}

	virtual ~vMFPhaseFunction() { }


	void serialize(Stream *stream, InstanceManager *manager) const {
		PhaseFunction::serialize(stream, manager);
		stream->writeFloat(m_kappa);
	}


	void configure() {
		PhaseFunction::configure();
		m_type = EAngleDependence;
        m_distr = VonMisesFisherDistr(m_kappa);
	}

	Float sample(PhaseFunctionSamplingRecord &pRec, 
			Sampler *sampler) const {

		Vector dir = m_distr.sample(sampler->next2D());
		pRec.wo = Frame(-pRec.wi).toWorld(dir);

		return 1.0f;
	}

	Float sample(PhaseFunctionSamplingRecord &pRec, 
			Float &pdf, Sampler *sampler) const {
		vMFPhaseFunction::sample(pRec, sampler);
		pdf = vMFPhaseFunction::eval(pRec);
		return 1.0f;
	}

	Float eval(const PhaseFunctionSamplingRecord &pRec) const {
		Float val = dot(-pRec.wi, pRec.wo);
        return m_distr.eval(val);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "vMFPhaseFunction[kappa=" << m_kappa << "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()

private:
	Float m_kappa;
    VonMisesFisherDistr m_distr;
};


MTS_IMPLEMENT_CLASS_S(vMFPhaseFunction, false, PhaseFunction)
MTS_EXPORT_PLUGIN(vMFPhaseFunction, "von Mises-Fisher phase function");
MTS_NAMESPACE_END
