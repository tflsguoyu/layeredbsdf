#include <algorithm>
#include <mitsuba/core/frame.h>

MTS_NAMESPACE_BEGIN


struct SphericalGaussian
{
    inline SphericalGaussian() {}

    inline SphericalGaussian(Float weight, Float kappa, const Vector &mu) {
        init(weight, kappa, mu);
    }

    inline SphericalGaussian(const SphericalGaussian &sg0, const SphericalGaussian &sg1) {
        init(sg0, sg1);
    }

    inline void init(Float weight, Float kappa, const Vector &mu) {
        this->weight = weight;
        this->kappa = kappa;
        this->mu = mu;
        this->norm = kappa/(1.0f - math::fastexp(-2.0f*kappa));
    }

    inline void init(const SphericalGaussian &sg0, const SphericalGaussian &sg1) {
        mu = sg0.kappa*sg0.mu + sg1.kappa*sg1.mu;
        kappa = mu.length();
        mu /= kappa;
        weight = sg0.weight*sg1.weight*math::fastexp(-sg0.kappa - sg1.kappa + kappa);
        norm = kappa/(1.0f - math::fastexp(-2.0f*kappa));
    }

    inline Float eval(const Vector &x) const {
        return weight*math::fastexp(kappa*(dot(mu, x) - 1.0f));
    }

    inline Vector sample(const Point2 &sample, Float &pdf) const {
        Float z = math::fastlog(kappa*sample.y/norm + math::fastexp(-2.0f*kappa))/kappa + 1.0f;
        Float r = math::safe_sqrt(1.0f - z*z);
        Float sinPhi, cosPhi;
        math::sincos(2.0f*M_PI*sample.x, &sinPhi, &cosPhi);
        pdf = INV_TWOPI*norm*math::fastexp(kappa*(z - 1.0f));
        return Frame(mu).toWorld(Vector(r*cosPhi, r*sinPhi, z));
    }

    inline Float pdf(const Vector &x) const {
        return INV_TWOPI*norm*math::fastexp(kappa*(dot(mu, x) - 1.0f));
    }

    Float weight, kappa, norm;
    Vector mu;
};


template <size_t nlobes>
class SphericalGaussianMixture
{
    template <size_t nlobes1> friend class SphericalGaussianMixture;

public:
    SphericalGaussianMixture(const Float *params = NULL) {
        if ( params ) {
            cdf[0] = 0.0f;
            for ( size_t i = 0; i < nlobes; ++i ) {
                components[i].init(params[i + i], params[i + i + 1], Vector(0.0f, 0.0f, 1.0f));
                cdf[i + 1] = cdf[i] + components[i].weight/components[i].norm;
            }
            for ( size_t i = 1; i <= nlobes; ++i )
                cdf[i] /= cdf[nlobes];
        }
    }

    template <size_t n0, size_t n1>
    SphericalGaussianMixture(const SphericalGaussianMixture<n0> &sgm0, const SphericalGaussianMixture<n1> &sgm1) {
        init(sgm0, sgm1);
    }

    template <size_t n0, size_t n1>
    void init(const SphericalGaussianMixture<n0> &sgm0, const SphericalGaussianMixture<n1> &sgm1) {
        BOOST_STATIC_ASSERT( nlobes == n0*n1 );

        for ( size_t i = 0, k = 0; i < n0; ++i )
            for ( size_t j = 0; j < n1; ++j, ++k )
                components[k].init(sgm0.components[i], sgm1.components[j]);
        cdf[0] = 0.0f;
        for ( size_t i = 0; i < nlobes; ++i )
            cdf[i + 1] = cdf[i] + components[i].weight/components[i].norm;
        for ( size_t i = 1; i <= nlobes; ++i )
            cdf[i] /= cdf[nlobes];
    }


    Float eval(const Vector &x) const {
        Float ret = 0.0;
        for ( size_t i = 0; i < nlobes; ++i )
            ret += components[i].eval(x);
        return ret;
    }

    Vector sample(const Point2 &_sample, Float &_pdf) const {
        Point2 sample(_sample);

        size_t idx = std::lower_bound(cdf + 1, cdf + nlobes + 1, sample.y) - cdf - 1;
        Float p0 = cdf[idx + 1] - cdf[idx];
        sample.y = (cdf[idx + 1] - sample.y)/p0;
        Vector ret = components[idx].sample(sample, _pdf);
        _pdf *= p0;
        for ( size_t i = 0; i < nlobes; ++i )
            if ( i != idx )
                _pdf += (cdf[i + 1] - cdf[i])*components[i].pdf(ret);
        return ret;
    }

    Float pdf(const Vector &x) const {
        Float ret = 0.0;
        for ( size_t i = 0; i < nlobes; ++i )
            ret += (cdf[i + 1] - cdf[i])*components[i].pdf(x);
        return ret;
    }

    void setMu(const Vector &mu) {
        for ( size_t i = 0; i < nlobes; ++i )
            components[i].mu = mu;
    }

protected:
    SphericalGaussian components[nlobes];
    Float cdf[nlobes + 1];
};


MTS_NAMESPACE_END
