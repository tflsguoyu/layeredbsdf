#include <mitsuba/render/texture.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/fresolver.h>

MTS_NAMESPACE_BEGIN

class InstancedTexture : public Texture2D {
public:
    InstancedTexture(const Properties &props) : Texture2D(props) {
        m_divideReso.x = props.getInteger("divideX");
        m_divideReso.y = props.getInteger("divideY");
        m_reso.x = props.getInteger("tileX");
        m_reso.y = props.getInteger("tileY");
        if ( props.hasProperty("fillBlockInfo") )
        {
            if ( props.hasProperty("blockInfo") || props.hasProperty("blockFile") )
                Log(EError, "Block information declared more than once!");
            m_blockID.resize(m_reso.x*m_reso.y, props.getInteger("fillBlockInfo"));

            m_blockFile = "";
        }
        else if ( props.hasProperty("blockInfo") )
        {
            if ( props.hasProperty("blockFile") )
                Log(EError, "Block information declared more than once!");

            std::istringstream iss(props.getString("blockInfo"));
            m_blockID.resize(m_reso.x*m_reso.y);
            for ( int i = 0; i < m_reso.x*m_reso.y; ++i )
            {
                if ( !(iss >> m_blockID[i]) )
                    Log(EError, "Failed to parse the information for block %d", i);
                if ( m_blockID[i] < 0 || m_blockID[i] >= m_divideReso.x*m_divideReso.y )
                    Log(EError, "Invalid block id: %d", m_blockID[i]);
            }

            m_blockFile = "";
        }
        else
            m_blockFile = props.getString("blockFile");

        float a00, a01, a10, a11;
        sscanf(props.getString("uvTransform", "1.0 0.0 0.0 1.0").c_str(), "%f %f %f %f", &a00, &a01, &a10, &a11);
        m_uvTransform = Matrix2x2(a00, a01, a10, a11);
    }

    InstancedTexture(Stream *stream, InstanceManager *manager) : Texture2D(stream, manager) {
        m_divideReso = Vector2i(stream);
        m_reso = Vector2i(stream);
        m_blockFile = stream->readString();
        m_uvTransform = Matrix2x2(stream);
    }
	
    void serialize(Stream *stream, InstanceManager *manager) const {
        Texture2D::serialize(stream, manager);
		
        m_divideReso.serialize(stream);
        m_reso.serialize(stream);
        stream->writeString(m_blockFile);
        m_uvTransform.serialize(stream);
    }
	
	void configure() {
        if (m_nested == NULL)
            Log(EError, "The scale plugin needs a nested texture!");
		
            if ( m_blockFile != "" ) {
				fs::path resolved = Thread::getThread()->getFileResolver()->resolve(m_blockFile);
				
				ref<FileStream> fs = new FileStream(resolved, FileStream::EReadOnly);
				int sz = fs->readInt();
				if ( sz != m_reso.x*m_reso.y )
					Log(EError, "Block information size mismatch: expected %d but got %d", m_reso.x*m_reso.y, sz);

				m_blockID.resize(sz);
				fs->readIntArray(&m_blockID[0], sz);
			}

            for ( int i = 0; i < m_reso.x*m_reso.y; ++i )
                if ( m_blockID[i] < 0 || m_blockID[i] >= m_divideReso.x*m_divideReso.y )
                    Log(EError, "Invalid block id: %d", m_blockID[i]);
    }
	
    void addChild(const std::string &name, ConfigurableObject *child) {
        if (child->getClass()->derivesFrom(MTS_CLASS(Texture)))
            m_nested = static_cast<Texture2D *>(child);
        else
            Texture::addChild(name, child);
    }

    inline Spectrum eval(const Point2 &_uv) const {
        Vector2 uv = m_uvTransform*Vector2(_uv);
		Point2 p(uv[0] - std::floor(uv[0]), uv[1] - std::floor(uv[1]));
		
		p.x *= m_reso.x;
		p.y *= m_reso.y;
		
		int bx = math::floorToInt(p.x), by = math::floorToInt(p.y);
		bx = math::clamp(bx, 0, m_reso.x - 1);
		by = math::clamp(by, 0, m_reso.y - 1);
		
        int bid = m_blockID[by*m_reso.x + bx];
        int tx = bid % m_divideReso.x, ty = bid/m_divideReso.x;

        p.x = (static_cast<Float>(tx) + p.x - std::floor(p.x))/static_cast<Float>(m_divideReso.x);
        p.y = (static_cast<Float>(ty) + p.y - std::floor(p.y))/static_cast<Float>(m_divideReso.y);
		return m_nested->eval(p);
    }

    Spectrum eval(const Point2 &uv, const Vector2 &d0, const Vector2 &d1) const {
        /* Filtering is currently not supported */
        return InstancedTexture::eval(uv);
    }

    bool usesRayDifferentials() const {
        return false;
    }
	
    Spectrum getAverage() const {
        return m_nested->getAverage();
    }

    Spectrum getMaximum() const {
        return m_nested->getMaximum();
    }

    Spectrum getMinimum() const {
        return m_nested->getMinimum();
    }	
	
    bool isConstant() const {
        return false;
    }

    bool isMonochromatic() const {
        return m_nested->isMonochromatic();
    }

    std::string toString() const {
        return "InstancedTexture[]";
    }

    MTS_DECLARE_CLASS()

protected:
	ref<Texture2D> m_nested;

    Vector2i m_divideReso;
    Vector2i m_reso;

    Matrix2x2 m_uvTransform;
	
    std::string m_blockFile;
    std::vector<int> m_blockID;
};

MTS_IMPLEMENT_CLASS_S(InstancedTexture, false, Texture2D)
MTS_EXPORT_PLUGIN(InstancedTexture, "Instanced texture");
MTS_NAMESPACE_END
