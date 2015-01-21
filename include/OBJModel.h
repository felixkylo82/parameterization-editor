#ifndef OBJ_MODEL_H
#define OBJ_MODEL_H


#include <vector>


class Transformation;
class Vector3D;


struct PolygonElement
{
	std::vector<size_t> vertexIndices;
	std::vector<size_t> normalIndices;
	std::vector<size_t> textureCoordinateIndices;
};

class OBJModel
{
public:
	OBJModel();
	OBJModel(const OBJModel &model);
	OBJModel(std::ifstream &is);
	OBJModel(const char filename[]);

	virtual ~OBJModel();

	bool loadFile(const char filename[]);
	bool loadFile(std::ifstream &is);
	void saveFile(const char filename[]) const;
	void saveFile(std::ofstream &os) const;

	size_t getNumberOfVertices() const;
	size_t getNumberOfNormal() const;
	size_t getNumberOfTextureCoordinates() const;
	size_t getNumberOfPolygons() const;

	const PolygonElement *getPolygon(size_t i) const;
	const Vector3D       &getVertex(size_t i) const;
	const Vector3D       &getNormal(size_t i) const;
	const Vector3D       &getTextureCoordinate(size_t i) const;

	size_t addVertex(const Vector3D &position);
	size_t addNormal(const Vector3D &normal);
	size_t addTextureCoordinate(const Vector3D &textureCoordinate);

	size_t addPolygon(const PolygonElement &polygon);

	void rotate(const double rotationMatrix[]);
	void translate(const Vector3D &displacement);
	void scale(const Vector3D &scalingFactor);

	void tranform(const Transformation& transformation);

	double normalize();
	void   reverseOrientation();
	void   collectCloseVertices(float epsilon = 1e-5f);

	void computeCentroid(double &area, Vector3D &centroid) const;

	void render(float epsilon = 0.0f) const;

	void clear();

	void genQuad(const Vector3D v[4], unsigned int resolution);
	void genQuad(const Vector3D v[4], const Vector3D n[4], unsigned int resolution);
	void genQuad(const Vector3D v[4], unsigned int resY, unsigned int resX);
	void genQuad(const Vector3D v[4], const Vector3D n[4], unsigned int resX, unsigned int resY);

	void genTriangle(const Vector3D v[3], unsigned int resolution);
	void genTriangle(const Vector3D v[3], const Vector3D n[3], unsigned int resolution);

	void genBox(const unsigned int resolution[3]);
	void genSphere(unsigned int resolution);
	void genRhombicDodechedron(const unsigned int resolution[4]);
	void genIcosahedron(unsigned int resolution);
	void genBrick();

protected:
	std::vector<Vector3D> vertexArray;
	std::vector<Vector3D> normalArray;
	std::vector<Vector3D> textureCoordinateArray;

	std::vector <PolygonElement *> polygons;
};

#endif // OBJ_MODEL_H
