#ifndef PATCH_H
#define PATCH_H

class Vector3D;
class OBJModel;
class DoublyConnectedEdgeList;

class Patch
{
public:

	Patch(unsigned int numberOfControlPoints, unsigned int resolution = 4);
	virtual ~Patch();

	const Vector3D* getControlPoints() const;
	const Vector3D* getControlNormal() const;

	unsigned int getNumberOfControlPoints() const;
	unsigned int getResolution() const;

	const OBJModel &getMesh() const;

	void resampling(const DoublyConnectedEdgeList &dcel);

	void updateCentroid();
	void getCentroid(double &area, Vector3D &centroid);

	void render() const;

	virtual void setControlPoints(const Vector3D controlPoints[], const Vector3D controlNormal[]) = 0;
	virtual void setResolution(int resolution) = 0;

protected:

	Vector3D*     controlPoints;
	Vector3D*     controlNormal;
	unsigned int  numberOfControlPoints;

	OBJModel*     mesh;
	unsigned int  resolution;

	Vector3D*     centroid;
	double        area;

	friend class PatchList;
};

class QuadPatch : public Patch
{
public:

	QuadPatch(const Vector3D controlPoints[], const Vector3D controlNormal[], unsigned int resolution = 4);
	virtual ~QuadPatch();

	void setControlPoints(const Vector3D controlPoints[], const Vector3D controlNormal[]);
	void setResolution(int resolution);
};

class TriPatch : public Patch
{
public:

	TriPatch(const Vector3D controlPoints[], const Vector3D controlNormal[], unsigned int resolution = 4);
	virtual ~TriPatch();

	void setControlPoints(const Vector3D controlPoints[], const Vector3D controlNormal[]);
	void setResolution(int resolution);
};


#endif //  PATCH_H
