#ifndef RESAMPLING_H
#define RESAMPLING_H

#include "Transformation.h"

class Vector3D;
class OBJModel;
class DoublyConnectedEdgeList;

class Resampling : public Transformation
{
public:

	static double timeUsed;

	Resampling();
	Resampling(const DoublyConnectedEdgeList &dcel);

	void setSurface(const DoublyConnectedEdgeList &dcel);
	void transform(OBJModel &mesh);

	virtual void execute(const Vector3D &position, Vector3D &nearestPoint, Vector3D &normal) const;

	void findNearestFace(const Vector3D &position, Vector3D &nearestPoint, Vector3D &normal) const;
	void rayTrace(const Vector3D &position, Vector3D &nearestPoint, Vector3D &normal) const;

private:

	const DoublyConnectedEdgeList *dcel;
};

#endif // RESAMPLING_H
