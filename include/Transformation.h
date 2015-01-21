#ifndef TRANSFORMATION_H
#define TRANSFORMATION_H


class Vector3D;


class Transformation
{
public:

	virtual void execute(const Vector3D &position, Vector3D &transformedPoint, Vector3D &transformedNormal) const = 0;
};


#endif // TRANSFORMATION_H