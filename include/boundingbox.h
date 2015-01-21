#ifndef BOUNDING_BOX_H
#define BOUNDING_BOX_H


#include <map>
#include <list>

#include "vector3d.h"


class Box;
class BoundingBox;
class LeafBoundingBox;
class InternalBoundingBox;


class Box
{
public:

	static bool rayTraceTwoSide;
	static void(*rayTraceFunction) (void *datum, const Vector3D &start, const Vector3D &direction, void **nearestDatum, Vector3D &hitPoint, double &minDist, bool twoSide);

	Box(const Vector3D &min, const Vector3D &max);

	const Vector3D& getMin() const;
	const Vector3D& getMax() const;

	bool intersectBox(const Box &b) const;
	bool intersectRay(const Vector3D &start, const Vector3D &direction) const;

	bool operator<(const Box &box) const;

protected:

	const Vector3D min;
	const Vector3D max;

};

class BoundingBox : public Box
{
public:

	enum Type
	{
		leafBox, internalBox
	};

	BoundingBox(const Vector3D &min, const Vector3D &max);
	virtual ~BoundingBox();

	virtual Type getType() const = 0;

	virtual void clear() = 0;

	virtual bool putDatum(const Box &box, void *datum) = 0;
	//virtual void *getDatum(const Vector3D &position) = 0;

	virtual void getData(const Box &box, std::list<void *> &data) const = 0;

	virtual void rayTrace(const Vector3D &start, const Vector3D &direction, void **nearestDatum, Vector3D &hitPoint, double &minDist) = 0;

	bool empty() const;

protected:

	bool emptyFlag;
};

class LeafBoundingBox : public BoundingBox
{
public:

	static unsigned int maxDataSize;

	LeafBoundingBox(const Vector3D &min, const Vector3D &max);
	virtual ~LeafBoundingBox();

	Type getType() const
	{
		return leafBox;
	}

	void clear();

	bool  putDatum(const Box &box, void *datum);
	//void *getDatum(const Vector3D &position);

	void getData(const Box &box, std::list<void *> &data) const;

	void rayTrace(const Vector3D &start, const Vector3D &direction, void **nearestDatum, Vector3D &hitPoint, double &minDist);

	bool isFull() const;

	const std::multimap<Box, void*> &getData() const;

	bool getSeparable() const;
	void setSeparable(bool separable);

private:

	std::multimap<Box, void*> data;

	bool separable;
};


class InternalBoundingBox : public BoundingBox
{
public:

	//InternalBoundingBox(const Vector3D &min, const Vector3D &max);
	InternalBoundingBox(const Vector3D &min, const Vector3D &max, int depth);
	virtual ~InternalBoundingBox();

	Type getType() const
	{
		return internalBox;
	}

	void clear();

	bool  putDatum(const Box &position, void *datum);
	bool  putDatum(const Box &position, void *datum, bool enableSplit);
	//void *getDatum(const Vector3D &position);

	void getData(const Box &box, std::list<void *> &data) const;

	void rayTrace(const Vector3D &start, const Vector3D &direction, void **nearestDatum, Vector3D &hitPoint, double &minDist);

protected:

	bool isChildrenBorn;

	BoundingBox *children[8];

	bool produceChildren(int depth);
	bool produceLeafChildren();
};


#endif // BOUNDING_BOX_H