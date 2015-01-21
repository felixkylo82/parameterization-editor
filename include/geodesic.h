#ifndef GEODESIC_H
#define GEODESIC_H


#include "dcel.h"


class Geodesic : public DoublyConnectedEdgeList::SearchPathList
{
public:

	static Geodesic *Geodesic::compute(
		const DoublyConnectedEdgeList &dcel, const Vector3D &from, const Vector3D &to, float epsilon = 1e-5f);

	static Vector3D findRelativePositionOnSurface(
		const DoublyConnectedEdgeList &dcel,
		const Vector3D &v1, const Vector3D &v2, const Vector3D &v3,
		const Vector3D &relativePosition);

	static Vector3D findRelativePositionOnSurface(
		const DoublyConnectedEdgeList &dcel,
		const Vector3D &v1, const Vector3D &v2, const Vector3D &v3,
		double r, double theta);

	static double length(const Geodesic::const_iterator &start, const Geodesic::const_iterator &end);

	//static bool   isOrdered(Geodesic *g1, Geodesic *g2, Geodesic *g3);

	static double timeForSearchingNearestFace;
	static double timeForGeodesicOptimization;

	bool validate(const DoublyConnectedEdgeList &dcel) const;

	double length() const;

	Vector3D findMidPoint(double distance) const;

	Vector3D forwardDirection() const;
	Vector3D backwardDirection() const;

	void optimize(const DoublyConnectedEdgeList &dcel);

	void print();
	void render();

protected:

	bool addPoint(const DoublyConnectedEdgeList &dcel, Geodesic::iterator &iter1);
	void movePoint(const DoublyConnectedEdgeList &dcel, Geodesic::iterator &iter1);
	bool mergePoint2(const DoublyConnectedEdgeList &dcel, Geodesic::iterator &iter1);
	void mergePoint3(const DoublyConnectedEdgeList &dcel, Geodesic::iterator &iter1);

#if 0
	bool mergePoint(const DoublyConnectedEdgeList &dcel, Geodesic::iterator &iter1, Geodesic::iterator &iter2, Geodesic::iterator &iter3);
#endif

private:

};


#endif // GEODESIC_H