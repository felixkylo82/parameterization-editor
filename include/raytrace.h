#ifndef RAYTRACE_H
#define RAYTRACE_H

#include <vector>

class Vector3D;

namespace RayTrace
{
	enum IntersectionType
	{
		none,
		point,
		lineSegment,
	};

	////////////////////////////////////////
	// Ray Tracing Functions

	IntersectionType rayTracePoint(const Vector3D &point,
		const Vector3D &start, const Vector3D &direction, Vector3D &hitPoint, double &minDist, bool twoSide = false);

	IntersectionType rayTraceLineSegment(const Vector3D &endPoint1, const Vector3D &endPoint2,
		const Vector3D &start, const Vector3D &direction, Vector3D &hitPoint, double &minDist, bool twoSide = false);

	IntersectionType rayTracePlane(const Vector3D &normal, double offset,
		const Vector3D &start, const Vector3D &direction, Vector3D &hitPoint, double &minDist, bool twoSide = false);

	IntersectionType rayTraceTriangle(const Vector3D &vertex1, const Vector3D &vertex2, const Vector3D &vertex3,
		const Vector3D &start, const Vector3D &direction, Vector3D &hitPoint, double &minDist, bool twoSide = false, bool straightlyInside = false);

	////////////////////////////////////////
	// Plane Intersection Tests

	IntersectionType planeIntersectPoint(const Vector3D &point,
		const Vector3D &planeNormal, double planeOffset);

	IntersectionType planeIntersectLineSegment(const Vector3D &endPoint1, const Vector3D &endPoint2,
		const Vector3D &planeNormal, double planeOffset, Vector3D &hitPoint);

	void planeIntersectTriangle(
		const Vector3D &endPoint1, const Vector3D &endPoint2, const Vector3D &endPoint3,
		const Vector3D &planeNormal, double planeOffset, std::vector<Vector3D> &frontTriangles, std::vector<Vector3D> &backTriangles);

	void planeIntersectTriangle(
		const Vector3D triangle[],
		const Vector3D &planeNormal, double planeOffset, std::vector<Vector3D> &frontTriangles, std::vector<Vector3D> &backTriangles);

	bool threePlanesIntesection(
		const Vector3D &normal1, double offset1,
		const Vector3D &normal2, double offset2,
		const Vector3D &normal3, double offset3,
		Vector3D &hitPoint);

	////////////////////////////////////////
	// Triangle-Triangle Intersection Tests

	// unstable
	bool triangleIntersectTriangle(
		const Vector3D &vertex11, const Vector3D &vertex12, const Vector3D &vertex13,
		const Vector3D &vertex21, const Vector3D &vertex22, const Vector3D &vertex23,
		double threshold = 5e-4, bool verbose = true);

	////////////////////////////////////////
	// Line Segment Intersection Tests

	IntersectionType lineSegmentIntersectPoint(const Vector3D &point,
		const Vector3D &endPoint1, const Vector3D &endPoint2, Vector3D &hitPoint);

	// unstable
	IntersectionType lineSegmentIntersectLineSegment(const Vector3D &endPoint11, const Vector3D &endPoint12,
		const Vector3D &endPoint21, const Vector3D &endPoint22, Vector3D &hitPoint);


	////////////////////////////////////////
	// Others

	//unstable
	bool insideTriangle(
		const Vector3D &vertex1, const Vector3D &vertex2, const Vector3D &vertex3,
		const Vector3D &testPoint);

};

#endif
