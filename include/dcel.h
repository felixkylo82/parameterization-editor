#ifndef DCEL_H
#define DCEL_H

#define _USE_DCEL_SEARCH_TREE

#ifdef _USE_DCEL_SEARCH_TREE
#include <cstdio>
#include "vptree.h"
#endif
#include "vector3D.h"

#include <list>
#include <set>
#include <map>
#include <cfloat>
#include <vector>
#include <fstream>

class OBJModel;
class InternalBoundingBox;

/***************************************************************************************************
*
* class DoublyConnectedEdgeList
*
* usage: to maintain and manipulate a doubly connected edge list (DCEL)
*        which is particular useful for efficient algorithms but not for rendering
*
***************************************************************************************************/

class DoublyConnectedEdgeList
{
public:

	class Primitive;

	class Vertex;
	class Edge;
	class Face;

	struct SearchPathNode;

	class CostComparison
	{
	public:

		bool operator()(const Vertex* x, const Vertex* y) const;
	};

	class VertexComparison
	{
	public:

		bool operator()(const Vertex* x, const Vertex* y) const;
	};

	class VertexReset
	{
	public:

		void operator()(Vertex* currentVertex) const;
	};

	class EdgeComparison
	{
	public:

		bool operator()(const Edge* x, const Edge* y) const;
	};

	typedef std::set<Vertex*, VertexComparison> VertexSet;
	typedef std::set<Edge  *, EdgeComparison  > EdgeSet;
	typedef std::set<Face *>                    FaceSet;

	typedef std::set<const Face *> ConstFaceSet;

	typedef std::list<SearchPathNode> SearchPathList;

	class Primitive
	{
	public:

		enum Type
		{
			vertex, edge, face,
		};

		virtual DoublyConnectedEdgeList::Primitive::Type getPrimitiveType() const = 0;

		virtual Vector3D computeNormal() const = 0;

		virtual void render(float epsilon = 0.0) const = 0;

	protected:

		size_t  id; // update only when outputting the VRML

		friend class DoublyConnectedEdgeList;
	};

	class Face : public Primitive
	{
	public:

		DoublyConnectedEdgeList::Primitive::Type getPrimitiveType() const
		{
			return DoublyConnectedEdgeList::Primitive::face;
		}

		const Edge *getComponent() const;	// for computing geodesic

		void getNeighboringFaces(ConstFaceSet &faces) const;
		void getVertices(std::vector<const Vertex *> &v) const;

		Vector3D computeNormal() const;

		void render(float epsilon = 0.0) const;

	private:

		Face(Edge *component);  // use DoublyConnectedEdgeList::addFace to create new instances

		void flipFaceOrientation();

		Edge *component;

		bool visited;

		friend class DoublyConnectedEdgeList;
		friend class FaceComparison;
	};

	class Edge : public Primitive
	{
	public:

		DoublyConnectedEdgeList::Primitive::Type getPrimitiveType() const
		{
			return DoublyConnectedEdgeList::Primitive::edge;
		}

		double getCost() const;
		//void  setCost(float cost);

		Vertex *getOrigin() const;
		//void    setVertex(Vertex *vertex);
		//Edge* setPreviousEdge(Edge* previous);
		//Edge* setTwinEdge    (Edge* twin    );
		//Face* setIncidentFace(Face* incidentFace);

		Edge* getNextEdge() const;
		Edge* getPreviousEdge() const;
		Edge* getTwinEdge() const;
		Face* getIncidentFace() const;

		Edge* setNextEdge(Edge* next);

		Vector3D computeNormal() const;

		void render(float epsilon = 0.0) const;

	private:

		Edge(Vertex *origin, double cost = 0.0f);

		Vertex *origin;

		double cost;

		Edge* next;
		Edge* previous;
		Edge* twin;

		Face* incidentFace;

		friend class DoublyConnectedEdgeList;
		friend class Face;
		friend class Vertex;
		friend class EdgeComparison;
	};

	class Vertex : public Primitive
	{
	public:

		virtual ~Vertex();

		DoublyConnectedEdgeList::Primitive::Type getPrimitiveType() const
		{
			return DoublyConnectedEdgeList::Primitive::vertex;
		}

		const Vector3D &getPosition() const;

		double getPathCost() const;

		const EdgeSet &getOutgoingEdges() const;

		void addNeighbor(Edge* edge);	// for computing geodesic

		size_t getNumberOfNeighbors() const;

		Vector3D computeNormal() const;
		Vector3D computeNormal(double edgeBias) const;

		void render(float epsilon = 0.0) const;

	protected:

		Vertex(const Vector3D &v);

	private:

		Vector3D position;
		Vector3D texCoord;

		EdgeSet outgoingEdges;

		// status
		bool    visited;
		Vertex *backtrace;
		double  pathCost;

		friend class DoublyConnectedEdgeList;
		friend class Face;
		friend class CostCopmparison;
		friend class VertexComparison;
		friend class VertexReset;
	};

	struct SearchPathNode
	{
		SearchPathNode(const DoublyConnectedEdgeList::Primitive *at);
		SearchPathNode(const DoublyConnectedEdgeList::Primitive *at, const Vector3D &position);

		const DoublyConnectedEdgeList::Primitive *at;
		Vector3D position;
	};


	static const unsigned int max_sides;


	DoublyConnectedEdgeList();
	DoublyConnectedEdgeList(const OBJModel &model);

	virtual ~DoublyConnectedEdgeList();


	void convertFromOBJModel(const OBJModel &model);

	void clear();

	size_t getNumberOfVertices() const;
	size_t getNumberOfEdges() const;
	size_t getNumberOfFaces() const;

	void getPolygonSide(unsigned int &min, unsigned int &max) const;

	bool flipFaceOrientation(Face *face);

	//void toCloseManifold();

	void fillMissingPolygons();

	bool validate() const;
	void render(float epsilon = 0.0f) const;
	void triangulate();

	void convertToOBJModel(OBJModel &model, bool enableEdgeBias, bool outputTexCoord) const;

	bool loadFile(const char filename[]);
	void saveFile(const char filename[]) const;

	bool loadFile(std::ifstream &fin);
	void saveFile(std::ofstream &fout) const;

	const VertexSet &getVertexSet() const;
	const EdgeSet   &getEdgeSet() const;
	const FaceSet   &getFaceSet() const;

	Vertex *addVertex(const Vector3D &position);
	void    addVertex(Vertex *vertex);
	Face   *removeVertex(const Vector3D &position);
	Face   *removeVertex(Vertex *vertex);
	void    moveVertex(Vertex *v1, const Vector3D position);
	Vertex *mergeVertex(Vertex *v1, Vertex *v2);
	Vertex *findVertex(const Vector3D &position) const;
	Vertex *findNearestVertex(const Vector3D &position) const;
	void    findVertex(const Vector3D &centre, double radius, std::list<Vertex *> &vertices) const;
	void    findCommonVertices(const Face *face1, const Face *face2, VertexSet &vertices) const;
	void    removeDisconnVert();
	Vertex *isCorner(const Face *face, const Vector3D &position, float epsilon = FLT_EPSILON) const;

	void   updateVertexIDs();
	size_t getVertexID(const Vertex *vertex) const;

	void   updateFaceIDs();
	size_t getFaceID(const Face *face) const;

	Edge *addEdge(Vertex *fromVertex, Vertex *toVertex, double cost = 0);
	Face *removeEdge(Edge *edge);

	Edge *findNearestEdge(const Vector3D &position, Vector3D &nearestPoint, Vertex * &nearestVertex, float epsilon = 0.0) const;

	Edge *findCommonEdge(const Face   *face1, const Face   *face2) const;
	Edge *findCommonEdge(const Vertex *vertex1, const Vertex *vertex2) const;

	Face *addFace(const std::list<Vector3D> &p);
	void  removeFace(Face *face, bool reserveEdgeConnection);

	Face      *findNearestFace(const Vector3D &position, Vector3D &nearestPoint, Vertex * &nearestVertex, Edge * &nearestEdge, float epsilon = 0.0) const;
	Primitive *findNearestFace(const Vector3D &position, Vector3D &nearestPoint, Vector3D &normal, float epsilon = 0.0) const;

	void findCommonFaces(const Edge   *edge1, const Edge   *edge2, ConstFaceSet &faces) const;
	void findCommonFaces(const Edge   *edge, const Vertex *vertex, ConstFaceSet &faces) const;
	void findCommonFaces(const Vertex *vertex1, const Vertex *vertex2, ConstFaceSet &faces) const;
	void findCommonFaces(const Face   *face1, const Face   *face2, ConstFaceSet &faces) const;
	void findCommonFaces(const Face   *face, const Edge   *edge, ConstFaceSet &faces) const;
	void findCommonFaces(const Face   *face, const Vertex *vertex, ConstFaceSet &faces) const;
	void findCommonFaces(const Primitive *p1, const Primitive *p2, ConstFaceSet &faces) const;

	bool haveCommonFaces(const Edge   *edge1, const Edge   *edge2) const;
	bool haveCommonFaces(const Edge   *edge, const Vertex *vertex) const;
	bool haveCommonFaces(const Vertex *vertex1, const Vertex *vertex2) const;
	bool haveCommonFaces(const Face   *face1, const Face   *face2) const;
	bool haveCommonFaces(const Face   *face, const Edge   *edge) const;
	bool haveCommonFaces(const Face   *face, const Vertex *vertex) const;
	bool haveCommonFaces(const Primitive *p1, const Primitive *p2) const;

#ifdef _USE_DCEL_SEARCH_TREE
	void constructSearchTrees(int depth = 3);
#endif

	Face *rayTrace(const Vector3D &start, const Vector3D &direction, Vector3D &hitPoint, bool twoSide) const;

	void transverseEdges(
		const Vertex *centre, const Primitive *primitive1, const Primitive *primitive2,
		std::list<const Edge *> &edgeList) const;


	/***************************************************************************************************
	*
	* void DoublyConnectedEdgeList::subdivide
	*
	* usage: to connect a vertex to all the vertices of a face. This will triangulate a face
	*        if the face contains the vertex
	*
	***************************************************************************************************/

	void subdivide(Face *face, Vertex *vertex);

	double breathFirstSearch(Vertex *fromVertex, Vertex *toVertex, SearchPathList &searchPath) const;
	double breathFirstSearch(const Vector3D &from, const Vector3D &to, SearchPathList &searchPath) const;

	double computeShortestPath(Vertex *fromVertex, Vertex *toVertex, SearchPathList &searchPath, long int maximumNumberOfVertex = -1) const;
	double computeShortestPath(const Vector3D &from, const Vector3D &to, SearchPathList &searchPath, long int maximumNumberOfVertex = -1) const;

protected:

	//Edge *mergeEdge      (Edge *v1, Edge *v2);

	Vertex *findNearestVertexWithoutSearchTree(const Vector3D &position) const;

	Face *addFace(Edge *component);

	void transverseEdges(
		const Vertex *centre, const Edge *edge1, const Edge *edge2,
		std::list<const Edge *> &edgeList) const;

	VertexSet vertices;
	EdgeSet   edges;
	FaceSet   faces;

	InternalBoundingBox *faceBoundingBoxes;

#ifdef _USE_DCEL_SEARCH_TREE
	VPTree*  vpTree;
	DataSet* ds;
	Vertex** vArray;
#endif

	bool isTriangularMesh;
	bool isFaceIDUpToDate;
	bool isVertexIDUpToDate;
};

#endif
