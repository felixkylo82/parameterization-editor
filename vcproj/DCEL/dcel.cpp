#define _USE_MATH_DEFINES

#include "dcel.h"

#include "raytrace.h"

#include "objmodel.h"
#include "boundingbox.h"

#include <set>
#include <queue>
#include <algorithm>

#include <cmath>

#include <windows.h>
#include <GL/gl.h>

const unsigned int DoublyConnectedEdgeList::max_sides = 10;
const float magnification = 50.0f;


DoublyConnectedEdgeList::Face::Face(DoublyConnectedEdgeList::Edge *newComponent)
{
	this->component = newComponent;
}

void DoublyConnectedEdgeList::Face::flipFaceOrientation()
{
	DoublyConnectedEdgeList::Edge   *startEdge, *currentEdge, *nextEdge;
	DoublyConnectedEdgeList::Vertex *startVertex;

	std::set<DoublyConnectedEdgeList::Face*> faceSet;
	if (!this->visited)
	{
		faceSet.insert(this);
	}

	while (!faceSet.empty())
	{
		std::set<DoublyConnectedEdgeList::Face*>::iterator fiter = faceSet.begin();
		Face* face = *fiter;
		faceSet.erase(fiter);

		nextEdge = currentEdge = startEdge = face->component;

		if (NULL != nextEdge)
		{
			nextEdge = nextEdge->getNextEdge();
		}

		startVertex = startEdge->origin;

		if (NULL != nextEdge)
		{
			do
			{
				currentEdge->origin = nextEdge->origin;

				currentEdge->next = currentEdge->previous;
				currentEdge->previous = nextEdge;

				currentEdge = nextEdge;
				nextEdge = nextEdge->next;
			} while (nextEdge && nextEdge != startEdge);

			if (nextEdge)
			{
				currentEdge->origin = startVertex;

				currentEdge->next = currentEdge->previous;
				currentEdge->previous = nextEdge;
			}
		}

		face->visited = true;

		currentEdge = startEdge = face->component;

		do
		{
			if (NULL != currentEdge->twin->incidentFace)
			{
				DoublyConnectedEdgeList::Face* incidentFace = currentEdge->twin->incidentFace;
				if (!incidentFace->visited)
				{
					faceSet.insert(incidentFace);
				}
			}
			else
			{
				currentEdge->twin->origin = currentEdge->next->origin;
			}

			currentEdge = currentEdge->previous;
		} while (currentEdge && currentEdge != startEdge);
	}
}

const DoublyConnectedEdgeList::Edge *DoublyConnectedEdgeList::Face::getComponent() const
{
	return this->component;
}

void DoublyConnectedEdgeList::Face::getNeighboringFaces(DoublyConnectedEdgeList::ConstFaceSet &faces) const
{
	DoublyConnectedEdgeList::Edge *startEdge, *currentEdge, *twinEdge;
	DoublyConnectedEdgeList::Face *face;

#ifdef DCEL_DEBUG
	assert(faces.size() == 0);
#endif

	currentEdge = startEdge = this->component;

	do
	{
		twinEdge = currentEdge->twin;

		face = twinEdge->incidentFace;
		if (face)
		{
			faces.insert(face);
		}

		currentEdge = currentEdge->next;
	} while (currentEdge != startEdge);
}

void DoublyConnectedEdgeList::Face::getVertices(std::vector<const Vertex *> &v) const
{
	DoublyConnectedEdgeList::Edge *currentEdge, *startEdge;


	currentEdge = startEdge = this->component;

	do
	{
		v.push_back(currentEdge->origin);

		currentEdge = currentEdge->next;
	} while (startEdge != currentEdge);
}

DoublyConnectedEdgeList::Edge::Edge(DoublyConnectedEdgeList::Vertex *newOrigin, double newCost)
{
	this->origin = newOrigin;
	this->cost = newCost;

	this->next = NULL;
	this->previous = NULL;
	this->twin = NULL;
	this->incidentFace = NULL;
}

double DoublyConnectedEdgeList::Edge::getCost() const
{
	return this->cost;
}

DoublyConnectedEdgeList::Vertex* DoublyConnectedEdgeList::Edge::getOrigin() const
{
	return this->origin;
}

DoublyConnectedEdgeList::Edge* DoublyConnectedEdgeList::Edge::setNextEdge(DoublyConnectedEdgeList::Edge* next)
{
	DoublyConnectedEdgeList::Edge *ret;


	ret = NULL;

	if (NULL != next)
	{
		if (NULL == this->next)
		{
			next->previous = this;
			ret = this->next = next;
		}
	}
	else
	{
		if (NULL != this->next)
		{
			this->next->previous = NULL;
			this->next = NULL;
		}
	}


	return ret;
}

DoublyConnectedEdgeList::Edge* DoublyConnectedEdgeList::Edge::getNextEdge() const
{
	return this->next;
}

DoublyConnectedEdgeList::Edge* DoublyConnectedEdgeList::Edge::getPreviousEdge() const
{
	return this->previous;
}

DoublyConnectedEdgeList::Edge* DoublyConnectedEdgeList::Edge::getTwinEdge() const
{
	return this->twin;
}

DoublyConnectedEdgeList::Face* DoublyConnectedEdgeList::Edge::getIncidentFace() const
{
	return this->incidentFace;
}

Vector3D DoublyConnectedEdgeList::Edge::computeNormal() const
{
	Vector3D n, t;

	if (this->incidentFace) n += this->incidentFace->computeNormal();
	if (this->twin->incidentFace) n += this->twin->incidentFace->computeNormal();

	return n;
}

void DoublyConnectedEdgeList::Edge::render(float epsilon) const
{
	GLdouble projection[16];

	if (0.0 != epsilon)
	{
		glGetDoublev(GL_PROJECTION_MATRIX, projection);

		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		glTranslated(0.0, 0.0, -epsilon);
		glMultMatrixd(projection);
		glMatrixMode(GL_MODELVIEW);
	}

	glBegin(GL_LINES);

	glNormal3dv(this->origin->computeNormal().normalize().data());
	glVertex3dv(this->origin->position.data());

	glNormal3dv(this->twin->origin->computeNormal().normalize().data());
	glVertex3dv(this->twin->origin->position.data());

	glEnd();

	if (0.0 != epsilon)
	{
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		glMatrixMode(GL_MODELVIEW);
	}
}

DoublyConnectedEdgeList::Vertex::Vertex(const Vector3D &v)
{
	this->position = v;

	DoublyConnectedEdgeList::VertexReset().operator ()(this);
}

DoublyConnectedEdgeList::Vertex::~Vertex()
{
}

const Vector3D &DoublyConnectedEdgeList::Vertex::getPosition() const
{
	return this->position;
}

void DoublyConnectedEdgeList::Vertex::addNeighbor(DoublyConnectedEdgeList::Edge *edge)
{
	this->outgoingEdges.insert(edge);
}

size_t DoublyConnectedEdgeList::Vertex::getNumberOfNeighbors() const
{
	return this->outgoingEdges.size();
}

Vector3D DoublyConnectedEdgeList::Vertex::computeNormal() const
{
	DoublyConnectedEdgeList::EdgeSet::const_iterator eiter;

	Vector3D normal;


	for (eiter = this->outgoingEdges.begin(); eiter != this->outgoingEdges.end(); ++eiter)
	{
		if ((*eiter)->incidentFace)
		{
			normal += (*eiter)->incidentFace->computeNormal();
		}
	}


	return normal;
}

Vector3D DoublyConnectedEdgeList::Vertex::computeNormal(double edgeBias) const
{
	DoublyConnectedEdgeList::Edge *e[7], *eswp;
	Vector3D normal, axis;
	double length;
	bool bias;
	unsigned int i;
	DoublyConnectedEdgeList::EdgeSet::const_iterator eiter;


	normal = this->computeNormal();
	length = normal.length();

	if (length > Vector3D::epsilon)
	{
		normal /= length;

		if (0.0 != edgeBias)
		{
			i = 0;
			eiter = this->outgoingEdges.begin();

			while (i < 7U && i < this->outgoingEdges.size())
			{
				e[i] = *eiter;

				++i;
				++eiter;
			}

			if (2 == this->outgoingEdges.size())
			{
				bias = true;
			}
			else if (3 <= this->outgoingEdges.size() && 7 >= this->outgoingEdges.size())
			{
				bias = true;

				for (i = 0; i < this->outgoingEdges.size(); ++i)
				{
					if (NULL == e[i]->incidentFace || NULL == e[i]->twin->incidentFace)
					{
						e[0] = e[i];
						break;
					}
				}

				//if(this->outgoingEdges.size() == i)
				//{
				//    bias = false;
				//}

				for (++i; i < this->outgoingEdges.size(); ++i)
				{
					if (NULL == e[i]->incidentFace || NULL == e[i]->twin->incidentFace)
					{
						e[1] = e[i];
						break;
					}
				}

				if (this->outgoingEdges.size() <= i)
				{
					bias = false;
				}
			}
			else
			{
				bias = false;
			}

			if (bias)
			{
				if (NULL == e[1]->incidentFace)
				{
					eswp = e[1];
					e[1] = e[0];
					e[0] = eswp;
				}

				//e[0] = e[0]->twin;

				axis = e[1]->twin->origin->position - e[0]->twin->origin->position;
				length = axis.length();
				if (length > Vector3D::epsilon)
				{
					axis /= length;
					normal = normal.rotate(axis, edgeBias);
				}
			}
		}
	}


	return normal;
}

void DoublyConnectedEdgeList::Vertex::render(float epsilon) const
{
	GLdouble projection[16];


	if (0.0 != epsilon)
	{
		glGetDoublev(GL_PROJECTION_MATRIX, projection);

		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		glTranslated(0.0, 0.0, -epsilon);
		glMultMatrixd(projection);
		glMatrixMode(GL_MODELVIEW);
	}

	glBegin(GL_POINTS);

	glNormal3dv(this->computeNormal().normalize().data());
	glVertex3dv(this->position.data());

	glEnd();

	if (0.0 != epsilon)
	{
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		glMatrixMode(GL_MODELVIEW);
	}
}

double DoublyConnectedEdgeList::Vertex::getPathCost() const
{
	return this->pathCost;
}

const DoublyConnectedEdgeList::EdgeSet &DoublyConnectedEdgeList::Vertex::getOutgoingEdges() const
{
	return this->outgoingEdges;
}

bool DoublyConnectedEdgeList::CostComparison::operator()(const DoublyConnectedEdgeList::Vertex* x, const DoublyConnectedEdgeList::Vertex* y) const
{
	double xc, yc;

	xc = x->getPathCost();
	yc = y->getPathCost();

	return (xc < yc || xc == yc && x < y);
}

bool DoublyConnectedEdgeList::VertexComparison::operator()(const DoublyConnectedEdgeList::Vertex* x, const DoublyConnectedEdgeList::Vertex* y) const
{
	const double *xp, *yp;


	xp = x->position.data();
	yp = y->position.data();

	if (xp[0] < yp[0])
		return true;

	if (yp[0] < xp[0])
		return false;

	if (xp[1] < yp[1])
		return true;

	if (yp[1] < xp[1])
		return false;

	if (xp[2] < yp[2])
		return true;

	return false;
}

void DoublyConnectedEdgeList::VertexReset::operator()(DoublyConnectedEdgeList::Vertex* currentVertex) const
{
	currentVertex->visited = false;
	currentVertex->backtrace = NULL;
	currentVertex->pathCost = FLT_MAX;
}

bool DoublyConnectedEdgeList::EdgeComparison::operator()(const DoublyConnectedEdgeList::Edge* x, const DoublyConnectedEdgeList::Edge* y) const
{
	return (x->origin < y->origin || x->origin == y->origin && x->twin->origin < y->twin->origin);
}

DoublyConnectedEdgeList::SearchPathNode::SearchPathNode(const DoublyConnectedEdgeList::Primitive *at) : at(at)
{
	if (DoublyConnectedEdgeList::Primitive::vertex == at->getPrimitiveType())
	{
		this->position = ((DoublyConnectedEdgeList::Vertex *)at)->position;
	}
}

DoublyConnectedEdgeList::SearchPathNode::SearchPathNode(const DoublyConnectedEdgeList::Primitive *at, const Vector3D &position) : at(at), position(position)
{
}

DoublyConnectedEdgeList::DoublyConnectedEdgeList() :

#ifdef _USE_DCEL_SEARCH_TREE

faceBoundingBoxes(NULL),
//vertexBoundingBoxes(NULL)

vpTree(NULL),
ds(NULL),
vArray(NULL),

#endif

isTriangularMesh(false),
isFaceIDUpToDate(true),
isVertexIDUpToDate(true)

{
}

DoublyConnectedEdgeList::DoublyConnectedEdgeList(const OBJModel &model) :

#ifdef _USE_DCEL_SEARCH_TREE

faceBoundingBoxes(NULL),

//vertexBoundingBoxes(NULL)

vpTree(NULL),
ds(NULL),
vArray(NULL),

#endif

isTriangularMesh(false),
isFaceIDUpToDate(true),
isVertexIDUpToDate(true)

{
	this->convertFromOBJModel(model);
}

DoublyConnectedEdgeList::~DoublyConnectedEdgeList()
{
	//destruction here

	this->clear();
}

void DoublyConnectedEdgeList::convertFromOBJModel(const OBJModel &model)
{
	const PolygonElement *polygon;
	DoublyConnectedEdgeList::Face *f;
	DoublyConnectedEdgeList::Edge *currentEdge, *firstEdge;
	unsigned int numberOfFaces;

	std::list<Vector3D> p;

	unsigned int i, j;


	this->clear();

	numberOfFaces = (unsigned int)model.getNumberOfPolygons();

	for (i = 0; i < numberOfFaces; ++i)
	{
		polygon = model.getPolygon(i);


		for (j = 0; j < polygon->vertexIndices.size(); ++j)
		{
			p.push_back(model.getVertex(polygon->vertexIndices[j]));
		}

		f = this->addFace(p);

		if (!polygon->textureCoordinateIndices.empty())
		{
			j = 0;
			currentEdge = firstEdge = f->component;
			do
			{
				currentEdge->origin->texCoord = model.getTextureCoordinate(polygon->textureCoordinateIndices[j]);

				++j;
				currentEdge = currentEdge->next;
			} while (currentEdge != firstEdge);
		}

		p.clear();
	}
}

void DoublyConnectedEdgeList::clear()
{
	DoublyConnectedEdgeList::VertexSet::iterator vitr;
	DoublyConnectedEdgeList::EdgeSet::iterator eitr;
	DoublyConnectedEdgeList::FaceSet::iterator fitr;


	for (vitr = this->vertices.begin(); vitr != this->vertices.end(); ++vitr)
	{
		delete *vitr;
	}
	this->vertices.clear();

	for (eitr = this->edges.begin(); eitr != this->edges.end(); ++eitr)
	{
		delete *eitr;
	}
	this->edges.clear();

	for (fitr = this->faces.begin(); fitr != this->faces.end(); ++fitr)
	{
		delete *fitr;
	}
	this->faces.clear();

#ifdef _USE_DCEL_SEARCH_TREE

	if (NULL != this->faceBoundingBoxes)
	{
		delete this->faceBoundingBoxes;
		this->faceBoundingBoxes = NULL;
	}

	//if(NULL != this->vertexBoundingBoxes)
	//{
	//    delete this->vertexBoundingBoxes;
	//    this->vertexBoundingBoxes = NULL;
	//}

	if (NULL != this->vpTree)
	{
		freeVPTree(this->vpTree);

		this->vpTree = NULL;
	}

	if (NULL != this->ds)
	{
		freeDataSet(this->ds);
		this->ds = NULL;
	}

	if (NULL != this->vArray)
	{
		delete this->vArray;
		this->vArray = NULL;
	}

#endif
}

void DoublyConnectedEdgeList::getPolygonSide(unsigned int &min, unsigned int &max) const
{
	DoublyConnectedEdgeList::Edge *startEdge, *currentEdge;

	unsigned int count;

	DoublyConnectedEdgeList::FaceSet::const_iterator fitr;


	min = 0xFFFFFFFF;
	max = 0;
	count = 0;

	for (fitr = this->faces.begin(); fitr != this->faces.end(); ++fitr)
	{
		currentEdge = startEdge = (*fitr)->component;
		count = 0;

		do
		{
			currentEdge = currentEdge->next;
			++count;
		} while (currentEdge != startEdge && count < max_sides);

		if (count > max_sides)
		{
			throw "An invalid face found!";
		}

		if (min > count)
		{
			min = count;
		}

		if (max < count)
		{
			max = count;
		}
	}
}

bool DoublyConnectedEdgeList::flipFaceOrientation(DoublyConnectedEdgeList::Face *face)
{
	DoublyConnectedEdgeList::EdgeSet originalEdges;

	bool flipped;

	DoublyConnectedEdgeList::FaceSet::iterator fiter;
	DoublyConnectedEdgeList::VertexSet::const_iterator viter;
	DoublyConnectedEdgeList::EdgeSet::const_iterator eiter;


	if (NULL == face || this->faces.find(face) == this->faces.end())
	{
		flipped = false;
	}
	else
	{
		originalEdges = this->edges;

		this->edges.clear();

		for (viter = this->vertices.begin(); viter != this->vertices.end(); ++viter)
		{
			(*viter)->outgoingEdges.clear();
		}

		for (fiter = this->faces.begin(); fiter != this->faces.end(); ++fiter)
		{
			(*fiter)->visited = false;
		}

		face->flipFaceOrientation();

		for (eiter = originalEdges.begin(); eiter != originalEdges.end(); ++eiter)
		{
			this->edges.insert(*eiter);
			(*eiter)->origin->outgoingEdges.insert(*eiter);
		}

		flipped = true;
	}


	return flipped;
}

#if 0
void DoublyConnectedEdgeList::toCloseManifold()
{
	DoublyConnectedEdgeList::EdgeSet::iterator eiter;
	DoublyConnectedEdgeList::Edge *currentEdge;

	int count;


	for (eiter = this->edges.begin(); eiter != this->edges.end(); ++eiter)
	{
		if (NULL == (*eiter)->next)
		{
			currentEdge = (*eiter)->twin;

			while (currentEdge->incidentFace)
			{
				currentEdge = currentEdge->previous->twin;
			}

			(*eiter)->setNextEdge(currentEdge);
		}
	}

	for (eiter = this->edges.begin(); eiter != this->edges.end(); ++eiter)
	{
		if (NULL == (*eiter)->incidentFace)
		{
			currentEdge = *eiter;

			if (currentEdge->next)
			{
				count = 0;

				do
				{
					if (currentEdge->next == currentEdge)
					{
						std::cerr << "currentEdge->next == currentEdge" << std::endl;
					}
					currentEdge = currentEdge->next;
					++count;
				} while (currentEdge->next && currentEdge != *eiter && 100 > count);

				std::cout << count << std::endl;

				if (currentEdge == *eiter)
				{
					this->addFace(*eiter);
					std::cout << "a face added" << std::endl;
				}
			}
		}
	}
}
#endif

void DoublyConnectedEdgeList::fillMissingPolygons()
{
	DoublyConnectedEdgeList::EdgeSet::iterator eiter;
	DoublyConnectedEdgeList::Edge *currentEdge;

	int count;


	for (eiter = this->edges.begin(); eiter != this->edges.end(); ++eiter)
	{
		if (NULL == (*eiter)->next)
		{
			currentEdge = (*eiter)->twin;

			while (currentEdge->incidentFace)
			{
				currentEdge = currentEdge->previous->twin;
			}

			(*eiter)->setNextEdge(currentEdge);
		}
	}

	for (eiter = this->edges.begin(); eiter != this->edges.end(); ++eiter)
	{
		if (NULL == (*eiter)->incidentFace)
		{
			currentEdge = *eiter;

			if (currentEdge->next)
			{
				count = 0;

				do
				{
					currentEdge = currentEdge->next;
					++count;
				} while (currentEdge && currentEdge != *eiter && DoublyConnectedEdgeList::max_sides > count);

				if (currentEdge == *eiter)
				{
					this->addFace(*eiter);
				}
			}
		}
	}

	for (eiter = this->edges.begin(); eiter != this->edges.end(); ++eiter)
	{
		if (NULL == (*eiter)->incidentFace && NULL != (*eiter)->next)
		{
			(*eiter)->setNextEdge(NULL);
		}
	}
}

bool DoublyConnectedEdgeList::validate() const
{
	DoublyConnectedEdgeList::VertexSet::const_iterator vitr;
	DoublyConnectedEdgeList::EdgeSet::const_iterator eitr;
	DoublyConnectedEdgeList::FaceSet::const_iterator fitr;

	DoublyConnectedEdgeList::Edge *startEdge, *currentEdge;

	unsigned int minSide, maxSide;
	unsigned int count;

	int genus;

	bool valid;


	valid = true;

	try
	{
		// check face integrity

		this->getPolygonSide(minSide, maxSide);

		for (fitr = this->faces.begin(); fitr != this->faces.end(); ++fitr)
		{
			currentEdge = startEdge = (*fitr)->component;

			do
			{
				if (NULL == currentEdge)
				{
					throw "An invalid face found [NULL == currentEdge] !";
				}

				if (currentEdge->incidentFace != *fitr)
				{
					throw "An invalid face found [currentEdge->incidentFace != *fitr] !";
				}

				if (this->edges.find(currentEdge) == this->edges.end())
				{
					throw "An invalid face found [this->edges.find(currentEdge) == this->edges.end()] !";
				}

				currentEdge = currentEdge->next;
			} while (currentEdge != startEdge);
		}

		// check edge integrity

		if (0 != this->edges.size() % 2)
		{
			throw "An invalid edge found [0 != this->edges.size() % 2] !";
		}

		for (eitr = this->edges.begin(); eitr != this->edges.end(); ++eitr)
		{
			if (NULL == (*eitr)->twin)
			{
				throw "An invalid edge found [NULL == (*eitr)->twin] !";
			}

			if ((*eitr)->twin->twin != *eitr)
			{
				throw "An invalid edge found [(*eitr)->twin->twin != *eitr] !";
			}

			if ((*eitr)->twin == (*eitr))
			{
				throw "An invalid edge found [(*eitr)->twin == (*eitr)] !";
			}

			if ((*eitr)->twin->origin == (*eitr)->origin)
			{
				throw "An invalid edge found [(*eitr)->twin->origin == (*eitr)->origin] !";
			}

			if (NULL == *eitr && NULL == (*eitr)->twin)
			{
				throw "An invalid edge found [NULL == *eitr && NULL == (*eitr)->twin] !";
			}

			if (NULL != (*eitr)->next)
			{
				if (NULL == (*eitr)->next->previous)
				{
					throw "An invalid edge found [NULL == (*eitr)->next->previous] !";
				}

				if ((*eitr)->next->previous != *eitr)
				{
					throw "An invalid edge found [(*eitr)->next->previous != *eitr] !";
				}
			}

			if (this->vertices.find((*eitr)->origin) == this->vertices.end())
			{
				throw "An invalid edge found [this->vertices.find((*eitr)->origin) == this->vertices.end()] !";
			}
		}

		// outgoing edge count

		count = 0;

		for (vitr = this->vertices.begin(); vitr != this->vertices.end(); ++vitr)
		{
			if (0 == (unsigned int)(*vitr)->outgoingEdges.size())
			{
				throw "The outgoing edge of a vertex is 0";
			}

			count += (unsigned int)(*vitr)->outgoingEdges.size();
		}

		if (count != this->edges.size())
		{
			throw "The outgoing edge count and the total edge count do not match";
		}

		for (vitr = this->vertices.begin(); vitr != this->vertices.end(); ++vitr)
		{
			for (eitr = (*vitr)->outgoingEdges.begin(); eitr != (*vitr)->outgoingEdges.end(); ++eitr)
			{
				if ((*eitr)->origin != *vitr)
				{
					throw "An invalid outgoing edge found [(*eitr)->origin != *vitr]";
				}
			}
		}
	}
	catch (const char err[])
	{
		std::cerr << err << std::endl;

		valid = false;
	}

	if (valid)
	{
		genus = (int)this->vertices.size() - (int)this->edges.size() / 2 + (int)this->faces.size();
		//std::cout << "V = " << this->vertices.size()  << std::endl;
		//std::cout << "E = " << this->edges.size() / 2 << std::endl;
		//std::cout << "F = " << this->faces.size()     << std::endl;
		if (2 != genus)
		{
			std::cout << "V - E + F = " << genus << std::endl;
			//valid = false;
		}
	}


	return valid;
}

void DoublyConnectedEdgeList::render(float epsilon) const
{
	GLdouble projection[16];

	DoublyConnectedEdgeList::FaceSet::const_iterator fitr;


	if (0.0 != epsilon)
	{
		glGetDoublev(GL_PROJECTION_MATRIX, projection);

		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		glTranslated(0.0, 0.0, -epsilon);
		glMultMatrixd(projection);
		glMatrixMode(GL_MODELVIEW);
	}

	for (fitr = this->faces.begin(); fitr != this->faces.end(); ++fitr)
	{
		(*fitr)->render();
	}

	if (0.0 != epsilon)
	{
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		glMatrixMode(GL_MODELVIEW);
	}
}

void DoublyConnectedEdgeList::triangulate()
{
#if 0
	DoublyConnectedEdgeList::FaceSet originalFaces;

	DoublyConnectedEdgeList::FaceSet::const_iterator fitr;
	DoublyConnectedEdgeList::Edge *startEdge, *currentEdge, *nextEdge, *newEdge, *newTwinEdge, *lastEdge, *incidentEdge;

	int count;


	originalFaces = this->faces;

	this->faces.clear();

	for (fitr = originalFaces.begin(); fitr != originalFaces.end(); ++fitr)
	{
		startEdge = (*fitr)->component;
		lastEdge = startEdge->previous;
		currentEdge = startEdge;
		count = 0;

		do
		{
			nextEdge = currentEdge->next;

			if (count >= 2)
			{
				incidentEdge = lastEdge->next;

				newEdge = this->addEdge(incidentEdge->origin, currentEdge->origin, Vector3D::distance(incidentEdge->origin->position, currentEdge->origin->position));
				newTwinEdge = newEdge->twin;

				newEdge->next = incidentEdge;
				incidentEdge->previous = newEdge;

				newEdge->previous = currentEdge->previous;
				currentEdge->previous->next = newEdge;

				newTwinEdge->next = currentEdge;
				currentEdge->previous = newTwinEdge;

				newTwinEdge->previous = lastEdge;
				lastEdge->next = newTwinEdge;

				this->addFace(incidentEdge);
			}

			currentEdge = nextEdge;

			++count;
		} while (currentEdge != lastEdge);

		this->addFace(lastEdge->next);

		delete *fitr;
	}

	//originalFaces.clear();
#endif

	DoublyConnectedEdgeList::FaceSet originalFaces;
	DoublyConnectedEdgeList::FaceSet::const_iterator fitr;

	DoublyConnectedEdgeList::Edge *startEdge, *currentEdge;

	Vector3D average;

	unsigned int minSide, maxSide;

	unsigned int count;


	if (!this->isTriangularMesh)
	{
		this->getPolygonSide(minSide, maxSide);

		if (maxSide != 3 || minSide != 3)
		{
			originalFaces = this->faces;

			for (fitr = originalFaces.begin(); fitr != originalFaces.end(); ++fitr)
			{
				currentEdge = startEdge = (*fitr)->component;

				average.setZero();
				count = 0;

				do
				{
					average += currentEdge->origin->position;
					currentEdge = currentEdge->next;
					++count;
				} while (currentEdge != startEdge);

				average /= (double)count;

				this->subdivide(*fitr, this->addVertex(average));
			}
		}

		isTriangularMesh = true;
	}
}

bool DoublyConnectedEdgeList::loadFile(const char filename[])
{
	std::ifstream fin;

	bool ret;


	fin.open(filename);

	ret = this->loadFile(fin);

	if (fin.is_open())
	{
		fin.close();
	}

	return ret;
}

bool DoublyConnectedEdgeList::loadFile(std::ifstream &fin)
{
	char line[256];

	std::vector<DoublyConnectedEdgeList::Vertex *> vertexList;

	char *currentCharacter;
	std::vector<char *> tokens, tokens2;

	Vector3D     vector;
	unsigned int vertexID[2];

	DoublyConnectedEdgeList::Edge *startEdge, *previousEdge, *currentEdge;

	unsigned int i;

	bool ret;


	ret = true;

	try
	{
		if (!fin.is_open())
		{
			throw "file cannot be opened";
		}

		this->clear();

		while (fin.good())
		{
			fin.getline(line, 256);

			tokens.clear();

			currentCharacter = line;

			while (true)
			{
				while (*currentCharacter == ' ' || *currentCharacter == '\t' || *currentCharacter == '\n')
				{
					++currentCharacter;
				}

				if (*currentCharacter != '\0')
				{
					tokens.push_back(currentCharacter);

					do
					{
						++currentCharacter;
					} while (*currentCharacter != ' ' && *currentCharacter != '\t' && *currentCharacter != '\n' && *currentCharacter != '\0');
				}

				if (*currentCharacter != '\0')
				{
					*currentCharacter = '\0';
					++currentCharacter;
				}
				else
				{
					break;
				}
			}

			if (tokens.size() != 0)
			{
				if (strcmp(tokens[0], "v") == 0)
				{
					//////////////////////////////////////////////////
					// 1a) Read vertex information
					//////////////////////////////////////////////////

					if (tokens.size() == 4)
					{
						for (i = 0; i < 3; ++i)
						{
							vector[i] = atof(tokens[i + 1]);

							if (vector[i] == HUGE_VAL)
							{
								throw "floating point overflow";
							}
						}

						vertexList.push_back(this->addVertex(vector));
					}
					else
					{
						throw "invalid format";
					}
				}
				//else if(strcmp(tokens[0], "e") == 0)
				//{
				//    //////////////////////////////////////////////////
				//    // 1b) Read edge information
				//    //////////////////////////////////////////////////
				//
				//    if(tokens.size() == 3)
				//    {
				//        for(i = 0; i < 2; ++i)
				//        {
				//            vertexID[i] = atoi(tokens[i + 1]);
				//
				//            if(vertexID[i] == 0)
				//            {
				//                throw "invalid integers";
				//            }
				//        }

				//        this->addEdge(vertexList[vertexID[0]], vertexList[vertexID[1]],
				//            Vector3D::distance(vertexList[vertexID[0]]->position, vertexList[vertexID[1]]->position));
				//    }
				//    else
				//    {
				//        throw "invalid format";
				//    }
				//}
				else if (strcmp(tokens[0], "f") == 0)
				{
					//////////////////////////////////////////////////
					// 1c) Read face information
					//////////////////////////////////////////////////

					if (tokens.size() >= 4)
					{
						vertexID[0] = atoi(tokens[1]);

						if (vertexID[0] == 0)
						{
							throw "invalid integers";
						}

						--vertexID[0];

						vertexID[1] = atoi(tokens[2]);

						if (vertexID[1] == 0)
						{
							throw "invalid integers";
						}

						--vertexID[1];

						startEdge = this->addEdge(vertexList[vertexID[0]], vertexList[vertexID[1]],
							Vector3D::distance(vertexList[vertexID[0]]->position, vertexList[vertexID[1]]->position));

						currentEdge = startEdge;

						for (i = 3; i < tokens.size(); ++i)
						{
							vertexID[0] = vertexID[1];
							vertexID[1] = atoi(tokens[i]);

							previousEdge = currentEdge;

							if (vertexID[1] == 0)
							{
								throw "invalid integers";
							}

							--vertexID[1];

							currentEdge = this->addEdge(vertexList[vertexID[0]], vertexList[vertexID[1]],
								Vector3D::distance(vertexList[vertexID[0]]->position, vertexList[vertexID[1]]->position));

							previousEdge->setNextEdge(currentEdge);
						}

						vertexID[0] = vertexID[1];
						vertexID[1] = atoi(tokens[1]) - 1;

						previousEdge = currentEdge;

						currentEdge = this->addEdge(vertexList[vertexID[0]], vertexList[vertexID[1]],
							Vector3D::distance(vertexList[vertexID[0]]->position, vertexList[vertexID[1]]->position));

						previousEdge->setNextEdge(currentEdge);

						currentEdge->setNextEdge(startEdge);

						this->addFace(startEdge);
					}
					else
					{
						throw "invalid format";
					}
				}
				//else if(strcmp(tokens[0], "END_DCEL") == 0)
				//{
				//    // skip below
				//    break;
				//}
				else if (tokens[0][0] == '#')
				{
					// skip comment
				}
				//else if(strcmp(tokens[0], "g") == 0)
				//{
				//}
				else// if(strcmp(token[0], "\0") != 0)
				{
					throw "unhandled token";
				}
			}
		}
	}
	catch (const char error[])
	{
		this->clear();

		std::cerr << error << std::endl;

		ret = false;
	}

	vertexList.clear();

	//std::cout << this->faces.size()    << std::endl;
	//std::cout << this->edges.size()    << std::endl;
	//std::cout << this->vertices.size() << std::endl;


	return ret;
}

void DoublyConnectedEdgeList::saveFile(const char filename[]) const
{
	std::ofstream fout;


	fout.open(filename);

	this->saveFile(fout);

	if (fout.is_open())
	{
		fout.close();
	}
}

void DoublyConnectedEdgeList::convertToOBJModel(OBJModel &model, bool enableEdgeBias, bool outputTexCoord) const
{
	DoublyConnectedEdgeList::VertexSet::const_iterator viter;
	DoublyConnectedEdgeList::EdgeSet::const_iterator eiter;
	DoublyConnectedEdgeList::FaceSet::const_iterator fiter;

	DoublyConnectedEdgeList::Edge *startEdge, *currentEdge;

	PolygonElement polygon;

	unsigned int i;


	i = 1;

	for (viter = this->vertices.begin(); viter != this->vertices.end(); ++viter)
	{
		model.addVertex((*viter)->position);
		if (enableEdgeBias)
		{
			model.addNormal((*viter)->computeNormal(M_PI / 4.0).data());
		}
		else
		{
			model.addNormal((*viter)->computeNormal(0.0).data());
		}
		if (outputTexCoord)
		{
			model.addTextureCoordinate((*viter)->texCoord);
		}

		(*viter)->id = i;

		++i;
	}

	for (fiter = this->faces.begin(); fiter != this->faces.end(); ++fiter)
	{
		currentEdge = startEdge = (*fiter)->component;

		do
		{
			polygon.vertexIndices.push_back(currentEdge->origin->id - 1);
			polygon.normalIndices.push_back(currentEdge->origin->id - 1);
			if (outputTexCoord)
			{
				polygon.textureCoordinateIndices.push_back(currentEdge->origin->id - 1);
			}

			currentEdge = currentEdge->next;
		} while (currentEdge != startEdge);

		model.addPolygon(polygon);

		polygon.textureCoordinateIndices.clear();
		polygon.normalIndices.clear();
		polygon.vertexIndices.clear();
	}
}

void DoublyConnectedEdgeList::saveFile(std::ofstream &fout) const
{
	DoublyConnectedEdgeList::VertexSet::const_iterator viter;
	DoublyConnectedEdgeList::EdgeSet::const_iterator eiter;
	DoublyConnectedEdgeList::FaceSet::const_iterator fiter;

	DoublyConnectedEdgeList::Edge *startEdge, *currentEdge;

	unsigned int i, j;


	if (fout.is_open())
	{
		i = 1;

		for (viter = this->vertices.begin(); viter != this->vertices.end(); ++viter)
		{
			fout << 'v';

			for (j = 0; j < 3; ++j)
			{
				fout << '\t' << (*viter)->position[j];
			}

			fout << '\n';

			(*viter)->id = i;

			++i;
		}

		//for(eiter = this->edges.begin(); eiter != this->edges.end(); ++eiter)
		//{
		//    fout << 'e' << '\t' << (*eiter)->origin->id << '\t' <<  (*eiter)->next->origin->id << '\n';
		//}

		for (fiter = this->faces.begin(); fiter != this->faces.end(); ++fiter)
		{
			fout << 'f';

			currentEdge = startEdge = (*fiter)->component;

			do
			{
				fout << '\t' << currentEdge->origin->id;

				currentEdge = currentEdge->next;
			} while (currentEdge != startEdge);

			fout << '\n';
		}

		//fout << "END_DCEL" << std::endl;
	}
}

const DoublyConnectedEdgeList::VertexSet &DoublyConnectedEdgeList::getVertexSet() const
{
	return this->vertices;
}

const DoublyConnectedEdgeList::EdgeSet &DoublyConnectedEdgeList::getEdgeSet() const
{
	return this->edges;
}

const DoublyConnectedEdgeList::FaceSet &DoublyConnectedEdgeList::getFaceSet() const
{
	return this->faces;
}

void DoublyConnectedEdgeList::addVertex(DoublyConnectedEdgeList::Vertex *v)
{
	if (this->vertices.find(v) == this->vertices.end())
	{
		this->vertices.insert(v);
	}
}

DoublyConnectedEdgeList::Vertex *DoublyConnectedEdgeList::addVertex(const Vector3D &position)
{
	DoublyConnectedEdgeList::Vertex *v;
	DoublyConnectedEdgeList::VertexSet::const_iterator cvitr;


	v = new DoublyConnectedEdgeList::Vertex(position);
	cvitr = this->vertices.find(v);

	if (cvitr != this->vertices.end())
	{
		delete v;
		v = (*cvitr);
	}
	else
	{
		this->vertices.insert(v);

		this->isVertexIDUpToDate = false;
	}


	return v;
}


DoublyConnectedEdgeList::Face *DoublyConnectedEdgeList::removeVertex(const Vector3D &position)
{
	DoublyConnectedEdgeList::Vertex *v = this->findVertex(position);
	return this->removeVertex(v);
}

DoublyConnectedEdgeList::Face *DoublyConnectedEdgeList::removeVertex(Vertex *v)
{
	DoublyConnectedEdgeList::EdgeSet::iterator eitr;
	DoublyConnectedEdgeList::Face *face;

#ifdef DCEL_DEBUG
	std::cout << "the number of faces before removing vertex: " << this->faces.size() << std::endl;
#endif

	face = NULL;

	while (!v->outgoingEdges.empty())
	{
		face = this->removeEdge(*v->outgoingEdges.begin());
		//v->edges.erase(*v->edges.begin());
	}

	this->vertices.erase(v);

	delete v;


	this->isVertexIDUpToDate = false;

#ifdef DCEL_DEBUG
	std::cout << "the number of faces after removing vertex: " << this->faces.size() << std::endl;
#endif

	return face;
}

void DoublyConnectedEdgeList::removeDisconnVert()
{
	DoublyConnectedEdgeList::VertexSet originalVertexSet;
	DoublyConnectedEdgeList::VertexSet::iterator viter;


	originalVertexSet = this->vertices;

	for (viter = originalVertexSet.begin(); viter != originalVertexSet.end(); ++viter)
	{
		if (0 == (*viter)->getNumberOfNeighbors())
		{
			this->removeVertex(*viter);
		}
	}
}

void DoublyConnectedEdgeList::moveVertex(Vertex *v1, const Vector3D newPosition)
{
	this->vertices.erase(v1);

	v1->position = newPosition;

	this->vertices.insert(v1);

	this->isVertexIDUpToDate = false;
}

DoublyConnectedEdgeList::Vertex *DoublyConnectedEdgeList::mergeVertex(Vertex *v1, Vertex *v2)
{
	DoublyConnectedEdgeList::Edge *edge, *edge2;

	DoublyConnectedEdgeList::EdgeSet::const_iterator eiter, eiter2;


	while (v1->outgoingEdges.size() != 0)
	{
		edge = *v1->outgoingEdges.begin();

		this->edges.erase(edge);
		this->edges.erase(edge->twin);
		edge->origin->outgoingEdges.erase(edge);
		edge->twin->origin->outgoingEdges.erase(edge->twin);

		edge->origin = v2;

		eiter2 = this->edges.find(edge);

		if (eiter2 == this->edges.end())
		{
			this->edges.insert(edge);
			this->edges.insert(edge->twin);

			edge->origin->outgoingEdges.insert(edge);
			edge->twin->origin->outgoingEdges.insert(edge->twin);
		}
		else
		{
			edge2 = *eiter2;

			if (edge->next)
			{
				edge->next->previous = edge2;
				edge2->next = edge->next;
			}

			if (edge->twin->next)
			{
				edge->twin->next->previous = edge2->twin;
				edge2->twin->next = edge->twin->next;
			}

			if (edge->previous)
			{
				edge->previous->next = edge2;
				edge2->previous = edge->previous;
			}

			if (edge->twin->previous)
			{
				edge->twin->previous->next = edge2->twin;
				edge2->twin->previous = edge->twin->previous;
			}

			if (edge->incidentFace)
			{
				edge2->incidentFace = edge->incidentFace;
				edge2->incidentFace->component = edge2;
			}

			if (edge->twin->incidentFace)
			{
				edge2->twin->incidentFace = edge->twin->incidentFace;
				edge2->twin->incidentFace->component = edge2->twin;
			}

			delete edge->twin;
			delete edge;

			edge = edge2;
		}
	}

	this->vertices.erase(v1);

	delete v1;

	this->isVertexIDUpToDate = false;


	return v2;
}

DoublyConnectedEdgeList::Face *DoublyConnectedEdgeList::removeEdge(Edge* edge)
{
	DoublyConnectedEdgeList::Edge *component;
	DoublyConnectedEdgeList::Face *face;

	DoublyConnectedEdgeList::EdgeSet::iterator edgeIterator;

	edgeIterator = this->edges.find(edge);

	face = NULL;

	if (edgeIterator != edges.end())
	{
		if (edge->twin->incidentFace) this->removeFace(edge->twin->incidentFace, false);	// 25/6/06
		if (edge->incidentFace) this->removeFace(edge->incidentFace, false);	// 25/6/06

		this->edges.erase(edge->twin);
		this->edges.erase(edgeIterator);

		edge->twin->origin->outgoingEdges.erase(edge->twin);
		edge->origin->outgoingEdges.erase(edge);

		if (0 == edge->origin->outgoingEdges.size())
		{
			this->removeVertex(edge->origin);
		}

		if (0 == edge->twin->origin->outgoingEdges.size())
		{
			this->removeVertex(edge->twin->origin);
		}

		component = edge->next;

		if (edge->twin->previous)
		{
			edge->twin->previous->next = edge->next;
		}

		if (edge->next)
		{
			edge->next->previous = edge->twin->previous;
		}

		if (edge->previous)
		{
			edge->previous->next = edge->twin->next;
		}

		if (edge->twin->next)
		{
			edge->twin->next->previous = edge->previous;
		}

		delete edge->twin;
		delete edge;

		if (NULL != component)
		{
			face = this->addFace(component);
		}
		else
		{
			face = NULL;
		}
	}

	return face;
}

DoublyConnectedEdgeList::Edge *DoublyConnectedEdgeList::findCommonEdge(
	const DoublyConnectedEdgeList::Face *face1,
	const DoublyConnectedEdgeList::Face *face2) const
{
	DoublyConnectedEdgeList::Edge *startEdge1, *currentEdge1;
	DoublyConnectedEdgeList::Edge *startEdge2, *currentEdge2;
	DoublyConnectedEdgeList::Edge *commonEdge;

#ifdef DCEL_DEBUG
	assert(face1 != NULL);
	assert(face2 != NULL);
	assert(face1 != face2);
#endif

	commonEdge = NULL;

	startEdge1 = face1->component;
	startEdge2 = face2->component;

	currentEdge1 = startEdge1;

	do
	{
		currentEdge2 = startEdge2;

		do
		{
			if (currentEdge1 == currentEdge2->twin)
			{
				commonEdge = currentEdge1;
				break;
			}

			currentEdge2 = currentEdge2->next;
		} while (currentEdge2 != startEdge2);

		currentEdge1 = currentEdge1->next;
	} while (!commonEdge && currentEdge1 != startEdge1);

	return commonEdge;
}

DoublyConnectedEdgeList::Edge *DoublyConnectedEdgeList::findCommonEdge(
	const DoublyConnectedEdgeList::Vertex *vertex1,
	const DoublyConnectedEdgeList::Vertex *vertex2) const
{
	DoublyConnectedEdgeList::Edge *commonEdge;

	DoublyConnectedEdgeList::EdgeSet::const_iterator eiter;

#ifdef DCEL_DEBUG
	assert(vertex1 != NULL);
	assert(vertex2 != NULL);
	assert(vertex1 != vertex2);
#endif

	commonEdge = NULL;

	for (eiter = vertex1->outgoingEdges.begin(); NULL == commonEdge && eiter != vertex1->outgoingEdges.end(); ++eiter)
	{
		if ((*eiter)->twin->origin == vertex2)
		{
			commonEdge = *eiter;
		}
	}

	return commonEdge;
}

DoublyConnectedEdgeList::Edge *DoublyConnectedEdgeList::addEdge(
	DoublyConnectedEdgeList::Vertex *fromVertex,
	DoublyConnectedEdgeList::Vertex *toVertex,
	double cost)
{
	DoublyConnectedEdgeList::Edge *edge, *twinEdge;
	DoublyConnectedEdgeList::EdgeSet::const_iterator ceitr, ceitr2;

	edge = new DoublyConnectedEdgeList::Edge(fromVertex, cost);
	twinEdge = new DoublyConnectedEdgeList::Edge(toVertex, cost);

	edge->twin = twinEdge;
	twinEdge->twin = edge;

	ceitr = this->edges.find(edge);

	if (ceitr != this->edges.end())
	{
		delete edge;
		delete twinEdge;

		edge = (*ceitr);
	}
	else
	{
		this->edges.insert(edge);
		this->edges.insert(twinEdge);

		fromVertex->addNeighbor(edge);
		toVertex->addNeighbor(twinEdge);
	}


	return edge;
}

DoublyConnectedEdgeList::Face *DoublyConnectedEdgeList::addFace(DoublyConnectedEdgeList::Edge *newComponent)
{
	DoublyConnectedEdgeList::Face *face;
	DoublyConnectedEdgeList::FaceSet::const_iterator cfitr;
	DoublyConnectedEdgeList::Edge *currentEdge;

	bool isConsistent;


	isConsistent = true;

	face = new DoublyConnectedEdgeList::Face(newComponent);

	cfitr = this->faces.find(face);
	if (cfitr != this->faces.end())
	{
		isConsistent = false;
	}
	else
	{
		currentEdge = newComponent;

		do
		{
			if (NULL != currentEdge && NULL == currentEdge->incidentFace)
			{
				currentEdge = currentEdge->next;
			}
			else
			{
				isConsistent = false;
			}
		} while (isConsistent && currentEdge != newComponent);
	}

	if (isConsistent)
	{
		currentEdge = newComponent;

		do
		{
			currentEdge->incidentFace = face;
			currentEdge = currentEdge->next;
		} while (currentEdge != newComponent);

		this->faces.insert(face);

		this->isFaceIDUpToDate = false;
	}
	else
	{
		delete face;
		face = NULL;
	}


	return face;
}

DoublyConnectedEdgeList::Face *DoublyConnectedEdgeList::addFace(const std::list<Vector3D> &p)
{
	std::vector<DoublyConnectedEdgeList::Vertex *>  v;
	std::vector<DoublyConnectedEdgeList::Edge   *>  e;
	DoublyConnectedEdgeList::Face                  *f;

	std::vector<double> edgeCost;

	std::list<Vector3D>::const_iterator piter;

	unsigned int j;


	f = NULL;

	for (piter = p.begin(); piter != p.end(); ++piter)
	{
		v.push_back(this->addVertex(*piter));
	}

	try
	{
		for (j = 0; j < v.size(); ++j)
		{
			edgeCost.push_back(Vector3D::distance(v[j]->getPosition(), v[(j + 1) % v.size()]->getPosition()));

			if (0.0 == edgeCost[j])
			{
				throw "Duplicated points in a polygon";
			}
		}

		for (j = 0; j < v.size(); ++j)
		{
			e.push_back(this->addEdge(v[j], v[(j + 1) % v.size()], edgeCost[j]));

			if (NULL != e.back()->incidentFace)
			{
				if (NULL != e.back()->twin->incidentFace)
				{
					throw "An invalid polygon found";
				}

				this->flipFaceOrientation(e[j]->incidentFace);

				e.back() = e.back()->twin;
				if (NULL != e.back()->incidentFace)
				{
					fprintf(stderr, "critical error at LINE %d in FILE %s!\n", __LINE__, __FILE__);
					exit(1);
				}
			}
		}

		for (j = 0; j < e.size(); ++j)
		{
			e[j]->setNextEdge(e[(j + 1) % e.size()]);
		}

		f = this->addFace(e[0]);

		if (NULL == f)
		{
			throw "NULL == f";
		}
	}
	catch (const char err[])
	{
		for (j = 0; j < e.size(); ++j)
		{
			if (NULL != e[j] && NULL == e[j]->incidentFace && NULL == e[j]->twin->incidentFace)
			{
				this->removeEdge(e[j]);
			}
		}

		for (j = 0; j < v.size(); ++j)
		{
			if (NULL != v[j] && v[j]->outgoingEdges.size() == 0)
			{
				this->removeVertex(v[j]);
			}
		}

		std::cout << err << std::endl;
	}


	return f;
}

void DoublyConnectedEdgeList::removeFace(DoublyConnectedEdgeList::Face *face, bool reserveEdgeConnection)
{
	DoublyConnectedEdgeList::FaceSet::iterator faceIterater;
	DoublyConnectedEdgeList::Edge *nextEdge, *currentEdge, *startEdge;


	faceIterater = this->faces.find(face);
	if (faceIterater != this->faces.end())
	{
		this->faces.erase(faceIterater);

		currentEdge = startEdge = face->component;

		do
		{
			nextEdge = currentEdge->next;

			currentEdge->incidentFace = NULL;

			if (!reserveEdgeConnection)
			{
				currentEdge->setNextEdge(NULL);

				if (NULL == currentEdge->incidentFace && NULL == currentEdge->twin->incidentFace)
				{
					this->removeEdge(currentEdge);
				}
			}

			currentEdge = nextEdge;
		} while (currentEdge && currentEdge != startEdge);

		delete face;

		this->isFaceIDUpToDate = false;
	}
}

void DoublyConnectedEdgeList::findCommonFaces(const DoublyConnectedEdgeList::Edge *edge1, const DoublyConnectedEdgeList::Edge *edge2, ConstFaceSet &faces) const
{
#ifdef DCEL_DEBUG
	assert(edge1->twin != NULL);
	assert(edge2->twin != NULL);
#endif

	const DoublyConnectedEdgeList::Edge *currentEdge;

	currentEdge = edge1;

	do
	{
		if (currentEdge == edge2 ||
			currentEdge == edge2->twin)
		{
			faces.insert(currentEdge->incidentFace);
			break;
		}

		currentEdge = currentEdge->next;
	} while (currentEdge && currentEdge != edge1);

	currentEdge = edge1->twin;

	do
	{
		if (currentEdge == edge2 ||
			currentEdge == edge2->twin)
		{
			faces.insert(currentEdge->incidentFace);
			break;
		}

		currentEdge = currentEdge->next;
	} while (currentEdge && currentEdge != edge1->twin);
}

void DoublyConnectedEdgeList::findCommonFaces(const DoublyConnectedEdgeList::Edge *edge, const DoublyConnectedEdgeList::Vertex *vertex, ConstFaceSet &faces) const
{
#ifdef DCEL_DEBUG
	assert(edge);
	assert(vertex);
#endif

	DoublyConnectedEdgeList::EdgeSet::const_iterator eiter;
	DoublyConnectedEdgeList::Edge *currentEdge;


	for (eiter = vertex->outgoingEdges.begin(); eiter != vertex->outgoingEdges.end(); ++eiter)
	{
		currentEdge = *eiter;

		do
		{
			if (currentEdge == edge ||
				currentEdge == edge->twin)
			{
				faces.insert(currentEdge->incidentFace);
				break;
			}

			currentEdge = currentEdge->next;
		} while (currentEdge && currentEdge != *eiter);
	}
}

void DoublyConnectedEdgeList::findCommonFaces(
	const DoublyConnectedEdgeList::Vertex *vertex1,
	const DoublyConnectedEdgeList::Vertex *vertex2,
	ConstFaceSet & faces) const
{

#ifdef DCEL_DEBUG
	assert(vertex1);
	assert(vertex2);
#endif

	DoublyConnectedEdgeList::EdgeSet::const_iterator endEdge1 = vertex1->outgoingEdges.end();
	DoublyConnectedEdgeList::EdgeSet::const_iterator startEdge2 = vertex2->outgoingEdges.begin();
	DoublyConnectedEdgeList::EdgeSet::const_iterator endEdge2 = vertex2->outgoingEdges.end();

	for (DoublyConnectedEdgeList::EdgeSet::const_iterator currentEdge1 = vertex1->outgoingEdges.begin(); currentEdge1 != endEdge1; ++currentEdge1)
	{
		for (DoublyConnectedEdgeList::EdgeSet::const_iterator currentEdge2 = startEdge2; currentEdge2 != endEdge2; ++currentEdge2)
		{
			this->findCommonFaces(*currentEdge1, *currentEdge2, faces);
		}
	}
}


void DoublyConnectedEdgeList::findCommonFaces(const DoublyConnectedEdgeList::Face *face1, const DoublyConnectedEdgeList::Face *face2, ConstFaceSet &faces) const
{
	ConstFaceSet facet1;
	ConstFaceSet facet2;

#ifdef DCEL_DEBUG
	assert(face1 != face2);
	assert(face1);
	assert(face2);
#endif

	if (face1 == face2)
	{
		faces.insert(face1);
	}
}

void DoublyConnectedEdgeList::findCommonFaces(const Face *face, const Edge *edge, ConstFaceSet &faces) const
{
	DoublyConnectedEdgeList::Edge *startEdge;
	DoublyConnectedEdgeList::Edge *currentEdge;
	DoublyConnectedEdgeList::Edge *edgeTwin;

#ifdef DCEL_DEBUG
	assert(edge);
	assert(edge->twin);
	assert(face);
	assert(face->component);
#endif

	//faces.clear();

	startEdge = face->component;
	currentEdge = startEdge;
	edgeTwin = edge->twin;

	do
	{
		if (currentEdge == edge || currentEdge == edgeTwin)
		{
			faces.insert(face);
			break;
		}

		currentEdge = currentEdge->next;
	} while (currentEdge != startEdge);
}

void DoublyConnectedEdgeList::findCommonFaces(const Face *face, const Vertex *vertex, ConstFaceSet &faces) const
{
	DoublyConnectedEdgeList::Edge *startEdge;
	DoublyConnectedEdgeList::Edge *currentEdge;

#ifdef DCEL_DEBUG
	assert(edge);
	assert(edge->twin);
	assert(face);
	assert(face->component);
#endif

	//faces.clear();

	startEdge = face->component;
	currentEdge = startEdge;

	do
	{
		if (currentEdge->origin == vertex)
		{
			faces.insert(face);
			break;
		}
		currentEdge = currentEdge->next;
	} while (currentEdge != startEdge);
}

void DoublyConnectedEdgeList::findCommonFaces(const DoublyConnectedEdgeList::Primitive *p1, const DoublyConnectedEdgeList::Primitive *p2, ConstFaceSet &faces) const
{
	switch (p1->getPrimitiveType())
	{
	case DoublyConnectedEdgeList::Primitive::vertex:
		switch (p2->getPrimitiveType())
		{
		case DoublyConnectedEdgeList::Primitive::vertex:
			this->findCommonFaces((DoublyConnectedEdgeList::Vertex *)p1, (DoublyConnectedEdgeList::Vertex *)p2, faces);
			break;

		case DoublyConnectedEdgeList::Primitive::edge:
			this->findCommonFaces((DoublyConnectedEdgeList::Edge *)p2, (DoublyConnectedEdgeList::Vertex *)p1, faces);
			break;

		case DoublyConnectedEdgeList::Primitive::face:
			this->findCommonFaces((DoublyConnectedEdgeList::Face *)p2, (DoublyConnectedEdgeList::Vertex *)p1, faces);
			break;
		}
		break;

	case DoublyConnectedEdgeList::Primitive::edge:
		switch (p2->getPrimitiveType())
		{
		case DoublyConnectedEdgeList::Primitive::vertex:
			this->findCommonFaces((DoublyConnectedEdgeList::Edge *)p1, (DoublyConnectedEdgeList::Vertex *)p2, faces);
			break;

		case DoublyConnectedEdgeList::Primitive::edge:
			this->findCommonFaces((DoublyConnectedEdgeList::Edge *)p1, (DoublyConnectedEdgeList::Edge *)p2, faces);
			break;

		case DoublyConnectedEdgeList::Primitive::face:
			this->findCommonFaces((DoublyConnectedEdgeList::Face *)p2, (DoublyConnectedEdgeList::Edge *)p1, faces);
			break;
		}
		break;

	case DoublyConnectedEdgeList::Primitive::face:
		switch (p2->getPrimitiveType())
		{
		case DoublyConnectedEdgeList::Primitive::vertex:
			this->findCommonFaces((DoublyConnectedEdgeList::Face *)p1, (DoublyConnectedEdgeList::Vertex *)p2, faces);
			break;

		case DoublyConnectedEdgeList::Primitive::edge:
			this->findCommonFaces((DoublyConnectedEdgeList::Face *)p1, (DoublyConnectedEdgeList::Edge *)p2, faces);
			break;

		case DoublyConnectedEdgeList::Primitive::face:
			this->findCommonFaces((DoublyConnectedEdgeList::Face *)p1, (DoublyConnectedEdgeList::Face *)p2, faces);
			break;
		}
		break;
	}
}

bool DoublyConnectedEdgeList::haveCommonFaces(const DoublyConnectedEdgeList::Edge *edge1, const DoublyConnectedEdgeList::Edge *edge2) const
{
#ifdef DCEL_DEBUG
	assert(edge1->twin != NULL);
	assert(edge2->twin != NULL);
#endif

	const DoublyConnectedEdgeList::Edge *currentEdge;

	currentEdge = edge1;

	do
	{
		if (currentEdge == edge2 ||
			currentEdge == edge2->twin)
		{
			return true;
		}

		currentEdge = currentEdge->next;
	} while (currentEdge && currentEdge != edge1);

	currentEdge = edge1->twin;

	do
	{
		if (currentEdge == edge2 ||
			currentEdge == edge2->twin)
		{
			return true;
		}

		currentEdge = currentEdge->next;
	} while (currentEdge && currentEdge != edge1->twin);


	return false;
}


bool DoublyConnectedEdgeList::haveCommonFaces(const DoublyConnectedEdgeList::Edge *edge, const DoublyConnectedEdgeList::Vertex *vertex) const
{
#ifdef DCEL_DEBUG
	assert(edge);
	assert(vertex);
#endif

	DoublyConnectedEdgeList::EdgeSet::const_iterator eiter;
	DoublyConnectedEdgeList::Edge *currentEdge;


	for (eiter = vertex->outgoingEdges.begin(); eiter != vertex->outgoingEdges.end(); ++eiter)
	{
		currentEdge = *eiter;

		do
		{
			if (currentEdge == edge ||
				currentEdge == edge->twin)
			{
				return true;
				break;
			}

			currentEdge = currentEdge->next;
		} while (currentEdge && currentEdge != *eiter);
	}


	return false;
}

bool DoublyConnectedEdgeList::haveCommonFaces(
	const DoublyConnectedEdgeList::Vertex *vertex1,
	const DoublyConnectedEdgeList::Vertex *vertex2) const
{
	DoublyConnectedEdgeList::EdgeSet::const_iterator currentEdge1, currentEdge2;


#ifdef DCEL_DEBUG
	assert(vertex1);
	assert(vertex2);
#endif

	for (currentEdge1 = vertex1->outgoingEdges.begin(); currentEdge1 != vertex1->outgoingEdges.end(); ++currentEdge1)
	{
		for (currentEdge2 = vertex2->outgoingEdges.begin(); currentEdge2 != vertex2->outgoingEdges.end(); ++currentEdge2)
		{
			if (this->haveCommonFaces(*currentEdge1, *currentEdge2))
			{
				return true;
			}
		}
	}


	return false;
}


bool DoublyConnectedEdgeList::haveCommonFaces(const DoublyConnectedEdgeList::Face *face1, const DoublyConnectedEdgeList::Face *face2) const
{
	ConstFaceSet facet1;
	ConstFaceSet facet2;

#ifdef DCEL_DEBUG
	assert(face1 != face2);
	assert(face1);
	assert(face2);
#endif

	//faces.clear();

	if (face1 == face2)
	{
		return true;
	}

	return false;
}

bool DoublyConnectedEdgeList::haveCommonFaces(const Face *face, const Edge *edge) const
{
	DoublyConnectedEdgeList::Edge *startEdge;
	DoublyConnectedEdgeList::Edge *currentEdge;
	DoublyConnectedEdgeList::Edge *edgeTwin;

#ifdef DCEL_DEBUG
	assert(edge);
	assert(edge->twin);
	assert(face);
	assert(face->component);
#endif

	//faces.clear();

	startEdge = face->component;
	currentEdge = startEdge;
	edgeTwin = edge->twin;

	do
	{
		if (currentEdge == edge || currentEdge == edgeTwin)
		{
			return true;
		}

		currentEdge = currentEdge->next;
	} while (currentEdge != startEdge);


	return false;
}

bool DoublyConnectedEdgeList::haveCommonFaces(const Face *face, const Vertex *vertex) const
{
	DoublyConnectedEdgeList::Edge *startEdge;
	DoublyConnectedEdgeList::Edge *currentEdge;

#ifdef DCEL_DEBUG
	assert(edge);
	assert(edge->twin);
	assert(face);
	assert(face->component);
#endif

	//faces.clear();

	startEdge = face->component;
	currentEdge = startEdge;

	do
	{
		if (currentEdge->origin == vertex)
		{
			return true;
		}
		currentEdge = currentEdge->next;
	} while (currentEdge != startEdge);


	return false;
}

bool DoublyConnectedEdgeList::haveCommonFaces(const DoublyConnectedEdgeList::Primitive *p1, const DoublyConnectedEdgeList::Primitive *p2) const
{
	bool have;


	have = false;

	switch (p1->getPrimitiveType())
	{
	case DoublyConnectedEdgeList::Primitive::vertex:
		switch (p2->getPrimitiveType())
		{
		case DoublyConnectedEdgeList::Primitive::vertex:
			have = this->haveCommonFaces((DoublyConnectedEdgeList::Vertex *)p1, (DoublyConnectedEdgeList::Vertex *)p2);
			break;

		case DoublyConnectedEdgeList::Primitive::edge:
			have = this->haveCommonFaces((DoublyConnectedEdgeList::Edge *)p2, (DoublyConnectedEdgeList::Vertex *)p1);
			break;

		case DoublyConnectedEdgeList::Primitive::face:
			have = this->haveCommonFaces((DoublyConnectedEdgeList::Face *)p2, (DoublyConnectedEdgeList::Vertex *)p1);
			break;
		}
		break;

	case DoublyConnectedEdgeList::Primitive::edge:
		switch (p2->getPrimitiveType())
		{
		case DoublyConnectedEdgeList::Primitive::vertex:
			have = this->haveCommonFaces((DoublyConnectedEdgeList::Edge *)p1, (DoublyConnectedEdgeList::Vertex *)p2);
			break;

		case DoublyConnectedEdgeList::Primitive::edge:
			have = this->haveCommonFaces((DoublyConnectedEdgeList::Edge *)p1, (DoublyConnectedEdgeList::Edge *)p2);
			break;

		case DoublyConnectedEdgeList::Primitive::face:
			have = this->haveCommonFaces((DoublyConnectedEdgeList::Face *)p2, (DoublyConnectedEdgeList::Edge *)p1);
			break;
		}
		break;

	case DoublyConnectedEdgeList::Primitive::face:
		switch (p2->getPrimitiveType())
		{
		case DoublyConnectedEdgeList::Primitive::vertex:
			have = this->haveCommonFaces((DoublyConnectedEdgeList::Face *)p1, (DoublyConnectedEdgeList::Vertex *)p2);
			break;

		case DoublyConnectedEdgeList::Primitive::edge:
			have = this->haveCommonFaces((DoublyConnectedEdgeList::Face *)p1, (DoublyConnectedEdgeList::Edge *)p2);
			break;

		case DoublyConnectedEdgeList::Primitive::face:
			have = this->haveCommonFaces((DoublyConnectedEdgeList::Face *)p1, (DoublyConnectedEdgeList::Face *)p2);
			break;
		}
		break;
	}


	return have;
}

#ifdef _USE_DCEL_SEARCH_TREE

void DoublyConnectedEdgeList::constructSearchTrees(int depth)
{
	Vector3D min, max;

	DoublyConnectedEdgeList::VertexSet::const_iterator viter;
	DoublyConnectedEdgeList::FaceSet::const_iterator fiter;

	DoublyConnectedEdgeList::Edge *e1, *e2, *e3;

	unsigned int i;


	if (this->isTriangularMesh && !this->vertices.empty())
	{
		viter = this->vertices.begin();

		min = max = (*viter)->position;

		if (viter != this->vertices.end())
		{
			++viter;
		}

		while (viter != this->vertices.end())
		{
			for (i = 0; i < 3; ++i)
			{
				if ((*viter)->position[i] > max[i])
				{
					max[i] = (*viter)->position[i];
				}
				else if ((*viter)->position[i] < min[i])
				{
					min[i] = (*viter)->position[i];
				}
			}

			++viter;
		}

		if (this->faceBoundingBoxes)
		{
			delete this->faceBoundingBoxes;
		}

		this->faceBoundingBoxes = new InternalBoundingBox(min, max, depth);

		for (fiter = this->faces.begin(); fiter != this->faces.end(); ++fiter)
		{
			e1 = (*fiter)->component;
			if (NULL == e1) continue;
			e2 = e1->next;
			if (NULL == e2) continue;
			e3 = e2->next;
			if (NULL == e3) continue;

			this->faceBoundingBoxes->putDatum(Box(
				Vector3D::minimum(Vector3D::minimum(e1->origin->position, e2->origin->position), e3->origin->position),
				Vector3D::maximum(Vector3D::maximum(e1->origin->position, e2->origin->position), e3->origin->position)), *fiter);
		}

		if (NULL != this->vpTree)
		{
			freeVPTree(this->vpTree);

			this->vpTree = NULL;
		}

		if (NULL != this->ds)
		{
			freeDataSet(this->ds);
			this->ds = NULL;
		}

		if (NULL != this->vArray)
		{
			delete this->vArray;

			this->vArray = NULL;
		}

		// TODO: build the vp tree here

		this->ds = mallocDataSet((unsigned int)this->vertices.size(), 3);
		this->vArray = new Vertex*[this->vertices.size()];

		try
		{
			if (NULL == this->ds || NULL == this->vArray)
			{
				throw "Not enough memory!";
			}

			i = 0;

			for (viter = this->vertices.begin(); viter != this->vertices.end(); ++viter)
			{
				memcpy(this->ds->points + 3 * i, (*viter)->position.data(), 3 * sizeof(double));

				this->vArray[i] = *viter;

				++i;
			}


			this->vpTree = buildVPTree(this->ds, 1); // use defualt branch
			if (NULL == this->vpTree)
			{
				throw "Not enough memory!";
			}
		}
		catch (const char err[])
		{
			if (NULL != this->vpTree)
			{
				freeVPTree(this->vpTree);

				this->vpTree = NULL;
			}

			if (NULL != this->ds)
			{
				freeDataSet(this->ds);

				this->ds = NULL;
			}

			if (NULL != this->vArray)
			{
				delete this->vArray;

				this->vArray = NULL;
			}

			std::cerr << err << std::endl;
		}
	}
}

#endif

void internalRayTrace(void *datum, const Vector3D &start, const Vector3D &direction, void **nearestDatum, Vector3D &hitPoint, double &minDist, bool twoSide)
{
	DoublyConnectedEdgeList::Face *face;

	const DoublyConnectedEdgeList::Edge *edge1, *edge2, *edge3;

	double currentDist;


	face = (DoublyConnectedEdgeList::Face *)datum;

	edge1 = face->getComponent();
	if (NULL == edge1) return;
	edge2 = edge1->getNextEdge();
	if (NULL == edge2) return;
	edge3 = edge2->getNextEdge();
	if (NULL == edge3) return;

	currentDist = minDist;

	RayTrace::rayTraceTriangle(
		edge1->getOrigin()->getPosition(), edge2->getOrigin()->getPosition(), edge3->getOrigin()->getPosition(),
		start, direction, hitPoint, minDist, twoSide);

	if (minDist < currentDist)
	{
		nearestDatum = &datum;
	}
}

DoublyConnectedEdgeList::Face *DoublyConnectedEdgeList::rayTrace(const Vector3D &start, const Vector3D &direction, Vector3D &hitPoint, bool twoSide) const
{
	DoublyConnectedEdgeList::Face *face;

	double minDist;


	face = NULL;

	if (this->isTriangularMesh && NULL != this->faceBoundingBoxes)
	{
		minDist = FLT_MAX;

		this->faceBoundingBoxes->rayTraceTwoSide = twoSide;
		this->faceBoundingBoxes->rayTraceFunction = internalRayTrace;

		this->faceBoundingBoxes->rayTrace(start, direction, (void **)&face, hitPoint, minDist);
	}


	return face;
}

void DoublyConnectedEdgeList::transverseEdges(
	const Vertex *centre, const Primitive *primitive1, const Primitive *primitive2,
	std::list<const Edge *> &edgeList) const
{
	const DoublyConnectedEdgeList::Edge *edge[2];

	const DoublyConnectedEdgeList::Primitive *currentPrimitive;

	const DoublyConnectedEdgeList::Edge *currentEdge, *startEdge;

	DoublyConnectedEdgeList::EdgeSet::const_iterator iter;

	unsigned int i, j;


	for (i = 0; i < 2; ++i)
	{
		switch (i)
		{
		case 0:
			currentPrimitive = primitive1;
			break;

		case 1:
			currentPrimitive = primitive2;
			break;
		}

		edge[i] = NULL;

		if (DoublyConnectedEdgeList::Primitive::vertex == currentPrimitive->getPrimitiveType())
		{
			for (iter = centre->outgoingEdges.begin(); NULL == edge[i] && iter != centre->outgoingEdges.end(); ++iter)
			{
				currentEdge = (*iter);

				if (currentEdge->twin->origin == currentPrimitive)
				{
					edge[i] = currentEdge;
				}
			}

			//std::cout << "Vertex" << std::endl;
		}
		else
		{
			for (j = 0; NULL == edge[i] && j < 2; ++j)
			{
				if (DoublyConnectedEdgeList::Primitive::edge == currentPrimitive->getPrimitiveType())
				{
#if 0
					switch (j)
					{
					case 0:
						startEdge = (DoublyConnectedEdgeList::Edge *)currentPrimitive;
						break;

					case 1:
						startEdge = (DoublyConnectedEdgeList::Edge *)currentPrimitive;
						startEdge = startEdge->twin;
						break;
					}
#endif
					startEdge = (DoublyConnectedEdgeList::Edge *)currentPrimitive;

					currentEdge = startEdge;

					do
					{
						if (currentEdge->origin == centre)
						{
							edge[i] = currentEdge;
							break;
						}

						currentEdge = currentEdge->next;
					} while (NULL != currentEdge && currentEdge != startEdge);

					if (NULL == edge[i])
					{
						startEdge = startEdge->twin;
					}

					//currentEdge = startEdge;

					//do
					//{
					//    if(currentEdge->origin == centre)
					//    {
					//        edge[i] = currentEdge;
					//        break;
					//    }

					//    currentEdge = currentEdge->next;
					//}
					//while(NULL != currentEdge && currentEdge != startEdge);

					//if(NULL == edge[i])
					//{
					//    std::cout << "error" << std::endl;
					//}

					//std::cout << "Edge" << std::endl;
				}
				else
				{
#if 0
					switch (j)
					{
					case 0:
						startEdge = ((DoublyConnectedEdgeList::Face *)currentPrimitive)->component;
						break;

					case 1:
						startEdge = ((DoublyConnectedEdgeList::Face *)currentPrimitive)->component;
						startEdge = startEdge->twin;
						break;
					}
#endif
					startEdge = ((DoublyConnectedEdgeList::Face *)currentPrimitive)->component;

					//std::cout << "Face" << std::endl;
				}

				if (NULL == edge[i])
				{
					currentEdge = startEdge;

					do
					{
						if (currentEdge->origin == centre)
						{
							edge[i] = currentEdge;
							break;
						}

						currentEdge = currentEdge->next;
					} while (NULL != currentEdge && currentEdge != startEdge);
				}
			}

			if (1 == i && NULL != edge[i])
			{
				if (NULL != edge[i]->previous)
				{
					edge[i] = edge[i]->previous->twin;
				}
				else
				{
					edge[i] = NULL;
				}
			}

			//unusedEdge[i] = false;
		}
	}

	if (NULL != edge[0] && NULL != edge[1])
	{
		this->transverseEdges(centre, edge[0], edge[1], edgeList);

		//if(unusedEdge[0] && edgeList.size() > 0)
		//{
		//    edgeList.pop_front();
		//}

		//if(unusedEdge[1] && edgeList.size() > 0)
		//{
		//    edgeList.pop_back();
		//}

		if (edgeList.size() > 2)
		{
			edgeList.pop_front();
			edgeList.pop_back();
		}
		else
		{
			edgeList.clear();
		}
	}
	else
	{
		//std::cout << edge[0] << std::endl;
		//std::cout << edge[1] << std::endl;

		//edge[0] = edge[0];
	}
}

void DoublyConnectedEdgeList::transverseEdges(
	const DoublyConnectedEdgeList::Vertex *centre,
	const DoublyConnectedEdgeList::Edge   *edge1,
	const DoublyConnectedEdgeList::Edge   *edge2,
	std::list<const DoublyConnectedEdgeList::Edge *> &edgeList) const
{
	const DoublyConnectedEdgeList::Edge *currentEdge;

	int count;


	currentEdge = edge1;
	count = 1;

	while (NULL != currentEdge->previous && currentEdge != edge2 && count < DoublyConnectedEdgeList::max_sides)
	{
		edgeList.push_back(currentEdge);

		currentEdge = currentEdge->previous->twin;
		++count;
	}

	edgeList.push_back(currentEdge);

	if (currentEdge != edge2)
	{
		edgeList.clear();
	}
}

DoublyConnectedEdgeList::Vertex* DoublyConnectedEdgeList::findVertex(const Vector3D &pos) const
{
	DoublyConnectedEdgeList::Vertex *v;
	DoublyConnectedEdgeList::VertexSet::const_iterator cvitr;

	DoublyConnectedEdgeList::Vertex target(pos);
	cvitr = this->vertices.find(&target);
	v = NULL;
	if (cvitr != this->vertices.end())
	{
		v = *cvitr;
	}

	return v;
}

DoublyConnectedEdgeList::Vertex* DoublyConnectedEdgeList::findNearestVertexWithoutSearchTree(const Vector3D &pos) const
{

	DoublyConnectedEdgeList::Vertex *v;
	double minDist, currentDist;

	DoublyConnectedEdgeList::Vertex* currentVertex;
	DoublyConnectedEdgeList::VertexSet::const_iterator currentVertexIterator, endVertex;//, startVertex;

	v = NULL;
	minDist = FLT_MAX;

	endVertex = this->vertices.end();

	for (currentVertexIterator = this->vertices.begin(); currentVertexIterator != endVertex; ++currentVertexIterator)
	{
		currentVertex = *currentVertexIterator;
		currentDist = Vector3D::distance(currentVertex->position, pos);

		if (currentDist < minDist)
		{
			v = currentVertex;
			minDist = currentDist;
		}
	}


	return v;
}

DoublyConnectedEdgeList::Vertex* DoublyConnectedEdgeList::findNearestVertex(const Vector3D &pos) const
{
	DoublyConnectedEdgeList::Vertex *v;

	double queryPt[3];
	int    resultID;
	double resultPt[3];


	v = NULL;

#ifdef _USE_DCEL_SEARCH_TREE

	if (NULL != this->vpTree)
	{
		// TODO: search the nearest point using the vp tree here

		memcpy(queryPt, pos.data(), 3 * sizeof(double));

		knnsearch(this->vpTree, queryPt, 1, resultPt, &resultID);

		v = this->vArray[resultID];
	}

#endif

	if (NULL == v)
	{
		//std::cerr << "NULL == v" << std::endl;
		v = this->findNearestVertexWithoutSearchTree(pos);
	}


	return v;
}

void DoublyConnectedEdgeList::findVertex(const Vector3D &centre, double radius, std::list<DoublyConnectedEdgeList::Vertex *> &results) const
{
	DoublyConnectedEdgeList::VertexSet::const_iterator viter;


	for (viter = this->vertices.begin(); viter != this->vertices.end(); ++viter)
	{
		if (Vector3D::distance((*viter)->position, centre) <= radius)
		{
			results.push_back(*viter);
		}
	}
}

DoublyConnectedEdgeList::Edge *DoublyConnectedEdgeList::findNearestEdge(const Vector3D &position, Vector3D &nearestPoint, DoublyConnectedEdgeList::Vertex * &nearestVertex, float epsilon) const
{
	DoublyConnectedEdgeList::Edge *nearestEdge;

	Vector3D t, v, w, w2, disp;

	double   dispt;
	double   minimumDist;

	DoublyConnectedEdgeList::Vertex *vertex;
	DoublyConnectedEdgeList::Edge   *currentEdge, *startEdge;

	DoublyConnectedEdgeList::EdgeSet::const_iterator eiter;


	nearestEdge = NULL;
	nearestVertex = this->findNearestVertex(position);

	if (nearestVertex)
	{
		nearestPoint = nearestVertex->position;

		minimumDist = (nearestPoint - position).length();

		if (minimumDist > epsilon)
		{
			vertex = nearestVertex;

			for (eiter = vertex->outgoingEdges.begin(); eiter != vertex->outgoingEdges.end(); ++eiter)
			{
				currentEdge = startEdge = *eiter;

				do
				{
					w = currentEdge->origin->position;
					w2 = currentEdge->twin->origin->position;
					t = (w2 - w).normalize();
					v = w - position;
					dispt = -v.dot(t);

					if (dispt > 0.0 && dispt < (w2 - w).length())
					{
						w += t * dispt;
						disp = w - position;

						if (disp.length() < minimumDist)
						{
							nearestPoint = w;
							nearestEdge = currentEdge;
							minimumDist = disp.length();
						}
					}

					currentEdge = currentEdge->next;
				} while (currentEdge && currentEdge->next != startEdge);
			}
		}
	}


	return nearestEdge;
}

DoublyConnectedEdgeList::Face *DoublyConnectedEdgeList::findNearestFace(const Vector3D &position, Vector3D &nearestPoint, Vertex * &nearestVertex, Edge * &nearestEdge, float epsilon) const
{
	DoublyConnectedEdgeList::Face   *nearestFace;

	DoublyConnectedEdgeList::Edge *e1, *e2, *e3;

	Vector3D v, w, n, disp;
	double   minimumDist;

	int i;


	nearestPoint = position;
	nearestVertex = NULL;
	nearestEdge = NULL;
	nearestFace = NULL;

	//if(this->isTriangularMesh)
	{
		nearestEdge = this->findNearestEdge(position, nearestPoint, nearestVertex, epsilon);

		if (nearestEdge)
		{
			v = nearestPoint - position;
			minimumDist = v.length();

			if (minimumDist > epsilon)
			{
				for (i = 0; i < 2; ++i)
				{
					e1 = (i ? nearestEdge->twin : nearestEdge);
					e2 = e1->next;
					if (NULL == e2) continue;
					e3 = e2->next;
					if (NULL == e3) continue;

					n = ((e2->origin->position - e1->origin->position).cross(e3->origin->position - e2->origin->position)).normalize();
					disp = n * n.dot(v);

					if (disp.length() < minimumDist)
					{
						w = position + disp;

						if ((e2->origin->position - e1->origin->position).cross(w - e2->origin->position).dot(n) > 0.0 &&
							(e3->origin->position - e2->origin->position).cross(w - e3->origin->position).dot(n) > 0.0 &&
							(e1->origin->position - e3->origin->position).cross(w - e1->origin->position).dot(n) > 0.0)
						{
							nearestPoint = w;
							nearestFace = e1->incidentFace;
							minimumDist = disp.length();
						}
					}
				}
			}
		}
	}


	return nearestFace;
}

DoublyConnectedEdgeList::Primitive *DoublyConnectedEdgeList::findNearestFace(const Vector3D &position, Vector3D &nearestPoint, Vector3D &normal, float epsilon) const
{
	Primitive *nearestPrimitive;
	Vertex    *nearestVertex;
	Edge      *nearestEdge;
	Face      *nearestFace;


	nearestFace = this->findNearestFace(position, nearestPoint, nearestVertex, nearestEdge, epsilon);

	if (NULL != nearestFace)
	{
		nearestPrimitive = nearestFace;
	}
	else if (NULL != nearestEdge)
	{
		nearestPrimitive = nearestEdge;
	}
	else if (NULL != nearestVertex)
	{
		nearestPrimitive = nearestVertex;
	}
	else
	{
		nearestPrimitive = NULL;
	}

	if (NULL != nearestPrimitive)
	{
		normal = nearestPrimitive->computeNormal();
	}


	return nearestPrimitive;
}

void DoublyConnectedEdgeList::findCommonVertices(const DoublyConnectedEdgeList::Face *face1, const DoublyConnectedEdgeList::Face *face2, DoublyConnectedEdgeList::VertexSet &vertices) const
{
#ifdef DCEL_DEBUG
	assert(face1 != face2);
	assert(vertices.size() == 0);
#endif

	DoublyConnectedEdgeList::Edge *startEdge1, *currentEdge1;
	DoublyConnectedEdgeList::Edge *startEdge2, *currentEdge2;

	startEdge1 = face1->component;
	startEdge2 = face2->component;

	currentEdge1 = startEdge1;

	do {

		currentEdge2 = startEdge2;

		do
		{
			if (currentEdge1->origin == currentEdge2->origin ||
				currentEdge1->origin == currentEdge2->twin->origin)
			{
				vertices.insert(currentEdge1->origin);
			}

			if (currentEdge1->twin->origin == currentEdge2->origin ||
				currentEdge1->twin->origin == currentEdge2->twin->origin)
			{
				vertices.insert(currentEdge1->twin->origin);
			}

			currentEdge2 = currentEdge2->next;

		} while (currentEdge2 != startEdge2);

		currentEdge1 = currentEdge1->next;

	} while (currentEdge1 != startEdge1);
}

DoublyConnectedEdgeList::Vertex *DoublyConnectedEdgeList::isCorner(const DoublyConnectedEdgeList::Face *face,
	//const float position[3],
	const Vector3D &position,
	float epsilon) const
{
	DoublyConnectedEdgeList::Vertex *corner, *candidate;
	DoublyConnectedEdgeList::Edge   *startEdge, *currentEdge;

	currentEdge = startEdge = face->component;

	corner = NULL;

	do
	{
		candidate = currentEdge->origin;
		if (Vector3D::distance(candidate->position, position) < epsilon)
		{
			corner = candidate;
		}

		currentEdge = currentEdge->next;

	} while (!corner && currentEdge != startEdge);

	return corner;
}

void DoublyConnectedEdgeList::updateVertexIDs()
{
	DoublyConnectedEdgeList::VertexSet::const_iterator viter;

	unsigned int i;


	if (!this->isVertexIDUpToDate)
	{
		i = 1;

		for (viter = this->vertices.begin(); viter != this->vertices.end(); ++viter)
		{
			(*viter)->id = i++;
		}

		this->isVertexIDUpToDate = true;
	}
}

size_t DoublyConnectedEdgeList::getVertexID(const DoublyConnectedEdgeList::Vertex *vertex) const
{
	if (this->isVertexIDUpToDate)
	{
		return vertex->id;
	}
	else
	{
		return 0;
	}
}

void DoublyConnectedEdgeList::updateFaceIDs()
{
	DoublyConnectedEdgeList::FaceSet::const_iterator fiter;

	unsigned int i;


	if (!this->isFaceIDUpToDate)
	{
		i = 1;

		for (fiter = this->faces.begin(); fiter != this->faces.end(); ++fiter)
		{
			(*fiter)->id = i++;
		}

		this->isFaceIDUpToDate = true;
	}
}

size_t DoublyConnectedEdgeList::getFaceID(const DoublyConnectedEdgeList::Face *face) const
{
	if (this->isFaceIDUpToDate)
	{
		return face->id;
	}
	else
	{
		return 0;
	}
}

void DoublyConnectedEdgeList::subdivide(Face *face, Vertex *vertex)
{
	DoublyConnectedEdgeList::Edge *startEdge, *currentEdge, *nextEdge;
	DoublyConnectedEdgeList::Edge *newEdge, *newTwinEdge;
	DoublyConnectedEdgeList::Edge *previousNewTwinEdge;
	DoublyConnectedEdgeList::Edge *firstNewEdge;

	DoublyConnectedEdgeList::Vertex *toVertex;

	double dist;

#ifdef DCEL_DEBUG
	assert(this->vertices.find(vertex) != this->vertices.end());
	assert(this->faces.find(face) != this->faces.end());
	std::cout << "the number of faces before subdivision: " << this->faces.size() << std::endl;
#endif

	previousNewTwinEdge = NULL;

	currentEdge = startEdge = face->component;

	this->removeFace(face, true);

	while ((nextEdge = currentEdge->next) != startEdge)
	{
		toVertex = currentEdge->twin->origin;

		dist = Vector3D::distance(toVertex->position, vertex->position);
		newEdge = this->addEdge(toVertex, vertex, dist);
		newTwinEdge = newEdge->twin;

		newTwinEdge->next = nextEdge;
		nextEdge->previous = newTwinEdge;

		newEdge->previous = currentEdge;
		currentEdge->next = newEdge;

		if (!previousNewTwinEdge)
		{
			firstNewEdge = newEdge;
			previousNewTwinEdge = newTwinEdge;
		}
		else
		{
			previousNewTwinEdge->previous = newEdge;
			newEdge->next = previousNewTwinEdge;

			this->addFace(newEdge);

			previousNewTwinEdge = newTwinEdge;
		}

		currentEdge = nextEdge;
	}

	toVertex = currentEdge->next->origin;

	newEdge = this->addEdge(toVertex, vertex, Vector3D::distance(toVertex->position, vertex->position));
	newTwinEdge = newEdge->twin;

	newTwinEdge->next = nextEdge;
	nextEdge->previous = newTwinEdge;

	newEdge->previous = currentEdge;
	currentEdge->next = newEdge;

	previousNewTwinEdge->previous = newEdge;
	newEdge->next = previousNewTwinEdge;

	this->addFace(newEdge);

	newTwinEdge->previous = firstNewEdge;
	firstNewEdge->next = newTwinEdge;

	this->addFace(firstNewEdge);

#ifdef DCEL_DEBUG
	std::cout << "the number of faces after subdivision: " << this->faces.size() << std::endl;
#endif
}

size_t DoublyConnectedEdgeList::getNumberOfVertices() const
{
	return this->vertices.size();
}

size_t DoublyConnectedEdgeList::getNumberOfEdges() const
{
	return this->edges.size();
}

size_t DoublyConnectedEdgeList::getNumberOfFaces() const
{
	return this->faces.size();
}

double DoublyConnectedEdgeList::breathFirstSearch(DoublyConnectedEdgeList::Vertex *fromVertex,
	DoublyConnectedEdgeList::Vertex *toVertex,
	SearchPathList &searchPath) const
{
	std::list<Vertex *> dirtyVertices;
	std::queue<DoublyConnectedEdgeList::Vertex*> verticesToBeVisited;

	bool found;

	DoublyConnectedEdgeList::VertexSet::const_iterator cvitr;
	DoublyConnectedEdgeList::Vertex * currentVertex;
	DoublyConnectedEdgeList::Vertex * neighboringVertex;

	DoublyConnectedEdgeList::EdgeSet::const_iterator ceitr;

	double pathDistance;

	pathDistance = FLT_MAX;

	// empty the serach path
#ifdef DCEL_DEBUG
	assert(searchPath.size() == 0);
#endif

	//searchPath.clear();

	if (!fromVertex || !toVertex)
	{
		return pathDistance;
	}


	// put the from vertex into a queue
	found = false;
	fromVertex->pathCost = 0.0f;
	verticesToBeVisited.push(fromVertex);
	dirtyVertices.push_back(fromVertex);

	// while the queue is not empty and the "to" vertex has not yet been found
	// pop(visit) the front vertex from the queue, set the flag on
	// add all its "unvisited" neighbors into the queue and
	// set the backtrace pointers of its neighbors to point to it
	while (!verticesToBeVisited.empty())
	{
		currentVertex = verticesToBeVisited.front();
		verticesToBeVisited.pop();
		currentVertex->visited = true;

		if (currentVertex == toVertex)
		{
			found = true;
			break;
		}

		for (ceitr = currentVertex->outgoingEdges.begin(); ceitr != currentVertex->outgoingEdges.end(); ++ceitr)
		{
			neighboringVertex = (*ceitr)->twin->origin;

			if (neighboringVertex->pathCost != 0.0f)
			{
				neighboringVertex->pathCost = 0.0f;
				neighboringVertex->backtrace = currentVertex;
				verticesToBeVisited.push(neighboringVertex);
				dirtyVertices.push_back(neighboringVertex);
			}
		}
	}

	if (found)
	{
		pathDistance = 0.0;

		currentVertex = toVertex;

		while (currentVertex != fromVertex)
		{
			searchPath.push_front(SearchPathNode(currentVertex));

			pathDistance += Vector3D::distance(currentVertex->position, currentVertex->backtrace->position);

			currentVertex = currentVertex->backtrace;
		}

		searchPath.push_front(SearchPathNode(currentVertex));
	}

	std::for_each(dirtyVertices.begin(), dirtyVertices.end(), DoublyConnectedEdgeList::VertexReset());

	return pathDistance;
}

double DoublyConnectedEdgeList::breathFirstSearch(const Vector3D &from, const Vector3D &to, SearchPathList &searchPath) const
{
	DoublyConnectedEdgeList::Vertex * fromVertex, *toVertex;

	fromVertex = this->findNearestVertex(from);
	toVertex = this->findNearestVertex(to);

	return breathFirstSearch(fromVertex, toVertex, searchPath);
}

double DoublyConnectedEdgeList::computeShortestPath(Vertex *fromVertex, Vertex *toVertex, SearchPathList &searchPath, long int maximumNumberOfVertex) const
{
	std::list<Vertex *> dirtyVertices;
	std::set<DoublyConnectedEdgeList::Vertex*, DoublyConnectedEdgeList::CostComparison> verticesToBeVisited;

	bool found;

	double newCost, pathCost;

	DoublyConnectedEdgeList::Vertex * currentVertex;
	DoublyConnectedEdgeList::Vertex * neighboringVertex;

	DoublyConnectedEdgeList::Edge * currentEdge;

	DoublyConnectedEdgeList::EdgeSet::const_iterator currentEdgeIterator, endEdgeIterator, startEdgeIterator;

	int i;


	// empty the serach path
	searchPath.clear();

	// find the "from" and "to" vertices
	if (!fromVertex || !toVertex)
	{
		std::cerr << "Unknown end points -- computeShortestPath cannot process!" << std::endl;
		return FLT_MAX;
	}

	if (fromVertex == toVertex)
	{
		//std::cerr << "Unknown end points -- computeShortestPath cannot process!" << std::endl;

		searchPath.push_back(DoublyConnectedEdgeList::SearchPathNode(fromVertex));

		return 0.0f;
	}

	// put the from vertex into a queue
	found = false;
	i = 0;
	fromVertex->pathCost = 0.0f;
	verticesToBeVisited.insert(fromVertex);
	dirtyVertices.push_back(fromVertex);

	// while the queue is not empty and the "to" vertex has not yet been found
	// pop(visit) the front vertex from the queue, set the flag on
	// add all its "unvisited" neighbors into the queue and
	// set the backtrace pointers of its neighbors to point to it
	while (!verticesToBeVisited.empty())
	{
		currentVertex = (*verticesToBeVisited.begin());
		verticesToBeVisited.erase(verticesToBeVisited.begin());

		if (currentVertex->visited)
			continue;

		currentVertex->visited = true;
		++i;

		//cout << "current path cost: " << currentVertex->pathCost << endl;
		if (currentVertex == toVertex)
		{
			found = true;
			break;
		}

		if (i == maximumNumberOfVertex)
		{
			break;
		}

		startEdgeIterator = currentVertex->outgoingEdges.begin();
		endEdgeIterator = currentVertex->outgoingEdges.end();

		for (currentEdgeIterator = startEdgeIterator; currentEdgeIterator != endEdgeIterator; ++currentEdgeIterator)
		{
			currentEdge = *currentEdgeIterator;
			neighboringVertex = currentEdge->twin->origin;

			if (!neighboringVertex->visited)
			{
				newCost = currentVertex->pathCost + currentEdge->cost;

				if (neighboringVertex->pathCost == FLT_MAX)
				{
					neighboringVertex->backtrace = currentVertex;
					neighboringVertex->pathCost = newCost;
					verticesToBeVisited.insert(neighboringVertex);
					dirtyVertices.push_back(neighboringVertex);
				}
				else if (newCost < neighboringVertex->pathCost)
				{
					verticesToBeVisited.erase(neighboringVertex);
					neighboringVertex->backtrace = currentVertex;
					neighboringVertex->pathCost = newCost;
					verticesToBeVisited.insert(neighboringVertex);
				}
			}
		}
	}

	if (found)
	{
		currentVertex = toVertex;

		while (currentVertex != fromVertex)
		{
			searchPath.push_front(DoublyConnectedEdgeList::SearchPathNode(currentVertex));
			currentVertex = currentVertex->backtrace;
		}

		searchPath.push_front(DoublyConnectedEdgeList::SearchPathNode(currentVertex));

		if (searchPath.size() == 0 || toVertex->pathCost == FLT_MAX)
		{
			std::cerr << "CRITICAL ERROR IN computeShortestPath! path with size 0 or distance = FLT_MAX" << std::endl;
		}
	}
	else
	{
		std::cerr << "CRITICAL ERROR IN computeShortestPath! path not found!" << std::endl;
	}

	pathCost = toVertex->pathCost;

	std::for_each(dirtyVertices.begin(), dirtyVertices.end(), DoublyConnectedEdgeList::VertexReset());

	return pathCost;
}

double DoublyConnectedEdgeList::computeShortestPath(const Vector3D &from, const Vector3D &to, SearchPathList &searchPath, long int maximumNumberOfVertex) const
{
	DoublyConnectedEdgeList::Vertex * fromVertex, *toVertex;

	fromVertex = this->findNearestVertex(from);
	toVertex = this->findNearestVertex(to);

	return this->computeShortestPath(fromVertex, toVertex, searchPath, maximumNumberOfVertex);
}

Vector3D DoublyConnectedEdgeList::Face::computeNormal() const
{
	DoublyConnectedEdgeList::Edge *startEdge, *currentEdge, *lastEdge;
	Vector3D v1, v2, normal;


	startEdge = this->component;
	currentEdge = startEdge;
	lastEdge = startEdge->previous;

	v1 = lastEdge->twin->origin->position - lastEdge->origin->position;

	do
	{
		v2 = currentEdge->twin->origin->position - currentEdge->origin->position;

		normal += v1.cross(v2);

		v1 = v2;
	} while (currentEdge != startEdge);

	return normal;
}

void DoublyConnectedEdgeList::Face::render(float epsilon) const
{
	GLdouble projection[16];

	DoublyConnectedEdgeList::Edge *currentEdge, *startEdge;


	if (0.0 != epsilon)
	{
		glGetDoublev(GL_PROJECTION_MATRIX, projection);

		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		glTranslated(0.0, 0.0, -epsilon);
		glMultMatrixd(projection);
		glMatrixMode(GL_MODELVIEW);
	}

	currentEdge = startEdge = this->component;

	glBegin(GL_POLYGON);

	do
	{
		glNormal3dv(currentEdge->origin->computeNormal().normalize().data());
		glVertex3dv(currentEdge->origin->position.data());

		currentEdge = currentEdge->next;
	} while (currentEdge != startEdge);

	glEnd();

	if (0.0 != epsilon)
	{
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		glMatrixMode(GL_MODELVIEW);
	}
}
