#define _USE_MATH_DEFINES

#include "patchlist.h"

#include "resampling.h"
#include "geodesic.h"
#include "patch.h"
#include "OBJModel.h"

#include <cmath>

#include <windows.h>
#include <GL/gl.h>


PatchList::PatchList(const DoublyConnectedEdgeList &newdcel) : dcel(newdcel), patchConnection(), patches(), controlPoints(), resolution(1)
{
}

PatchList::~PatchList()
{
	this->clear();
}

void PatchList::clear()
{
	FaceToPatchMap::iterator piter;


	this->controlPoints.clear();

	for (piter = this->patches.begin(); piter != this->patches.end(); ++piter)
	{
		delete piter->second;
	}

	this->patches.clear();
	this->patchConnection.clear();
}

bool PatchList::loadFile(const char filename[])
{
	DoublyConnectedEdgeList::FaceSet::const_iterator fiter;
	const DoublyConnectedEdgeList::Edge *startEdge, *currentEdge;

	Vector3D controlPoints[4], controlNormal[4];

	int count;

	bool ret;


	ret = true;

	this->clear();

	try
	{
		if (!this->patchConnection.loadFile(filename))
		{
			throw "Fail to load the file!";
		}

		if (!this->patchConnection.validate())
		{
			throw "Patch connection error!";
		}

		for (fiter = this->patchConnection.getFaceSet().begin(); fiter != this->patchConnection.getFaceSet().end(); ++fiter)
		{
			startEdge = (*fiter)->getComponent();

			currentEdge = startEdge;
			count = 0;

			do
			{
				controlPoints[count] = currentEdge->getOrigin()->getPosition();
				controlNormal[count] = currentEdge->getOrigin()->computeNormal().normalize();

				currentEdge = currentEdge->getNextEdge();
				++count;
			} while (currentEdge != startEdge && count < 4);

			if (count < 3 || currentEdge != startEdge)
			{
				throw "Invalid patch file format!";
			}

			//if(count == 4)
			{
				this->patches[*fiter] = new QuadPatch(controlPoints, controlNormal, this->resolution);
			}
		}
	}
	catch (const char err[])
	{
		std::cerr << err << std::endl;

		this->clear();

		ret = false;
	}


	return ret;
}

void PatchList::saveFile(const char filename[]) const
{
	this->patchConnection.saveFile(filename);
}

const DoublyConnectedEdgeList &PatchList::getPatchConnection() const
{
	return this->patchConnection;
}

const PatchList::FaceToPatchMap &PatchList::getFaceToPatchMap() const
{
	return this->patches;
}

const DoublyConnectedEdgeList &PatchList::getOriginalSurface() const
{
	return this->dcel;
}

void PatchList::addPatch(Patch *patch/*, float epsilon*/)
{
	//DoublyConnectedEdgeList::Vertex *v[4];
	//DoublyConnectedEdgeList::Edge   *e[4];
	DoublyConnectedEdgeList::Face   *f;

	std::list<Vector3D> p;

	unsigned int i;


#if 0
	for (i = 0; i < patch->numberOfControlPoints; ++i)
	{
		//v[i] = this->patchConnection.findNearestVertex(patch->controlPoints[i]);

		//if(NULL == v[i] || (v[i]->getPosition() - patch->controlPoints[i]).length() > epsilon)
		//{
		v[i] = this->patchConnection.addVertex(patch->controlPoints[i]);
		//}
		//else
		//{
		//    patch->controlPoints[i] = v[i]->getPosition();

		//    std::cout << "used duplicated points!" << std::endl;
		//}
	}

	for (i = 0; i < patch->numberOfControlPoints; ++i)
	{
		e[i] = this->patchConnection.addEdge(
			v[i], v[(i + 1) % patch->numberOfControlPoints],
			Vector3D::distance(v[i]->getPosition(), v[(i + 1) % patch->numberOfControlPoints]->getPosition()));
	}

	for (i = 0; i < patch->numberOfControlPoints; ++i)
	{
		e[i]->setNextEdge(e[(i + 1) % patch->numberOfControlPoints]);
	}

	f = this->patchConnection.addFace(e[0]);

	//if(NULL == f)
	//{
	//    // reverse orientation

	//    for(i = 0; i < patch->numberOfControlPoints; ++i)
	//    {
	//        e[i]->setNextEdge(NULL);
	//        e[i]->getTwinEdge()->setNextEdge(e[(i + 1) % patch->numberOfControlPoints]->getTwinEdge());
	//    }

	//    f = this->patchConnection.addFace(e[0]->getTwinEdge());
	//}

#endif

	for (i = 0; i < patch->numberOfControlPoints; ++i)
	{
		p.push_back(patch->controlPoints[i]);
	}

	f = this->patchConnection.addFace(p);

	this->patches[f] = patch;

	//for(i = 0; i < patch->numberOfControlPoints; ++i)
	//{
	//    this->syncMeshPoint(dcel, v[i]);
	//}

	patch->setResolution(this->resolution);
}

int PatchList::getPatchResolution() const
{
	return this->resolution;
}

void PatchList::setPatchResolution(int newResolution)
{
	FaceToPatchMap::const_iterator piter;


	if (this->resolution != newResolution)
	{
		for (piter = this->patches.begin(); piter != this->patches.end(); ++piter)
		{
			piter->second->setResolution(newResolution);
		}

		this->resolution = newResolution;
	}
}

void PatchList::subdividePatchList(PatchList &subdividedPatchList) const
{
	FaceToPatchMap::const_iterator piter;

	const Vector3D *originalControlPoints;
	const Vector3D *originalControlNormal;

	Vector3D p[4], n[4];

	double t[2], s[2];

	Patch *patch;

	unsigned int y[2], x[2], i, j;


	try
	{
		if (&subdividedPatchList == this)
		{
			throw "&subdividedPatchList == this";
		}

		if (&subdividedPatchList.dcel != &this->dcel)
		{
			throw "&subdividedPatchList.dcel != &this->dcel";
		}

		subdividedPatchList.clear();
		subdividedPatchList.setPatchResolution((this->resolution + 1) / 2);

		for (piter = this->patches.begin(); piter != this->patches.end(); ++piter)
		{
			switch (piter->second->getNumberOfControlPoints())
			{
			case 3:

				break;

			case 4:

				originalControlPoints = piter->second->getControlPoints();
				originalControlNormal = piter->second->getControlNormal();

				y[1] = 0;
				t[1] = 0.0;

				while (y[1] < 2)
				{
					y[0] = y[1];
					t[0] = t[1];

					++y[1];
					t[1] = (double)y[1] / 2.0;

					x[1] = 0;
					s[1] = 0.0;

					while (x[1] < 2)
					{
						x[0] = x[1];
						s[0] = s[1];

						++x[1];
						s[1] = (double)x[1] / 2.0;

						for (j = 0; j < 2; ++j)
						{
							for (i = 0; i < 2; ++i)
							{
								p[2 * j + i] =
									originalControlPoints[0] * ((1.0 - s[i]) * (1.0 - t[j])) +
									originalControlPoints[1] * ((s[i]) * (1.0 - t[j])) +
									originalControlPoints[2] * ((s[i]) * (t[j])) +
									originalControlPoints[3] * ((1.0 - s[i]) * (t[j]));
							}
						}

						p[3].swap(p[2]);

						for (j = 0; j < 2; ++j)
						{
							for (i = 0; i < 2; ++i)
							{
								n[2 * j + i] =
									originalControlNormal[0] * ((1.0 - s[i]) * (1.0 - t[j])) +
									originalControlNormal[1] * ((s[i]) * (1.0 - t[j])) +
									originalControlNormal[2] * ((s[i]) * (t[j])) +
									originalControlNormal[3] * ((1.0 - s[i]) * (t[j]));
							}
						}

						n[3].swap(n[2]);

						//std::cout << p[0] << p[1] << p[2] << p[3] << std::endl;

						patch = new QuadPatch(p, n, subdividedPatchList.getPatchResolution());

						if (NULL == patch)
						{
							throw "NULL == patch";
						}

						patch->resampling(this->dcel);

						subdividedPatchList.addPatch(patch);
					}
				}

				break;

			default:

				break;
			}
		}

		if (!subdividedPatchList.patchConnection.validate())
		{
			throw "!subdividedPatchList.patchConnection.validate()";
		}
	}
	catch (const char err[])
	{
		subdividedPatchList.clear();

		std::cerr << err << std::endl;
	}
}

void PatchList::resampling()
{
	FaceToPatchMap::const_iterator piter;


	for (piter = this->patches.begin(); piter != this->patches.end(); ++piter)
	{
		piter->second->resampling(dcel);
	}
}

void PatchList::syncMeshPoint(const DoublyConnectedEdgeList::Vertex *v1)
{
	DoublyConnectedEdgeList::Face   *f;

	Patch *patch;

	unsigned int count;

	Vector3D controlPoints[4], controlNormal[4];

	const DoublyConnectedEdgeList::Edge *startEdge, *currentEdge;

	DoublyConnectedEdgeList::VertexSet::const_iterator viter;
	DoublyConnectedEdgeList::EdgeSet::const_iterator eiter;


	for (eiter = v1->getOutgoingEdges().begin(); eiter != v1->getOutgoingEdges().end(); ++eiter)
	{
		f = (*eiter)->getIncidentFace();

		if (f)
		{
			patch = patches[f];

			count = 0;
			currentEdge = startEdge = f->getComponent();

			do
			{
				controlPoints[count] = currentEdge->getOrigin()->getPosition();
				controlNormal[count] = currentEdge->getOrigin()->computeNormal().normalize();

				currentEdge = currentEdge->getNextEdge();
				++count;
			} while (currentEdge != startEdge);

			patch->setControlPoints(controlPoints, controlNormal);
			patch->resampling(dcel);
		}
	}
}

void PatchList::syncMeshPoints()
{
	Vector3D controlPoints[4], controlNormal[4];

	DoublyConnectedEdgeList::VertexSet::const_iterator viter;
	DoublyConnectedEdgeList::EdgeSet::const_iterator eiter;


	for (viter = this->patchConnection.getVertexSet().begin(); viter != this->patchConnection.getVertexSet().end(); ++viter)
	{
		this->syncMeshPoint(*viter);
	}
}

void PatchList::outputMesh(OBJModel &model) const
{
	Patch *patch;

	PolygonElement newPolygon;

	const PolygonElement *currentPolygon;

	size_t voffset, noffset, toffset;

	FaceToPatchMap::const_iterator piter;

	size_t i, j;


	for (piter = this->patches.begin(); piter != this->patches.end(); ++piter)
	{
		patch = piter->second;

		voffset = model.getNumberOfVertices();
		noffset = model.getNumberOfNormal();
		toffset = model.getNumberOfTextureCoordinates();

		for (i = 0; i < patch->mesh->getNumberOfVertices(); ++i)
		{
			model.addVertex(patch->mesh->getVertex(i));
		}

		for (i = 0; i < patch->mesh->getNumberOfNormal(); ++i)
		{
			model.addNormal(patch->mesh->getNormal(i));
		}

		for (i = 0; i < patch->mesh->getNumberOfTextureCoordinates(); ++i)
		{
			model.addTextureCoordinate(patch->mesh->getTextureCoordinate(i));
		}

		for (i = 0; i < patch->mesh->getNumberOfPolygons(); ++i)
		{
			newPolygon.vertexIndices.clear();
			newPolygon.normalIndices.clear();
			newPolygon.textureCoordinateIndices.clear();

			currentPolygon = patch->mesh->getPolygon(i);

			for (j = 0; j < currentPolygon->vertexIndices.size(); ++j)
			{
				newPolygon.vertexIndices.push_back(currentPolygon->vertexIndices[j] + voffset);

				if (!currentPolygon->normalIndices.empty())
				{
					newPolygon.normalIndices.push_back(currentPolygon->normalIndices[j] + noffset);
				}

				if (!currentPolygon->textureCoordinateIndices.empty())
				{
					newPolygon.textureCoordinateIndices.push_back(currentPolygon->textureCoordinateIndices[j] + toffset);
				}
			}

			model.addPolygon(newPolygon);
		}
	}
}

void PatchList::moveControlPoint(DoublyConnectedEdgeList::Vertex * v1, const Vector3D &to, float epsilon)
{
	DoublyConnectedEdgeList::Vertex *v2;

	Vector3D controlPoints[4], controlNormal[4];

	DoublyConnectedEdgeList::EdgeSet::const_iterator eiter;


	if (NULL != v1)
	{
		v2 = this->patchConnection.findNearestVertex(to);

		if (NULL != v2 && v1 != v2 && (v2->getPosition() - to).length() < epsilon)
		{
			//std::cout
			//    << "merge the point\n" << v2->getPosition() << '\n'
			//    << "with\n"            << v1->getPosition() << std::endl;

			std::cout << "MERGED two points\n" << std::endl;

			this->patchConnection.mergeVertex(v2, v1);
		}

		//std::cout
		//    << "move the point from\n" << v1->getPosition() << '\n'
		//    << "to\n"                  << to                << std::endl;

		std::cout << "MOVED a point" << std::endl;

		this->patchConnection.moveVertex(v1, to);

#if 0
		for (eiter = v1->getOutgoingEdges().begin(); eiter != v1->getOutgoingEdges().end(); ++eiter)
		{
			f = (*eiter)->getIncidentFace();

			if (f)
			{
				patch = patches[f];

				count = 0;
				currentEdge = startEdge = f->getComponent();

				do
				{
					controlPoints[count] = currentEdge->getOrigin()->getPosition();
					controlNormal[count] = currentEdge->getOrigin()->computeNormal().normalize();

					currentEdge = currentEdge->getNextEdge();
					++count;
				} while (currentEdge != startEdge);

				patch->setControlPoints(controlPoints, controlNormal);
				patch->resampling(dcel);
			}
		}
#endif

		this->syncMeshPoint(v1);
	}
}

void PatchList::upateIDs()
{
	this->patchConnection.updateVertexIDs();
	this->patchConnection.updateFaceIDs();
}

void PatchList::removePatch(DoublyConnectedEdgeList::Face *face)
{
	Patch *patch;


	patch = this->patches[face];

	if (NULL != patch)
	{
		this->patches.erase(face);
		delete patch;
		this->patchConnection.removeFace(face, false);
	}
}

void PatchList::projectParameterization(OBJModel &initialParameterization, unsigned int gridDensity)
{
	DoublyConnectedEdgeList::Face *face;
	const DoublyConnectedEdgeList::Edge *currentEdge, *startEdge;

	double length;

	Patch *patch;

	int count;

	Vector3D controlPoints[4], controlNormal[4];

	DoublyConnectedEdgeList::FaceSet::const_iterator iter;


	this->clear();

	Resampling(dcel).transform(initialParameterization);

	this->patchConnection.convertFromOBJModel(initialParameterization);
	this->patchConnection.validate();

	for (iter = this->patchConnection.getFaceSet().begin(); iter != this->patchConnection.getFaceSet().end(); ++iter)
	{
		face = *iter;

		count = 0;
		currentEdge = startEdge = face->getComponent();

		do
		{
			controlPoints[count] = currentEdge->getOrigin()->getPosition();
			controlNormal[count] = currentEdge->getOrigin()->computeNormal();

			length = controlNormal[count].length();

			if (length > Vector3D::epsilon)
			{
				controlNormal[count] /= length;
			}
			else
			{
				controlNormal[count].setZero();
			}

			currentEdge = currentEdge->getNextEdge();
			++count;
		} while (count < 4 && currentEdge != startEdge);

		if (currentEdge == startEdge)
		{
			this->patches.erase(face);

			switch (count)
			{
			case 3:
				patch = new TriPatch(controlPoints, controlNormal, gridDensity);
				break;
			case 4:
				patch = new QuadPatch(controlPoints, controlNormal, gridDensity);
				break;
			}

			patch->resampling(dcel);

			this->patches[face] = patch;
		}
	}
}

void PatchList::updatePatchConnectionVertices(const std::list<Vector3D> &newPositions)
{
	DoublyConnectedEdgeList::VertexSet originalVertexSet;
	DoublyConnectedEdgeList::VertexSet::const_iterator viter;

	std::list<Vector3D>::const_iterator piter;


	if (this->patchConnection.getVertexSet().size() == newPositions.size())
	{
		originalVertexSet = this->patchConnection.getVertexSet();

		piter = newPositions.begin();

		for (viter = originalVertexSet.begin(); viter != originalVertexSet.end(); ++viter)
		{
			this->patchConnection.moveVertex(*viter, *piter);

			++piter;
		}
	}
	else
	{
		std::cout << "Cannot update vertices on the patch list DCEL" << std::endl;
	}
}

#if (PATCHLIST_MAJOR_VERSION == 1)

void PatchList::minimizeStretch(const DoublyConnectedEdgeList &dcel)
{
	std::list<Geodesic *> geodesics;
	std::list<Vector3D> newPositions;

	Geodesic *geodesic;

	Vector3D start, direction, newPosition, sum;

	Vector3D dir1, dir2, normal;

	DoublyConnectedEdgeList::Vertex *v1, *v2, *v3;

	DoublyConnectedEdgeList::Edge *currentEdge, *startEdge;

	DoublyConnectedEdgeList::VertexSet::const_iterator viter;
	DoublyConnectedEdgeList::EdgeSet::const_iterator eiter;

	std::list<Geodesic *>::iterator giter, giter2;

	unsigned int minSide, maxSide;

	unsigned int count;


	this->patchConnection.getPolygonSide(minSide, maxSide);

	if (4 == minSide && 4 == maxSide)
	{
		count = 0;

		for (viter = this->patchConnection.getVertexSet().begin(); viter != this->patchConnection.getVertexSet().end(); ++viter)
		{
			sum.setZero();

			v1 = *viter;

			try
			{
				currentEdge = startEdge = *v1->getOutgoingEdges().begin();

				do
				{
					v2 = currentEdge->getNextEdge()->getOrigin();

					geodesic = Geodesic::compute(dcel, v1->getPosition(), v2->getPosition());

					if (NULL == geodesic)
					{
						throw "Fail to find a geodesic!";
					}

					sum += geodesic->forwardDirection() * geodesic->length();

					geodesics.push_back(geodesic);

					v3 = currentEdge->getNextEdge()->getNextEdge()->getOrigin();

					geodesic = Geodesic::compute(dcel, v1->getPosition(), v3->getPosition());

					if (NULL == geodesic)
					{
						throw "Fail to find a geodesic!";
					}

					sum += geodesic->forwardDirection() * (geodesic->length() / M_SQRT2);

					geodesics.push_back(geodesic);

					currentEdge = currentEdge->getPreviousEdge()->getTwinEdge();
				} while (currentEdge != startEdge);

				if (v1->getOutgoingEdges().empty())
				{
					throw "No edge!";
				}

				sum /= (4.0 * (double)v1->getOutgoingEdges().size());

				giter = geodesics.end();
				--giter;

				giter2 = geodesics.begin();

				while (giter2 != geodesics.end())
				{
					dir1 = (*giter)->forwardDirection();
					dir2 = (*giter2)->forwardDirection();

					normal = dir1.cross(dir2);

					if (dir1.cross(sum).dot(normal) >= 0.0 &&
						dir2.cross(sum).dot(normal) <= 0.0)
					{
						break;
					}

					//ori1 = dir1.cross(sum );
					//ori2 = sum .cross(dir2);

					//if( ori1.length() <= FLT_EPSILON ||
					//    ori2.length() <= FLT_EPSILON ||
					//    ori1.dot(ori2) > 0.0)
					//{
					//    break;
					//}

					giter = giter2;
					++giter2;
				}

				if (giter2 == geodesics.end())
				{
					newPositions.push_back(v1->getPosition());
					throw "Error 4!";
				}

				if (sum.length() <= FLT_EPSILON)
				{
					newPosition = v1->getPosition();
				}
				else
				{
					newPosition = Geodesic::findRelativePositionOnSurface(
						dcel,
						v1->getPosition(), (*giter)->rbegin()->position, (*giter2)->rbegin()->position,
						sum);
				}


				newPositions.push_back(newPosition);
			}
			catch (const char err[])
			{
				std::cerr << err << std::endl;
			}

			for (giter = geodesics.begin(); giter != geodesics.end(); ++giter)
			{
				delete *giter;
			}
			geodesics.clear();

			if (++count % 20 == 0)
			{
				printf(".");
				fflush(stdout);
			}
		}

		this->updatePatchConnectionVertices(newPositions);

		printf("\n");
	}
	else
	{
		std::cout << "The stretch minization operator only applies to quad-based parameterization" << std::endl;
	}
}

#if 0

void PatchList::minimizeAreaDistortion(const DoublyConnectedEdgeList &dcel)
{
	DoublyConnectedEdgeList::Vertex *v1, *v2, *v3;

	Geodesic *g1, *g2;
	Vector3D start, direction, newPosition, sum, dir1, dir2;

	double areaSum, area;

	DoublyConnectedEdgeList::VertexSet originalVertexSet;
	std::list<Vector3D> newPositions;

	DoublyConnectedEdgeList::VertexSet::const_iterator viter;
	DoublyConnectedEdgeList::EdgeSet::const_iterator eiter;

	std::list<Vector3D>::const_iterator piter;


	for (viter = this->patchConnection.getVertexSet().begin(); viter != this->patchConnection.getVertexSet().end(); ++viter)
	{
		areaSum = 0.0;
		sum.setZero();

		v1 = *viter;

		for (eiter = v1->getOutgoingEdges().begin(); eiter != v1->getOutgoingEdges().end(); ++eiter)
		{
			v2 = (*eiter)->getNextEdge()->getOrigin();
			v3 = (*eiter)->getPreviousEdge()->getOrigin();

			g1 = Geodesic::compute(dcel, v1->getPosition(), v3->getPosition());
			g2 = Geodesic::compute(dcel, v1->getPosition(), v2->getPosition());

			if (NULL != g1 && NULL != g2)
			{
				dir1 = g1->forwardDirection() * g1->length();
				dir2 = g2->forwardDirection() * g2->length();

				area = dir1.cross(dir2).length();

				areaSum += area;
				sum += (dir1 + dir2) * (area / 3.0);

				delete g1;
				delete g2;
			}
			else
			{
				if (NULL != g1)
				{
					delete g1;
				}

				if (NULL != g2)
				{
					delete g2;
				}

				throw "Fail to find a geodesic!";
			}
		}


		if (!v1->getOutgoingEdges().empty())
		{
			sum /= areaSum;

			start = v1->getPosition() - sum / 5.0;
			direction = -v1->computeNormal().normalize();

			dcel.rayTrace(start - direction * 1e-1f, direction, newPosition, false);

			//this->patchConnection.moveVertex(v1, newPosition);
			newPositions.push_back(newPosition);
		}
		else
		{
			newPositions.push_back(v1->getPosition());
		}
	}

	this->updatePatchConnectionVertices(newPositions);
}

#endif

void PatchList::minimizeAreaDistortionAccurately(const DoublyConnectedEdgeList &dcel)
{
	std::list<Geodesic *> geodesics;
	std::list<Vector3D> newPositions;

	Geodesic *g1, *g2;

	Vector3D start, direction, newPosition;

	Vector3D dir1, dir2, normal;

	Vector3D centroidSum, centroid;

	double areaSum, area;

	Patch *patch;

	DoublyConnectedEdgeList::Face   *face;
	DoublyConnectedEdgeList::Vertex *v1, *v2, *v3;
	DoublyConnectedEdgeList::Edge   *currentEdge, *nextEdge, *startEdge;

	DoublyConnectedEdgeList::VertexSet::const_iterator viter;
	DoublyConnectedEdgeList::EdgeSet::const_iterator eiter;

	std::list<Geodesic *>::iterator giter, giter2;

	unsigned int minSide, maxSide;

	unsigned int count;


	this->patchConnection.getPolygonSide(minSide, maxSide);

	if (3 <= minSide && 4 >= maxSide)
	{
		count = 0;

		for (viter = this->patchConnection.getVertexSet().begin(); viter != this->patchConnection.getVertexSet().end(); ++viter)
		{
			centroidSum.setZero();

			areaSum = 0.0;

			v1 = *viter;

			try
			{
				currentEdge = startEdge = *v1->getOutgoingEdges().begin();

				do
				{
					g1 = Geodesic::compute(dcel, v1->getPosition(), currentEdge->getTwinEdge()->getOrigin()->getPosition());

					if (NULL == g1)
					{
						throw "Fail to find a geodesic!";
					}

					geodesics.push_back(g1);

					face = currentEdge->getIncidentFace();

					if (NULL == face)
					{
						throw "Null face";
					}

					if (NULL == patch)
					{
						throw "Null patch";
					}

					patch = this->patches[face];

					patch->getCentroid(area, centroid);

					areaSum += area;
					centroidSum += centroid * area;

					if (NULL == currentEdge->getPreviousEdge())
					{
						currentEdge = NULL;
					}
					else
					{
						currentEdge = currentEdge->getPreviousEdge()->getTwinEdge();
					}
				} while (currentEdge && currentEdge != startEdge);

				if (0.0 == areaSum)
				{
					throw "Zero Area!";
				}

				if (currentEdge != startEdge)
				{
					throw "Open manifold!";
				}

				centroidSum /= areaSum;
				centroidSum -= v1->getPosition();

				centroidSum /= 2.0;

				giter = geodesics.end();
				--giter;

				giter2 = geodesics.begin();

				while (giter2 != geodesics.end())
				{
					dir1 = (*giter)->forwardDirection();
					dir2 = (*giter2)->forwardDirection();

					normal = dir1.cross(dir2);

					if (dir1.cross(centroidSum).dot(normal) >= 0.0 &&
						dir2.cross(centroidSum).dot(normal) <= 0.0)
					{
						break;
					}

					//ori1 = dir1.cross(sum );
					//ori2 = sum .cross(dir2);

					//if( ori1.length() <= FLT_EPSILON ||
					//    ori2.length() <= FLT_EPSILON ||
					//    ori1.dot(ori2) > 0.0)
					//{
					//    break;
					//}

					giter = giter2;
					++giter2;
				}

				if (giter2 == geodesics.end())
				{
					newPositions.push_back(v1->getPosition());
					throw "Error 4!";
				}

				if (centroidSum.length() <= FLT_EPSILON)
				{
					newPosition = v1->getPosition();
				}
				else
				{
					newPosition = Geodesic::findRelativePositionOnSurface(
						dcel,
						v1->getPosition(), (*giter)->rbegin()->position, (*giter2)->rbegin()->position,
						centroidSum);
				}

				newPositions.push_back(newPosition);
			}
			catch (const char err[])
			{
				std::cerr << err << std::endl;
			}

			for (giter = geodesics.begin(); giter != geodesics.end(); ++giter)
			{
				delete *giter;
			}
			geodesics.clear();

			if (++count % 20 == 0)
			{
				printf(".");
				fflush(stdout);
			}
		}

		this->updatePatchConnectionVertices(newPositions);

		this->syncMeshPoints(dcel);

		printf("\n");
	}
	else
	{
		std::cout << "The stretch minization operator only applies to quad-based parameterization" << std::endl;
	}
}

void PatchList::minimizeAreaDistortion(const DoublyConnectedEdgeList &dcel)
{
	std::list<Geodesic *> geodesics;
	std::list<Vector3D> newPositions;

	std::vector<double> rs, thetas;

	Geodesic *g1, *g2;

	Vector3D start, direction, newPosition;//, sum;

	Vector3D dir1, dir2, dir3, normal;

	double angleSum, angleSum2, angle, angle2, areaSum, area;
	double Cr, Ctheta, Cx, Cy;

	DoublyConnectedEdgeList::Vertex *v1, *v2, *v3;

	DoublyConnectedEdgeList::Edge *currentEdge, *nextEdge, *startEdge;

	DoublyConnectedEdgeList::VertexSet::const_iterator viter;
	DoublyConnectedEdgeList::EdgeSet::const_iterator eiter;

	std::list<Geodesic *>::iterator giter, giter2;

	unsigned int minSide, maxSide;

	unsigned int count;

	unsigned int i;


	this->patchConnection.getPolygonSide(minSide, maxSide);

	if (3 <= minSide && 4 >= maxSide)
	{
		count = 0;

		for (viter = this->patchConnection.getVertexSet().begin(); viter != this->patchConnection.getVertexSet().end(); ++viter)
		{
			v1 = *viter;

			if (!this->isControlPoint(v1))
			{
				newPosition = v1->getPosition();

				rs.clear();
				thetas.clear();

				angleSum = 0.0;
				areaSum = 0.0;

				try
				{
					nextEdge = currentEdge = startEdge = *v1->getOutgoingEdges().begin();

					v3 = nextEdge->getTwinEdge()->getOrigin();
					g2 = Geodesic::compute(dcel, v1->getPosition(), v3->getPosition());

					if (NULL == g2)
					{
						throw "Fail to find a geodesic!";
					}

					nextEdge = nextEdge->getPreviousEdge()->getTwinEdge();

					do
					{
						v2 = v3;
						v3 = nextEdge->getTwinEdge()->getOrigin();

						g1 = g2;
						g2 = Geodesic::compute(dcel, v1->getPosition(), v3->getPosition());

						if (NULL == g2)
						{
							throw "Fail to find a geodesic!";
						}

						geodesics.push_back(g1);

						dir1 = g1->forwardDirection() * g1->length();
						dir2 = g2->forwardDirection() * g2->length();

						area = dir1.cross(dir2).length();
						angle = dir1.angle(dir2);

						dir3 = (dir1 + dir2) * (area / 3.0);

						rs.push_back(dir3.length());
						thetas.push_back(dir1.angle(dir3) + angleSum);

						areaSum += area;
						angleSum += angle;

						currentEdge = nextEdge;
						nextEdge = nextEdge->getPreviousEdge()->getTwinEdge();
					} while (currentEdge != startEdge);

					if (0.0 == areaSum)
					{
						throw "Zero Area!";
					}

					Cx = 0.0;
					Cy = 0.0;

					for (i = 0; i < rs.size(); ++i)
					{
						Cx += rs[i] * cos(2.0 * M_PI * thetas[i] / angleSum);
						Cy += rs[i] * sin(2.0 * M_PI * thetas[i] / angleSum);
					}

					Cx /= (2.0 * areaSum);
					Cy /= (2.0 * areaSum);

					Cr = sqrt(Cx * Cx + Cy * Cy);
					//Ctheta = (fmod(atan2(Cy, Cx) + 2.0 * M_PI, 2.0 * M_PI)) / 2.0 / M_PI * angleSum;
					Ctheta = (fmod(atan2(Cy, Cx) / 2.0 / M_PI + 1.0, 1.0)) * angleSum;

					//sum /= (2.0 * areaSum);

					if (Cr > FLT_EPSILON)
					{
						angleSum2 = 0.0;

						giter2 = giter = geodesics.begin();

						if (giter2 != geodesics.end())
						{
							++giter2;
						}

						while (giter2 != geodesics.end())
						{
							dir1 = (*giter)->forwardDirection();
							dir2 = (*giter2)->forwardDirection();

							angle2 = dir1.angle(dir2);
							angleSum2 += angle2;

							if (angleSum2 >= Ctheta)
							{
								break;
							}

							giter = giter2;
							++giter2;
						}

						if (giter2 == geodesics.end())
						{
							giter2 = geodesics.begin();

							dir1 = (*giter)->forwardDirection();
							dir2 = (*giter2)->forwardDirection();

							angle2 = dir1.angle(dir2);
						}

						newPosition = Geodesic::findRelativePositionOnSurface(
							dcel,
							v1->getPosition(), (*giter)->rbegin()->position, (*giter2)->rbegin()->position,
							Cr, angle2 - angleSum2 + Ctheta);
					}

					newPositions.push_back(newPosition);
				}
				catch (const char err[])
				{
					std::cerr << err << std::endl;

					newPositions.push_back(newPosition);
				}

				for (giter = geodesics.begin(); giter != geodesics.end(); ++giter)
				{
					delete *giter;
				}
				geodesics.clear();

				if (++count % 20 == 0)
				{
					printf(".");
					fflush(stdout);
				}
			}
		}

		this->updatePatchConnectionVertices(newPositions);

		printf("\n");
	}
	else
	{
		std::cout << "The stretch minization operator only applies to quad-based parameterization" << std::endl;
	}
}

void PatchList::minimizeAreaDistortionWithConstraint(const DoublyConnectedEdgeList &dcel)
{
	DoublyConnectedEdgeList::VertexSet originalVertexSet;

	std::vector<Geodesic *> geodesics;

	std::vector<double> rs, thetas;

	Geodesic *g1, *g2;

	Vector3D start, direction, newPosition;//, sum;

	Vector3D dir1, dir2, dir3, normal;

	double angleSum, angleSum2, angle, angle2, areaSum, area;
	double Cr, Ctheta, Cx, Cy;

	DoublyConnectedEdgeList::Vertex *v1, *v2, *v3;

	DoublyConnectedEdgeList::Edge *currentEdge, *nextEdge, *startEdge;

	DoublyConnectedEdgeList::VertexSet::const_iterator         viter;
	DoublyConnectedEdgeList::VertexSet::const_reverse_iterator viter2;
	DoublyConnectedEdgeList::EdgeSet::const_iterator         eiter;

	std::vector<Geodesic *>::iterator giter, giter2;

	unsigned int minSide, maxSide;

	unsigned int count;

	unsigned int iteration, i;


	this->patchConnection.getPolygonSide(minSide, maxSide);

	if (3 <= minSide && 4 >= maxSide)
	{
		count = 0;

		for (iteration = 0; iteration < 2; ++iteration)
		{
			originalVertexSet = this->patchConnection.getVertexSet();

			viter = originalVertexSet.begin();
			viter2 = originalVertexSet.rbegin();

			while (viter != originalVertexSet.end())
			{
				v1 = (iteration == 0 ? *viter : *viter2);

				if (!this->isControlPoint(v1))
				{
					newPosition = v1->getPosition();

					rs.clear();
					thetas.clear();

					angleSum = 0.0;
					areaSum = 0.0;

					try
					{
						nextEdge = currentEdge = startEdge = *v1->getOutgoingEdges().begin();

						v3 = nextEdge->getTwinEdge()->getOrigin();
						g2 = Geodesic::compute(dcel, v1->getPosition(), v3->getPosition());

						if (NULL == g2)
						{
							throw "Fail to find a geodesic!";
						}

						if (nextEdge->getPreviousEdge())
						{
							nextEdge = nextEdge->getPreviousEdge()->getTwinEdge();
						}
						else
						{
							nextEdge = NULL;
						}

						do
						{
							if (NULL == nextEdge)
							{
								throw "Open manifold!";
							}

							v2 = v3;
							v3 = nextEdge->getTwinEdge()->getOrigin();

							g1 = g2;
							g2 = Geodesic::compute(dcel, v1->getPosition(), v3->getPosition());

							if (NULL == g2)
							{
								throw "Fail to find a geodesic!";
							}

							geodesics.push_back(g1);

							dir1 = g1->forwardDirection() * g1->length();
							dir2 = g2->forwardDirection() * g2->length();

							area = dir1.cross(dir2).length();

							if (area > FLT_EPSILON)
							{
								angle = dir1.angle(dir2);

								dir3 = (dir1 + dir2) * (area / 3.0);

								rs.push_back(dir3.length());
								thetas.push_back(dir1.angle(dir3) + angleSum);

								areaSum += area;
								angleSum += angle;
							}

							currentEdge = nextEdge;

							if (nextEdge->getPreviousEdge())
							{
								nextEdge = nextEdge->getPreviousEdge()->getTwinEdge();
							}
						} while (currentEdge != startEdge);

						if (0.0 == areaSum)
						{
							throw "Zero Area!";
						}

						Cx = 0.0;
						Cy = 0.0;

						for (i = 0; i < rs.size(); ++i)
						{
							Cx += rs[i] * cos(2.0 * M_PI * thetas[i] / angleSum);
							Cy += rs[i] * sin(2.0 * M_PI * thetas[i] / angleSum);
						}

						Cx /= (2.0 * areaSum);
						Cy /= (2.0 * areaSum);

						Cr = sqrt(Cx * Cx + Cy * Cy);
						//Ctheta = (fmod(atan2(Cy, Cx) + 2.0 * M_PI, 2.0 * M_PI)) / 2.0 / M_PI * angleSum;
						Ctheta = (fmod(atan2(Cy, Cx) / 2.0 / M_PI + 1.0, 1.0)) * angleSum;

						//sum /= (2.0 * areaSum);

						if (Cr > FLT_EPSILON)
						{
							angleSum2 = 0.0;

							giter2 = giter = geodesics.begin();

							if (giter2 != geodesics.end())
							{
								++giter2;
							}

							while (giter2 != geodesics.end())
							{
								dir1 = (*giter)->forwardDirection();
								dir2 = (*giter2)->forwardDirection();

								angle2 = dir1.angle(dir2);
								angleSum2 += angle2;

								if (angleSum2 >= Ctheta)
								{
									break;
								}

								giter = giter2;
								++giter2;
							}

							if (giter2 == geodesics.end())
							{
								giter2 = geodesics.begin();

								dir1 = (*giter)->forwardDirection();
								dir2 = (*giter2)->forwardDirection();

								angle2 = dir1.angle(dir2);
							}

							newPosition = Geodesic::findRelativePositionOnSurface(
								dcel,
								v1->getPosition(), (*giter)->rbegin()->position, (*giter2)->rbegin()->position,
								Cr, angle2 - angleSum2 + Ctheta);

							this->patchConnection.moveVertex(v1, newPosition);
						}
					}
					catch (const char err[])
					{
						std::cerr << err << std::endl;
					}

					for (giter = geodesics.begin(); giter != geodesics.end(); ++giter)
					{
						delete *giter;
					}
					geodesics.clear();

					if (++count % 20 == 0)
					{
						printf(".");
						fflush(stdout);
					}
				}
				++viter;
				++viter2;
			}
		}

		printf("\n");
	}
	else
	{
		std::cout << "The stretch minization operator only applies to quad-based parameterization" << std::endl;
	}
}

#else

Vector3D PatchList::minimizeStretch(const std::list<Geodesic *> &geodesics, DoublyConnectedEdgeList::Vertex *v1)
{
	Geodesic *geodesic;

	Vector3D start, direction, newPosition, sum;

	Vector3D dir1, dir2, normal;

	DoublyConnectedEdgeList::VertexSet::const_iterator viter;
	DoublyConnectedEdgeList::EdgeSet::const_iterator eiter;

	std::list<Geodesic *>::const_iterator giter, giter2;

	unsigned int i;


	if (v1->getOutgoingEdges().empty())
	{
		throw "A vertex does not have any outgoing edge!";
	}

	i = 0;

	for (giter = geodesics.begin(); giter != geodesics.end(); ++giter)
	{
		geodesic = *giter;

		if (0 == i)
		{
			sum += geodesic->forwardDirection() * geodesic->length();
		}
		else
		{
			sum += geodesic->forwardDirection() * (geodesic->length() / M_SQRT2);
		}

		i = (i + 1) % 2;
	}

	sum /= (4.0 * (double)v1->getOutgoingEdges().size());


	return sum;
}

Vector3D PatchList::minimizeAreaDistortion(DoublyConnectedEdgeList::Vertex *v1)
{
	Vector3D start, direction, newPosition;

	Vector3D dir1, dir2, normal;

	Vector3D centroidSum, centroid;

	double areaSum, area;

	Patch *patch;

	Geodesic *geodesic;

	DoublyConnectedEdgeList::Face   *face;
	DoublyConnectedEdgeList::Edge   *currentEdge, *startEdge; //, *nextEdge;

	DoublyConnectedEdgeList::VertexSet::const_iterator viter;
	DoublyConnectedEdgeList::EdgeSet::const_iterator eiter;

	std::list<Geodesic *>::iterator giter, giter2;


	areaSum = 0.0;

	currentEdge = startEdge = *v1->getOutgoingEdges().begin();

	do
	{
		face = currentEdge->getIncidentFace();

		if (NULL == face)
		{
			throw "Null face";
		}

		patch = this->patches[face];

		if (NULL == patch)
		{
			throw "Null patch";
		}

		//patch->mesh->computeCentroid(area, centroid);
		patch->getCentroid(area, centroid);

		geodesic = Geodesic::compute(dcel, v1->getPosition(), centroid);

		if (NULL == geodesic)
		{
			throw "Cannot compute geodesic!";
		}

		areaSum += area;
		centroidSum += geodesic->forwardDirection() * (geodesic->length() * area);

		delete geodesic;

		if (NULL == currentEdge->getPreviousEdge())
		{
			currentEdge = NULL;
		}
		else
		{
			currentEdge = currentEdge->getPreviousEdge()->getTwinEdge();
		}
	} while (currentEdge && currentEdge != startEdge);

	if (0.0 == areaSum)
	{
		throw "Zero Area!";
	}

	if (currentEdge != startEdge)
	{
		throw "Open manifold!";
	}

	centroidSum /= (2.0 * areaSum);
	//centroidSum -= v1->getPosition();
	//centroidSum /= 2.0;


	return centroidSum;
}

void PatchList::optimizeParamerization(const double weight[])
{
	DoublyConnectedEdgeList::VertexSet originalVertexSet;

	std::list<Geodesic *> geodesics;

	Geodesic *geodesic;

	DoublyConnectedEdgeList::Vertex *v1, *v2, *v3;

	Vector3D sum, newPosition;

	Vector3D dir1, dir2, normal;

	unsigned int minSide, maxSide;

	DoublyConnectedEdgeList::VertexSet::const_iterator viter;
	DoublyConnectedEdgeList::VertexSet::const_reverse_iterator viter2;

	DoublyConnectedEdgeList::Edge   *currentEdge, *startEdge; //, *nextEdge;

	std::list<Geodesic *>::iterator giter, giter2;

	unsigned int count;

	unsigned int iteration;


	this->patchConnection.getPolygonSide(minSide, maxSide);

	originalVertexSet = this->patchConnection.getVertexSet();

	count = 0;

	for (iteration = 0; iteration < 2; ++iteration)
	{
		for (viter = originalVertexSet.begin(),
			viter2 = originalVertexSet.rbegin();
			viter != originalVertexSet.end(); ++viter, ++viter2)
		{
			try
			{
				v1 = (iteration == 0 ? *viter : *viter2);

				if (!this->isControlPoint(v1))
				{
					currentEdge = startEdge = *v1->getOutgoingEdges().begin();

					do
					{
						v2 = currentEdge->getNextEdge()->getOrigin();

						geodesic = Geodesic::compute(dcel, v1->getPosition(), v2->getPosition());

						if (NULL == geodesic)
						{
							throw "Fail to find a geodesic!";
						}

						geodesics.push_back(geodesic);

						v3 = currentEdge->getNextEdge()->getNextEdge()->getOrigin();

						geodesic = Geodesic::compute(dcel, v1->getPosition(), v3->getPosition());

						if (NULL == geodesic)
						{
							throw "Fail to find a geodesic!";
						}

						geodesics.push_back(geodesic);

						if (NULL == currentEdge->getPreviousEdge())
						{
							currentEdge = NULL;
						}
						else
						{
							currentEdge = currentEdge->getPreviousEdge()->getTwinEdge();
						}
					} while (currentEdge && currentEdge != startEdge);

					if (currentEdge != startEdge)
					{
						throw "Open manifold";
					}

					sum.setZero();

					if (3 <= minSide && 4 >= maxSide)
					{
						sum += this->minimizeStretch(geodesics, v1) * weight[0];
					}

					if (4 == minSide && 4 == maxSide)
					{
						sum += this->minimizeAreaDistortion(v1) * weight[1];

					}

					//std::cout << sum << std::endl;

					if (sum.length() <= FLT_EPSILON)
					{
						newPosition = v1->getPosition();
					}
					else
					{
						//sum = sum - v1->computeNormal() * sum.dot(v1->computeNormal());

						giter = geodesics.end();
						--giter;

						giter2 = geodesics.begin();

						while (giter2 != geodesics.end())
						{
							dir1 = (*giter)->forwardDirection();
							dir2 = (*giter2)->forwardDirection();

							normal = dir1.cross(dir2);

							if (dir1.cross(sum).dot(normal) >= 0.0 &&
								dir2.cross(sum).dot(normal) <= 0.0)
							{
								break;
							}

							giter = giter2;
							++giter2;
						}

						if (giter2 == geodesics.end())
						{
							//newPositions.push_back(v1->getPosition());
							throw "Error 4!";
						}

						newPosition = Geodesic::findRelativePositionOnSurface(
							dcel,
							v1->getPosition(), (*giter)->rbegin()->position, (*giter2)->rbegin()->position,
							sum);
					}

					//newPositions.push_back(newPosition);

					//std::cout << v1->getPosition() << '\t' << newPosition << std::endl;

					this->patchConnection.moveVertex(v1, newPosition);

					this->syncMeshPoint(v1);

					if (++count % 20 == 0)
					{
						printf(".");
						fflush(stdout);
					}
				}
			}
			catch (const char err[])
			{
				std::cerr << err << std::endl;
			}

			for (giter = geodesics.begin(); giter != geodesics.end(); ++giter)
			{
				delete *giter;
			}
			geodesics.clear();
		}
	}

	printf("\n");
}

#endif

Patch *PatchList::getPatch(const DoublyConnectedEdgeList::Face *f) const
{
	PatchList::FaceToPatchMap::const_iterator iterator = this->patches.find(f);

	if (iterator == this->patches.end())
	{
		return NULL;
	}
	else
	{
		return iterator->second;
	}
}

bool PatchList::isControlPoint(DoublyConnectedEdgeList::Vertex *v) const
{
	return this->controlPoints.find(v) != this->controlPoints.end();
}

bool PatchList::markControlPoint(DoublyConnectedEdgeList::Vertex *v)
{
	bool marked;


	if (!this->isControlPoint(v))
	{
		this->controlPoints.insert(v);

		marked = true;
	}
	else
	{
		marked = false;
	}


	return marked;
}

bool PatchList::unmarkControlPoint(DoublyConnectedEdgeList::Vertex *v)
{
	bool unmarked;


	if (this->isControlPoint(v))
	{
		this->controlPoints.erase(v);

		unmarked = true;
	}
	else
	{
		unmarked = false;
	}


	return unmarked;
}

void PatchList::renderPatches(float epsilon) const
{
	GLdouble projection[16];

	FaceToPatchMap::const_iterator piter;


#if 0
	{
		static Geodesic *g[434][8];

		static bool first = true;

		glGetDoublev(GL_PROJECTION_MATRIX, projection);

		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		glTranslated(0.0, 0.0, -1e-2f);
		glMultMatrixd(projection);
		glMatrixMode(GL_MODELVIEW);

		if (first)
		{
			memset(g, NULL, 434 * 8 * sizeof(void *));

			if (!this->patchConnection.getVertexSet().empty())
			{
				unsigned int i, j;

				DoublyConnectedEdgeList::VertexSet::const_iterator viter = this->patchConnection.getVertexSet().begin();

				for (i = 0; i < this->patchConnection.getVertexSet().size(); ++i)
				{
					//if(i != 51)
					//{
					//    ++viter;
					//    continue;
					//}

					DoublyConnectedEdgeList::EdgeSet edges = (*viter)->getOutgoingEdges();

					j = 0;

					for (DoublyConnectedEdgeList::EdgeSet::const_iterator iter = edges.begin(); iter != edges.end(); ++iter)
					{
						extern DoublyConnectedEdgeList dcel;

						//if(j != 6)
						//{
						//    j += 2;
						//    continue;
						//}

						g[i][j++] = Geodesic::compute(dcel, (*iter)->getOrigin()->getPosition(), (*iter)->getTwinEdge()->getOrigin()->getPosition());
						g[i][j++] = Geodesic::compute(dcel, (*iter)->getOrigin()->getPosition(), (*iter)->getNextEdge()->getNextEdge()->getOrigin()->getPosition());
						j++;
					}

					++viter;
				}

				first = false;
			}
		}

		for (unsigned int i = 0; i < this->patchConnection.getVertexSet().size(); ++i)
		{
			for (unsigned int j = 0; j < 8; ++j)
			{
				if (g[i][j]) g[i][j]->render();
			}
		}

		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		glMatrixMode(GL_MODELVIEW);
	}
#endif

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

	for (piter = this->patches.begin(); piter != this->patches.end(); ++piter)
	{
		piter->second->render();
	}

	if (0.0 != epsilon)
	{
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		glMatrixMode(GL_MODELVIEW);
	}
}


void PatchList::renderControlPoints(float epsilon) const
{
	GLdouble projection[16];

	DoublyConnectedEdgeList::VertexSet::const_iterator viter;


	glPushAttrib(GL_POINT_BIT);

	glPointSize(7.0f);

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

	for (viter = this->controlPoints.begin(); viter != this->controlPoints.end(); ++viter)
	{
		glNormal3dv((*viter)->computeNormal().data());
		glVertex3dv((*viter)->getPosition().data());
	}

	glEnd();

	if (0.0 != epsilon)
	{
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		glMatrixMode(GL_MODELVIEW);
	}

	glPopAttrib();
}
