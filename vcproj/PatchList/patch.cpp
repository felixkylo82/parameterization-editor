#include "patch.h"

#include "resampling.h"
#include "dcel.h"
#include "OBJModel.h"
#include "vector3D.h"

#include <fstream>

#include <windows.h>
#include <GL/gl.h>


Patch::Patch(unsigned int newNumberOfControlPoints, unsigned int newResolution) :

mesh(new OBJModel()),
numberOfControlPoints(newNumberOfControlPoints),
resolution(newResolution),
centroid(new Vector3D())

{
	this->controlPoints = new Vector3D[this->numberOfControlPoints];
	this->controlNormal = new Vector3D[this->numberOfControlPoints];

	//this->setControlPoints(anchorPoints, anchorNormal);
}

Patch::~Patch()
{
	delete[]this->controlNormal;
	delete[]this->controlPoints;
}

const Vector3D* Patch::getControlPoints() const
{
	return this->controlPoints;
}

const Vector3D* Patch::getControlNormal() const
{
	return this->controlNormal;
}

unsigned int Patch::getNumberOfControlPoints() const
{
	return this->numberOfControlPoints;
}

unsigned int Patch::getResolution() const
{
	return this->resolution;
}

const OBJModel &Patch::getMesh() const
{
	return *this->mesh;
}

void Patch::resampling(const DoublyConnectedEdgeList &dcel)
{
	Resampling(dcel).transform(*this->mesh);

	this->updateCentroid();
}

void Patch::updateCentroid()
{
	this->mesh->computeCentroid(this->area, *this->centroid);

	//std::cout << this->centroid << std::endl;
}

void Patch::getCentroid(double &area, Vector3D &centroid)
{
	area = this->area;
	centroid = *this->centroid;
}

void Patch::render() const
{
	unsigned int i;


	glPushAttrib(GL_POINT_BIT /*| GL_CURRENT_BIT*/);

	//glColor3ub (150, 110, 20);
	glPointSize(7.0f);

	glBegin(GL_POINTS);

	for (i = 0; i < this->numberOfControlPoints; ++i)
	{
		glNormal3dv(this->controlNormal[i].data());
		glVertex3dv(this->controlPoints[i].data());
	}

	glEnd();

	//glColor3ub (20, 150, 20);
	this->mesh->render();

	glPopAttrib();
}


QuadPatch::QuadPatch(const Vector3D newControlPoints[4], const Vector3D newControlNormal[4], unsigned int newResolution) :

Patch(4, newResolution)

{
	//memcpy(this->controlPoints, newControlPoints, 4 * sizeof(Vector3D));

	//this->mesh->genQuad(this->controlPoints, this->resolution);

	this->setControlPoints(newControlPoints, newControlNormal);
}

QuadPatch::~QuadPatch()
{
}

void QuadPatch::setResolution(int newResolution)
{
	if (this->resolution != newResolution)
	{
		this->mesh->clear();
		this->mesh->genQuad(this->controlPoints, this->controlNormal, newResolution);

		this->resolution = newResolution;

		this->updateCentroid();
	}
}

void QuadPatch::setControlPoints(const Vector3D newControlPoints[4], const Vector3D newControlNormal[4])
{
	//unsigned int i;


	//for(i = 0; i < 4; ++i)
	//{
	//    this->controlPoints[i] = anchorPoints[i];
	//}

	memcpy(this->controlPoints, newControlPoints, 4 * sizeof(Vector3D));
	memcpy(this->controlNormal, newControlNormal, 4 * sizeof(Vector3D));

	this->mesh->clear();
	this->mesh->genQuad(this->controlPoints, this->controlNormal, this->resolution);

	this->updateCentroid();
}

TriPatch::TriPatch(const Vector3D newControlPoints[3], const Vector3D newControlNormal[3], unsigned int newResolution) :

Patch(3, newResolution)

{
	//memcpy(this->controlPoints, newControlPoints, 3 * sizeof(Vector3D));

	//this->mesh->genTriangle(this->controlPoints, this->resolution);

	this->setControlPoints(newControlPoints, newControlNormal);
}

TriPatch::~TriPatch()
{
}

void TriPatch::setResolution(int newResolution)
{
	if (this->resolution != newResolution)
	{
		this->mesh->clear();
		this->mesh->genTriangle(this->controlPoints, this->controlNormal, newResolution);

		this->resolution = newResolution;

		this->updateCentroid();
	}
}

void TriPatch::setControlPoints(const Vector3D newControlPoints[3], const Vector3D newControlNormal[3])
{
	memcpy(this->controlPoints, newControlPoints, 3 * sizeof(Vector3D));
	memcpy(this->controlNormal, newControlNormal, 3 * sizeof(Vector3D));

	this->mesh->clear();
	this->mesh->genTriangle(this->controlPoints, this->controlNormal, this->resolution);

	this->updateCentroid();
}
