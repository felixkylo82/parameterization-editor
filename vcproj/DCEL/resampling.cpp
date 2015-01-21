#include "resampling.h"

#include "dcel.h"
#include "OBJModel.h"

#include "common.h"

double Resampling::timeUsed = 0.0;

Resampling::Resampling()
{
}

Resampling::Resampling(const DoublyConnectedEdgeList &dcel)
{
	this->setSurface(dcel);
}

void Resampling::setSurface(const DoublyConnectedEdgeList &dcel)
{
	this->dcel = &dcel;
}

void Resampling::transform(OBJModel &model)
{
	double startTime, endTime;

	startTime = myTime();

	if (this->dcel)
	{
		model.tranform(*this);
		model.collectCloseVertices();
	}

	endTime = myTime();

	Resampling::timeUsed += endTime - startTime;
}

void Resampling::execute(const Vector3D &position, Vector3D &nearestPoint, Vector3D &normal) const
{
	this->rayTrace(position, nearestPoint, normal);
}

void Resampling::findNearestFace(const Vector3D &position, Vector3D &nearestPoint, Vector3D &normal) const
{
	this->dcel->findNearestFace(position, nearestPoint, normal);
}

void Resampling::rayTrace(const Vector3D &position, Vector3D &nearestPoint, Vector3D &normal) const
{
	DoublyConnectedEdgeList::Face* face;

	nearestPoint = position;

	face = this->dcel->rayTrace(position, -normal, nearestPoint, true);

	if (face)
	{
		normal = face->computeNormal().normalize();
	}
}
