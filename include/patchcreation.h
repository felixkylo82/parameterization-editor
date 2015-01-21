#ifndef PATCH_CREATION_H
#define PATCH_CREATION_H


#include "vector3D.h"


class OBJModel;
class PatchList;


class PatchCreation
{
public:

	enum State
	{
		ready, triangle, quad, triangleStrip, quadStrip
	};

	PatchCreation(PatchList &patchList);

	State getCurrentState() const;
	bool  isReady() const;

	void begin(enum State state);
	void end();

	void vertex(const Vector3D &position);
	bool devertex();

	OBJModel *generateParameterization(unsigned int type, const unsigned int patchDensity[]);
	void      projectParameterization(OBJModel &initialParameterization);

	bool isOrientationConsistent(const Vector3D &position, const Vector3D &normal) const;

	bool getAutoReverseOrientation() const;
	void setAutoReverseOrientation(bool enable);

	void renderPendingPoints() const;

	bool getAutoSave() const;
	void setAutoSave(bool enable);

	void setUndoPoint();
	bool undo();
	bool redo();

private:

	enum State  currentState;
	PatchList  *patchList;

	Vector3D anchorPoints[4];
	unsigned int numberOfAnchorPoints;

	bool autoReverseOrientation;


	bool autoSave;

	const int maximumUndo;
	int       possibleLevel;
	int       undoLevel;
	int       latestLevel;
};


#endif // PATCH_CREATION_H
