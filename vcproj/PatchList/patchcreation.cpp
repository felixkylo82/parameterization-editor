#include "patchcreation.h"

#include "patchlist.h"
#include "patch.h"
#include "OBJModel.h"

#include "windows.h"
#include "GL/gl.h"


PatchCreation::PatchCreation(PatchList &newPatchList) :

currentState          (ready        ),
patchList             (&newPatchList),
numberOfAnchorPoints  (0            ),
autoReverseOrientation(true         ),
//resolution            (1            ),

autoSave(true),

maximumUndo  (10),
undoLevel    (-1),
latestLevel  (-1),
possibleLevel(-1)

//epsilon             (5e-2f        )
{
    this->setUndoPoint();
}

PatchCreation::State PatchCreation::getCurrentState() const
{
    return this->currentState;
}

bool PatchCreation::isReady() const
{
    return (PatchCreation::ready == this->currentState);
}

void PatchCreation::begin(enum State state)
{
    if(PatchCreation::ready != this->currentState)
    {
        throw "Inconsistent state!";
    }

    this->currentState         = state;
    this->numberOfAnchorPoints = 0    ;
}

void PatchCreation::end()
{
    if(/*0 != this->numberOfAnchorPoints || */PatchCreation::ready == this->currentState)
    {
        this->currentState = PatchCreation::ready;
        throw "Inconsistent state!";
    }

    this->currentState = PatchCreation::ready;
}

void PatchCreation::vertex(const Vector3D &position)
{
    bool reversed;

    Patch *patch;

    Vector3D normal;
    Vector3D tmp;


    if(PatchCreation::ready ==this->currentState)
    {
        throw "Inconsistent state!";
    }

    this->anchorPoints[this->numberOfAnchorPoints] = position;

    //std::cout << "Anchor point position = " << position << std::endl;

    ++this->numberOfAnchorPoints;

    switch(this->currentState)
    {
    case PatchCreation::quad:

        if(4 == this->numberOfAnchorPoints)
        {
            if(this->autoReverseOrientation)
            {
                normal = Vector3D::computeNormal(anchorPoints[0], anchorPoints[1], anchorPoints[2]);

                if(!normal.isZero() && !this->isOrientationConsistent(anchorPoints[1], normal.normalize()))
                {
                    anchorPoints[0].swap(anchorPoints[3]);
                    anchorPoints[1].swap(anchorPoints[2]);
                }
            }

            patch = new QuadPatch(anchorPoints, anchorPoints, this->patchList->getPatchResolution());

            patch->resampling(this->patchList->getOriginalSurface());
            this->patchList->addPatch(patch/*, epsilon*/);

            this->numberOfAnchorPoints = 0;

            std::cout << "A quad created" << std::endl;

            if(this->autoSave)
            {
                this->setUndoPoint();
            }
        }

        break;

    case PatchCreation::quadStrip:

        if(4 == this->numberOfAnchorPoints)
        {
            tmp = anchorPoints[2];
            anchorPoints[2] = anchorPoints[3];
            anchorPoints[3] = tmp;

            if(this->autoReverseOrientation)
            {
                normal = Vector3D::computeNormal(anchorPoints[0], anchorPoints[1], anchorPoints[2]);

                if(!normal.isZero() && !this->isOrientationConsistent(anchorPoints[1], normal.normalize()))
                {
                    anchorPoints[0].swap(anchorPoints[3]);
                    anchorPoints[1].swap(anchorPoints[2]);

                    reversed = true;
                }
                else
                {
                    reversed = false;
                }
            }

            patch = new QuadPatch(anchorPoints, anchorPoints, this->patchList->getPatchResolution());

            if(this->autoReverseOrientation)
            {
                if(reversed)
                {
                    anchorPoints[0].swap(anchorPoints[3]);
                    anchorPoints[1].swap(anchorPoints[2]);
                }
            }

            patch->resampling(this->patchList->getOriginalSurface());
            this->patchList->addPatch(patch/*, epsilon*/);

            anchorPoints[0] = anchorPoints[3];
            anchorPoints[1] = anchorPoints[2];

            this->numberOfAnchorPoints = 2;

            std::cout << "A quad created" << std::endl;

            if(this->autoSave)
            {
                this->setUndoPoint();
            }
        }

        break;

    case PatchCreation::triangle:

        //this->anchorPoints[this->numberOfAnchorPoints] = position;
        //++this->numberOfAnchorPoints;

        if(3 == this->numberOfAnchorPoints)
        {
            if(this->autoReverseOrientation)
            {
                normal = Vector3D::computeNormal(anchorPoints[0], anchorPoints[1], anchorPoints[2]);

                if(!normal.isZero() && !this->isOrientationConsistent(anchorPoints[1], normal.normalize()))
                {
                    anchorPoints[1].swap(anchorPoints[2]);
                }
            }

            patch = new TriPatch(anchorPoints, anchorPoints, this->patchList->getPatchResolution());

            patch->resampling(this->patchList->getOriginalSurface());
            this->patchList->addPatch(patch/*, epsilon*/);

            this->numberOfAnchorPoints = 0;

            std::cout << "A triangle created" << std::endl;

            if(this->autoSave)
            {
                this->setUndoPoint();
            }
        }

        break;

    case PatchCreation::triangleStrip:

        switch(this->numberOfAnchorPoints)
        {
        case 3:

            anchorPoints[1].swap(anchorPoints[2]);

            if(this->autoReverseOrientation)
            {
                normal = Vector3D::computeNormal(anchorPoints[0], anchorPoints[1], anchorPoints[2]);

                if(!normal.isZero() && !this->isOrientationConsistent(anchorPoints[1], normal.normalize()))
                {
                    anchorPoints[1].swap(anchorPoints[2]);

                    reversed = true;
                }
                else
                {
                    reversed = false;
                }
            }

            patch = new TriPatch(anchorPoints, anchorPoints, this->patchList->getPatchResolution());

            if(this->autoReverseOrientation)
            {
                if(reversed)
                {
                    anchorPoints[1].swap(anchorPoints[2]);
                }
            }

            patch->resampling(this->patchList->getOriginalSurface());
            this->patchList->addPatch(patch/*, epsilon*/);

            anchorPoints[1].swap(anchorPoints[2]);

            std::cout << "A triangle created" << std::endl;

            if(this->autoSave)
            {
                this->setUndoPoint();
            }

            break;

        case 4:

            if(this->autoReverseOrientation)
            {
                normal = Vector3D::computeNormal(anchorPoints[1], anchorPoints[2], anchorPoints[3]);

                if(!normal.isZero() && !this->isOrientationConsistent(anchorPoints[2], normal.normalize()))
                {
                    anchorPoints[2].swap(anchorPoints[3]);

                    reversed = true;
                }
                else
                {
                    reversed = false;
                }
            }

            patch = new TriPatch(anchorPoints + 1, anchorPoints + 1, this->patchList->getPatchResolution());

            if(this->autoReverseOrientation)
            {
                if(reversed)
                {
                    anchorPoints[2].swap(anchorPoints[3]);
                }
            }

            patch->resampling(this->patchList->getOriginalSurface());
            this->patchList->addPatch(patch/*, epsilon*/);

            anchorPoints[0] = anchorPoints[2];
            anchorPoints[1] = anchorPoints[3];

            this->numberOfAnchorPoints = 2;

            std::cout << "A triangle created" << std::endl;

            if(this->autoSave)
            {
                this->setUndoPoint();
            }

            break;
        }
    }

    std::cout << "The number of anchor points = " << this->numberOfAnchorPoints << std::endl;
}

bool PatchCreation::devertex()
{
    bool remove;


    remove = (this->numberOfAnchorPoints > 0);

    if(remove)
    {
        --this->numberOfAnchorPoints;
    }


    return remove;
}

OBJModel *PatchCreation::generateParameterization(unsigned int type, const unsigned int resolution[])
{
    OBJModel *initialParameterization;


    initialParameterization = new OBJModel();

    switch(type)
    {
    case 0:
        initialParameterization->genRhombicDodechedron(resolution);
        break;
    case 1:
        initialParameterization->genIcosahedron(resolution[0]);
        break;
    case 2:
        initialParameterization->genBox(resolution);
        break;
    }

    //initialParameterization->normalize();
    //initialParameterization->scale(Vector3D(M_SQRT2, M_SQRT2, M_SQRT2));
    //initialParameterization->scale(Vector3D(8.0, 8.0, 8.0));

    initialParameterization->scale(Vector3D(2.0, 2.0, 2.0));


    return initialParameterization;
}

void PatchCreation::projectParameterization(OBJModel &initialParameterization)
{
    //this->patchList->clear();
    this->patchList->projectParameterization(initialParameterization, this->patchList->getPatchResolution());
}

//float PatchCreation::getEpsilon() const
//{
//    return this->epsilon;
//}
//
//void PatchCreation::setEpsilon(float newEpsilon)
//{
//    this->epsilon = newEpsilon;
//}

bool PatchCreation::isOrientationConsistent(const Vector3D &position, const Vector3D &normal) const
{
    Vector3D nearestPoint, relatedNormal;


    this->patchList->getOriginalSurface().findNearestFace(position, nearestPoint, relatedNormal);


    return (normal.dot(relatedNormal) >= 0);
}

bool PatchCreation::getAutoReverseOrientation() const
{
    return this->autoReverseOrientation;
}

void PatchCreation::setAutoReverseOrientation(bool enable)
{
    this->autoReverseOrientation = enable;
}

#if 0

int PatchCreation::getResolution() const
{
    return this->resolution;
}

void PatchCreation::setResolution(int newResolution)
{
    if(this->resolution != newResolution)
    {
        this->patchList->setPatchResolution(newResolution );
        this->patchList->resampling        (*this->surface);

        this->resolution = newResolution;
    }
}

#endif

void PatchCreation::renderPendingPoints() const
{
    unsigned int i;


    if(!this->ready)
    {
        glBegin(GL_POINTS);

        for(i = 0; i < this->numberOfAnchorPoints; ++i)
        {
            glVertex3dv(this->anchorPoints[i].data());
        }

        glEnd();
    }
}

bool PatchCreation::getAutoSave() const
{
    return this->autoSave;
}

void PatchCreation::setAutoSave(bool enable)
{
    this->autoSave = enable;
}

void PatchCreation::setUndoPoint()
{
    char filename[128];


    this->undoLevel = this->latestLevel = (this->undoLevel + 1) % this->maximumUndo;

    if(this->latestLevel > this->possibleLevel)
    {
        this->possibleLevel = this->latestLevel;
    }

    sprintf_s(filename, "temp/undo.%d.patch", this->latestLevel);

    this->patchList->saveFile(filename);
}

bool PatchCreation::undo()
{
    char filename[128];

    bool ret;


    ret = false;

    if( (this->undoLevel - this->latestLevel + this->maximumUndo - 1) % this->maximumUndo > 0 &&
        (this->undoLevel + this->maximumUndo - 1) % this->maximumUndo <= this->possibleLevel)
    {
        this->undoLevel = (this->undoLevel + this->maximumUndo - 1) % this->maximumUndo;

        //printf("%d\n", this->undoLevel);

        sprintf_s(filename, "temp/undo.%d.patch", this->undoLevel);

        if(this->patchList->loadFile(filename))
        {
            //this->patchList->setResolution(this->resolution);
            this->patchList->resampling();

            ret = true;
        }
    }


    return ret;
}

bool PatchCreation::redo()
{
    char filename[128];

    bool ret;


    ret = false;

    if((this->latestLevel - this->undoLevel + this->maximumUndo) % this->maximumUndo > 0 &&
        (this->undoLevel + 1) % this->maximumUndo <= this->possibleLevel)
    {
        this->undoLevel = (this->undoLevel + 1) % this->maximumUndo;

        //printf("%d\n", this->undoLevel);

        sprintf_s(filename, "temp/undo.%d.patch", this->undoLevel);

        if(this->patchList->loadFile(filename))
        {
            //this->patchList->setResolution(this->resolution);
            this->patchList->resampling();

            ret = true;
        }
    }


    return ret;
}
