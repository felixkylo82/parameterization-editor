#ifndef PATCH_LIST_H
#define PATCH_LIST_H


#define PATCHLIST_MAJOR_VERSION 2


#include "dcel.h"


class Patch;
class Geodesic;


class PatchList
{
public:

    typedef std::map<const DoublyConnectedEdgeList::Face *, Patch *> FaceToPatchMap;

    PatchList (const DoublyConnectedEdgeList &dcel);
    ~PatchList();

    //void 
    void clear ();

    bool loadFile(const char filename[])      ;
    void saveFile(const char filename[]) const;

    const DoublyConnectedEdgeList &getPatchConnection() const;
    const FaceToPatchMap          &getFaceToPatchMap () const;
    const DoublyConnectedEdgeList &getOriginalSurface() const;

    void addPatch   (Patch *patch);
    void removePatch(DoublyConnectedEdgeList::Face *face);

    int  getPatchResolution() const;
    void setPatchResolution(int resolution);

    void subdividePatchList(PatchList &subdividedPatchList) const;

    void resampling();

    void syncMeshPoint (const DoublyConnectedEdgeList::Vertex *v);
    void syncMeshPoints();

    void outputMesh(OBJModel &model) const;

    void moveControlPoint(DoublyConnectedEdgeList::Vertex * vertex, const Vector3D &to, float epsilon = 0.0f);

    void upateIDs();

    void projectParameterization(OBJModel &initialParameterization, unsigned int gridDensity);

#if (PATCHLIST_MAJOR_VERSION == 1)

    void minimizeStretch       (const DoublyConnectedEdgeList &dcel);

    void minimizeAreaDistortion              (const DoublyConnectedEdgeList &dcel);
    void minimizeAreaDistortionAccurately    (const DoublyConnectedEdgeList &dcel);
    void minimizeAreaDistortionWithConstraint(const DoublyConnectedEdgeList &dcel);

#elif (PATCHLIST_MAJOR_VERSION == 2)

    Vector3D minimizeStretch       (const std::list<Geodesic *> &geodesics, DoublyConnectedEdgeList::Vertex *v1);
    Vector3D minimizeAreaDistortion(DoublyConnectedEdgeList::Vertex *v1);

    void optimizeParamerization(const double weight[]);

#endif

    Patch *getPatch(const DoublyConnectedEdgeList::Face *) const;

    bool isControlPoint    (DoublyConnectedEdgeList::Vertex *) const;
    bool markControlPoint  (DoublyConnectedEdgeList::Vertex *);
    bool unmarkControlPoint(DoublyConnectedEdgeList::Vertex *);

    void renderPatches      (float epsilon = 0.0f) const;
    void renderControlPoints(float epsilon = 0.0f) const;

protected:

    void updatePatchConnectionVertices(const std::list<Vector3D> &newPositions);

private:

    const DoublyConnectedEdgeList &dcel;

    DoublyConnectedEdgeList patchConnection;

    FaceToPatchMap patches;

    DoublyConnectedEdgeList::VertexSet controlPoints;

    int resolution;
};


#endif // PATCH_LIST_H
