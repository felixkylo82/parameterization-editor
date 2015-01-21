# parameterization-editor

I. Directories
========================================

1.  bin\
output for executables
2.  data\models
OBJ models
3.  data\patchlists
files storing lists of patches in OBJ format
4.  include\
header files shared by vc projects
5.  lib\
output for libraries
6.  temp\
intermediate output of parametarization editor
7.  vcproj\DCEL
a library project for doubly connected edge lists and geodesics
8.  vcproj\OBJModel
a library project for simple OBJ model manipulation
9.  vcproj\ParameterizationEditor
an executable project for parameterization editor
10. vcproj\PatchList
a library project for patch list manipulation
11. vcproj\RayTrace
a library project for ray tracing and intersection tests
12. vcproj\utilites and vcproj\cUtilites
library projects for commonly used functions


II. Source Files
========================================

1.  vcproj\DCEL\dcel.cpp
doubly connected edge list implemetation
2.  vcproj\DCEL\geodesic.cpp
compute geodesics on a surface represented in a doubly connected edge list
3.  vcproj\DCEL\resampling.cpp
move points of a OBJ model to a surface represented in a doubly connected edge list
4.  vcproj\OBJModel\OBJModel.cpp
manipluate obj models
5.  vcproj\ParameterizationEditor\ParameterizationEditor.cpp
a simple GLUT program for editing a patch list and performing parameterization optimization
6.  vcproj\PatchList\patch.cpp
manipulate a patch
7.  vcproj\PatchList\patchList.cpp
manipulate a patch list
8.  vcproj\PatchList\patchcreation.cpp
help create patches on a surface represented in a doubly connected edge list
9.  vcproj\RayTrace\raytrace.cpp
commonly used functions for ray tracing and intersection tests
10. vcproj\RayTrace\boundingbox.cpp
naive implementation of hierarchical bounding box
11. vcproj\utilities\jacobi.cpp
diagonalize symmetric matrices
12. vcproj\utilities\vector3D.cpp
operation on 3D vectors


III. Header Files
========================================

1.  include\boundingbox.h
2.  include\common.h
3.  include\dcel.h
4.  include\dumpscreen.h
5.  include\geodesic.h
6.  include\jacobi.h
7.  include\OBJModel.h
8.  include\patch.h
9.  include\patchcreation.h
10. include\patchlist.h
11. include\raytrace.h
12. include\resampling.h
13. include\standard.h
14. include\vector3D.h
15. include\vptree.h


IV. External Libraries
========================================

1.  GLUT


V. Globals
========================================

1.  OBJModel                  model
an input surface
2.  DoublyConnectedEdgeList   dcel
the DCEL of the input surface
3.  PatchList                 patchList
patches which lie on the input surface
4.  PatchCreation             pc
A helper for creating patches and modifying patches
5.  OBJModel                 *initialParameterization
an initial parameterization in OBJ format
6.  GUI states x 12


VI. Patches & Patch Lists
========================================

1.  A patch is a polygon, typically a quadrilateral, whose vertices lie on a given surface represented in a doubly connected edge list
2.  A patch list is a list of patches. Connected patches may form another mesh and the connectivity is, again, repersented in a doubly connected edge list


VII. Creating Parameterization
========================================

1. Start up with arguments <surface in OBJ format> <patch list in OBJ format>
With the function, int main(int argc, char** argv), the program
normalizes the surface,
converts the surface into a DCEL,
removes points without any connection in the DCEL,
converts the DCEL into a triangular mesh if it is not,
validates the DCEL,
loads the patch file, and
initializes openGL states and create a GLUT window

2. Edit a parameterization with the program

3. Keep minimize stretch
With the function, void PatchList::optimizeParamerization(const double weight[]), the program,
for each vertex in the vertex set of the patch connection, 
computes all the geodesics connecting it to its neighbors,
computes a tangential vector, T, which indicate the magnitude and the direction of movement of this vertex,
computes a geodesic whose first derivative matches T, and
moves this vertex along the geodesics such that the magnitude of movement matches ||T||

* Project parameterization for recovery
With the function, void PatchList::projectParameterization(OBJModel &initialParameterization, unsigned int gridDensity), the program
takes input points on an intial parameterization in OBJ format,
projects the initial parameterization onto a given surface represented in a doubly connected list, and
outputs the final parameterization in a doubly connected edge list where each face is a patch and the list represent patch connection as a whole


VIII. GUI Key Menu
========================================

Esc    exit
1      generate a rhombic dodechedron parameterization
2      generate an icosahedron parameterization
3      generate a cube parameterization
+/-    Increase (/ decrease) the resampling resolution
[/]    Undo (/ redo) geometry
A      Enable (/ disable) auto-save
a      Enable (/ disable) automatic flipping triangle (/ quad) orientation
d      Remove an inserted point
D      Remove a selected face
e      Show model edges
f      Show model faces
F      Flip model face orientation
G      Project parameterization onto the surface
H/h    Manual
k      Mark (/ unmark) the selected point as a anchor point
L/l    Load a patch list
m      Reduce parameterization stretch
Q      Begin (/ end) a quad strip
q      Begin (/ end) quads
r      Reset view
S      Save a model dcel
s      Save a patch list
T      Begin (/ end) a triangle strip
t      Begin (/ end) triangles
V      Validate the model dcel
v      Validate the patch list dcel
u      Update search trees


IX. KNOWN BUGS OR PROBLEMS
========================================

1.  In some rare cases, geodesics cannot pass the validate() function which checks a list of properties of valid geodesics
2.  OBJModel.cpp has only handled a limited number of tokens
3.  When converting a OBJ model into doubly connected using convertFromOBJModel(), vertices of each polygon have to be ordered in anti-clockwise order. Automatic ordering has been implemented but has not been tested completely.
