#include "patchcreation.h"
#include "patchlist.h"
#include "patch.h"
#include "geodesic.h"
#include "resampling.h"
#include "dcel.h"
#include "OBJModel.h"

//#include <cmath>
#include <iostream>
#include <cstring>

#include <windows.h>
#include <GL/glut.h>


//////////////////////////////////////////////////////////////////////////////////////////
//
// Working Directory: $(SolutionDir)
//
// Arguments:
//
// data/models/squirrel.obj data/patchlists/squirrel.patch
// data/models/duck.obj     data/patchlists/duck.patch
// data/models/bimba.obj    data/patchlists/bimba.patch
//
//////////////////////////////////////////////////////////////////////////////////////////


char modelFilePrefix[512];
int  patch_sequence_no;

OBJModel                model;
DoublyConnectedEdgeList dcel;

PatchList     patchList(dcel);
PatchCreation pc(patchList);

OBJModel *initialParameterization = NULL;

GLuint model_display_list = 0;

Vector3D mouse_obj_pos;

int mouseButton = 0;
int previousX;
int previousY;

//int selectMode  = 0;
//int selectState = 0;

DoublyConnectedEdgeList::Primitive *selectedPrimitive;

//int draw_patch  = 1;
int draw_points = 0;
int draw_edges = 0;
int draw_faces = 1;
int draw_param_edges = 1;
int draw_param_faces = 0;

GLsizei windowWidth = 1024;
GLsizei windowHeight = 1024;



float epsilon = 2e-2f;


void resetView()
{
	glLoadIdentity();
	gluLookAt(
		0.0, 0.0, 10.0,
		//5.0,  5.0,  5.0,
		0.0, 0.0, 0.0,
		0.0, 1.0, 0.0);
}

void init()
{
	const GLfloat material_ambient[] = { 0.25f, 0.25f, 0.25f, 1.0f };
	const GLfloat material_diffuse[] = { 0.70f, 0.70f, 0.70f, 1.0f };
	const GLfloat material_specular[] = { 0.05f, 0.05f, 0.05f, 1.0f };
	const GLfloat material_shininess = 10.0f;

	const GLfloat light_position[] = { 0.0f, 1.0f, 1.0f, 0.0f };
	const GLfloat light_ambient[] = { 0.8f, 0.8f, 0.8f, 0.0f };
	const GLfloat light_diffuse[] = { 1.0f, 1.0f, 1.0f, 0.0f };
	const GLfloat light_specular[] = { 1.0f, 1.0f, 1.0f, 0.0f };


	glMaterialfv(GL_FRONT, GL_AMBIENT, material_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, material_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, material_specular);
	glMaterialf(GL_FRONT, GL_SHININESS, material_shininess);

	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);

	glEnable(GL_LIGHT0);
	//glDisable(GL_LIGHT0);
	glEnable(GL_LIGHTING);
	//glDisable(GL_LIGHTING);

	glEnable(GL_COLOR_MATERIAL);

	glEnable(GL_NORMALIZE);

	glEnable(GL_CULL_FACE);
	//glCullFace(GL_FRONT);

	glEnable(GL_DEPTH_TEST);

	glShadeModel(GL_SMOOTH);

	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_POINT_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_LINE_SMOOTH);

	glEnable(GL_BLEND);

	//glClearColor(1.0, 1.0, 1.0, 0.0);
	glClearColor(0.0, 0.0, 0.0, 0.0);

	glPointSize(3.0f);

	resetView();
}

void reshape(int width, int height)
{
	windowWidth = width;
	windowHeight = height;

	glViewport(0, 0, (GLsizei)width, (GLsizei)height);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(25.0, (GLdouble)width / (GLdouble)height, 1.0, 50.0);
	glMatrixMode(GL_MODELVIEW);
}

bool getMousePosition(int x, int y, int snap, Vector3D &mouse_obj_pos)
{
	GLint    viewport[4];
	GLdouble projmatrix[16];
	GLdouble mvmatrix[16];

	GLint w, h;

	GLfloat *depth_buffer;
	GLfloat  depth;

	GLdouble obj_pos[3];

	Vector3D nearestPoint;

	DoublyConnectedEdgeList::Vertex    *v;
	DoublyConnectedEdgeList::Primitive *p;

	bool ret;


	ret = true;

	glGetIntegerv(GL_VIEWPORT, viewport);
	glGetDoublev(GL_PROJECTION_MATRIX, projmatrix);
	glGetDoublev(GL_MODELVIEW_MATRIX, mvmatrix);

	w = viewport[2];
	h = viewport[3];


	///////////////////////////
	// Read the depth buffer

	depth_buffer = new GLfloat[w * h * sizeof(GLfloat)];

	if (depth_buffer)
	{
		glClear(GL_DEPTH_BUFFER_BIT);

		//dcel.render();
		glCallList(model_display_list);

		glReadBuffer(GL_BACK);
		glPixelStorei(GL_PACK_ALIGNMENT, 1);
		glDepthRange(0.0, 1.0);
		glReadPixels(0, 0, w, h, GL_DEPTH_COMPONENT, GL_FLOAT, depth_buffer);

		///////////////////////////
		// Read depth value at (x,y)

		depth = depth_buffer[x + (h - y - 1) * w];

		///////////////////////////
		// Call gluUnProject

		ret = depth != 1.0 && gluUnProject(
			(GLdouble)x, (GLdouble)(h - y - 1), depth,
			mvmatrix, projmatrix, viewport,
			obj_pos, obj_pos + 1, obj_pos + 2);

		//if(!ret)
		//{
		//    fprintf ( stderr , "Warning (gluUnProject returns FALSE)...\n" ) ;
		//}

		delete[]depth_buffer;
	}
	else
	{
		fprintf(stderr, "readDepth : Not enough memory to allocate for depth!\n");

		ret = false;
	}

	mouse_obj_pos = Vector3D(obj_pos);

	if (ret)
	{
		Vector3D mouse_normal;

		switch (snap)
		{
		case 1:

			v = patchList.getPatchConnection().findNearestVertex(mouse_obj_pos);
			//p = dcel.findNearestFace(mouse_obj_pos, nearestPoint, mouse_normal);

			if (NULL != v)
			{
				nearestPoint = v->getPosition();

				if ((nearestPoint - mouse_obj_pos).length() < epsilon)
				{
					selectedPrimitive = v;

					mouse_obj_pos = nearestPoint;
					//std::cout << "SELCTED a point on a patch at " << mouse_obj_pos << std::endl;
					std::cout << "SELCTED a point" << std::endl;
				}
			}

			break;

		case 2:

			p = patchList.getPatchConnection().findNearestFace(mouse_obj_pos, nearestPoint, mouse_normal, epsilon);
			//p = dcel.findNearestFace(mouse_obj_pos, nearestPoint, mouse_normal);

			if (NULL != p)
			{
				//if((nearestPoint - mouse_obj_pos).length() < 10.0 * epsilon)
				{
					selectedPrimitive = p;

					if (DoublyConnectedEdgeList::Primitive::vertex == p->getPrimitiveType())
					{
						mouse_obj_pos = nearestPoint;
						//std::cout << "SELCTED a point on a patch at " << mouse_obj_pos << std::endl;
						std::cout << "SELCTED a point" << std::endl;
					}
					else if (DoublyConnectedEdgeList::Primitive::edge == p->getPrimitiveType())
					{
						std::cout << "SELCTED a edge" << std::endl;
					}
					else if (DoublyConnectedEdgeList::Primitive::face == p->getPrimitiveType())
					{
						std::cout << "SELCTED a face" << std::endl;
					}
				}
			}

			break;

		default:
			break;
		}
	}


	return ret;
}

void mouse(int button, int state, int x, int y)
{
	if (state == GLUT_DOWN)
	{
		mouseButton |= (button << 1 | (button == GLUT_LEFT_BUTTON));

		switch (button)
		{
		case GLUT_LEFT_BUTTON:
			break;

		case GLUT_MIDDLE_BUTTON:

			getMousePosition(x, y, 2, mouse_obj_pos);

			break;

		case GLUT_RIGHT_BUTTON:

			getMousePosition(x, y, 1, mouse_obj_pos);

			if (!pc.isReady())
			{
				pc.vertex(Vector3D(mouse_obj_pos));
			}

			break;
		}
	}
	else
	{
		switch (button)
		{
		case GLUT_RIGHT_BUTTON:

			if (pc.isReady())
			{
				if (selectedPrimitive &&
					DoublyConnectedEdgeList::Primitive::vertex == selectedPrimitive->getPrimitiveType() &&
					getMousePosition(x, y, 0, mouse_obj_pos))
				{
					patchList.moveControlPoint((DoublyConnectedEdgeList::Vertex *) selectedPrimitive, mouse_obj_pos, epsilon);
					pc.setUndoPoint();
					selectedPrimitive = NULL;
				}
			}
		};

		mouseButton &= ~(button << 1 | (button == GLUT_LEFT_BUTTON));
	}

	previousX = x;
	previousY = y;

	glutPostRedisplay();
}

void motion(int x, int y)
{
	const float rotateScale = 80.0f;
	const float translateScale = 10.0f;

	GLfloat modelview[16];

	double rotation[9];

	float dx, dy;
	float nx, ny;
	float magnitude;
	float angle;


	if (x != previousX || y != previousY)
	{
		switch (mouseButton)
		{
		case 1:

			dx = ((float)x - (float)previousX) / (float)windowHeight;
			dy = ((float)previousY - (float)y) / (float)windowWidth;

			magnitude = sqrt(dx * dx + dy * dy);
			dx /= magnitude;
			dy /= magnitude;

			angle = rotateScale * magnitude;

			nx = -dy;
			ny = dx;

			if (NULL == initialParameterization)
			{
				glGetFloatv(GL_MODELVIEW_MATRIX, modelview);

				glLoadIdentity();
				glTranslatef(modelview[12], modelview[13], modelview[14]);
				glRotatef(angle, nx, ny, 0.0f);
				glTranslatef(-modelview[12], -modelview[13], -modelview[14]);
				glMultMatrixf(modelview);
			}
			else
			{
				glPushMatrix();
				glLoadIdentity();
				glRotatef(angle, nx, ny, 0.0f);

				glGetFloatv(GL_MODELVIEW_MATRIX, modelview);

				rotation[0] = modelview[0]; rotation[1] = modelview[4]; rotation[2] = modelview[8];
				rotation[3] = modelview[1]; rotation[4] = modelview[5]; rotation[5] = modelview[9];
				rotation[6] = modelview[2]; rotation[7] = modelview[6]; rotation[8] = modelview[10];

				initialParameterization->rotate(rotation);

				glPopMatrix();
			}

			break;

		case 2:

			dx = translateScale * ((float)x - (float)previousX) / (float)windowHeight;
			dy = translateScale * ((float)previousY - (float)y) / (float)windowWidth;

			if (NULL == initialParameterization)
			{
				glGetFloatv(GL_MODELVIEW_MATRIX, modelview);

				glLoadIdentity();
				glTranslatef(dx, dy, 0.0f);
				glMultMatrixf(modelview);
			}
			else
			{
				initialParameterization->translate(Vector3D(dx, dy, 0.0f));
			}

			break;

		case 3:

			dx = translateScale * ((float)x - (float)previousX) / (float)windowHeight;
			dy = translateScale * ((float)previousY - (float)y) / (float)windowWidth;

			if (NULL == initialParameterization)
			{
				glGetFloatv(GL_MODELVIEW_MATRIX, modelview);

				glLoadIdentity();
				glTranslatef(0.0f, 0.0f, -dy);
				glMultMatrixf(modelview);
			}
			else
			{
				initialParameterization->translate(Vector3D(0.0f, 0.0f, -dy));
			}

			break;

		case 4:

			if (pc.isReady())
			{
				getMousePosition(x, y, 0, mouse_obj_pos);
			}

			break;
		}

		previousX = x;
		previousY = y;

		glutPostRedisplay();
	}
}

void freeMemory()
{
	if (0 != model_display_list)
	{
		glDeleteLists(model_display_list, 1);
	}
}

void keyboard(unsigned char key, int x, int y)
{
	const Vector3D meshPoints[4] =
	{
		Vector3D(1.2, -0.2, -0.2),
		Vector3D(1.2, 0.2, -0.2),
		Vector3D(1.2, 0.2, 0.2),
		Vector3D(1.2, -0.2, 0.2),
	};

	char filename[512];

	DoublyConnectedEdgeList::Vertex *v;

	unsigned int resolution[4] = { 4, 4, 4, 4, };


	//printf("%d\n", key);

	switch (key)
	{
	case 27:

		freeMemory();
		exit(0);

		break;

	case '1':
	case '2':
	case '3':

		if (NULL != initialParameterization)
		{
			delete initialParameterization;
		}

		switch (key)
		{
		case '1':
			scanf_s("%d %d %d %d", resolution, resolution + 1, resolution + 2, resolution + 3);
			initialParameterization = pc.generateParameterization(0, resolution);
			break;
		case '2':
			scanf_s("%d", resolution);
			initialParameterization = pc.generateParameterization(1, resolution);
			break;
		case '3':
			scanf_s("%d %d %d", resolution, resolution + 1, resolution + 2);
			initialParameterization = pc.generateParameterization(2, resolution);
			break;
		}

		if (NULL == initialParameterization)
		{
			printf("Not enough memory!\n");
		}

		break;

	case '+':

		patchList.setPatchResolution(2 * patchList.getPatchResolution());

		patchList.resampling();

		printf("Resampling resolution = %d\n", patchList.getPatchResolution());

		break;

	case '-':

		if (patchList.getPatchResolution() >= 2)
		{
			patchList.setPatchResolution(patchList.getPatchResolution() / 2);

			patchList.resampling();
		}

		printf("Resampling resolution = %d\n", patchList.getPatchResolution());

		break;

	case '[':

		if (pc.undo())
		{
			selectedPrimitive = NULL;

			printf("undo geometry\n");
		}
		else
		{
			printf("fail to undo geometry\n");
		}

		break;

	case ']':

		if (pc.redo())
		{
			selectedPrimitive = NULL;

			printf("redo geometry\n");
		}
		else
		{
			printf("fail to redo geometry\n");
		}

		break;

	case 'A':

		pc.setAutoSave(!pc.getAutoSave());

		if (pc.getAutoSave())
		{
			printf("Enable autosave\n");
		}
		else
		{
			printf("Disable autosave\n");
		}

		break;

	case 'a':

		pc.setAutoReverseOrientation(!pc.getAutoReverseOrientation());

		if (pc.getAutoReverseOrientation())
		{
			printf("Enable automatic flipping triangle (/ quad) orientation\n");
		}
		else
		{
			printf("Disable automatic flipping triangle (/ quad) orientation\n");
		}

		break;

	case 'D':

		if (NULL != selectedPrimitive && DoublyConnectedEdgeList::Primitive::face == selectedPrimitive->getPrimitiveType())
		{
			patchList.removePatch((DoublyConnectedEdgeList::Face *)selectedPrimitive);
			pc.setUndoPoint();
			selectedPrimitive = NULL;

			printf("A face removed\n");
		}

		break;

	case 'd':

		if (pc.devertex())
		{
			printf("The last inserted point removed\n");
		}
		else
		{
			printf("Fail to remove a point\n");
		}

		break;

	case 'E':

		draw_param_edges = (draw_param_edges + 1) % 2;

		break;

	case 'e':

		draw_edges = (draw_edges + 1) % 2;

		if (draw_edges)
		{
			printf("Show edges\n");
		}
		else
		{
			printf("Hide edges\n");
		}

		break;

	case 'f':

		draw_faces = (draw_faces + 1) % 2;

		if (draw_faces)
		{
			printf("Show faces\n");
		}
		else
		{
			printf("Hide faces\n");
		}

		break;

	case 'F':

		dcel.flipFaceOrientation(*dcel.getFaceSet().begin());
		printf("Flip model face orientation\n");

		if (0 == model_display_list)
		{
			model_display_list = glGenLists(1);
		}

		glNewList(model_display_list, GL_COMPILE);
		dcel.render();
		glEndList();

		break;

	case 'G':

		if (NULL != initialParameterization)
		{
			pc.projectParameterization(*initialParameterization);

			delete initialParameterization;
			initialParameterization = NULL;

			pc.setUndoPoint();
		}

		break;

	case 'H':
	case 'h':

		printf("\n");
		printf("========== Manual ==========\n");
		printf("Esc    exit\n");
		printf("1      generate a rhombic dodechedron parameterization\n");
		printf("2      generate an icosahedron parameterization\n");
		printf("3      generate a cube parameterization\n");
		printf("+/-    Increase (/ decrease) the resampling resolution\n");
		printf("[/]    Undo (/ redo) geometry\n");
		printf("A      Enable (/ disable) auto-save\n");
		printf("a      Enable (/ disable) automatic flipping triangle (/ quad) orientation\n");
		printf("d      Remove an inserted point\n");
		printf("D      Remove a selected face\n");
		printf("e      Show model edges\n");
		printf("f      Show model faces\n");
		printf("F      Flip model face orientation\n");
		printf("G      Project parameterization onto the surface\n");
		printf("H/h    Manual\n");
		printf("k      Mark (/ unmark) the selected point as a anchor point\n");
		printf("L/l    Load a patch list\n");
		printf("m      Reduce parameterization stretch\n");
		//printf("p      Show model points\n");
		printf("Q      Begin (/ end) a quad strip\n");
		printf("q      Begin (/ end) quads\n");
		printf("r      Reset view\n");
		printf("S      Save a model dcel\n");
		printf("s      Save a patch list\n");
		printf("T      Begin (/ end) a triangle strip\n");
		printf("t      Begin (/ end) triangles\n");
		printf("V      Validate the model dcel\n");
		printf("v      Validate the patch list dcel\n");
		printf("u      Update search trees\n");
		printf("============================\n");

		break;

	case 'k':

		if (NULL != selectedPrimitive)
		{
			if (DoublyConnectedEdgeList::Primitive::vertex == selectedPrimitive->getPrimitiveType())
			{
				v = (DoublyConnectedEdgeList::Vertex *) selectedPrimitive;

				if (patchList.isControlPoint(v))
				{
					patchList.unmarkControlPoint(v);

					printf("Unarked an anchor point\n");
				}
				else
				{
					patchList.markControlPoint(v);

					printf("Marked an anchor point\n");
				}
			}
		}

		break;

	case 'L':
	case 'l':

		sprintf_s(filename, "%s.%d.patch", modelFilePrefix, patch_sequence_no - 1);

		printf("Loading the patch file [%s]\n", filename);

		patchList.loadFile(filename);
		patchList.resampling();

		break;

	case 'M':

		//printf("Reduce parameterization area distortion\n");
		printf("Reduce parameterization area distortion (x50)\n");

		for (unsigned int i = 0; i < 50; ++i)
		{
			printf("iteration %d", i);

			//patchList.minimizeAreaDistortionAccurately(dcel);

			const double weight[] = { 0.25, 0.75, };

			patchList.optimizeParamerization(weight);

			printf("Time for geodesic optimization   = %f\n", Geodesic::timeForGeodesicOptimization);
			printf("Time for searching nearest faces = %f\n", Geodesic::timeForSearchingNearestFace);
			printf("Time for resampling              = %f\n", Resampling::timeUsed);

			extern void display();
			display();

			sprintf_s(filename, "temp/int.%d.patch", i);

			patchList.saveFile(filename);
		}
		//patchList.syncMeshPoints(dcel);

		pc.setUndoPoint();

		break;

#if 0

	case 'm':

		//printf("Reduce parameterization stretch\n");
		printf("Reduce parameterization stretch (x50)\n");

		for (unsigned int i = 0; i < 50; ++i)
		{
			printf("iteration %d", i);

			patchList.minimizeStretch(dcel);

			printf("Time for geodesic optimization   = %f\n", Geodesic::timeForGeodesicOptimization);
			printf("Time for searching nearest faces = %f\n", Geodesic::timeForSearchingNearestFace);
			printf("Time for resampling              = %f\n", Resampling::timeUsed);

			//if(i % 5 == 4)
			//{
			//    extern void display();
			//    patchList.syncMeshPoints(dcel);
			//    display();
			//}

			sprintf_s(filename, "temp/int.%d.patch", i);

			patchList.saveFile(filename);
		}
		patchList.syncMeshPoints(dcel);

		pc.setUndoPoint();

		break;

#endif

	case 'p':

		draw_points = (draw_points + 1) % 2;

		if (draw_points)
		{
			printf("Show points\n");
		}
		else
		{
			printf("Hide points\n");
		}

		break;

	case 'Q':

		if (pc.isReady())
		{
			printf("Begin a quad strip\n");
			pc.begin(PatchCreation::quadStrip);
		}
		else
		{
			printf("End a quad strip\n");
			pc.end();
		}

		break;

	case 'q':

		if (pc.isReady())
		{
			printf("Begin quads\n");
			pc.begin(PatchCreation::quad);
		}
		else
		{
			printf("End quads\n");
			pc.end();
		}

		break;

		//case 'R':
		//    break;

	case 'r':

		printf("Reset view\n");
		resetView();

		break;

	case 'S':

		sprintf_s(filename, "%s.model.obj", modelFilePrefix);

		printf("Saving the model dcel[%s]\n", filename);

		dcel.saveFile(filename);

		break;

	case 's':

		sprintf_s(filename, "%s.%d.patch", modelFilePrefix, patch_sequence_no++);

		printf("Saving the patch list[%s]\n", filename);

		patchList.saveFile(filename);

		break;

	case 'T':

		if (pc.isReady())
		{
			printf("Begin a triangle strip\n");
			pc.begin(PatchCreation::triangleStrip);
		}
		else
		{
			printf("End a triangle strip\n");
			pc.end();
		}

		break;

	case 't':

		if (pc.isReady())
		{
			printf("Begin triangles\n");
			pc.begin(PatchCreation::triangle);
		}
		else
		{
			printf("End triangles\n");
			pc.end();
		}

		break;

	case 'V':

		printf("Validating the model DCEL\n");
		dcel.validate();

		break;

	case 'v':

		printf("Validating the patch list DCEL\n");
		patchList.getPatchConnection().validate();

		break;

	case 'u':

		std::cout << "Update search trees" << std::endl;
		dcel.constructSearchTrees(1);

		break;

	case 'w':
		break;
	}

	glutPostRedisplay();
}

void display()
{
	const float white[] = { 1.0f, 1.0f, 1.0f, 1.0f, };

	GLdouble projection[16];


	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


	if (0 == model_display_list)
	{
		model_display_list = glGenLists(1);
		glNewList(model_display_list, GL_COMPILE);
		dcel.render();
		glEndList();
	}

	{
		glColor3ub(130, 40, 40);
		glPointSize(7.0f);

		glGetDoublev(GL_PROJECTION_MATRIX, projection);

		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		glTranslated(0.0, 0.0, -6e-3f);
		glMultMatrixd(projection);
		glMatrixMode(GL_MODELVIEW);

#ifdef RENDER_CORNERS
		{
			DoublyConnectedEdgeList::VertexSet::const_iterator viter;

			for (viter = patchList.getPatchConnection().getVertexSet().begin(); viter != patchList.getPatchConnection().getVertexSet().end(); ++viter)
			{
				if (4 != (*viter)->getNumberOfNeighbors())
				{
					(*viter)->render();
				}
			}
		}
#endif

		if (!pc.isReady())
		{
			pc.renderPendingPoints();
		}
		else if (NULL != selectedPrimitive)
		{
			//glBegin    (GL_POINTS);
			//glVertex3dv(mouse_obj_pos.data());
			//glEnd      ();
			switch (selectedPrimitive->getPrimitiveType())
			{
			case DoublyConnectedEdgeList::Primitive::vertex:

				selectedPrimitive->render();

				break;

			case DoublyConnectedEdgeList::Primitive::face:

				//patchList.getPatch((DoublyConnectedEdgeList::Face *)selectedPrimitive)->updateCentroid();
				patchList.getPatch((DoublyConnectedEdgeList::Face *)selectedPrimitive)->render();

				break;

			default:

				break;
			}
		}

		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		glMatrixMode(GL_MODELVIEW);

		glPointSize(1.0f);
	}

	if (draw_param_edges)
	{
		glDisable(GL_LIGHTING);

		glColor3ub(140, 100, 40);

		patchList.renderControlPoints(7e-3f);

		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

		//glColor3ub (150, 60, 60);
		//glLineWidth(2.0);

		//patchList.getPatchConnection().render(7e-3f);

		//glPointSize(2.0f);
		glColor3ub(50, 50, 140);

		if (NULL != initialParameterization)
		{
			initialParameterization->render(6e-3f);
		}
		else
		{
			patchList.renderPatches(6e-3f);
		}
	}

	if (draw_param_faces)
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		//glPointSize(1.0f);

		glEnable(GL_LIGHTING);

		glColor3ub(120, 120, 150);


		if (NULL != initialParameterization)
		{
			initialParameterization->render(5e-3f);
		}
		else
		{
			patchList.renderPatches(5e-3f);
		}
	}

	if (draw_points)
	{
	}

	if (draw_edges)
	{
		glDisable(GL_LIGHTING);
		glColor3ub(20, 20, 20);
		glLineWidth(2.0f);

		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

		glGetDoublev(GL_PROJECTION_MATRIX, projection);

		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		//glTranslated  (0.0, 0.0, -5e-5f);
		glTranslated(0.0, 0.0, -7e-3f);
		glMultMatrixd(projection);
		glMatrixMode(GL_MODELVIEW);

		glCallList(model_display_list);

		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		glMatrixMode(GL_MODELVIEW);

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		glLineWidth(1.0f);
		glEnable(GL_LIGHTING);
	}

	if (draw_faces)
	{
		glColor3ub(130, 130, 130);
		//model.render();
		//dcel.render();
		glCallList(model_display_list);
	}

	glutSwapBuffers();

	glutReportErrors();
}


int main(int argc, char** argv)
{
	patch_sequence_no = 0;

	if (argc == 1)
	{
		strcpy_s(modelFilePrefix, "sphere");

		unsigned int resolution[] = { 8, 8, 8, 8, };

		//model.genBrick ( );
		//model.genBox   (1);
		//model.genIcosahedron(8);
		model.genRhombicDodechedron(resolution);

		//std::cout << "Loading an OBJ model" << std::endl;
		//model.loadFile(argv[1]);

		std::cout << "Converting the model into DCEL" << std::endl;
		dcel.convertFromOBJModel(model);

		//dcel.validate();
	}
	else if (argc >= 2)
	{
		strcpy_s(modelFilePrefix, argv[1]);
		modelFilePrefix[strlen(modelFilePrefix) - 4] = '\0';

		std::cout << "Loading an OBJ model" << std::endl;
		model.loadFile(argv[1]);

		//model.reverseOrientation();

		std::cout << "Normalizing the model" << std::endl;
		model.normalize();

		//std::cout << "Collecting close points" << std::endl;
		//model.collectCloseVertices();

		std::cout << "Converting the model into DCEL" << std::endl;
		dcel.convertFromOBJModel(model);

		std::cout << "Removing points without any connection in the DCEL" << std::endl;
		dcel.removeDisconnVert();

		//std::cout << "Converting the model to a close manifold" << std::endl;
		//dcel.toCloseManifold();
	}

	std::cout << "Filling missing polygons" << std::endl;
	dcel.fillMissingPolygons();

	std::cout << "Converting the model into a triangular mesh if it is not" << std::endl;
	dcel.triangulate();

	std::cout << "Construct search trees" << std::endl;
	dcel.constructSearchTrees(1);

	std::cout << "Validating the model DCEL" << std::endl;
	if (dcel.validate())
	{
		//if(argc <= 2)
		//{
		//    std::cout << "Building a spherical parameterization" << std::endl;
		//    patchList.sphericalParameterization(dcel, 4, 1);
		//}
		//else
		if (argc > 2)
		{
			std::cout << "Loading the patch file" << std::endl;
			patchList.loadFile(argv[2]);
			patchList.resampling();
		}

		//glutInit(&argc, argv);

		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);

		glutInitWindowSize(windowWidth, windowHeight);
		glutCreateWindow(argv[0]);

		init();

		glutDisplayFunc(display);
		glutReshapeFunc(reshape);
		glutMouseFunc(mouse);
		glutMotionFunc(motion);
		glutKeyboardFunc(keyboard);

		glutMainLoop();
	}
	else
	{
		std::cerr << "The DCEL is corrupted!" << std::endl;
	}

	freeMemory();


	return 0;
}
