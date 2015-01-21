#include "OBJModel.h"

#include "Transformation.h"
#include "vector3D.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <cfloat>

#include <windows.h>
#include <GL/gl.h>


OBJModel::OBJModel()
{
}

OBJModel::OBJModel(const OBJModel &model) :

vertexArray(model.vertexArray),
normalArray(model.normalArray),
textureCoordinateArray(model.textureCoordinateArray)

{
	std::vector <PolygonElement *>::const_iterator pitr;


	for (pitr = model.polygons.begin(); pitr != model.polygons.end(); ++pitr)
	{
		this->polygons.push_back(new PolygonElement(**pitr));

		//std::cout << this->polygons.back()->vertexIndices.size() << std::endl;
	}
}

OBJModel::OBJModel(std::ifstream &is)
{
	this->loadFile(is);
}

OBJModel::OBJModel(const char filename[])
{
	this->loadFile(filename);
}

OBJModel::~OBJModel()
{
	this->clear();
}

bool OBJModel::loadFile(const char filename[])
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

bool OBJModel::loadFile(std::ifstream &fin)
{
	char line[256];

	char *currentCharacter;
	std::vector<char *> tokens, tokens2;

	Vector3D vector;
	PolygonElement *polygon;

	std::vector<PolygonElement *>::const_iterator piter;

	int index;

	unsigned int i, j;

	bool ret;


	//////////////////////////////////////////////////
	// 1) Read OBJ model
	//////////////////////////////////////////////////

	polygon = NULL;
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

						this->vertexArray.push_back(vector);
					}
					else
					{
						throw "invalid format";
					}
				}
				else if (strcmp(tokens[0], "vt") == 0)
				{
					//////////////////////////////////////////////////
					// 1b) Read edge information
					//////////////////////////////////////////////////

					if (tokens.size() == 3)
					{
						for (i = 0; i < 2; ++i)
						{
							vector[i] = atof(tokens[i + 1]);

							if (vector[i] == HUGE_VAL)
							{
								throw "floating point overflow";
							}
						}

						this->textureCoordinateArray.push_back(vector);
					}
					else
					{
						throw "invalid format";
					}
				}
				else if (strcmp(tokens[0], "vn") == 0)
				{
					//////////////////////////////////////////////////
					// 1b) Read normal information
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

						this->normalArray.push_back(vector);
					}
					else
					{
						throw "invalid format";
					}
				}
				else if (strcmp(tokens[0], "f") == 0)
				{
					//////////////////////////////////////////////////
					// 1c) Read face information
					//////////////////////////////////////////////////

					polygon = new PolygonElement;

					for (i = 1; i < tokens.size(); ++i)
					{
						tokens2.clear();

						currentCharacter = tokens[i];

						if (*currentCharacter == '/')
						{
							throw "invalid format";
						}

						while (true)
						{
							tokens2.push_back(currentCharacter);

							while (*currentCharacter != '/' && *currentCharacter != '\0')
							{
								++currentCharacter;
							}

							if (*currentCharacter != '\0')
							{
								*currentCharacter = '\0';
								++currentCharacter;

								while (*currentCharacter != '\0' && (*currentCharacter < '0' || *currentCharacter > '9'))
								{
									++currentCharacter;
								}
							}
							else
							{
								break;
							}
						}

						for (j = 0; j < tokens2.size(); ++j)
						{
							index = atoi(tokens2[j]);

							if (index > 0)
							{
								--index;

								switch (j)
								{
								case 0:
									polygon->vertexIndices.push_back(index);
									break;

								case 1:
									polygon->textureCoordinateIndices.push_back(index);
									break;

								case 2:
									polygon->normalIndices.push_back(index);
									break;

								default:
									throw "invalid file format";
								}
							}
							else
							{
								printf("%s\n", tokens2[j]);

								throw "non positive index not handled";
							}
						}
					}

					if (polygon->vertexIndices.size() >= 3 &&
						(polygon->normalIndices.size() == polygon->vertexIndices.size() || polygon->normalIndices.size() == 0) &&
						(polygon->textureCoordinateIndices.size() == polygon->vertexIndices.size() || polygon->textureCoordinateIndices.size() == 0))
					{
						this->polygons.push_back(polygon);
					}
					else
					{
						delete polygon;
					}

					polygon = NULL;
				}
				else if (tokens[0][0] == '#')
				{
					// skip comment
				}
				else if (tokens[0][0] == 'g')
				{
					// skip group token
				}
				else if (tokens[0][0] == 's')
				{
					// skip 's' token
				}
				else// if(strcmp(token[0], "\0") != 0)
				{
					// skip unhandled token
					throw "unhandled token";
				}
			}
		}
	}
	catch (const char error[])
	{
		if (polygon)
		{
			delete polygon;
		}

		this->clear();

		std::cerr << error << std::endl;

		ret = false;
	}

	return ret;
}

void OBJModel::saveFile(const char filename[]) const
{
	std::ofstream fout;


	fout.open(filename);

	this->saveFile(fout);

	if (fout.is_open())
	{
		fout.close();
	}
}

void OBJModel::saveFile(std::ofstream &fout) const
{
	std::vector<Vector3D        >::const_iterator viter;
	std::vector<PolygonElement *>::const_iterator piter;

	PolygonElement *polygon;

	unsigned int i;


	try
	{
		if (!fout.is_open())
		{
			throw "file cannot be opened";
		}

		for (viter = this->vertexArray.begin(); viter != this->vertexArray.end(); ++viter)
		{
			fout << "v " << (*viter)[0] << ' ' << (*viter)[1] << ' ' << (*viter)[2] << '\n';
		}

		if (this->vertexArray.size() != 0)
		{
			fout << std::endl;
		}

		for (viter = this->textureCoordinateArray.begin(); viter != this->textureCoordinateArray.end(); ++viter)
		{
			fout << "vt " << (*viter)[0] << ' ' << (*viter)[1] << '\n';
		}

		if (this->textureCoordinateArray.size() != 0)
		{
			fout << std::endl;
		}

		for (viter = this->normalArray.begin(); viter != this->normalArray.end(); ++viter)
		{
			fout << "vn " << (*viter)[0] << ' ' << (*viter)[1] << ' ' << (*viter)[2] << '\n';
		}

		if (this->normalArray.size() != 0)
		{
			fout << std::endl;
		}

		for (piter = this->polygons.begin(); piter != this->polygons.end(); ++piter)
		{
			fout << "f ";

			polygon = *piter;

			for (i = 0; i < polygon->vertexIndices.size(); ++i)
			{
				fout << polygon->vertexIndices[i] + 1 << '/';

				if (polygon->vertexIndices.size() == polygon->textureCoordinateIndices.size())
				{
					fout << polygon->textureCoordinateIndices[i] + 1;
				}
				fout << '/';

				if (polygon->vertexIndices.size() == polygon->normalIndices.size())
				{
					fout << polygon->normalIndices[i] + 1;
				}
				fout << ' ';
			}

			fout << '\n';
		}
	}
	catch (const char error[])
	{
		std::cout << error << std::endl;
	}
}

void OBJModel::clear()
{
	std::vector<PolygonElement *>::const_iterator piter;

	this->vertexArray.clear();
	this->normalArray.clear();
	this->textureCoordinateArray.clear();

	for (piter = this->polygons.begin(); piter != this->polygons.end(); ++piter)
	{
		delete (*piter);
	}

	this->polygons.clear();
}

size_t OBJModel::getNumberOfVertices() const
{
	return this->vertexArray.size();
}

size_t OBJModel::getNumberOfNormal() const
{
	return this->normalArray.size();
}

size_t OBJModel::getNumberOfTextureCoordinates() const
{
	return this->textureCoordinateArray.size();
}

size_t OBJModel::getNumberOfPolygons() const
{
	return this->polygons.size();
}

const PolygonElement *OBJModel::getPolygon(size_t i) const
{
	return this->polygons[i];
}

const Vector3D &OBJModel::getVertex(size_t i) const
{
	return this->vertexArray[i];
}

const Vector3D &OBJModel::getNormal(size_t i) const
{
	return this->normalArray[i];
}

const Vector3D &OBJModel::getTextureCoordinate(size_t i) const
{
	return this->textureCoordinateArray[i];
}

void OBJModel::reverseOrientation()
{
	std::vector <PolygonElement *>::iterator iter;


	for (iter = this->polygons.begin(); iter != this->polygons.end(); ++iter)
	{
		std::reverse((*iter)->normalIndices.begin(), (*iter)->normalIndices.end());
		std::reverse((*iter)->textureCoordinateIndices.begin(), (*iter)->textureCoordinateIndices.end());
		std::reverse((*iter)->vertexIndices.begin(), (*iter)->vertexIndices.end());
	}
}

size_t OBJModel::addVertex(const Vector3D &position)
{
	size_t index;


	index = this->vertexArray.size();

	this->vertexArray.push_back(position);


	return index;
}

size_t OBJModel::addNormal(const Vector3D &normal)
{
	size_t index;


	index = this->normalArray.size();

	this->normalArray.push_back(normal);


	return index;
}

size_t OBJModel::addTextureCoordinate(const Vector3D &textureCoordinate)
{
	size_t index;


	index = this->textureCoordinateArray.size();

	this->textureCoordinateArray.push_back(textureCoordinate);


	return index;
}

size_t OBJModel::addPolygon(const PolygonElement &polygon)
{
	size_t index;


	index = this->polygons.size();


	this->polygons.push_back(new PolygonElement(polygon));


	return index;
}

void OBJModel::rotate(const double rotationMatrix[])
{
	Vector3D position, normal;

	Vector3D x, y, z;

	unsigned int i;


	x = Vector3D(rotationMatrix + 0);
	y = Vector3D(rotationMatrix + 3);
	z = Vector3D(rotationMatrix + 6);

	for (i = 0; i < this->vertexArray.size(); ++i)
	{
		position[0] = x.dot(this->vertexArray[i]);
		position[1] = y.dot(this->vertexArray[i]);
		position[2] = z.dot(this->vertexArray[i]);

		this->vertexArray[i] = position;
	}

	for (i = 0; i < this->normalArray.size(); ++i)
	{
		normal[0] = x.dot(this->normalArray[i]);
		normal[1] = y.dot(this->normalArray[i]);
		normal[2] = z.dot(this->normalArray[i]);

		this->normalArray[i] = normal;
	}
}

void OBJModel::translate(const Vector3D &displacement)
{
	std::vector<Vector3D>::iterator viter;


	for (viter = this->vertexArray.begin(); viter != this->vertexArray.end(); ++viter)
	{
		(*viter) += displacement;
	}
}

void OBJModel::scale(const Vector3D &scalingFactor)
{
	std::vector<Vector3D>::iterator viter;


	for (viter = this->vertexArray.begin(); viter != this->vertexArray.end(); ++viter)
	{
		(*viter)[0] = (*viter)[0] * scalingFactor[0];
		(*viter)[1] = (*viter)[1] * scalingFactor[1];
		(*viter)[2] = (*viter)[2] * scalingFactor[2];
	}
}

#if 0

void OBJModel::normalize()
{
	std::vector<Vector3D>::iterator viter, viter2;

	Vector3D min, max, centre;
	double scale_x, scale_y, scale;


	/*
	for(viter = this->vertexArray.begin(); viter != this->vertexArray.end(); ++viter)
	{
	for(viter2 = viter; viter2 != this->vertexArray.end(); ++viter2)
	{
	if(dist3D((*viter), viter2->coordinate) <= 5e-3)
	{
	add3Df ((*viter), viter2->coordinate, viter2->coordinate);
	mult3Df(viter2->coordinate, 0.5);
	copy3Df((*viter), viter2->coordinate);
	}
	}
	}
	*/

	if (this->vertexArray.size() > 0)
	{
		viter = this->vertexArray.begin();

		//copy3Df(min, (*viter));
		//copy3Df(max, (*viter));

		max = min = *viter;

		++viter;

		while (viter != this->vertexArray.end())
		{
			if (min[0] > (*viter)[0])
			{
				min[0] = (*viter)[0];
			}
			else if (max[0] < (*viter)[0])
			{
				max[0] = (*viter)[0];
			}

			if (min[1] > (*viter)[1])
			{
				min[1] = (*viter)[1];
			}
			else if (max[1] < (*viter)[1])
			{
				max[1] = (*viter)[1];
			}

			if (min[2] > (*viter)[2])
			{
				min[2] = (*viter)[2];
			}
			else if (max[2] < (*viter)[2])
			{
				max[2] = (*viter)[2];
			}

			++viter;
		}

		//add3Df (min, max, centre);
		//mult3Df(centre, 0.5);

		centre = (min + max) / 2.0;

		scale_x = max[0] - min[0];
		scale_y = max[1] - min[1];

		scale = (scale_x > scale_y ? scale_x : scale_y);
		scale /= 2.0;

		for (viter = this->vertexArray.begin(); viter != this->vertexArray.end(); ++viter)
		{
			(*viter)[0] = ((*viter)[0] - centre[0]) / scale;
			(*viter)[1] = ((*viter)[1] - centre[1]) / scale;
			(*viter)[2] = ((*viter)[2] - centre[2]) / scale;
		}

		for (viter = this->normalArray.begin(); viter != this->normalArray.end(); ++viter)
		{
			//unify3Df(*viter);
			*viter = viter->normalize();
		}
	}
}

#endif

double OBJModel::normalize()
{
	Vector3D centroid;

	Vector3D min, max;
	double scale_x, scale_y, scale_z, scale;
	double length, area;

	std::vector<PolygonElement *>::const_iterator piter;

	std::vector<Vector3D>::iterator viter, viter2;


	if (!this->vertexArray.empty())
	{
		this->computeCentroid(area, centroid);

		for (viter = this->vertexArray.begin(); viter != this->vertexArray.end(); ++viter)
		{
			*viter -= centroid;
		}

		viter = this->vertexArray.begin();

		max = min = *viter;

		++viter;

		while (viter != this->vertexArray.end())
		{
			if (min[0] > (*viter)[0])
			{
				min[0] = (*viter)[0];
			}
			else if (max[0] < (*viter)[0])
			{
				max[0] = (*viter)[0];
			}

			if (min[1] > (*viter)[1])
			{
				min[1] = (*viter)[1];
			}
			else if (max[1] < (*viter)[1])
			{
				max[1] = (*viter)[1];
			}

			if (min[2] > (*viter)[2])
			{
				min[2] = (*viter)[2];
			}
			else if (max[2] < (*viter)[2])
			{
				max[2] = (*viter)[2];
			}

			++viter;
		}

		scale_x = max[0] - min[0];
		scale_y = max[1] - min[1];
		scale_z = max[2] - min[2];

		scale = scale_x;
		if (scale_y > scale) scale = scale_y;
		if (scale_z > scale) scale = scale_z;
		scale /= 2.0;

		this->scale(Vector3D(1.0 / scale, 1.0 / scale, 1.0 / scale));

		for (viter = this->normalArray.begin(); viter != this->normalArray.end(); ++viter)
		{
			length = viter->length();

			if (length <= FLT_EPSILON)
			{
				viter->setZero();
			}
			else
			{
				*viter /= length;
			}
		}
	}

	return scale;
}

void OBJModel::collectCloseVertices(float epsilon)
{
	unsigned int i, j;


	for (i = 0; i < this->vertexArray.size(); ++i)
	{
		for (j = i + 1; j < this->vertexArray.size(); ++j)
		{
			if ((this->vertexArray[i] - this->vertexArray[j]).length() <= epsilon)
			{
				this->vertexArray[j] = this->vertexArray[i];
				this->normalArray[j] = this->normalArray[i];
			}
		}
	}
}

void OBJModel::computeCentroid(double &area, Vector3D &centroid) const
{
	const PolygonElement *polygon;

	std::vector<PolygonElement *>::const_iterator piter;

	double triArea;

	unsigned int i;


	area = 0.0;
	centroid.setZero();

	for (piter = this->polygons.begin(); piter != this->polygons.end(); ++piter)
	{
		polygon = *piter;

		for (i = 2; i < polygon->vertexIndices.size(); ++i)
		{
			triArea = (
				this->vertexArray[polygon->vertexIndices[i - 1]] - this->vertexArray[polygon->vertexIndices[0]]).cross(
				this->vertexArray[polygon->vertexIndices[i]] - this->vertexArray[polygon->vertexIndices[i - 1]]).length() / 2.0;

			area += triArea;

			centroid += (
				this->vertexArray[polygon->vertexIndices[0]] +
				this->vertexArray[polygon->vertexIndices[i - 1]] +
				this->vertexArray[polygon->vertexIndices[i]]) * (triArea / 3.0);
		}
	}

	if (0.0 != area)
	{
		centroid /= area;
	}
}

void OBJModel::render(float epsilon) const
{
	const PolygonElement *polygon;

	GLdouble projection[16];

	std::vector<PolygonElement *>::const_iterator piter;

	unsigned int i;


	glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT);

	if (this->normalArray.size() == 0)
	{
		glDisable(GL_LIGHTING);
	}

	//glDisable(GL_CULL_FACE);

	//glEnable(GL_COLOR_MATERIAL);

	//glColor3ub(230, 230, 230);

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

	for (piter = this->polygons.begin(); piter != this->polygons.end(); ++piter)
	{
		polygon = *piter;

		glBegin(GL_POLYGON);

		//glNormal3dv(polygon->faceNormal);

		for (i = 0; i < polygon->vertexIndices.size(); ++i)
		{
			if (this->normalArray.size())
			{
				glNormal3dv(this->normalArray[polygon->normalIndices[i]].data());
			}
			glVertex3dv(this->vertexArray[polygon->vertexIndices[i]].data());
		}

		glEnd();
	}

	if (0.0 != epsilon)
	{
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		glMatrixMode(GL_MODELVIEW);
	}

	glPopAttrib();
}

#if 0

void OBJModel::genQuad(unsigned int orientation, unsigned int resolution)
{
	PolygonElement *polygon;

	unsigned int offset;
	unsigned int i, j;


	offset = (unsigned int)this->vertexArray.size();

	for (i = 0; i < resolution + 1; ++i)
	{
		for (j = 0; j < resolution + 1; ++j)
		{
			switch (orientation)
			{
			case 0:
				this->vertexArray.push_back(Vector3D(0.5f, (float)i / (float)resolution - 0.5f, (float)j / (float)resolution - 0.5f));
				this->normalArray.push_back(Vector3D(1.0f, 0.0f, 0.0f));
				break;
			case 1:
				this->vertexArray.push_back(Vector3D(-0.5f, (float)j / (float)resolution - 0.5f, (float)i / (float)resolution - 0.5f));
				this->normalArray.push_back(Vector3D(-1.0f, 0.0f, 0.0f));
				break;
			case 2:
				this->vertexArray.push_back(Vector3D((float)j / (float)resolution - 0.5f, 0.5f, (float)i / (float)resolution - 0.5f));
				this->normalArray.push_back(Vector3D(0.0f, 1.0f, 0.0f));
				break;
			case 3:
				this->vertexArray.push_back(Vector3D((float)i / (float)resolution - 0.5f, -0.5f, (float)j / (float)resolution - 0.5f));
				this->normalArray.push_back(Vector3D(0.0f, -1.0f, 0.0f));
				break;
			case 4:
				this->vertexArray.push_back(Vector3D((float)i / (float)resolution - 0.5f, (float)j / (float)resolution - 0.5f, 0.5f));
				this->normalArray.push_back(Vector3D(0.0f, 0.0f, 1.0f));
				break;
			case 5:
				this->vertexArray.push_back(Vector3D((float)j / (float)resolution - 0.5f, (float)i / (float)resolution - 0.5f, -0.5f));
				this->normalArray.push_back(Vector3D(0.0f, 0.0f, -1.0f));
				break;
			}
		}
	}

	for (i = 0; i < resolution; ++i)
	{
		for (j = 0; j < resolution; ++j)
		{
			polygon = new PolygonElement;

			polygon->vertexIndices.push_back(offset + (resolution + 1) * (i + 1) + (j + 0));
			polygon->vertexIndices.push_back(offset + (resolution + 1) * (i + 1) + (j + 1));
			polygon->vertexIndices.push_back(offset + (resolution + 1) * (i + 0) + (j + 1));
			polygon->vertexIndices.push_back(offset + (resolution + 1) * (i + 0) + (j + 0));

			polygon->normalIndices.push_back(offset + (resolution + 1) * (i + 1) + (j + 0));
			polygon->normalIndices.push_back(offset + (resolution + 1) * (i + 1) + (j + 1));
			polygon->normalIndices.push_back(offset + (resolution + 1) * (i + 0) + (j + 1));
			polygon->normalIndices.push_back(offset + (resolution + 1) * (i + 0) + (j + 0));

			this->polygons.push_back(polygon);
		}
	}
}

#endif

#if 0

void OBJModel::genQuad(const Vector3D &t, const Vector3D &b, double disp, unsigned int resolution)
{
	PolygonElement *polygon;

	Vector3D n;

	unsigned int offset;
	unsigned int i, j;


	offset = (unsigned int)this->vertexArray.size();

	n = t.cross(b).normalize();

	for (i = 0; i < resolution + 1; ++i)
	{
		for (j = 0; j < resolution + 1; ++j)
		{
			this->vertexArray.push_back(n * disp + t * ((double)i / (double)resolution - 0.5) + b * ((double)j / (double)resolution - 0.5));
			this->normalArray.push_back(n);
		}
	}

	for (i = 0; i < resolution; ++i)
	{
		for (j = 0; j < resolution; ++j)
		{
			polygon = new PolygonElement;

			polygon->vertexIndices.push_back(offset + (resolution + 1) * (i + 1) + (j + 0));
			polygon->vertexIndices.push_back(offset + (resolution + 1) * (i + 1) + (j + 1));
			polygon->vertexIndices.push_back(offset + (resolution + 1) * (i + 0) + (j + 1));
			polygon->vertexIndices.push_back(offset + (resolution + 1) * (i + 0) + (j + 0));

			polygon->normalIndices.push_back(offset + (resolution + 1) * (i + 1) + (j + 0));
			polygon->normalIndices.push_back(offset + (resolution + 1) * (i + 1) + (j + 1));
			polygon->normalIndices.push_back(offset + (resolution + 1) * (i + 0) + (j + 1));
			polygon->normalIndices.push_back(offset + (resolution + 1) * (i + 0) + (j + 0));

			this->polygons.push_back(polygon);
		}
	}
}

#endif

void OBJModel::genQuad(const Vector3D v[4], unsigned int resolution)
{
	this->genQuad(v, resolution, resolution);
}

void OBJModel::genQuad(const Vector3D v[4], const Vector3D n[4], unsigned int resolution)
{
	this->genQuad(v, n, resolution, resolution);
}

void OBJModel::genQuad(const Vector3D v[4], unsigned int resY, unsigned int resX)
{
	Vector3D n[4];

	unsigned int i;


	for (i = 0; i < 4; ++i)
	{
		n[i] = (v[i] - v[(i + 3) % 4]).cross(v[(i + 1) % 4] - v[i]);
	}

	this->genQuad(v, n, resY, resX);
}

void OBJModel::genQuad(const Vector3D v[4], const Vector3D n[4], unsigned int resX, unsigned int resY)
{
	PolygonElement *polygon;

	double s, t;
	double w[4];

	Vector3D normal;
	double   length;

	unsigned int offset;
	unsigned int i, j;


	offset = (unsigned int)this->vertexArray.size();

	for (i = 0; i < resY + 1; ++i)
	{
		s = (double)i / (double)resY;

		for (j = 0; j < resX + 1; ++j)
		{
			t = (double)j / (double)resX;

			w[0] = (1 - s) * (1 - t);
			w[1] = (1 - s) * (t);
			w[2] = (s)* (t);
			w[3] = (s)* (1 - t);

			this->vertexArray.push_back(v[0] * w[0] + v[1] * w[1] + v[2] * w[2] + v[3] * w[3]);
			normal = n[0] * w[0] + n[1] * w[1] + n[2] * w[2] + n[3] * w[3];
			length = normal.length();
			if (length > Vector3D::epsilon)
			{
				normal /= length;
			}
			else
			{
				normal.setZero();
			}
			this->normalArray.push_back(normal);
		}
	}

	for (i = 0; i < resY; ++i)
	{
		for (j = 0; j < resX; ++j)
		{
			polygon = new PolygonElement;

			polygon->vertexIndices.push_back(offset + (resX + 1) * (i + 0) + (j + 0));
			polygon->vertexIndices.push_back(offset + (resX + 1) * (i + 0) + (j + 1));
			polygon->vertexIndices.push_back(offset + (resX + 1) * (i + 1) + (j + 1));
			polygon->vertexIndices.push_back(offset + (resX + 1) * (i + 1) + (j + 0));

			polygon->normalIndices.push_back(offset + (resX + 1) * (i + 0) + (j + 0));
			polygon->normalIndices.push_back(offset + (resX + 1) * (i + 0) + (j + 1));
			polygon->normalIndices.push_back(offset + (resX + 1) * (i + 1) + (j + 1));
			polygon->normalIndices.push_back(offset + (resX + 1) * (i + 1) + (j + 0));

			this->polygons.push_back(polygon);
		}
	}
}

void OBJModel::genTriangle(const Vector3D v[3], unsigned int resolution)
{
	Vector3D n[3];

	unsigned int i;


	for (i = 0; i < 3; ++i)
	{
		n[i] = (v[i] - v[(i + 2) % 3]).cross(v[(i + 1) % 3] - v[i]).normalize();
	}

	this->genTriangle(v, n, resolution);
}

void OBJModel::genTriangle(const Vector3D v[3], const Vector3D n[3], unsigned int resolution)
{
	PolygonElement *polygon;

	Vector3D v1, v2, n1, n2;

	double s, t;

	unsigned int offset;
	unsigned int c1, c2;
	unsigned int i, j;


	offset = (unsigned int)this->vertexArray.size();

	for (i = 0; i <= resolution; ++i)
	{
		s = (double)i / (double)resolution;

		v1 = v[0] * (1.0 - s) + v[1] * s;
		v2 = v[0] * (1.0 - s) + v[2] * s;
		n1 = n[0] * (1.0 - s) + n[1] * s;
		n2 = n[0] * (1.0 - s) + n[2] * s;

		for (j = 0; j <= i; ++j)
		{
			t = (i == 0 ? 0 : (double)j / (double)i);

			this->vertexArray.push_back(v1 * (1.0 - t) + v2 * t);
			this->normalArray.push_back(n1 * (1.0 - t) + n2 * t);
		}
	}

	c2 = 0;

	for (i = 0; i < resolution; ++i)
	{
		c1 = c2;
		c2 += (i + 1);

		for (j = 0; j < 2 * i + 1; ++j)
		{
			polygon = new PolygonElement;

			switch (j % 2)
			{
			case 0:

				polygon->vertexIndices.push_back(offset + c1 + (j / 2));
				polygon->vertexIndices.push_back(offset + c2 + (j / 2));
				polygon->vertexIndices.push_back(offset + c2 + (j / 2) + 1);

				polygon->normalIndices.push_back(offset + c1 + (j / 2));
				polygon->normalIndices.push_back(offset + c2 + (j / 2));
				polygon->normalIndices.push_back(offset + c2 + (j / 2) + 1);

				break;

			case 1:

				polygon->vertexIndices.push_back(offset + c1 + (j / 2));
				polygon->vertexIndices.push_back(offset + c2 + (j / 2) + 1);
				polygon->vertexIndices.push_back(offset + c1 + (j / 2) + 1);

				polygon->normalIndices.push_back(offset + c1 + (j / 2));
				polygon->normalIndices.push_back(offset + c2 + (j / 2) + 1);
				polygon->normalIndices.push_back(offset + c1 + (j / 2) + 1);

				break;
			}

			this->polygons.push_back(polygon);
		}
	}
}

void OBJModel::genBox(const unsigned int resolution[3])
{
	//this->genQuad(Vector3D(1.0, 0.0, 0.0), Vector3D(0.0, 1.0, 0.0), 0.5, resolution);
	//this->genQuad(Vector3D(0.0, 1.0, 0.0), Vector3D(0.0, 0.0, 1.0), 0.5, resolution);
	//this->genQuad(Vector3D(0.0, 0.0, 1.0), Vector3D(1.0, 0.0, 0.0), 0.5, resolution);

	//this->genQuad(Vector3D(0.0, 1.0, 0.0), Vector3D(1.0, 0.0, 0.0), 0.5, resolution);
	//this->genQuad(Vector3D(0.0, 0.0, 1.0), Vector3D(0.0, 1.0, 0.0), 0.5, resolution);
	//this->genQuad(Vector3D(1.0, 0.0, 0.0), Vector3D(0.0, 0.0, 1.0), 0.5, resolution);

	Vector3D quads[][4] =
	{
		{ Vector3D(-0.5, -0.5, 0.5), Vector3D(0.5, -0.5, 0.5), Vector3D(0.5, 0.5, 0.5), Vector3D(-0.5, 0.5, 0.5), },
		{ Vector3D(0.5, -0.5, 0.5), Vector3D(0.5, -0.5, -0.5), Vector3D(0.5, 0.5, -0.5), Vector3D(0.5, 0.5, 0.5), },
		{ Vector3D(0.5, -0.5, -0.5), Vector3D(-0.5, -0.5, -0.5), Vector3D(-0.5, 0.5, -0.5), Vector3D(0.5, 0.5, -0.5), },
		{ Vector3D(-0.5, -0.5, -0.5), Vector3D(-0.5, -0.5, 0.5), Vector3D(-0.5, 0.5, 0.5), Vector3D(-0.5, 0.5, -0.5), },
		{ Vector3D(-0.5, 0.5, 0.5), Vector3D(0.5, 0.5, 0.5), Vector3D(0.5, 0.5, -0.5), Vector3D(-0.5, 0.5, -0.5), },
		{ Vector3D(-0.5, -0.5, -0.5), Vector3D(0.5, -0.5, -0.5), Vector3D(0.5, -0.5, 0.5), Vector3D(-0.5, -0.5, 0.5), },
	};

	this->genQuad(quads[0], quads[0], resolution[1], resolution[0]);
	this->genQuad(quads[1], quads[1], resolution[2], resolution[0]);
	this->genQuad(quads[2], quads[2], resolution[1], resolution[0]);
	this->genQuad(quads[3], quads[3], resolution[2], resolution[0]);
	this->genQuad(quads[4], quads[4], resolution[1], resolution[2]);
	this->genQuad(quads[5], quads[5], resolution[1], resolution[2]);

	this->collectCloseVertices();
	//this->normalize();
}

void OBJModel::genSphere(unsigned int resolution)
{
	this->genIcosahedron(resolution);
}

#if 0

void OBJModel::genRhombicDodechedron(unsigned int resolution)
{
	/////////////////////////////////////////////
	// 2d) Rhombic Dodecahedron
	//     - created by adding six points on the face centers of a cube

	const double sin45 = 1.0 / sqrt(2.0);
	const double rhom_dodec_x = 1.0 / sqrt(3.0);
	const double rhom_dodec_y = sqrt(2.0 / 3.0);

#define ALIGN_WITH_OTHER_GEOM

#ifndef ALIGN_WITH_OTHER_GEOM

	const double rhom_dodec_basePt[][3] =
	{
		// This base mesh can align with XYZ axes

		{ 0.0, 0.0, 1.0 },		// 0: tip

		{ rhom_dodec_x, -rhom_dodec_x, rhom_dodec_x },	// 1
		{ rhom_dodec_x, rhom_dodec_x, rhom_dodec_x },	// 2
		{ -rhom_dodec_x, rhom_dodec_x, rhom_dodec_x },	// 3
		{ -rhom_dodec_x, -rhom_dodec_x, rhom_dodec_x },	// 4

		{ 1.0, 0.0, 0.0 },
		{ 0.0, 1.0, 0.0 },
		{ -1.0, 0.0, 0.0 },
		{ 0.0, -1.0, 0.0 },

		{ rhom_dodec_x, -rhom_dodec_x, -rhom_dodec_x },
		{ rhom_dodec_x, rhom_dodec_x, -rhom_dodec_x },
		{ -rhom_dodec_x, rhom_dodec_x, -rhom_dodec_x },
		{ -rhom_dodec_x, -rhom_dodec_x, -rhom_dodec_x },

		{ 0.0, 0.0, -1.0 },
	};

#else

	const double rhom_dodec_basePt[][3] =
	{
		// This base mesh can align its 1st base polygon with that of
		// the icosa. and cube for easier comparison

		{ 0.0, sin45, sin45 },				// 0: tip

		{ -rhom_dodec_x, 0.0, rhom_dodec_y },	// 1
		{ rhom_dodec_x, 0.0, rhom_dodec_y },	// 2
		{ rhom_dodec_x, rhom_dodec_y, 0.0 },	// 3
		{ -rhom_dodec_x, rhom_dodec_y, 0.0 },	// 4

		{ 0.0, -sin45, sin45 },				// 5
		{ 1.0, 0.0, 0.0 },				// 6
		{ 0.0, sin45, -sin45 },				// 7
		{ -1.0, 0.0, 0.0 },				// 8

		{ -rhom_dodec_x, -rhom_dodec_y, 0.0 },	// 9
		{ rhom_dodec_x, -rhom_dodec_y, 0.0 },	// 10
		{ rhom_dodec_x, 0.0, -rhom_dodec_y },	// 11
		{ -rhom_dodec_x, 0.0, -rhom_dodec_y },	// 12

		{ 0.0, -sin45, -sin45 },				// 13
	};

#endif

	const unsigned int rhom_dodec_baseIndex[][4] =
	{
		{ 1, 5, 2, 0 },
		{ 2, 6, 3, 0 },
		{ 3, 7, 4, 0 },
		{ 4, 8, 1, 0 },

		{ 8, 9, 5, 1 },
		{ 5, 10, 6, 2 },
		{ 6, 11, 7, 3 },
		{ 7, 12, 8, 4 },

		{ 9, 13, 10, 5 },
		{ 10, 13, 11, 6 },
		{ 11, 13, 12, 7 },
		{ 12, 13, 9, 8 },
	};

	const unsigned int rhom_dodec_numBaseQuad = 12;
	const unsigned int rhom_dodec_numBaseVertex = 14;


	unsigned int i, j;

	Vector3D v[4];


	for (i = 0; i < rhom_dodec_numBaseQuad; ++i)
	{
		for (j = 0; j < 4; ++j)
		{
			v[j] = Vector3D(rhom_dodec_basePt[rhom_dodec_baseIndex[i][j]]);
		}

		this->genQuad(v, resolution);
	}

	for (i = 0; i < this->vertexArray.size(); ++i)
	{
		this->normalArray[i] = this->vertexArray[i] = this->vertexArray[i].normalize();
	}

	//this->collectCloseVertices();
}

#endif

//////////////////////////////////////////////////

static void
build_cycleMap_RhombicDodechedron(int cycle_map[14][14])
{
	int cycle0_frm[6] = { 0, 2, 6, 11, 7, 4 };
	int cycle0_to[6] = { 1, 5, 10, 13, 12, 8 };
	int cycle1_frm[6] = { 0, 3, 7, 12, 8, 1 };
	int cycle1_to[6] = { 2, 6, 11, 13, 9, 5 };
	int cycle2_frm[6] = { 0, 2, 5, 9, 8, 4 };
	int cycle2_to[6] = { 3, 6, 10, 13, 12, 7 };
	int cycle3_frm[6] = { 0, 1, 5, 10, 6, 3 };
	int cycle3_to[6] = { 4, 8, 9, 13, 11, 7 };

	int i, j;

	// unuse -> -1
	for (i = 0; i < 14; i++)
		for (j = 0; j < 14; j++)
			cycle_map[i][j] = -1;

	// cycle ID: 0 to 3
	for (i = 0; i < 6; i++)
	{
		cycle_map[cycle0_frm[i]][cycle0_to[i]] = 0;
		cycle_map[cycle0_to[i]][cycle0_frm[i]] = 0;
	}
	for (i = 0; i < 6; i++)
	{
		cycle_map[cycle1_frm[i]][cycle1_to[i]] = 1;
		cycle_map[cycle1_to[i]][cycle1_frm[i]] = 1;
	}
	for (i = 0; i < 6; i++)
	{
		cycle_map[cycle2_frm[i]][cycle2_to[i]] = 2;
		cycle_map[cycle2_to[i]][cycle2_frm[i]] = 2;
	}
	for (i = 0; i < 6; i++)
	{
		cycle_map[cycle3_frm[i]][cycle3_to[i]] = 3;
		cycle_map[cycle3_to[i]][cycle3_frm[i]] = 3;
	}
}


void OBJModel::genRhombicDodechedron(const unsigned int resolution[4])
{
	/////////////////////////////////////////////
	// 2d) Rhombic Dodecahedron
	//     - created by adding six points on the face centers of a cube

	const double sin45 = 1.0 / sqrt(2.0);
	const double rhom_dodec_x = 1.0 / sqrt(3.0);
	const double rhom_dodec_y = sqrt(2.0 / 3.0);

#define ALIGN_WITH_OTHER_GEOM
#ifndef ALIGN_WITH_OTHER_GEOM

	const double rhom_dodec_basePt[][3] =
	{
		// This base mesh can align with XYZ axes

		{ 0.0, 0.0, 1.0 },		// 0: tip

		{ rhom_dodec_x, -rhom_dodec_x, rhom_dodec_x },	// 1
		{ rhom_dodec_x, rhom_dodec_x, rhom_dodec_x },	// 2
		{ -rhom_dodec_x, rhom_dodec_x, rhom_dodec_x },	// 3
		{ -rhom_dodec_x, -rhom_dodec_x, rhom_dodec_x },	// 4

		{ 1.0, 0.0, 0.0 },
		{ 0.0, 1.0, 0.0 },
		{ -1.0, 0.0, 0.0 },
		{ 0.0, -1.0, 0.0 },

		{ rhom_dodec_x, -rhom_dodec_x, -rhom_dodec_x },
		{ rhom_dodec_x, rhom_dodec_x, -rhom_dodec_x },
		{ -rhom_dodec_x, rhom_dodec_x, -rhom_dodec_x },
		{ -rhom_dodec_x, -rhom_dodec_x, -rhom_dodec_x },

		{ 0.0, 0.0, -1.0 },
	};

#else

	const double rhom_dodec_basePt[][3] =
	{
		// This base mesh can align its 1st base polygon with that of
		// the icosa. and cube for easier comparison

		{ 0.0, sin45, sin45 },				// 0: tip

		{ -rhom_dodec_x, 0.0, rhom_dodec_y },	// 1
		{ rhom_dodec_x, 0.0, rhom_dodec_y },	// 2
		{ rhom_dodec_x, rhom_dodec_y, 0.0 },	// 3
		{ -rhom_dodec_x, rhom_dodec_y, 0.0 },	// 4

		{ 0.0, -sin45, sin45 },				// 5
		{ 1.0, 0.0, 0.0 },				// 6
		{ 0.0, sin45, -sin45 },				// 7
		{ -1.0, 0.0, 0.0 },				// 8

		{ -rhom_dodec_x, -rhom_dodec_y, 0.0 },	// 9
		{ rhom_dodec_x, -rhom_dodec_y, 0.0 },	// 10
		{ rhom_dodec_x, 0.0, -rhom_dodec_y },	// 11
		{ -rhom_dodec_x, 0.0, -rhom_dodec_y },	// 12

		{ 0.0, -sin45, -sin45 },				// 13
	};

#endif

	const unsigned int rhom_dodec_baseIndex[][4] =
	{
		{ 1, 5, 2, 0 },
		{ 2, 6, 3, 0 },
		{ 3, 7, 4, 0 },
		{ 4, 8, 1, 0 },

		{ 8, 9, 5, 1 },
		{ 5, 10, 6, 2 },
		{ 6, 11, 7, 3 },
		{ 7, 12, 8, 4 },

		{ 9, 13, 10, 5 },
		{ 10, 13, 11, 6 },
		{ 11, 13, 12, 7 },
		{ 12, 13, 9, 8 },
	};

	const unsigned int rhom_dodec_numBaseQuad = 12;
	const unsigned int rhom_dodec_numBaseVertex = 14;

	unsigned int i, j;

	// Build: cycle map matrix
	int cycle_map[14][14];
	build_cycleMap_RhombicDodechedron(cycle_map);

	// Generate quads
	Vector3D v[4];
	for (i = 0; i < rhom_dodec_numBaseQuad; ++i)
	{
		int cycle1, cycle2;

		for (j = 0; j < 4; ++j)
		{
			v[j] = Vector3D(rhom_dodec_basePt[rhom_dodec_baseIndex[i][j]]);
		}

		cycle1 = cycle_map[rhom_dodec_baseIndex[i][0]][rhom_dodec_baseIndex[i][1]];
		cycle2 = cycle_map[rhom_dodec_baseIndex[i][1]][rhom_dodec_baseIndex[i][2]];

		this->genQuad(v, resolution[cycle1], resolution[cycle2]);
	}

	for (i = 0; i < this->vertexArray.size(); ++i)
	{
		this->normalArray[i] = this->vertexArray[i] = this->vertexArray[i].normalize();
	}

	this->collectCloseVertices();
}

//////////////////////////////////////////////////

void OBJModel::genIcosahedron(unsigned int resolution)
{
	enum INDEX
	{
		ZA, ZB, ZC, ZD,
		YA, YB, YC, YD,
		XA, XB, XC, XD,
	};

	const double icos_tau = 0.8506508084;		/* t=(1+sqrt(5))/2, tau=t/sqrt(1+t^2)  */
	const double icos_one = 0.5257311121;		/* one=1/sqrt(1+t^2) , unit sphere     */

	const double icosahedron_basePt[][3] =
	{
		{ icos_tau, icos_one, 0 },        // ZA
		{ -icos_tau, icos_one, 0 },        // ZB
		{ -icos_tau, -icos_one, 0 },        // ZC
		{ icos_tau, -icos_one, 0 },        // ZD
		{ icos_one, 0, icos_tau },        // YA
		{ icos_one, 0, -icos_tau },        // YB
		{ -icos_one, 0, -icos_tau },        // YC
		{ -icos_one, 0, icos_tau },        // YD
		{ 0, icos_tau, icos_one },        // XA
		{ 0, -icos_tau, icos_one },        // XB
		{ 0, -icos_tau, -icos_one },        // XC
		{ 0, icos_tau, -icos_one },        // XD
	};

#define ICOSA_ORDER 0

#ifdef ICOSA_ORDER

#if (ICOSA_ORDER == 0) // Given order in original source code

	const unsigned int icosahedron_baseIndex[][3] =
	{
		{ YA, XA, YD },		// 0
		{ YA, YD, XB },		// 1
		{ YB, YC, XD },		// 2
		{ YB, XC, YC },		// 3
		{ ZA, YA, ZD },		// 4
		{ ZA, ZD, YB },		// 5
		{ ZC, YD, ZB },		// 6
		{ ZC, ZB, YC },		// 7
		{ XA, ZA, XD },		// 8
		{ XA, XD, ZB },		// 9
		{ XB, XC, ZD },		// 10
		{ XB, ZC, XC },		// 11
		{ XA, YA, ZA },		// 12
		{ XD, ZA, YB },		// 13
		{ YA, XB, ZD },		// 14
		{ YB, ZD, XC },		// 15
		{ YD, XA, ZB },		// 16
		{ YC, ZB, XD },		// 17
		{ YD, ZC, XB },		// 18
		{ YC, XC, ZC }		// 19
	};

#elif (ICOSA_ORDER == 1) // Triangle strip order

	const unsigned int icosahedron_baseIndex[][3] =
	{
		{ YA, XA, YD },		// 0
		{ YA, YD, XB },		// 1
		{ YD, ZC, XB },		// 18
		{ ZC, YD, ZB },		// 6
		{ YD, XA, ZB },		// 16
		{ XA, XD, ZB },		// 9
		{ XA, ZA, XD },		// 8
		{ XA, YA, ZA },		// 12
		{ ZA, YA, ZD },		// 4
		{ YA, XB, ZD },		// 14
		{ XB, XC, ZD },		// 10
		{ XB, ZC, XC },		// 11
		{ YC, XC, ZC },		// 19
		{ ZC, ZB, YC },		// 7
		{ YC, ZB, XD },		// 17
		{ YB, YC, XD },		// 2
		{ XD, ZA, YB },		// 13
		{ ZA, ZD, YB },		// 5
		{ YB, XC, YC },		// 3
		{ YB, ZD, XC }		// 15
	};

#elif (ICOSA_ORDER == 2) // "Proxmity to zero" order

	const unsigned int icosahedron_baseIndex[][3] =
	{
		{ YA, XA, YD },		// 0
		{ YA, YD, XB },		// 1
		{ XA, YA, ZA },		// 12
		{ YD, XA, ZB },		// 16
		{ YA, XB, ZD },		// 14
		{ ZA, YA, ZD },		// 4
		{ XA, ZA, XD },		// 8
		{ XA, XD, ZB },		// 9
		{ ZC, YD, ZB },		// 6
		{ YD, ZC, XB },		// 18
		{ XB, ZC, XC },		// 11
		{ XB, XC, ZD },		// 10
		{ ZA, ZD, YB },		// 5
		{ XD, ZA, YB },		// 13
		{ YC, ZB, XD },		// 17
		{ ZC, ZB, YC },		// 7
		{ YB, ZD, XC },		// 15
		{ YB, YC, XD },		// 2
		{ YC, XC, ZC },		// 19
		{ YB, XC, YC }		// 3
	};

#endif

#endif

	const unsigned int icosahedron_baseIndexQuad[][4] =
	{
		{ YA, XA, YD, XB },    // 0
		{ YD, ZB, ZC, XB },    // 1
		{ XA, XD, ZB, YD },    // 2
		{ XA, YA, ZA, XD },    // 3
		{ YA, XB, ZD, ZA },    // 4
		{ XB, ZC, XC, ZD },    // 5
		{ YC, XC, ZC, ZB },    // 6
		{ YC, ZB, XD, YB },    // 7
		{ ZA, ZD, YB, XD },    // 8
		{ YB, ZD, XC, YC },    // 9
	};

	const unsigned int icosahedron_numBaseTri = 20;
	const unsigned int icosahedron_numBaseQuad = 10;
	const unsigned int icosahedron_numBaseVertex = 12;


	unsigned int i, j;

	Vector3D v[4];

#if 0

	for (i = 0; i < icosahedron_numBaseTri; ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			v[j] = Vector3D(icosahedron_basePt[icosahedron_baseIndex[i][j]]);
		}

		this->genTriangle(v, resolution);
	}

	for (i = 0; i < this->vertexArray.size(); ++i)
	{
		this->normalArray[i] = this->vertexArray[i] = this->vertexArray[i].normalize();
	}

#endif

	for (i = 0; i < icosahedron_numBaseQuad; ++i)
	{
		for (j = 0; j < 4; ++j)
		{
			v[j] = Vector3D(icosahedron_basePt[icosahedron_baseIndexQuad[i][j]]);
		}

		this->genQuad(v, resolution);
	}

	for (i = 0; i < this->vertexArray.size(); ++i)
	{
		this->normalArray[i] = this->vertexArray[i] = this->vertexArray[i].normalize();
	}

	this->collectCloseVertices();
}

void OBJModel::genBrick()
{
	PolygonElement *polygon;


	this->vertexArray.push_back(Vector3D(-1.0f, 0.333333f, -0.666667f));
	this->vertexArray.push_back(Vector3D(-1.0f, 0.333333f, 0.0f));
	this->vertexArray.push_back(Vector3D(-0.333333f, 0.333333f, 0.0f));
	this->vertexArray.push_back(Vector3D(-0.333333f, 0.333333f, 0.666667f));
	this->vertexArray.push_back(Vector3D(0.333333f, 0.333333f, 0.666667f));
	this->vertexArray.push_back(Vector3D(0.333333f, 0.333333f, 0.0f));
	this->vertexArray.push_back(Vector3D(1.0f, 0.333333f, 0.0f));
	this->vertexArray.push_back(Vector3D(1.0f, 0.333333f, -0.666667f));
	this->vertexArray.push_back(Vector3D(0.333333f, 0.333333f, -0.666667f));
	this->vertexArray.push_back(Vector3D(-0.333333f, 0.333333f, -0.666667f));

	this->vertexArray.push_back(Vector3D(1.0f, -0.333333f, -0.666667f));
	this->vertexArray.push_back(Vector3D(1.0f, -0.333333f, 0.0f));
	this->vertexArray.push_back(Vector3D(0.333333f, -0.333333f, 0.0f));
	this->vertexArray.push_back(Vector3D(0.333333f, -0.333333f, 0.666667f));
	this->vertexArray.push_back(Vector3D(-0.333333f, -0.333333f, 0.666667f));
	this->vertexArray.push_back(Vector3D(-0.333333f, -0.333333f, 0.0f));
	this->vertexArray.push_back(Vector3D(-1.0f, -0.333333f, 0.0f));
	this->vertexArray.push_back(Vector3D(-1.0f, -0.333333f, -0.666667f));
	this->vertexArray.push_back(Vector3D(-0.333333f, -0.333333f, -0.666667f));
	this->vertexArray.push_back(Vector3D(0.333333f, -0.333333f, -0.666667f));

	polygon = new PolygonElement;

	polygon->vertexIndices.push_back(0);
	polygon->vertexIndices.push_back(1);
	polygon->vertexIndices.push_back(2);
	polygon->vertexIndices.push_back(9);

	this->polygons.push_back(polygon);

	polygon = new PolygonElement;

	polygon->vertexIndices.push_back(2);
	polygon->vertexIndices.push_back(3);
	polygon->vertexIndices.push_back(4);
	polygon->vertexIndices.push_back(5);

	this->polygons.push_back(polygon);

	polygon = new PolygonElement;

	polygon->vertexIndices.push_back(5);
	polygon->vertexIndices.push_back(6);
	polygon->vertexIndices.push_back(7);
	polygon->vertexIndices.push_back(8);


	this->polygons.push_back(polygon);

	polygon = new PolygonElement;

	polygon->vertexIndices.push_back(8);
	polygon->vertexIndices.push_back(9);
	polygon->vertexIndices.push_back(2);
	polygon->vertexIndices.push_back(5);

	this->polygons.push_back(polygon);

	polygon = new PolygonElement;

	polygon->vertexIndices.push_back(10);
	polygon->vertexIndices.push_back(11);
	polygon->vertexIndices.push_back(12);
	polygon->vertexIndices.push_back(19);

	this->polygons.push_back(polygon);

	polygon = new PolygonElement;

	polygon->vertexIndices.push_back(12);
	polygon->vertexIndices.push_back(13);
	polygon->vertexIndices.push_back(14);
	polygon->vertexIndices.push_back(15);

	this->polygons.push_back(polygon);

	polygon = new PolygonElement;

	polygon->vertexIndices.push_back(15);
	polygon->vertexIndices.push_back(16);
	polygon->vertexIndices.push_back(17);
	polygon->vertexIndices.push_back(18);

	this->polygons.push_back(polygon);

	polygon = new PolygonElement;

	polygon->vertexIndices.push_back(18);
	polygon->vertexIndices.push_back(19);
	polygon->vertexIndices.push_back(12);
	polygon->vertexIndices.push_back(15);

	this->polygons.push_back(polygon);

	polygon = new PolygonElement;

	polygon->vertexIndices.push_back(1);
	polygon->vertexIndices.push_back(0);
	polygon->vertexIndices.push_back(17);
	polygon->vertexIndices.push_back(16);

	this->polygons.push_back(polygon);

	polygon = new PolygonElement;

	polygon->vertexIndices.push_back(2);
	polygon->vertexIndices.push_back(1);
	polygon->vertexIndices.push_back(16);
	polygon->vertexIndices.push_back(15);

	this->polygons.push_back(polygon);

	polygon = new PolygonElement;

	polygon->vertexIndices.push_back(3);
	polygon->vertexIndices.push_back(2);
	polygon->vertexIndices.push_back(15);
	polygon->vertexIndices.push_back(14);

	this->polygons.push_back(polygon);

	polygon = new PolygonElement;

	polygon->vertexIndices.push_back(4);
	polygon->vertexIndices.push_back(3);
	polygon->vertexIndices.push_back(14);
	polygon->vertexIndices.push_back(13);

	this->polygons.push_back(polygon);

	polygon = new PolygonElement;

	polygon->vertexIndices.push_back(5);
	polygon->vertexIndices.push_back(4);
	polygon->vertexIndices.push_back(13);
	polygon->vertexIndices.push_back(12);

	this->polygons.push_back(polygon);

	polygon = new PolygonElement;

	polygon->vertexIndices.push_back(6);
	polygon->vertexIndices.push_back(5);
	polygon->vertexIndices.push_back(12);
	polygon->vertexIndices.push_back(11);

	this->polygons.push_back(polygon);

	polygon = new PolygonElement;

	polygon->vertexIndices.push_back(7);
	polygon->vertexIndices.push_back(6);
	polygon->vertexIndices.push_back(11);
	polygon->vertexIndices.push_back(10);

	this->polygons.push_back(polygon);

	polygon = new PolygonElement;

	polygon->vertexIndices.push_back(8);
	polygon->vertexIndices.push_back(7);
	polygon->vertexIndices.push_back(10);
	polygon->vertexIndices.push_back(19);

	this->polygons.push_back(polygon);

	polygon = new PolygonElement;

	polygon->vertexIndices.push_back(9);
	polygon->vertexIndices.push_back(8);
	polygon->vertexIndices.push_back(19);
	polygon->vertexIndices.push_back(18);

	this->polygons.push_back(polygon);

	polygon = new PolygonElement;

	polygon->vertexIndices.push_back(0);
	polygon->vertexIndices.push_back(9);
	polygon->vertexIndices.push_back(18);
	polygon->vertexIndices.push_back(17);

	this->polygons.push_back(polygon);
}

void OBJModel::tranform(const Transformation& transformation)
{
	Vector3D transformedPoint;
	Vector3D transformedNormal;

	unsigned int i;


	for (i = 0; i < this->vertexArray.size(); ++i)
	{
		transformation.execute(this->vertexArray[i], transformedPoint, this->normalArray[i]);
		this->vertexArray[i] = transformedPoint;
	}
}
