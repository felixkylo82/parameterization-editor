#include "vector3D.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <cfloat>
#include <cstring>



//#ifdef FLT_EPSILON
//const double Vector3D::epsilon = FLT_EPSILON;
//#else
//const double Vector3D::epsilon = 1.2e-7;
//#endif

const double Vector3D::epsilon = 1e-10;

#ifdef M_PI
const double Vector3D::pi = M_PI;
#else
const double Vector3D::pi = 3.1415926535897932384626433832795;
#endif


std::ostream &operator<<(std::ostream &os, const Vector3D &v)
{
    return (os << '(' << v.data()[0] << ", " << v.data()[1] << ", " << v.data()[2] << ')');
}


const Vector3D Vector3D::zero  = Vector3D(0.0, 0.0, 0.0);
const Vector3D Vector3D::xAxis = Vector3D(1.0, 0.0, 0.0);
const Vector3D Vector3D::yAxis = Vector3D(0.0, 1.0, 0.0);
const Vector3D Vector3D::zAxis = Vector3D(0.0, 0.0, 1.0);


double Vector3D::distance(const Vector3D &v1, const Vector3D &v2)
{
    return (v1 - v2).length();
}

Vector3D Vector3D::computeNormal(const Vector3D &v1, const Vector3D &v2, const Vector3D &v3, bool normalize)
{
    Vector3D w1 = v1 - v3;
    Vector3D w2 = v2 - v1;

    if(normalize)
    {
        double d1 = w1.length();
        double d2 = w2.length();

        if(d1 <= Vector3D::epsilon && d2 <= Vector3D::epsilon)
        {
            throw "Degenerated triangle";
        }

        w1 /= d1;
        w2 /= d2;

        Vector3D w3 = w1.cross(w2);
        double   d3 = w3.length();

        if(d3 <= Vector3D::epsilon)
        {
            throw "Degenerated triangle";
        }

        return (w3 / d3);
    }
    else
    {
        return w1.cross(w2);
    }
}

Vector3D Vector3D::minimum(const Vector3D &v1, const Vector3D &v2)
{
    return Vector3D(
        Vector3D::minimum(v1.coordinate[0], v2.coordinate[0]),
        Vector3D::minimum(v1.coordinate[1], v2.coordinate[1]),
        Vector3D::minimum(v1.coordinate[2], v2.coordinate[2]));
}

Vector3D Vector3D::maximum(const Vector3D &v1, const Vector3D &v2)
{
    return Vector3D(
        Vector3D::maximum(v1.coordinate[0], v2.coordinate[0]),
        Vector3D::maximum(v1.coordinate[1], v2.coordinate[1]),
        Vector3D::maximum(v1.coordinate[2], v2.coordinate[2]));
}

double Vector3D::minimum(double x1, double x2)
{
    return (x1 <= x2 ? x1 : x2);
}

double Vector3D::maximum(double x1, double x2)
{
    return (x1 >= x2 ? x1 : x2);
}

double Vector3D::area (const Vector3D &v1, const Vector3D &v2, const Vector3D &v3)
{
    return (v2 - v1).cross(v3 - v2).length();
}

double Vector3D::angle(const Vector3D &v1, const Vector3D &v2, const Vector3D &v3)
{
    return (v2 - v1).angle(v3 - v2);
}

Vector3D::Vector3D()
{
    this->setZero();
}

Vector3D::Vector3D(const Vector3D &v)
{
    memcpy(this->coordinate, v.coordinate, 3 * sizeof(double));
}

Vector3D::Vector3D(double x, double y, double z)
{
    this->coordinate[0] = x;
    this->coordinate[1] = y;
    this->coordinate[2] = z;
}

Vector3D::Vector3D(float x, float y, float z)
{
    this->coordinate[0] = x;
    this->coordinate[1] = y;
    this->coordinate[2] = z;
}

Vector3D::Vector3D(const double newCoordinate[3])
{
    memcpy(this->coordinate, newCoordinate, 3 * sizeof(double));
}

Vector3D::Vector3D(const float newCoordinate[3])
{
    this->coordinate[0] = newCoordinate[0];
    this->coordinate[1] = newCoordinate[1];
    this->coordinate[2] = newCoordinate[2];
}

const double &Vector3D::operator[](int index) const
{
    return this->coordinate[index];
}

double &Vector3D::operator[](int index)
{
    return this->coordinate[index];
}

void Vector3D::toArray(double newCoordinate[]) const
{
    memcpy(newCoordinate, this->coordinate, 3 * sizeof(double));
}

void Vector3D::toArray(float newCoordinate[]) const
{
    newCoordinate[0] = (float)this->coordinate[0];
    newCoordinate[1] = (float)this->coordinate[1];
    newCoordinate[2] = (float)this->coordinate[2];
}

const double *Vector3D::data() const
{
    return this->coordinate;
}

double Vector3D::length() const
{
    return sqrt(this->dot(*this));
}

void Vector3D::setZero()
{
    memcpy(this->coordinate, zero.coordinate, 3 * sizeof(double));
}

Vector3D Vector3D::normalize() const
{
    return *this / this->length();
}

Vector3D Vector3D::roundOff() const
{
    return Vector3D(
        floor(this->coordinate[0] + 0.5 - Vector3D::epsilon),
        floor(this->coordinate[1] + 0.5 - Vector3D::epsilon),
        floor(this->coordinate[2] + 0.5 - Vector3D::epsilon));
}

Vector3D Vector3D::operator+(const Vector3D &v) const
{
    return Vector3D(
        this->coordinate[0] + v.coordinate[0],
        this->coordinate[1] + v.coordinate[1],
        this->coordinate[2] + v.coordinate[2]);
}

Vector3D Vector3D::operator-(const Vector3D &v) const
{
    return Vector3D(
        this->coordinate[0] - v.coordinate[0],
        this->coordinate[1] - v.coordinate[1],
        this->coordinate[2] - v.coordinate[2]);
}

Vector3D Vector3D::operator-() const
{
    Vector3D v;

    v.coordinate[0] = -this->coordinate[0];
    v.coordinate[1] = -this->coordinate[1];
    v.coordinate[2] = -this->coordinate[2];

    return v;
}

Vector3D Vector3D::operator*(double scalar) const
{
    return Vector3D(
        this->coordinate[0] * scalar,
        this->coordinate[1] * scalar,
        this->coordinate[2] * scalar);
}

Vector3D Vector3D::operator/(double scalar) const
{
    if(fabs(scalar) <= Vector3D::epsilon)
    {
        std::cerr << "Vector operation exception occured!" << std::endl;

        throw "Vector divided by zero!";
    }

    return Vector3D(
        this->coordinate[0] / scalar,
        this->coordinate[1] / scalar,
        this->coordinate[2] / scalar);
}

Vector3D Vector3D::cross(const Vector3D &v) const
{
    return Vector3D(
        this->coordinate[1] * v.coordinate[2] - this->coordinate[2] * v.coordinate[1],
        this->coordinate[2] * v.coordinate[0] - this->coordinate[0] * v.coordinate[2],
        this->coordinate[0] * v.coordinate[1] - this->coordinate[1] * v.coordinate[0]);
}

double Vector3D::dot(const Vector3D &v) const
{
    return
        this->coordinate[0] * v.coordinate[0] +
        this->coordinate[1] * v.coordinate[1] +
        this->coordinate[2] * v.coordinate[2] ;
}

double Vector3D::angle(const Vector3D &v) const
{
    double cosineAngle;


    //cosineAngle = this->dot(v) / this->length() / v.length();

    cosineAngle = this->normalize().dot(v.normalize());

    if(cosineAngle >= 1.0)
    {
        return 0.0;
    }
    else if(cosineAngle <= -1.0)
    {
        return Vector3D::pi;
    }
    else
    {
        return acos(cosineAngle);
    }
}

Vector3D Vector3D::rotate(const Vector3D &axis, double angle) const
{
    double c, s;
    double mat[9];
    double n1, n2, n3;

    //Vector3D outVec;


    // form the rotation matrix
    c = cos(angle);
    s = sin(angle);

    n1 = axis[0];
    n2 = axis[1];
    n3 = axis[2];

    mat[0] = c + n1*n1*(1-c);
    mat[1] = n1*n2*(1-c) - s*n3;
    mat[2] = n3*n1*(1-c) + s*n2;

    mat[3] = n1*n2*(1-c) + s*n3;
    mat[4] = c + n2*n2*(1-c);
    mat[5] = n3*n2*(1-c) - s*n1;

    mat[6] = n3*n1*(1-c) - s*n2;
    mat[7] = n3*n2*(1-c) + s*n1;
    mat[8] = c + n3*n3*(1-c);

    // perform rotation
    //outVec[0] = this->dot(mat + 0);
    //outVec[1] = this->dot(mat + 3);
    //outVec[2] = this->dot(mat + 6);


    //return outVec;

    return Vector3D(
        this->dot(mat + 0),
        this->dot(mat + 3),
        this->dot(mat + 6));
}

void Vector3D::swap(Vector3D &v)
{
    Vector3D swp(v);

    v     = *this;
    *this = swp  ;
}

Vector3D &Vector3D::operator=(const Vector3D &v)
{
    memcpy(this->coordinate, v.coordinate, 3 * sizeof(double));

    return *this;
}

Vector3D &Vector3D::operator+=(const Vector3D &v)
{
    this->coordinate[0] += v.coordinate[0];
    this->coordinate[1] += v.coordinate[1];
    this->coordinate[2] += v.coordinate[2];

    return *this;
}

Vector3D &Vector3D::operator-=(const Vector3D &v)
{
    this->coordinate[0] -= v.coordinate[0];
    this->coordinate[1] -= v.coordinate[1];
    this->coordinate[2] -= v.coordinate[2];

    return *this;
}

Vector3D &Vector3D::operator*=(double scalar)
{
    this->coordinate[0] *= scalar;
    this->coordinate[1] *= scalar;
    this->coordinate[2] *= scalar;

    return *this;
}

Vector3D &Vector3D::operator/=(double scalar)
{
    if(fabs(scalar) <= Vector3D::epsilon)
    {
        std::cerr << "Vector operation exception occured!" << std::endl;

        throw "Vector divided by zero!";
    }

    this->coordinate[0] /= scalar;
    this->coordinate[1] /= scalar;
    this->coordinate[2] /= scalar;

    return *this;
}

bool Vector3D::operator==(const Vector3D &v2) const
{
    return
        this->coordinate[0] == v2.coordinate[0] &&
        this->coordinate[1] == v2.coordinate[1] &&
        this->coordinate[2] == v2.coordinate[2] ;
}

bool Vector3D::operator!=(const Vector3D &v2) const
{
    return !this->operator ==(v2);
}

bool Vector3D::operator<(const Vector3D &v2) const
{
    unsigned int i;

    for(i = 0; i < 3; ++i)
    {
        if(this->coordinate[i] < v2.coordinate[i])
        {
            return true;
        }

        if(this->coordinate[i] > v2.coordinate[i])
        {
            return false;
        }
    }

    return false;
}

bool Vector3D::operator<=(const Vector3D &v2) const
{
    unsigned int i;

    for(i = 0; i < 3; ++i)
    {
        if(this->coordinate[i] < v2.coordinate[i])
        {
            return true;
        }

        if(this->coordinate[i] > v2.coordinate[i])
        {
            return false;
        }
    }

    return true;
}

bool Vector3D::operator>(const Vector3D &v2) const
{
    return !(*this <= v2);
}

bool Vector3D::operator>=(const Vector3D &v2) const
{
    return !(*this < v2);
}

bool Vector3D::isZero() const
{
    return *this == Vector3D::zero;
}
