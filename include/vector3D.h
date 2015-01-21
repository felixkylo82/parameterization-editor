#ifndef VECTOR_H
#define VECTOR_H


#include <iostream>


class Vector3D;


std::ostream &operator<<(std::ostream &os, const Vector3D &v);


class Vector3D
{
public:

    static double distance(const Vector3D &v1, const Vector3D &v2);

    static Vector3D computeNormal(const Vector3D &v1, const Vector3D &v2, const Vector3D &v3, bool normalize = false);

    static Vector3D minimum(const Vector3D &v1, const Vector3D &v2);
    static Vector3D maximum(const Vector3D &v1, const Vector3D &v2);

    static double minimum(double x1, double x2);
    static double maximum(double x1, double x2);

    static double area (const Vector3D &v1, const Vector3D &v2, const Vector3D &v3);
    static double angle(const Vector3D &v1, const Vector3D &v2, const Vector3D &v3);

    static const Vector3D zero ;
    static const Vector3D xAxis;
    static const Vector3D yAxis;
    static const Vector3D zAxis;

    static const double epsilon;
    static const double pi     ;


    Vector3D();
    Vector3D(const Vector3D &v);
    Vector3D(double x, double y, double z);
    Vector3D(float  x, float  y, float  z);
    Vector3D(const double newCoordinate[3]);
    Vector3D(const float  newCoordinate[3]);

    const double &operator[](int index) const;
    double       &operator[](int index)      ;

    void toArray(double newCoordinate[]) const;
    void toArray(float  newCoordinate[]) const;

    const double *data() const;

    double length() const;

    void     setZero  ();
    Vector3D normalize() const;
    Vector3D roundOff () const;

    Vector3D operator+(const Vector3D &v) const;
    Vector3D operator-(const Vector3D &v) const;

    Vector3D operator-() const;

    Vector3D operator*(double scalar) const;
    Vector3D operator/(double scalar) const;

    Vector3D cross (const Vector3D &v) const;
    double   dot   (const Vector3D &v) const;
    double   angle (const Vector3D &v) const;
    Vector3D rotate(const Vector3D &axis, double angle) const;

    void     swap (Vector3D &v);

    Vector3D &operator =(const Vector3D &v);
    Vector3D &operator+=(const Vector3D &v);
    Vector3D &operator-=(const Vector3D &v);

    Vector3D &operator*=(double scalar);
    Vector3D &operator/=(double scalar);

    bool operator==(const Vector3D &v2) const;
    bool operator!=(const Vector3D &v2) const;

    bool operator< (const Vector3D &v2) const;
    bool operator<=(const Vector3D &v2) const;

    bool operator> (const Vector3D &v2) const;
    bool operator>=(const Vector3D &v2) const;

    bool isZero() const;

private:

    double coordinate[3];

    friend std::ostream &operator<<(std::ostream &os, const Vector3D &v);
};

#endif // VECTOR_H
