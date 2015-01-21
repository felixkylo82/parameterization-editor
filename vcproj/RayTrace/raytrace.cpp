#include "raytrace.h"

#include "vector3D.h"

#include <cfloat>
#include <cmath>

namespace RayTrace
{
	const static double rayEpsilon = Vector3D::epsilon;

	RayTrace::IntersectionType rayTracePoint(const Vector3D &point, const Vector3D &start, const Vector3D &direction, Vector3D &hitPoint, double &minDist, bool twoSide)
    {
        Vector3D v;

        double m, dist;


        v    = (point - start);
        dist = v.length();

        if(dist <= rayEpsilon)
        {
            hitPoint = point;
            //minDist  = 0.0;
            if(dist < minDist) minDist  = dist;

            return RayTrace::point;
        }

        //m  = v.dot(direction);

        //if(fabs(1.0 - m  / dist) <= rayEpsilon)
        //if(v.cross(direction).length() / dist <= rayEpsilon)
        if(v.cross(direction).length() <= rayEpsilon)
        {
            m  = v.dot(direction);

            if(twoSide)
            {
                m = fabs(m);
            }

            if(m >= 0.0)
            {
                if(m < minDist)
                {
                    hitPoint = point;
                    minDist  = m;
                }

                return RayTrace::point;
            }
        }


        return RayTrace::none;
    }

    RayTrace::IntersectionType rayTraceLineSegment(const Vector3D &endPoint1, const Vector3D &endPoint2, const Vector3D &start, const Vector3D &direction, Vector3D &hitPoint, double &minDist, bool twoSide)
    {
        double a1, a2, a;
        double d1, d2, e;

        double m;

        IntersectionType t1, t2;


        t1 = rayTracePoint(endPoint1, start, direction, hitPoint, minDist); 
        t2 = rayTracePoint(endPoint2, start, direction, hitPoint, minDist); 

        if(RayTrace::point == t1 || RayTrace::point == t2)
        {
            return RayTrace::point;
        }

        Vector3D v1, v2;
        //Vector3D w1, w2;
        Vector3D w2;
        Vector3D u;
        Vector3D n1, n2, n;
        Vector3D p, v;


        v1 = endPoint1 - start;
        v2 = endPoint2 - start;

        d1 = v1.length();
        d2 = v2.length();

        if(d1 <= rayEpsilon || d2 <= rayEpsilon)
        {
            std::cerr << "hit points ??" << std::endl;
            exit(1);
        }

        if(!twoSide)
        {
            if( direction.dot(v1) < 0.0 &&
                direction.dot(v2) < 0.0)
            {
                return RayTrace::none;
            }
        }

        w2 = endPoint2 - endPoint1;

        d2 = w2.length();

        if(d2 <= rayEpsilon)
        {
            //std::cerr << "Degenerated line segment!" << std::endl;
            return RayTrace::none;
        }

        w2 /= d2;

        u = w2.cross(direction);

        e = u.length();

        if(e <= rayEpsilon || fabs(u.dot(v1) / e) > rayEpsilon)
        {
            return RayTrace::none;
        }

        n1 = v1.cross(direction);
        n2 = direction.cross(v2);

        if(n1.dot(n2) >= 0.0)
        {
            a1 = n1.length();
            a2 = n2.length();

            a = a1 + a2;

            if(a <= rayEpsilon)
            {
                e = -w2.dot(v1);

                if(e > 0.0 && e < d2)
                {
                    p = w2 * e + endPoint1;
                }
                else
                {
                    return RayTrace::none;
                }
            }
            else
            {
                p = (endPoint1 * a2 + endPoint2 * a1) / a;
            }

            v = p - start;
            m = v.dot(direction);

            if(twoSide)
            {
                m = fabs(m);
            }

            if(m >= 0.0)
            {
                if(m < minDist)
                {
                    hitPoint = p;
                    minDist  = m;
                }

                return RayTrace::point;
            }
        }


        return RayTrace::none;
    }

    bool insideTriangle(
        const Vector3D &vertex1,
        const Vector3D &vertex2,
        const Vector3D &vertex3,
        const Vector3D &testPoint)
    {
        Vector3D hitPoint;


        Vector3D v1 = vertex2 - vertex1;
        Vector3D v2 = vertex3 - vertex2;
        Vector3D v3 = vertex1 - vertex3;

#if 1
        double d1 = v1.length();
        double d2 = v2.length();
        double d3 = v3.length();

        if(d1 <= rayEpsilon || d2 <= rayEpsilon || d3 <= rayEpsilon)
        {
            fprintf(stderr, "A unhandled case found [degenerated triangle]!\n");
            return false;
        }

        v1 /= d1; v2 /= d2; v3 /= d3;
#endif

        Vector3D w1 = testPoint - vertex1;
        Vector3D w2 = testPoint - vertex2;
        Vector3D w3 = testPoint - vertex3;

#if 0
        double e1 = w1.length();
        double e2 = w2.length();
        double e3 = w3.length();

        if(e1 <= rayEpsilon || e2 <= rayEpsilon || e3 <= rayEpsilon)
        {
            return true;
        }

        w1 /= e1; w2 /= e2; w3 /= e3;
#endif

        Vector3D normal = v1.cross(v2);

#if 1
        double len = normal.length();

        if(len <= rayEpsilon)
        {
            fprintf(stderr, "A unhandled case found [degenerated triangle]!\n");
            return false;
        }

        normal /= len;
#endif

        return (
            v1.cross(w1).dot(normal) >= 0.0 &&
            v2.cross(w2).dot(normal) >= 0.0 &&
            v3.cross(w3).dot(normal) >= 0.0 );
    }

    RayTrace::IntersectionType rayTracePlane(const Vector3D &normal, double offset,
        const Vector3D &start, const Vector3D &direction, Vector3D &hitPoint, double &minDist, bool twoSide)
    {
        double r, s, t;


        r = offset - normal.dot(start);
        s = normal.dot(direction);

        if(r == 0.0)
        {
            hitPoint = start;
            minDist  = 0.0;
            return RayTrace::point;
        }

        if(s <= rayEpsilon)
        {
            return RayTrace::none;
        }

        t = r / s;

        if(twoSide || t >= 0)
        {
            hitPoint = start + direction * t;
            minDist = fabs(t);
            return RayTrace::point;
        }

        return RayTrace::none;
    }

    RayTrace::IntersectionType rayTraceTriangle(
        const Vector3D &vertex1,
        const Vector3D &vertex2,
        const Vector3D &vertex3,
        const Vector3D &start,
        const Vector3D &direction,
        Vector3D &hitPoint,
        double &minDist,
        bool twoSide,
        bool straightlyInside)
    {
        Vector3D r, myHitPoint;

        double myMinDist;

        double d1, d2, d3;
        double m, n;

        RayTrace::IntersectionType t1, t2, t3;


        myMinDist = minDist;
        
        t1 = rayTraceLineSegment(vertex1, vertex2, start, direction, myHitPoint, myMinDist, twoSide);
        t2 = rayTraceLineSegment(vertex2, vertex3, start, direction, myHitPoint, myMinDist, twoSide);
        t3 = rayTraceLineSegment(vertex3, vertex1, start, direction, myHitPoint, myMinDist, twoSide);

        if(t1 || t2 || t3)
        {
            if(straightlyInside)
            {
                return RayTrace::none;
            }
            else
            {
                if(myMinDist < minDist)
                {
                    hitPoint = myHitPoint;
                    minDist  = myMinDist ;
                }
                return RayTrace::point;
            }
        }

        if(twoSide)
        {
            r = vertex1 - start;
        }
        else
        {
            // facing the triangle
            r = vertex1 - start;
            d1 = direction.dot(r);
            if (d1 <= 0.0)
            {
                r = vertex2 - start;
                d2 = direction.dot(r);
                if (d2 <= 0.0)
                {
                    r = vertex3 - start;
                    d3 = direction.dot(r);
                    if (d3 <= 0.0) 
                    {
                        return RayTrace::none;
                    }
                }
            }
        }

        // alignment the normal (two faced triangle)

        Vector3D v1 = vertex2 - vertex1;
        Vector3D v2 = vertex3 - vertex2;
        Vector3D v3 = vertex1 - vertex3;

        d1 = v1.length();
        d2 = v2.length();
        d3 = v3.length();

        if(d1 <= rayEpsilon || d2 <= rayEpsilon || d3 <= rayEpsilon)
        {
            //std::cerr << "Degenerated triangle!" << std::endl;
            return RayTrace::none;
        }

        v1 /= d1;
        v2 /= d2;

        Vector3D hitNormal = v1.cross(v2);

        d1 = hitNormal.length();

        if(d1 <= rayEpsilon)
        {
            //std::cerr << "Degenerated triangle!" << std::endl;
            return RayTrace::none;
        }

        hitNormal /= d1;

        n = direction.dot(hitNormal);

        if(fabs(n) <= rayEpsilon)
        {
            return RayTrace::none;
        }

        // hit the triangle?
        m = hitNormal.dot(r) / n;

        if(m < 0.0 && !twoSide)
        {
            return RayTrace::none;
        }

        /* check if hitNormal need to reverse in direction */
        //if (n > 0.0)
        //{
        //    hitNormal = -hitNormal;
        //}

        Vector3D p = start + direction * m;

        /* check if hitCoord hit the boundary or it is inside the triangle*/

        Vector3D s1 = p - vertex1;
        Vector3D s2 = p - vertex2;
        Vector3D s3 = p - vertex3;

        Vector3D y1 = s1.cross(v1);
        Vector3D y2 = s2.cross(v2);
        Vector3D y3 = s3.cross(v3);

        if( y1.dot(y2) > 0.0 && y2.dot(y3) > 0.0 && y3.dot(y1) > 0.0)
        {
            m = fabs(m);

            if(m < minDist)
            {
                hitPoint = p;
                minDist  = m;
            }

            return RayTrace::point;
        }


        return RayTrace::none;
    }

    RayTrace::IntersectionType planeIntersectPoint(const Vector3D &point,
        const Vector3D &planeNormal, double planeOffset)
    {
        if(fabs(point.dot(planeNormal) - planeOffset) <= rayEpsilon)
        {
            return RayTrace::point;
        }
        else
        {
            return RayTrace::none;
        }
    }

    RayTrace::IntersectionType planeIntersectLineSegment(const Vector3D &endPoint1, const Vector3D &endPoint2,
        const Vector3D &planeNormal, double planeOffset, Vector3D &hitPoint)
    {
        double offset1, offset2;

        IntersectionType t1, t2;


        t1 = planeIntersectPoint(endPoint1, planeNormal, planeOffset);
        t2 = planeIntersectPoint(endPoint2, planeNormal, planeOffset);

        if(RayTrace::point == t1)
        {
            if(RayTrace::point == t2)
            {
                //printf("Intersected a line segment\n");

                return RayTrace::lineSegment;
            }

            hitPoint = endPoint1;

            return RayTrace::point;
        }
        else if(RayTrace::point == t2)
        {
            hitPoint = endPoint2;

            return RayTrace::point;
        }

        offset1 = endPoint1.dot(planeNormal) - planeOffset;
        offset2 = endPoint2.dot(planeNormal) - planeOffset;

        if(offset1 * offset2 <= 0.0)
        {
            offset1 = fabs(offset1);
            offset2 = fabs(offset2);

            if(offset1 + offset2 > rayEpsilon)
            {
                hitPoint = (endPoint1 * offset2 + endPoint2 * offset1) / (offset1 + offset2);
                return RayTrace::point;
            }
            else
            {
                std::cerr << "critical error\n" << std::endl;
                exit(1);
            }
        }


        return none;
    }

    void planeIntersectTriangle(
        const Vector3D &endPoint1, const Vector3D &endPoint2, const Vector3D &endPoint3,
        const Vector3D &planeNormal, double planeOffset,
        std::vector<Vector3D> &frontTriangles, std::vector<Vector3D> &backTriangles)
    {
        Vector3D triangle[3] = {endPoint1, endPoint2, endPoint3, };
        RayTrace::planeIntersectTriangle(triangle, planeNormal, planeOffset, frontTriangles, backTriangles);
    }

    void planeIntersectTriangle(
        const Vector3D triangle[],
        const Vector3D &planeNormal, double planeOffset,
        std::vector<Vector3D> &frontTriangles,
        std::vector<Vector3D> &backTriangles)
    {
        Vector3D hitPoints[3];
        Vector3D v[4], p[4];

        double d[3], e;

        IntersectionType t[3];

        std::vector<Vector3D> *myTriangles[2];

        double theta1, theta2;
        
        int numVertices;

        int i, j;


        for(i = 0; i < 3; ++i)
        {
            t[i] = planeIntersectLineSegment(triangle[i], triangle[(i + 1) % 3], planeNormal, planeOffset, hitPoints[i]);
        }

        e = 0.0;

        for(i = 0; i < 3; ++i)
        {
            e += d[i] = triangle[i].dot(planeNormal) - planeOffset;
        }

        if(RayTrace::lineSegment == t[0] && RayTrace::lineSegment == t[1] && RayTrace::lineSegment == t[2])
        {
            // Intersected a triangle
            if(e >= 0.0)
            {
                frontTriangles.push_back(triangle[0]);
                frontTriangles.push_back(triangle[1]);
                frontTriangles.push_back(triangle[2]);
            }
            else
            {
                backTriangles.push_back(triangle[0]);
                backTriangles.push_back(triangle[1]);
                backTriangles.push_back(triangle[2]);
            }
            return;
        }

        for(i = 0; i < 3; ++i)
        {
            if(RayTrace::lineSegment == t[i])
            {
                if(d[(i + 2) % 3] >= 0.0)
                {
                    frontTriangles.push_back(triangle[0]);
                    frontTriangles.push_back(triangle[1]);
                    frontTriangles.push_back(triangle[2]);
                }
                else
                {
                    backTriangles.push_back(triangle[0]);
                    backTriangles.push_back(triangle[1]);
                    backTriangles.push_back(triangle[2]);
                }
                
                return;
            }
        }

        if(RayTrace::none == t[0] && RayTrace::none == t[1] && RayTrace::none == t[2])
        {
            if(e >= 0.0)
            {
                frontTriangles.push_back(triangle[0]);
                frontTriangles.push_back(triangle[1]);
                frontTriangles.push_back(triangle[2]);
            }
            else
            {
                backTriangles.push_back(triangle[0]);
                backTriangles.push_back(triangle[1]);
                backTriangles.push_back(triangle[2]);
            }
            return;
        }

        for(i = 0; i < 3; ++i)
        {
            if( RayTrace::point == t[(i + 0) % 3] &&
                RayTrace::point == t[(i + 1) % 3] &&
                Vector3D::distance(hitPoints[(i + 0) % 3], hitPoints[(i + 1) % 3]) > rayEpsilon)
            {
                if(d[(i + 1) % 3] >= 0.0)
                {
                    myTriangles[0] = &frontTriangles;
                    myTriangles[1] = &backTriangles ;
                }
                else
                {
                    myTriangles[0] = &backTriangles ;
                    myTriangles[1] = &frontTriangles;
                }

                if( Vector3D::distance(hitPoints[(i + 0) % 3], triangle [(i + 1) % 3]) > rayEpsilon &&
                    Vector3D::distance(hitPoints[(i + 1) % 3], triangle [(i + 1) % 3]) > rayEpsilon)
                {
                    myTriangles[0]->push_back(hitPoints[(i + 1) % 3]);
                    myTriangles[0]->push_back(hitPoints[(i + 0) % 3]);
                    myTriangles[0]->push_back(triangle [(i + 1) % 3]);
                }

                p[0] = hitPoints[(i + 0) % 3];
                p[1] = hitPoints[(i + 1) % 3];
                p[2] = triangle [(i + 2) % 3];
                p[3] = triangle [(i + 0) % 3];

                numVertices = 4;
                j           = 0;

                while(j < numVertices)
                {
                    v[j] = p[(j + 1) % 4] - p[j];

                    if(v[j].length() <= rayEpsilon)
                    {
                        --numVertices;

                        for(int k = j; k < numVertices; ++k)
                        {
                            p[k] = p[k + 1];
                        }
                    }
                    else
                    {
                        ++j;
                    }
                }

                switch(j)
                {
                case 4:

                    theta1 = v[0].angle(v[1]) + v[2].angle(v[3]);
                    theta2 = v[1].angle(v[2]) + v[3].angle(v[0]);

                    if(theta1 > theta2)
                    {
                        myTriangles[1]->push_back(p[0]);
                        myTriangles[1]->push_back(p[1]);
                        myTriangles[1]->push_back(p[2]);

                        myTriangles[1]->push_back(p[2]);
                        myTriangles[1]->push_back(p[3]);
                        myTriangles[1]->push_back(p[0]);
                    }
                    else
                    {
                        myTriangles[1]->push_back(p[1]);
                        myTriangles[1]->push_back(p[2]);
                        myTriangles[1]->push_back(p[3]);

                        myTriangles[1]->push_back(p[3]);
                        myTriangles[1]->push_back(p[0]);
                        myTriangles[1]->push_back(p[1]);
                    }

                    break;

                //case 3:
                default:

                    //std::cout << "a degenerated quad" << std::endl;

                    myTriangles[1]->push_back(p[0]);
                    myTriangles[1]->push_back(p[1]);
                    myTriangles[1]->push_back(p[2]);

                    break;

                //default:

                //    std::cerr << "a unhandled case found!" << std::endl;
                //    std::cout << "case = " << j << std::endl;
                //    std::cout << t[0] << '\t' << t[1] << '\t' << t[2] << std::endl;

                    break;
                }

                return;
            }
        }

        //for(i = 0; i < 3; ++i)
        {
            //if(RayTrace::point == t[i])
            {
                //if(d[(i + 2) % 3] >= 0.0)
                if(e >= 0.0)
                {
                    frontTriangles.push_back(triangle[0]);
                    frontTriangles.push_back(triangle[1]);
                    frontTriangles.push_back(triangle[2]);
                }
                else
                {
                    backTriangles.push_back(triangle[0]);
                    backTriangles.push_back(triangle[1]);
                    backTriangles.push_back(triangle[2]);
                }
                return;
            }
        }

        //fprintf(stderr, "a unhandled case found!\n");
    }

    bool threePlanesIntesection(
        const Vector3D &normal1, double offset1,
        const Vector3D &normal2, double offset2,
        const Vector3D &normal3, double offset3,
        Vector3D &hitPoint)
    {
        const double *A[3] =
        {
            normal1.data(), normal2.data(), normal3.data(),
        };
        const double y[3] = {offset1, offset2, offset3, };
        
        double B[3][3];
        double det;

        int i, j;


        for(i = 0; i < 3; ++i)
        {
            for(j = 0; j < 3; ++j)
            {
                B[j][i] = (
                    A[(i + 1) % 3][(j + 1) % 3] * A[(i + 2) % 3][(j + 2) % 3] -
                    A[(i + 2) % 3][(j + 1) % 3] * A[(i + 1) % 3][(j + 2) % 3] );
            }
        }

        det = A[0][0] * B[0][0] + A[0][1] * B[1][0] + A[0][2] * B[2][0];

        if(det == 0)
        {
            return false;
        }

        for(i = 0; i < 3; ++i)
        {
            for(j = 0; j < 3; ++j)
            {
                B[i][j] /= det;
            }
        }

        for(i = 0; i < 3; ++i)
        {
            hitPoint[i] = B[i][0] * y[0] + B[i][1] * y[1] + B[i][2] * y[2];
        }

        return true;
    }

    bool triangleIntersectTriangle(
        const Vector3D &vertex11, const Vector3D &vertex12, const Vector3D &vertex13,
        const Vector3D &vertex21, const Vector3D &vertex22, const Vector3D &vertex23,
        double threshold, bool verbose)
    {
        Vector3D x[2][3] =
        {
            {vertex11, vertex12, vertex13, },
            {vertex21, vertex22, vertex23, },
        };
        
        Vector3D p[2][3] =
        {
            {vertex11 - vertex21, vertex12 - vertex21, vertex13 - vertex21, },
            {vertex21 - vertex11, vertex22 - vertex11, vertex23 - vertex11, },
        };
        
        Vector3D v[2][3] =
        {
            {vertex12 - vertex11, vertex13 - vertex12, vertex11 - vertex13, },
            {vertex22 - vertex21, vertex23 - vertex22, vertex21 - vertex23, },
        };
        
        //Vector3D n[2] =
        //{
        //    v[0][0].cross(v[0][1]), v[1][0].cross(v[1][1]),
        //};

        Vector3D n[2];

        Vector3D u[2][3];

        double d[2][3], e, f, g, min;

        int count;

        int h, i, j, k;


        for(i = 0; i < 2; ++i)
        {
            for(j = 0; j < 3; ++j)
            {
                e = v[i][j].length();

                if(e > rayEpsilon)
                {
                    u[i][j] = v[i][j] / e;
                }
            }
        }

        //for(i = 0; i < 2; ++i)
        //{
        //    e = n[i].length();

        //    if(e <= rayEpsilon)
        //    {
        //        n[i].setZero();
        //    }
        //    else
        //    {
        //        n[i] /= e;
        //    }
        //}

        for(i = 0; i < 2; ++i)
        {
            n[i] = v[i][0].cross(v[i][1]);
            e    = n[i].length();
            if(e > rayEpsilon)
            {
                n[i] /= e;
            }
            else
            {
                n[i].setZero();
            }
        }

        h = 1;

        for(i = 0; i < 2; ++i)
        {
            for(j = 0; j < 3; ++j)
            {
                d[i][j] = n[h].dot(p[i][j]);
            }

            h = i;
        }

        min   = FLT_MAX;
        count = 0;

        h = 1;

        for(i = 0; i < 2; ++i)
        {
            j = 2;

            for(k = 0; k < 3; ++k)
            {
                if(d[i][j] * d[i][k] < 0.0)
                {
                    f = fabs(d[i][j]);
                    g = fabs(d[i][k]);
                    e = f + g;
                    
                    if(e > rayEpsilon)
                    {
                        if(f > threshold && g > threshold)
                        {
                            Vector3D q = x[i][j] + v[i][j] * (f / e);

                            Vector3D w[3] =
                            {
                                q - x[h][0], q - x[h][1], q - x[h][2], 
                            };

                            //for(l = 0; l < 3; ++l)
                            //{
                            //    if(w[l].length() <= rayEpsilon)
                            //    {
                            //        w[l].setZero();
                            //    }
                            //}

                            double disp1 = u[h][0].cross(w[0]).dot(n[h]);
                            double disp2 = u[h][1].cross(w[1]).dot(n[h]);
                            double disp3 = u[h][2].cross(w[2]).dot(n[h]);

                            if(disp1 > threshold && disp2 > threshold && disp3 > threshold)
                            {
                                min = std::fmin(std::fmin(disp1, disp2), std::fmin(min, disp3));
                                min = std::fmin(min, std::fmin(f, g));

                                if(2 == ++count)
                                {
                                    if(verbose)
                                    {
                                        printf("MIN = %f\n", min);
                                    }
                                    return true;
                                }
                            }
                        }
                    }
                }

                j = k;
            }

            h = i;
        }


        return false;
    }

    IntersectionType lineSegmentIntersectPoint(const Vector3D &point,
        const Vector3D &endPoint1, const Vector3D &endPoint2, Vector3D &hitPoint)
    {
        Vector3D v1, v2, w, p;//n;

        double d1, d2, e;


        v1 = point - endPoint1; 
        v2 = point - endPoint2;

        d1 = v1.length();
        d2 = v2.length();

        if(d1 <= rayEpsilon)
        {
            hitPoint = endPoint1;
            return RayTrace::point;
        }
        
        if(d2 <= rayEpsilon)
        {
            hitPoint = endPoint2;
            return RayTrace::point;
        }

        if(v1.dot(v2) <= 0.0)
        {
            //v1 /= d1;
            //v2 /= d2;

            w = endPoint2 - endPoint1;
            e = w.length();

            if(e > rayEpsilon)
            {
                w /= e;
                p  = endPoint1 + w * v1.dot(w);
                if((p - point).length() <= rayEpsilon)
                {
                    hitPoint = p;
                    return RayTrace::point;
                }
            }
#if 0
            n = v1.cross (v2);
            e = n .length(  );

            if(e <= rayEpsilon)
            {
                hitPoint = point; // ..
                return RayTrace::point;
            }
#endif
        }


        return RayTrace::none;
    }

    IntersectionType lineSegmentIntersectLineSegment(const Vector3D &endPoint11, const Vector3D &endPoint12,
        const Vector3D &endPoint21, const Vector3D &endPoint22, Vector3D &hitPoint)
    {
        IntersectionType t1, t2;

        Vector3D v1, v2, w2;

        Vector3D myHitPoint;

        double d1, d2;

        double minDist;


        v1 = endPoint12 - endPoint11;
        v2 = endPoint22 - endPoint21;

        d1 = v1.length();
        d2 = v2.length();

        if(d1 <= rayEpsilon)
        {
            // degenerated line segment
            return lineSegmentIntersectPoint((endPoint11 + endPoint12) * 0.5, endPoint21, endPoint22, hitPoint);
        }

        if(d2 <= rayEpsilon)
        {
            // degenerated line segment
            return lineSegmentIntersectPoint((endPoint21 + endPoint22) * 0.5, endPoint11, endPoint12, hitPoint);
        }

        t1 = lineSegmentIntersectPoint(endPoint11, endPoint21, endPoint22, hitPoint);
        t2 = lineSegmentIntersectPoint(endPoint12, endPoint21, endPoint22, hitPoint);

        if(t1 == RayTrace::point && t2 == RayTrace::point)
        {
            std::cerr << "Intersected a line segment" << std::endl;
            return RayTrace::lineSegment;
        }
        else if(t1 == RayTrace::point || t2 == RayTrace::point)
        {
            return RayTrace::point;
        }

        w2 = v2 / d2;

        minDist = DBL_MAX;

        if(RayTrace::none != RayTrace::rayTraceLineSegment(endPoint11, endPoint12, endPoint21, w2, myHitPoint, minDist, false))
        {
            if(minDist < d2)
            {
                hitPoint = myHitPoint;
                return RayTrace::point;
            }
        }

        return RayTrace::none;
    }
}
