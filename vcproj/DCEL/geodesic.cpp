#define _USE_MATH_DEFINES
#include "geodesic.h"

#include "common.h"

#include <cmath>

#include <windows.h>
#include <GL/gl.h>

Geodesic *Geodesic::compute(const DoublyConnectedEdgeList &dcel, const Vector3D &from, const Vector3D &to, float epsilon)
{
    DoublyConnectedEdgeList::Primitive *prim  [2];
    DoublyConnectedEdgeList::Vertex    *vertex[2];
    DoublyConnectedEdgeList::Edge      *edge  [2];
    DoublyConnectedEdgeList::Face      *face  [2];

    Vector3D point[2];

    Geodesic *geodesic;

    double previousLength, currentLength;

    Geodesic::iterator giter;

    double startTime, endTime;

    unsigned int i;


    geodesic = NULL;

    startTime = myTime();

    face[0] = dcel.findNearestFace(from, point[0], vertex[0], edge[0]);
    face[1] = dcel.findNearestFace(to  , point[1], vertex[1], edge[1]);

    endTime = myTime();

    Geodesic::timeForSearchingNearestFace += endTime - startTime;

    if(NULL != vertex[0] && NULL != vertex[1])
    {
        geodesic = new Geodesic;

        if(NULL != geodesic)
        {
            //if(vertex[0] == vertex[1])
            //{
            //    geodesic->push_back(DoublyConnectedEdgeList::SearchPathNode(vertex[0], point[0]));
            //}
            //else if(NULL != edge[0] && NULL != edge[1] && edge[0] == edge[1])
            //{
            //    geodesic->push_back(DoublyConnectedEdgeList::SearchPathNode(edge[0], point[0]));
            //    geodesic->push_back(DoublyConnectedEdgeList::SearchPathNode(edge[1], point[1]));
            //}
            //else if(NULL != face[0] && NULL != face[1] && face[0] == face[1])
            //{
            //    geodesic->push_back(DoublyConnectedEdgeList::SearchPathNode(face[0], point[0]));
            //    geodesic->push_back(DoublyConnectedEdgeList::SearchPathNode(face[1], point[1]));
            //}
            //else
            {
                for(i = 0; i < 2; ++i)
                {
                    prim[i] = face[i];

                    if(NULL == prim[i])
                    {
                        prim[i] = edge[i];
                    }

                    if(NULL == prim[i])
                    {
                        prim[i] = vertex[i];
                    }
                }

                //DoublyConnectedEdgeList::ConstFaceSet fList;

                //dcel.findCommonFaces(prim[0], prim[1], fList);

                //if(!fList.empty())
                if(dcel.haveCommonFaces(prim[0], prim[1]))
                {
                    geodesic->push_back(DoublyConnectedEdgeList::SearchPathNode(prim[0], point[0]));
                    geodesic->push_back(DoublyConnectedEdgeList::SearchPathNode(prim[1], point[1]));
                }
                else
                {
                    dcel.computeShortestPath(vertex[0], vertex[1], *geodesic);

                    for(i = 0; i < 2; ++i)
                    {
                        if(NULL != face[i])
                        {
                            if(0 == i)
                            {
                                geodesic->push_front(DoublyConnectedEdgeList::SearchPathNode(face[i], point[i]));
                            }
                            else
                            {
                                geodesic->push_back(DoublyConnectedEdgeList::SearchPathNode(face[i], point[i]));
                            }
                        }
                        else if(NULL != edge[i])
                        {
                            if(0 == i)
                            {
                                geodesic->push_front(DoublyConnectedEdgeList::SearchPathNode(edge[i], point[i]));
                            }
                            else
                            {
                                geodesic->push_back(DoublyConnectedEdgeList::SearchPathNode(edge[i], point[i]));
                            }
                        }
                    }
                }
            }
        }

        startTime = myTime();

        for(giter = geodesic->begin(); giter != geodesic->end(); ++giter)
        {
            geodesic->mergePoint3(dcel, giter);
        }

#ifdef _DEBUG

        geodesic->validate(dcel);
        //geodesic->print();

#endif

        currentLength = geodesic->length();

        do
        {
            previousLength = currentLength;

            geodesic->optimize(dcel);

            currentLength = geodesic->length();
        }
        while(previousLength - currentLength > epsilon);

        endTime = myTime();

        Geodesic::timeForGeodesicOptimization += endTime - startTime;
    }

    return geodesic;
}

#if 0

bool Geodesic::isOrdered(Geodesic *g1, Geodesic *g2, Geodesic *g3)
{
    if( NULL != g1 &&
        NULL != g2 &&
        NULL != g3 )
    {
        Vector3D d1 = g1->forwardDirection();
        Vector3D d2 = g2->forwardDirection();
        Vector3D d3 = g3->forwardDirection();

        return d1.cross(d2).dot(d2.cross(d3)) > 0.0;
    }
}

#endif

bool Geodesic::validate(const DoublyConnectedEdgeList &dcel) const
{
    //DoublyConnectedEdgeList::ConstFaceSet faces;

    const DoublyConnectedEdgeList::Primitive *p;
    const DoublyConnectedEdgeList::Vertex    *v;

    Geodesic::const_iterator iter, iter2;

    bool ret;


    ret = true;

    try
    {
        for(iter = this->begin(); iter != this->end(); ++iter)
        {
            p = iter->at;
            if(DoublyConnectedEdgeList::Primitive::vertex == p->getPrimitiveType())
            {
                v = (const DoublyConnectedEdgeList::Vertex *)p;

                if(v->getPosition() != iter->position)
                {
                    throw "Inconsistent node position!";
                }
            }
        }

        if(this->size() >= 2)
        {
            //iter = this->begin();
            //iter2 = iter;
            //++iter2;

            //while(iter2 != this->end())
            //{
            //    if(iter->at == iter2->at)
            //    {
            //        throw "Consecutive nodes are located at the same dcel primitive!";
            //    }

            //    iter = iter2;
            //    ++iter2;
            //}

            //iter = this->begin();
            //iter2 = iter;
            //++iter2;

            //while(iter2 != this->end())
            //{
            //    if((iter->position - iter2->position).length() < FLT_EPSILON)
            //    {
            //        throw "Consecutive nodes are located at the same location!";
            //    }

            //    iter = iter2;
            //    ++iter2;
            //}

            iter = this->begin();
            iter2 = iter;
            ++iter2;

            while(iter2 != this->end())
            {
                //faces.clear();
                //dcel.findCommonFaces(iter->at, iter2->at, faces);

                //if(faces.empty())
                if(!dcel.haveCommonFaces(iter->at, iter2->at))
                {
                    throw "Consecutive dcel primitives do not have a common face!";
                }

                iter = iter2;
                ++iter2;
            }

            for(iter = this->begin(); iter != this->end(); ++iter)
            {
                ++(iter2 = iter);

                while(iter2 != this->end())
                {
                    if(iter->at == iter2->at)
                    {
                        throw "Two non-consecutive nodes are located at the same dcel primitive!";
                    }

                    ++iter2;
                }
            }

            for(iter = this->begin(); iter != this->end(); ++iter)
            {
                ++(iter2 = iter);

                if(iter2 != this->end())
                {
                    ++iter2;
                }

                while(iter2 != this->end())
                {
                    //faces.clear();
                    //dcel.findCommonFaces(iter->at, iter2->at, faces);

                    //if(!faces.empty())
                    if(dcel.haveCommonFaces(iter->at, iter2->at))
                    {
                        throw "Non consecutive dcel primitives have a common face";
                    }

                    ++iter2;
                }
            }
        }
    }
    catch(const char err[])
    {
        std::cerr << err << std::endl;

        ret = false;
    }


    return ret;
}

Vector3D Geodesic::findRelativePositionOnSurface(
    const DoublyConnectedEdgeList &dcel,
    const Vector3D &v1, const Vector3D &v2, const Vector3D &v3,
    const Vector3D &relativePosition)
{
    Geodesic *g;

    double a1, a2;

    Vector3D position;


    g = NULL;

    Vector3D dir1 = v2 - v1;
    Vector3D dir2 = v3 - v1;

    a1 = dir1.cross(relativePosition).length();
    a2 = dir2.cross(relativePosition).length();

    if(a1 + a2 <= FLT_EPSILON)
    {
        g = Geodesic::compute(dcel, v1, (v2 + v3) / 2.0);

        position = g->findMidPoint(relativePosition.length());
    }
    else
    {
        g = Geodesic::compute(dcel, v2, v3);

        if(g)
        {
            position = g->findMidPoint(g->length() * (a1 / (a1 + a2)));
            
            delete g;

            g = Geodesic::compute(dcel, v1, position);

            if(g)
            {
                position = g->findMidPoint(relativePosition.length());
            }
        }
    }

    if(g)
    {
        delete g;
    }


    return position;
}

Vector3D Geodesic::findRelativePositionOnSurface(
    const DoublyConnectedEdgeList &dcel,
    const Vector3D &v1, const Vector3D &v2, const Vector3D &v3,
    double r, double theta)
{
    Geodesic *g;

    double a1, a2, x, l;

    Vector3D v4, position;


    g = NULL;

    Vector3D dir1 = v2 - v1;
    Vector3D dir2 = v3 - v1;

    double angle = dir1.angle(dir2);

    a1 = dir1.length() * sin(theta        );
    a2 = dir2.length() * sin(angle - theta);

    if(a1 + a2 <= FLT_EPSILON)
    {
        g = Geodesic::compute(dcel, v1, (v2 + v3) / 2.0);

        l = g->length();

        if(r / l > 1 - FLT_EPSILON)
        {
            r = l * (1 - FLT_EPSILON);
        }

        position = g->findMidPoint(r);
    }
    else
    {
        g = Geodesic::compute(dcel, v2, v3);

        if(g)
        {
            l = g->length();

            x = l * a1 / (a1 + a2);

            if(x / l < FLT_EPSILON)
            {
                x = l * FLT_EPSILON;
            }
            else if(x / l > 1 - FLT_EPSILON)
            {
                x = l * (1 - FLT_EPSILON);
            }

            position = g->findMidPoint(x);
            
            delete g;

            g = Geodesic::compute(dcel, v1, position);

            if(g)
            {
                l = g->length();

                if(l <= 2.0 * FLT_EPSILON)
                {
                    position = v1; // keep unchanged
                }
                else
                {
                    if(r / l > 1 - FLT_EPSILON)
                    {
                        r = l * (1 - FLT_EPSILON);
                    }

                    position = g->findMidPoint(r);
                }
            }
        }
    }

    if(g)
    {
        delete g;
    }


    return position;
}

double Geodesic::timeForSearchingNearestFace = 0.0;
double Geodesic::timeForGeodesicOptimization = 0.0;

double Geodesic::length() const
{
    double sum;


    if(this->size() >= 2)
    {
        sum = Geodesic::length(this->begin(), this->end());
    }
    else
    {
        sum = 0.0;
    }


    return sum;
}

Vector3D Geodesic::findMidPoint(double distance) const
{
    Geodesic::const_iterator current;
    Geodesic::const_iterator next   ;

    Vector3D direction, position;

    double currentLength, m;


    if(!this->empty())
    {
        if(distance <= 0.0)
        {
            distance = 0.0;
        }

        currentLength = 0.0;

        next = this->begin();

        if(next == this->end())
        {
            position = current->position;
        }
        else
        {
            do
            {
                current = next;
                ++next;

                if(next != this->end())
                {
                    direction = next->position - current->position;
                    currentLength += direction.length();
                }
            }
            while(currentLength < distance && next != this->end());

            if(next == this->end())
            {
                position = current->position;
            }
            else
            {
                m = direction.length();

                if(m > FLT_EPSILON)
                {
                    position = next->position - direction * ((currentLength - distance) / m);
                }
                else
                {
                    position = (current->position + next->position) / 2.0;
                }
            }
        }
    }


    return position;
}

Vector3D Geodesic::forwardDirection () const
{
    Geodesic::const_iterator next;

    Vector3D direction;

    double length;


    length = 0.0;
    ++(next = this->begin());

    do
    {
        direction = next->position - this->begin()->position;

        length = direction.length();
    }
    while(length <= FLT_EPSILON && ++next != this->end());

    direction /= length;

    return direction;
}

Vector3D Geodesic::backwardDirection() const
{
    Geodesic::const_reverse_iterator next;

    Vector3D direction;

    double length;


    length = 0.0;
    ++(next = this->rbegin());

    do
    {
        direction = next->position - this->begin()->position;

        length = direction.length();
    }
    while(length <= FLT_EPSILON && ++next != this->rend());

    direction /= length;

    return direction;
}

double Geodesic::length(const Geodesic::const_iterator &start, const Geodesic::const_iterator &end)
{
    double sum;
    Geodesic::const_iterator iter1, iter2;


    sum = 0.0;

    iter1 = start;
    iter2 = iter1;
    ++iter2;

    while(iter2 != end)
    {
        sum += Vector3D::distance(iter1->position, iter2->position);

        iter1 = iter2;
        ++iter2;
    }


    return sum;
}

void Geodesic::optimize(const DoublyConnectedEdgeList &dcel)
{
    Geodesic::iterator iter1;

    bool stable;


    stable = false;

    //while(!stable)
    {
        stable = true;

        for(iter1 = this->begin(); iter1 != this->end(); ++iter1)
        {
            if(this->mergePoint2(dcel, iter1))
            {
                stable = false;
                continue;
            }

            if(this->addPoint(dcel, iter1))
            {
                stable = false;
                continue;
            }

            this->movePoint(dcel, iter1);
        }

        for(iter1 = this->begin(); iter1 != this->end(); ++iter1)
        {
            this->mergePoint3(dcel, iter1);
        }
    }

#ifdef _DEBUG

    this->validate(dcel);

#endif
}

bool Geodesic::addPoint(const DoublyConnectedEdgeList &dcel, Geodesic::iterator &iter1)
{
    std::list<const DoublyConnectedEdgeList::Edge *> edgeList;
    std::list<const DoublyConnectedEdgeList::Edge *>::const_iterator eiter;

    const DoublyConnectedEdgeList::Vertex *centre;

    double originalLength, currentLength, minimumLength;
    double sumOfAngle;

    Vector3D dir[2];

    double r1, r2, R3, r3, v1[2], v2[2], M[4], Minv[4], Mdet, alpha;

    Geodesic geodesic[2];
    int minGeodesicIndex;

    unsigned int i;

    Geodesic::iterator iter2, iter3;


    ++(iter2 = iter1);
    if(iter2 == this->end()) return false;

    ++(iter3 = iter2);
    if(iter3 == this->end()) return false;

    minGeodesicIndex = -1;

    if(DoublyConnectedEdgeList::Primitive::vertex == iter2->at->getPrimitiveType())
    {
        centre = (DoublyConnectedEdgeList::Vertex *)iter2->at;

        minimumLength = originalLength =
            Vector3D::distance(iter1->position, iter2->position) +
            Vector3D::distance(iter2->position, iter3->position) ;

        for(i = 0; i < 2; ++i)
        {
            edgeList.clear();

            switch(i)
            {
            case 0:
                dcel.transverseEdges(centre, iter1->at, iter3->at, edgeList);
                break;

            case 1:
                dcel.transverseEdges(centre, iter3->at, iter1->at, edgeList);
                edgeList.reverse();
                break;
            }

            if(!edgeList.empty())
            {
                sumOfAngle = 0.0;

                eiter  = edgeList.begin();

                dir[1] = (iter1->position - centre->getPosition());

                //DoublyConnectedEdgeList::ConstFaceSet fList;

                //dcel.findCommonFaces(iter1->at, *eiter, fList);

                while(sumOfAngle < M_PI && eiter != edgeList.end())
                {
                    dir[0] = dir[1];
                    dir[1] = ((*eiter)->getTwinEdge()->getOrigin()->getPosition() - centre->getPosition());

                    //sumOfAngle += acos(dir[1].dot(dir[0]));
                    sumOfAngle += dir[1].angle(dir[0]);

                    ++eiter;
                }

                if(sumOfAngle < M_PI)
                {
                    dir[0] = dir[1];
                    dir[1] = (iter3->position - centre->getPosition()).normalize();

                    //sumOfAngle += acos(dir[1].dot(dir[0]));
                    sumOfAngle += dir[1].angle(dir[0]);
                }

                if(sumOfAngle < M_PI)
                {
                    r1 = (iter1->position - centre->getPosition()).length();
                    r2 = (iter3->position - centre->getPosition()).length();

                    v1[0] = r1 ;
                    v1[1] = 0.0;

                    v2[0] = r2 * cos(sumOfAngle);
                    v2[1] = r2 * sin(sumOfAngle);

                    M[0] = v1[0] - v2[0];
                    M[1] = v1[1] - v2[1];

                    geodesic[i].clear();

                    alpha = 0.0;

                    eiter = edgeList.begin();

                    dir[1] = (iter1->position - centre->getPosition()).normalize();

                    while(eiter != edgeList.end())
                    {
                        dir[0]  = dir[1];
                        dir[1]  = ((*eiter)->getTwinEdge()->getOrigin()->getPosition() - centre->getPosition());
                        R3      = dir[1].length();
                        dir[1] /= R3;

                        //alpha += acos(dir[1].dot(dir[0]));
                        alpha += dir[0].angle(dir[1]);

                        M[2] = cos(alpha);
                        M[3] = sin(alpha);

                        Mdet = M[0] * M[3] - M[1] * M[2];

                        //Minv[0] =  M[3] / Mdet;
                        Minv[1] = -M[1] / Mdet;
                        //Minv[2] = -M[2] / Mdet;
                        Minv[3] =  M[0] / Mdet;

                        //v3[0] = v1[0] * (1.0 + M[0] * (Minv[0]));
                        //v3[1] = v1[1] * (1.0 + M[1] * (Minv[2]));

                        //r3 = sqrt(v3[0] * v3[0] + v3[1] * v3[1]);

                        r3 = Minv[1] * v1[0] + Minv[3] * v1[1];

                        if(r3 < R3)
                        {
                            geodesic[i].push_back(DoublyConnectedEdgeList::SearchPathNode(*eiter, centre->getPosition() + dir[1] * r3));
                        }
                        else
                        {
                            geodesic[i].push_back(DoublyConnectedEdgeList::SearchPathNode((*eiter)->getTwinEdge()->getOrigin()));
                        }

                        ++eiter;
                    }

                    geodesic[i].mergePoint3(dcel, geodesic[i].begin());

                    currentLength = geodesic[i].length() +
                        Vector3D::distance(iter1->position, geodesic[i].front().position) +
                        Vector3D::distance(iter3->position, geodesic[i].back ().position) ;

                    if(currentLength < minimumLength)
                    {
                        minGeodesicIndex = i;
                        minimumLength    = currentLength;
                    }
                }
            }
        }

        if(-1 != minGeodesicIndex)
        {
            this->insert(iter2, geodesic[minGeodesicIndex].begin(), geodesic[minGeodesicIndex].end());
            this->erase (iter2);
        }
    }


    return (-1 != minGeodesicIndex);
}

void Geodesic::movePoint(const DoublyConnectedEdgeList &dcel, Geodesic::iterator &iter1)
{
    DoublyConnectedEdgeList::Edge *edge;

    Vector3D dir    ;
    Vector3D v1 , v2;
    Vector3D w1 , w2;

    double d1, d2;
    double e1, e2;
    double r , R ;

    Geodesic::iterator iter2, iter3;


    ++(iter2 = iter1);
    if(iter2 == this->end()) return;

    ++(iter3 = iter2);
    if(iter3 == this->end()) return;

    if(DoublyConnectedEdgeList::Primitive::edge == iter2->at->getPrimitiveType())
    {
        edge = (DoublyConnectedEdgeList::Edge *)iter2->at;

        dir  = edge->getTwinEdge()->getOrigin()->getPosition() - edge->getOrigin()->getPosition();
        R    = dir.length();
        dir /= R;

        v1 = iter1->position - edge->getOrigin()->getPosition();

        d1 = v1.dot(dir);
        w1 = dir * d1;
        e1 = Vector3D::distance(v1, w1);

        v2 = iter3->position - edge->getOrigin()->getPosition();

        d2 = v2.dot(dir);
        w2 = dir * d2;
        e2 = Vector3D::distance(v2, w2);

        if(e1 + e2 > FLT_EPSILON)
        {
            r = (d1 * e2 + d2 * e1) / (e1 + e2);

            if(0 < r && r < R)
            {
                iter2->position = edge->getOrigin()->getPosition() + dir * r;
            }
            else
            {
                if(r <= 0.0)
                {
                    iter2->at       = edge->getOrigin();
                    iter2->position = edge->getOrigin()->getPosition();
                }
                else
                {
                    iter2->at       = edge->getTwinEdge()->getOrigin();
                    iter2->position = edge->getTwinEdge()->getOrigin()->getPosition();
                }

                if( iter2->at == iter1->at ||
                    iter2->at == iter3->at )
                {
                    this->erase(iter2);
                }
            }
        }
    }
}

bool Geodesic::mergePoint2(const DoublyConnectedEdgeList &dcel, Geodesic::iterator &iter1)
{
    const DoublyConnectedEdgeList::Vertex *v;
    Vector3D x;

    DoublyConnectedEdgeList::Edge *e1, *e2, *e3;

    double length1, length2;

    bool merged;

    Geodesic::iterator iter2, iter3, iter4;


    merged = false;

    ++(iter2 = iter1);
    if(iter2 == this->end()) return false;

    ++(iter3 = iter2);
    if(iter3 == this->end()) return false;

#if 0

    if( iter1->at == iter2->at ||
        iter3->at == iter2->at )
    {
        this->erase(iter2);

        merged = true;
    }
    else if((iter1->position - iter2->position).length() <= FLT_EPSILON)
    {
        //x = (iter1->position + iter2->position) / 2.0;
        //v = dcel.findnearestvertex(x);

        v = dcel.findNearestVertex(iter2->position);

        if((v->getPosition() - iter1->position).length() <= FLT_EPSILON)
        {
            iter1->position = v->getPosition();
            iter1->at       = v;

            this->erase(iter2);

            merged = true;
        }
    }
    else if((iter3->position - iter2->position).length() <= FLT_EPSILON)
    {
        //x = (iter3->position + iter2->position) / 2.0;
        //v = dcel.findNearestVertex(x);

        v = dcel.findNearestVertex(iter2->position);

        if((v->getPosition() - iter3->position).length() <= FLT_EPSILON)
        {
            iter3->position = v->getPosition();
            iter3->at       = v;

            this->erase(iter2);

            merged = true;
        }
    }

#endif
    
    if(
        DoublyConnectedEdgeList::Primitive::edge == iter2->at->getPrimitiveType() &&
        DoublyConnectedEdgeList::Primitive::edge == iter3->at->getPrimitiveType() )
    {
        e1 = (DoublyConnectedEdgeList::Edge *)iter2->at;
        e2 = (DoublyConnectedEdgeList::Edge *)iter3->at;

        if( e1->getOrigin() == e2               ->getOrigin() ||
            e1->getOrigin() == e2->getTwinEdge()->getOrigin() )
        {
            v = e1->getOrigin();
        }
        else if(
            e1->getTwinEdge()->getOrigin() == e2               ->getOrigin() ||
            e1->getTwinEdge()->getOrigin() == e2->getTwinEdge()->getOrigin() )
        {
            v = e1->getTwinEdge()->getOrigin();
        }
        else
        {
            v = NULL;
        }

        if(NULL != v)
        {
            ++(iter4 = iter3);

            if(iter4 == this->end()) return false;

            while(iter4 != this->end())
            {
                if(DoublyConnectedEdgeList::Primitive::edge != iter4->at->getPrimitiveType())
                {
                    break;
                }

                e3 = (DoublyConnectedEdgeList::Edge *)iter4->at;

                if( e3               ->getOrigin() != v &&
                    e3->getTwinEdge()->getOrigin() != v )
                {
                    break;
                }

                ++iter4;
            }

            if(iter4 == this->end())
            {
                --iter4;
            }

            --(iter3 = iter4);

            if(iter2 != iter3)
            {
                length1 =
                    (iter1->position - iter2->position).length() +
                    (iter3->position - iter4->position).length() +
                    Geodesic::length(iter2, iter4) ;

                length2 =
                    (iter1->position - v->getPosition()).length() +
                    (iter4->position - v->getPosition()).length() ;

                if(length2 < length1)
                {
                    if(iter1->at == v)
                    {
                        this->erase(iter2, iter4);
                    }
                    else
                    {
                        iter2->at       = v;
                        iter2->position = v->getPosition();

                        this->erase(++iter2, iter4);
                    }

                }
            }
        }
    }


    return merged;
}

void Geodesic::mergePoint3(const DoublyConnectedEdgeList &dcel, Geodesic::iterator &iter1)
{
    Geodesic::iterator iter2, iter3, iter4;


    ++(iter2 = iter1);
    if(iter2 == this->end()) return;

    ++(iter4 = iter2);
    if(iter4 == this->end()) return;

    iter3 = iter2;

    while(iter4 != this->end())
    {
        //DoublyConnectedEdgeList::ConstFaceSet fList;
        //dcel.findCommonFaces(iter1->at, iter4->at, fList);

        //if(!fList.empty())
        if(dcel.haveCommonFaces(iter1->at, iter4->at))
        {
            iter3 = iter4;
        }

        ++iter4;
    }


    if(iter3 != iter2)
    {
        this->erase(iter2, iter3);
    }
}

#if 0

bool Geodesic::mergePoint(const DoublyConnectedEdgeList &dcel, Geodesic::iterator &iter1, Geodesic::iterator &iter2, Geodesic::iterator &iter3)
{
    bool ret;


    //ret  = this->mergePoint2(dcel, iter2, iter3);
    //ret |= this->mergePoint2(dcel, iter1, iter2);

    ret = this->mergePoint2(dcel, iter1, iter2);

    if(!ret)
    {
        ret = this->mergePoint2(dcel, iter2, iter3);
    }

    ++(iter2 = iter1);
    
    iter3 = iter2;
    if(iter2 != this->end())
    {
        ++iter3;
    }


    return ret;
}

#endif

void Geodesic::print()
{
    Geodesic::iterator iter1;


    for(iter1 = this->begin(); iter1 != this->end(); ++iter1)
    {
        switch(iter1->at->getPrimitiveType())
        {
        case DoublyConnectedEdgeList::Primitive::vertex:
                std::cout << "V\t";
                break;

        case DoublyConnectedEdgeList::Primitive::edge:
                std::cout << "E\t";
                break;

        case DoublyConnectedEdgeList::Primitive::face:
                std::cout << "F\t";
                break;
        }

        std::cout << iter1->at << '\t' << iter1->position << '\n';
    }

    std::cout << std::endl;
}

void Geodesic::render()
{
    glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_LINE_BIT);

    glDisable  (GL_LIGHTING);
    glColor3ub (140, 20, 20);
    glLineWidth(2.0f);

    glBegin(GL_LINE_STRIP);

    for(Geodesic::const_iterator piter = this->begin(); piter != this->end(); ++piter)
    {
        glVertex3dv((*piter).position.data());
    }

    glEnd();

    glPopAttrib();
}
