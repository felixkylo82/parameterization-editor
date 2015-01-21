#include "boundingbox.h"

#include <cfloat>


bool Box::rayTraceTwoSide = false;

Box::Box(const Vector3D &min, const Vector3D &max) :

min(min), max(max)

{
}

const Vector3D& Box::getMin() const
{
    return this->min;
}

const Vector3D& Box::getMax() const
{
    return this->max;
}

bool Box::intersectBox(const Box &b) const
{
    bool intercepted;


    if( this->getMax()[0] < b.getMin()[0] ||
        this->getMin()[0] > b.getMax()[0] ||
        this->getMax()[1] < b.getMin()[1] ||
        this->getMin()[1] > b.getMax()[1] ||
        this->getMax()[2] < b.getMin()[2] ||
        this->getMin()[2] > b.getMax()[2] )
    {
        intercepted = false;
    }
    else
    {
        intercepted = true;
    }


    return intercepted;
}

bool Box::intersectRay(const Vector3D &start, const Vector3D &direction) const
{
    double Tnear = -FLT_MAX;
    double Tfar  =  FLT_MAX;

    double T1, T2, tmp;

    bool hit;

    unsigned int i;


    hit = true;

    for(i = 0; hit && i < 3; ++i)
    {
        if(0.0 == direction[i])
        {
            if(start[i] < this->min[i] || start[i] >= this->max[i])
            {
                hit = false;
            }
        }
        else
        {
            T1 = (this->min[i] - start[i]) / direction[i];
            T2 = (this->max[i] - start[i]) / direction[i];

            if(T1 > T2)
            {
                tmp = T1 ;
                T1  = T2 ;
                T2  = tmp;
            }

            if(T1 > Tnear)
            {
                Tnear = T1;
            }

            if(T2 < Tfar)
            {
                Tfar = T2;
            }

            if(Tnear > Tfar)
            {
                hit = false;
            }

            if(!Box::rayTraceTwoSide && Tfar < 0.0)
            {
                hit = false;
            }
        }
    }


    return hit;
}

bool Box::operator<(const Box &box) const
{
    if(this->min < box.min)
    {
        return true;
    }
    else if(this->min != box.min)
    {
        return false;
    }

    if(this->max < box.max)
    {
        return true;
    }
    else if(this->max != box.max)
    {
        return false;
    }


    return false;
}

BoundingBox::BoundingBox(const Vector3D &min, const Vector3D &max) :

Box(min, max), emptyFlag(true)

{
}


BoundingBox::~BoundingBox()
{
}

bool BoundingBox::empty() const
{
    return this->emptyFlag;
}


//void LeafBoundingBox::setRayTraceFunntion(void (*newRayTraceFunction)(void *datum, const Vector3D &start, const Vector3D &direction, void **nearestDatum, Vector3D &hitPoint, double &minDist))
//{
//    LeafBoundingBox::rayTraceFunction = newRayTraceFunction;
//}

unsigned int LeafBoundingBox::maxDataSize = 10;

void (*Box::rayTraceFunction) (void *datum, const Vector3D &start, const Vector3D &direction, void **nearestDatum, Vector3D &hitPoint, double &minDist, bool twoSide) = NULL;

LeafBoundingBox::LeafBoundingBox(const Vector3D &min, const Vector3D &max) :

BoundingBox(min, max), separable(true)

{
}


LeafBoundingBox::~LeafBoundingBox()
{
    this->clear();
}

void LeafBoundingBox::clear()
{
    this->data.clear();
}

bool LeafBoundingBox::putDatum(const Box &box, void *datum)
{
    if(this->intersectBox(box))
    {
        if(datum == NULL)
        {
            this->data.erase(box);
        }
        else
        {
            this->data.insert(std::multimap<Box, void*>::value_type(box, datum));
        }
    }

    return !(this->emptyFlag = this->data.empty());
}

#if 0
void *LeafBoundingBox::getDatum(const Vector3D &position)
{
    std::map<Vector3D, void*>::iterator iter;


    iter = this->data.find(position);

    if(iter == this->data.end())
    {
        return NULL;
    }
    else
    {
        return iter->second;
    }
}
#endif

void LeafBoundingBox::getData(const Box &box, std::list<void *> &data) const
{
    std::multimap<Box, void *>::const_iterator iter;


    if(!this->emptyFlag && this->intersectBox(box))
    {
        for(iter = this->data.begin(); iter != this->data.end(); ++iter)
        {
            if(iter->first.intersectBox(box))
            {
                data.push_back(iter->second);
            }
        }
    }
}

void LeafBoundingBox::rayTrace(const Vector3D &start, const Vector3D &direction, void **nearestDatum, Vector3D &hitPoint, double &minDist)
{
    std::multimap<Box, void *>::const_iterator iter;


    if(!this->emptyFlag && this->intersectRay(start, direction))
    {
        for(iter = this->data.begin(); iter != this->data.end(); ++iter)
        {
            if(iter->first.intersectRay(start, direction))
            {
                if(LeafBoundingBox::rayTraceFunction)
                {
                    LeafBoundingBox::rayTraceFunction(iter->second, start, direction, nearestDatum, hitPoint, minDist, Box::rayTraceTwoSide);
                }
            }
        }
    }
}

bool LeafBoundingBox::isFull() const
{
    return this->data.size() >= LeafBoundingBox::maxDataSize;
}

const std::multimap<Box, void*> &LeafBoundingBox::getData() const
{
    return this->data;
}

bool LeafBoundingBox::getSeparable() const
{
    return this->separable;
}

void LeafBoundingBox::setSeparable(bool separable)
{
    this->separable = separable;
}

//InternalBoundingBox::InternalBoundingBox(const Vector3D &min, const Vector3D &max) :
//
//BoundingBox(min, max), isChildrenBorn(false)
//
//{
//    memset(this->children, NULL, 8 * sizeof(void *));
//}

InternalBoundingBox::InternalBoundingBox(const Vector3D &min, const Vector3D &max, int depth) :

BoundingBox(min, max), isChildrenBorn(false)

{
    memset(this->children, NULL, 8 * sizeof(void *));

    this->produceChildren(depth);
}

InternalBoundingBox::~InternalBoundingBox()
{
    this->clear();
}

void InternalBoundingBox::clear()
{
    unsigned int i;

    for(i = 0; i < 8; ++i)
    {
        delete this->children[i];
        this->children[i]    = NULL;
        this->isChildrenBorn = false;
    }
}


bool InternalBoundingBox::produceChildren(int depth)
{
    unsigned int i;

    Vector3D start, offset;

    bool ret;


    if(isChildrenBorn)
    {
        ret = false;
    }
    else
    {
        if(depth == 1)
        {
            this->produceLeafChildren();
        }
        else if(depth > 1)
        {
            offset = (max - min) / 2.0;

            for(i = 0; i < 8; ++i)
            {
                start = min + Vector3D(((i & 1) == 0) ? 0.0 : offset[0], ((i & 2) == 0) ? 0.0 : offset[1], ((i & 4) == 0) ? 0.0 : offset[2]);

                //this->children[i] = new InternalBoundingBox(start, start + offset);
                //((InternalBoundingBox *)this->children[i])->produceChildren(depth - 1);

                this->children[i] = new InternalBoundingBox(start, start + offset, depth - 1);

                if(NULL == this->children[i])
                {
                    throw "not enough memory";
                }
            }
        }

        this->isChildrenBorn = true;
        ret = true;
    }


    return ret;
}


bool InternalBoundingBox::produceLeafChildren()
{
    unsigned int i;

    Vector3D start, offset;

    bool ret;


    if(isChildrenBorn)
    {
        ret = false;
    }
    else
    {
        offset = (this->max - this->min) / 2.0;

        for(i = 0; i < 8; ++i)
        {
            start = this->min + Vector3D(((i & 1) == 0) ? 0.0 : offset[0], ((i & 2) == 0) ? 0.0 : offset[1], ((i & 4) == 0) ? 0.0 : offset[2]);
            this->children[i] = new LeafBoundingBox(start, start + offset);

            if(NULL == this->children[i])
            {
                throw "not enough memory";
            }
        }

        this->isChildrenBorn = true;
        ret = true;
    }


    return ret;
}

bool InternalBoundingBox::putDatum(const Box &box, void *datum)
{
    return this->putDatum(box, datum, true);
}

bool InternalBoundingBox::putDatum(const Box &box, void *datum, bool enableSplit)
{
    Vector3D avg;

    LeafBoundingBox     *leafBox;
    InternalBoundingBox *newBox ;

    std::multimap<Box, void *>::const_iterator iter;

    bool split;//, intercepted;

    unsigned int index, j;


    if(!this->isChildrenBorn)
    {
        this->produceLeafChildren();
    }

    this->emptyFlag = true;

    avg   = (this->max + this->min) / 2.0;
    //index = ((position[0] < avg[0]) ? 0 : 1) | ((position[1] < avg[1]) ? 0 : 2) | ((position[2] < avg[2]) ? 0 : 4);

    for(index = 0; index < 8; ++index)
    {
#if 0
        if( box.getMax()[0] <= this->children[index]->getMin()[0] ||
            box.getMin()[0] >= this->children[index]->getMax()[0] ||
            box.getMax()[1] <= this->children[index]->getMin()[1] ||
            box.getMin()[1] >= this->children[index]->getMax()[1] ||
            box.getMax()[2] <= this->children[index]->getMin()[2] ||
            box.getMin()[2] >= this->children[index]->getMax()[2] )
        {
            intercepted = false;
        }
        else
        {
            intercepted = true;
        }
#endif

        if(this->children[index]->intersectBox(box))
        {
            split = false;

            if(enableSplit && BoundingBox::leafBox == this->children[index]->getType())
            {
                leafBox = (LeafBoundingBox *)this->children[index];

                if(leafBox->isFull() && leafBox->getSeparable())
                {
                    newBox = new InternalBoundingBox(this->children[index]->getMin(), this->children[index]->getMax(), 1);

                    if(NULL == newBox)
                    {
                        throw "not enough memory";
                    }

                    //newBox->produceLeafChildren();

                    for(iter = leafBox->getData().begin(); iter != leafBox->getData().end(); ++iter)
                    {
                        newBox->putDatum(iter->first, iter->second, false);
                    }

                    split = true;

                    for(j = 0; split && j < 8; ++j)
                    {
                        if(((LeafBoundingBox *)newBox->children[j])->getData().size() >= LeafBoundingBox::maxDataSize)
                        {
                            split = false;
                        }
                    }

                    if(split)
                    {
                        //std::cout << "Split box!" << std::endl;
                        this->children[index] = newBox;
                        delete leafBox;
                    }
                    else
                    {
                        leafBox->setSeparable(false);

                        //std::cout << "Cannot split box!" << std::endl;
                        delete newBox;
                    }
                }
            }

            if(this->children[index]->putDatum(box, datum))
            {
                this->emptyFlag = false;
            }
        }
    }

    for(index = 0; this->emptyFlag && index < 8; ++index)
    {
        if(!this->children[index]->empty())
        {
            this->emptyFlag = false;
        }
    }


    return !this->emptyFlag;
}

void InternalBoundingBox::getData(const Box &box, std::list<void *> &data) const
{
    unsigned int i;


    if(!this->emptyFlag && this->intersectBox(box))
    {
        for(i = 0; i < 8; ++i)
        {
            this->children[i]->getData(box, data);
        }
    }
}

#if 0
void *InternalBoundingBox::getDatum(const Vector3D &position)
{
    Vector3D avg;


    avg = (max + min) / 2.0;

    return (this->containsData ? this->children[((position[0] < avg[0]) ? 0 : 1) | ((position[1] < avg[1]) ? 0 : 2) | ((position[2] < avg[2]) ? 0 : 4)]->getDatum(position) : false);
}
#endif

void InternalBoundingBox::rayTrace(const Vector3D &start, const Vector3D &direction, void **nearestDatum, Vector3D &hitPoint, double &minDist)
{
    unsigned int i;


    if(!this->emptyFlag && this->intersectRay(start, direction))
    {
        for(i = 0; i < 8; ++i)
        {
            this->children[i]->rayTrace(start, direction, nearestDatum, hitPoint, minDist);
        }
    }
}
