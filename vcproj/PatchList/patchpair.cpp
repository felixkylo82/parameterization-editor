#include "patchpair.h"


PatchPairList::PatchPairList()
{
}

PatchPairList::~PatchPairList()
{
    this->clear();
}

void PatchPairList::clear()
{
    std::vector<PatchPair>::iterator ppiter;


    this->shadowList.clear();

    for(ppiter = this->patchPairs.begin(); ppiter != this->patchPairs.end(); ++ppiter)
    {
        delete ppiter->second;
    }

    this->patchPairs.clear();
}

#if 0

PatchPairList::PatchPairList(const PatchList &list1, const PatchList &list2)
{
    this->createPairList(list1, list2);
}

void PatchPairList::createPairList(const PatchList &list1, const PatchList &list2)
{
    const PatchList::FaceToPatchMap &map1 = list1.getFaceToPatchMap();
    const PatchList::FaceToPatchMap &map2 = list2.getFaceToPatchMap();


    PatchList::FaceToPatchMap::const_iterator miter1;
    PatchList::FaceToPatchMap::const_iterator miter2;


    if(map1.size() != map2.size())
    {
        throw "The sizes of maps do not match!";
    }

    this->clear();

    miter1 = map1.begin();
    miter2 = map2.begin();

    while(miter1 != map1.end())
    {
        this->push_back(PatchPair(miter1->second, miter2->second));

        ++miter1;
        ++miter2;
    }
}
#endif

void PatchPairList::createPairList(const PatchList &list1)
{
    const PatchList::FaceToPatchMap &map1 = list1.getFaceToPatchMap();

    Patch *newPatch;

    PatchList::FaceToPatchMap::const_iterator miter1;


    this->clear();

    miter1 = map1.begin();

    while(miter1 != map1.end())
    {
        newPatch = NULL;

        switch(miter1->second->getNumberOfControlPoints())
        {
        case 3:

            newPatch = new TriPatch(miter1->second->getControlPoints(), miter1->second->getControlNormal(), miter1->second->getResolution());

            break;

        case 4:

            newPatch = new QuadPatch(miter1->second->getControlPoints(), miter1->second->getControlNormal(), miter1->second->getResolution());

            break;
        }

        this->shadowList.push_back(newPatch);
        this->patchPairs.push_back(PatchPair(miter1->second, newPatch));

        ++miter1;
    }
}
