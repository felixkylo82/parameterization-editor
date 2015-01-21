/*
 * vptree.h
 *
 * VP-tree Definition
 *
 * By Philip Fu
 *
 * Sun Jan 14 20:53:26 EST 2001
 *
 */




#ifndef _VPTREE
#define _VPTREE



#ifdef __cplusplus
extern "C" {
#endif

////////////////////////////////////////////////////////////////
// Global Definitions
////////////////////////////////////////////////////////////////


#define _NUM_PT_RANDOM_SUBSET	1000	// Minimum number of points in Random subset
#define _NUM_CAND_VP			100		// Number of candidate VPs

#define _DEFAULT_BRANCH			2		// Branching factor of the vp-tree (must > 1)

#define _MAX_POINTS_NODE		10		// Maximum number of points in a node
#define _MIN_POINTS_NODE		1		// Minimum number of points in a node

#define _VPTREE_ID				"VPTREE-FILE"
#define _VPTREE_LEAF			"LEAF"
#define _VPTREE_INTERNAL		"INTERNAL"




////////////////////////////////////////////////////////////////
// Data Structure
////////////////////////////////////////////////////////////////


typedef struct
{
    int  nPoints;
    int  dimension;

    double *points;		// size : numPoints * dimension

} DataSet;


typedef struct VPNodeStru
{
    int isLeaf;			// _TRUE or _FALSE

    struct VPNodeStru **child;	// children (array of numBranch pointers)
    double *medians;		// array of (numBranch-1) distances 
				// separating each group such that the
				// (i)th group is bounded by medians[i]

    int nPoints;		// number of data points

    double *points;		// if (isLeaf)
				//    this is the data points
				// else
				//    this is the vantage Point(s)
				//    (note that vp is also data point)

    int *dataID;		// dataIDs' keep track of where the 
                                // points goes in the VP tree

} VPNode;


typedef struct
{
    int numBranch;
    int dimension;

    VPNode *root;

} VPTree;




////////////////////////////////////////////////////////////////
// Function Declarations
////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////
// (1) Reading the dataset
//     (return NULL on error)
//
//     Format of dataFile :
//     <number of points>
//     <dimension>
//     <x1> <y1> <z1>
//     ......
//     ......

DataSet *readDataSet(char *dataFile);



/////////////////////////////////////////////////////////
// (2) Building the VP Tree
//     - if numBranch    <= 1, use _DEFAULT_BRANCH
//     - if numDimension <= 0, use _DEFAULT_DIMENSION
//     - note that ordering of points in dataSet is changed after calling
//       (return the root node)

VPTree *
buildVPTree(DataSet *dataSet, int numBranch);

// For internal use only!
VPNode *
buildVPNode(double *points, int *dataID, int nPoints, 
            int numBranch, int dimension);



/////////////////////////////////////////////////////////
// (3) Read/Write VPTree file

// on error, return NULL
VPTree *readVPTree(char *treeFile, int *nPoints);

// on error, return -1, else number of data points written
int writeVPTree(char *treeFile, VPTree *vpTree);



/////////////////////////////////////////////////////////
// (4) K-Nearest Neighbor Search


// Perform K-Nearest Neighbor Search
// Input :
// - vpTree
// - queryPt
// - k
// Return :
// - list of points -> resultPt is returned in argument
// - list of dataID -> resultID is returned in argument
// Note : - memory for resultID and resultPt should have been
//        - allocated before calling.

void
knnsearch(VPTree *vpTree, double *queryPt, int k,
          double *resultPt, int *resultID);



/////////////////////////////////////////////////////////
// (5) Free the resource and simple utility

void freeDataSet(DataSet *dataSet);
void freeVPTree(VPTree *vpTree);

// For Internal use only!
void freeVPNode(VPNode *vpNode, int numBranch);

// Print point
void printPt(FILE *out, double *pt, int dim);


// [Felix] 

DataSet *
mallocDataSet(unsigned int nPoints, unsigned int dimension);

#ifdef __cplusplus
}
#endif


#endif
