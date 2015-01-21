/*
 * VPTree
 *
 * standard.h - header file for pre-defined declarations 
 *              and general data structure
 *
 * By Philip Fu
 *
 * 1/13/2001
 *
 */


#ifndef _STANDARD
#define _STANDARD

#include <math.h>



///////////////////////////////////////////////////////////////////////
//
//  Global debug identifier
//

//#define _DEBUG



///////////////////////////////////////////////////////////////////////
//  
//  General Definition
//

#define _MAX_BUF_SIZE	255
#define _MAX_STR_SIZE	255
#define _MAX_PATH_SIZE	255

#define _OKAY		1
#define _ERROR		-1
#define _TRUE		1
#define _FALSE		0

#ifdef FLT_EPSILON
#define _EPSILON FLT_EPSILON
#else
#define _EPSILON	1e-7
#endif

#ifndef M_PI
#define M_PI	3.1415926535897932384626433832795
#endif

#define TWO_PI (2.0*M_PI)
#define HALF_PI (0.5*M_PI)

#define toDegree(X)     ((X)*180.0/M_PI)
#define toRadian(X)     ((X)*M_PI/180.0)

#define errexit(msg)			{ printf(msg);           exit(1); }
#define errexit2(msg,arg1)		{ printf(msg,arg1);      exit(1); }
#define errexit3(msg,arg1,arg2)	{ printf(msg,arg1,arg2); exit(1); }

#ifndef max
#define max(x,y)		( (x) > (y) ? (x) : (y) )
#endif

#ifndef min
#define min(x,y)		( (x) < (y) ? (x) : (y) )
#endif

#define toggle(bool)		( (bool) == _TRUE ? _FALSE : _TRUE )



///////////////////////////////////////////////////////////////////////
// All Path Information should reference to _TOP_DIR

#define _TOP_DIR	"/nfs/italy/home/cwfu/research/vpnav/"



#endif
