//
// common.h
//
// By Tien-Tsin Wong 
//

#ifndef __COMMON_H
#define __COMMON_H

#include <stdio.h>
#include <stdlib.h>

#define MAX(a,b)		((a>b)? a : b)
#define MIN(a,b)		((a<b)? a : b)
#define ERREXIT(s)		{fprintf(stderr,s); exit(1);}
#define ERREXIT2(s1,s2) {fprintf(stderr,s1,s2); exit(1);}
#define ERRMSG(s)		{fprintf(stderr,s);}

#define DEL_PTR(p)		{if (p) {delete p; p=NULL;}}
#define DEL_ARR(p)		{if (p) {delete [] p; p=NULL;}}

#define FALSE  0
#define TRUE   1

#define FAIL   0
#define OKAY   1


//#ifdef WIN32
//
//// These values are defined in value.h on Unix platform,
//// but not on the Windows, so we define them here.
//#ifndef M_PI
//#define	M_PI		3.14159265358979323846
//#define	M_2_PI		M_PI * 2
//#define	M_PI_2		M_PI / 2
//#endif
//
//#define MAXFLOAT	3.402823466e38
//#define MINFLOAT	1.175494351e-38
//
//#endif

// Simple window coordinate calculation routines.
void SpanReorder(int *interval);
void WindowIntersect(int *xwin1, int *ywin1, int *xwin2,
		     int *ywin2, int *xout,  int *yout);

// Some utility routines to handle cross platform issues
void SwapByte(unsigned char *buf, size_t elemsize, size_t nelem);
size_t cp_fread (void *buf, size_t size, size_t n, FILE* file);
size_t cp_fwrite(void *buf, size_t size, size_t n, FILE* file);

// Get the file size of binary file
long getFileSize(FILE *fptr);

// It returns seconds elapsed since start of program, as a double
double myTime();


#endif	// __COMMON_H

