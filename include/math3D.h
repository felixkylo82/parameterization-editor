///////////////////////////////////////////////////////////////////////
//
//  Manipulating Movie Path
//
//  math3D.h
//
//  written by Philip Fu (cwfu@acm.org)
//
//  - copyright
//  - please notify me for any modification
//
//  Mon Jan 22 7:54:04 EST 2001
//
///////////////////////////////////////////////////////////////////////




#ifndef __MATH3D_
#define __MATH3D_


//#define _DEBUG


#ifndef _EPSILON
#define _EPSILON 1e-7
#endif

#ifndef _TRUE
#define _TRUE 1
#endif

#ifndef _FALSE
#define _FALSE 0
#endif


#ifdef __cplusplus
extern "C" {
#endif


void printVec(const double *v, int c);


// check if trangle v1,v2,v3 is anticlockwise w.r.t. +z
extern int anticlockwise2D(const double *v1, const double *v2, const double *v3);

extern double clamp(const double val1);

extern double dist2D(const double *pt1, const double *pt2);
extern double dist3D(const double *pt1, const double *pt2);

extern double length2D(const double *pt1);
extern double length3D(const double *pt1);

extern double sqrLength2D(const double *pt1);
extern double sqrLength3D(const double *pt1);

extern void add2D(const double *vect1, const double *vect2, double *vect3);
extern void add3D(const double *vect1, const double *vect2, double *vect3);

extern void minus2D(const double *vect1, const double *vect2, double *vect3);
extern void minus3D(const double *vect1, const double *vect2, double *vect3);

extern void add_mul_3D(const double *vect1, const double t, const double *vect2, double *v1_add_t_mul_v2);
extern void add_mul_2D(const double *vect1, const double t, const double *vect2, double *v1_add_t_mul_v2);

extern double dot2D(const double *vect1, const double *vect2);
extern double dot3D(const double *vect1, const double *vect2);

// v cross w = r
extern void cross3D(const double *v, const double *w, double *r);

extern void copy2D(double *vect1, const double *vect2);
extern void copy3D(double *vect1, const double *vect2);

extern int equalVec2D(const double *vect1, const double *vect2);
extern int equalVec3D(const double *vect1, const double *vect2);

extern void mult2D(double *vect, const double scale);
extern void mult3D(double *vect, const double scale);

extern void div2D(double *vect, const double scale);
extern void div3D(double *vect, const double scale);

extern int unify2D(double *vect);
extern int unify3D(double *vect);

extern void findNormal3D(const double *pt1, const double *pt2, const double *pt3, double *normal);
extern void average3D(const double *plane1, const double *plane2, const double *plane3, double *normal);

extern double hitTriangle(const double *ver1,   const double *ver2, const double *ver3,
                          const double *rayOrg, const double *rayDir,
                          double *hitCoord, double *hitNormal);
extern double hitSphere(const double *center, const double radius,
                        const double *rayOrg, const double *rayDir,
                        double *hitCoord, double *hitNormal);

extern int anticlockwise2D(const double *v1, const double *v2, const double *v3);

extern int areParallelVec2D(const double *vect1, const double *vect2);
extern int areParallelVec3D(const double *vect1, const double *vect2);

extern int insideTri2D ( float testX  , float testY  ,
                         float triV1X , float triV1Y ,
                         float triV2X , float triV2Y ,
                         float triV3X , float triV3Y ,
                         float epsilon ) ;
extern int insideTri3D ( double * testPt , double * triV1 , double * triV2 , double * triV3 , double tolerance ) ;

extern double area2D ( double x0 , double y0 ,
                       double x1 , double y1 ,
                       double x2 , double y2 ) ;
extern double area3D ( double * triV1 , double * triV2 , double * triV3 ) ;


// angle is in radian
extern void rotate2D(double *outVec, const double *inVec, const double angle);
extern void rotate3D_3(double *inOutVec, const double *axis, const double angle);


extern void printMat( FILE *stream, const char *format, const double *D, const int row, const int column );
extern void transposeMat( double *D, const int n );
extern void multMat( const double *A, const double *B, double *result, const int n );

// cast 4x4 matrix to be identity and 4x4 matrix inversion
extern void identity(double m[16]);
extern int  invert4by4(double src[16]);

// invert the square matrix where D is a n by n row-wise matrix
// return 1 on error else 0
extern int invertMat( double *D, const int n );


// Float-based functions

void  rotVec      ( const float * ptrV , const float * rotMat , float * ptrVOut ) ;
void  invRotVec   ( const float * ptrV , const float * rotMat , float * ptrVOut ) ;
void  minus3Df    ( const float * v1 , const float * v2 , float * v_result ) ;
int   equalVec3Df ( const float * v1 , const float * v2 , float epsilon ) ;
void  add3Df      ( const float * v1 , const float * v2 , float * v_result ) ;
void  add_mul_3Df ( const float * vect1 , const float t , const float * vect2 , float * v1_add_t_mul_v2 ) ;
void  mult3Df     ( float * v_result , float multvalue ) ;
void  copy3Df     ( float * v_result , const float * v_src ) ;
void  cross3Df    ( const float * v , const float * w , float * r ) ;
float dist3Df     ( const float * pt1 , const float * pt2 ) ;
void  unify3Df    ( float * vect ) ;
int   equalVecf   ( float * v1 , float * v2 , double eps ) ;
float dot3Df      ( const float *vect1, const float *vect2);
float length3Df   ( const float *pt1);
float area3Df     ( float * triV1 , float * triV2 , float * triV3 ) ;


#ifdef __cplusplus
}
#endif


#endif
