/*
* sort.h
*
* utility for sorting!
*
* Sat Jan 13 19:34:22 EST 2001
*
*/


#ifndef _SORT
#define _SORT


#ifdef __cplusplus
extern "C" {
#endif

	// Quick sort D[] and the corresponding comp[] in ascending order of
	// values of comp[].

	void q_sort(int *D, double *comp, int low, int high, int size);


#ifdef __cplusplus
}
#endif


#endif
