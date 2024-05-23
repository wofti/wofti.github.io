/* bam.h */
/* Bernd Bruegmann, 12/99 */


#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <malloc.h>
#include <time.h>



#ifdef REDUCEORDERTO2
#define REDUCEORDERTO4
#endif
#ifdef REDUCEORDERTO4
#define REDUCEORDERTO6
#endif
#ifdef REDUCEORDERTO6
#define REDUCEORDERTO8
#endif


//mth: NO SUPPORT FOR NOT USING THESE TWO FLAGS ANY LONGER
#ifndef BOX
#define BOX
#endif

#ifndef BOXES
#define BOXES
#endif




#include "bam_automatic_include.h"
