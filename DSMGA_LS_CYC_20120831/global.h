/***************************************************************************
 *   Copyright (C) 2005 by Tian-Li Yu,,,                                   *
 *   tianliyu@fishlaptop.ytgroup                                           *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#ifndef _GLOBAL_H
#define _GLOBAL_H

#define NDEBUG // for assert
#include <cassert>
#include <cmath>

#include "UFO.h"
#include "fitOC.h"

#include "spin-glass.h"
#include "myrand.h"
#include "bitwisedistance.h"
#include "tablelookup.h"
#include "Tsuji.h"
class UFO;
class FitOverConflict;
//class Tsuji;
//#include "spin-glass.h"

#define CHI
#define TRAP_K 5
#define MAX_K 5
#define DIS_M 3
#define GROUP_SIZE 5
#define MAX_PARENT 15
//#define MAX_OVERLAP 5


extern int maxMemory;

extern bool SHOW_HC;
extern bool SHOW_DSM;
extern bool SHOW_LINKAGE;
extern bool SHOW_POPULATION;
extern bool SHORT_HAND;
extern bool SHOW_REPLACEMENT;
extern bool SHOW_SELECTION_INDEX;
extern bool SHOW_MAPPING;
extern bool SHOW_COMPRESSION;
extern bool SHOW_BUSINESSMEN;
extern bool SHOW_BISECTION;

extern char outputFilename[100];
extern void gstop ();
extern void outputErrMsg (char *errMsg);
extern int pow2 (int x);
extern void findMax (int number, int index[], int values[], int size);
extern void findMax (int number, int index[], double values[], int size);
extern void findMin (int number, int index[], int values[], int size);
extern void findMin (int number, int index[], double values[], int size);
extern int compareForBR (const void *a, const void *b);
extern int myDoubleCompare (const void *a, const void *b);
extern double lambertw(double w);

extern MyRand myRand;
extern BitwiseDistance myBD;
extern TableLookUp myTableLookUp;

extern char bbiErrorOutFile[100];
extern int* decisionV;
extern int V_size;
extern double dataStd;
extern int dataMin,dataMax;

extern UFO ufo;
extern Tsuji tsuji;
extern FitOverConflict fitoc;
extern SpinGlassParams mySpinGlassParams;

inline int quotientLong(int a)
{
    return (a / (sizeof(unsigned long) * 8) );
}


inline int remainderLong(int a)
{
    return (a & (sizeof(unsigned long) * 8 - 1));
}

inline double jointEntropy(double p00, double p01, double p10, double p11)
{
    double result = 0.0;
    result -= p00 * log(p00);
    result -= p01 * log(p01);
    result -= p10 * log(p10);
    result -= p11 * log(p11);

    return result;
}

inline double mutualInformation(double p00, double p01, double p10, double p11)
{
    double result = 0.0;

    double p0x = p00+p01;
    double p1x = p10+p11;
    double px0 = p00+p10;
    double px1 = p01+p11;
    
    // Bhattacharyya Distance
//    result += -log(sqrt(p0x*p0y) + sqrt(p1x*p1y));
    //result += (sqrt(p0x*p0y) + sqrt(p1x*p1y));

    // Hellinger Distance
    //result = 0.5 * sqrt(2 + 2 * result);

    
    result += (p00 < 1e-6) ? 0.0 : p00 * log (p00 / p0x / px0);
    result += (p01 < 1e-6) ? 0.0 : p01 * log (p01 / p0x / px1);
    result += (p10 < 1e-6) ? 0.0 : p10 * log (p10 / p1x / px0);
    result += (p11 < 1e-6) ? 0.0 : p11 * log (p11 / p1x / px1);
    
    return result;
}

inline double metric(double p00, double p01, double p10, double p11)
{
   return mutualInformation(p00,p01,p10,p11)/jointEntropy(p00,p01,p10,p11);
}

inline double square(double a)
{
    return a*a;
}

inline void swapInt(int *i, int *j)
{
    int temp;
    temp = *i;
    *i = *j;
    *j = temp;
}


int Partition(int * list,  const double * p, const int & start, const int& last ) ;
void QuickSort(int * list, const double * p, const int & start, const int& last ) ;
void generateMoneTest(const int& ell, const int& std,int repeat,int maxOverlap);
#endif
