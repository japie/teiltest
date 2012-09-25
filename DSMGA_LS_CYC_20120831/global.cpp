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



#include <cstdio>
#include <cstdlib>
#include "global.h"
#include "myrand.h"
#include "statistics.h"
#include "Tsuji.h"
#include <cstring>
#include <iostream>
#include <fstream>

using namespace std;

int  maxMemory = 0;
bool SHOW_HC = false;
bool SHOW_DSM = false;
bool SHOW_LINKAGE = false;
bool SHOW_POPULATION = false;
bool SHORT_HAND = false;
bool SHOW_SELECTION_INDEX = false;
bool SHOW_REPLACEMENT = false;
bool SHOW_MAPPING = false;
bool SHOW_COMPRESSION = false;
bool SHOW_BUSINESSMEN = false;
bool SHOW_BISECTION = true;

char bbiErrorOutFile[100];
int* decisionV;
int V_size = 0;
double dataStd = 0;
int dataMax = 0;
int dataMin = INT_MAX;



SpinGlassParams mySpinGlassParams;

char outputFilename[100];
MyRand myRand;
BitwiseDistance myBD;
TableLookUp myTableLookUp;
Tsuji tsuji;
UFO ufo;
FitOverConflict fitoc;

void outputErrMsg (char *errMsg)
{
  printf ("%s\n", errMsg);
  exit (1);
}

int
pow2 (int x)
{
  return (1 << x);
}

void
findMax (int number, int index[], double values[], int size)
{
  int i, j, k;
  double *max = new double[number];


  for (i = 0; i < number; i++)
    {
      index[i] = -1;
      max[i] = -1;
    }

  for (i = 0; i < size; i++)
    for (j = 0; j < number; j++)

      if (values[i] > max[j])
	{
	  for (k = number - 1; k > j; k--)
	    {
	      index[k] = index[k - 1];
	      max[k] = max[k - 1];
	    }
	  index[j] = i;
	  max[j] = values[i];

	  break;
	}

  delete[]max;
}

/** find number of maxima and return them in index, only for positive number*/
void
findMax (int number, int index[], int values[], int size)
{
  int i, j, k;
  int *max = new int[number];


  for (i = 0; i < number; i++)
    {
      index[i] = -1;
      max[i] = -1;
    }

  for (i = 0; i < size; i++)
    for (j = 0; j < number; j++)

      if (values[i] > max[j])
	{
	  for (k = number - 1; k > j; k--)
	    {
	      index[k] = index[k - 1];
	      max[k] = max[k - 1];
	    }
	  index[j] = i;
	  max[j] = values[i];

	  break;
	}

  delete[]max;
}

/** findMin--int */
void
findMin (int number, int index[], int values[], int size)
{
  int i, j, k;
  int *min = new int[number];


  for (i = 0; i < number; i++)
    {
      index[i] = -1;
      min[i] = (int)INF * 2;
    }

  for (i = 0; i < size; i++)
    for (j = 0; j < number; j++)

      if (values[i] < min[j])
	{
	  for (k = number - 1; k > j; k--)
	    {
	      index[k] = index[k - 1];
	      min[k] = min[k - 1];
	    }
	  index[j] = i;
	  min[j] = values[i];

	  break;
	}

  delete[]min;
}


/** findMin--double */
void
findMin (int number, int index[], double values[], int size)
{
  int i, j, k;
  double *min = new double[number];


  for (i = 0; i < number; i++)
    {
      index[i] = -1;
      min[i] = INF * 2;
    }

  for (i = 0; i < size; i++)
    for (j = 0; j < number; j++)

      if (values[i] < min[j])
	{
	  for (k = number - 1; k > j; k--)
	    {
	      index[k] = index[k - 1];
	      min[k] = min[k - 1];
	    }
	  index[j] = i;
	  min[j] = values[i];

	  break;
	}

  delete[]min;
}

int myDoubleCompare(const void *a, const void *b)
{
	if (*((double *)a) > *((double *)b))
		return 1;
	else if (*((double *)a) < *((double *)b))
		return -1;
	else 
		return 0;
}

// return the x which x*exp(x) = w
double lambertw(double w)
{
    double upper = 100.0;
    double lower = 0.0;

    assert(upper * exp(upper) > w);
    assert(w > 0);

    double middle;

    while ((upper - lower) / upper > 0.01) {
	middle = (upper + lower) / 2;
	
	double z = middle * exp(middle);
	if (z < w)
	    lower = middle;
	else
	    upper = middle;
    }

    return (upper + lower) / 2;
}


   //PARTITION(A, p, r)
    //1. x „] A[p] /* break up A wrt x */
    //2. i „] p -1
    //3. j „] r +1
    //4. while TRUE do
    //        5. repeat j „] j -1
    //        6. until A[j] <= x
    //        7. repeat i „] i +1
    //        8. until A[i] >= x
    //        9. if i < j
    //            10. then exchange A[i] „\ A[j]
    //            11. else return j

int Partition(int * list,  const double * p, const int & start, const int& last ) {
    double x = p[list[start]];
    int i = start-1;
    int j = last+1;
    while (1)
    {
        while( p[list[--j]] < x )
        {
        }
        while( p[list[++i]] > x )
        {
        }
        if (i < j)
            swapInt(&list[i], &list[j]);
        else
            return j;
    }
}

//QUICKSORT(A, p, r)
///* Call QUICKSORT(A, 1, length[A]) to sort an entire array */
//1. if p < r then
    //    2. q „] PARTITION(A, p, r)
    //    3. QUICKSORT(A, p, q)
    //    4. QUICKSORT(A, q+1, r)

void QuickSort(int * list, const double * p, const int & start, const int& last ) {
    if ( start < last )
    {
        int q = Partition( list, p, start, last );
        QuickSort(list, p , start, q);
        QuickSort(list, p , q+1, last);
    }
}


void generateMoneTest(const int& ell, const int& std,int repeat,int maxOverlap) {
    V_size = ell/DIS_M*GROUP_SIZE;
    decisionV = new int[V_size];
    int* cache = new int[ell]();
    for ( int i = 0 ; i < ell ; i++ ) cache[i] = 0;
    if ( std != -1 )
    {
         
         int counter = 0;
         char filename[200];
         string token;
         sprintf(filename, "NormalSet/NSet_%d_%d_%d_%d_%d.dat", ell, ell/DIS_M ,maxOverlap,std,repeat);
         printf("Loading: %s\n", filename);
         //FILE *fp = fopen(filename, "r");
         ifstream input;
         input.open(filename,ifstream::in);
         if(!input.is_open())
         {
            cout << "wrong file path, enter new file path:"<<endl;
            cin >> filename;
            input.open(filename,ifstream::in);
         }
         while(input)
         {
            input >> token ;
            counter++;
            if ( counter > V_size)
              break;
            decisionV[counter-1] = (atoi(token.c_str()));
            token.clear();
            cache[decisionV[counter-1]]++;
            
         }
         /*for( int i = 0; i < V_size ;i++)
         {
         cout << decisionV[i] << " ";
         if(i%GROUP_SIZE == (GROUP_SIZE-1) )
         cout << endl;
         }*/
         Statistics geneSt;
         for( int i = 0 ; i < ell ; i++)
	 {
	    geneSt.record(cache[i]);
	 }
         (geneSt.getMax() > dataMax)?dataMax = geneSt.getMax():dataMax = dataMax;
	 dataStd = geneSt.getStdev();
	 (dataMin < geneSt.getMin())?dataMin = geneSt.getMin():dataMin = dataMin;
	   
    }

    

}
