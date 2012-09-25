 /***************************************************************************
  *   Copyright (C) 2004 by Tian-Li Yu                                      *
  *   tianliyu@illigal.ge.uiuc.edu                                          *
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


#include <math.h>
#include <iostream>
#include <fstream>

#include "statistics.h"
#include "dsmga.h"
#include "global.h"
#include "fitOC.h"
//#define MAX_GEN 500
#define EXP_TIME 10
#define SHOW_BISECTION true
Statistics POPU_ST;
Statistics NFE_ST;
using namespace std;

#define mincutPlusXO 0
#define mincutXO 5
#define SBSXO 1

bool dsmgaOneRun(int ell, int runSize, int maxGen, int& usedGen){
    DSMGA dsmga(ell, runSize, 2, 1, 0, maxGen, -1);
                        //output, mode
    int gen = dsmga.doIt(false, SBSXO);    
    if (!dsmga.foundOptima()) {
      if (SHOW_BISECTION) {
        //printf("-");
        //fflush(NULL);
      }
      usedGen = 0;
      return false;
    }
    else{ 
      if (SHOW_BISECTION) {
        //printf("+");
        //fflush(NULL);
      }
      usedGen = gen;
      return true;
    }
}


bool experimentOneRun(int ell, int bbsize, float oRatio, float cRatio, int lower, int upper, 
                      int convergenceCombo, int maxGen){
    int nInitial = (int) (0.5 * ell * log((double)ell) / log(2.71828));	
    int left, right, middle;   
    myTableLookUp.readTable(hIFF);
    int populationSize = nInitial/2;
    bool foundOptima = false;
    int flag = ocFitness.createProblem(ell, oRatio, cRatio, bbsize);
    if (flag != 0) return false;     
    if (lower < 0 || upper < 0) {
     if (SHOW_BISECTION) printf("Bisection phase 1\n");
      
      do {
	      populationSize *= 2;
	      if(populationSize >= 50000){
	        cout << "pupulation exceed 50000, use 50000" << endl;
	        populationSize = 50000;
	        break;
	      }
	      if (SHOW_BISECTION) printf("[%d]\n", populationSize);
	      foundOptima = true;
	      int usedGen = 0;
	      for (int j = 0; j < convergenceCombo; j++) {
		      foundOptima = dsmgaOneRun(ell, populationSize, maxGen, usedGen);
		      if (!foundOptima)
		        break;
                      else
                        cout << "+" << endl;
              }
              cout << endl;
	      //if (SHOW_BISECTION) printf("\n");
      } 
      while (!foundOptima);
      
      left = populationSize/2;
      right = populationSize;
    }
    else {
	   left = lower;
	   right = upper;
    }
    if (SHOW_BISECTION) printf("\n");    
/////////////////////////////////////////////////////////////////    
    middle = (left + right)/2;
    if (SHOW_BISECTION) printf("Bisection phase 2\n");
    while ((right > 1.05 * left) && right > left + 2) {	   
      middle = (left + right) / 2;
	    if (SHOW_BISECTION) printf("[%d] left=%d right=%d\n", middle, left, right);
	    foundOptima = true;
	    int usedGen = 0;
	    for (int j = 0; j < convergenceCombo; j++){
	      foundOptima = dsmgaOneRun(ell, middle, maxGen, usedGen);
	      if (!foundOptima)
	        break;
        else
          cout << "+" << endl;  
	    }
      cout << endl;

	    if (foundOptima) 
	      right = middle;
	    else 
	      left = middle;
      //if (SHOW_BISECTION) printf("\n");
    };
    middle = (left + right) / 2;
    if (SHOW_BISECTION) printf("\n");
    //printf("===============\n%d\n", middle);
    if (SHOW_BISECTION) printf("Bisection phase 3\n");
    printf("using popultion size:[%d]\n", middle);
    int nfe = 0;
    int successTime = 0;
    for (int j = 0; j < convergenceCombo; j++){
      int usedGen = 0;
      foundOptima = dsmgaOneRun(ell, middle, maxGen, usedGen);
      if (foundOptima){
        cout << "+" << endl;
        nfe += usedGen*middle;
        successTime += 1;
      }
      cout << endl;
    }    
    //printf("\n===============\n%d\n", middle);
    cout << "Population: " << middle << endl;
    cout << "Nfe: " << nfe/(float)successTime << endl;
    POPU_ST.record(middle);
    NFE_ST.record(nfe/(float)successTime);
    return true;
}


int main (int argc, char *argv[]){

    if (argc != 9)
    {
      printf ("bisection ell convergenceCombo lower upper bbsize oRatio cRatio maxGen\n");
      return -1;
    }

    int ell = atoi (argv[1]);
    int convergenceCombo = atoi (argv[2]);	// problem size
    int lower = atoi(argv[3]);
    int upper = atoi(argv[4]);
    //int std = atoi(argv[5]);
    //int maxOverlap = atoi(argv[6]);
    int bbsize = atoi (argv[5]);
    double oRatio = atoi (argv[6]);
    double cRatio = atoi(argv[7]);
    oRatio /= 100.0;
    cRatio /= 100.0;
    int maxGen = atoi(argv[8]);
    // initial population size
    bool settingOk = true;
    POPU_ST.reset();
    NFE_ST.reset();
    for (int i = 0; i < EXP_TIME; ++i){
      cout << "-------" << i+1 << "-------" << endl;
      settingOk = experimentOneRun(ell, bbsize, oRatio, cRatio, 
                  lower, upper, convergenceCombo, maxGen);
      if(!settingOk){
        cout << "this setting is not OK" << endl;
        return -1;
      }
    }
    cout << "===============" << endl;
    cout << POPU_ST.getMean() << endl;
    cout << NFE_ST.getMean() << endl;
    cout << POPU_ST.getVariance() << endl;
    cout << NFE_ST.getVariance() << endl;
}

