/***************************************************************************
 *   Copyright (C) 2011 by Tian-Li Yu                                      *
 *   tianliyu@cc.ee.ntu.edu.tw                                             *
 *                                                                         *
 ***************************************************************************/


#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>

#include "statistics.h"
#include "dsmga.h"
#include "global.h"
#include "fitOC.h"
#define MAX_GEN 500

#define mincutPlusXO 0
#define mincutXO 5
#define SBSXO 1

using namespace std;


int
main (int argc, char *argv[])
{

    if (argc != 5 && argc !=7 ) {
        printf ("bisection ell numConvergence lower upper [spin #x#] [spin problem #]\n");
        return -1;
    }

    int ell = atoi (argv[1]);
    int numConvergence = atoi (argv[2]); // problem size
    int lower = atoi(argv[3]);
    int upper = atoi(argv[4]);

    int problemSize;
    int problemNum;
    if (argc == 7) {
        problemSize = atoi (argv[5]);
        problemNum = atoi (argv[6]);
    }

    //int nInitial = (int) (0.5 * ell * log((double)ell) / log(2.71828)); // initial population size
    int nInitial = 30; // initial population size
    if (lower > 0)
        nInitial = lower;

    //int nInitial = 10;

    //ufo.createProblem(ell, TRAP_K, 1.5, 0);
    //tsuji.sampling(0.6, TRAP_K, ell, 10000, 200);

    //myRand.seed(123);

    int j;

    Statistics st;

    int left, right, middle;

	/*
    if (argc == 7) {
        char filename[200];

        if (problemSize < 17)
            sprintf(filename, "./2D/%dx%d/%d_%d.%d", problemSize, problemSize, problemSize, problemSize*problemSize, problemNum);
        else
            sprintf(filename, "./2D/%dx%d/s_%dx%d.%d", problemSize, problemSize, problemSize, problemSize, problemNum);

        printf("Loading: %s\n", filename);
        FILE *fp = fopen(filename, "r");

        if (problemSize < 17)
            loadSimonsSpinGlassInstance(fp, &mySpinGlassParams);
        else
            loadSpinGlassInstance(fp, &mySpinGlassParams);
    }*/


    //myTableLookUp.readTable(hXOR);

    int populationSize = nInitial/2;
    bool foundOptima;
		double oRatio = 0.2;
		double cRatio = 0.1;
		int bbsize = 5;
	//
    int flag = fitoc.createProblem(ell, oRatio, cRatio, bbsize);
    if (flag != 0) return false;     	
	//
    if (lower < 0 || upper < 0) {

        if (SHOW_BISECTION) printf("Bisection phase 1\n");

        do {

            populationSize *= 2;

            if (SHOW_BISECTION) printf("[%d]: ", populationSize);

            foundOptima = true;


            for (j=0; j<numConvergence; j++) {

                DSMGA ga(ell, populationSize, 2, 1, 0, MAX_GEN, -1);
                ga.doIt(false, mincutPlusXO);


                if (!ga.foundOptima()) {

                    foundOptima = false;

                    if (SHOW_BISECTION) {
                        printf("-");
                        fflush(NULL);
                    }
                    break;
                }

                if (SHOW_BISECTION) {
                    printf("+");
                    fflush(NULL);
                }
            }

            if (SHOW_BISECTION) printf("\n");

        } while (!foundOptima);

        left = populationSize/2;
        right = populationSize;
    }

    else {
        left = lower;
        right = upper;
    }


    middle = (left + right)/2;
    Statistics stGen, stLS;

    if (SHOW_BISECTION) printf("Bisection phase 2\n");

    while ((right > 1.05 * left) && right > left + 2) {

        middle = (left + right) / 2;

        if (SHOW_BISECTION) printf("[%d]: ", middle);

        foundOptima = true;

        for (j=0; j<numConvergence; j++) {

            DSMGA ga(ell, middle, 2, 1, 0, MAX_GEN, -1);
            ga.doIt(false, mincutPlusXO);

            if (!ga.foundOptima()) {
                foundOptima = false;
                if (SHOW_BISECTION) {
                    printf("-");
                    fflush(NULL);
                }
                break;
            }
            if (SHOW_BISECTION) {
                printf("+");
                fflush(NULL);
            }
            if (j==0){
							stGen.reset();
							stLS.reset();
						}
            stGen.record(ga.getGeneration());
            stLS.record(ga.lsNFE);	
        }

        if (foundOptima)
            right = middle;
        else
            left = middle;


        if (SHOW_BISECTION) printf("\n");
    };



    //if (argc==7)
      //  freeSpinGlassInstance(&mySpinGlassParams);

    middle = (left + right) / 2;

    printf("%f %d %f %f\n", stLS.getMean(), middle, stGen.getMean(), middle*stGen.getMean());

    return EXIT_SUCCESS;

}

