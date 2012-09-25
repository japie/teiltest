/* 
 * File:   testUFO.cpp
 * Author: Rick
 *
 * Created on January 18, 2011, 1:36 AM
 */

#include <math.h>
#include <iostream>
#include <cstdlib>
#include "statistics.h"
#include "dsmga.h"
#include "global.h"
#include "UFO.h"
#define REPEAT 10
#define MAX_GEN 200

using namespace std;

Statistics *st_bb;


int main (int argc, char *argv[])
{

    //myTableLookUp.readTable(hXOR);
    //myTableLookUp.readTable(hIFF);


    //maxMemory=0;

    if (argc != 8)
    {
        printf
            ("DSMGA nInitial totalBB bbSize omega conflict display mode \n");
        return -1;
    }

    int nInitial = atoi (argv[1]);
    int totalBB = atoi (argv[2]);	// problem size
    int bbSize = atoi (argv[3]);
    double omega = atof (argv[4]);
    int conflict = atoi (argv[5]);
    int display = atoi (argv[6]);
    int mode = atoi(argv[7]);



    //int selectionPressure = atoi (argv[3]);	// selection pressure
    int selectionPressure = 2;
    //double pc = atof (argv[4]);	// pc
    double pc = 1;
    //double pm = atof (argv[5]);	// pm
    double pm = 0;
    //int maxGen = atoi (argv[6]);	// max generation
    int maxGen = MAX_GEN;
    //int maxFe = atoi (argv[7]);	// max fe
    int maxFe = -1;
    //int repeat = atoi (argv[8]);	// how many time to repeat
    int repeat = REPEAT;
    //int display = atoi (argv[9]); // display each generation or not
    //int rand_seed = atoi (argv[10]); 	// rand seed
    int rand_seed = -1;
    //int std = atoi(argv[11]);

    if (rand_seed != -1)  // time
        myRand.seed((unsigned long)rand_seed);

    int i;

    Statistics stGen;
    int usedGen;
    int failNum = 0;
    int maxMemoryUsage = 0;

    int ell = 0;

    for (i = 0; i < repeat; i++)
    {

        ufo.createProblem(totalBB,bbSize, omega, conflict);
        ell = ufo.getNumOfGenes();

        ufo.printStatistics();

        printf("conflict : %d\n", ufo.getZeros());
        printf("max fitness : %f\n", ufo.getMaxFitness());


        DSMGA dsmga (ell, nInitial, selectionPressure, pc, pm, maxGen, maxFe);

        if (display == 1)
            usedGen = dsmga.doIt (true,mode);
        else
            usedGen = dsmga.doIt (false,mode);


        if (!dsmga.foundOptima()) {
            failNum++;
            printf ("-");
        }
        else {
            stGen.record (usedGen);
            printf ("+");
        }

        maxMemoryUsage += maxMemory;
        fflush (NULL);

    }

    cout << endl;
    cout  << "Max DSM memory usage:" << (double)maxMemoryUsage/(double)repeat << " bytes." << endl;
//    cout  << "Memory usage of DSM + population: " << (double)maxMemory/(double)repeat + nInitial*ell/8 << " bytes." << endl;
    printf ("\n");
    printf ("generation :  min/mean/max/std  %f / %f / %f / %f\n", stGen.getMin(), stGen.getMean(), stGen.getMax(), stGen.getStdev() );

    return EXIT_SUCCESS;
}

