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
#include <cstdlib>
#include <fstream>


#include "statistics.h"
#include "dsmga.h"
#include "global.h"

using namespace std;

Statistics *st_bb;


int
main (int argc, char *argv[])
{

    //myTableLookUp.readTable(hXOR);
    //myTableLookUp.readTable(hIFF);


    //maxMemory=0;

    if (argc != 12)
    {
        printf
            ("DSMGA ell nInitial selectionPressure pc pm maxGen maxFe repeat display rand_seed std\n");
        return -1;
    }

    int ell = atoi (argv[1]);	// problem size
    int nInitial = atoi (argv[2]);	// initial population size
    int selectionPressure = atoi (argv[3]);	// selection pressure
    double pc = atof (argv[4]);	// pc
    double pm = atof (argv[5]);	// pm
    int maxGen = atoi (argv[6]);	// max generation
    int maxFe = atoi (argv[7]);	// max fe
    int repeat = atoi (argv[8]);	// how many time to repeat
    int display = atoi (argv[9]); // display each generation or not
    int rand_seed = atoi (argv[10]); 	// rand seed
    int std = atoi(argv[11]);

    if (rand_seed != -1)  // time
        myRand.seed((unsigned long)rand_seed);

    int i;

    Statistics stGen;
    int usedGen;

    int failNum = 0;

    int maxMemoryUsage = 0;
    if ( std != -1 )
      generateMoneTest(ell,std);
      
   
    
    
    for (i = 0; i < repeat; i++) {

        DSMGA dsmga (ell, nInitial, selectionPressure, pc, pm, maxGen, maxFe);
       
        if (display == 1)
            usedGen = dsmga.doIt (true);
        else
            usedGen = dsmga.doIt (false);


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
    cout  << "Memory usage of DSM + population: " << (double)maxMemory/(double)repeat + nInitial*ell/8 << " bytes." << endl;
    printf ("\n");
    printf ("%f  %d\n", stGen.getMean (), failNum);

    return EXIT_SUCCESS;
}
