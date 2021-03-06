

#include "dsmga.h"
/***************************************************************************
  *   Copyright (C) 2010 by Rick Chang                                      *
  *   r99921037@ntu.edu.tw                                       *
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
#include "global.h"
#include "UFO.h"

using namespace std;
#define MAX_GEN 200

 char * xoMethod = "mune";

 int main (int argc, char *argv[] )  {

   if (argc != 10)
     {
       printf ("repeated-bisection numConvergence lower ell bbSize omega std totalTests fileCount mode\n");
       return -1;
     }

   int numConvergence = atoi (argv[1]);
   int lower = atoi(argv[2]);
   int ell = atoi(argv[3]);
   int bbSize = atoi(argv[4]);
   double omega = atof(argv[5]);	    // use how many problems tested
   double std = atof(argv[6]);
   int totalTests = atoi(argv[7]);
   int fileCount = atoi(argv[8]);
   int mode = atoi(argv[9]);



   myTableLookUp.readTable(hIFF);

   char totalFileName[100];

   int totalBB = 0;


    if (mode == 0)
        xoMethod = "mune";
    else if (mode == 1)
        xoMethod = "tlxo";
    else if (mode == 2)
        xoMethod = "gibbs";
    else if (mode == 3)
        xoMethod = "onePointXO";
    else if (mode == 4)
        xoMethod = "uniformXO";
    else if (mode == 5)
        xoMethod = "simpleMinCut";
    else
        xoMethod = "bbUniformXO";

    Statistics nfeSt;
    char outFileName[100];

    Statistics ufoOmegaSt;
    Statistics ufoConflictSt;
    Statistics popSt;

   for ( int d =0; d < totalTests; ++d )
   {

       lower = lower / 8;
       /*
       if ( ! ufo.createProblem(ell,bbSize,omega, conflict) )
       {
           printf("Can not create the demanded problem. \n");
           continue;
       }*/
       tsuji.verbose();
       tsuji.getBBI(omega,bbSize,ell,std,d);
       //totalBB = ufo.getNumOfBBs();
       totalBB = tsuji.getM();
       cout<< d <<"th,  m = " << totalBB << "  k = " << bbSize << "  omega = " << omega << "  ell = " << ell << " std = " << std << "  \n";
       //ufo.print();
       //ufo.printStatistics();
       //printf("conflict : %d\n", ufo.getTotalConflictedBB());

       //ufoOmegaSt.record(ufo.getMeanOmega());
       //ufoConflictSt.record((double)ufo.getTotalConflictedBB());

       
       //sprintf(totalFileName, "./result/totalResult-%s-m_%03d-k_%d-w_%f-ell_%d-c_%d-total_%02d", xoMethod, totalBB, bbSize, omega, ell, ufo.getTotalConflictedBB(), totalTests );

       sprintf(totalFileName, "./result/totalResult-%s-m_%03d-k_%d-w_%f-ell_%d-std_%f-total_%02d", xoMethod, totalBB, bbSize, omega, ell, std, totalTests );


       sprintf(bbiErrorOutFile, "./bbi_error/%s-m_%03d-k_%d-w_%f-ell_%d-std_%f", xoMethod, totalBB, bbSize, omega, ell, std );






       sprintf(outFileName, "./result/tsuji-%s-m_%03d-k_%d-w_%f-ell_%d-std_%f-total_%02d-%d", xoMethod, totalBB, bbSize, omega, ell, std, totalTests, fileCount);
       ofstream outFile(outFileName, ios_base::app);
        if (!outFile.is_open())
        {
            cout << "Error opening file " << outFileName <<endl;
            return -1;
        }

       //if ( std != -1 )
          //generateMoneTest(ell,std);

       int nInitial = (int) (0.5 * ell * log((double)ell) / log(2.71828));	// initial population size
       if (nInitial < lower * 1.05) nInitial = (int)(lower * 1.05);

       outFile<< d <<"th,  m = " << totalBB << "  k = " << bbSize << "  omega = " << omega << "  ell = " << ell << " std = " << std << "  \n";

   
       




       int j;

       Statistics stGen;
       int populationSize = nInitial;
       int numOptima = 0;
       
//       if(std != -1)
//           generateMoneTest(ell,std,d,maxOverlap);

       do {
               if (SHOW_BISECTION)
                   if (stGen.getNumber() <= 0)
                       printf("[%d]: ", populationSize);
           //if(std != -1)
           //  generateMoneTest(ell,std,numOptima,maxOverlap);
             
           DSMGA dsmga(ell, populationSize, 2, 1, 0, MAX_GEN, -1);
           dsmga.doIt(false,mode);

           if (!dsmga.foundOptima()) {

                   if (SHOW_BISECTION) {
                       printf("- [gen:%f]\n", stGen.getMean());
                       fflush(NULL);
                   }

               lower = populationSize;
               numOptima = 0;
               populationSize *= 2;
               stGen.reset();

           }
           else {

                   if (SHOW_BISECTION) {
                       printf("+");
                       fflush(NULL);
                   }

               if ((populationSize < lower * 1.05 ) || (populationSize < lower + 2)) {
                   numOptima ++;
                   stGen.record(dsmga.getGeneration());
               }
               else {
                   numOptima = 0;
                   populationSize = (populationSize + lower) / 2;
                   printf(" [gen:%d]\n", dsmga.getGeneration());
                   stGen.reset();
               }

           }


       } while (numOptima < numConvergence);



       printf("\n===============\n%d %f\n", populationSize, stGen.getMean());
       double nfe = populationSize*stGen.getMean();
       nfeSt.record(nfe);
       outFile << "populationSize : " << populationSize << " MeanGeneration : " << stGen.getMean() << " nfe : " << nfe <<endl;
       outFile.close();
       popSt.record(populationSize);
   }
   

   ofstream totalOutFile(totalFileName, ios_base::app);
   totalOutFile << " fileCount = " <<fileCount  << "  meanOmega = " << omega  << "  std = " << std << "  nfe =  (min/mean/max/std)  "
           << nfeSt.getMin() <<" / " << nfeSt.getMean() << " / " << nfeSt.getMax() << " / "<< nfeSt.getStdev() << "  mean population = " << popSt.getMean()  <<endl;
   fflush(NULL);

   totalOutFile.close();
   return 1;

 }



