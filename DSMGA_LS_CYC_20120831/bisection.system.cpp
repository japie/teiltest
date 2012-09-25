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
 #include "hdsmga.h"
 #include "global.h"

 using namespace std;


 int
 main (int argc, char *argv[])
 {

   if (argc != 7)
     {
       printf ("GA ell numConvergence lower upper repeat nodeNumber\n");
       return -1;
     }

   int ell = atoi (argv[1]);
   int numConvergence = atoi (argv[2]);	// problem size
   int lower = atoi(argv[3]);
   int upper = atoi(argv[4]);
   int bisectionRepeat = atoi (argv[5]);	// how many time to repeat
   int nodeNum = atoi(argv[6]);	    // node number

   char tmpFilename[100];
   sprintf(tmpFilename, "dsmga-result-%04d", nodeNum);

   int nInitial = (int) (0.5 * ell * log((double)ell) / log(2.71828));	// initial population size

   //int nInitial = 10;


   int i, j;

   Statistics st;

   int left, right, middle;


      
   for (i = 0; i < bisectionRepeat; i++) {
   

       int populationSize = nInitial/2;
       bool foundOptima;
       
       if (lower < 0 || upper < 0) {

	   if (SHOW_BISECTION) printf("Bisection phase 1\n");

	   do {
	  
	       populationSize *= 2;

	       if (SHOW_BISECTION) printf("[%d]: ", populationSize);

	       foundOptima = true;


	       for (j=0; j<numConvergence; j++) {

		   char command[200];
		   sprintf(command, "nohup ./DSMGA++ %d %d 2 1 0 200 -1 1 0 -1 > %s", ell, populationSize, tmpFilename);
		   system(command);

		   char ch;
		   std::ifstream outfile(tmpFilename);
		   outfile >> ch;
                   outfile.close();


		   if (ch != '+') {
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

       if (SHOW_BISECTION) printf("Bisection phase 2\n");

       while ((right > 1.05 * left) && right > left + 2) {
	   
	   middle = (left + right) / 2;

	   if (SHOW_BISECTION) printf("[%d]: ", middle);

	   foundOptima = true;

	   for (j=0; j<numConvergence; j++) {

	       char command[200];
	       sprintf(command, "nohup ./DSMGA++ %d %d 2 1 0 200 -1 1 0 -1 > %s", ell, middle, tmpFilename);
	       system(command);
	       char ch;
	       std::ifstream outfile(tmpFilename);
	       outfile >> ch;
               outfile.close();

	       if (ch != '+') {
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

	   if (foundOptima) 
	       right = middle;
	   else 
	       left = middle;

	   if (SHOW_BISECTION) printf("\n");
       };


       st.record(middle);
   }

   printf("===============\n%f\n", st.getMean());

   return EXIT_SUCCESS;

 }

