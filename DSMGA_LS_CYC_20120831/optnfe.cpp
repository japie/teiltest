#include <math.h>
#include <iostream>
#include <fstream>

#include "statistics.h"
#include "dsmga.h"
#include "global.h"

#define MAX_GEN 300

// bisection
#define NUMCONV 10
#define BI_REPEAT 50
#define RELIABLE 50

// optnfe
#define PRECISION 0.01
#define GA_REPEAT 50
#define NFE_REPEAT 50
#define COUNTER 1

bool SHOW_DETAIL = false;

using namespace std;

int
main (int argc, char *argv[])
{

	if (argc != 7)
	{
		printf ("GA ell lower upper std mode maxOverlap\n");
		return -1;
	}

	int ell = atoi(argv[1]);
	int lower = atoi(argv[2]);
	int upper = atoi(argv[3]);
        int std = atoi(argv[4]);
        int mode = atoi(argv[5]);
        int maxOverlap = atoi(argv[6]);
        //if ( std != -1 )
        //    generateMoneTest(ell,std);


	int round;
	int failround;

	int i, j;
	int populationSize;
	int minpopulationSize = 0;

//
// bisection
//

	if (SHOW_DETAIL)
	{
		printf("***************\n");
		printf("Bisection Phase\n");
		printf("***************\n");
	}

	failround = 0;
	for(round=0; round<BI_REPEAT; round++)
	{
//		populationSize = (int) (0.5 * ell * log((double)ell) / log(2.71828)) / 2;
		populationSize = lower/2;

                if ( populationSize <= 0 )
                    populationSize = 1;

		int left, right, middle;
		bool foundOptima = false;

		// bisection phase 1

		if (lower < 0 || upper < 0)
		{
			if (SHOW_DETAIL)
			{
				printf("round %d phase 1\n", round+1);
				printf("---------------\n");
			}

			do
			{
				populationSize *= 2;

				if (SHOW_DETAIL)
					printf("[%d]: ", populationSize);

				foundOptima = true;

				for (j=0; j<NUMCONV; j++)
				{
				        if(std!=-1)
					    generateMoneTest(ell,std,j,maxOverlap);
					DSMGA dsmga(ell, populationSize, 2, 1, 0, MAX_GEN, -1);
					//dsmga.doIt(false);
                                        dsmga.doIt(false,mode);

					if (!dsmga.foundOptima())
					{
						foundOptima = false;

						if (SHOW_DETAIL)
						{
							printf("-");
							fflush(NULL);
						}

						break;
					}

					if (SHOW_DETAIL)
					{
						printf("+");
						fflush(NULL);
					}
				}

				if (SHOW_DETAIL)
					printf("\n");
			} while (!foundOptima);

			left = populationSize/2;
			right = populationSize;
		}
		else
		{
			left = lower;
			right = upper;
		}

		if (SHOW_DETAIL)
			printf("===============\n");

	// bisection phase 2
	
		if (SHOW_DETAIL)
		{
			printf("round %d phase 2\n", round+1);
			printf("---------------\n");
		}

		while ((right > 1.05 * left) && right > left + 2)
		{
			middle = (left + right) / 2;

			if (SHOW_DETAIL)
				printf("[%d]: ", middle);

			foundOptima = true;

			for (j=0; j<NUMCONV; j++)
			{
				DSMGA dsmga(ell, middle, 2, 1, 0, MAX_GEN, -1);
				//dsmga.doIt(false);
                                dsmga.doIt(false,mode);

				if (!dsmga.foundOptima())
				{
					foundOptima = false;

					if (SHOW_DETAIL)
					{
						printf("-");
						fflush(NULL);
					}

					break;
				}

				if (SHOW_DETAIL)
				{
					printf("+");
					fflush(NULL);
				}
			}

			if (foundOptima)
				right = middle;
			else
				left = middle;

			if (SHOW_DETAIL)
				printf("\n");
		};

		middle = (left + right) / 2;

		if (SHOW_DETAIL)
		{
			printf("===============\n");
			printf("%d\n", middle);
			printf("===============\n\n");
		}

		// if can find reliable population size
		if (!foundOptima && right == upper)
		{
			failround++;
		}
		else
			minpopulationSize += middle;
	}

	if (failround != BI_REPEAT)
		minpopulationSize /= (BI_REPEAT - failround);

	if (SHOW_DETAIL)
	{
		double failrate = ( 100.0 * (double)failround ) / (double)BI_REPEAT;
		if ((BI_REPEAT - failround) == RELIABLE)
		{
			printf("minimum population size: %d\n\n", minpopulationSize);
			printf("failed rate: %f %c\n\n", failrate, '%');
		}
		else
		{
			printf("Bisection cannot find reliable population size in [%d, %d]\n", lower, upper);
			printf("failed rate: %f %c\n\n", failrate, '%');
			return EXIT_SUCCESS;
		}
	}
//
// optnfe
//

	if (SHOW_DETAIL)
	{
		printf("***************\n");
		printf("Opt. nfe Phase\n");
		printf("***************\n");
	}

	double avgoptnfe[3] = {0.0, 0.0, 0.0};
	failround = 0;
	for(round=0; round<NFE_REPEAT; round++)
	{
		double optimalnfe[3];	// nfe, population size, convergence time
		double history[5][3];
		for(i=0; i<5; i++)
			for(j=0; j<3;j++)
				history[i][j] = 0.0;

	// nfe phase 1

		populationSize = minpopulationSize;
		bool foundLower = false;

		double memo[3];
		if (SHOW_DETAIL) 
		{
			printf("round %d phase 1\n", round+1);
			printf("---------------\n");
		}	

		while (1)
		{
			if (SHOW_DETAIL)
			{
				printf("[%d]: ", populationSize);
				fflush(NULL);
			}

			Statistics stGen;
	    	int usedGen;
			for (j=0; j<GA_REPEAT; j++)
			{
				DSMGA dsmga(ell, populationSize, 2, 1, 0, MAX_GEN, -1);
				usedGen = dsmga.doIt(false);
                                usedGen = dsmga.doIt(false, mode);

				if (dsmga.foundOptima())
					stGen.record (usedGen);

				if (SHOW_DETAIL && ((j % COUNTER) == 0))
				{
					printf(">");
					fflush(NULL);
				}
			}
			memo[1] = (double)populationSize;
			memo[2] = stGen.getMean();
			memo[0] = memo[1] * memo[2];
			if (SHOW_DETAIL)
				printf(" %f = %d x %f\n", memo[0], (int)memo[1], memo[2]);

			if (!foundLower)
			{
				history[0][0] = memo[0];
				history[0][1] = memo[1];
				history[0][2] = memo[2];
				foundLower = true;
			}
			else
			{
				if ( (memo[0] >= history[0][0]) || ((memo[0] >= history[4][0]) && (history[4][0] != 0.0)) )
				{
					double memo_swap = memo[0];
					memo[0] = history[4][0];
					history[4][0] = memo_swap;

					memo_swap = memo[1];
					memo[1] = history[4][1];
					history[4][1] = memo_swap;

					memo_swap = memo[2];
					memo[2] = history[4][2];
					history[4][2] = memo_swap;
					break;
				}
				else
				{
					if (history[4][0] >= memo[0])
					{
						history[0][0] = history[4][0];
						history[0][1] = history[4][1];
						history[0][2] = history[4][2];
					}
					history[4][0] = memo[0];
					history[4][1] = memo[1];
					history[4][2] = memo[2];

					if (history[0][0] == history[4][0])
						break;
				}
			}

			populationSize += (int)(minpopulationSize * 0.1);
		}

		if (SHOW_DETAIL)
			printf("===============\n");

	// nfe phase 2
   
		if (SHOW_DETAIL)
			printf("round %d phase 2\n", round+1);

		double precision = (double)lower * PRECISION;
		if (precision < 1.0)
			precision = 1.0;

		double history_temp[3] = {0.0, 0.0, 0.0};
		int skip = 0;
		bool flag = false;
		while (1)
		{
			if ( (history[4][1] - history[0][1]) <= precision )
			{
				double temp = history[0][0];
				int index = 0;

				for (i=1; i<=4; i++)
				{
					if ( (history[i][0]<temp) && (history[i][0]!=0.0) && (history[i][1]>=history[0][1]) && (history[i][1]<=history[4][1]) )
					{
						temp = history[i][0];
						index = i;
					}
				}
				optimalnfe[0] = history[index][0];
				optimalnfe[1] = history[index][1];
				optimalnfe[2] = history[index][2];

				break;
			}
			else
			{
				if (SHOW_DETAIL)				
				{
					printf("---------------\n");
					printf("[%d <-> %d]\n", (int)history[0][1], (int)history[4][1]);
				}

				for (i=0; i<=4; i++)
				{
					populationSize = (int)( ( (double)(4-i)*history[0][1] + (double)i*history[4][1] ) / 4.0 );

					if (SHOW_DETAIL)
					{
						printf("[%d]: ", populationSize);
						fflush(NULL);
					}

					// skip
					if (i==0 || i==4)
						skip = 1;
					else
					{
						if (populationSize == (int)memo[1])
							skip = 2;
						else if (populationSize == (int)history_temp[1])
							skip = 3;
						else if (populationSize == (int)history[i-1][1])
							skip = 4;
						else
							skip = 0;
					}

					if (skip != 0)
					{
						if (skip == 2)
						{
							history[i][0] = memo[0];
							history[i][1] = memo[1];
							history[i][2] = memo[2];
						}
						else if (skip == 3)
						{
							history[i][0] = history_temp[0];
							history[i][1] = history_temp[1];
							history[i][2] = history_temp[2];
						}
						else if (skip == 4)
						{
							history[i][0] = history[i-1][0];
							history[i][1] = history[i-1][1];
							history[i][2] = history[i-1][2];
						}

						if (SHOW_DETAIL)
						{
							printf("skipped ");
							for (j=0; j<(GA_REPEAT/COUNTER-8); j++)
								printf(">");
						}
					}
					else
					{
						Statistics stGen;
    				    int usedGen;
						for (j=0; j<GA_REPEAT; j++) 
						{
							DSMGA dsmga(ell, populationSize, 2, 1, 0, MAX_GEN, -1);
	    	    		    //usedGen = dsmga.doIt(false);
                                    usedGen = dsmga.doIt(false, mode);
	
				            if (dsmga.foundOptima())
    			        	    stGen.record (usedGen);

							if (SHOW_DETAIL && ((j % COUNTER) == 0))
							{
								printf(">");
								fflush(NULL);
							}	
						}
	
						history[i][1] = (double)populationSize;
			    	    history[i][2] = stGen.getMean();
						history[i][0] = history[i][1] * history[i][2];
					}
	
					if (SHOW_DETAIL)
						printf(" %f = %d x %f\n", history[i][0], (int)history[i][1], history[i][2]);

					if ( (i != 0) && ((history[i][0]-history[i-1][0]) >= 0.0) && (skip != 4) )
					{
						if (history[i][0]-history[i-1][0]==0)
							flag = true;
						else
							flag = false;

						break;
					}
				}

				int a, b, c;
				switch (i)
				{
					case 1:
						a=0; b=1; c=0;
						break;
					case 2:
						a=0; b=2; c=1;
						if (flag) {	a=1; b=2; c=0; }
						break;
					case 3:
						a=1; b=3; c=2;
						if (flag) {	a=2; b=3; c=0; }
						break;
					case 4:
						a=2; b=4; c=3;
						if (flag) {	a=3; b=4; c=0; }
						break;
					default:
						a=3; b=4; c=0;
				}
				history[0][0] = history[a][0];
				history[0][1] = history[a][1];
				history[0][2] = history[a][2];
				history[4][0] = history[b][0];
				history[4][1] = history[b][1];
				history[4][2] = history[b][2];
				if (c != 0)
				{
					history_temp[0] = history[c][0];
					history_temp[1] = history[c][1];
					history_temp[2] = history[c][2];
				}
			}
		}

		if (SHOW_DETAIL)
		{
			printf("===============\n");
			printf("%f  =  %d  x  %f\n", optimalnfe[0], (int)optimalnfe[1], optimalnfe[2]);
			printf("===============\n\n");
		}
		else
			printf("%f %d %f\n", optimalnfe[0], (int)optimalnfe[1], optimalnfe[2]);

		avgoptnfe[0] += optimalnfe[0];
		avgoptnfe[1] += optimalnfe[1];
		avgoptnfe[2] += optimalnfe[2];
	}

	avgoptnfe[0] /= (double)NFE_REPEAT;
	avgoptnfe[1] /= (double)NFE_REPEAT;
	avgoptnfe[2] /= (double)NFE_REPEAT;

	if (SHOW_DETAIL)
		printf("optimal nfe:  %f  =  %f  x  %f  (population size  x  convergence time)\n", avgoptnfe[0], avgoptnfe[1], avgoptnfe[2]);
	else
		printf("%f %f %f\n", avgoptnfe[0], avgoptnfe[1], avgoptnfe[2]);

	return EXIT_SUCCESS;
}

