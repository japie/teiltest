/*
 * File:   testGlobal.cpp
 * Author: Rick
 *
 * Created on January 20, 2011, 1:18 AM
 */


#include <stdlib.h>
#include "UFO.h"
#include "statistics.h"
#include "chromosome.h"

using namespace std;

/*
 *
 */
int main(int argc, char** argv) {

    if ( argc != 5 )
    {
        printf ("testConflict ell bbSize omega conflict\n");
        return 0;
    }

    int ell = atoi(argv[1]);
    int bbSize = atoi(argv[2]);
    double omega = atof(argv[3]);
    int  conflict = atof(argv[4]);

//    UFO ufo;
//    ufo.createProblem(20,5,4,3);
//    ufo.print();



//    char fileName[100];
//    sprintf(fileName, "./conflict/m_%d-k_%d-w_%f", totalBB, bbSize, omega);
//
//    FILE * outfile;
//    outfile = fopen (fileName, "a");
//    if (outfile == NULL)
//    {
//        printf("Can not open file %s and write.", fileName);
//    }
    printf("ell = %d\n",ell);
    UFO ufo1;
    ufo1.createProblem(ell, bbSize, omega, conflict);
    printf("ConflictedBB = %d\n", ufo1.getTotalConflictedBB());
//    ufo1.print();
//    ufo1.checkGlobalOptimum();


//    for ( int j = 0; j < 10000; ++j )
//    {
//
//        UFO ufo1;
//        ufo1.createProblem(ell, bbSize, omega, conflict);
//        if ( ! ufo1.checkGlobalOptimum() )
//        {
//            printf("%d\n", j);
//            break;
//        }
//
////        int ell = ufo1.getNumOfGenes();
////        printf("ell = %d\n", ell);
////        printf("conflict = %d\n", ufo1.getZeros());
////
////        Chromosome c(ell);
////        for ( int i = 0; i < ell; ++i )
////        {
////            c.setVal(i,0);
////        }
////
////
////        for ( int i = 0; i < pow(2.0,(double)ell); ++i )
////        {
////            for ( int q = 0; q < ell; ++q )
////            {
////                c.setVal(q,0);
////            }
////            int still = i;
////            int index = 0;
////            while (still > 0)
////            {
////                c.setVal(index, still%2);
////                still = still /2;
////                ++index;
////            }
////            c.printOut();
////            if ( ufo1.getFitness(c) > ufo1.getMaxFitness() && i != pow(2.0,(double)ell) -1 )
////            {
////                printf("here!!!! ");
////                system("pause");
////            }
////            printf("\n");
////        }
//    }





    return 0;
}

