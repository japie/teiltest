/* 
 * File:   testConflict.cpp
 * Author: Rick
 *
 * Created on January 20, 2011, 1:18 AM
 */


#include <stdlib.h>
#include "UFO.h"
#include "statistics.h"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    if ( argc != 4 )
    {
        printf (" testConflict totalBB bbSize omega\n");
        return 0;
    }

    int totalBB = atoi(argv[1]);
    int bbSize = atoi(argv[2]);
    double omega = atof(argv[3]);

    char fileName[100];
    sprintf(fileName, "./conflict/m_%d-k_%d-w_%f", totalBB, bbSize, omega);

    FILE * outfile;
    outfile = fopen (fileName, "a");
    if (outfile == NULL)
    {
        printf("Can not open file %s and write.", fileName);
    }

    UFO ufo1;

//    Statistics st;
//
//    for ( double i = 0; i < 100; ++i )
//    {
//        ufo1.createProblem(totalBB, bbSize, 2.221,0);
//        st.record(ufo.getUnusedGenes());
//        ufo1.printStatistics();
//    }
//
//    printf("mean Unused = %f\n",st.getMean());
//
//    ufo1.print();



//    for ( double omega = 1; omega < 10.1; omega+=0.2)
//    {
//        Statistics st0;
//        Statistics st1;
//        Statistics st2;
//        Statistics st3;
//        Statistics st4;
//        Statistics st5;
//        int max = -1;
//        for ( int i = 0; i < 1000; ++i )
//        {
//
//            ufo1.createProblem(totalBB, bbSize, omega, 0);
//            if ( ufo1.getG2BMaxDifference() > max )
//                max = ufo1.getG2BMaxDifference();
//            //fprintf(outfile, "%f   %d\n", ufo.getMeanOmega(), ufo.getZeros());
//            //st.record((double)ufo.getZeros());
//            st0.record((double)ufo1.getMeanOmega());
//            st1.record((double)ufo1.getG2BStd());
//            st2.record((double)ufo1.getB2BStd());
//            st3.record((double)ufo1.getMaxOmega());
//            st4.record((double)ufo1.getB2BMean());
//            st5.record((double)ufo1.getMinOmega());
//            //fprintf(stdout, "%f   %d\n", ufo.getMeanOmega(), ufo.getZeros());
//        }
//
//        fprintf(outfile, "omega = %f   mean = %f   g2bStd = %f   b2BMean = %f   b2bStd = %f   maxOmega = %f  minOmega = %f\n", omega, st0.getMean(), st1.getMean(), st4.getMean(), st2.getMean(), st3.getMean(), st5.getMean() );
//
//        fprintf(stdout, "omega = %f   mean = %f   g2bStd = %f   b2BMean = %f   b2bStd = %f   maxOmega = %f  minOmega = %f\n", omega, st0.getMean(), st1.getMean(), st4.getMean(), st2.getMean(), st3.getMean(), st5.getMean() );
//        fprintf(stdout, "max difference = %d\n", max);
//    }



    for ( double omega = 1; omega < 10.1; omega+=0.2)
    {
        Statistics st0;
        for ( int i = 0; i < 10; ++i )
        {

            ufo1.createProblem(totalBB, bbSize, omega, -1);
            //fprintf(outfile, "%f   %d\n", ufo.getMeanOmega(), ufo.getZeros());
            st0.record((double)ufo1.getZeros());
            //st0.record((double)ufo1.getMeanOmega());

            //fprintf(stdout, "%f   %d\n", ufo.getMeanOmega(), ufo.getZeros());
        }

        fprintf(outfile, "omega = %f   mean = %f\n", omega, st0.getMean());
        fprintf(stdout, "omega = %f   mean = %f\n", omega, st0.getMean());
    }

    

    fclose (outfile);


    return 0;
}

