/***************************************************************************
 *   Copyright (C) 2005 by Tian-Li Yu                                      *
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
#ifndef _DSMCLUSTERINGCHROMOSOME_H
#define _DSMCLUSTERINGCHROMOSOME_H

/**
  @author Tian-Li Yu
 */

#include "chromosome.h"
#include "bbi.h"
#include "twodarray.h"
#include "statistics.h"
#include "triMatrix.h"

class DSMClusteringChromosome
{
    public:
        DSMClusteringChromosome ();

        ~DSMClusteringChromosome ();

        void init (int n_part_n, double w1, double w2);
        void initialGuess ();
        void initialGuess (int *guess);

        bool partInCluster (int p, int c);
        bool partInSomeCluster (int p);
        void invert (int p, int c);
        void addPartToCluster (int p, int c);
        void removePartFromCluster (int p, int c);

        void createDSMV2 (Chromosome * population, int populationSize, int selectionPressure, int *selectionIndex);
        void createDSMCC (Chromosome * population, int populationSize, int selectionPressure, int *selectionIndex);

        void createDSMNonlinearity (Chromosome * population, int populationSize, int selectionPressure, int *selectionIndex);

        double getInitialFitness ();
        double getFitnessIfInvertPC (double originalFitness, int p, int c);
        double getDeltaFitnessIfInvertPC (int p, int c);

        void speedyHillClimbingV3 ();             // modify only row and column

        BBI getBBI ();

        int *gene;
        bool hasInitialGuess;

    private:

        BBI clusterInformation;
        TwoDArray < int > clusterNumberArray;


        bool inCluster (int clusterSize, int start[], int end[], int x, int y);
        bool inSomeCluster (int clusterSize, TwoDArray < int >&fastIndex, int partI,
            int partJ);

        bool inPrePre (TwoDArray < int >&fastIndex, int part, int current,
            int clusterSize);
        bool inPre (TwoDArray < int >&fastIndex, int part, int current,
            int clusterSize);
        bool inPost (TwoDArray < int >&fastIndex, int part, int current,
            int clusterSize);

        void swapInt (int *a, int *b);
        int getValue (int p, int c);
        void setValue (int p, int c, int val);

	double gaussianSplit (double mean, double var);
	double gammaSplit (double mean, double var);
	double betaSplit (double mean, double var);
	double gammaLn(double x);
	double incompleteGamma(double a, double x);
	double incompleteBeta(double a, double b, double x);
	double betacf(double a, double b, double x);
	void computeIMeanVar(double *mean, double *variance, double q00, double q01, double q10, double q11, double mi, int n);


        static int partN;
        static double w1;
        static double w2;
        static double w3;
        static TriMatrix<double> dsm;

};
#endif
