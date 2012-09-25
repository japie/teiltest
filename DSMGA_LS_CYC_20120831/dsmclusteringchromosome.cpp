
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
#include "dsmclusteringchromosome.h"
#include "global.h"
#include "fastcounting.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

//extern int maxMemory;

int DSMClusteringChromosome::partN;
double DSMClusteringChromosome::w1;
double DSMClusteringChromosome::w2;
double DSMClusteringChromosome::w3;
TriMatrix <double> DSMClusteringChromosome::dsm;

DSMClusteringChromosome::DSMClusteringChromosome ()
{
    gene = NULL;
    hasInitialGuess = false;
}


DSMClusteringChromosome::~DSMClusteringChromosome ()
{
    delete[]gene;
}


void
DSMClusteringChromosome::init (int n_partN, double n_w1, double n_w2)
{
    partN = n_partN;
    w1 = n_w1;
    w2 = n_w2;
    w3 = 1.0 - w1 - w2;

    gene = new int[partN * partN];

    for (int i=0; i<partN * partN; i++)
        gene[i] = 0;

    dsm.init(partN);
    clusterNumberArray.init (partN, partN);
    clusterInformation.init (partN);
    clusterInformation.bbNum = partN;
}


void
DSMClusteringChromosome::initialGuess (int *guess)
{
    if (hasInitialGuess)
        return;

    hasInitialGuess = true;

    int i,j;

    for (i = 0; i < partN; i++)
        for  (j = 0; j < partN; j++)
    {
        if (guess[i*partN+j] == 1)
            addPartToCluster (i, j);
        else
            removePartFromCluster (i, j);
    }

}


void
DSMClusteringChromosome::initialGuess ()
{
    if (hasInitialGuess)
        return;

    hasInitialGuess = true;

    int i, j;
    for (i = 0; i < partN; i++)
        for (j = 0; j < partN; j++)
            if (i == j)
                addPartToCluster (i, j);
    else
        removePartFromCluster (i, j);

}


int
DSMClusteringChromosome::getValue (int p, int c)
{
    return gene[c * partN + p];
}


void
DSMClusteringChromosome::setValue (int p, int c, int val)
{
    int i;

    int originalValue = getValue (p, c);

    if (originalValue == val)
        return;

    if (val == 0)
    {

        for (i = 1; i <= clusterInformation.bb[c][0]; i++)
        {
            int part = clusterInformation.bb[c][i];
            if (part != p)
            {
                clusterNumberArray (part, p)--;
                clusterNumberArray (p, part)--;
            }
        }

        for (i = 1; i <= clusterInformation.bb[c][0]; i++)
            if (clusterInformation.bb[c][i] == p)
        {
            // set it to the last one
            clusterInformation.bb[c][i] =
                clusterInformation.bb[c][clusterInformation.bb[c][0]];
            break;
        }
        clusterInformation.bb[c][0]--;

    }
    else
    {

        for (i = 1; i <= clusterInformation.bb[c][0]; i++)
        {
            int part = clusterInformation.bb[c][i];
            if (part != p)
            {
                clusterNumberArray (part, p)++;
                clusterNumberArray (p, part)++;
            }
        }

        clusterInformation.bb[c][clusterInformation.bb[c][0] + 1] = p;
        clusterInformation.bb[c][0]++;
    }

    // Last thing: set the value!
    gene[c * partN + p] = val;
}


bool
DSMClusteringChromosome::partInCluster (int p, int c)
{
    return (getValue (p, c) == 1);
}


bool
DSMClusteringChromosome::partInSomeCluster (int p)
{
    int c;
    for (c = 0; c < partN; c++)
        if (p != c)
            if (partInCluster (p, c))
                return true;

    return false;

}


void
DSMClusteringChromosome::addPartToCluster (int p, int c)
{
    setValue (p, c, 1);
}


void
DSMClusteringChromosome::removePartFromCluster (int p, int c)
{
    setValue (p, c, 0);
}


void
DSMClusteringChromosome::invert (int p, int c)
{
    setValue (p, c, 1 - getValue (p, c));
}


double
DSMClusteringChromosome::getInitialFitness ()
{
    int cluster_n = partN;
    int *cl = new int[cluster_n];
    int clusterSize = 0;

    int i, j;
    double fitness;
    double fit1, fit2, fit3;

    for (i = 0; i < cluster_n; i++)
    {
        cl[i] = 0;
        for (j = 0; j < partN; j++)
            cl[i] += partInCluster (j, i) ? 1 : 0;
        if (cl[i] > 1)
            clusterSize++;
    }

    // MODEL description
    fit1 = 0.0;
    for (i = 0; i < cluster_n; i++)
    {
        if (cl[i] > 1)
            fit1 += log2 (partN) * (1 + cl[i]);
    }

    fit2 = fit3 = 0.0;

    for (i = 0; i < partN; i++)
    {
        for (j = 0; j < partN; j++)
        {

            if (i == j) continue;

            if (dsm (i, j) > 1e-10)
            {

                fit3 += (2 * log2 (partN) + 1);
            }
        }
    }

    double w3 = (1 - w1 - w2);

    fitness = w1 * fit1 + w2 * fit2 + w3 * fit3;

    delete []cl;

    return fitness;

}


/** This is a speedy version of
 *  invert(p,c);
 *  getFitness();
 *  invert(p,c);
 */
double
DSMClusteringChromosome::getFitnessIfInvertPC (double originalFitness, int p, int c)
{
    return originalFitness + getDeltaFitnessIfInvertPC(p, c);
}


double
DSMClusteringChromosome::getDeltaFitnessIfInvertPC (int p, int c)
{
    int i;
    double fit1 = 0.0;
    double fit2 = 0.0;
    double fit3 = 0.0;
    double log2partN = log2 (partN);

    int val = getValue (p, c);
    int oldCL = clusterInformation.bb[c][0];

    int newCL = (val == 0) ? (oldCL + 1) : (oldCL - 1);

    if (newCL > MAX_K) return INF;

    if (oldCL == 1)                               // a new cluster is created
        fit1 = log2partN * 2;
    else if (newCL == 1)                          // a cluster is deleted
        fit1 = -log2partN * 2;
    else
        fit1 = log2partN * (newCL - oldCL);

    if (val == 0)
    {                                             // DsmOneLine increases

        for (i = 1; i <= oldCL; i++)
        {
            int part = clusterInformation.bb[c][i];

            if (part != p && clusterNumberArray(part,p) == 0)
            {                                     // dsm1(part,p): 0->1
                if (dsm (part, p) == 1)
                    fit3 -= (2 * log2partN + 1);
                else
                    fit2 += (2 * log2partN + 1);
            }

        }
    }
    else
    {                                             // DsmOneLine decreases

        for (i = 1; i <= oldCL; i++)
        {
            int part = clusterInformation.bb[c][i];

            if (part != p && clusterNumberArray(part,p) == 1)
            {                                     // dsm1(part,p): 1->0
                if (dsm (part, p) == 1)
                    fit3 += (2 * log2partN + 1);
                else
                    fit2 -= (2 * log2partN + 1);
            }

        }

    }

    double w3 = (1 - w1 - w2);

    double fitness = w1 * fit1 + 2 * (w2 * fit2 + w3 * fit3);

    return fitness;

}


void
DSMClusteringChromosome::speedyHillClimbingV3 ()
{
    int i, j;

    initialGuess ();
    double bestFitness = getInitialFitness ();

    int generation = 0;

    if (SHOW_HC)
        printf ("%d --> %f\n", generation, bestFitness);

    bool improvement;
    int *index = new int[partN * partN];
    double *cache = new double[partN * partN];

    typedef struct tagEntry
    {
        int part;
        int cluster;
        double deltaFitness;
    } Entry;

    int indexP = 0;
    int indexC = 0;

    TwoDArray <bool> entryUsed(partN, partN);
    Entry *entry = new Entry[partN * partN];
    int numEntry = 0;

    for (i=0; i<partN; i++)
        for (j=0; j<partN; j++)
    {
        if (i==j)
            cache[i*partN+j] = INF;
        else if (i<j)
            cache[i*partN+j] = getDeltaFitnessIfInvertPC (i, j);
        else
            cache[i*partN+j] = cache[j*partN+i];  // symetrc at the beginning

        if (cache[i*partN+j] < 0)
        {
            entry[numEntry].deltaFitness =  cache[i*partN+j];
            entry[numEntry].part = i;
            entry[numEntry].cluster = j;
            numEntry++;
            entryUsed(i, j) = true;
        }
        else
            entryUsed(i, j) = false;
    }

    while (numEntry > 0)
    {
        // Find the deepest scent

        double minDeltaFitness = INF;
        int minIndex = 0;

        for (i=0; i<numEntry; i++)
            if (minDeltaFitness > entry[i].deltaFitness)
        {
            minDeltaFitness = entry[i].deltaFitness;
            minIndex = i;
        }

        indexP = entry[minIndex].part;
        indexC = entry[minIndex].cluster;

        invert (indexP, indexC);
        improvement = true;

        bestFitness += entry[minIndex].deltaFitness;

        generation++;
        if (SHOW_HC)
            printf ("%d --> %f\n", generation, bestFitness);

        for (i=0; i<numEntry; i++)
        {
            double deltaFitness = getDeltaFitnessIfInvertPC(entry[i].part, entry[i].cluster);
            if (deltaFitness < 0)
            {
                entry[i].deltaFitness = deltaFitness;
            }
            else                                  // remove this entry
            {
                entryUsed(entry[i].part, entry[i].cluster) = false;
                entry[i].part = entry[numEntry-1].part;
                entry[i].cluster = entry[numEntry-1].cluster;
                i--;
                numEntry--;
            }
        }

        // recheck all part with indexC
        for (i=0; i<partN; i++)
        {
            if ( !entryUsed(i, indexC) && i!= indexC)
            {
                double deltaFitness = getDeltaFitnessIfInvertPC(i, indexC);
                if (deltaFitness < 0)
                {
                    entry[numEntry].part = i;
                    entry[numEntry].cluster = indexC;
                    entry[numEntry].deltaFitness = deltaFitness;
                    numEntry++;
                    entryUsed(i, indexC) = true;
                }
            }
        }

    }

    delete[]entry;
    delete[]index;
    delete[]cache;
}


inline void
DSMClusteringChromosome::swapInt (int *i, int *j)
{
    int temp;
    temp = *i;
    *i = *j;
    *j = temp;
}


bool DSMClusteringChromosome::inPrePre (TwoDArray < int >&fastIndex, int part, int current, int clusterSize)
{
    if (current < 2)
        return false;
    return (fastIndex (current - 2, part) == 1);
}


bool DSMClusteringChromosome::inPre (TwoDArray < int >&fastIndex, int part, int current, int clusterSize)
{
    if (current < 1)
        return false;
    return (fastIndex (current - 1, part) == 1);
}


bool DSMClusteringChromosome::inPost (TwoDArray < int >&fastIndex, int part, int current, int clusterSize)
{
    if (current >= clusterSize - 1)
        return false;
    return (fastIndex (current + 1, part) == 1);
}


bool DSMClusteringChromosome::inCluster (int clusterSize, int start[], int end[], int x, int y)
{

    int i;
    for (i = 0; i < clusterSize; i++)
    {
        if (start[i] <= x && x <= end[i] && start[i] <= y && y <= end[i])
            return true;
    }

    return false;
}


bool DSMClusteringChromosome::inSomeCluster (int clusterSize, TwoDArray < int >&fastIndex, int part_i, int part_j)
{
    for (int i = 0; i < clusterSize; i++)
        if (fastIndex (i, part_i) == 1 && fastIndex (i, part_j) == 1)
            return true;
    return false;
}


BBI DSMClusteringChromosome::getBBI ()
{
    BBI
        bbi (partN);
    int
        p,
        c;
    int
        clusterCount = 0;
    int
        partCount = 0;
    bool
        used[partN];
    for (p = 0; p < partN; p++)
        used[p] = false;
    for (c = 0; c < partN; c++)
    {
        partCount = 0;
        for (p = 0; p < partN; p++)
        {
            if (partInCluster (p, c))
            {
                bbi.bb[clusterCount][partCount + 1] = p;
                partCount++;
            }
        }
        if (partCount > 1)
        {
            bbi.bb[clusterCount][0] = partCount;
            clusterCount++;
        }
    }

    for (c = 0; c < clusterCount; c++)
        for (p = 1; p <= bbi.bb[c][0]; p++)
            used[bbi.bb[c][p]] = true;
    for (p = 0; p < partN; p++)
    {
        if (!used[p])
        {
            bbi.bb[clusterCount][0] = 1;
            bbi.bb[clusterCount][1] = p;
            clusterCount++;
        }
    }

    bbi.bbNum = clusterCount;
    return bbi;
}


void DSMClusteringChromosome::createDSMCC (Chromosome* population, int populationSize, int selectionPressure, int *selectionIndex)
{
    int i, j, k, m;
    Statistics stLow, stHigh;

    FastCounting *fastCounting =  new FastCounting[partN];

    for (i = 0; i < partN; i++)
        fastCounting[i].init(populationSize);

    for (i = 0; i < populationSize; i++)
        for (j = 0; j < partN; j++)
            fastCounting[j].setVal(i, population[selectionIndex[i]].getVal(j));

    // Randomly pick 5 MI and choose the median
    double mi[5];
    double q00[5],q01[5],q10[5],q11[5];
    int randI[5], randJ[5];
    // Can be too long for small ell
    assert(partN > 10);
    myRand.uniformArray(randI, 5, 0, partN-1);
    myRand.uniformArray(randJ, 5, 0, partN-1);
    for (m=0; m<5; m++)
    {
        i = randI[m];
        j = randJ[m];
        int n00, n01, n10, n11;
        n00 = n01 = n10 = n11 = 0;

        for (k = 0; k < fastCounting[0].lengthLong; k++)
        {
            unsigned long val1 = fastCounting[i].gene[k];
            unsigned long val2 = fastCounting[j].gene[k];

            unsigned long long01 = (~val1) & (val2);
            unsigned long long10 = (val1) & (~val2);
            unsigned long long11 = (val1) & (val2);

            n01 += myBD.countOne(long01);
            n10 += myBD.countOne(long10);
            n11 += myBD.countOne(long11);
        }

        n00 = populationSize - n01 - n10 - n11;

        q00[m] = (double)n00/(double)populationSize;
        q01[m] = (double)n01/(double)populationSize;
        q10[m] = (double)n10/(double)populationSize;
        q11[m] = (double)n11/(double)populationSize;

        mi[m] = mutualInformation(q00[m], q01[m], q10[m], q11[m]);
    }

    // O(n^2) for finding median, Okay since n = 5
    int median = 0;
    for (i=0; i<5; i++)
    {
        int smallerCount = 0;
        int greaterCount = 0;
        for (j=0; j<5; j++)
        {
            if (i==j) continue;
            if (mi[j] < mi[i]) smallerCount++;
            else if (mi[j] > mi[i]) greaterCount++;
        }
        if (abs(smallerCount-greaterCount) <= 1)
        {
            median = i;
            break;
        }
    }

    double independentMean, independentVariance;
    computeIMeanVar(&independentMean, &independentVariance, q00[median], q01[median], q10[median], q11[median], mi[median], populationSize);
    double split = betaSplit(independentMean, independentVariance);
    //printf("pass 0: %e\n", split);

    for (i = 0; i < partN; i++)
    {
        for (j = i+1; j < partN; j++)
        {

            int n00, n01, n10, n11;
            n00 = n01 = n10 = n11 = 0;

            for (k = 0; k < fastCounting[0].lengthLong; k++)
            {
                unsigned long val1 = fastCounting[i].gene[k];
                unsigned long val2 = fastCounting[j].gene[k];

                unsigned long long01 = (~val1) & (val2);
                unsigned long long10 = (val1) & (~val2);
                unsigned long long11 = (val1) & (val2);

                n01 += myBD.countOne(long01);
                n10 += myBD.countOne(long10);
                n11 += myBD.countOne(long11);

            }

            n00 = populationSize - n01 - n10 - n11;

            double p00, p01, p10, p11;
            p00 = (double)n00/(double)populationSize;
            p01 = (double)n01/(double)populationSize;
            p10 = (double)n10/(double)populationSize;
            p11 = (double)n11/(double)populationSize;

            double temp = mutualInformation(p00, p01, p10, p11);

            dsm.write(i,j,temp);

        }

    }
    delete []fastCounting;

/*
    for (int pass=1; pass<3; pass++)
    {
        Statistics q00, q01, q10, q11;
        Statistics mi;

        for (i = 0; i < partN; i++)
        {
            for (j = i+1; j < partN; j++)
            {

		double temp = dsm(i,j);
                if (temp > 0 && temp < split)
                    mi.record(temp);

            }

        }
        split = betaSplit(mi.getMean(), mi.getVariance());
        printf("pass %i: %e\n", pass, split);
    }
*/


    if (SHOW_DSM)
    {
        for (i = 0; i < partN; i++)
        {
            for (j = 0; j < partN; j++)
                if (i==j)
                    printf ("%0.3f ", 1.0);
            else if (i > j)
                printf ("%0.3f ", dsm (j, i));
            else
                printf ("%0.3f ", dsm (i, j));
            printf ("\n");
        }
        printf ("===========================\n");
        getchar();
    }

    for (i=0; i<partN; i++)
    {
        for (j=i+1; j<partN; j++)
        {
            double x = dsm(i,j);

            if (x > split)
                dsm.write(i,j,1.0);
            else
                dsm.write(i,j,0.0);

        }
    }

    if (SHOW_DSM)
    {
        for (i = 0; i < partN; i++)
        {
            for (j = 0; j < partN; j++)
                printf ("%0.1f ", dsm (i, j));
            printf ("\n");
        }
        getchar();
    }

}


void DSMClusteringChromosome::computeIMeanVar(double *mean, double *variance, double q00, double q01, double q10, double q11, double mi, int n)
{
    double J = mi;
    double q0x = q00+q01;
    double q1x = q10+q11;
    double qx0 = q00+q10;
    double qx1 = q01+q11;
    double K = 0.0;
    if (q00 > 1e-10) K+=q00 * square(log(q00/q0x/qx0));
    if (q01 > 1e-10) K+=q01 * square(log(q01/q0x/qx1));
    if (q10 > 1e-10) K+=q10 * square(log(q10/q1x/qx0));
    if (q11 > 1e-10) K+=q11 * square(log(q11/q1x/qx1));

    double M = 0.0;
    if (q00 > 1e-10) M+=(1-q00/q0x-q00/qx0+q00) * log(q00/q0x/qx0);
    if (q01 > 1e-10) M+=(1-q01/q0x-q01/qx1+q01) * log(q01/q0x/qx1);
    if (q10 > 1e-10) M+=(1-q10/q1x-q10/qx0+q10) * log(q10/q1x/qx0);
    if (q11 > 1e-10) M+=(1-q11/q1x-q11/qx1+q11) * log(q11/q1x/qx1);

    double Q = 1.0;
    if (q00 > 1e-10) Q-=q00*q00/q0x/qx0;
    if (q01 > 1e-10) Q-=q01*q01/q0x/qx1;
    if (q10 > 1e-10) Q-=q10*q10/q1x/qx0;
    if (q11 > 1e-10) Q-=q11*q11/q1x/qx1;

    *mean = 0.5/(n+1) + J;
    *variance = (K-J*J)/(n+1) + (M+0.5-J-Q)/(n+1)/(n+2);
}


double DSMClusteringChromosome::gaussianSplit(double mean, double variance)
{
    // epsilon <= 2 / ell^3
    // epsilon = 0.25 * exp(-z^2/2)

    double epsilon = 2.0 / partN / partN / partN;
    double z = sqrt( -2.0 * log(4 * epsilon));
    return mean + z * sqrt(variance);
}


double DSMClusteringChromosome::gammaSplit(double mean, double variance)
{

    double theta = variance / mean;
    double k = mean / theta;

    double gammaK = exp(gammaLn(k));

    double epsilon = 2.0 / partN / partN / partN;

    double leftZ = 0.0;
    double rightZ = 1.0;
    double middleZ;

    do
    {
        middleZ = (leftZ + rightZ ) / 2;
        double error = 1.0 - incompleteGamma(k, middleZ/theta) / gammaK;
        if (error > epsilon)
            leftZ = middleZ;
        else
            rightZ = middleZ;
    } while ( (rightZ - leftZ) > 1e-10* leftZ);

    return middleZ;
}


// gamma function
double DSMClusteringChromosome::gammaLn(double xx)
{
    double x,y,tmp,ser;
    static double cof[6]=
    {
        76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,
        0.1208650973866179e-2,-0.5395239384953e-5
    };
    int j;
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);

}


// gamma function
double DSMClusteringChromosome::incompleteGamma(double a, double x)
{
    assert(a > 0.0);
    assert(x >= 0.0);

    double sum = 0.0;
    double entity = 1.0 / a;

    for (int i=0; i<1000; i++)
    {
        sum += entity;
        entity *= x / (a+i+1);
    }

    return pow(x,a) * exp(-x) * sum;
}


double DSMClusteringChromosome::incompleteBeta(double a, double b, double x)
{
    float bt;
    assert (x >= 0.0 && x <= 1.0);
    if (x == 0.0 || x == 1.0) bt=0.0;
    else
        bt=exp(gammaLn(a+b)-gammaLn(a)-gammaLn(b)+a*log(x)+b*log(1.0-x));
    if (x < (a+1.0)/(a+b+2.0))
        return bt*betacf(a,b,x)/a;
    else
        return 1.0-bt*betacf(b,a,1.0-x)/b;
}


#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
double DSMClusteringChromosome::betacf(double a, double b, double x)
{
    int m,m2;
    float aa,c,d,del,h,qab,qam,qap;
    qab=a+b;
    qap=a+1.0;
    qam=a-1.0;
    c=1.0;
    d=1.0-qab*x/qap;
    if (fabs(d) < FPMIN) d=FPMIN;
    d=1.0/d;
    h=d;
    for (m=1;m<=MAXIT;m++)
    {
        m2=2*m;
        aa=m*(b-m)*x/((qam+m2)*(a+m2));
        d=1.0+aa*d;
        if (fabs(d) < FPMIN) d=FPMIN;
        c=1.0+aa/c;
        if (fabs(c) < FPMIN) c=FPMIN;
        d=1.0/d;
        h *= d*c;
        aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
        d=1.0+aa*d;
        if (fabs(d) < FPMIN) d=FPMIN;
        c=1.0+aa/c;
        if (fabs(c) < FPMIN) c=FPMIN;
        d=1.0/d;
        del=d*c;
        h *= del;
        if (fabs(del-1.0) < EPS) break;
    }
    assert (m <= MAXIT);
    return h;
}


double DSMClusteringChromosome::betaSplit(double mean, double variance)
{
    // alpha / beta

//    mean /= (1-log(2));
//    variance /= (1-log(2)) * (1-log(2));
    
    //modified by k
    mean /= log(2);
    variance /= log(2)*log(2);

    double ratio = mean / (1-mean);
    double beta = (ratio - variance*pow(ratio+1, 2)) / variance / pow(ratio+1, 3);
    double alpha = ratio * beta;
    double epsilon = 1.0 / partN / partN / partN;

    //printf("a, b = %e %e\n", alpha, beta);
    double leftZ = 0.0;
    double rightZ = 1.0;
    double middleZ;

    do
    {
        middleZ = (leftZ + rightZ ) / 2;
        // 1-I_x(a,b) = I_{1-x}(b,a);
        double error = incompleteBeta(beta, alpha, 1.0-middleZ);
        if (error > epsilon)
            leftZ = middleZ;
        else
            rightZ = middleZ;
        //printf ("=== z:%e  error:%e\n", middleZ, error);
    } while ( (rightZ - leftZ) > 1e-10* leftZ);

    return middleZ * (1-log(2)) ;                 // / log(2.0);

}
