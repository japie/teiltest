/*
 * File:   UFO.h
 * Author: Rick
 *
 * Created on January 17, 2011, 7:45 PM
 */

#ifndef UFO_H
#define UFO_H


#include <cstdlib>
#include <cstdio>
#include <vector>
#include <list>
#include "myrand.h"
#include "statistics.h"
#include "bbi.h"
//#include "chromosome.h"


//#ifdef CHECK_GLOBAL_WHEN_CALCULATATION


using std::vector;
using std::list;

class UFO
{
public:
    UFO();
    UFO(const int& m, const int& k, const double& desiredW);
    UFO(const UFO& orig);
    virtual ~UFO();

    bool createProblem(const int& n, const int& k, const double& desiredW, int conflict);


    void print();
    void printStatistics();
    BBI getBBI();
    int getNumOfGenes() const;
    int getNumOfBBs() const;

    double getFitness( int const * chromosome ) const;
    
    double getMaxFitness() const;

    int getZeros() const;

    double getMeanOmega();
    double getMaxOmega();
    double getMinOmega();
    double getG2BStd();
    double getB2BStd();
    double getB2BMean();
    int getG2BMaxDifference();

    int getUnusedGenes();

    bool checkGlobalOptimum();
    bool checkAndSetGlobalOptimum();

    int getTotalConflictedBB();



private:
    void reset();
    bool initial(const int& n, const int& k, const double& desiredW);
    bool construct(const int& n, const int& k, const double& desiredW);
    void constructGraph();
    bool BFSandUpdate();
    void traceBackAndUpdate(int* parent);
    void evaluate();

    double kTrap( const int& k, const int& totalOne, const int& type) const;
    double kTrap9(const int& k, const int& totalOne, const int& type) const;
    double kMax(const int& k, const int& totalOne, const int& type) const;
    double ikTrap( const int& k, const int& totalOne, const int& type) const;


    int conflict();
    int conflict2( const int& totalConflict);
    int conflict3( const int& totalConflict);
    int conflict4( const int& totalConflict);

    int randomConflict( const int& totalConflict);

    int findConflictedBB();


//    int findSmallestBBAndUpdate( vector<vector<int> >& list);
//    int findSmallestGeneAndUpdate( vector<vector<int> >& list);

private:
    int** bb2gene;
    int** gene2bb;
    int* bbSize;
    int* geneSize;
    int* geneOmega;
    int m;
    int n;
    int k;
    int w;
    double desiredW;

    list<int>* graph;
    int** room;
    MyRand r;
    Statistics g2B;
    Statistics b2B;
    bool evaluated;
    BBI bbi;

    int* bbType;
    int zeros;
    int unUsedGenes;

    double globalOptimum;
    int* bestChromosome;
    int conflictedBB;

};

#endif /* UFO_H */

