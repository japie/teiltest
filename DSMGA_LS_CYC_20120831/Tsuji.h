#ifndef TSUJI_H
#define TSUJI_H


#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include "statistics.h"
#include "bbi.h"
#include "myrand.h"

using namespace std;

#define SHOW_GEN2BB true

class Tsuji{
  friend class Chromosome;
  public:
    Tsuji();
    ~Tsuji();

    void sampling(double w, int k, int ell, double std, int fCount);
    void verbose();
    BBI getBBI(double w, int k, int ell, double std, int fCount);
    BBI getBBI();
    //double getFitness( const Chromosome& c) const;
    double getMaxFitness();
    void setParams(double w, int k, int ell, double std);
    int getM();
  protected:
    bool loadSample(double w, int k, int ell, double std, int fCount);
    void reset();
    double checkRange(double d, int ell);
    void mySplit(const string &line,const string &delimiter);
    BBI     bbi;
    bool    ver;
    MyRand  myRand;
    int**   gen2BB;
    vector<int> tokens;
    double  myw;
    int     myk;
    int     myell;
    double  mystd;
    int     mym;
    

};

#endif
