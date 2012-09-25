/***************************************************************************
 *   Copyright (C) 2011 by TEIL                                        *
 *                                                                         *
 ***************************************************************************/
#ifndef _CHROMOSOME_H
#define _CHROMOSOME_H

#include "global.h"
#include "UFO.h"

class Chromosome
{

    friend class GA;
    friend class DSMGA;
    friend class TableLookUp;

public:

    static enum Function {
        ONEMAX,
        MKTRAP,
        HIFF,
        HXOR,
        HTRAP,
        TRAP3339,
				MKCYCLIC,
        SPINGLASS,
        TSUJI,
        UFO,
        FITOC
    } function;


    Chromosome ();
    Chromosome (int n_ell);

    ~Chromosome ();

    friend void crossover (Chromosome & p1, Chromosome & p2, Chromosome & c1, Chromosome & c2);

    void init (int n_ell);

    int getVal (int index) const {
        assert (index >= 0 && index < length);

        int q = quotientLong(index);
        int r = remainderLong(index);

        if ( (gene[q] & (1lu << r)) == 0 )
            return 0;
        else
            return 1;
    }

    void setVal (int index, int val) {
        assert (index >= 0 && index < length);

        //outputErrMsg ("Index overrange in Chromosome::operator[]");

        int q = quotientLong(index);
        int r = remainderLong(index);

        if (val == 1)
            gene[q] |= (1lu<<r);
        else
            gene[q] &= ~(1lu<<r);

        evaluated = false;
    }

    /** real evaluator */
    double evaluate ();

    bool isEvaluated () const;

    bool operator== (const Chromosome & c) const;
    Chromosome & operator= (const Chromosome & c);

    int makeInt (int *bb) const;
	long int makeInt (const vector<int>& bb) const;
    int steepestDescent();

    double getFitness ();
    double oneMax () const;
    double MKTrap0 (double high, double low) const;
    double MKTrap1 (double high, double low) const;
    double MKCyclic(double fHigh, double fLow) const;
    double trap0 (int u, double high, double low, int trapK) const;
    double trap1 (int u, double high, double low, int trapK) const;
    double hTrap () const;
    double MKHC () const;
    double Trap3339 () const;
    double spinGlass () const;
    double TsujiFitness () const;
    double UFOFitness () const;
    double OCFitness() const;
    
    double hIFF () const;
    int hIFFBigH (int, int, bool) const;
    bool hIFFSmallH (int, int, bool) const;

    double hXOR () const;
    int hXORBigH (int, int, bool) const;
    bool hXORSmallH (int, int, bool) const;

    int getCorrectBBsMKTrap () const;

    void printOut () const;
    void shortPrintOut () const;

    int getLength () const;
    void setLength ();

    double getMaxFitness () const;
protected:

    unsigned int getInt(int start, int length) const;

    unsigned long *gene;
    int length;
    int lengthLong;
    double fitness;
    bool evaluated;

};
#endif
