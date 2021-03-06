/***************************************************************************
 *   Copyright (C) 2011 by TEIL                                            *
 ***************************************************************************/

#include <cstdio>
#include <cstring>
#include "spin-glass.h"
#include "chromosome.h"


//Chromosome::Function Chromosome::function = UFO;
//Chromosome::Function Chromosome::function = FITOC;
//Chromosome::Function Chromosome::function = MKTRAP;
//Chromosome::Function Chromosome::function = SPINGLASS;
//Chromosome::Function Chromosome::function = FITOC;
Chromosome::Function Chromosome::function = MKCYCLIC;


Chromosome::Chromosome () {
    length = 0;
    lengthLong = 0;
    gene = NULL;
    evaluated = false;
}

Chromosome::Chromosome (int n_length) {
    gene = NULL;

    init (n_length);
}


Chromosome::~Chromosome () {
    if (gene != NULL) delete []gene;
}

void
Chromosome::init (int n_length) {
    length = n_length;
    lengthLong = quotientLong(length)+1;

    if (gene != NULL)
        delete []gene;

    gene = new unsigned long [lengthLong];
    gene[lengthLong-1] = 0;

    evaluated = false;
}

unsigned int Chromosome::getInt (int start, int length) const {
    assert (length < (int) sizeof(unsigned int) * 8);

    int q = quotientLong(start);
    int r = remainderLong(start);


    if ( r + length <= (int) sizeof(unsigned long) * 8) {
        // with in one long
        unsigned long mask = (~(0lu)) >> (sizeof(unsigned long) * 8 - length);
        mask <<= r;

        return (unsigned int) ((gene[q] & mask) >> r);
    } else {
        unsigned long mask = (~(0lu)) << r;
        int part1 = ((gene[q] & mask) >> r);

        int length2 = (r + length) - sizeof(unsigned long) * 8;
        int length1 = length - length2;

        mask = (~(0lu)) >> (sizeof(unsigned long) * 8 - length2);
        int part2 = (gene[q+1] & mask);

        return ((part2 << length1) & part1);
    }

}

double
Chromosome::getFitness () {
    if (evaluated)
        return fitness;
    else
        return (fitness = evaluate ());
}

int Chromosome::steepestDescent() {

    int index;
    double origF;
    int nfe = 0;

    do {
        origF = getFitness();
        double maxF = origF;
        index = -1;

        for (int i=0; i<length; ++i) {
            setVal(i, 1-getVal(i));  // try
            double modified = getFitness();
            if (modified > maxF) {
                maxF = modified;
                index = i;
            }
            setVal(i, 1-getVal(i));  // undo
        }

        nfe += length;

        if (index != -1) {
            setVal(index, 1-getVal(index));
            evaluated = true;
            fitness = maxF;
        }

    } while (index != -1);

    fitness = origF;
    evaluated = true;
    return nfe;
}

bool Chromosome::isEvaluated () const {
    return evaluated;
}

double
Chromosome::evaluate () {
    evaluated = true;
    double accum = 0.0;

    switch (function) {
    case ONEMAX:
        accum = oneMax();
        break;
    case MKTRAP:
        accum = MKTrap1(1, 0.9);
        break;
    case MKCYCLIC:
				accum = MKCyclic(1, 0.9);
				break;
    case HIFF:
        accum = hIFF();
        break;
    case HXOR:
        accum = hXOR();
        break;
    case HTRAP:
        accum = hTrap();
        break;
    case SPINGLASS:
        accum = spinGlass();
        break;
    case TRAP3339:
        accum = Trap3339();
        break;
    case TSUJI:
        accum = TsujiFitness();
        break;
    case UFO:
        accum = UFOFitness();
        break;
    case FITOC:
        accum = OCFitness();
        break;
    default:
        accum = MKTrap1(1, 0.9);
        break;
    }


    return accum;

}



double
Chromosome::spinGlass () const {

    int *x = new int[length];
    double result;

    for (int i=0; i<length; i++)
        if (getVal(i) == 1)
            x[i] = 1;
        else
            x[i] = -1;

    result = spinGlassValue(x, &mySpinGlassParams);

    delete []x;

    return result;
}

double
Chromosome::TsujiFitness() const {
    /*if( length != ell )
     {
     cout << "Length does not match!" << endl;
     return 0.0;
     } */
    double result = 0.0;
    int unitary = 0;
    for ( int i = 0; i < tsuji.mym; ++i ) {
        unitary = 0;
        for ( int j = 0; j < tsuji.myk; ++j ) {
            if ( getVal(tsuji.bbi.bb[i][j+1]) == 1 )
                ++unitary;
        }
        result += trap1(unitary,1.0,0.9,TRAP_K);
    }
    return result;
}

double Chromosome::OCFitness () const {

    int *x = new int[length];

    for ( int i = 0; i < length; ++i ) {
        x[i] = getVal(i);
    }
    double result = fitoc( x );
    delete []x;
    return result;
}


double Chromosome::UFOFitness () const {

    int *x = new int[length];

    for ( int i = 0; i < length; ++i ) {
        x[i] = getVal(i);
    }
    double result = ufo.getFitness( x );
    delete []x;
    return result;
}


// OneMax
double
Chromosome::oneMax () const {
    int i;
    double result = 0;

    for (i = 0; i < length; i++)
        result += getVal(i);

    return result;
}

bool Chromosome::operator== (const Chromosome& c) const {
    if (length != c.length)
        return false;

    for (int i=0; i<lengthLong; i++)
        if (gene[i] != c.gene[i])
            return false;

    return true;
}

Chromosome& Chromosome::operator= (const Chromosome& c) {

    if (length != c.length) {
        length = c.length;
        init (length);
    }

    evaluated = c.evaluated;
    fitness = c.fitness;
    lengthLong = c.lengthLong;

    memcpy(gene, c.gene, sizeof(long) * lengthLong);

    return *this;
}

int
Chromosome::getCorrectBBsMKTrap () const {
    if (length % TRAP_K != 0) {
        outputErrMsg ("TRAP_K doesn't divide length");
    }

    int i, j;
    int u;
    int correctBBNum = 0;

    int TRAP_M = length / TRAP_K;

    for (i = 0; i < TRAP_M; i++) {
        u = 0;
        for (j = 0; j < TRAP_K; j++)
            u += getVal(i * TRAP_K + j);

        if (u == TRAP_K)
            correctBBNum++;
    }

    return correctBBNum;

}

double Chromosome::trap0 (int unitary, double fHigh, double fLow, int trapK) const {
    if (unitary > trapK)
        return 0;

    return trap1 (trapK - unitary, fHigh, fLow, trapK);
}

double Chromosome::trap1 (int unitary, double fHigh, double fLow, int trapK) const {
    if (unitary > trapK)
        return 0;

    if (unitary == trapK)
        return fHigh;
    else
        return fLow - unitary * fLow / (trapK-1);
}


/** Currently two levels */
// Need be rewrite
double Chromosome::hTrap () const {
    int i, j;
    int u;

    if (length % TRAP_K != 0)
        outputErrMsg ("TRAP_K doesn't divide length");


    double result = 0.0;


    int *localCh = new int[length];
    int *nextLevel = new int[length / TRAP_K];

    int localLength = length;

    for (i = 0; i < length; i++)
        localCh[i] = getVal(i);

    int factor = 1;
    while (localLength >= TRAP_K) {

        int TRAP_M = localLength / TRAP_K;

        for (i = 0; i < TRAP_M; i++) {
            u = 0;

            for (j = 0; j < TRAP_K; j++)
                u += localCh[i * TRAP_K + j];

            if (u == TRAP_K)
                nextLevel[i] = 1;
            else if (u == 0)
                nextLevel[i] = 0;
            else
                nextLevel[i] = 20; // trap0 and trap1 won't count this

            if (localLength >= TRAP_K * TRAP_K)
                result += factor * trap1(u, 1.0, 1.0, TRAP_K);
            else // last level
                result += factor * trap1(u, 1.0, 0.9, TRAP_K);
        }

        factor *= TRAP_K;
        localLength /= TRAP_K;

        for (i = 0; i < localLength; i++)
            localCh[i] = nextLevel[i];
    }

    delete []localCh;
    delete []nextLevel;

    return result;

}

double Chromosome::Trap3339() const {
    double result = 0.0;

    int m;

    m = length / 3;
    for (int i=0; i<m; i++) {
        int u=0;
        for (int j=0; j<3; j++)
            u += getVal(i*3 + j);
        result += trap1 (u, 1, 0.9, 3);
    }

    m = length / 9;
    for (int i=0; i<m; i++) {
        int u=0;
        for (int j=0; j<9; j++)
            u += getVal(i*9 + j);
        result += trap1 (u, 1, 0.9, 9);
    }

    return result;

}

double Chromosome::MKTrap0 (double fHigh, double fLow) const {
    int i, j;
    int u;

    int trapM = length / TRAP_K;

    if (length % TRAP_K != 0)
        outputErrMsg ("TRAP_K doesn't divide length");

    double result = 0;

    for (i = 0; i < trapM; i++) {
        u = 0;
        for (j = 0; j < TRAP_K; j++)
            u += getVal(i * TRAP_K + j);

        result += trap0 (u, fHigh, fLow, TRAP_K);
    }

    return result;
}

double
Chromosome::MKCyclic(double fHigh, double fLow) const {
    int i, j;
    int u;
    int TRAP_M = length / (TRAP_K-1);
    if (length % (TRAP_K-1) != 0)
      outputErrMsg ("TRAP_k doesn't divide length for Cyclic Setting");
    double result = 0;
    for (i = 0; i < TRAP_M; i++) {
      u = 0;
      int idx = i * TRAP_K - i;
      for (j = 0; j < TRAP_K; j++){
	int pos = idx + j;
        if (pos == length)
	  pos = 0;
        else if (pos > length)
          outputErrMsg ("CYCLIC BUG");
        //
	u += getVal(pos);
      }
      result += trap1 (u, fHigh, fLow, TRAP_K);
    }
    return result;
}





double
Chromosome::MKTrap1 (double fHigh, double fLow) const {
    int i, j;
    int u;

    int TRAP_M = length / TRAP_K;

    if (length % TRAP_K != 0)
        outputErrMsg ("TRAP_K doesn't divide length");

    double result = 0;

    for (i = 0; i < TRAP_M; i++) {
        u = 0;
        for (j = 0; j < TRAP_K; j++)
            u += getVal(i * TRAP_K + j);

        result += trap1 (u, fHigh, fLow, TRAP_K);
    }

    return result;
}

double
Chromosome::MKHC () const {
    int i, j;
    int u;

    int TRAP_M = length / TRAP_K;

    if (length % TRAP_K != 0)
        outputErrMsg ("TRAP_K doesn't divide length");

    double result = 0;

    for (i = 0; i < TRAP_M; i++) {
        u = 0;
        for (j = 0; j < TRAP_K; j++)
            u += getVal(i * TRAP_K + j);

        if (u == 4)
            result += 1.0;
        else if (u == 3)
            result += 0.0;
        else if (u == 2)
            result += 0.0;
        else if (u == 1)
            result += 0.9;
        else
            result += 0.0;

    }

    return result;
}

void
Chromosome::printOut () const {
    int i;
    for (i = 0; i < length; i++)
        printf ("%d", getVal(i));
}

void Chromosome::shortPrintOut () const {
    int i, j;
    int u;

    int TRAP_M = length / TRAP_K;

    for (i = 0; i < TRAP_M; i++) {
        u = 0;
        for (j = 0; j < TRAP_K; j++)
            u += getVal(i * TRAP_K + j);

        if (u == TRAP_K)
            ::printf ("1");
        else if (u == 0)
            ::printf ("0");
        else
            ::printf ("*");

    }
}

int
Chromosome::makeInt (int *bb) const {
    int value = 0;

    for (int i = 1; i <= bb[0]; i++)
        if (getVal(bb[i]) == 1)
            value |= (1 << (i - 1));

    return value;
}

long int Chromosome::makeInt (const vector<int>& bb) const
{
	long int value = 0;
	/*
	 cout <<endl;
	 for ( int i = 0; i < bb.size(); ++i )
	 {
	 cout << bb[i] << "  ";
	 }
	 cout << endl;
	 */
	for ( int i = 0; (int)i < bb.size(); ++i )
	{
            if (getVal(bb[i]) == 1)
            {
                value |= (1 << (i));
                //cout << "value = " << value << endl;
            }
	}

	return value;
}

int Chromosome::getLength () const {
    return length;
}

double Chromosome::getMaxFitness () const {

    switch (function) {
    case ONEMAX:
        return length - 1e-6;
    case MKTRAP:
        return (length/TRAP_K-1e-6);
    case MKCYCLIC:
				return (length/(TRAP_K - 1) - 1e-6);
    case HTRAP:
        return (length/TRAP_K) * int(log(length)/log(TRAP_K)) - 1e-6;
    case HXOR:
    case HIFF:
        return (length * (log2(length) + 1) - 1e-6);
    case TRAP3339:
        return (length/9)*4-1e-6;
    case SPINGLASS:
        return mySpinGlassParams.optimalValue - 1e-6;
    case TSUJI:
        return tsuji.getMaxFitness() - 1e-6;
    case UFO:
        return ufo.getMaxFitness() -1e-6;
    case FITOC:
        return fitoc.getMaxFitness() -1e-6;
    default:
        // Never converge
        return 1e6-1e-6;
    }
}

double Chromosome::hIFF () const {
    return (double) hIFFBigH(0, length, true);
}

int Chromosome::hIFFBigH (int start, int hIFFLength, bool lookUpTable) const {

    if (hIFFLength == 1) return 1;

    if (lookUpTable && hIFFLength <=16)
        return myTableLookUp.lookUpInt(getInt(start, hIFFLength));

    int startL = start;
    int startR = start + hIFFLength/2;

    if (hIFFSmallH(start, hIFFLength, lookUpTable))
        return hIFFBigH(startL, hIFFLength/2, lookUpTable) + hIFFBigH(startR, hIFFLength/2, lookUpTable) + hIFFLength;
    else
        return hIFFBigH(startL, hIFFLength/2, lookUpTable) + hIFFBigH(startR, hIFFLength/2, lookUpTable);
}

bool Chromosome::hIFFSmallH (int start, int hIFFLength, bool lookUpTable) const {
    int i;

    if (hIFFLength == 1) return true;

    if (lookUpTable && hIFFLength <=16)
        return myTableLookUp.lookUpBool(getInt(start, hIFFLength));

    int startL = start;
    int startR = start + hIFFLength/2;

    bool iff = true;
    for (i=0; i<hIFFLength/2; i++)
        if (getVal(startL+i) != getVal(startR+i)) iff = false;

    if ( iff && hIFFSmallH(startL, hIFFLength/2, lookUpTable) && hIFFSmallH(startR, hIFFLength/2, lookUpTable))
        return true;

    return false;
}



double Chromosome::hXOR () const {
    return (double) hXORBigH(0, length, true);
}

int Chromosome::hXORBigH (int start, int hXORLength, bool lookUpTable) const {

    if (lookUpTable && hXORLength <= 16)
        return myTableLookUp.lookUpInt(getInt(start, hXORLength));

    if (hXORLength == 1) return 1;

    int startL = start;
    int startR = start + hXORLength/2;

    if (hXORSmallH(start, hXORLength, lookUpTable))
        return hXORBigH(startL, hXORLength/2, lookUpTable) + hXORBigH(startR, hXORLength/2, lookUpTable) + hXORLength;
    else
        return hXORBigH(startL, hXORLength/2, lookUpTable) + hXORBigH(startR, hXORLength/2, lookUpTable);
}

bool Chromosome::hXORSmallH (int start, int hXORLength, bool lookUpTable) const {
    int i;

    if (lookUpTable && hXORLength <= 16)
        return myTableLookUp.lookUpBool(getInt(start, hXORLength));

    if (hXORLength == 1) return true;

    int startL = start;
    int startR = start + hXORLength/2;

    bool isComplement = true;
    for (i=0; i<hXORLength/2; i++)
        if (getVal(startL+i) == getVal(startR+i)) isComplement = false;

    if ( isComplement && hXORSmallH(startL, hXORLength/2, lookUpTable) && hXORSmallH(startR, hXORLength/2, lookUpTable))
        return true;

    return false;
}



