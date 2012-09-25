/***************************************************************************
 *   Copyright (C) 2004 by tianliyu                                        *
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

#include <stdio.h>
#include <cstring>
#include "chromosome.h"
#include "UFO.h"
#include "Tsuji.h"
#include "fitOC.h"

double 
Chromosome::getOCFitness() const {
    int *x = new int[length];
    
    for ( int i = 0; i < length; ++i )
    {
        x[i] = getVal(i);
    }
    double result = ocFitness( x );
    delete []x;
    return result;
}

double 
Chromosome::getOCMaxFitness() const {
    return ocFitness.getMaxFitness();
}


double 
Chromosome::getMaxFitness () const
{
  // Never converge
  //return 1e6-1e-6;

  // For Trap
  //return (length/TRAP_K-0.11);	
  // For Trap5
  //return (length/DIS_M - 0.11);
  // For new Tsuji
  //return getTsujiMaxFitness();
  // for UFO
  //return getUFOMaxFitness();
  // For hTrap
  //return (length/TRAP_K) * int(log(length)/log(TRAP_K)) - 1e-6;	
  // For hXOR hIFF
  //return (length * (log2(length) + 1) - 1e-6);
  // For ONEMAX
  //return length - 1e-6;
  // For Over Conflict
  return getOCMaxFitness();  
  
}

double
Chromosome::evaluate ()
{
  evaluated = true;

  //double accum = hIFF ();
  //double accum = hXOR ();
  //double accum = hTrap ();
  //double accum = MKTrap1 (1, 0.9);
  
  //double accum = MKtrap5();   // for mune test problem

  //double accum = getUFOFitness();
  //double accum = getTsujiFitness();
  //return accum;
  
  // For Over Conflict
  return getOCFitness();
  //return oneMax ();
  
  //return MKHC ();

}



Chromosome::Chromosome ()
{
  length = 0;
  lengthLong = 0;
  gene = NULL;
  evaluated = false;
}

Chromosome::Chromosome (int n_length)
{
  gene = NULL;

  init (n_length);
}


Chromosome::~Chromosome ()
{
  if (gene != NULL) delete []gene;
}

void
Chromosome::init (int n_length)
{
  length = n_length;
  lengthLong = quotientLong(length)+1;

  if (gene != NULL)
      delete []gene;

  gene = new unsigned long [lengthLong];
  gene[lengthLong-1] = 0;

  evaluated = false;
}

unsigned int Chromosome::getInt (int start, int length) const
{
    assert (length < (int) sizeof(unsigned int) * 8);

    int q = quotientLong(start);
    int r = remainderLong(start);


    if ( r + length <= (int) sizeof(unsigned long) * 8) {
        // with in one long
        unsigned long mask = (~(0lu)) >> (sizeof(unsigned long) * 8 - length);
        mask <<= r;

        return (unsigned int) ((gene[q] & mask) >> r);
    }
    else {
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
Chromosome::getFitness ()
{
 
  if (evaluated)
    return fitness;
  else
    return (fitness = evaluate ());
}

bool Chromosome::isEvaluated () const
{
  return evaluated;
}

// OneMax
double
Chromosome::oneMax () const 
{
  int i;
  double result = 0;

  for (i = 0; i < length; i++)
    result += getVal(i);

  return result;
}

Chromosome& Chromosome::operator= (const Chromosome& c)
{

  if (length != c.length)
    {
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
Chromosome::getCorrectBBsMKTrap () const 
{
  if (length % TRAP_K != 0)
    {
      outputErrMsg ("TRAP_K doesn't divide length");
    }

  int i, j;
  int u;
  int correctBBNum = 0;

  int TRAP_M = length / TRAP_K;

  for (i = 0; i < TRAP_M; i++)
    {
      u = 0;
      for (j = 0; j < TRAP_K; j++)
	u += getVal(i * TRAP_K + j);

      if (u == TRAP_K)
	correctBBNum++;
    }

  return correctBBNum;

}

double
Chromosome::trap0 (int unitary, double fHigh, double fLow) const
{
  if (unitary > TRAP_K)
    return 0;

  return trap1 (TRAP_K - unitary, fHigh, fLow);
}

double Chromosome::trap1 (int unitary, double fHigh, double fLow) const
{
  if (unitary > TRAP_K)
    return 0;

  if (unitary == TRAP_K)
    return fHigh;
  else
    return fLow - unitary * fLow / (TRAP_K);
}

double 
Chromosome::trap5(int unitary) const
{
  if (unitary > 5)
  return 0;

  if (unitary == 5)
  return 1;
  else 
  return double((4-unitary)/double(5));


}

/** Currently two levels */
// Need be rewrite
double Chromosome::hTrap () const
{
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
		nextLevel[i] = 20;	// trap0 and trap won't count this

	      if (localLength >= TRAP_K * TRAP_K) 
		      result += factor * trap1(u, 1.0, 1.0);
	      else // last level
		      result += factor * trap1(u, 1.0, 0.9);
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

double
Chromosome::MKTrap0 (double fHigh, double fLow) const
{
  int i, j;
  int u;

  int TRAP_M = length / TRAP_K;

  if (length % TRAP_K != 0)
    outputErrMsg ("TRAP_K doesn't divide length");

  double result = 0;

  for (i = 0; i < TRAP_M; i++)
    {
      u = 0;
      for (j = 0; j < TRAP_K; j++)
	u += getVal(i * TRAP_K + j);

      result += trap0 (u, fHigh, fLow);
    }

  return result;
}

double
Chromosome::MKTrap1 (double fHigh, double fLow) const
{
  int i, j;
  int u;

  int TRAP_M = length / TRAP_K;

  if (length % TRAP_K != 0)
    outputErrMsg ("TRAP_K doesn't divide length");

  double result = 0;

  for (i = 0; i < TRAP_M; i++)
    {
      u = 0;
      for (j = 0; j < TRAP_K; j++)
	u += getVal(i * TRAP_K + j);

      result += trap1 (u, fHigh, fLow);
    }

  return result;
}

double 
Chromosome::MKtrap5 () const
{
  int u;
  int TRAP_M = length / DIS_M;

  if (length % DIS_M != 0)
    outputErrMsg("TRAP_K does not divide length");

  double result = 0;
  
  
  for (int i = 0; i < TRAP_M ; i++)
  {
    u = 0;
    //getchar();
    for (int j = 0; j < GROUP_SIZE ;j++)
    {
      u += getVal(decisionV[i*GROUP_SIZE+j]);
      //printf("i,j,idx = %d %d %d\n",i,j,decisionV[i*5+j]);
    }
    
    result += trap5(u);
    //printf("(i,result) = %d %f \n",i,result);
  }
  
  return result;
}


double 
Chromosome::getTsujiFitness() const
{
     /*if( length != ell )
     {
       cout << "Length does not match!" << endl;  
       return 0.0;                        
     } */                 
     double result = 0.0;
     int unitary = 0;
     for ( int i = 0; i < tsuji.mym; ++i )
     {
        unitary = 0;
        for ( int j = 0; j < tsuji.myk; ++j )
        {
            if ( getVal(tsuji.bbi.bb[i][j+1]) == 1 )
                ++unitary;
        }
        result += trap5(unitary);
     }
     return result;
}



double
Chromosome::MKHC () const
{
  int i, j;
  int u;

  int TRAP_M = length / TRAP_K;

  if (length % TRAP_K != 0)
    outputErrMsg ("TRAP_K doesn't divide length");

  double result = 0;

  for (i = 0; i < TRAP_M; i++)
    {
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
Chromosome::printOut () const
{
  int i;
  for (i = 0; i < length; i++)
    printf ("%d", getVal(i));
}

void Chromosome::shortPrintOut () const
{
  int i, j;
  int u;

  int TRAP_M = length / TRAP_K;

  for (i = 0; i < TRAP_M; i++)
    {
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
Chromosome::makeInt (int *bb) const
{
  int value = 0;

  for (int i = 1; i <= bb[0]; i++)
    if (getVal(bb[i]) == 1)
      value |= (1 << (i - 1));

  return value;
}

int Chromosome::getLength () const
{
  return length;
}

double Chromosome::hIFF () const
{
  return (double) hIFFBigH(0, length, true);
}

int Chromosome::hIFFBigH (int start, int hIFFLength, bool lookUpTable) const
{

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

bool Chromosome::hIFFSmallH (int start, int hIFFLength, bool lookUpTable) const
{
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



double Chromosome::hXOR () const
{
  return (double) hXORBigH(0, length, true);
}

int Chromosome::hXORBigH (int start, int hXORLength, bool lookUpTable) const
{

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

bool Chromosome::hXORSmallH (int start, int hXORLength, bool lookUpTable) const
{
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


double Chromosome::getUFOFitness() const {
    return ufo.getFitness(*this);
}

double Chromosome::getTsujiMaxFitness() const {
    return tsuji.getMaxFitness();
}



double Chromosome::getUFOMaxFitness() const {
    return ufo.getMaxFitness();
}
