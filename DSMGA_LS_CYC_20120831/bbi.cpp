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
#include "bbi.h"
#include <vector>
#include <cstdlib>

using std::vector;

BBI::BBI ()
{
  bb = NULL;
  ell = 0;
  prepareBB = 0;
}

BBI::BBI (int n_ell)
{
  init (n_ell);
}

BBI::BBI (const BBI & bbi)
{
  int i, j;

//    //printf ("1why");
//    if ( bb != NULL )
//    {
//        for (i = 0; i < prepareBB; i++)
//            delete[](bb[i]);
//        delete[]bb;
//    }

   init (bbi.ell);

    //printf ("2why");
    for (i = 0; i < ell; i++)
        for (j = 0; j < ell + 1; j++)
            bb[i][j] = bbi.bb[i][j];

    //printf ("3why");
    bbNum = bbi.bbNum;


}

void
BBI::init (int n_ell)
{
  int i;

  ell = n_ell;

  prepareBB = 4*ell;

  bb = new int *[prepareBB];

  for (i = 0; i < prepareBB; i++)
    bb[i] = new int[ell + 1];

  for (i = 0; i < prepareBB; i++)
    bb[i][0] = 0;

  bbNum = 0;
}

BBI & BBI::operator= (const BBI & bbi)
{
  int
    i,
    j;

  if (ell != bbi.ell)
    {
      for (i = 0; i < ell; i++)
	delete[](bb[i]);

      delete[]bb;
      init (bbi.ell);
    }

  for (i = 0; i < ell; i++)
    for (j = 0; j < ell + 1; j++)
      bb[i][j] = bbi.bb[i][j];

  bbNum = bbi.bbNum;

  return *this;
}

BBI::~BBI ()
{
  int i;

  if (bb != NULL)
    {
      for (i = 0; i < 4*ell; i++)
	delete[](bb[i]);

      delete[]bb;
    }

  ell = 0;
  bb = NULL;
}

void
BBI::printOut ()
{
  int i, j;

  for (i = 0; i < ell; i++)
    {
      if (bb[i][0] != 0)
	{

	  printf ("[ %i", bb[i][1]);

	  for (j = 2; j <= bb[i][0]; j++)
	    printf ("-%i", bb[i][j]);

	  printf (" ]");
	}
    }

  printf ("\n");

}

bool
BBI::operator!= (const BBI& bbi) const
{
  return (!((*this) == bbi));
}

bool
BBI::operator== (const BBI& bbi) const
{

  if (bbNum != bbi.bbNum) return false;

  bool *matched = new bool [bbNum];

  int i,j,k;


  for (i=0; i<bbNum; i++)
    matched[i] = false;


  for (i=0; i<bbNum; i++) {

    bool foundMatch = false;

    for (j=0; j<bbNum; j++) {
      if (bb[i][0] != bbi.bb[j][0]) continue;
      if (matched[j]) continue;

      foundMatch = true;
      for (k=1; k<=bb[i][0]; k++)
	if (bb[i][k] != bbi.bb[j][k]) {
	  foundMatch = false;
	  break;
	}

      if (foundMatch) {
        matched[j] = true;
        break;
      }
    }

    if (!foundMatch) {
	    delete []matched;
	    return false;
    }

  }

  delete []matched;
  return true;

}

int BBI::getNumofMatchedBB (const BBI& bbi) const
{
  bool *matched = new bool [bbi.bbNum];

  int i,j,k;
  int counter = 0;


  for (i=0; i<bbi.bbNum; i++)
    matched[i] = false;

    
  for (i=0; i<bbNum; i++) {
    
    bool foundMatch = false;

    for (j=0; j<bbi.bbNum; j++) {
      if (bb[i][0] != bbi.bb[j][0]) continue;
      if (matched[j]) continue;

      foundMatch = true;
      for (k=1; k<=bb[i][0]; k++)
	if (bb[i][k] != bbi.bb[j][k]) {
	  foundMatch = false;
	  break;
	}

      if (foundMatch) {
        matched[j] = true;
        break;
      }
    }

    if (foundMatch) counter++;

  }

  delete []matched;
  return counter;
}


int BBI::calculateBBIError(const BBI & bbi, int & correct, int& miss ) const {
    
    int ** geneMatrix = new int*[ell];
    int ** correctMatrix = new int*[ell];
    int ** correctCount = new int*[ell];
    int ** errorCount = new int*[ell];
    for ( int i = 0; i < ell; ++i )
    {
        geneMatrix[i] = new int[ell];
        correctMatrix[i]= new int[ell];
        correctCount[i] = new int[ell];
        errorCount[i] = new int[ell];
    }
    
    for ( int i = 0; i < ell; ++i )
    {
        for ( int j = 0; j < ell; ++j )
        {
            geneMatrix[i][j] = 0;
            correctMatrix[i][j]= 0;
            correctCount[i][j] = 0;
            errorCount[i][j] = 0;
        }
    }
    
    
    for ( int i = 0; i < bbNum; ++i )
    {
        for ( int j = 1; j < bb[i][0] +1; ++j )
        {
            int thisGene = bb[i][j];
            int neiGene = 0;
            for ( int k = j +1; k < bb[i][0] +1; ++k )
            {
                neiGene = bb[i][k];
                ++geneMatrix[thisGene][neiGene];
                ++geneMatrix[neiGene][thisGene];
            }
        }
    }
    
    for ( int i = 0; i < bbi.bbNum; ++i )
    {
        for ( int j = 1; j < bbi.bb[i][0] +1; ++j )
        {
            int thisGene = bbi.bb[i][j];
            int neiGene = 0;
            for ( int k = j +1; k < bbi.bb[i][0] +1; ++k )
            {
                neiGene = bbi.bb[i][k];
                ++geneMatrix[thisGene][neiGene];
                ++geneMatrix[neiGene][thisGene];
            }
        }
    }
    
    int total = 0;
    correct = 0;
    miss = 0;
    for ( int i = 0; i < ell; ++i )
    {
        for ( int j = 0; j < ell; ++j )
        {
            total = total + geneMatrix[i][j];
            if (geneMatrix[i][j] >= correctMatrix[i][j])
            {
                correctCount[i][j] = correctMatrix[i][j];
            }
            else
            {
                correctCount[i][j] = geneMatrix[i][j];
            }
            correct = correct + correctCount[i][j];

            if ( correctMatrix[i][j] > geneMatrix[i][j] )
            {
                miss = miss + correctMatrix[i][j] - geneMatrix[i][j];
            }
        }
    }
    
    
    
    
    for ( int i = 0; i < ell; ++i )
    {
        delete [] geneMatrix[i];
        delete [] correctMatrix[i];
        delete [] correctCount[i];
        delete [] errorCount[i];
    }
    delete [] geneMatrix;
    delete [] correctMatrix;
    delete [] correctCount;
    delete [] errorCount;

    return total;
    
//    vector<int>* gene2GG = new vector<int>[ell] ;
//    vector<int>* correctGene2GG = new vector<int>[ell] ;
//
//    for ( int i = 0; i < bbNum; ++i )
//    {
//        for ( int j = 1; j < bb[i][0] +1; ++j )
//        {
//            int thisGene = bb[i][j];
//            int neiGene = 0;
//            for ( int k = j +1; k < bb[i][0] +1; ++k )
//            {
//                neiGene = bb[i][k];
//                bool found = false;
//                for ( int p = 0 ; p < gene2GG[thisGene].size(); ++p )
//                {
//                    if ( gene2GG[thisGene][p] == neiGene )
//                    {
//                        found = true;
//                        break;
//                    }
//                }
//                if ( ! found )
//                {
//                    gene2GG[thisGene].push_back(neiGene);
//                    gene2GG[neiGene].push_back(thisGene);
//                }
//            }
//        }
//    }
//
//
//    for ( int i = 0; i < bbi.bbNum; ++i )
//    {
//        for ( int j = 1; j < bbi.bb[i][0] +1; ++j )
//        {
//            int thisGene = bbi.bb[i][j];
//            int neiGene = 0;
//            for ( int k = j +1; k < bbi.bb[i][0] +1; ++k )
//            {
//                neiGene = bbi.bb[i][k];
//                bool found = false;
//                for ( int p = 0 ; p < correctGene2GG[thisGene].size(); ++p )
//                {
//                    if ( correctGene2GG[thisGene][p] == neiGene )
//                    {
//                        found = true;
//                        break;
//                    }
//                }
//                if ( ! found )
//                {
//                    correctGene2GG[thisGene].push_back(neiGene);
//                    correctGene2GG[neiGene].push_back(thisGene);
//                }
//            }
//        }
//    }
//
//
//    // count error gene linkages
//
//    int error = 0;
//    correct = 0;
//    int * gene = new int[ell];
//    for ( int i = 0; i < ell; ++i )
//    {
//        for ( int j = 0; j < ell; ++j )
//        {
//            gene[j] = 0;
//        }
//        for ( int j = 0 ; j < gene2GG[i].size(); ++j )
//        {
//            ++gene[gene2GG[i][j]];
//        }
//        for ( int j = 0 ; j < correctGene2GG[i].size(); ++j )
//        {
//            --gene[correctGene2GG[i][j]];            
//        }
//
//        for ( int j = 0; j < ell; ++j )
//        {
//            error = error + abs(gene[j]);
//        }
//
//    }
//    delete [] gene;
//
//
//    delete [] gene2GG;
//    delete [] correctGene2GG;
//
//    return error;
}
