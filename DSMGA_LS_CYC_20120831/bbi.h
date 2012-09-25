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
#ifndef BBI_H
#define BBI_H

/**
@author tianliyu
*/
class BBI
{
public:
  BBI ();
  BBI (int n_ell);
  BBI (const BBI & bbi);
  void init (int n_ell);
  void printOut ();
  BBI & operator = (const BBI & bbi);
  bool operator == (const BBI & bbi) const;
  bool operator != (const BBI & bbi) const;

  int getNumofMatchedBB(const BBI & bbi) const;
  int calculateBBIError(const BBI & bbi, int& correct, int& miss) const;

   ~BBI ();

  int bbNum;
  int prepareBB;
  int ell;
  int **bb;
};

#endif
