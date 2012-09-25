/***************************************************************************
 *   Copyright (C) 2005 by Tian-Li Yu,,,                                   *
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
#ifndef TWODARRAY_H
#define TWODARRAY_H

/**
@author Tian-Li Yu
*/

#include <stdio.h>
#include <stdlib.h>

#define NDEBUG
#include <assert.h>

template < class T > class TwoDArray
{
public:
  TwoDArray ()
  {
    lengthX = 0;
    lengthY = 0;
    array = NULL;
  }

  TwoDArray (int n_lengthX, int n_lengthY)
  {
    lengthX = 0;
    lengthY = 0;
    array = NULL;
    init (n_lengthX, n_lengthY);
  }

  /** initialize TwoDArray with n_lengrhX by n_lengthY */
  void init (int n_lengthX, int n_lengthY)
  {
    release();
    
    int i,j;
    lengthX = n_lengthX;
    lengthY = n_lengthY;

    array = new T *[lengthX];
    for (i = 0; i < lengthX; i++)
      array[i] = new T[lengthY];

   
    for (i=0; i<lengthX; i++)
	for (j=0; j<lengthY; j++)
	    array[i][j] = T(0);
		   
  }

  ~TwoDArray ()
  {
      release();
  }

  void release ()
  {
    int i;

    if (array == NULL)
      return;
    else {
	for (i = 0; i < lengthX; i++)
            delete[](array[i]);

	delete[]array;
    }

    array = NULL;
    lengthX = lengthY = 0;

  }

  T & operator() (int i, int j) {
      
      assert (i >= 0 && i < lengthX && j >= 0 && j < lengthY);
      return array[i][j];
  }

  T *operator[] (int i) {
      
      assert(i>=0 && i< lengthX);
      return array[i];
  }

  int size() {
    
    int _size = 0;
    for(int i = 0 ; i < lengthX ; ++i)
      for(int j = 0 ; j < lengthY ; ++j)
        _size += sizeof(array[i][j]);

    _size += sizeof(lengthX);
    _size += sizeof(lengthX);

    return _size;
  }

private:
  int lengthX;
  int lengthY;
  T **array;

};

#endif
