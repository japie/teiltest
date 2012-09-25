
#include <stdio.h>
#include "bitwisedistance.h"

BitwiseDistance::BitwiseDistance()
{
  // 2^16 = 65536
  distanceMap = new char[(1 << 16)];
  readMap();
}

BitwiseDistance::~BitwiseDistance()
{
  delete []distanceMap;
}

void BitwiseDistance::readMap()
{
    FILE *fp = fopen("distance.map", "rb");
    
    if (fp == NULL) {
	createArray();
	fp = fopen("distance.map", "wb");
	fwrite(distanceMap, sizeof(char), (1<<16), fp);
	fclose(fp);
    }
    else {
	fread(distanceMap, sizeof(char), (1<<16), fp);
	fclose(fp);
    }	
	
}

void BitwiseDistance::createArray()
{
    unsigned int i, j;
    for (i=0; i<(1<<16); i++) {
	char ones = 0;
	for (j=0; j<16; j++) {
	    if ( ((i >> j) & 1) == 1 )
		ones ++;
	}
	distanceMap[i] = ones;
    }
	
}

int BitwiseDistance::countOne(unsigned long distance)
{
  size_t i;
  int result = 0;
  for (i=0; i<sizeof(unsigned long)*8/16; i++) {
    int partialDistance = (int) (distance & 0xffff);
    result += distanceMap[partialDistance];
    distance >>= 16;
  }

  return result;
}

int BitwiseDistance::getHammingDistance(unsigned long a, unsigned long b)
{
  return countOne(a ^ b);
}



