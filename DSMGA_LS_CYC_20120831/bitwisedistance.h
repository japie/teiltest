
#ifndef _BITWISEDISTANCE_H_
#define _BITWISEDISTANCE_H_

class BitwiseDistance {

public:
  BitwiseDistance();
  ~BitwiseDistance();

  void readMap();
  void createArray();
  int countOne(unsigned long);
  int getHammingDistance(unsigned long, unsigned long);


 private:
  char *distanceMap;

};


#endif
