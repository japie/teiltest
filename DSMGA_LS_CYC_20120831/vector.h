#ifndef VECTOR_H
#define VECTOR_H

class Vector{
 public:
  int  n_elem;  // the size of elem[]
  int  n;       // the number of elements
  int  *elem;   // elements
  char *clem;   // flag for elements
  
  Vector(int n1){
    n      = 0;
    n_elem = n1;
    elem   = new int  [n_elem];
    clem   = new char [n_elem];
  }
  ~Vector(){
    delete [] elem;
    delete [] clem;
  };
  char is_in(int i){
    return clem[i];
  }

void ccopy(char *from, char *to, int n)
{
  int i;
  for(i = 0; i < n; i++)
    to[i] = from[i];
}

void icopy(int *from, int *to, int n)
{
  int i;
  for(i = 0; i < n; i++)
    to[i] = from[i];
}

  void add(int i){
    if(!is_in(i) && n < n_elem){
      elem[n] = i;
      clem[i] = 1;
      n++;
    }
  }
  void del(int i){
    if(is_in(i)){
      clem[i] = 0;
      int pos = where(i);
      elem[pos] = elem[n-1];
      n--;
    }
  }
  int where(int i){
    int j;
    for(j = 0; j < n; j++)
      if(elem[j] == i)	return j;
  }
  void init(){
    int i;
    n = 0;
    for(i = 0; i < n_elem; i++)
      clem[i] = 0;
  }
  void copy(Vector &v){
    n = v.n; 
    if(n_elem != v.n_elem){
      delete [] elem; 
      delete [] clem;
      elem = new int  [v.n_elem];
      clem = new char [v.n_elem];
    }
    n_elem = v.n_elem;
    icopy(v.elem, elem, n     );
    ccopy(v.clem, clem, n_elem);
  }
  void print(){
    int i;
    for(i = 0; i < n; i++)
      printf("%3d%c",elem[i],(i!=n-1)?' ':'\n');
  }
  void printc(){
    int i;
    for(i = 0; i < n; i++)
      printf("%d",clem[i]);
    if(n > 0)
      printf("\n");
  }
};

#endif
