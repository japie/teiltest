#include "myrand.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include "statistics.h"



#define  LEN_SUB 5
//#define  MAX_OVERLAP 3
using namespace std;




double checkRange(double d,int l)
{
   double token = d;
   if ( d*(-1) > 0 )
     token += l*(1+int(d)/l*(-1));
   else if ( int(d) >= l)
     token -= (int(d)/l)*l;
   else
     return d;//do nothing
   
   return token;
}


int main(int argc,char* argv[] )
{
    if (argc != 6)
    {
      printf("invalid usage, argvs = ell numGroup maxOverlap std numFiles\n");
      return 0;
    }
    //numGroup == ell / groupLength == j; default = 3 ( in Monetomo test )
    int ell      = atoi(argv[1]);
    int numGroup = atoi(argv[2]);
    int maxOverlap = atoi(argv[3]);
    double std   = atoi(argv[4]);
    int intstd   = int (std);
    int numFiles = atoi(argv[5]);
    MyRand rand;
    char buf[50];
    int* cache = new int[ell]();
    Statistics StGene,StStd,StMean,StMax;

    for(int fileCount = 0; fileCount < numFiles ; fileCount++)
    {    
        for ( int i = 0 ; i < ell ; i++ ) cache[i] = 0;

        sprintf(buf,"NormalSet/NSet_%d_%d_%d_%d_%d.dat",ell,numGroup,maxOverlap,intstd,fileCount);
        
        string ofileName(buf);
        ofstream ofile(ofileName.c_str());
        int ** r = new int*[numGroup];
		for ( int p= 0; p <numGroup; ++p )
			r[p] = new int[LEN_SUB];
       
        for ( int i = 0 ; i < numGroup ; i++)
        {
           double mean = 3*i;
           for ( int j = 0 ; j < LEN_SUB ; j++)
           {
               double temp  = rand.normal(mean,std);
    	   //check range
               temp = checkRange(temp,ell);
    	   //check occurance
    	       bool occur = true;
               int t = int(temp);
    	   
    	   while(occur)
    	   {
             int k = 0;
	     while(cache[t] >= maxOverlap)
	     {
                 t = int(checkRange(rand.normal(mean,std),ell));
		 //cout << t << "," << cache[t] << endl;
	     }
	     while(k<j)
    	     {
    	       while ( t == r[i][k] || cache[t] >= (maxOverlap) )
    	       {
	         //if(cache[t] > MAX_OVERLAP) cout << cache[t] << endl;
    	         temp = checkRange(rand.normal(mean,std),ell);
    		 t = int(temp);
		 
    		 occur = true;
    		 //printf("i,k,t,r = %d %d %d %d\n",i,k,t,r[i][k]);
    		 k = 0;
    		 //continue;
    	       }
	       k++;
	       
	       
	     }
	     
    	     r[i][j] = t;
             cache[t]++;
	     if(cache[t] > maxOverlap)
             cout << t<< ","<<cache[t]<< endl;
    	     occur = false;
	     
    	   }
    	   
    	   ofile << r[i][j] << " ";
           }
           ofile << endl;
        
        }
        for ( int i = 0 ; i < ell ; i++ ) StGene.record(cache[i]);
        StStd.record(StGene.getStdev());
        StMax.record(StGene.getMax());
        StMean.record(StGene.getMean());
	StGene.reset();
		
		for ( int p= 0; p <numGroup; ++p )
			delete [] r[p]; 
		delete [] r;

    }
    printf("overlap mean:  -Max- %f -Mean- %f  -Stdev- %f over %d times \n" , StMax.getMean(),StMean.getMean(),StStd.getMean(),numFiles);
    printf("maxOverlap in the sampled data: %f\n",StMax.getMax());
    //system("pause");
    return 0;
}

