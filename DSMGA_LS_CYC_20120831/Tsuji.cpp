#ifndef TSUJI_CPP
#define TSUJI_CPP

#include "Tsuji.h"

Tsuji::Tsuji()
{
    reset();
}
Tsuji::~Tsuji()
{
    for(int i = 0 ; i < mym ; i++)  
    delete[] gen2BB[i]; 
    delete[] gen2BB; 
}

void
Tsuji::sampling(double w, int k, int ell, double std, int fCount){
    if(system("cd Tsuji"))
      system("mkdir Tsuji");
    //init cache for single gene  
    int* cache = new int[ell](); 
    for ( int i = 0 ; i < ell ; i++ ) cache[i] = 0;
                         
    Statistics StGene;
    //init Tsuji params
    setParams(w,k,ell,std); 
    double mu = double(ell / mym) ;   
    //init ouput file                 
    char buf[128];
    sprintf(buf,"Tsuji/Tsuji-w_%f-k_%d-ell_%d-std_%f-%d.dat",w,k,ell,std,fCount);
    string ofilename = string(buf);
    ofstream ofile;
    ofile.open(ofilename.c_str());
    gen2BB = new int*[mym];
    for ( int i = 0 ; i < mym ; i++ )
    {
   	   gen2BB[i] = new int[k];
       for ( int j = 0 ; j < k ; j++ )
	   {
	   	   double token = checkRange(myRand.normal(double(mu*i),std),ell);
  	       
	   	   bool occur = true;
	   	   while(occur)
		   {
		   	 
		   	 if( j >= 1 )
	   	     for ( int p = j-1 ; p >= 0 ; p-- )
	   	     {
			  	   
	   	           while( int(token) == gen2BB[i][p] )
	        	   {
				   	   token = checkRange(myRand.normal(double(mu*i),std),ell);	        	   
				   	   p = j-1;
   	   	               
	   			   }
	         }
	         occur = false;
		   }
	       gen2BB[i][j] = int(token);
		   cache[int(token)]++;	   
	   }	   
   	   
    }

    int BBNum = mym;
    vector<int> singleGene;
    for ( int i = 0 ; i < ell ; i++ )
    {
       StGene.record(cache[i]);
       if(cache[i] == 0)
       {
         BBNum++;
         singleGene.push_back(i);
        }
    }
    bbi.init(ell);
    bbi.bbNum = BBNum;
    for ( int i = 0 ; i < BBNum ; i++ )
    {
   	   int BBsize;
       if( i < mym )
	   {
	       BBsize = k;
	       bbi.bb[i][0] = BBsize;
	       ofile << bbi.bb[i][0] << ":" ;
	       for ( int j = 0 ; j < BBsize ; j++)
           {
	           bbi.bb[i][j+1] = gen2BB[i][j];
			   ofile << bbi.bb[i][j+1] << " " ;
	   	   }
		   ofile << "\n";
	   }	     
   	   else 
   	   {
		   
   	       BBsize = 1;
   	       bbi.bb[i][0] = BBsize;
   	       bbi.bb[i][1] = singleGene[i-mym];
   	       ofile << bbi.bb[i][0] << ":" ;
   	       ofile << bbi.bb[i][1] << endl;
	   }          
	   
    }
    if(ver) bbi.printOut();
    if(SHOW_GEN2BB) printf("gen2BB: (min/mean/max/std) = %f / %f / %f / %f \n",StGene.getMin(),StGene.getMean(),StGene.getMax(),StGene.getStdev());
   
   
    ofile.close();
    delete cache;
                       
}

void 
Tsuji::reset()
{
    tokens.clear();
    ver = false;          
    myw = 0;                        
    myk = 0;
    myell = 0;
    mystd = 0;
    mym = 0;          
              
}


double 
Tsuji::checkRange(double d,int l)
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

void
Tsuji::setParams(double iw, int ik, int iell, double istd)
{
    myw = iw;                        
    myk = ik;
    myell = iell;
    mystd = istd;
    mym = int(iw*iell/ik);
}

BBI 
Tsuji::getBBI(double w, int k, int ell, double std, int fCount)
{    
     setParams(w,k,ell,std); 
     if( loadSample( w,  k,  ell,  std,  fCount) == true)
       return bbi;
     else
       sampling( w,  k,  ell,  std,  fCount);
     return bbi; 			   
}

BBI
Tsuji::getBBI()
{
    return bbi;
}

int
Tsuji::getM()
{
    return mym;
}
bool
Tsuji::loadSample(double w, int k, int ell, double std, int fCount)
{
 	if(ver) cout << "loading file..." << endl;					 
 	char buf[128];					 
    sprintf(buf,"Tsuji/Tsuji-w_%f-k_%d-ell_%d-std_%f-%d.dat",w,k,ell,std,fCount);
    string ifilename = string(buf);					 
    ifstream infile;
	infile.open(ifilename.c_str(),ifstream::in);
	if( !infile.is_open()) return false;
	//succ loading
	setParams(w,k,ell,std); 
	int numBB = 0 ;
	int bbSize = 0;	
	
	bbi.init(ell);
	bbi.bbNum = ell;
	while(infile)
	{
	 	
	    string aline;
		getline(infile,aline);
		//if(ver) cout << aline << endl;
        mySplit(aline,": "); 			 
        if( (infile.eof()) )
        break;
		//cout << aline << endl;
		bbSize = tokens[0];
		if( bbSize == 0 )
		break;
		//cout << "bbsize = " << bbSize << endl;
		bbi.bb[numBB][0] = bbSize;
		
		for( int i = 0 ; i < bbSize ; i ++ )
		{
	 		//cout << tokens[i+1] << endl; 
	        bbi.bb[numBB][i+1] = (tokens[i+1]);
		}
		numBB++;		 		  
		tokens.clear();
    }
	bbi.bbNum = numBB;
	
    if (ver) bbi.printOut();
	return true;					 
}
/*
double 
Tsuji::getFitness(const Chromosome& c) const
{
     if( c.getLength() != ell )
     {
       cout << "Length does not match!" << endl;  
       return 0.0;                        
     }                  
     double fitness = 0.0;
     int unitary = 0;
     for ( int i = 0; i < mym; ++i )
     {
        unitary = 0;
        for ( int j = 0; j < myk; ++j )
        {
            if ( c.getVal(bbi.bb[i][j+1]) == 1 )
                ++unitary;
        }
        fitness = fitness + c.trap5(unitary);
     }
     return fitness;
}*/

void 
Tsuji::verbose()
{
    ver = true;
}

double
Tsuji::getMaxFitness()
{
    return mym - 0.11;
}

void
Tsuji::mySplit(const string &line,const string &delimiter)
{                                           
               
    char* token;
    char tempL[1024];
    line.copy(tempL,line.size(),0);
    token = strtok(tempL,delimiter.c_str());
    //cout << token << endl;
    //tokens.clear();
    while( token != NULL )
    {
  	 if ( token != "\n" )	
  	 {
          tokens.push_back(atoi(token));  
          token = strtok(NULL,delimiter.c_str());
	 }
    }
    //if(ver) {for(int i = 0 ; i < tokens.size(); i++) cout << tokens[i] << ",";}                                                       
}

#endif
