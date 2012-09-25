#include <stdlib.h>
#include <stdio.h>
   
#include "random.h"
#include "nk.h"

//#define DEBUG

#define MULTIPLE_LOCAL_SEARCH_RESTARTS

/* #define COUNT_NODES */

double nk_opt_f;

#ifdef COUNT_NODES
long long num_cuts;
#endif

// ==============================================================================

void generate_nk(int n, int k, NK_instance *nk)
{
  int i,j, two2k1;

  int *generated=(int*) calloc(n,sizeof(int));

  if (k<=0)
    {
      printf("ERROR: Expecting k>0, got k=%d. Exiting.\n",k);
      exit(-1);
    }

  nk->n = n;
  nk->k = k;
  two2k1 = (1<<(k+1));

  nk->neighbors   = (int**) calloc(n,sizeof(int*));
  nk->friends     = (int**) calloc(n,sizeof(int*));
  nk->num_friends = (int*)  calloc(n,sizeof(int));
  nk->f           = (double**) calloc(n,sizeof(double*));
  nk->is_neighbor = (int**) calloc(n,sizeof(int*));

  nk->num_subproblems = (int*) calloc(n,sizeof(int));
  nk->subproblems = (int**) calloc(n,sizeof(int*));

  nk->num_subproblems2 = (int*) calloc(n,sizeof(int));
  nk->subproblems2 = (int**) calloc(n,sizeof(int*));

  // ---------------------------------------------

  for (i=0; i<n; i++)
    {
      nk->neighbors[i] = (int*) calloc(k+1,sizeof(int));
      nk->f[i] = (double*) calloc(two2k1,sizeof(double));
      nk->num_friends[i] = 0;
      nk->is_neighbor[i] = (int*) calloc(n,sizeof(int));

      for (j=0; j<n; j++)
	nk->is_neighbor[i][j]=-1;

      nk->num_subproblems[i]=0;
      nk->num_subproblems2[i]=0;
    }

  // ---------------------------------------------

  for (i=0; i<n; i++)
    {
      for (j=0; j<n; j++)
	generated[j]=0;
      generated[i]=1;

      nk->num_subproblems[i]++;

      for (j=0; j<k; j++)
	{
	  int neighbor;

	  do {
	    neighbor=intRand(n);
	  } while (generated[neighbor]);

	  generated[neighbor]=1;
	  nk->neighbors[i][j]=neighbor;
	  nk->num_friends[neighbor]++;
	  
	  nk->num_subproblems[neighbor]++;
	}

      nk->neighbors[i][k]=i;

      sort_int_array(nk->neighbors[i],k+1);

      for (j=0; j<k; j++)
	nk->num_subproblems2[nk->neighbors[i][j]]++;

      for (j=0; j<two2k1; j++)
	nk->f[i][j]=drand();
    }

  // ---------------------------------------------

  for (i=0; i<n; i++)
    {
      if (nk->num_friends[i]>0)
	{
	  nk->friends[i] = (int*) calloc(nk->num_friends[i],sizeof(int));
	  nk->num_friends[i]=0;
	}
      else
	nk->friends[i]=NULL;

      nk->subproblems[i] = (int*) calloc(nk->num_subproblems[i],sizeof(int));
      nk->num_subproblems[i] = 0;

      if (nk->num_subproblems2[i]>0)
	{
	  nk->subproblems2[i] = (int*) calloc(nk->num_subproblems2[i],sizeof(int));
	  nk->num_subproblems2[i] = 0;
	}
      else
	nk->subproblems2[i]=NULL;
    }

  // ---------------------------------------------

  for (i=0; i<n; i++)
    {
      for (j=0; j<=k; j++)
	{
	  int neighbor=nk->neighbors[i][j];
	  
	  if (neighbor!=i)
	    nk->friends[neighbor][nk->num_friends[neighbor]++]=i;
	  
	  nk->is_neighbor[i][neighbor]=j;
	  nk->subproblems[neighbor][nk->num_subproblems[neighbor]++]=i;
	  
	  if (j<k)
	    nk->subproblems2[neighbor][nk->num_subproblems2[neighbor]++]=i;
	}
    }
  
  free(generated);
}

// ==============================================================================

void free_nk(NK_instance *nk)
{
  int i;

  for (i=0; i<nk->n; i++)
    {
      free(nk->neighbors[i]);
      if (nk->friends[i])
	free(nk->friends[i]);
      free(nk->is_neighbor[i]);
      if (nk->subproblems[i])
	free(nk->subproblems[i]);
      if (nk->subproblems2[i])
	free(nk->subproblems2[i]);
      free(nk->f[i]);
    }

  free(nk->neighbors);
  free(nk->num_friends);
  free(nk->subproblems);
  free(nk->num_subproblems);
  free(nk->subproblems2);
  free(nk->num_subproblems2);
  free(nk->is_neighbor);
  free(nk->f);
  free(nk->friends);
}

// ==============================================================================

void prepare_for_solve_nk(NK_instance *nk, int **index, int *index_size, double *max_contrib, double ***best)
{
  int i,j;
  int two2k1 = (1<<(nk->k+1));

  for (i=0; i<nk->n; i++)
    index_size[i]=0;

  for (i=0; i<nk->n; i++)
    {
      int max_idx=i;
      for (j=0; j<=nk->k; j++)
	{
	  int neighbor=nk->neighbors[i][j];

	  if (neighbor>max_idx)
	    max_idx=neighbor;
	}

      index_size[max_idx]++;

      best[i]=(double**)calloc(nk->k,sizeof(double*));

      for (j=0; j<nk->k; j++)
	{
	  int ii,jj,num=1<<(j+1);
	  best[i][j]=(double*)calloc(num,sizeof(double));

	  for (ii=0; ii<num; ii++)
	    {
	      int base=(ii<<(nk->k-j));
	      int num2=((int)1)<<(nk->k-j);

	      double max_f=-1E30;

	      for (jj=0; jj<num2; jj++)
		if (nk->f[i][base+jj]>max_f)
		  max_f=nk->f[i][base+jj];
	      
	      best[i][j][ii]=max_f;
	    }
	}
    }

  for (i=0; i<nk->n; i++)
    {
      if (index_size[i]>0)
	index[i]=(int*) calloc(index_size[i],sizeof(int));
      else
	index[i]=NULL;

      index_size[i]=0;
    }

  for (i=0; i<nk->n; i++)
    {
      int max_idx=i;
      for (j=0; j<=nk->k; j++)
	{
	  int neighbor=nk->neighbors[i][j];
	  
	  if (neighbor>max_idx)
	    max_idx=neighbor;
	}
      
      index[max_idx][index_size[max_idx]++]=i;
    }

  for (i=0; i<nk->n; i++)
    {
      double max=-1E10;
      max_contrib[i]=0;

      for (j=0; j<two2k1; j++)
	if (nk->f[i][j]>max)
	  max=nk->f[i][j];

      max_contrib[i]=max;
    }
}

// ==============================================================================

double solve_nk(NK_instance *nk)
{
  int i,j;
  int *x = (int*) calloc(nk->n,sizeof(int));
  int **index;        // index[i] stores indices of subproblems that contain variable i and others <i
  int *index_size;   // index_size[i] stores the number of entries in index[i]
  double *max_contrib = (double*) calloc(nk->n,sizeof(double));
  double ***best;

  index=(int**) calloc(nk->n,sizeof(int*));
  index_size=(int*) calloc(nk->n,sizeof(int));
  best=(double***) calloc(nk->n,sizeof(double**));

  prepare_for_solve_nk(nk,index,index_size,max_contrib,best);

  double max_remain=0;
  
  for (i=0; i<nk->n; i++)
    max_remain+=max_contrib[i];

/*   nk_opt_f=-10E10; */

  nk_opt_f=local_search_nk(nk);
  
  printf("\nRunning BB\n");

#ifdef COUNT_NODES
  num_cuts=0;
#endif

  bb_nk(x,0,0,index,index_size,max_contrib,max_remain,nk,best);

#ifdef COUNT_NODES
  printf("Number of cuts: %llu (out of %llu)\n",num_cuts,((long long)1)<<(nk->n));
#endif
  
  free(x);
  free(max_contrib);

  for (i=0; i<nk->n; i++)
    {
      if (index[i]!=NULL)
	free(index[i]);
      for (j=0; j<nk->k; j++)
	free(best[i][j]);
      free(best[i]);
    }

  free(index);
  free(index_size);
  free(best);

  return nk_opt_f;
}

// ==============================================================================

double evaluate_nk(char *x, NK_instance *nk)
{
  int i,j;
  double f=0;

  for (i=0; i<nk->n; i++)
  {
    int idx=0;
  
    for (j=0; j<=nk->k; j++)
	  idx=((idx<<1)+x[nk->neighbors[i][j]]);
      
      f+=nk->f[i][idx];
  }

  return f;
}

// ==============================================================================

void bb_nk(int *x, int current, double current_f, int **index, int *index_size, double *max_contrib, double max_remain, NK_instance *nk, double ***best)
{
  int i,j,ii,n=nk->n,k=nk->k,num;
  double d;

  double inc0=0, inc1=0;
  
  for (i=0,num=index_size[current]; i<num; i++)
    {
      int idx=0;
      int which=index[current][i];
      int *neighbors=nk->neighbors[which];
      
      for (j=0; j<k; j++, neighbors++)
	idx=((idx<<1)+x[*neighbors]);
      
      idx<<=1;
      
      inc0+=nk->f[which][idx];
      inc1+=nk->f[which][idx+1];

      max_remain-=max_contrib[which];
    }
  
  if (current==n-1)
    {
      if (inc0>inc1)
	{
	  x[current]=0;
	  current_f+=inc0;
	}
      else
	{
	  x[current]=1;
	  current_f+=inc1;
	}

      if (current_f>nk_opt_f)
	{
	  nk_opt_f=current_f;
	  printf("Found better optimum: %lf\n",nk_opt_f);
	  printf("    (");
	  for (i=0; i<n; i++)
	    printf("%u",x[i]);
	  printf(")\n");
	}

      return;
    }
  
  // -------------------------------------------------------

  double max_remain0=0;
  for (i=0,num=nk->num_subproblems2[current]; i<num; i++)
    {
      int which=nk->subproblems2[current][i];

      int idx=0;
      int neighbor;

      for (ii=0; ((neighbor=nk->neighbors[which][ii])<current); ii++)
	idx=((idx<<1)+x[neighbor]);

/*       double max_f=best[which][ii][idx<<1]; */
      
/*       idx<<=(k-ii+1); */
      
/*       int max=(1<<(k-ii)); */
/*       double max_f=-10E10; */
      
/*       for (ii=0; ii<max; ii++) */
/* 	if ((d=nk->f[which][idx+ii])>max_f) */
/* 	  max_f=d; */
      
/*       double max_f=best[which][ii][idx<<1]; */

      max_remain-=max_contrib[which];
      max_remain0+=(max_contrib[which]=best[which][ii][idx<<1]);
/*       max_remain0+=max_f; */
    }

  max_remain+=max_remain0;
  
  if (current_f+inc0+max_remain>nk_opt_f)
    {
      x[current]=0;
      bb_nk(x,current+1,current_f+inc0,index,index_size,max_contrib,max_remain,nk,best);
    }
#ifdef COUNT_NODES
  else 
    num_cuts+=(1L<<(nk->n-current-2));
#endif

  max_remain-=max_remain0;

  // -------------------------------------------------------

  double max_remain1=0;
  for (i=0,num=nk->num_subproblems2[current]; i<num; i++)
    {
      int which=nk->subproblems2[current][i];
      
      int idx=0;
      int neighbor;
      
      for (ii=0; ((neighbor=nk->neighbors[which][ii])<current); ii++)
	idx=((idx<<1)+x[neighbor]);

/*       double max_f=best[which][ii][(idx<<1)+1]; */
      
/*       idx=(idx<<1)+1; */
/*       idx<<=(k-ii); */
      
/*       int max=(1<<(k-ii)); */
/*       double max_f=-10E10; */
      
/*       for (ii=0; ii<max; ii++) */
/* 	if ((d=nk->f[which][idx+ii])>max_f) */
/* 	  max_f=d; */
      
      max_remain1+=(max_contrib[which]=best[which][ii][(idx<<1)+1]);
/*       max_remain1+=max_f; */
    }

  max_remain+=max_remain1;

  if (current_f+inc1+max_remain>nk_opt_f)
    {
      x[current]=1;
      bb_nk(x,current+1,current_f+inc1,index,index_size,max_contrib,max_remain,nk,best);
    }
#ifdef COUNT_NODES
  else 
    num_cuts+=(1L<<(nk->n-current-2));
#endif

  // -------------------------------------------------------

  for (i=0,num=nk->num_subproblems[current]; i<num; i++)
    {
      int which=nk->subproblems[current][i];
      
      int idx=0;
      int neighbor;

      for (ii=0; ((neighbor=nk->neighbors[which][ii])<current); ii++)
	idx=((idx<<1)+x[neighbor]);
      
      idx<<=(k-ii+1);
      
      int max=(1<<(k-ii+1));
      double max_f=-10E10;
      
      for (ii=0; ii<max; ii++)
	if ((d=nk->f[which][idx+ii])>max_f)
	  max_f=d;
      
      max_contrib[which]=max_f;
    }
}

// ==============================================================================

int num_subproblems_nk(NK_instance *nk, int low, int high)
{
  int i,j;
  int num=0;

  for (i=low; i<high; i++)
    {
      int count_this=1;

      for (j=0; j<nk->k; j++)
	if ((nk->neighbors[i][j]<low)||(nk->neighbors[i][j]>=high))
	  count_this=0;

      num+=count_this;
    }

  return num;
}

// ==============================================================================

void sort_int_array(int *x, int n)
{
  int i;

/*   printf("before sort: "); */
/*   for (i=0; i<n; i++) */
/*     printf("%u ",x[i]); */
/*   printf("\n"); */

  for (i=1; i<n; i++)
    {
      int tmp=x[i];
      int j=i;

      while ((j>0)&&(x[j-1]>tmp))
	{
	  x[j]=x[j-1];
	  j--;
	}
      
      x[j]=tmp;
    }

/*   printf("after sort: "); */
/*   for (i=0; i<n; i++) */
/*     printf("%u ",x[i]); */
/*   printf("\n\n"); */
}

// ==============================================================================

double local_search_nk(NK_instance *nk)
{
  int i;
  int n=nk->n;
  double current_f;
  char    *x = (char*) malloc(n);
  double best_f=-10E30;
  int run;
  int num_runs=n*n;
  long num_flips=0;

  printf("Local search for i=%u (%u restarts)\n",n,num_runs);

#ifdef MULTIPLE_LOCAL_SEARCH_RESTARTS
  for (run=0; run<num_runs; run++)
#endif
    {
      for (i=0; i<n; i++)
	x[i]=(drand()<0.5)? 0:1;

      current_f=run_local_search_nk(x,nk,&num_flips);
            
      if (current_f>best_f)
	{
	  printf("  -> local optimum = %lf\n",current_f);
	  printf("    (");
	  for (i=0; i<n; i++)
	    printf("%u",x[i]);
	  printf(")\n");
	  best_f=current_f;
	}
    }
  
  free(x);

  return best_f;
}

// ==============================================================================

void save_nk(char *fname, NK_instance *nk)
{
  int i,j;
  FILE *f=fopen(fname,"w");

  fprintf(f,"%u %u\n",nk->n,nk->k);

  for (i=0; i<nk->n; i++)
    for (j=0; j<=nk->k; j++)
      fprintf(f,"%u\n",nk->neighbors[i][j]);

  for (i=0; i<nk->n; j++)
    for (j=0; j<(1<<(nk->k+1)); j++)
      fprintf(f,"%lf\n",nk->f[i][j]);

  fprintf(f,"Optimum=?");

  fclose(f);
}

// ==============================================================================

void load_nk(char *fname, NK_instance *nk)
{
  int i,j;
  int n,k;
  FILE *f=fopen(fname,"r");
  char s[250];

  fscanf(f,"%u %u\n",&nk->n,&nk->k);

  n=nk->n;
  k=nk->k;

  nk->neighbors   = (int**) calloc(n,sizeof(int*));
  nk->friends     = (int**) calloc(n,sizeof(int*));
  nk->num_friends = (int*)  calloc(n,sizeof(int));
  nk->f           = (double**) calloc(n,sizeof(double*));
  nk->is_neighbor = (int**) calloc(n,sizeof(int*));

  nk->num_subproblems = (int*) calloc(n,sizeof(int));
  nk->subproblems = (int**) calloc(n,sizeof(int*));

  nk->num_subproblems2 = (int*) calloc(n,sizeof(int));
  nk->subproblems2 = (int**) calloc(n,sizeof(int*));

  // ---------------------------------------------

  for (i=0; i<nk->n; i++)
    {
      nk->neighbors[i] = (int*) calloc(nk->k+1,sizeof(int));

      for (j=0; j<=nk->k; j++)
	{
	  int neighbor;
	  fscanf(f,"%u\n",&neighbor);
	  nk->neighbors[i][j]=neighbor;
	  nk->num_subproblems[neighbor]++;
	  if (j<nk->k)
	    nk->num_subproblems2[nk->neighbors[i][j]]++;
	  nk->num_friends[neighbor]++;
	}

      nk->is_neighbor[i] = (int*) calloc(nk->n,sizeof(int));
      for (j=0; j<nk->n; j++)
	nk->is_neighbor[i][j]=0;
    }

  // ---------------------------------------------

  for (i=0; i<nk->n; i++)
    {
      int max=(1<<(nk->k+1));

      nk->f[i] = (double*) calloc(max,sizeof(double));

      for (j=0; j<max; j++)
	fscanf(f,"%lf\n",&nk->f[i][j]);
    }

  // ---------------------------------------------

  for (i=0; i<nk->n; i++)
    {
      if (nk->num_friends[i]>0)
	{
	  nk->friends[i] = (int*) calloc(nk->num_friends[i],sizeof(int));
	  nk->num_friends[i]=0;
	}
      else
	nk->friends[i]=NULL;
      
      nk->subproblems[i] = (int*) calloc(nk->num_subproblems[i],sizeof(int));
      nk->num_subproblems[i] = 0;
      
      if (nk->num_subproblems2[i]>0)
	{
	  nk->subproblems2[i] = (int*) calloc(nk->num_subproblems2[i],sizeof(int));
	  nk->num_subproblems2[i] = 0;
	}
      else
	nk->subproblems2[i]=NULL;
    }

  // ---------------------------------------------  

  for (i=0; i<n; i++)
    {
      for (j=0; j<=k; j++)
	{
	  int neighbor=nk->neighbors[i][j];
	  
	  if (neighbor!=i)
	    nk->friends[neighbor][nk->num_friends[neighbor]++]=i;
	  
	  nk->is_neighbor[i][neighbor]=j;
	  nk->subproblems[neighbor][nk->num_subproblems[neighbor]++]=i;
	  
	  if (j<k)
	    nk->subproblems2[neighbor][nk->num_subproblems2[neighbor]++]=i;
	}
    }


  // ---------------------------------------------

  fgets(s,250,f);
/*   printf("sk_string=%s",s+8); */
  if (s[8]=='?')
    nk->optimum=10E+10;
  else
    sscanf(s,"Optimum: %lf\n",&nk->optimum);

  printf("nk_optimum = %lf\n",nk->optimum);

  // ---------------------------------------------

  fclose(f);
}
 
// ==============================================================================

double evaluate_flip(char *x, int i, NK_instance *nk)
{
  int ii,j;
  double f=0;

  for (ii=0; ii<nk->num_subproblems[i]; ii++)
    {
      int which=nk->subproblems[i][ii];
      int idx=0;
      int pos_i=-1000;

/*       printf("positions: "); */
      for (j=0; j<=nk->k; j++)
	{
	  int pos=nk->neighbors[which][j];
	 
/* 	  printf("%u(%u) ",pos,x[pos]); */
 
	  idx=((idx<<1)+x[pos]);

	  if (pos==i)
	    pos_i=j;
	}

/*       printf("\n"); */
      
      f-=nk->f[which][idx];

/*       printf("pos_i=%u\n",pos_i); */
/*       printf("idx before: %u\n",idx); */
      idx^=(int) (((int)1)<<(nk->k-pos_i));
/*       printf("idx after:  %u (%u)\n",idx,1<<(nk->k-pos_i)); */

/*       x[i]=1-x[i]; */
/*       idx=0; */
/*       for (j=0; j<=nk->k; j++) */
/* 	{ */
/* 	  int pos=nk->neighbors[which][j]; */
	  
/* 	  idx=((idx<<1)+x[pos]); */

/* 	  if (pos==i) */
/* 	    pos_i=j; */
/* 	} */
/*       x[i]=1-x[i]; */

/*       printf("idx correct:  %u\n",idx); */

      f+=nk->f[which][idx];
    }
  
  return f;
}

// ==============================================================================

double run_local_search_nk(char *x, NK_instance *nk, long *num_flips)
{
  int i;
  int n=nk->n;
  double current_f;
  int num_flips_local=0;

  double max_improvement;
  current_f=evaluate_nk(x,nk);
  
  // !!! can still optimize but not the bottleneck so what
  
  do {
    max_improvement=-1;
    int change=-1;
    for (i=0; i<n; i++)
      {
	double improvement;

	improvement=evaluate_flip(x,i,nk);
	
/* 	// just for debug!!! */
/* 	x[i]=1-x[i]; */
/* 	double new_f=evaluate_nk(x,nk); */
/* 	x[i]=1-x[i]; */
/* 	if (fabs(improvement-(new_f-current_f))>0.00000000001) */
/* 	  printf("Problem: Improvement: %lf (should be %lf)\n",improvement,new_f-current_f); */
/* 	else */
/* 	  printf("Improvement computed ok (%lf=%lf)\n",improvement,new_f-current_f); */
/* 	// !!! */

	if (improvement>max_improvement)
	  {
	    max_improvement=improvement;
	    change=i;
	  }
      }
    
    if (max_improvement>0)
      {
	x[change]=1-x[change];
	current_f+=max_improvement;
	num_flips_local++;
      }
    
    /* 	printf("  improvement = %e\n",max_improvement); */
  } while (max_improvement>10E-15);
  
  (*num_flips)+=num_flips_local;

  return current_f;
}
