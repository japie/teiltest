#include <iostream>
#include <string.h>
#include <math.h>
#include "global.h"
#include "fitOC.h"
using namespace std;

int
FitOverConflict::createProblem(int length, double over, double conf, int s)
{
  //cout << "bb size  " << s << endl;
  //cout << "overlap  " << over << endl;
  //cout << "conflict " << conf << endl;
	tr_value = 0.9;
	//cout << "START" << endl;
	flag = CalcGraphInfo(length, over, conf, s);
	//cout << "CALC FINISH" << endl;
	if(flag == 0) 
	{
		MakeGraph();
		//cout << "MAKE GRAPH FINISH" << endl;
		JudgeBBs();
	}
	//cout << "END" << endl;
	//cout << "fitOC problem create done" << endl;
	/*
	for(int i = 0; i < ell; ++i)
		printf("%d", bitanswer[i]);
	printf(" %f\n==========\n", getMaxFitness());
	*/
	return flag;
}



double
FitOverConflict::getMaxFitness()
{
  return real_score;
}

int
FitOverConflict::findb(int x)
{
	if(big[x]==x)
	{
		return x;
	}
	else
	{
		return big[x]=findb(big[x]);
	}
}

int
FitOverConflict::JudgeIfPointCanMerge(int Blk1, int Blk2)
{
	if(Blk1>=totalBB && Blk2>=totalBB) return 0;
	if(Blk1==Blk2 || used[Blk1*totalBB+Blk2]==1) return 0;

	if(findb(tmp[Blk1][BitForLink[Blk1]])==findb(tmp[Blk2][BitForLink[Blk2]]))
	{
		BitForLink[Blk1]=(BitForLink[Blk1]+myRand.flip())%size;
		BitForLink[Blk2]=(BitForLink[Blk2]+myRand.flip())%size;
		return 0;
	}
	for(int m=0;m<size;m++)
	{
		if(findb(tmp[Blk1][BitForLink[Blk1]])==findb(tmp[Blk2][m])) return 0;
		if(findb(tmp[Blk2][BitForLink[Blk2]])==findb(tmp[Blk1][m])) return 0;
	}
	return 1;
}

void
FitOverConflict::MergePoints(int Blk1, int Blk2)
{//cout << "link: " << Blk1 << " " << Blk2 << endl;
	used[Blk1*(totalBB)+Blk2]=1;
	used[Blk2*(totalBB)+Blk1]=1;
	big[findb(tmp[Blk1][BitForLink[Blk1]])]=findb(tmp[Blk2][BitForLink[Blk2]]);
	BitForLink[Blk1]=(BitForLink[Blk1]+1)%size;
	BitForLink[Blk2]=(BitForLink[Blk2]+1)%size;

	degree[Blk1]++;
	degree[Blk2]++;

	return;
}

int 
FitOverConflict::CalcGraphInfo(int length, double over, double conf, int s)
{
	//cout << "bbsize: " << s << "\noverlap ratio: " << over << "\nconflict ratio: " << conf << endl;
	flag=0;
  over+=1;
	if(conf>=0.5+1E-8)
	{
		cout << "too many conf" << endl;
		return 1;
	}
	if(over<=1.0-1E-8)
	{
		cout << "too less overlap" << endl;
		return 1;
	}

    	//count the number of blocks, linkages..
	totalBB=(int)(1.0*length*over/s+0.5);
	LinkageNum=totalBB*s-length;
	ConBB=(int)(conf*totalBB+0.4999);
	DomBB=totalBB-ConBB;

	// save data
	real_score=DomBB+tr_value*(ConBB);
	totalBB=DomBB+ConBB;
	size=s;
	ell=length;
	
	if(LinkageNum<ConBB)
	{
		cout << "too less overlap" << endl;
		return 1;
	}
	else if(LinkageNum>totalBB*totalBB/4)
	{
		cout << "too much overlap" << endl;
		return 1;
	}
	
	return 0;
}

void
FitOverConflict::MakeGraph()
{
	int i,j,k,t;
	k = 0;
	//int BB1,BB2;
	
	// reset the data
	for(i=0;i<size*totalBB;i++)
	{
		big[i]=i;
		j=i/size;
		k=i%size;
		tmp[j][k]=i;
		BitForLink[j]=0;
	}

	memset(used,0,sizeof(used));
	memset(degree,0,sizeof(degree));
	//cout << totalBB << " " << LinkageNum << " " << ConBB << endl;

	// link SCC
	/*int round = LinkageNum / totalBB;
	for(i = 1 ; i <= round ; i++)
	{
		for(j = 0 ; j < totalBB ; j++)
		{
			BB1 = j;
			BB2 = (BB1 + i) % totalBB;
			for(;;)
			{
				if(JudgeIfPointCanMerge(BB1, BB2) == 1)
				{
					LinkageNum--;
					break;
				}
			}
			MergePoints(BB1, BB2);
		}
	}	*/
//cout << "linking" << endl;	
	// link the rest of linkages
	j=0;
	for(i=0;i<LinkageNum;i++)
	{
		for(j=0;j<totalBB;j++)
		{
			sortde[j]=j;
		}
		for(j=0;j<10*totalBB;j++)
		{
			int r1=rand()%totalBB;
			int r2=rand()%totalBB;
			k=sortde[r1];
			sortde[r1]=sortde[r2];
			sortde[r2]=k;
		}

		//sort by degree
		for(j=1;j<=totalBB;j++)
		{
			for(k=1;k<totalBB;k++)
			{
				if(degree[sortde[k]]<degree[sortde[k-1]])
				{
					t=sortde[k];
					sortde[k]=sortde[k-1];
					sortde[k-1]=t;
				}
			}
		}

		// choose one
		for(j=1;j<totalBB;j++)
		{
			for(k=0;k<j;k++)
			{
				if(JudgeIfPointCanMerge(sortde[j],sortde[k])==1) goto successfind;
			}
		}

		//cout << "bad!!!!!!!!!!!!!!!1" << endl;
		
		
		successfind:;
		MergePoints(sortde[j],sortde[k]);
		if(i<ConBB)
		{
			tmpdata[i][0]=sortde[j];
			tmpdata[i][1]=sortde[k];
		}
	}
	return;
}

void
FitOverConflict::JudgeBBs()
{
	int i, j;//, k;
	// judge the block which share the same cell
	j=0;
	for(i=0;i<size*totalBB;i++)
	{
		if(big[i]==i)
		{
			tran[i]=j;//cout << i << " " << j << endl;
			j++;
		}
		else
		{
			tran[i]=-1;
		}
	}
//cout << "judge ceil" << endl;
	for(i=0;i<totalBB;i++)
	{
		BBs[i][0]=1;
		for(j=1;j<=size;j++)
		{//cout << findb(tmp[i][j-1]) << " ";
			BBs[i][j]=tran[findb(tmp[i][j-1])];
		}//cout << endl;
	}
	if(ConBB!=0)
	{
		for(j=0;j<ConBB;j++)
		{//cout << tmpdata[j][1] << endl;
			BBs[tmpdata[j][1]][0] = 0;
		}
	}

	for(i=0;i<ell;i++)
	{
		bitanswer[i]=rand()%2;
	}

	// copy the data into graph
	for(i=0;i<totalBB;i++)
	{
		graph[i][0]=BBs[i][0];//cout << graph[i][0];
		for(j=1;j<=size;j++)
		{
			graph[i][j]=BBs[i][j];//cout << " " << graph[i][j];
		}//cout << endl;
	}

	return;
}



FitOverConflict::FitOverConflict()
{
  ;//do nothing
}

double 
FitOverConflict::operator()(const int* const chromosome)
{
	double sc=0;
	int t=0;

	for(int i=0;i<totalBB;i++)
	{
		t=0;
		for(int j=1;j<=size;j++)
 		{
			assert (graph[i][j] >= 0 && graph[i][j] < ell);
			if(bitanswer[graph[i][j]] == chromosome[graph[i][j]]) t++;
		}

		if(graph[i][0]==1)
		{//cout << i << " " << t << " " << block_score(t,1) << endl;
			sc+=block_score(t,1);
		}
		else
		{//cout << i << " " << t << " " << block_score(t,0) << endl;
			sc+=block_score(t,0);
		}
	}

    return sc;
}

double 
FitOverConflict::block_score(int t, int type)
{
	if(type==1)
	{
		if(t==5) return 1;
		else if(t==4) return 0.5;
		else if(t==3) return 0;
		else if(t==2) return tr_value/3;
		else if(t==1) return tr_value/3*2;
		else return tr_value;
	}
	else
	{
		if(t==5) return tr_value;
		else if(t==4) return tr_value/2;
		else if(t==3) return 0;
		else if(t==2) return 0.33;
		else if(t==1) return 0.66;
		else return 0.99;
	}
	return 0;
}
