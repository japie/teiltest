#ifndef FITOC_H
#define FITOC_H

class FitOverConflict
{
    public:        
        FitOverConflict();
        int createProblem(int length, double over, double conf, int s); //s = bbsize
        double operator()(const int* const chromosome);
        double getMaxFitness();
        int bitanswer[505];
    private:
	void MergePoints(int, int);
	int JudgeIfPointCanMerge(int, int);
	int CalcGraphInfo(int, double, double, int);
	void MakeGraph();
	void JudgeBBs();
	double block_score(int, int);

        int flag;
        double real_score;
        int size, ell, totalBB;
	int DomBB, ConBB, LinkageNum;
	double tr_value;
        
        int graph[1005][9];
        int BBs[505][9];
	int tmpdata[505][2];

	// graph tmp
        int tmp[505][9];
	int BitForLink[505];

	// for mode1
	int degree[505];
	int sortde[505];

	// for disjoint set        
	int big[5005];
        int tran[5005];
        bool used[250005];

        int findb(int x);
};

#endif
