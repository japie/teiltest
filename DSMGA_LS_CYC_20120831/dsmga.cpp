#include "dsmga.h"
#include "chromosome.h"
#include "UFO.h"
#include "fitOC.h"
#include "fastcounting.h"

#include <iostream>
#include <list>
#include <vector>
#include <bits/stl_vector.h>



using namespace std;

bool USE_LOCAL_SEARCH_AT_INIT = true;

DSMGA::DSMGA(int n_ell, int n_nInitial, int n_selectionPressure, double n_pc,
        double n_pm, int n_maxGen, int n_maxFe) {
    previousFitnessMean = -INF;
    fe = 0;
    ell = n_ell;
    nInitial = n_nInitial;
    nCurrent = nInitial;
    selectionPressure = n_selectionPressure;
    pc = n_pc;
    pm = n_pm;
    maxGen = n_maxGen;
    maxFe = n_maxFe;
    bbiCalculated = false;
		lsNFE = 0;
    population = new Chromosome[nInitial];
    offspring = new Chromosome[nInitial];

    for (int i = 0; i < nInitial; i++) {
        population[i].init(ell);
        offspring[i].init(ell);
    }
    selectionIndex = new int[nInitial];

    bbi.init(ell);

    initializePopulation();

}

DSMGA::~DSMGA() {
    delete[]population;
    delete[]offspring;
    delete[]selectionIndex;
}

void DSMGA::getBBI() {
    if (!bbiCalculated) {
        //getBBIOneMax ();
        //bbi = getBBIMune();
        
        //ufo.print();

        //bbi = ufo.getBBI();
        //bbi = tsuji.getBBI();
        DSMClusteringChromosome dsmcc;

        dsmcc.init(ell, 0.3333, 0.3333);


        dsmcc.createDSMCC(population, nCurrent, selectionPressure, selectionIndex);

        dsmcc.speedyHillClimbingV3();
        bbi = dsmcc.getBBI ();
        //calculateBBIError();
    }
    bbiCalculated = true;
}

int DSMGA::getDistance(Chromosome & c1, Chromosome & c2) {

    int i;

    int distance = 0;

    for (i = 0; i < c1.lengthLong; i++)
        distance += myBD.getHammingDistance(c1.gene[i], c2.gene[i]);

    return distance;
}

int DSMGA::getBBDistance(Chromosome & c1, Chromosome & c2) {
    getBBI();

    int i, j;

    int distance = 0;

    for (i = 0; i < bbi.bbNum; i++) {

        bool same = true;
        for (j = 1; j <= bbi.bb[i][0]; j++)
            if (c1.getVal(bbi.bb[i][j]) != c2.getVal(bbi.bb[i][j]))
                same = false;

        if (!same)
            distance++;
    }

    return distance;
}

void DSMGA::replacePopulation() {
    int i;

    if (SHOW_POPULATION) {
        printf("===== Current population =====\n");
        for (i = 0; i < nCurrent; i++) {
            printf("(%d)", i);
            if (SHORT_HAND)
                population[i].shortPrintOut();
            else
                population[i].printOut();
            printf("\n");
        }
    }

    if (SHOW_SELECTION_INDEX) {
        printf("===== Selection index =====\n");
        for (i = 0; i < nCurrent; i++)
            printf("%d ", selectionIndex[i]);
        printf("\n");
    }

    //fullReplace();
    RTR();
}

void DSMGA::fullReplace() {
    int i;

    for (i = 0; i < nCurrent; i++)
        population[i] = offspring[i];

    bbiCalculated = false;
}

/** RTR with window size = N */
void DSMGA::RTR() {

    int i, j;

    int windowSize = (nCurrent < ell) ? nCurrent : ell;
    int *randArray = new int[windowSize];

    for (i = 0; i < nCurrent; i++) {
        int index = -1;
        int distance;
        int minDistance = ell + 1; // max distance
        double minFitness = INF;

        myRand.uniformArray(randArray, windowSize, 0, nCurrent - 1);

        for (j = 0; j < windowSize; j++) {
            //distance = getBBDistance (offspring[i], population[randArray[j]]);
            distance = getDistance(offspring[i], population[randArray[j]]);
            if (distance < minDistance) {
                index = randArray[j];
                minDistance = distance;
                minFitness = population[index].getFitness();
            } else if (distance == minDistance
                    && minFitness > population[randArray[j]].getFitness()) {
                index = randArray[j];
                minFitness = population[index].getFitness();
            }
        }

        if (offspring[i].getFitness() > population[index].getFitness()) {

            population[index] = offspring[i];

            if (SHOW_REPLACEMENT) {
                printf("Replacing (%d) with ", index);
                if (SHORT_HAND)
                    offspring[i].shortPrintOut();
                else
                    offspring[i].printOut();
                printf("\n");
            }
        } else if (SHOW_REPLACEMENT) {
            printf("Not replacing (%d) with ", index);
            if (SHORT_HAND)
                offspring[i].shortPrintOut();
            else
                offspring[i].printOut();

            printf("\n");
        }
    }

    bbiCalculated = false;

    delete []randArray;

}

bool DSMGA::isSteadyState() {

    if (stFitness.getNumber() <= 0)
        return false;

    if (previousFitnessMean < stFitness.getMean()) {
        previousFitnessMean = stFitness.getMean() + 1e-6;
        return false;
    }

    return true;
}

bool DSMGA::findCorrectModel(int requiredMatchedM) {
    getBBIMKTrap();
    BBI perfectBBI(ell);
    perfectBBI = bbi;

    bbiCalculated = false;
    selection();
    getBBI();

    int numofMatchedBB = perfectBBI.getNumofMatchedBB(bbi);

    if (numofMatchedBB >= requiredMatchedM) return true;
    else return false;
}

int DSMGA::findInitialPopulation() {

    int i, j;

    getBBIMKTrap();
    BBI perfectBBI(ell);
    perfectBBI = bbi;

    bbiCalculated = false;
    selection();
    getBBI();

    int oldPopulationSize;
    int populationSize = nCurrent;
    int numofMatchedBB;

    do {

        oldPopulationSize = populationSize;
        populationSize *= 2;

        Chromosome * newPopulation = new Chromosome[populationSize];

        for (i = 0; i < oldPopulationSize; i++) {
            newPopulation[i].init(ell);
            newPopulation[i] = population[i];
        }

        for (i = oldPopulationSize; i < populationSize; i++) {

            for (j = 0; j < ell; j++)
                if (myRand.flip())
                    newPopulation[i].setVal(j, 0);
                else
                    newPopulation[i].setVal(j, 1);
        }

        delete []population;
        delete []selectionIndex;

        population = newPopulation;
        selectionIndex = new int[populationSize];

        nCurrent = populationSize;

        bbiCalculated = false;
        selection();
        getBBI();

        numofMatchedBB = perfectBBI.getNumofMatchedBB(bbi);

        if (populationSize > 10000) exit(-1);

    } while (numofMatchedBB < perfectBBI.bbNum - 1);

    int left = populationSize / 2;
    int right = populationSize;
    int middle;

    do {
        middle = (left + right) / 2;
        nCurrent = middle;

        bbiCalculated = false;
        selection();

        getBBI();

        /*
        bbi.printOut();
        printf("\n%d\n=========\n", middle);
         */
        numofMatchedBB = perfectBBI.getNumofMatchedBB(bbi);
        if (numofMatchedBB < perfectBBI.bbNum - 1)
            left = middle;
        else
            right = middle;
    } while (left + 2 < right);

    printf("\n");

    return middle;

}

int DSMGA::doIt(bool output, int mode) {
    generation = 0;
    while (!shouldTerminate()) {
        oneRun(output, mode);
        //bbi.printOut();
        //getchar();
    }
    //bbi.printOut();
    //getchar();

    return generation;
}

 
 void DSMGA::oneRun(bool output, int mode) {
    int i;
    
    selection();
    if (mode == 0)
        crossOver();
    else if (mode == 1)
        strengthBasedBayesianXO_TL();
    else if (mode == 2)
        gibbsXO();
    else if (mode == 3)
        onePointXO();
    else if (mode == 4)
        uniformXO();
    else if (mode == 5)
	    simpleMinCut();
    else
        bbUniformXO();


    mutation();

    replacePopulation();


    double max = -INF;
    stFitness.reset();
    for (i = 0; i < nCurrent; i++) {
        double fitness = population[i].getFitness();
        if (fitness > max) {
            max = fitness;
            bestIndex = i;
        }
        stFitness.record(fitness);
    }

    if (output)
        showStatistics();

    generation++;
}

bool DSMGA::shouldTerminate() {
    bool
    termination = false;

    if (maxFe != -1) {
        if (fe > maxFe)
            termination = true;
    }

    if (maxGen != -1) {
        if (generation > maxGen)
            termination = true;
    }
    //getMaxFitness will get the right answer 
    if (population[0].getMaxFitness() <= stFitness.getMax() + 1E-6)
        termination = true;

    return termination;

}

bool DSMGA::foundOptima() {
    return (stFitness.getMax() + 1E-6 >= population[0].getMaxFitness());
}

void DSMGA::showStatistics() {

    printf("Gen:%d  Fitness:(Max/Mean/Min):%f/%f/%f Chromsome Length:%d\n",
            generation, stFitness.getMax(), stFitness.getMean(),
            stFitness.getMin(), population[0].getLength());
    printf("best chromosome:");
    population[bestIndex].printOut();
    printf("\n");

    if (SHOW_LINKAGE)
        bbi.printOut();

    fflush(NULL);
}

void
DSMGA::getBBIOneMax() {
    int i;

    for (i = 0; i < ell; i++) {
        bbi.bb[i][0] = 1;
        bbi.bb[i][1] = i;
    }

    bbi.bbNum = ell;
}

void DSMGA::getBBIMKTrap() {
    int i, j;

    for (i = 0; i < ell / TRAP_K; i++) {
        bbi.bb[i][0] = TRAP_K;
        for (j = 1; j <= TRAP_K; j++)
            bbi.bb[i][j] = TRAP_K * i + (j - 1);
    }

    bbi.bbNum = ell / TRAP_K;
}

BBI DSMGA::getBBIMune() {

    BBI correctBBI(ell);

    //getchar();
    int i, j, idx;

    int* cache = new int[ell];

    for (int i = 0; i < ell; i++)
        cache[i] = 0;

    int numSingles = 0;
    for (i = 0; i < ell / DIS_M; i++) {
        correctBBI.bb[i][0] = GROUP_SIZE;
        for (j = 1; j <= GROUP_SIZE; j++) {
            idx = GROUP_SIZE * i + (j - 1);
            correctBBI.bb[i][j] = decisionV[idx];

            cache[decisionV[idx]]++;
            //printf("cache-%d ++ \n",decisionV[idx]);
        }
    }

    for (i = 0; i < ell; i++) {
        if (cache[i] == 0) {

            //printf("cache-%d is empty \n",i);
            correctBBI.bb[ell / DIS_M + numSingles][0] = 1;
            correctBBI.bb[ell / DIS_M + numSingles][1] = i;
            numSingles++;
        }

    }

    correctBBI.bbNum = ell / DIS_M + numSingles;
    //bbi.printOut();
    //getchar();
    delete[] cache;

    return correctBBI;
}


void DSMGA::selection() {

    tournamentSelection();
}


// tournamentSelection without replacement

void DSMGA::tournamentSelection() {
    int i, j;
    int randArrSize = selectionPressure * nCurrent;

    int randArray[randArrSize];

    for (i = 0; i < selectionPressure; i++)
        myRand.uniformArray(randArray + (i * nCurrent), nCurrent, 0,
            nCurrent - 1);

    for (i = 0; i < nCurrent; i++) {

        int winner = 0;
        double winnerFitness = -INF;

        for (j = 0; j < selectionPressure; j++) {
            int challenger = randArray[selectionPressure * i + j];
            double challengerFitness = population[challenger].getFitness();

            if (challengerFitness > winnerFitness) {
                winner = challenger;
                winnerFitness = challengerFitness;
            }

        }
        selectionIndex[i] = winner;
    }
}

void DSMGA::initializePopulation() {
    int i, j;

    for (i = 0; i < nInitial; i++)
        for (j = 0; j < ell; j++)
            if (myRand.flip())
                population[i].setVal(j, 1);
            else
                population[i].setVal(j, 0);
		if(USE_LOCAL_SEARCH_AT_INIT)
			localSearch();
		/*
		for (i = 0; i < nInitial; i++){
                    for (j = 0; j < ell; j++){
					printf("%d", population[i].getVal(j));
				}
				printf(" %f\n", population[i].getFitness());
		}
		printf("\n");
		*/
}


void DSMGA::localSearch() {
	int i, j;
	//printf("%s", "ls start\n");
	lsNFE = 0;
	for (i = 0; i < nInitial; i++){
		//
		double currentf, maxf;
		int maxIdx, lastFlip;
		lastFlip = maxIdx = -1;
		currentf = population[i].getFitness();
		do{
			maxf = currentf - 1;
			maxIdx = -1;
			for(j = 0; j < ell; j++){
				if(j != lastFlip){
					int bit = population[i].getVal(j);
					bit = 1 - bit;
					population[i].setVal(j, bit);
					lsNFE++;
					double newf = population[i].getFitness();
					if(newf > maxf){
						maxf = newf;
						maxIdx = j;
					}
					bit = 1 - bit;
					population[i].setVal(j, bit);
				}//end of if j != lastFlip
			}//end of iter of ell
			if(maxf > currentf){
				int bit = population[i].getVal(maxIdx);
				bit = 1 - bit;
				population[i].setVal(maxIdx, bit);
				currentf = maxf;
				lastFlip = maxIdx;
			}
			else break;
		}while(true);
	}//end of iter of population
	//printf("%d\n", "lsNFE\n");
}



int DSMGA::countOverlap(Chromosome & p1, Chromosome & p2) {
    int bbNum = bbi.bbNum;
    int j;

    int edgeCount = 0;
    for (j = 0; j < bbNum; j++) {
        bool flag = false;
        for (int k = 1; k <= bbi.bb[j][0]; k++) {
            if (p1.getVal(bbi.bb[j][k]) != p2.getVal(bbi.bb[j][k])) {
                flag = true;
                break;
            }
        }

        if (flag) {
            for (int k = 0; k < bbNum; k++) {
                if (j != k && isOverlapping(bbi.bb[j], bbi.bb[k], p1, p2)) {
                    edgeCount++;
                }
            }
        }
    }
    return edgeCount;
}

/*
void DSMGA::crossOver ()
{
    int i;
    getBBI ();
	
    if (nCurrent & 0x1 == 0)
    {                                             // nCurrent is even
        for (i = 0; i < nCurrent; i += 2)
        {
            crossOver (population[selectionIndex[i]],
                population[selectionIndex[i + 1]], offspring[i],
                offspring[i + 1]);
        }
    }
    else
    {
        for (i = 0; i < nCurrent - 1; i += 2)
        {
            crossOver (population[selectionIndex[i]],
                population[selectionIndex[i + 1]], offspring[i],
                offspring[i + 1]);
        }
        offspring[nCurrent - 1] = population[selectionIndex[nCurrent - 1]];
    }

}
 */

void DSMGA::mutation() {

    int i;
    getBBI();

    for (i = 0; i < nCurrent; i++) {
        bbMutation(offspring[i]);
    }
}

void DSMGA::bbMutation(Chromosome & ch) {
    int i, j;

    for (i = 0; i < bbi.bbNum; i++) {
        if (myRand.uniform() < pm) {
            for (j = 1; j <= bbi.bb[i][0]; j++) {
                ch.setVal(bbi.bb[i][j], myRand.flip());
            }
        }
    }
}

void DSMGA::bbUniformXO() {

    int i;

    if ((nCurrent & 0x1) == 0) {
    	// nNextGeneration is even

        for (i = 0; i < nCurrent; i += 2)
            bbUniformXO (population[selectionIndex[i]], population[selectionIndex[i + 1]],
                offspring[i], offspring[i + 1], 0.5);

    }
    else {
        for (i = 0; i < nCurrent - 1; i += 2) {
            bbUniformXO (population[selectionIndex[i]], population[selectionIndex[i + 1]],
                offspring[i], offspring[i + 1] , 0.5);
        }
        offspring[nCurrent - 1] =
            population[selectionIndex[nCurrent - 1]];
    }
}


void DSMGA::bbUniformXO(Chromosome & p1, Chromosome & p2, Chromosome & c1, Chromosome & c2, double prob) {
    int i, j;

    for (i = 0; i < bbi.bbNum; i++) {
        if (myRand.flip(prob)) {
            for (j = 1; j <= bbi.bb[i][0]; j++) {
                c1.setVal(bbi.bb[i][j], p1.getVal(bbi.bb[i][j]));
                c2.setVal(bbi.bb[i][j], p2.getVal(bbi.bb[i][j]));
            }
        } else {
            for (j = 1; j <= bbi.bb[i][0]; j++) {
                c1.setVal(bbi.bb[i][j], p2.getVal(bbi.bb[i][j]));
                c2.setVal(bbi.bb[i][j], p1.getVal(bbi.bb[i][j]));
            }
        }
    }

}


void DSMGA::calculateBBIError() {

    BBI correctBBI = getBBIMune();

    int match = bbi.getNumofMatchedBB(correctBBI);

    int correct = 0;
    int miss = 0;
    int total = bbi.calculateBBIError(correctBBI, correct , miss);


    FILE* outfile = fopen(bbiErrorOutFile,"a");
    fprintf (outfile, "population : %d   generation : %d   match : %d  bbNum : %d  totalL : %d  correct :  %d  miss : %d  \n",nCurrent, generation, match, bbi.bbNum, total, correct, miss);

    fclose(outfile);
    outfile=NULL;

//    vector<int>* gene2BB;           // gene by index, BB by strength big -> small
//    try {
//        gene2BB = new vector<int>[ell];                // gene by index, BB by strength big -> small
//    }
//    catch (std::bad_alloc)
//    {
//        printf("cannot new array gene2BB\n");
//        bbi.printOut();
//        exit(-1);
//    }
//    for ( int i = 0; i < bbi.bbNum; ++i )
//    {
//        for ( int j = 1; j < bbi.bb[i][0]+1 ; ++j )
//        {
//            gene2BB[bbi.bb[i][j]].push_back(i);
//        }
//    }
//
//
//    delete [] gene2BB;
}



void DSMGA::strengthBasedBayesianXO_TL() {


    getBBI ();
    

    //bbi.printOut();
    /*
	 //cout << "population after selection: \n";
	 for ( int i = 0; i < nCurrent; ++i )
	 {
	 population[selectionIndex[i]].printOut();
	 cout << endl;
	 }
	 */

    //system("pause");

    //  BBI should be obtained before this line
    if(bbi.bbNum==0)
    {
        cout << "No BB information, do normal crossover.\n";
        return crossOver();
    }

    const double threshold = 0.0;


    double*  bbStrength;                    // by index
    int*         sortedBB;                        // by strength
    bool*      bbFinished;                    // by index
    int*         undecidedGeneCount ;  // by index
    try {
        bbStrength = new double[bbi.bbNum];           // by index
        sortedBB = new int[bbi.bbNum];                      // by strength
        bbFinished = new bool[bbi.bbNum];               // by index
        undecidedGeneCount = new int[bbi.bbNum];  // by index
    }
    catch (std::bad_alloc)
    {
        printf("cannot new array bbStrength, sortedBB, bbFinished, or undecidedGeneCount\n");
        printf("bbi.bbNum = %d\n", bbi.bbNum);
        bbi.printOut();
        fflush(NULL);
        exit(-1);
    }

    vector<int>* gene2BB;           // gene by index, BB by strength big -> small
    try {
        gene2BB = new vector<int>[ell];                // gene by index, BB by strength big -> small
    }
    catch (std::bad_alloc)
    {
        printf("cannot new array gene2BB\n");
        bbi.printOut();
        exit(-1);
    }

    bool* geneDecided;              // by index
    try {
        geneDecided = new bool[ell];                               // by index

    }
    catch (std::bad_alloc)
    {
        printf("cannot new array geneDecided\n");
        bbi.printOut();
        exit(-1);
    }


    for ( int i = 0; i < ell; ++i )
    {
        geneDecided[i] = false;
    }


    for ( int i = 0; i < bbi.bbNum; ++i )
    {
        bbStrength[i] = getBBStrength(bbi.bb[i]);   // getBBStrength  needs to be modified better
        sortedBB[i] = i;
        undecidedGeneCount[i] = bbi.bb[i][0];
        bbFinished[i] = false;
    }

    // each gene2BB is not sorted
//    for ( int i = 0; i < bbi.bbNum; ++i )
//    {
//        int index = sortedBB[i];
//        for ( int j = 1; j < bbi.bb[index][0]+1 ; ++j )
//        {
//            gene2BB[bbi.bb[index][j]].push_back(index);
//        }
//    }


    QuickSort(sortedBB , bbStrength , 0 , bbi.bbNum-1);



//	     //for test
//	    cout << "BB num:   strength \n";
//	    for (int i = 0; i < bbi.bbNum; ++i )
//	    {
//	        cout << sortedBB[i] << " : " << bbStrength[sortedBB[i]] << endl;
//	    }
//            //system("pause");
//	     //end for test


    for ( int i = 0; i < bbi.bbNum; ++i )
    {
        int index = sortedBB[i];
        for ( int j = 1; j < bbi.bb[index][0]+1 ; ++j )
        {
            gene2BB[bbi.bb[index][j]].push_back(index);
        }
    }

//    	     //for test
//	    cout << "gene2BB :  \n";
//	    for (int i = 0; i < ell; ++i )
//	    {
//                for ( int j = 0; j < gene2BB[i].size(); ++j )
//                {
//                    cout << gene2BB[i][j] << "  ";
//                }
//                cout <<endl;
//	    }
//            system("pause");
//	     //end for test




    int decidedGene = 0;
    vector<int> candidateGene;
    candidateGene.reserve(5);
    vector<int> parentGene;
    parentGene.reserve(5);


    int thisBB = 0;
    int thisGene = 0;


    list<int> unfinishedBBList(bbi.bbNum);
    int temp = 0;
    for ( list<int>::iterator iter = unfinishedBBList.begin(); iter != unfinishedBBList.end(); ++iter)
    {
        (*iter) = sortedBB[temp];
        ++temp;
    }




    while ( decidedGene != ell )
    {

//        // for test
//        cout << "Unfinished BB : \n";
//        for ( list<int>::iterator iter = unfinishedBBList.begin(); iter != unfinishedBBList.end(); ++iter )
//        {
//            cout << *iter << "  ";
//        }
//        cout <<endl;
//        system("pause");
//        //end


        //  choose strongest BB
        if ( unfinishedBBList.empty() )
        {
            ;
//            cout << "Error happens, all BBs are finished, but not all genes\n";
//            assert(false);
        }
        else  // ! unfinishedBB.empty()
        {
            double maxStrength = -1.0;
            thisBB = unfinishedBBList.front();
            for (list<int>::iterator iter=unfinishedBBList.begin(); iter!=unfinishedBBList.end(); ++iter)
            {
                int bb = *iter;
                bool overlapping = false;

                for (int j=1; j<=bbi.bb[bb][0]; ++j)
                {
                    if (geneDecided[bbi.bb[bb][j]])
                    {
                        overlapping = true;
                        break;
                    }
                }

                if (overlapping)
                {
                    if ( bbStrength[*iter] > maxStrength)
                    {
                        thisBB = bb;
                        maxStrength = bbStrength[*iter];
                    }
                }
            }
        }

        // choose this Gene
        for ( int i = 1; i < bbi.bb[thisBB][0]+1 ; ++i )
        {
            if ( geneDecided[bbi.bb[thisBB][i]] == false )
                candidateGene.push_back(bbi.bb[thisBB][i]);
        }

        //cout << "candidateGene.size:  " << candidateGene.size()  <<endl;
        assert(candidateGene.size() > 0 );
        thisGene = candidateGene[myRand.uniformInt(0,candidateGene.size()-1)];

//        cout << "thisGene: " << thisGene << endl;
//        cout << "thisBB: " << thisBB << endl;

//  find thisGene's parentGenes
        bool* hasBeenParent = new bool[ell];
        for ( int r = 0; r < ell; ++r )
        {
            hasBeenParent[r] = false;
        }
        hasBeenParent[thisGene] = true;

        for ( unsigned int i = 0; i < gene2BB[thisGene].size(); ++i )
        {
            int checkBB = gene2BB[thisGene][i];
            if (bbStrength[checkBB] <= threshold )
            {
                continue;
            }


            for ( int j = 1; j < bbi.bb[checkBB][0]+1 ; ++j )
            {
                if ( geneDecided[bbi.bb[checkBB][j]] == true
                        && hasBeenParent[bbi.bb[checkBB][j]] == false )
                {
                    parentGene.push_back(bbi.bb[checkBB][j]);
                    hasBeenParent[bbi.bb[checkBB][j]] = true;
                }
            }
        }
        delete [] hasBeenParent;




        // dealing with parentGene.size() > MAX_PARENT
        // select those with smaller entropy
        if ( parentGene.size() > MAX_PARENT )
        {
            int * list = new int[parentGene.size()];
            int * tempParentGene = new int[parentGene.size()];
            double * negEntropy = new double[parentGene.size()];
						int parentGeneSize = parentGene.size();
            for ( int i = 0; i < parentGeneSize; ++i )
            {
                list[i] = i;
                tempParentGene[i] = parentGene[i];
                negEntropy[i] = 0;
                int tempGene = parentGene[i];
                int freq0 = 0;
                for (int j=0; j<nCurrent;  ++j)
                {
                    if (population[selectionIndex[j]].getVal(tempGene)==0)
                        freq0++;
                }
                if (freq0 != 0)
                {
                    double p0 = double(freq0) / double(nCurrent);
                    negEntropy[i] = p0 * log2(p0);
                }
                if (freq0 != nCurrent)
                {
                    double p1 = double(nCurrent-freq0) / double(nCurrent);
                    negEntropy[i] = negEntropy[i] + p1 * log2(p1);
                }

            }

            QuickSort( list, negEntropy, 0, parentGene.size()-1 );

            parentGene.clear();
            for ( int i = 0; i < MAX_PARENT; ++i )
            {              
                parentGene.push_back(tempParentGene[list[i]]);
            }
            delete [] list;
            delete [] negEntropy;
            delete [] tempParentGene;
        }




        // count every parent values combinition
        if (parentGene.empty() )
        {
            double p = 0.0;
            for ( int i = 0; i < nCurrent; ++i)
            {
                    if (population[selectionIndex[i]].getVal(thisGene) == 1)
                            p += 1.0;
            }

            p /= double(nCurrent);

            //cout << "No parent for gene :  "<< thisGene << "  , use marginal P = " << p <<endl;

            for ( int i=0; i < nCurrent; ++i )
            {
                    if (myRand.flip(p))
                            offspring[i].setVal(thisGene, 1);
                    else
                            offspring[i].setVal(thisGene, 0);
            }
        }
        else
        {
            long int numOfPattern = (int)pow(double(2),double(parentGene.size()));
            const long int maxRecordSize = (int)pow(double(2), double(25));
            long int patternOffset = 0;            
            
            long int patternLeft = numOfPattern;
            long int recordSize = 0;
            
            while ( (patternLeft = numOfPattern - patternOffset) > 0 )
            {
                long int fromPattern = patternOffset;
                
                if ( patternLeft > maxRecordSize )
                {
                    recordSize = maxRecordSize;
                    patternLeft = patternLeft - maxRecordSize;
                    patternOffset = patternOffset + maxRecordSize;                    
                }
                else
                {
                    recordSize = patternLeft;
                    patternLeft = 0;
                    patternOffset = patternOffset + recordSize;
                }                
                long int toPattern =  patternOffset -1;

                double p = 0.0;
                int *count;
                int *total;
                try {
                    count = new int[recordSize];
                    total = new int[recordSize];
                }
                catch (std::bad_alloc)
                {
                    printf("Could not allocate memory at int *count = new int[recordSize]; int *total = new int[recordSize];\n");
                    printf("numOfPattern = %ld\n", numOfPattern);
                    fflush(NULL);
                    exit(-1);
                }
                for ( int i = 0; i < recordSize; ++i )
                {
                    count[i] = 0;
                    total[i] = 0;
                }



                for ( int i=0; i < nCurrent; ++i )
                {
                    long int pattern = population[selectionIndex[i]].makeInt(parentGene);
                    if ( pattern < fromPattern || pattern > toPattern )
                        continue;
                    
    //                cout << "numOfPattern = " << numOfPattern <<"  pattern = " << pattern <<endl;
    //                cout << "population[selectionIndex[["<<i<<"]] : \n";
    //                population[selectionIndex[i]].printOut();
    //                cout <<endl;
                    ++total[pattern - fromPattern];
                    if (population[selectionIndex[i]].getVal(thisGene) == 1)
                    {
                        p = p + 1.0;
                        ++count[pattern - fromPattern];
                    }
                }

                p = p / double(nCurrent);

                // check every offspring if he satisfies parentGene's pattern
                for (int i=0; i<nCurrent; i++)
                {
                    long int pattern = offspring[i].makeInt(parentGene);
                    if ( pattern < fromPattern || pattern > toPattern )
                        continue;
                    if (total[pattern - fromPattern] == 0)  // new pattern, use solo p
                    {
                        //cout <<"Offspring "<< i << " No parent pattern for gene " << thisGene <<" in population"<< "  , use marginal P = " << p <<endl;
                        if (myRand.flip(p))
                            offspring[i].setVal(thisGene, 1);
                        else
                            offspring[i].setVal(thisGene, 0);
                    }
                    else
                    {
                        long int index =  pattern - fromPattern;
                        //cout << "Offspring "<< i <<"  Gene:  "<< thisGene << "  , use conditional P = " << double(count[pattern]) / double(total[pattern]) <<endl;
                        if ( double(count[index]) / double(total[index]) < -1e-6 ||  double(count[index]) / double(total[index]) > 1)
                        {
                            cout << "Offspring "<< i <<"  Gene:  "<< thisGene << "  , use conditional P = " << double(count[index]) / double(total[index]) <<endl;
                            cout << "pattern = " << pattern <<endl;
                            cout << "count[pattern]  = "<< count[index] <<endl;
                            cout << "total[pattern] = "<< total[index]<<endl;
                        }

                        assert(double(count[index]) / double(total[index]) >= 0 && double(count[index]) / double(total[index]) <= 1);
                        if (myRand.flip(double(count[index]) / double(total[index])))
                            offspring[i].setVal(thisGene, 1);
                        else
                            offspring[i].setVal(thisGene, 0);
                    }
                }

                delete [] total;
                delete [] count;

            }


//            long int numOfPattern = (int)pow(2,parentGene.size());
//            printf("numOfPattern = %d   parentGene.size() = %d\n", numOfPattern, parentGene.size());
//            fflush(NULL);
//	    double p = 0.0;
//
//            int *count;
//            int *total;
//            try {
//                count = new int[numOfPattern];
//                total = new int[numOfPattern];
//            }
//            catch (std::bad_alloc)
//            {
//                printf("Could not allocate memory at int *count = new int[numOfPattern]; int *total = new int[numOfPattern];\n");
//                printf("numOfPattern = %d\n", numOfPattern);
//                fflush(NULL);
//                exit(-1);
//            }
//            for ( int i = 0; i < numOfPattern; ++i )
//            {
//                count[i] = 0;
//                total[i] = 0;
//            }
//
//
//
//            for ( int i=0; i < nCurrent; ++i )
//            {
//                long int pattern = population[selectionIndex[i]].makeInt(parentGene);
//                assert(pattern < numOfPattern && pattern >= 0);
////                cout << "numOfPattern = " << numOfPattern <<"  pattern = " << pattern <<endl;
////                cout << "population[selectionIndex[["<<i<<"]] : \n";
////                population[selectionIndex[i]].printOut();
////                cout <<endl;
//                ++total[pattern];
//                if (population[selectionIndex[i]].getVal(thisGene) == 1)
//                {
//                    p = p + 1.0;
//                    ++count[pattern];
//                }
//            }
//
//            p = p / double(nCurrent);
//
//            for (int i=0; i<nCurrent; i++)
//            {
//                long int pattern = offspring[i].makeInt(parentGene);
//		assert(pattern < numOfPattern && pattern >= 0);
//                if (total[pattern] == 0)  // new pattern, use solo p
//                {
//                    //cout <<"Offspring "<< i << " No parent pattern for gene " << thisGene <<" in population"<< "  , use marginal P = " << p <<endl;
//                    if (myRand.flip(p))
//                        offspring[i].setVal(thisGene, 1);
//                    else
//                        offspring[i].setVal(thisGene, 0);
//                }
//                else
//                {
//                    //cout << "Offspring "<< i <<"  Gene:  "<< thisGene << "  , use conditional P = " << double(count[pattern]) / double(total[pattern]) <<endl;
//                    if ( double(count[pattern]) / double(total[pattern]) < -1e-6 ||  double(count[pattern]) / double(total[pattern]) > 1)
//                    {
//                        cout << "Offspring "<< i <<"  Gene:  "<< thisGene << "  , use conditional P = " << double(count[pattern]) / double(total[pattern]) <<endl;
//                        cout << "pattern = " << pattern <<endl;
//                        cout << "count[pattern]  = "<< count[pattern] <<endl;
//                        cout << "total[pattern] = "<< total[pattern]<<endl;
//                    }
//
//                    assert(double(count[pattern]) / double(total[pattern]) >= 0 && double(count[pattern]) / double(total[pattern]) <= 1);
//                    if (myRand.flip(double(count[pattern]) / double(total[pattern])))
//                        offspring[i].setVal(thisGene, 1);
//                    else
//                        offspring[i].setVal(thisGene, 0);
//                }
//            }
//
//            delete [] total;
//            delete [] count;
        }

        ++decidedGene;
        geneDecided[thisGene] = true;


        for ( unsigned int i = 0; i < gene2BB[thisGene].size(); ++i )
        {
            int checkBB = gene2BB[thisGene][i];
            assert(undecidedGeneCount[checkBB] > 0);
            if ( (--undecidedGeneCount[checkBB] ) == 0 )
            {
                bbFinished[checkBB] = true;
            }
        }

        list<int>::iterator iter = unfinishedBBList.begin();
        list<int>::iterator iter_end = unfinishedBBList.end();
        while (iter != iter_end)
        {
            if (undecidedGeneCount[*iter] == 0)
            {
                iter = unfinishedBBList.erase(iter);
                iter_end = unfinishedBBList.end();
            }
            else
            {
                ++iter;
            }
        }

        thisBB = 0;
        thisGene = 0;

        candidateGene.clear();
        candidateGene.reserve(5);
        parentGene.clear();
        parentGene.reserve(5);
    }

    delete []        bbStrength;
    delete []        sortedBB;
    delete []        bbFinished;
    delete []        undecidedGeneCount;
    delete []        gene2BB;
    delete []        geneDecided;
    //cout << "Finish all offspring\n";
}

double DSMGA::getBBStrength(int *bb)
{
    // we cannot deal with too large BB.
    // should use a hash table when possible
    if ( bb[0] >= 20 )
    {
        printf("the size of bb is %d\n", bb[0]);
    }
    assert (bb[0] < 20);

    double entropyAll = 0.0;
    double entropyIndividual = 0.0;

    int i,j;
    int pow2bb = pow2(bb[0]);
    int *freq = (int *) calloc (sizeof(int), pow2bb);

    for (i=0; i<nCurrent; i++) {
	int pattern = population[selectionIndex[i]].makeInt(bb);
	freq[pattern]++;
    }

    for (i=0; i<pow2bb; i++)
	if (freq[i] > 0 && freq[i] < nCurrent) {
	    double p = double(freq[i]) / double(nCurrent);
	    entropyAll -= p * log2(p);
	}

    free(freq);


    for (i=1; i<=bb[0]; i++) {
        int freq0 = 0;
	for (j=0; j<nCurrent; j++) {
	    if (population[selectionIndex[j]].getVal(bb[i])==0)
		freq0++;
	}
	if (freq0 != 0) {
	    double p0 = double(freq0) / double(nCurrent);
	    entropyIndividual -= p0 * log2(p0);
	}
	if (freq0 != nCurrent) {
	    double p1 = double(nCurrent-freq0) / double(nCurrent);
	    entropyIndividual -= p1 * log2(p1);
	}

    }


    if (entropyIndividual > entropyAll)
    	return (entropyIndividual - entropyAll);
    else
	return 0.0;
}

// Edited Version. Implementing Monetomo's crossover. by Hsuan Lee

bool DSMGA::isOverlapping(int *bb1, int *bb2, Chromosome &p1, Chromosome &p2) {
    int i;

    char *genes = (char*) calloc(ell, sizeof (char));

    for (i = 1; i <= bb1[0]; i++)
        genes[bb1[i]] = 1;

    for (i = 1; i <= bb2[0]; i++)
        if (genes[bb2[i]] == 1) {
            if (p1.getVal(bb2[i]) != p2.getVal(bb2[i])) {
                free(genes);
                return true;
            }
        }

    free(genes);
    return false;
}

bool isOverlap(int *bb1, int *bb2) {

    for (int i = 1; i <= bb1[0]; i++)
        for (int j = 1; j <= bb2[0]; j++)
            if (bb1[i] == bb2[j])
                return true;
    return false;
}

void DSMGA::construct_the_original_graph(Network &onet) {
    int i, j, n, m;

    n = bbi.bbNum;

    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            // add an edge if l[i] and l[j] overlap
            if (isOverlap(bbi.bb[i], bbi.bb[j])) {
                onet.add_edge(i, j);
                onet.add_edge(j, i);
            }
        }
    }

    onet.max_incomes = 0;
    for (i = 0; i < n; i++) {
        m = 0;
        for (j = 0; j < n; j++) {
            m += onet.E[i][j];
        }
        onet.max_incomes = (onet.max_incomes < m) ? m : onet.max_incomes;
    }

    onet.edge_list = new int * [n];
    for (i = 0; i < n; i++) {
        onet.edge_list[i] = new int [onet.max_incomes];
        m = 0;
        for (j = 0; j < n; j++) {
            if (onet.E[i][j]) onet.edge_list[i][m++] = j;
        }
        if (m < onet.max_incomes) onet.edge_list[i][m] = -1;

        /* Check that edges are added in edge_list correctly
        printf("%d : ",i);
        for(j = 0; j < onet.max_incomes && onet.edge_list[i][j]!=-1; j++){
          printf("%d ",onet.edge_list[i][j]);
        }
        printf("\n");
         */
    }

}

char do_disrupt_bb(Chromosome &p1, Chromosome &p2, int *bb1, int *bb2) {
    for (int i = 1; i <= bb1[0]; i++)
        for (int j = 1; j <= bb2[0]; j++)
            if (bb1[i] == bb2[j])
                if (p1.getVal(bb1[i]) != p2.getVal(bb1[i]))
                    return 1;

    return 0;
}

bool DSMGA::sameSubstrings(Chromosome &p1, Chromosome &p2, int *bb) {
    for (int k = 1; k <= bb[0]; k++) {
        if (p1.getVal(bb[k]) != p2.getVal(bb[k])) {
            return false;
        }
    }
    return true;
}

void DSMGA::reconstruct_graph(Chromosome &p1, Chromosome &p2, Network &oNet, Network &nNet) {
    int i, j;
    int node1, node2;
    // remove nodes whose BBs are identical
    nNet.alloc_map(oNet.n);
    nNet.n = 0;

    for (i = 0; i < oNet.n; i++) {
        
        if (!sameSubstrings(p1, p2, bbi.bb[i])) {
            nNet.map[nNet.n] = i;
            oNet.map[i] = nNet.n;
            //cout << nNet.n << ": Node " << nNet.map[nNet.n] << endl ;
            nNet.n++;
        } else {
            oNet.map[i] = -1;
        }
    }

    // add appropriate edges where edges which don't cause BB disruptions are removed
    for (i = 0; i < oNet.n; i++) {
        for (j = 0; j < oNet.max_incomes && oNet.edge_list[i][j] != -1; j++) {
            node1 = i;
            node2 = oNet.edge_list[i][j];
            if (oNet.map[node1] == -1 || oNet.map[node2] == -1) continue;
            if (nNet.E[oNet.map[node1]][oNet.map[node2]]) continue;
            if (do_disrupt_bb(p1, p2, bbi.bb[node1], bbi.bb[node2]) == 1) {
                //add edges where overlapping relationship doesn't cause BB disrutipn
                nNet.add_edge(oNet.map[node1], oNet.map[node2]);
                nNet.add_edge(oNet.map[node2], oNet.map[node1]);
            }
        }
    }
    //cout << "test4" << endl;
}

void DSMGA::reconstruct_graph_s(Chromosome &p1, Chromosome &p2, Network &oNet, Network &nNet) {
    int i, j;
    int node1, node2;
    // remove nodes whose BBs are identical
    nNet.alloc_map(oNet.n);
    nNet.n = 0;

    for (i = 0; i < oNet.n; i++) {
        
        //if (!sameSubstrings(p1, p2, bbi.bb[i])) {
            nNet.map[nNet.n] = i;
            oNet.map[i] = nNet.n;
            //cout << nNet.n << ": Node " << nNet.map[nNet.n] << endl ;
            nNet.n++;
        //} else {
        //    oNet.map[i] = -1;
        //}
    }

    // add appropriate edges where edges which don't cause BB disruptions are removed
    for (i = 0; i < oNet.n; i++) {
        for (j = 0; j < oNet.max_incomes && oNet.edge_list[i][j] != -1; j++) {
            node1 = i;
            node2 = oNet.edge_list[i][j];
            if (oNet.map[node1] == -1 || oNet.map[node2] == -1) continue;
            if (nNet.E[oNet.map[node1]][oNet.map[node2]]) continue;
            //if (do_disrupt_bb(p1, p2, bbi.bb[node1], bbi.bb[node2]) == 1) {
                //add edges where overlapping relationship doesn't cause BB disrutipn
                nNet.add_edge(oNet.map[node1], oNet.map[node2]);
                nNet.add_edge(oNet.map[node2], oNet.map[node1]);
            //}
        }
    }
    //cout << "test4" << endl;
}

void get2rands(int *v1, int *v2, int n) {
    *v1 = myRand.uniformInt(0, 10000)%n;
    do {
        *v2 = myRand.uniformInt(0, 10000)%n;
    } while (*v1 == *v2);
}

class Linkage {
public:
    int n_elem; // memory size of *elem
    int n; // number of elements
    int *elem;

    Linkage() {
        n = 0;
        elem = 0;
        n_elem = 0;

        elem = NULL;
    }

    ~Linkage() {
        free();
    }
    void add(int i);
    void allocate(int n_elem1);
    char is_in(int i);
    void free();

    void init() {
        n = 0;
    }

    void addset(Linkage &from) {
        int i;
        for (i = 0; i < from.n; i++) {
            add(from.elem[i]);
        }
    }

    void addset(int *bb) {
        for (int i = 1; i <= bb[0]; i++)
            add(bb[i]);
    }

    void exchange(Chromosome &p1, Chromosome &p2, Chromosome &c1, Chromosome &c2, int len) {
        int i;

        for (int i = 0; i < len; i++) {
            c1.setVal(i, p1.getVal(i));
            c2.setVal(i, p2.getVal(i));
        }

        for (i = 0; i < n; i++) {
            c1.setVal(elem[i], p2.getVal(elem[i]));
            c2.setVal(elem[i], p1.getVal(elem[i]));
        }
    }
};

void Linkage::allocate(int n_elem1) {
    if (elem) delete [] elem;
    n_elem = n_elem1;
    elem = new int [n_elem];
}

char Linkage::is_in(int i) {
    int j;
    for (j = 0; j < n; j++) {
        if (elem[j] == i)
            return 1;
    }
    return 0;
}

void Linkage::add(int i) {
    int j;

    if (!is_in(i)) {
        if (n >= n_elem) {
            n_elem = n + 1;
            int *tmpElem = new int [n_elem];
            for (j = 0; j < n; j++) {
                tmpElem[j] = elem[j];
            }
            tmpElem[n] = i;
            n++;
            if (elem)
                delete [] elem;
            elem = tmpElem;
        } else {
            elem[n] = i;
            n++;
        }
    }
}

void Linkage::free() {
    if (elem) delete [] elem;

    elem = NULL;
    n_elem = 0;
    n = 0;
}

void DSMGA::cdc_for_1pair_s(Chromosome &p1, Chromosome &p2, Chromosome &c1, Chromosome &c2, Network &onet) {

    int vs, vt, i;
    int ncut;
    Network net(bbi.bbNum);

    // reconstuct the graph onet according to the pair of strings (p1, p2)
    reconstruct_graph_s(p1, p2, onet, net);

    Linkage exchange_set;
    exchange_set.allocate(ell);

    if (net.n >= 3) {
        // find a mincut, make a exchange set and exchange strings
        Vector v(net.n);

        // -- get 2 random nodes
        get2rands(&vs, &vt, net.n);
        //cout << "vs = " << vs << ", vt = " << vt << endl;

        // -- find a mincut
        ncut = net.mincut(vs, vt, v);
        
        // -- make a exchange set
        for (i = 0; i < v.n; i++) {
            exchange_set.addset(bbi.bb[ net.map[v.elem[i]] ]);
        }
        // -- exchange strings
        exchange_set.exchange(p1, p2, c1, c2, ell);
    } else if (net.n == 2) {
        int a = myRand.uniformInt(0, 100)%2;
        exchange_set.addset(bbi.bb[net.map[a]]);
        exchange_set.exchange(p1, p2, c1, c2, ell);
    } else {
        // exchange_set = {}
        exchange_set.exchange(p1, p2, c1, c2, ell);
    }
}

void DSMGA::cdc_for_1pair(Chromosome &p1, Chromosome &p2, Chromosome &c1, Chromosome &c2, Network &onet) {

    int vs, vt, i;
    int ncut;
    Network net(bbi.bbNum);

    // reconstuct the graph onet according to the pair of strings (p1, p2)
    reconstruct_graph(p1, p2, onet, net);

    Linkage exchange_set;
    exchange_set.allocate(ell);

    if (net.n >= 3) {
        // find a mincut, make a exchange set and exchange strings
        Vector v(net.n);

        // -- get 2 random nodes
        get2rands(&vs, &vt, net.n);
        //cout << "vs = " << vs << ", vt = " << vt << endl;

        // -- find a mincut
        ncut = net.mincut(vs, vt, v);
        
        // -- make a exchange set
        for (i = 0; i < v.n; i++) {
            exchange_set.addset(bbi.bb[ net.map[v.elem[i]] ]);
        }
        // -- exchange strings
        exchange_set.exchange(p1, p2, c1, c2, ell);
    } else if (net.n == 2) {
        int a = myRand.uniformInt(0, 100)%2;
        exchange_set.addset(bbi.bb[net.map[a]]);
        exchange_set.exchange(p1, p2, c1, c2, ell);
    } else {
        // exchange_set = {}
        exchange_set.exchange(p1, p2, c1, c2, ell);
    }
}

int Network::mincut(int vs, int vt, Vector &v) {
    int i, j, k, ii;
    int srt = -1;
    int *prm = new int [n];
    int *delta = new int [n + 1];
    Vector L(n + 1), T(n + 1), M(n + 1), M2(n + 1); // T = L-S, M=^L, M2 temporary use
    M2.init();

    myRand.uniformArray(prm, n, 0, n - 1);

    while (1) {
        //- STEP1: labeling method to find flow-ZOUKARO

        // labeling method 0: Let L = {vs}, S = {}
        L.init();
        L.add(vs);
        T.init();
        T.add(vs); // T = L-S
        M.init();

        for (ii = 0; ii < n; ii++) {
            i = prm[ii];
            if (L.is_in(i) == 0) {
                M.add(i);
            } // M = ^L
        }
        
        while (!L.is_in(vt) && T.n > 0) {
            // labeling method 1: Stop if vt \in L or L==S
            //                    otherwise select i (\in L-S=T) and let S=S+{i}
            i = T.elem[myRand.uniformInt(0, 10000) % T.n];
            T.del(i);
            srt = L.n;
			
            //  labeling method 2: if i\notin L then L=L+{j} (s.t Eedge(i,j) exist)
            //                     and delta[j] = i
            //         	     return labeling method 1
            M2.copy(M);
            for (k = 0; k < M2.n; k++) {
                j = M2.elem[k];
                //cout << i << " " << j << endl;
                if (E[i][j] && (F[i][j] < E[i][j])) {
                    //cout << "flow pushed" << endl;
                    L.add(j);
                    M.del(j); // M = ^L
                    T.add(j);
                    delta[j] = i;

                }
            }
        }
        //- BREAK if L doesn't increase
        if (srt == L.n) {
            break;
        }
        //- STEP 2: add flow according to the delta[]
        //          and find new flow F, then return STEP 1
        j = vt;
        while (j != vs) {
            i = delta[j];
            F[i][j] = 1;
            j = i;
        }
    }

    if (L.n > M.n) {
        v.copy(M);
    } else {
        v.copy(L);
    }

    // add nodes without any edges randomly
    for (i = 0; i < n; i++) {
        k = 0;
        for (j = 0; j < n; j++) {
            k += E[i][j];
        }
        if (k == 0 && myRand.uniform() < 0.5)
            v.add(i);
    }

    int ncuts = 0;
    for (i = 0; i < L.n; i++) {
        for (j = 0; j < M.n; j++) {
            if (E[L.elem[i]][M.elem[j]]) ncuts++;
        }
    }

    delete [] prm;
    delete [] delta;

    return ncuts;
}

void DSMGA::simpleMinCut() 
{
	int i;
    getBBI();

    Network onet;
    onet.alloc(bbi.bbNum);
    onet.init();
    onet.alloc_map(bbi.bbNum);
    construct_the_original_graph(onet);

    if (nCurrent % 2 == 0) { // nCurrent is even
        for (i = 0; i < nCurrent; i += 2) {
            cdc_for_1pair_s(population[selectionIndex[i]],
                    population[selectionIndex[i + 1]], offspring[i],
                    offspring[i + 1], onet);
        }
    } else {
        for (i = 0; i < nCurrent - 1; i += 2) {
            cdc_for_1pair_s(population[selectionIndex[i]],
                    population[selectionIndex[i + 1]], offspring[i],
                    offspring[i + 1], onet);
        }
        offspring[nCurrent - 1] = population[selectionIndex[nCurrent - 1]];
    }
    //cout << "A generation of crossovers are done" << endl;
}

void DSMGA::crossOver() {
    int i;
    getBBI();

    Network onet;
    onet.alloc(bbi.bbNum);
    onet.init();
    onet.alloc_map(bbi.bbNum);
    construct_the_original_graph(onet);

    if (nCurrent % 2 == 0) { // nCurrent is even
        for (i = 0; i < nCurrent; i += 2) {
            cdc_for_1pair(population[selectionIndex[i]],
                    population[selectionIndex[i + 1]], offspring[i],
                    offspring[i + 1], onet);
        }
    } else {
        for (i = 0; i < nCurrent - 1; i += 2) {
            cdc_for_1pair(population[selectionIndex[i]],
                    population[selectionIndex[i + 1]], offspring[i],
                    offspring[i + 1], onet);
        }
        offspring[nCurrent - 1] = population[selectionIndex[nCurrent - 1]];
    }
    //cout << "A generation of crossovers are done" << endl;
}


void DSMGA::onePointXO ()
{
    int i;

    if ((nCurrent & 0x1) == 0) {
    	// nNextGeneration is even

        for (i = 0; i < nCurrent; i += 2)
            onePointXO (population[selectionIndex[i]], population[selectionIndex[i + 1]],
                offspring[i], offspring[i + 1]);

    }
    else {
        for (i = 0; i < nCurrent - 1; i += 2) {
            onePointXO (population[selectionIndex[i]], population[selectionIndex[i + 1]],
                offspring[i], offspring[i + 1]);
        }
        offspring[nCurrent - 1] =
            population[selectionIndex[nCurrent - 1]];
    }

}


void DSMGA::onePointXO (const Chromosome & p1, const Chromosome & p2, Chromosome & c1, Chromosome & c2)
{
    int i;
    int crossSite = myRand.uniformInt(1, ell-1);

    for (i = 0; i < crossSite; i++) {
            c1.setVal (i, p1.getVal(i));
            c2.setVal (i, p2.getVal(i));
    }

    for (i = crossSite; i < ell; i++) {
            c1.setVal (i, p2.getVal(i));
            c2.setVal (i, p1.getVal(i));
    }
}


void DSMGA::uniformXO ()
{
    int i;

    if ((nCurrent & 0x1) == 0) {
    	// nNextGeneration is even

        for (i = 0; i < nCurrent; i += 2)
            uniformXO (population[selectionIndex[i]], population[selectionIndex[i + 1]],
                offspring[i], offspring[i + 1], 0.5);

    }
    else {
        for (i = 0; i < nCurrent - 1; i += 2) {
            uniformXO (population[selectionIndex[i]], population[selectionIndex[i + 1]],
                offspring[i], offspring[i + 1] , 0.5);
        }
        offspring[nCurrent - 1] =
            population[selectionIndex[nCurrent - 1]];
    }

}

void DSMGA::uniformXO (const Chromosome & p1, const Chromosome & p2, Chromosome & c1, Chromosome & c2, double prob)
{
    int i;

    for (i = 0; i < ell; i++) {
        if (myRand.flip (prob)) {
            c1.setVal (i, p1.getVal(i));
            c2.setVal (i, p2.getVal(i));
        }
        else {
            c1.setVal (i, p2.getVal(i));
            c2.setVal (i, p1.getVal(i));
        }
    }
}



void DSMGA::getMutualInformation ( double** mi) {

    FastCounting *fastCounting =  new FastCounting[ell];

    for (int i = 0; i < ell; i++)
        fastCounting[i].init(nCurrent);

    for (int i = 0; i < nCurrent; i++)
        for (int j = 0; j < ell; j++)
            fastCounting[j].setVal(i, population[selectionIndex[i]].getVal(j));

    for ( int x = 0; x < ell ; ++x )
    {
        for ( int y = x+1; y < ell; ++ y)
        {
            int k;

            double q00,q01,q10,q11;

            int n00, n01, n10, n11;
            n00 = n01 = n10 = n11 = 0;

            for (k = 0; k < fastCounting[0].lengthLong; k++)
            {
                unsigned long val1 = fastCounting[x].gene[k];
                unsigned long val2 = fastCounting[y].gene[k];

                unsigned long long01 = (~val1) & (val2);
                unsigned long long10 = (val1) & (~val2);
                unsigned long long11 = (val1) & (val2);

                n01 += myBD.countOne(long01);
                n10 += myBD.countOne(long10);
                n11 += myBD.countOne(long11);
            }

            n00 = nCurrent - n01 - n10 - n11;

            q00 = (double)n00/(double)nCurrent;
            q01 = (double)n01/(double)nCurrent;
            q10 = (double)n10/(double)nCurrent;
            q11 = (double)n11/(double)nCurrent;

            mi[y][x]= mi[x][y] = mutualInformation(q00, q01, q10, q11);
        }
    }

    delete [] fastCounting;
}


void DSMGA::gibbsXO() {

    const int totalIter = 300;
    const int maxNeighbors = 13;
    //const double entropyThreshold = 0.8;
    const int groupNum = 4;

    //strengthBasedBayesianXO_TL();
    //crossOver();

    //getBBI ();


    //bbi.printOut();
    /*
	 //cout << "population after selection: \n";
	 for ( int i = 0; i < nCurrent; ++i )
	 {
	 population[selectionIndex[i]].printOut();
	 cout << endl;
	 }
	 */

    //system("pause");

    //  BBI should be obtained before this line
    /*
    if(bbi.bbNum==0)
    {
        cout << "No BB information, do normal crossover.\n";
        return crossOver();
    }
    */

    vector<int>* gene2Gene = new vector<int>[ell];            // gene to neighbors

    bool** adjMatrix = new bool*[ell];
    for ( int i = 0; i < ell; ++i )
    {
        adjMatrix[i] = new bool[ell];
    }
    for ( int i = 0; i < ell; ++i )
    {
        for ( int j = 0; j < ell; ++j )
        {
            adjMatrix[i][j] = true;
        }
        adjMatrix[i][i] = false;
    }
/*
    for ( int i = 0; i < ell; ++i )
    {
        for ( int j = 0; j < ell; ++j )
        {
            adjMatrix[i][j] = false;
        }
    }

    for ( int i = 0; i < bbi.bbNum; ++i )
    {
        for ( int j = 1; j < bbi.bb[i][0]+1 ; ++j )
        {
            for ( int k = j+1 ; k<bbi.bb[i][0]+1; ++k )
            {
                adjMatrix[ bbi.bb[i][j] ][ bbi.bb[i][k] ]  = true;
                adjMatrix[ bbi.bb[i][k] ][ bbi.bb[i][j] ]  = true;
            }
        }
    }
*/
    for ( int i = 0 ; i < ell; ++i )
    {
        for ( int j = 0; j < ell; ++j )
        {
            if ( adjMatrix[i][j] )
                gene2Gene[i].push_back(j);
        }
    }

    for ( int i = 0; i < ell; ++i )
    {
        delete [] adjMatrix[i];
    }

    delete [] adjMatrix;





//    	     //for test
//	    cout << "gene2BB :  \n";
//	    for (int i = 0; i < ell; ++i )
//	    {
//                for ( int j = 0; j < gene2BB[i].size(); ++j )
//                {
//                    cout << gene2BB[i][j] << "  ";
//                }
//                cout <<endl;
//	    }
//            system("pause");
//	     //end for test



    // for test
//    printf("before\n");
//
//    for ( int i = 0; i < ell ; ++i )
//    {
//        printf("total : %d    ", gene2Gene[i].size());
//        for ( int j = 0; j < gene2Gene[i].size(); ++j )
//        {
//            printf( "%d   ", gene2Gene[i][j]);
//        }
//        printf("\n");
//    }



    // gene2Gene size constraint  due to memory size limit
    double * entropy = new double[ell];

    for ( int i = 0; i < ell ; ++i )
    {
        entropy[i] = 0;
        int freq0 = 0;
        for ( int j=0; j < nCurrent;  ++j )
        {
            if (population[selectionIndex[j]].getVal(i)==0)
                freq0++;
        }
        if (freq0 != 0)
        {
            double p0 = double(freq0) / double(nCurrent);
            entropy[i] = -1 *p0 * log2(p0);
        }
        if (freq0 != nCurrent)
        {
            double p1 = double(nCurrent-freq0) / double(nCurrent);
            entropy[i] = entropy[i] - p1 * log2(p1);
        }
    }

    // for test
//    for ( int i = 0; i < ell; ++i )
//    {
//        printf("gene %d   entropy = %f \n", i, entropy[i]);
//    }


    // calculate mutual information

    double ** mi = new double*[ell];
    for ( int i = 0; i < ell; ++i )
    {
        mi[i] = new double[ell];
    }

    getMutualInformation(mi);


    int neiSize = 0;
    for ( int i = 0 ; i < ell ; ++i )
    {
        if  ( gene2Gene[i].empty() )
            continue;

        neiSize = gene2Gene[i].size();

//        double* tempEntropy = new double[neiSize] ;
//        int * bigFirst = new int [neiSize];
//        for ( int j = 0 ; j < neiSize; ++j )
//        {
//            tempEntropy[j] = entropy[gene2Gene[i][j]];
//            bigFirst[j] = j;
//        }
//
//        QuickSort(bigFirst, tempEntropy, 0, neiSize-1);
//        vector<int> temp;
//        for ( int j = neiSize-1; j > neiSize-1- maxNeighbors; --j )
//        {
//            if ( tempEntropy[bigFirst[j]] < entropyThreshold )
//                temp.push_back( gene2Gene[i][bigFirst[j]] );
//            else
//                break;
//            if ( j == 0 )
//                break;
//        }

        double* tempEntropy = new double[neiSize] ;
        int * bigFirst = new int [neiSize];
        for ( int j = 0 ; j < neiSize; ++j )
        {
            tempEntropy[j] = mi[i][ gene2Gene[i][j] ];
            bigFirst[j] = j;
        }

        QuickSort(bigFirst, tempEntropy, 0, neiSize-1);
        vector<int> temp;
        for ( int j = 0; j < maxNeighbors; ++j )
        {
            if ( j == neiSize )
                break;
            temp.push_back( gene2Gene[i][bigFirst[j]] );
        }



        gene2Gene[i] = temp;
        delete [] tempEntropy;
        delete [] bigFirst;

    }

    delete [] entropy;

    // for test
//    printf("after :\n");
//
//    for ( int i = 0; i < ell ; ++i )
//    {
//        printf("total : %d    ", gene2Gene[i].size());
//        for ( int j = 0; j < gene2Gene[i].size(); ++j )
//        {
//            printf( "%d   ", gene2Gene[i][j]);
//        }
//        printf("\n");
//    }
//
//    getchar();



//
//
    // calculate genes' potentials
    double** potentials = new double*[ell];     // potential for gene == 1
    long int totalPattern = 0;


    for ( int i = 0; i < ell; ++i )
    {

        totalPattern = (long int) pow(2.0, (double)gene2Gene[i].size());
        potentials[i] = new double[totalPattern];
//        cout << "totalPattern : " << totalPattern <<endl;
        double* total = new double[totalPattern];
        for ( long int j = 0; j < totalPattern; ++j )
        {
            potentials[i][j] = 1/ ( generation+1) /2 + 1e-6;
            total[j] = potentials[i][j] *2 ;
        }

        for ( int j=0; j < nCurrent; ++j )
        {
            long int pattern = population[selectionIndex[j]].makeInt(gene2Gene[i]);
//            cout << "pattern : " << pattern <<endl;
            if (population[selectionIndex[j]].getVal(i) == 1)
            {
                ++potentials[i][pattern];
            }
            ++total[pattern];
        }

//        printf("potential  of %d    ", i);
        for ( long int j = 0; j < totalPattern; ++j )
        {
            potentials[i][j] = potentials[i][j] / total[j];
//            printf("%f  ",  potentials[i][j] );
        }
//        printf("\n");

        delete [] total;

    }
//        cout <<"1\n";
//    for ( int i = 0; i < ell; ++i )
//    {
//        cout << i << "  ";
//        delete [] potentials[i];
//    }


    // random initialize

    for (int i = 0 ; i < groupNum; ++i)
    {
        for ( int j = 0; j < ell;  ++j)
        {
            if (myRand.flip())
                offspring[i].setVal(j, 1);
            else
                offspring[i].setVal(j, 0);
        }
    }


    for (int k = 0 ; k < groupNum; ++k )
    {
        for ( int iter = 0 ; iter < totalIter; ++iter )
        {
            for ( int i = 0; i < ell; ++i )
            {
                long int pattern = offspring[k].makeInt(gene2Gene[i]);
                if (myRand.flip(potentials[i][pattern]))
                    offspring[k].setVal(i, 1);
                else
                    offspring[k].setVal(i, 0);
            }
        }
    }

    for ( int p = groupNum ; p < nCurrent; ++p )
    {
        for ( int i = 0; i < ell; ++i )
        {
            long int pattern = offspring[p%groupNum].makeInt(gene2Gene[i]);
            if (myRand.flip(potentials[i][pattern]))
                offspring[p%groupNum].setVal(i, 1);
            else
                offspring[p%groupNum].setVal(i, 0);
        }
        offspring[p] = offspring[p%groupNum];
    }

    for ( int i = 0; i < ell; ++i )
    {
        delete [] mi[i];
        delete [] potentials[i];
    }
    delete [] mi;
    delete [] gene2Gene;
    delete [] potentials;


}

void DSMGA::minCut () {
  int i;
  getBBI ();

  if ((nCurrent & 0x1) == 0)
    {				// nNextGeneration is even
      for (i = 0; i < nCurrent; i += 2)
	{
	  minCut (population[selectionIndex[i]],
		     population[selectionIndex[i + 1]], offspring[i],
		     offspring[i + 1]);
	}
    }
  else
    {
      for (i = 0; i < nCurrent - 1; i += 2)
	{
	  minCut (population[selectionIndex[i]],
		     population[selectionIndex[i + 1]], offspring[i],
		     offspring[i + 1]);
	}
      offspring[nCurrent - 1] =
	population[selectionIndex[nCurrent - 1]];
    }

}

void DSMGA::minCut (Chromosome & p1, Chromosome & p2, Chromosome & c1, Chromosome & c2) {
    if (myRand.uniform () < pc)
    {
        //minimalCut1(p1,p2,c1,c2);
        minimalCut2(p1,p2,c1,c2);
        /*
        double prob = 0.5;
        bbUniformXO (p1, p2, c1, c2, prob);
        */
    }
    else
    {
        c1 = p1;
        c2 = p2;
    }
}


void DSMGA::minimalCut2 (Chromosome &p1, Chromosome &p2, Chromosome &c1, Chromosome &c2) {
  int i;
  int edge = (int) sqrt (ell + 0.5);
  int r;

  double ppp = double (edge % 4) / 4.0;

  if (myRand.flip(ppp))
	  r = (edge / 4) + 1;
  else
	  r = (edge / 4);


  int x = myRand.uniformInt(0, edge-1);
  int y = myRand.uniformInt(0, edge-1);

  for (i=0; i<ell; i++) {
	  int xx = i / edge;
	  int yy = i % edge;
	  if (((x <= xx && xx <= x+r-1) || (x+r-1 > edge && xx <= x+r-1-edge)) &&
	     ((y <= yy && yy <= y+r-1) || (y+r-1 > edge && yy <= y+r-1-edge))) {
				  c1.setVal (i, p1.getVal(i));
				  c2.setVal (i, p2.getVal(i));
	  }
	  else {
				  c1.setVal (i, p2.getVal(i));
				  c2.setVal (i, p1.getVal(i));
	  }

  }
}

