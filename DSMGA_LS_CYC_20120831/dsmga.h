
#ifndef _DSMGA_H
#define _DSMGA_H

/**
@author Tian-Li Yu
 */

#include "global.h"
#include "chromosome.h"
#include "bbi.h"
#include "statistics.h"
#include "dsmclusteringchromosome.h"
#include "vector.h"
#include <iostream>

class Network {
public:
    int n; // the number of nodes
    int n_initial;
    int **E; // edges
    int **F; // flows
    int **edge_list;
    int max_incomes; // the maximum number of edges in a node
    int *map; // mappings from nodes to linkage sets

    Network() {
        E = NULL;
        F = NULL;
        map = NULL;
        edge_list = NULL;
    }

    ~Network() {
        free();
    }

    Network(int n1) {
        n = n1;
        E = NULL;
        F = NULL;
        map = NULL;
        edge_list = NULL;
        alloc();
        init();
    }

    void alloc() {
        n_initial = n;
        int i;

        if (E) {
            for (i = 0; i < n; i++)
                if (E[i]) delete[] E[i];
            delete [] E;
        }
        if (F) {
            for (i = 0; i < n; i++)
                if (F[i]) delete[] F[i];
            delete [] F;
        }

        E = new int * [n];
        F = new int * [n];
        for (i = 0; i < n; i++) {
            E[i] = new int [n];
            F[i] = new int [n];
        }
    }

    void alloc(int n1) {
        n = n1;
        alloc();
    }

    void free() {
        int i;
        for (i = 0; i < n_initial; i++)
            delete[] E[i];
        delete [] E;
        for (i = 0; i < n_initial; i++)
            delete[] F[i];
        delete [] F;
		if (edge_list) {
            for (i = 0; i < n_initial; i++)
                delete[] edge_list[i];
            delete [] edge_list;
        }
        delete[] map;
    }

    void init() {
        int i, j;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                E[i][j] = 0;
                F[i][j] = 0;
            }
        }
    }

    void add_edge(int i, int j) {
        E[i][j] = 1;
    }

    void alloc_map(int n) {
        if (map != NULL) {
            delete [] map;
        }
        map = new int [n];
    }
    int mincut(int vs, int vt, Vector &v);

    void print_edges() {
        int i, j;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++)
                printf("%d", E[i][j]);
            printf("\n");
        }
    }
};

class DSMGA {
public:
    DSMGA(int n_ell, int n_nInitial, int n_selectionPressure, double n_pc,
            double n_pm, int n_maxGen, int n_maxFe);

    ~DSMGA();

    void init(int n_ell, int n_nInitial, int n_selectionPressure, double n_pc,
            double n_pm, int n_maxGen, int n_maxFe);

    void initializePopulation();
		void localSearch();	
    void evaluate();

    void selection();

    /** tournament selection without replacement*/
    void tournamentSelection();

    int countOverlap(Chromosome &, Chromosome &);
    void crossOver();
    void onePointXO ();
    void onePointXO (const Chromosome & p1, const Chromosome & p2, Chromosome & c1, Chromosome & c2);
    void uniformXO();
    void uniformXO (const Chromosome & p1, const Chromosome & p2, Chromosome & c1, Chromosome & c2, double prob);
    void crossOver(Chromosome &, Chromosome &, Chromosome &, Chromosome &);
    void bbUniformXO();
    void bbUniformXO(Chromosome &, Chromosome &, Chromosome &, Chromosome &,
            double);
    void bfs(int **cap, int **flow, bool *splitter, int **edges, int node);
    void crossOverMone(Chromosome &p1, Chromosome &p2, Chromosome &c1, Chromosome &c2);
    bool isOverlapping(int *bb1, int *bb2, Chromosome &p1, Chromosome &p2);
    void mutation();
    void bbMutation(Chromosome &);

    /** get BBI via DSM clustering */
    void getBBI();
    void getBBIOneMax();
    void getBBIMKTrap();
    BBI getBBIMune();
    //void oneRun (bool output = true);
    void oneRun(bool output = true, int mode = 1);
    //int doIt (bool output = true);
    //int doIt(bool output = true, int mode = 1);
    int doIt(bool output = true, int mode = 1);

    bool shouldTerminate();
    int getNextPopulation();

    bool foundOptima();

    int getGeneration() const {
        return generation;
    }

    void replacePopulation();
    void fullReplace();
    void RTR();

    void showStatistics();

    int findInitialPopulation();
    bool findCorrectModel(int requiredMatchedM);

    bool isSteadyState();

    void strengthBasedBayesianXO_TL();
    double getBBStrength(int *bb);
    void calculateBBIError() ;

    void minCut ();
    void minimalCut2 (Chromosome &p1, Chromosome &p2, Chromosome &c1, Chromosome &c2);
    void minCut (Chromosome & p1, Chromosome & p2, Chromosome & c1, Chromosome & c2);


    void gibbsXO();
    void getMutualInformation ( double** mi);


    void construct_the_original_graph(Network &onet);
    void reconstruct_graph(Chromosome &p1, Chromosome &p2, Network &oNet, Network &nNet);
	void reconstruct_graph_s(Chromosome &p1, Chromosome &p2, Network &oNet, Network &nNet);
    void cdc_for_1pair(Chromosome &p1, Chromosome &p2, Chromosome &c1, Chromosome &c2, Network &onet);
	void cdc_for_1pair_s(Chromosome &p1, Chromosome &p2, Chromosome &c1, Chromosome &c2, Network &onet);
	void simpleMinCut();
    bool sameSubstrings(Chromosome &p1, Chromosome &p2, int *bb);
	double getMaxFitness() {return stFitness.getMax();}
	int getbest(int i) {return population[bestIndex].getVal(i);}
	double lsNFE;
protected:

    BBI bbi; // BB Information
    bool bbiCalculated;

    int ell; // chromosome length
    int nInitial; // initial population size
    int nCurrent; // current population size
    int selectionPressure;
    
    double pc; // prob of XO
    double pm; // prob of Mutation
    Chromosome *population;
    Chromosome *offspring;
    int *selectionIndex;
    int maxGen;
    int maxFe;
    int repeat;
    int fe;
    int generation;
    int bestIndex;
    
    Statistics stFitness;
    int getDistance(Chromosome & c1, Chromosome & c2);
    int getBBDistance(Chromosome & c1, Chromosome & c2);

    double previousFitnessMean;

};
#endif
