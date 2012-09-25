/*
 * File:   UFO.cpp
 * Author: Rick
 *
 * Created on January 17, 2011, 7:45 PM
 */

#include "UFO.h"
#include "math.h"
#include "float.h"
#include <limits.h>
#include <cassert>

UFO::UFO()
{
    m = k = n = w = 0;
    bb2gene = NULL;
    gene2bb = NULL;
    graph = NULL;
    room = NULL;
    bbSize = NULL;
    geneSize = NULL;
    geneOmega = NULL;
    bbType = NULL;
    g2B.reset();
    b2B.reset();
    zeros = 0;
    unUsedGenes = 0;
    globalOptimum = 0;
    conflictedBB = 0;
    bestChromosome = NULL;
}

UFO::UFO(const UFO& orig)
{
    reset();
    m = orig.m;
    k = orig.k;
    n = orig.n;
    w = orig.w;
    zeros = orig.zeros;
    
        
}

UFO::UFO( const int& m, const int& k, const double& desiredW)
{
    bb2gene = NULL;
    gene2bb = NULL;
    graph = NULL;
    room = NULL;
    bbSize = NULL;
    geneSize = NULL;
    geneOmega = NULL;
    evaluated = false;
    bbType = NULL;
    g2B.reset();
    b2B.reset();
    zeros = 0;
    unUsedGenes = 0;
    globalOptimum = 0;
    conflictedBB = 0;


    construct(m,k,desiredW);

}


UFO::~UFO()
{
    reset();
}

void UFO::reset()
{
    if ( bb2gene != NULL ) {
        for ( int i = 0; i < m; ++i ) {
            delete [] bb2gene[i];
        }

        for ( int i = 0; i < n; ++i ) {
            delete [] gene2bb[i];
        }

        for ( int i = 0; i < m+n+2; ++i ) {
            delete [] room[i];
        }
        delete [] room;
        delete [] bb2gene;
        delete [] gene2bb;

        delete [] graph;

        delete [] bbSize;
        delete [] geneSize;
        delete [] geneOmega;
        delete [] bbType;
    }
    m = k = n = w = 0;
    bb2gene = NULL;
    gene2bb = NULL;
    graph = NULL;
    room = NULL;
    bbSize = NULL;
    geneSize = NULL;
    geneOmega = NULL;
    evaluated = false;
    bbType = NULL;
    g2B.reset();
    b2B.reset();
    zeros = 0;
    unUsedGenes = 0;
    globalOptimum = 0;
    conflictedBB = 0;
    
    delete [] bestChromosome;
    bestChromosome = NULL;
}


bool UFO::createProblem(const int& n, const int& k, const double& desiredW, int confliction)
{
    reset();

    bool ok = construct(n,k,desiredW);
    
    if ( !ok )
        printf("NOT OK\n");
    
    if ( confliction != 0 && ok) {
        for ( int i = 0; i < 10; ++i ) {
            if ( randomConflict(confliction) == confliction ||  confliction == -1 ) {
                checkAndSetGlobalOptimum();
                return true;
            } else {
                reset();
                construct(n,k,desiredW);
            }
        }
        return false;
    } else if ( ok) {
        checkAndSetGlobalOptimum();
        return ok;
    } else
        return ok;
}


bool UFO::initial(const int& tempN, const int& tempK, const double& tempDesiredW)
{

    if ( tempDesiredW == 0) {
        printf("Desired omega can't be 0. \n");
        return false;
    }

    reset();

    desiredW = tempDesiredW;
    double tempM= (double)tempN*tempDesiredW /(double)tempK;

    n = tempN;
    //printf("hello , n = %d\n",n);
    m = (int)tempM;
    k = tempK;


    bestChromosome = new int[n];


    if ( k > n ) {
        printf("k = %d   n = %d  \n", k, n);
        printf("The setting is wrong!!\n");
        return false;
    }

    bb2gene = new int*[m];
    bbSize = new int[m];
    for ( int i = 0; i < m; ++i ) {
        bb2gene[i] = new int[k];
        bbSize[i] = 0;
    }


    geneSize = new int[n];
    geneOmega = new int[n];
    for ( int i = 0; i < n; ++i ) {
        geneSize[i] = 0;
        geneOmega[i] = 0;
    }


    // decide every gene's omega
    int baseW = m*k/n;
    int leftFlow = m*k - baseW*n;
    for ( int i = 0; i < n; ++i ) {
        geneOmega[i] = baseW;
    }
    for ( int i = 0; i < leftFlow; ++i ) {
        geneOmega[i] = geneOmega[i] + 1;
    }

    gene2bb = new int*[n];
    for ( int i = 0; i < n; ++i ) {
        gene2bb[i] = new int[geneOmega[i]];
    }

    constructGraph();

    bbType = new int[m];
    for ( int i = 0; i < m; ++i ) {
        bbType[i] = 1; // One-Trap
    }

    return true;
}

bool UFO::construct(const int& tempN, const int& tempK, const double& tempDesiredW)
{

    if ( ! initial(tempN, tempK, tempDesiredW) ) {
        return false;
    }

    while (this->BFSandUpdate());



    for ( int i = 1; i < 1+m; ++i ) {
        int bbIndex = i -1;
        for ( int j = m+1; j < m+1+n; ++j ) {
            int geneIndex = j - (m+1);
            if ( room[i][j] == 0 ) {
                bb2gene[bbIndex][bbSize[bbIndex]] = geneIndex;
                gene2bb[geneIndex][geneSize[geneIndex]] = bbIndex;
                ++bbSize[bbIndex];
                ++geneSize[geneIndex];
            }
        }
    }


    return true;
}


// the graph is a flow graph.   index of s = 0,   bb = 1 ~ m,   gene = m +1 ~ m+n,  t = m+n+1
void UFO::constructGraph()
{
    graph = new list<int>[2+m+n];

    // initial room matrix;
    int total = m+n+2;
    room = new int*[total];
    for ( int i = 0; i < total; ++i ) {
        room[i] = new int[total];
    }
    for ( int i = 0; i < total; ++i ) {
        for ( int j = 0; j < total; ++j ) {
            room[i][j] = 0;
        }
    }



    // from source to every bb
    int* temp = new int[m];
    r.uniformArray(temp, m, 1, m);
    for ( int i = 0; i < m; ++i ) {
        //printf("temp %d = %d\n", i, temp[i]);
        graph[0].push_back(temp[i]);
    }



    for ( int i = 1; i < 1+m; ++i ) {
        room[0][i] = k;
    }
    delete [] temp;

    // from bb to every gene
    for ( int i = 1; i < 1+m; ++i ) {
        int* temp = new int[n];
        r.uniformArray(temp, n, m+1, m+n);
        for ( int j = 0; j < n; ++j ) {
            //printf("temp %d = %d\n", j, temp[j]);
            graph[i].push_back(temp[j]);
        }
        for ( int j = m+1; j < m+1+n; ++j ) {
            room[i][j] = 1;
        }
        delete [] temp;
    }

    // from gene to target
    for ( int i = m+1; i <m+1+n; ++i ) {
        graph[i].push_back(m+n+1);
        room[i][m+n+1] = geneOmega[i-m-1];
    }

}


bool UFO::BFSandUpdate()
{
    int total = m+n+2;

    int* color = new int[total];    // white = 0, gray = 1, black = -1
    int* parent = new int[total];

    for ( int i = 0; i < total; ++i ) {
        parent[i] = -1;
        color[i] = 0;
    }

    color[0] = 1;

    list<int> queue;
    queue.push_back(0);

    while ( ! queue.empty() ) {
        int u = queue.front();
        int v = 0;
        for ( list<int>::iterator iter = graph[u].begin(); iter != graph[u].end(); ++iter ) {
            v = *iter;
            if ( color[v] == 0 ) {
                color[v] = 1;
                parent[v] = u;
                queue.push_back(v);
            }
        }
        queue.pop_front();
        color[u] = -1;
    }

    if ( color[m+n+1] != -1 ) {
        delete [] color;
        delete [] parent;
        return false;
    }

    // trace back and update edge
    traceBackAndUpdate(parent);




    delete [] color;
    delete [] parent;
    return true;
}


void UFO::traceBackAndUpdate(int* parent)
{

    // find path
    list<int> stack;
    int index = m+n+1;
    int minRoom = INT_MAX;

    stack.push_front(index);
    int p = 0;
    while ( index != 0 ) {
        p = parent[index];
        if ( room[p][index] < minRoom)
            minRoom = room[p][index];
        index = p;
        stack.push_front(index);
    }

    // update room matrix and graph
    // index shound = 0
    stack.pop_front();
    int child = 0;

    while ( index != m+n+1 ) {
        child = stack.front();
        room[index][child] = room[index][child] - minRoom;
        assert(room[index][child]>= 0);
        room[child][index] = room[child][index] + minRoom;

        if ( room[index][child] == 0 ) {
            graph[index].remove(child);
            //graph[child].push_back(index);
            // for random
            int at = r.uniformInt(0,graph[child].size());
            list<int>::iterator iter = graph[child].begin();
            for ( int i = 0; i < at; ++i )
                ++iter;
            graph[child].insert(iter, index);
        }
        index = child;
        stack.pop_front();
    }


}


void UFO::print()
{
    if ( !this->evaluated )
        this->evaluate();

    printf("Desired omega = %f\n", desiredW);
    printf("m = %d\n", m);
    printf("k = %d\n", k);
    printf("n = %d\n", n);
    printf("w = %d\n", w);

    printStatistics();

    printf("BB : \n");
    for ( int i = 0; i < m; ++i ) {
        printf("%d  (%d)[%d]   ", i, bbType[i], bbSize[i]);
        for ( int j = 0; j < bbSize[i]; ++j ) {
            printf("%d  ", bb2gene[i][j]);
        }
        printf("\n");
    }

    printf("Gene to BB : \n");
    for ( int i = 0; i < n; ++i ) {
        int total = 0;
        printf("%d  [%d]   ", i, geneSize[i]);
        for ( int j = 0; j < geneSize[i]; ++j ) {
            printf("%d  ", gene2bb[i][j]);
            if ( bbType[gene2bb[i][j]] == 0)
                ++total;
        }
        printf("  (%d)", total);
        printf("\n");
    }

    printf("Zeros : %d\n", this->zeros);
    printf("ConflictedBBs : %d\n", this->conflictedBB);
    printf("Global optimum : %f\n", this->globalOptimum);
    printf("Best Chromosome : ");
    for ( int i = 0; i < n; ++i )
        printf("%d",bestChromosome[i]);
    printf("\n");
}

void UFO::printStatistics()
{
    if ( ! evaluated ) {
        evaluate();
    }

    printf("g2B:  min/mean/max =  %f / %f / %f \n", g2B.getMin(),g2B.getMean(),g2B.getMax());
    printf("b2B:  min/mean/max =  %f / %f / %f \n", b2B.getMin(),b2B.getMean(),b2B.getMax());
    printf("Unused genes = %d\n", unUsedGenes);

}


void UFO::evaluate()
{

    // handle g2B b2B

    for ( int i = 0; i < n ; ++i ) {
        if (geneSize[i] != 0)
            g2B.record((double)geneSize[i]);
        else
            ++unUsedGenes;
    }

    bool** bbGraph = new bool*[m];
    for ( int i = 0; i < m; ++i ) {
        bbGraph[i] = new bool[m];
    }
    for ( int i = 0; i < m; ++i ) {
        for ( int j = 0; j < m; ++j ) {
            bbGraph[i][j] = false;
        }
    }

    for ( int i = 0; i < n; ++i ) {
        for ( int j =0; j < geneSize[i]; ++j ) {
            int thisBB = gene2bb[i][j];
            for ( int k = j+1; k < geneSize[i]; ++k ) {
                bbGraph[thisBB][gene2bb[i][k]] = true;
                bbGraph[gene2bb[i][k]][thisBB] = true;
            }
        }
    }

    int * count =  (int*) calloc (m,sizeof(int));
    for ( int i = 0; i < m ; ++i ) {
        for ( int j = 0; j < m; ++j ) {
            if (bbGraph[i][j])
                ++count[i];
        }
    }

    for ( int i = 0; i < m; ++i ) {
        //printf("BB %d  overlaps with %d  BBs\n", i, count[i]);
        b2B.record((double)count[i]);
    }

    for ( int i = 0; i < m; ++i )
        delete [] bbGraph[i];
    delete [] bbGraph;

    free(count);


    // create bbi

    BBI tempBBI(n);
    tempBBI.bbNum = m;
    for ( int i = 0; i < m; ++i ) {
        tempBBI.bb[i][0] = bbSize[i];
        for ( int j = 0; j < bbSize[i]; ++j ) {
            tempBBI.bb[i][j+1] = bb2gene[i][j];
        }
    }

    bbi = tempBBI;



    // calculate conflictedBB
    this->findConflictedBB();

    evaluated = true;
}

BBI UFO::getBBI()
{
    if ( ! evaluated ) {
        evaluate();
    }
    return bbi;
}


double UFO::kTrap( const int& k, const int& totalOne, const int& type) const
{
    if (totalOne > k)
        return 0.0;

    int u = totalOne;

    if ( type == 1 ) {
        if (u == k)
            return 1.0;
        else
            return 0.8/(double)(k-1)*(double)(k-1-u);
    } else { // type == 0
        if (u == k)
            return 0.8;
        else
            return 1/(double)(k-1)*(double)(k-1-u);
    }

}


double UFO::kTrap9( const int& k, const int& totalOne, const int& type) const
{
    if (totalOne > k)
        return 0.0;

    int u = totalOne;

    if ( type == 1 ) {
        if (u == k)
            return 1.0;
        else
            return 0.9/(double)(k-1)*(double)(k-1-u);
    } else { // type == 0
        if (u == k)
            return 0.9;
        else
            return 1/(double)(k-1)*(double)(k-1-u);
    }
}

double UFO::kMax( const int& k, const int& totalOne, const int& type) const
{
    if (totalOne > k)
        return 0.0;

    int u = totalOne;

    if ( type == 1 ) {
        return 1*u/k;
    } else { // type == 0
        return 1-1*u/k;
    }
}

double UFO::ikTrap( const int& k, const int& totalOne, const int& type) const
{
    if (totalOne > k)
        return 0.0;

    int u = totalOne;

    if ( type == 1 ) {
        if (u == k)
            return 1.0;
        else
            return 0.8/(double)(k-1)*(double)(k-1-u);
    } else { // type == 0
        if (u == k)
            return 1.0-1.0;
        else
            return 1.0-0.8/(double)(k-1)*(double)(k-1-u);

    }

}




double UFO::getFitness( int const * chromosome ) const
{
    
    double fitness = 0;
    int totalOne = 0;
    for ( int i = 0; i < m; ++i ) {
        totalOne = 0;
        for ( int j = 0; j < bbSize[i]; ++j ) {
            if ( chromosome[bb2gene[i][j]] == 1 )
                ++totalOne;
        }
        
        fitness = fitness + kTrap9(bbSize[i], totalOne, bbType[i]);
        //fitness = fitness + kMax(bbSize[i],totalOne,bbType[i]);
        //fitness = fitness + ikTrap(bbSize[i], totalOne, bbType[i]);

    }

#ifdef CHECK_GLOBAL_WHEN_CALCULATATION
    if ( fitness > this->getMaxFitness() + 5e-6 ) {
        int ones = 0;
        for (int i = 0; i < n; ++i ) {
            if (chromosome[i] == 1)
                ++ones;
        }
        if ( ones != n ) {
            printf("\nOhoh~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
            printf("fitness = %f   ", fitness );
            printf("chromosome = ");
            for  ( int i = 0; i < n; ++i )
            {
                printf("%d", chromosome[i]);
            }
            printf("\n");
        }
    }
#endif


    return fitness;
}

int UFO::getNumOfGenes() const
{
    return n;
}

int UFO::getNumOfBBs() const
{
    return m;
}

double UFO::getMaxFitness() const
{

    return globalOptimum - 1e-6;
//    if ( ! use9 )
//        //return (double) ( m - zeros ) + (double) ( zeros ) *0.8 -1e-6;
//        return (double) ( m - zeros ) -1e-6;
//    else
//        return (double) ( m - zeros ) + (double) ( zeros ) *0.9 -1e-6;
//    //return m - 1e-6;
}


// must have been constructed
int UFO::conflict()
{
    list<int>* bg = new list<int>[m];
    list<int>* gb = new list<int>[n];

    for ( int i = 0; i < m; ++i ) {
        for ( int j = 0; j < bbSize[i]; ++j ) {
            bg[i].push_back(bb2gene[i][j]);
        }
    }
    for ( int i = 0; i < n; ++i ) {
        for ( int j = 0; j < geneSize[i]; ++j ) {
            gb[i].push_back(gene2bb[i][j]);
        }
    }

    int * cost = new int[m];
    for ( int i = 0; i < m; ++i ) {
        cost[i] = 2;
    }


    list<int> One;
    int * tempIndex = new int[m];
    r.uniformArray(tempIndex,m,0,m-1);
    for ( int i = 0; i < m; ++i ) {
        One.push_back(tempIndex[i]);
    }
    delete [] tempIndex;

    zeros = 0;



    for ( list<int>::iterator iter = One.begin(); iter != One.end(); ++iter ) {
        int thisBB = *iter;
        int count = 0;

        bool * view =  (bool*) calloc (m,sizeof(bool));
        for ( list<int>::iterator geneIter = bg[thisBB].begin(); geneIter != bg[thisBB].end(); ++geneIter ) {
            int thisGene = *geneIter;
            for ( list<int>::iterator bbIter = gb[thisGene].begin(); bbIter != gb[thisGene].end(); ++bbIter ) {
                if ( ! view[*bbIter]  ) {
                    if ( bbType[*bbIter] == 1 )
                        ++count;
                    view[*bbIter] = true;
                }
            }
        }

        if ( count > cost[thisBB] ) {
            bbType[thisBB] = 0;
            ++zeros;
            for ( int i = 0 ; i < m; ++i ) {
                if ( view[i] ) {
                    ++cost[i];
                }
            }

            for ( int i = 0; i < n; ++i ) {
                gb[i].remove(thisBB);
            }
//            iter = One.erase(iter);
//            --iter;
        }

        free(view);
    }

    delete [] cost;
    delete [] bg;
    delete [] gb;

    return zeros;
}

int UFO::getZeros() const
{
    return zeros;
}



int UFO::conflict2( const int& totalConflict)
{

    int* cost = new int[n];
    for ( int i = 0; i < n; ++i ) {
        cost[i] = 3;
    }
    zeros = 0;

    for ( int i = 0; i < m; ++i ) {
        bool canConflict = true;

        for ( int j = 0; j < bbSize[i]; ++j ) {
            int thisGene = bb2gene[i][j];
            int totalOne = 0;
            for ( int k = 0; k < geneSize[thisGene]; ++k ) {
                if ( bbType[gene2bb[thisGene][k]] ==1 )
                    ++totalOne;
            }
            if ( totalOne < cost[thisGene] ) {
                canConflict = false;
                break;
            }
        }

        if ( canConflict ) {
            bbType[i] = 0;
            ++zeros;
            if ( zeros == totalConflict )
                break;

            for ( int j = 0; j < bbSize[i]; ++j ) {
                ++cost[bb2gene[i][j]];
            }
        }

    }

    delete [] cost;

    return zeros;
}


double UFO::getMeanOmega()
{
    if ( ! evaluated )
        evaluate();
    return g2B.getMean();
}

double UFO::getMaxOmega()
{
    if ( ! evaluated )
        evaluate();
    return g2B.getMax();
}

double UFO::getMinOmega()
{
    if ( ! evaluated )
        evaluate();
    return g2B.getMin();
}

double UFO::getG2BStd()
{
    if ( ! evaluated )
        evaluate();
    return g2B.getStdev();
}

double UFO::getB2BStd()
{
    if ( ! evaluated )
        evaluate();
    return b2B.getStdev();
}


double UFO::getB2BMean()
{
    if ( ! evaluated )
        evaluate();
    return b2B.getMean();
}


int UFO::getUnusedGenes()
{
    if ( ! evaluated )
        evaluate();
    return unUsedGenes;
}

int UFO::getTotalConflictedBB()
{
    if ( ! evaluated )
        evaluate();
    return conflictedBB;
}

int UFO::getG2BMaxDifference()
{
    if ( ! evaluated )
        evaluate();
    return (int)g2B.getMax()-(int)g2B.getMin();
}


int UFO::conflict3( const int& totalConflict)
{


    zeros = 0;

    for ( int i = 0; i < m; ++i ) {
        //bool canConflict = true;

        bool* checked = (bool*)calloc(m,sizeof(bool));
        int totalOne = 0;
        int totalZero = 0;

        for ( int j = 0; j < bbSize[i]; ++j ) {
            int thisGene = bb2gene[i][j];
            for ( int k = 0; k < geneSize[thisGene]; ++k ) {
                int thisBB = gene2bb[thisGene][k];
                if ( !checked[thisBB] ) {
                    if ( bbType[thisBB] ==1  )
                        ++totalOne;
                    else
                        ++totalZero;
                    checked[thisBB] = true;
                }
            }
        }

        --totalOne;
        ++totalZero;

        if ( 0.2*totalOne > totalZero ) {
            bbType[i] = 0;
            ++zeros;
            if ( zeros == totalConflict ) {
                free(checked);
                break;
            }
        }
        free(checked);
    }

    return zeros;
}

bool UFO::checkGlobalOptimum()
{

    int c[n];
    double maxFitness = 0;
    int maxC[n];
    for ( int i = 0; i < n; ++i ) {
        c[i] = 0;
    }
    printf("check if when m = %d fmax  = %f \n", m, this->getMaxFitness());

    for ( int i = 0; i < pow(2.0,(double)n); ++i ) {
        for ( int q = 0; q < n; ++q ) {
            c[q] = 0;
        }
        int still = i;
        int index = 0;
        while (still > 0) {
            c[index] =  still%2;
            still = still /2;
            ++index;
        }
        //c.printOut();
        if ( getFitness(c) > maxFitness  ) {
            maxFitness = getFitness(c);
            for ( int i = 0; i < n; ++i )
                maxC[i] = c[i];
        }
    }


    int total = 0;
    for ( int i = 0; i < n; ++i ) {
        if ( maxC[i] == 1 )
            ++total;
    }

    if (maxFitness > this->getMaxFitness() && total != n) {
        printf("maxFitness = %f > %f \n",maxFitness,this->getMaxFitness());
        for ( int i = 0; i < n; ++i )
        {
            printf("%d", maxC[i]);
        }
        printf("\n");
        return false;
    }

    printf("matched! f = %f \n",maxFitness);
    for ( int i = 0; i < n; ++i )
    {
        printf("%d", maxC[i]);
    }
    printf("\n");

    return true;
}

int UFO::conflict4( const int& totalConflict)
{


    zeros = 0;

    int* cost = (int*) calloc(m, sizeof(int));

    for ( int i = 0; i < m; ++i ) {
        //bool canConflict = true;

        bool* checked = (bool*)calloc(m,sizeof(bool));
        int totalOne = 0;
        int totalZero = 0;

        int pseudo = 0;

        for ( int j = 0; j < bbSize[i]; ++j ) {
            int thisGene = bb2gene[i][j];
            for ( int k = 0; k < geneSize[thisGene]; ++k ) {
                int thisBB = gene2bb[thisGene][k];
                if ( !checked[thisBB] ) {
                    if ( bbType[thisBB] ==1  ) {
                        ++totalOne;
                    } else {
                        ++totalZero;
                        pseudo = pseudo + cost[thisBB];
                    }
                    --totalOne;
                    ++totalZero;
                    checked[thisBB] = true;
                }
            }
        }

        pseudo = pseudo + cost[i];

        if ( 1*totalOne > 0.8*totalOne + 1*totalZero + pseudo ) {
            bbType[i] = 0;
            ++zeros;
            for ( int j = 0; j < bbSize[i]; ++j ) {
                int thisGene = bb2gene[i][j];
                for ( int k = 0; k < geneSize[thisGene]; ++k ) {
                    int thisBB = gene2bb[thisGene][k];
                    if ( bbType[thisBB] ==1  ) {
                        ++cost[thisBB];
                    }
                }
            }
            if ( zeros == totalConflict ) {
                free(checked);
                break;
            }
        }
        free(checked);
    }

    free(cost);
    return zeros;
}


bool UFO::checkAndSetGlobalOptimum()
{
    if ( zeros == 0 ) {
        this->globalOptimum = m;
        for ( int i = 0; i < n; ++i ) {
            bestChromosome[i] = 1;
        }
        return true;
    }


    printf("checkAndSet : \n");
    int c[n];
    double maxFitness = 0;
    int maxC[n];
    for ( int i = 0; i < n; ++i ) {
        c[i]  = 0;
    }

    for ( int i = 0; i < pow(2.0,(double)n); ++i ) {
        int still = i;
        int index = 0;
        while (still > 0) {
            c[index] = still%2;
            still = still /2;
            ++index;
        }


        double fit = getFitness(c);
        if ( fit > maxFitness  ) {
            maxFitness = fit;
            for ( int i = 0; i < n; ++i )
                maxC[i] = c[i];
        }
    }
    
    printf( "global optimum : ");
    for ( int i = 0; i < n; ++i )
    {
        printf("%d", maxC[i]);
    }
    printf("\n");
    printf("   fitness =  %f\n", maxFitness);
    this->globalOptimum = maxFitness;
    for ( int i = 0; i < n; ++i )
        bestChromosome[i] = maxC[i];

    return true;
}

int UFO::randomConflict( const int& totalConflict)
{
    int total = totalConflict;
    if ( totalConflict == -1 )
        total = m/2;

    int* index = (int*) malloc(sizeof(int)*m);
    r.uniformArray(index,m,0,m-1);
    for ( int i = 0; i < total; ++i  ) {
        bbType[index[i]]=0;
    }
    free(index);
    zeros = total;
    return zeros;
}

int UFO::findConflictedBB()
{
    int totalConflicted = 0;
    int totalOne = 0;
    for ( int i = 0; i < m; ++i ) {
        totalOne = 0;
        for ( int j = 0; j < bbSize[i]; ++j ) {
            if ( bestChromosome[bb2gene[i][j]] == 1 )
                ++totalOne;
        }

        if (  kTrap9(bbSize[i], totalOne, bbType[i]) < 1 - 5e-6 )
            ++totalConflicted;

    }
    this->conflictedBB = totalConflicted;
    return totalConflicted;
}
