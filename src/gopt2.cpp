/*
 ============================================================================
 Name        : gopt2.c
 Author      : Leonid Tkachenko
 Version     :
 Copyright   : Your copyright notice
 Description : Hello OpenMP World in C
 ============================================================================
 */
#include <omp.h>
#include <time.h>
#include <stdio.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include <exception>
#include "../include/CPUgopt.h"
#include "../include/testFuncOwnImpl.h"
/**
 * Hello OpenMP World prints the number of threads and the current thread id
 */
int main (int argc, char *argv[]) {


	int rank = 13;

    double *box = new double[rank * 2];
    double *argmin = new double[rank * 2];
    double outMin;
    double eps = 0.001;
    int splitCoeff = 2;
    int status = -1;

    box[0] = -100;
    box[1] = 100;
    box[2] = -100;
    box[3] = 100;
    box[4] = -100;
    box[5] = 100;
    box[6] = -100;
    box[7] = 100;
    box[8] = -100;
    box[9] = 100;

    box[10] = -100;
    box[11] = 100;

    box[12] = -100;
    box[13] = 100;

    box[14] = -100;
    box[15] = 100;

    box[16] = -100;
    box[17] = 100;

    box[18] = -100;
    box[19] = 100;



    box[20] = -100;
    box[21] = 100;



    box[22] = -100;
    box[23] = 100;


    box[24] = -100;
    box[25] = 100;

/*
    box[26] = -100;
    box[27] = 100;



    box[28] = -100;
    box[29] = 100;
*/
    auto start = std::chrono::high_resolution_clock::now();

    calcOptValueOnCPU(box, 1, rank, fnCalcFunLimitsRozenbroke, eps, &outMin, &status, argmin);

    auto end = std::chrono::high_resolution_clock::now();


    std::cout << "\n";
    //std::cout << "min = " << outMin << "\t";
    printf("[%.8lf]\n",outMin);
    for(int i = 0; i < rank; i++)
    {
    	printf("[%.8lf; %8lf]", argmin[2*i],argmin[2*i+1]);
    }
    std::cout << "\n";
    std::cout << "time in millisecs: " << ((std::chrono::duration_cast<std::chrono::milliseconds>(end - start)).count())/1 << "\t";
    std::cout << "\n";


    delete [] box;



	//fnGetOptValueOnCPUSort(inBox, inRank, inNumBoxesSplitCoeff, inEps, inMaxIter, fnCalcFunLimitsMultiple2, outBox,&outMin, &outEps, &status);
	//fnGetOptValueOnCPUSort(inBox, inRank, inNumBoxesSplitCoeff, inEps, inMaxIter, fnCalcFunLimitsHypebolic2, outBox, &outMin, &outEps, &status);
    //fnGetOptValueOnCPUSort(inBox, inRank, inNumBoxesSplitCoeff, inEps, inMaxIter, fnCalcFunLimitsAluffiPentini2, outBox, &outMin, &outEps, &status);
    //fnGetOptValueOnCPUSort(inBox, inRank, inNumBoxesSplitCoeff, inEps, inMaxIter, fnCalcFunLimitsRozenbroke, outBox, &outMin, &outEps, &status);





 return 0;
}


