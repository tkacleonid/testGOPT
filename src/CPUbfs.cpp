/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File: CPU_GO_FULL_SPLIT.cpp
 * Author: Leodid, Tkachenko
 *
 *
 */


#include "../include/CPUgopt.h"
#include <sys/mman.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <omp.h>
//#include "interval.h"




/**
*	Calculus minimum value for function on CPU
*	@param inbox pointer to Box
*	@param inRank number of variables
*	@param inNumBoxesSplitCoeff number of parts for each dimension
*	@param inEps required accuracy
*	@param inMaxIter maximum count of iterations
*	@param inFun pointer to optimazing function
*	@param outBox pointer to optimal box
*	@param outMin pointer to optimal value
*	@param outEps pointer to reached accuracy
*	@param outEps pointer to status of solving optimization problem
*/


#define FILEPATH "mmaped.bin"

void calcOptValueOnCPU(const double *_boxes, int _numBoxes, int _rank, void (*_fun)(const double *, int, double *), double _eps, double *_min, int *_status, double *_argmin)
{
	const int maxArrayLen = 5000000;

	double *boxes = new double[_rank*maxArrayLen*2];
	double *restBoxes = new double[_rank*maxArrayLen*2];
	double *funBounds = new double[3*maxArrayLen];
	const int splitCoeff = 2;

	const int nt = 8;
	const int size_b = maxArrayLen/ nt + 1;
	double *ar[nt];
	double *argmin[nt];

	for(int i= 0; i < nt; i++)
	{
		ar[i] = new double[size_b];
		argmin[i] = new double[size_b];
	}


	//copy Input Boxes in work set #1
	memcpy(restBoxes, _boxes, _rank*2*_numBoxes*sizeof(double));

	//calculate initial function record
	_fun(_boxes,_rank,funBounds);
	double funRecord = funBounds[2];
	double funLB;

	*_status = 1;

	int numBoxes = _numBoxes;

	int fd;


	fd = open(FILEPATH,O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);
	if(fd == -1)
	{
		perror("Error opening file for writing");
		exit(EXIT_FAILURE);

	}

	int result = lseek(fd, 10000000000,SEEK_SET);
	if(result == -1)
	{
		close(fd);
		perror("Error calling lseel");
		exit(EXIT_FAILURE);
	}

	result = write(fd,"",1);
	if(result != 1)
	{
		close(fd);
		perror("Error calling write");
		exit(EXIT_FAILURE);
	}




	int numBoxesInFile = 0;
	int s;
	double *map;
	off_t pa_offset,offset;

	while(true)
	{
	//	printf("\t\t-------count = %d\n\n",numBoxes);

		if(numBoxes*splitCoeff >= maxArrayLen)
		{
	//		printf("\nnumBoxesInFile = %d\n",numBoxesInFile);
			s = numBoxes/2;
			//printf("\nssssssssssssssssss = %d\n",numBoxesInFile*_rank*2*sizeof(double));
			offset = numBoxesInFile*_rank*2*sizeof(double);
			pa_offset = offset & ~(sysconf(_SC_PAGE_SIZE) - 1);
			map = (double *)mmap(0,s*_rank*2*sizeof(double),PROT_READ | PROT_WRITE, MAP_SHARED, fd, pa_offset);
			if(map == MAP_FAILED)
			{
				perror("MAP FAILED");
				exit(EXIT_FAILURE);
			}
			numBoxesInFile += s;
			memcpy(map,restBoxes+(numBoxes - s)*_rank*2,s*_rank*2*sizeof(double));
			if(munmap(map,s*_rank*2*sizeof(double)) == -1)
			{
				perror("Error un-mapping the file");
				exit(EXIT_FAILURE);
			}
			numBoxes -= s;

		}
		else if(numBoxes*splitCoeff <= maxArrayLen/4 && numBoxesInFile > 0)
		{
			s = maxArrayLen/4;
			if(numBoxesInFile <= s)  s = numBoxesInFile;

			offset = numBoxesInFile-s > 0? numBoxesInFile-s : 0;

			offset = offset*_rank*2*sizeof(double);
			pa_offset = offset & ~(sysconf(_SC_PAGE_SIZE) - 1);

	//		printf("\nssssssssssssssssss = %d\n",s);
			map = (double *)mmap(0,s*_rank*2*sizeof(double),PROT_READ | PROT_WRITE, MAP_SHARED, fd, pa_offset);
			if(map == MAP_FAILED)
			{
				perror("MAP FAILED");
				exit(EXIT_FAILURE);
			}
			numBoxesInFile -= s;
			memcpy(restBoxes+numBoxes*_rank*2,map,s*_rank*2*sizeof(double));
			if(munmap(map,s*_rank*2*sizeof(double)) == -1)
			{
				perror("Error un-mapping the file");
				exit(EXIT_FAILURE);
			}
			numBoxes += s;
		}

		//split boxes in work set #1
		int k = 0;
		int thid;




		#pragma omp parallel  for num_threads(nt)

		for(k = 0; k < numBoxes; k++)
		{

			//Find dimension with maximum length
			int maxDimensionIndex = 0;
			double maxDimension = restBoxes[(k*_rank)*2 + 1] - restBoxes[(k*_rank)*2];
			double h;
			for(int i = 0; i < _rank; i++)
			{
				h = (restBoxes[(k*_rank+i)*2 + 1] - restBoxes[(k*_rank+i)*2]);
				if (maxDimension < h)
				{
					maxDimension = h;
					maxDimensionIndex = i;
				}

			}
			h = maxDimension/splitCoeff;

			for(int n = 0; n < splitCoeff; n++)
			{
				for(int i = 0; i < _rank; i++)
				{
					if (i==maxDimensionIndex)
					{
						boxes[((k*splitCoeff + n)*_rank+i)*2] = restBoxes[(k*_rank+i)*2] + h*n;
						boxes[((k*splitCoeff + n)*_rank+i)*2 + 1] = restBoxes[(k*_rank+i)*2] + h*(n+1);
					} else
					{
						boxes[((k*splitCoeff + n)*_rank+i)*2] = restBoxes[(k*_rank+i)*2];
						boxes[((k*splitCoeff + n)*_rank+i)*2 + 1] = restBoxes[(k*_rank+i)*2 + 1];
					}

					//printf("[%f; %f]\t",boxes[((k*splitCoeff + n)*_rank+i)*2],boxes[((k*splitCoeff + n)*_rank+i)*2 + 1]);
				}

				_fun(&boxes[((k*splitCoeff + n)*_rank)*2],_rank,&funBounds[(k*splitCoeff + n)*3]);
				/*

				 if(funRecord > funBounds[(k*splitCoeff + n)*3 + 2] )
				{
					funRecord = funBounds[(k*splitCoeff + n)*3+2];
					memcpy(_argmin,&boxes[((k*splitCoeff + n)*_rank)*2],2*_rank*sizeof(double));
				}
				if( k == 0 && n == 0) funLB = funBounds[(k*splitCoeff + n)*3];
				else if (funLB > funBounds[(k*splitCoeff + n)*3]) funLB = funBounds[(k*splitCoeff + n)*3];

*/
				//printf("%f\t %f",funBounds[(k*splitCoeff + n)*3],funBounds[(k*splitCoeff + n)*3+2]);
				//printf("\n");
			}
		}

		funLB = funBounds[0];
		double funLBs[nt],funRecords[nt];

		for(int i = 0; i < nt; i++)
		{
			funLBs[i] = funLB;
			funRecords[i] = funRecord;
		}

#pragma omp parallel  for num_threads(nt)
		for(int i = 0; i < numBoxes*splitCoeff; i++)
		{
			thid = omp_get_thread_num();
			if(funRecords[thid] > funBounds[i*3 + 2] )
			{
				funRecords[thid] = funBounds[i*3+2];
				memcpy(argmin[thid],&boxes[i*_rank*2],2*_rank*sizeof(double));
			}
			if (funLBs[thid] > funBounds[i*3]) funLBs[thid] = funBounds[i*3];

			/*
			if(funRecord > funBounds[i*3 + 2] )
			{
					funRecord = funBounds[i*3+2];
					memcpy(_argmin,&boxes[i*_rank*2],2*_rank*sizeof(double));
			}
			if (funLB > funBounds[i*3]) funLB = funBounds[i*3];
			*/
		}

		for(int i = 0; i < nt; i++)
		{
			if(funRecord > funRecords[i] )
			{
					funRecord = funRecords[i];
					memcpy(_argmin,argmin[i],2*_rank*sizeof(double));
			}
			if (funLB > funLBs[i]) funLB = funLBs[i];
		}




/*
		funLB = funBounds[0];
		for(int i = 0; i < numBoxes*splitCoeff; i++)
		{
			if(funLB > funBounds[i*3])	funLB = funBounds[i*3];
		}
		*/

		//printf("\n\n-----%f\t %f-----\n\n",funLB,funRecord);


		double curEps = funRecord - funLB < 0 ? -(funRecord - funLB) : funRecord - funLB;

		if(curEps < _eps && numBoxesInFile == 0)
		{
			*_min = funRecord;
			*_status = 0;

			delete [] restBoxes;
			delete [] boxes;
			delete [] funBounds;

			return;
		}

		int cnt = 0;
		for(int i = 0; i < numBoxes*splitCoeff; i++)
		{
			if(funBounds[i*3] <= funRecord - _eps)
			{
				for(int j = 0; j < _rank; j++)
				{
					restBoxes[(cnt*_rank+j)*2] = boxes[(i*_rank+j)*2];
					restBoxes[(cnt*_rank+j)*2+1] = boxes[(i*_rank+j)*2+1];
				}
				cnt++;
			}
		}

		numBoxes = cnt;


		//if(numBoxes > 0) printf("\n\n%d: %f",numBoxes,funRecord);


		//countIter++;

	}
}







