#pragma once

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <chrono>




void calcOptValueOnCPU(const double *_boxes, int _numBoxes, int _rank, void (*_fun)(const double *, int, double *), double _eps, double *_min, int *_status, double *_argmin);


