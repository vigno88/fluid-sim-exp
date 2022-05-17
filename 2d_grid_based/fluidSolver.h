#ifndef FLUID_H
#define FLUID_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "util.h"

#define DIFFUSION  0.5
#define VISCOSITY  0.5
#define DT 1

// Pass a grid of points to the solver, the solver will update those points 
void solver_init(float* grid, int size);

// This will update the grid passed as a pointer at the init
void solver_step();

// This will delete the ressources acquired for the grid
void solver_delete();

#endif