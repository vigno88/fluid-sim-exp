#ifndef FLUID_H
#define FLUID_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "util.h"

#define DIFFUSION  0.1
#define VISCOSITY  0.1
#define DT 0.1

// Pass a grid of points to the solver, the solver will update those points 
void solver_init(float* grid, int size);

// This will update the grid passed as a pointer at the init
void solver_step();

// This will delete the ressources acquired for the grid
void solver_delete();



// Implementations 

struct FluidGrid {
    int size;
    float dt;
    float diff;
    float visc;
    
    float *s;
    float *density;
    
    float *Vx;
    float *Vy;

    float *Vx0;
    float *Vy0;
};
typedef struct FluidGrid FluidGrid;

FluidGrid* grid;

static void advect(int b, float *d, float *d0,  float *velocX, float *velocY, float dt, int N);
FluidGrid *FluidGridCreate(int size, float* grid_d);
void FluidGridFree(FluidGrid *grid);
void FluidGridAddDensity(FluidGrid* grid, int x, int y, float amount);
void FluidGridAddVelocity(FluidGrid *grid, int x, int y, float amountX, float amountY);
void FluidGridStep(FluidGrid *grid);
static void set_bnd(int b, float *x, int N);
static void lin_solve(int b, float *x, float *x0, float a, float c, int iter, int N);
static void diffuse (int b, float *x, float *x0, float diff, float dt, int iter, int N);
static void project(float *velocX, float *velocY, float *p, float *div, int iter, int N);
static void advect(int b, float *d, float *d0,  float *velocX, float *velocY, float dt, int N);

float sum(float *array, int size) {
    float sum = 0;
    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
            sum += array[index(j,i,size)];
        }
    }
    return sum;
}

void solver_init(float* grid_d,int size) {
    grid = FluidGridCreate(size, grid_d);

    for (int i = 5; i < 15; i++) {
        for (int j = 5; j < 15;j++) {
            FluidGridAddDensity(grid, i, j, 1);
        }
    }
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size;j++) {
            FluidGridAddVelocity(grid, j,i,2,0);
        }
    }
}

void solver_step() {
    FluidGridStep(grid);
  //  printf("sum grid: %f\n", sum(grid->density,grid->size));
}

FluidGrid *FluidGridCreate(int size, float* grid_d) {
    FluidGrid *grid = new FluidGrid(); 
    int N = size;
    
    grid->size = size;
    grid->dt = DT;
    grid->diff = DIFFUSION;
    grid->visc = VISCOSITY;
    
    grid->s = (float*) calloc(N * N, sizeof(float));
    grid->density = grid_d;
    
    grid->Vx = (float*) calloc(N * N, sizeof(float));
    grid->Vy = (float*) calloc(N * N, sizeof(float));
    
    grid->Vx0 = (float*) calloc(N * N, sizeof(float));
    grid->Vy0 = (float*) calloc(N * N, sizeof(float));
    
    return grid;
}

void FluidGridFree(FluidGrid *grid) {
    delete(grid->s);
    
    delete(grid->Vx);
    delete(grid->Vy);
    
    delete(grid->Vx0);
    delete(grid->Vy0);
    
    delete(grid);
}

void FluidGridAddDensity(FluidGrid* grid, int x, int y, float amount) {
    int N = grid->size;
    grid->density[index(x, y, N)] += amount;
}
void FluidGridAddVelocity(FluidGrid *grid, int x, int y, float amountX, float amountY) {
    int N =  grid->size;
    int index = index(x, y, N);
    
    grid->Vx[index] += amountX;
    grid->Vy[index] += amountY;
}

void FluidGridStep(FluidGrid *grid) {
    int N          = grid->size;
    float visc     = grid->visc;
    float diff     = grid->diff;
    float dt       = grid->dt;
    float *Vx      = grid->Vx;
    float *Vy      = grid->Vy;
    float *Vx0     = grid->Vx0;
    float *Vy0     = grid->Vy0;
    float *s       = grid->s;
    float *density = grid->density;
    
    diffuse(1, Vx0, Vx, visc, dt, 4, N);
    diffuse(2, Vy0, Vy, visc, dt, 4, N);
    
    project(Vx0, Vy0, Vx, Vy, 4, N);
    
    advect(1, Vx, Vx0, Vx0, Vy0, dt, N);
    advect(2, Vy, Vy0, Vx0, Vy0, dt, N);
    
    project(Vx, Vy, Vx0, Vy0, 4, N);
    
    diffuse(0, s, density, diff, dt, 4, N);
    advect(0, density, s, Vx, Vy, dt, N);
}
static void set_bnd(int b, float *x, int N) {
    for(int i = 1; i < N - 1; i++) {
        x[index(i, 0  , N)] = b == 2 ? -x[index(i, 1  , N)] : x[index(i, 1  , N)];
        x[index(i, N-1, N)] = b == 2 ? -x[index(i, N-2, N)] : x[index(i, N-2, N)];
    }

    for(int j = 1; j < N - 1; j++) {
        x[index(0  , j, N)] = b == 1 ? -x[index(1  , j, N)] : x[index(1  , j, N)];
        x[index(N-1, j, N)] = b == 1 ? -x[index(N-2, j, N)] : x[index(N-2, j, N)];
    }
    
    x[index(0, 0, N)]       = 0.33f * (x[index(1, 0, 0)] + x[index(0, 1, 0)]);
    x[index(0, N-1, N)]     = 0.33f * (x[index(1, N-1, 0)] + x[index(0, N-2, 0)]);
    x[index(N-1, 0, N)]     = 0.33f * (x[index(N-2, 0, 0)] + x[index(N-1, 1, 0)]);
    x[index(N-1, N-1, N)]   = 0.33f * (x[index(N-2, N-1, 0)] + x[index(N-1, N-2, 0)]);
}

static void lin_solve(int b, float *x, float *x0, float a, float c, int iter, int N) {
    float cRecip = 1.0 / c;
    for (int k = 0; k < iter; k++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                x[index(i, j, N)] =
                    (x0[index(i, j, N)]
                        + a*(    x[index(i+1, j  , N )]
                                +x[index(i-1, j  , N )]
                                +x[index(i  , j+1, N )]
                                +x[index(i  , j-1, N )]
                        )) * cRecip;
                }
            }
        set_bnd(b, x, N);
    }
}

static void diffuse (int b, float *x, float *x0, float diff, float dt, int iter, int N) {
    float a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a, iter, N);
}

static void project(float *velocX, float *velocY, float *p, float *div, int iter, int N) {
    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            div[index(i, j, N)] = -0.5f*(
                     velocX[index(i+1, j  , N  )]
                    -velocX[index(i-1, j  , N  )]
                    +velocY[index(i  , j+1, N  )]
                    -velocY[index(i  , j-1, N  )]
                )/N;
            p[index(i, j, N)] = 0;
        }
    }
    set_bnd(0, div, N); 
    set_bnd(0, p, N);
    lin_solve(0, p, div, 1, 6, iter, N);
    
    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            velocX[index(i, j, N)] -= 0.5f * (  p[index(i+1, j, N)] - p[index(i-1, j, N)]) * N;
            velocY[index(i, j, N)] -= 0.5f * (  p[index(i, j+1, N)] - p[index(i, j-1, N)]) * N;
        }
    }
    set_bnd(1, velocX, N);
    set_bnd(2, velocY, N);
}

static void advect(int b, float *d, float *d0,  float *velocX, float *velocY, float dt, int N){
    float i0, i1, j0, j1, k0, k1;
    
    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);
    float dtz = dt * (N - 2);
    
    float s0, s1, t0, t1;
    float tmp1, tmp2, x, y;
    
    float Nfloat = N;
    float ifloat, jfloat, kfloat;
    int i, j;
    
    for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
        for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
            tmp1 = dtx * velocX[index(i, j, N)];
            tmp2 = dty * velocY[index(i, j, N)];
            x    = ifloat - tmp1; 
            y    = jfloat - tmp2;
            
            if(x < 0.5f) x = 0.5f; 
            if(x > Nfloat + 0.5f) x = Nfloat + 0.5f; 
            i0 = floorf(x); 
            i1 = i0 + 1.0f;
            if(y < 0.5f) y = 0.5f; 
            if(y > Nfloat + 0.5f) y = Nfloat + 0.5f; 
            j0 = floorf(y);
            j1 = j0 + 1.0f; 
            
            s1 = x - i0; 
            s0 = 1.0f - s1; 
            t1 = y - j0; 
            t0 = 1.0f - t1;
            
            int i0i = i0;
            int i1i = i1;
            int j0i = j0;
            int j1i = j1;
            
            d[index(i, j, N)] = 
                s0 * ( t0 * d0[index(i0i, j0i, N)]
                    +( t1 * d0[index(i0i, j1i, N)]))
                +s1 * ( t0 * d0[index(i1i, j0i, N)]
                    +( t1 * d0[index(i1i, j1i, N)])); 
        }
    }
    set_bnd(b, d, N);
}

#endif