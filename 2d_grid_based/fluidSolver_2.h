// Source: https://www.dgp.toronto.edu/public_user/stam/reality/Research/pdf/GDC03.pdf
#ifndef FLUID2_H
#define FLUID2_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>


#define DIFFUSION  0.0f
#define VISCOSITY  0.0f
#define DT 1.0f

#define IX(i,j) ((i)+(N+2)*(j))
#define SWAP(x0,x) {float *tmp=x0; x0=x; x=tmp;}

// Pass a grid of points to the solver, the solver will update those points 
void solver_init(float* grid, int size);

// This will update the grid passed as a pointer at the init
void solver_step();

// This will delete the ressources acquired for the grid
void solver_delete();


// solver related
void project(int N, float *u, float *v, float *p, float *div);
void set_bnd(int N, int b, float *x);

// Implementations 
struct FluidGrid {
    int N;
    float dt;
    float diff;
    float visc;
    
    float *dens;
    float *dens_prev;
    
    float *u;
    float *v;

    float *u_prev;
    float *v_prev;
};
typedef struct FluidGrid FluidGrid;
FluidGrid* grid;
FluidGrid *FluidGridCreate(int size, float* grid_d);
void FluidGridFree(FluidGrid *grid);


FluidGrid *FluidGridCreate(int size, float* grid_d) {
    FluidGrid *grid = new FluidGrid(); 
    grid->N = size;
    int N = size;

    grid->dt = DT;
    grid->diff = DIFFUSION;
    grid->visc = VISCOSITY;

    grid->dens_prev = (float*) calloc((N+2) * (N+2), sizeof(float));                                 
    grid->dens = grid_d;
    
    grid->u = (float*) calloc((N+2) * (N+2), sizeof(float));
    grid->v = (float*) calloc((N+2) * (N+2), sizeof(float));
    
    grid->u_prev = (float*) calloc((N+2) * (N+2), sizeof(float));
    grid->v_prev = (float*) calloc((N+2) * (N+2), sizeof(float));
    
    return grid;
}

void FluidGridFree(FluidGrid *grid) {
    delete(grid->dens_prev);
    
    delete(grid->u);
    delete(grid->v);
    
    delete(grid->u_prev);
    delete(grid->u_prev);
    
    delete(grid);
}

float sum(float *array, int N) {
    float sum = 0;
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            sum += array[IX(i,j)];
        }
    }
    return sum;
}

// increase the array x by the product of dt and the array s
void add_source(int N, float *x, float *s, float dt) {
    for(int i = 0; i < (N+2)*(N+2); i++) x[i] += dt*s[i];
    // set_bnd(N,1,grid->u);
    // set_bnd(N,2,grid->v);
    // set_bnd(N,0,grid->dens);
}

// diffuse takes an array x0, simulate diffusion from it and put it into x
void diffuse(int N, int b, float *x, float *x0, float diff, float dt ) {
    float a = dt*diff*N*N; // diffusion rate, it is scaled by the number of cells in simulation

    // Simulate flow using averaging (20 iter) - Gauss-Seidel method
    for(int k = 0; k < 20; k++) {
        for(int i = 1; i <= N; i++) {
            for(int j = 1; j <= N; j++) {
                x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)] + x[IX(i+1,j)] + 
                                x[IX(i,j-1)] + x[IX(i,j+1)]))/(1+4*a);
            }
        }
        set_bnd(N,b,x);
    }
}

// The advection is implemented by getting previous possible position of a cell
// By using the negative of the velocities, then interpolate the new value from
// the values of the 4 closest cell (linear interpolation)
void advect(int N, int b, float *d, float *d0, float *u, float *v, float dt) {
    float dtN = dt*N; // scale dt by size grid
    
    for(int i = 1; i <= N; i++) {
        for(int j = 1; j <= N; j++) {
            // coords of previous cell (current cell - vel of the cell)
            float x =  i - dtN * u[IX(i,j)];
            float y =  j - dtN * v[IX(i,j)];

            // Bound values to stay in the grid
            if(x < 0.5) x = 0.5;
            if(x > N+0.5) x = N+0.5;
            if(y < 0.5) y = 0.5;
            if(y > N+0.5) y = N+0.5;
            
            // Get index of 4 surrounding cell 
            int i0 = (int) x;
            int i1 = i0+1;
            int j0 = (int) y;
            int j1 = j0+1;
            // Linear interpolation of the 4 surrounding cell
            float s1 = x - i0; float s0 = 1 - s1; float t1 = y - j0; float t0 = 1-t0;
            d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)])+s1*(t0*d0[i1,j0]+t1*d0[IX(i1,j0)]);
        }
    }
    set_bnd(N,b,d);
}

void dens_step(int N, float *x, float *x0, float *u, float *v, float diff, float dt) {
 //   add_source(N, x, x0, dt);
    SWAP(x0,x); 
    diffuse(N, 0, x, x0, diff, dt);
    SWAP(x0,x); 
    advect(N, 0, x, x0, u, v, dt);
}

void vel_step(int N, float *u, float *v, float *u0, float *v0, float visc, float dt) {
//    add_source(N, u, u0, dt);
//    add_source(N, v, v0, dt);
    
    SWAP(u0,u); 
    SWAP(v0,v); 
    diffuse(N, 1, u, u0, visc, dt);
    diffuse(N, 2, v, v0, visc, dt);
    
    project(N, u, v, u0, v0);
    SWAP(u0,u); 
    SWAP(v0,v); 
    advect(N, 1, u, u0, u0, v0, dt);
    advect(N, 2, v, v0, u0, v0, dt);
    project(N, u, v, u0, v0);
}

void project(int N, float *u, float *v, float *p, float *div) {
    float h = 1.0/N;

    // compute the amount of noncompressibility in each cell
    // The goal is to find if the velocities in the grid introduces compressibility, we do that by checking
    // if the continuity equatio is valid for each cell
    for(int i = 1; i <= N; i++) {
        for(int j = 1; j<= N; j++) {
            div[IX(i,j)] = -0.5*h*(u[IX(i+1,j)] - u[IX(i-1,j)] + v[IX(i,j+1)] - v[IX(i, j-1)]);
            p[IX(i,j)] = 0;
        }
    }
    set_bnd(N,0,div);
    set_bnd(N,0,p);

    // Gauss-seidel to solve the pressure at each cell
    for(int k = 0; k < 20; k++) {
        for(int i = 1; i<=N; i++) {
            for (int j = 0; j <= N; j++ ) {
                p[IX(i,j)] = (div[IX(i,j)] + p[IX(i-1,j)] + p[IX(i+1,j)] +  p[IX(i,j-1)] + p[IX(i,j+1)])/4;
            }
        }
        set_bnd(N,0,p);
    }
    // Substract the pressure gradient from each cell
    for(int i = 1; i<=N; i++) {
        for (int j = 1; j <= N; j++ ) {
            u[IX(i,j)] -= 0.5*(p[IX(i+1,j)] - p[IX(i-1,j)])/h;
            v[IX(i,j)] -= 0.5*(p[IX(i,j+1)] - p[IX(i,j-1)])/h;
        }
    }
    set_bnd(N,1,u);
    set_bnd(N,2,v);
}

void set_bnd(int N, int b, float *x)  {
    for (int i = 1; i<=N; i++) {
        x[IX(0  ,i)] = b == 1 ? -x[IX(1,i)] : x[IX(1,i)];
        x[IX(N+1,i)] = b == 1 ? -x[IX(N,i)] : x[IX(N,i)];
        x[IX(i,0  )] = b == 2 ? -x[IX(i,1)] : x[IX(i,1)];
        x[IX(i,N+1)] = b == 2 ? -x[IX(i,N)] : x[IX(i,N)];
    }
    x[IX(0  ,  0)] = 0.5*(x[IX(1,0  )] + x[IX(0  ,1)]);
    x[IX(0  ,N+1)] = 0.5*(x[IX(1,N+1)] + x[IX(0  ,N)]);
    x[IX(N+1,0  )] = 0.5*(x[IX(N,0  )] + x[IX(N+1,1)]);
    x[IX(N+1,N+1)] = 0.5*(x[IX(N,N+1)] + x[IX(N+1,N)]);
}




void solver_init(float* grid_d,int N) {
    grid = FluidGridCreate(N, grid_d);
}

void solver_step() {
    vel_step(grid->N, grid->u, grid->v, grid->u_prev, grid->v_prev, grid->visc, grid->dt);
    dens_step(grid->N, grid->dens, grid->dens_prev, grid->u, grid->v, grid->diff, grid->dt);
}
#endif