#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "utilities.h"
#include "minimize.h"

# define GRIDSIZE 250

/* Parameters */
#define w 1.3712
#define r 0.0129
#define beta 0.96
#define mu 2.0
#define sigma 0.4
/* The L shock has 3 states: {smin, 1, smax}, and is uniform.  We
 * assume that |smin -1| = |smax - 1| = d, so that E L = 1, and 
 * Var L = (2/3) * d**2.  We choose d so that Var L = sigma**2 */
#define d 0.49  /* sqrt(3.0 / 2.0) * sigma  */
#define smin 0.51  /* 1 - d */
#define smax 1.49  /* 1 + d */

/* Some constants derived from the parameters */
double m = 1 - mu;
double betad3 = beta / 3.0;
double R = 1 + r;
double shock_vals[] = {smin, 1, smax};
double gridmin = w * smin;
double gridmax = 14.0;

/* The grid */
double grid[GRIDSIZE]; 
/* First diffs used for linear interpolation */
double grid_diffs[GRIDSIZE-1]; 
/* The approximate optimal policy and successive iterates of the Bellman 
 * operator are represented as piecewise linear functions, with values
 * at the interpotation points stored in the following array */
double vals[GRIDSIZE];  
/* First diffs used for linear interpolation */
double vals_diffs[GRIDSIZE-1];  
double vals_temp[GRIDSIZE];  /* Used when recomputing vals */
/* And the function itself */
double lp(double x) 
{
    return lininterp(grid, vals, grid_diffs, vals_diffs, GRIDSIZE, x);
}

/* Utility function and its first derivative */
double U(double c) 
{
    return (pow(c + 1e-10, m) - 1) / m;
}

double Up(double c) 
{
    return pow(c + 1e-10, -mu);
}

/* The objective function for computing v-greedy policies.  The v 
 * used in the computation is that specified by the current value of
 * vals and vals_diffs.  The argument z is passed as a void pointer
 * in order to create a gsl_function.  The function is multiplied by
 * -1 because we are minimizing. */
double objective(double a, void *params) 
{
    double *z = (double *) params;
    return - U(*z - a) - betad3 * (lp(w * smin + R * a) +
                lp(w + R * a) + lp(w * smax + R * a));
}

void bellman(void) 
{
    int i;
    double result;
    double lower_bound = 0.0;
    double epsilon = 0.000001;
    gsl_function F;
    F.function = &objective;
    for (i = 0; i < GRIDSIZE; i++)
    {
        F.params = &grid[i];
        int exit_val = minimize_convex(&F, lower_bound, grid[i], &result, epsilon);
        vals_temp[i] = - objective(result, &grid[i]);  /* need to take negative */
        lower_bound = result;
    }
    for (i = 0; i < GRIDSIZE; i++)
    {
        vals[i] = vals_temp[i];
    }
}


int main(void)
{
    int i;
    /* Initialize grid and grid_diffs */
    linspace(grid, gridmin, gridmax, GRIDSIZE);
    first_diffs(grid_diffs, grid, GRIDSIZE);
    for (i = 0; i < GRIDSIZE; i++)
    {
        vals[i] = U(grid[i]);
    }
    first_diffs(vals_diffs, vals, GRIDSIZE);
    for (i= 0; i < 20; i++)
    {
        bellman();
    }
    /*
    for (i = 0; i < GRIDSIZE / 2; i++)
    {
        printf("%g %g\n", grid[i], objective(grid[i], &grid[GRIDSIZE / 2]));
    }
    */

    double lower_bound = 0.0;
    double result;
    double epsilon = 0.000001;
    gsl_function F;
    F.function = &objective;
    for (i = 0; i < GRIDSIZE; i++)
    {
        F.params = &grid[i];
        int exit_val = minimize_convex(&F, lower_bound, grid[i], &result, epsilon);
        vals_temp[i] = result;
        lower_bound = result;
    }
    for (i = 0; i < GRIDSIZE; i++)
    {
        vals[i] = vals_temp[i];
    }
    for (i = 0; i < GRIDSIZE; i++)
    {
        printf("%g %g\n", grid[i], w * smax + R * vals[i]);
    }
    return 0;
}
