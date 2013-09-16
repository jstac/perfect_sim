#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "utilities.h"

#define GRIDSIZE 150
#define MAX_DEPTH 500

/* Parameters */
#define w 1.3712
#define R 1.0129  /* 1 + r */
/* The L shock has 3 states: {smin, 1, smax}, and is uniform.  We
 * assume that |smin -1| = |smax - 1| = d, so that E L = 1, and 
 * Var L = (2/3) * d**2.  We choose d so that Var L = sigma**2
 * = 0.4**2, which implies that d = 0.49 */
#define smin 0.51  /* 1 - d */
#define smax 1.49  /* 1 + d */
#define zb 0.952744  /* Value of zb corresponding to parameters */

/* Possible shock values */
double shock_vals[] = {smin, 1, smax};

/* The grid */
double gridmin = w * smin;
double gridmax = 14.0;
double grid[GRIDSIZE]; 
double gridstep;
/* Y values for interpolation (values of policy on grid) */
double vals[GRIDSIZE];  
/* First diffs used for linear interpolation */
double vals_diffs[GRIDSIZE-1];  

/* To store the shock path */
double shock_path[MAX_DEPTH];

/* The linear interpolation function */
double lp(double x) 
{
    return lininterp2(grid, gridstep, vals, vals_diffs, GRIDSIZE, x);
}

void initialize_pol_func(void)
{
    /* Read in the grid and values on the grid */
    FILE *f_grid;
    FILE *f_pol;
    int line_length = 100;
    char line[line_length];
    double temp;
    int i;

    f_grid = fopen("grid_dat.txt", "r");
    f_pol = fopen("pol_dat.txt", "r");

    i = 0;
    while (fgets(line, line_length, f_grid) != NULL)
    {
        sscanf(line, "%lf", &temp);
        grid[i] = temp;
        i++;
    }

    i = 0;
    while (fgets(line, line_length, f_pol) != NULL)
    {
        sscanf(line, "%lf", &temp);
        vals[i] = temp;
        i++;
    }

    fclose(f_grid);
    fclose(f_pol);

    /* Compute the diffs */
    first_diffs(vals_diffs, vals, GRIDSIZE);
    gridstep = grid[1] - grid[0];
}

/* For debugging */
void print_path(void)
{
    int i, n = 10000;
    double mean = 0.0, z = gridmax;
    srand(time(NULL));  
    for (i = 0; i < n; i++)
    {
        mean += z;
        z = w * shock_vals[rand() % 3] + R * lp(z);
        printf("%d %g\n", i, z);
    }
    mean = mean / n;
    /* printf("mean: %g\n", mean); */
}


/* If a perfect sample is generated then returns 0, 
 * else returns 1 */
int perf_sample(double *sample)
{
    double z;
    int i, regenerated, current_depth;

    /* Initialize the entire shock path */
    for (i = 0; i < MAX_DEPTH; i++)
    {
        shock_path[i] = shock_vals[rand() % 3];
    }

    current_depth = -1;
    regenerated = 0;
    while (regenerated == 0)
    {
        current_depth += 100;
        if (current_depth > MAX_DEPTH)
        {
            return 1;  /* Failed */
        }
        z = gridmax;
        for (i = current_depth; i >= 0; i--)
        {
            z = w * shock_path[i] + R * lp(z);
            if (z < zb && i > 0)
            {
                regenerated = 1;
            }
        }
    }

    *sample = z;
    return 0;
}


int main()
{
    int status;
    int i = 0, number_of_samples = 100000;
    double temp;
    double K = 0.0;

    srand(time(NULL));  

    initialize_pol_func();
    while (i < number_of_samples)
    {
        status = perf_sample(&temp);
        if (status == 0)
        {
            /* Uncomment the next line to print out the perfect samples
             * printf("%g\n", temp); */
            i++;
        }
    }
    return 0;
}
