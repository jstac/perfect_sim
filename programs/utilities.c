
void first_diffs(double *x_diffs, double *x, int n) {
    int i;
    for (i = 0; i < n - 1; i++) 
    {
        x_diffs[i] = x[i+1] - x[i];
    }
}


/* A la NumPy linspace.  *ls is a pointer to the target array */
void linspace(double *ls, double lower, double upper, int n)
{
    double step = (upper - lower) / (n - 1);
    int i;
    for (i = 0; i < n; i++) 
    {
        ls[i] = lower;
        lower += step;
    }
}

/* Linear interpolation */
double lininterp(double *x_vals, double *y_vals, double *x_diffs, 
        double *y_diffs, int n, double x) 
{
    int i;
    if (x < x_vals[0]) 
    {
        return y_vals[0];
    }
    for (i = 0; i < n; i++)
    {
        if (x < x_vals[i+1]) 
        {
            return y_vals[i] + (x - x_vals[i]) * y_diffs[i] / x_diffs[i];
        }
    }
    /* x >= x[n-1], so */
    return y_vals[n-1];
}

/* Linear interpolation on an even grid, where x_step is the step
 * between points. */
double lininterp2(double *x_vals, double x_step, double *y_vals, 
        double *y_diffs, int n, double x) 
{
    int i;
    if (x <= x_vals[0])
    {
        return y_vals[0];
    }
    if (x >= x_vals[n-1])
    {
        return y_vals[n-1];
    }
    i = (int) ((x - x_vals[0]) / x_step);
    return y_vals[i] + (x - x_vals[i]) * y_diffs[i] / x_step;
}



