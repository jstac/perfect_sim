
#include <gsl/gsl_math.h>

/* gsl_function, lower, upper, guess, return value, tolerance (e.g. 0.001) */
double minimizeF(gsl_function *, double, double, double, double *, double);
