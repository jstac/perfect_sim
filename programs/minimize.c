
#include <gsl/gsl_errno.h>  /* Defines GSL_SUCCESS, etc. */
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>

int minimize_convex(gsl_function *F,
        double a, double b, double *x_min, double tol)
{
    int status;
    double h = (b - a) * .0000001;   /* Used to test slope at boundaries */
    /* First deal with the special cases */
    if (b - a < tol)
    {
        *x_min = b;
        status = GSL_SUCCESS;
    }
    /* If the min is at a, then the derivative at a is >= 0.  Test for
     * this case. */
    else if (GSL_FN_EVAL(F, a + h) - GSL_FN_EVAL(F, a) >= 0)
    {
        *x_min = a;
        status = GSL_SUCCESS;
    }
    /* If the min is at b, then the derivative at b is >= 0.  Test for
     * this case. */
    else if (GSL_FN_EVAL(F, b - h) - GSL_FN_EVAL(F, b) >= 0)
    {
        *x_min = b;
        status = GSL_SUCCESS;
    }
    else
    {
        /* Choose x_guess so that it's value is less than either of the two
         * endpoint values. Since we've got this far, we know that at least
         * of of F(a + h) and F(b - h) has this property. */
        double x_guess;
        x_guess = (GSL_FN_EVAL(F, a + h) < GSL_FN_EVAL(F, b - h)) ? 
            a + h : b - h;
        int iter = 0, max_iter = 200;
        const gsl_min_fminimizer_type *T;
        gsl_min_fminimizer *s;
        T = gsl_min_fminimizer_goldensection;
        s = gsl_min_fminimizer_alloc(T);
        gsl_min_fminimizer_set(s, F, x_guess, a, b);

        do
        {
           iter++;
           status = gsl_min_fminimizer_iterate(s);  /* perform iteration */
           status = 
               gsl_min_test_interval(a, b, tol, 0.0); /* |a - b| < tol? */

           a = gsl_min_fminimizer_x_lower(s);
           b = gsl_min_fminimizer_x_upper(s);

           if (status == GSL_SUCCESS)
           {
               *x_min = gsl_min_fminimizer_x_minimum(s);  /* current est */
           }
        }
        while (status == GSL_CONTINUE && iter < max_iter);

        gsl_min_fminimizer_free(s);
    }
    return status;
}

/* An example of useage:

double f(double x, void *params)
{
    double *p = (double *) params;
    return cos(x) + *p;
}

double C = 3.0;

int main (void)
{
    double m = 2.0, result;
    double a = 0.0, b = 6.0;
    double epsilon = 0.001;
    int exit_val;

    gsl_function F;

    F.function = &f;
    F.params = &C;

    exit_val = minimize_convex(&F, a, b, m, &result, epsilon);
    printf("Minimizer: %g\n", result);
    printf("Function value: %g\n", f(result, &C));
    printf("%d\n", exit_val);

    return 0;
}

*/

