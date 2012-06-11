#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include "bhi.h"

double drand(double a, double b)
{
    return a + drand48() * (b - a);
}

void seedrand()
{
    struct timeval tmptv;
    static int init = 1;

    if(init)
    {
        init = 0;
        gettimeofday (&tmptv, (struct timezone *)NULL);
        long seed = tmptv.tv_usec + 1000000 * tmptv.tv_sec;
        srand48(seed);
        //printf("### +++ SEED: %ld +++ \t%e\t%e\n", seed, drand(.0,1.), drand(.0,1.));
    }
}

double gauss()
{
    static int iset = 0;
    static double gset;
    double fac, rsq, v1, v2;

    seedrand();

    if(iset == 0)
    {
        do
        {
            v1  = 2. * drand48() - 1.;
            v2  = 2. * drand48() - 1.;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1. || rsq == .0);

        fac  = sqrt(-2. * log(rsq) / rsq);
        gset = v1 * fac;
        iset = 1;

        return v2 * fac;
    }

    else
    {
        iset = 0;
        return gset;
    }
}

void setup_imf(double *m, double *alpha, int n, double *r, double *f)
{
    int i;
    r[0] = 1.;
    for(i = 0; i < n - 1; i++)
    {
        r[i+1] = r[i] * pow(m[i+1], alpha[i] - alpha[i+1]);
        printf("############################## %e\t%e\t%e\t%e\n", r[i], m[i], m[i+1], alpha[i]);
    }

    printf("############################## %e\t%e\t%e\t%e\n", r[i], m[i], m[i+1], alpha[i]);

    double x = .0;
    for(i = 0; i < n; i++)
    {
        x += r[i] * (pow(m[i+1], alpha[i]+1.) - pow(m[i], alpha[i]+1.)) / (alpha[i] + 1.);
    }

    for(i = 0; i < n; i++)
    {
        r[i] /= x;
    }

    x = .0;
    f[0] = .0;

    for(i = 0; i < n; i++)
    {
        f[i+1] = f[i] + r[i] * (pow(m[i+1], alpha[i]+1.) - pow(m[i], alpha[i]+1.)) / (alpha[i] + 1.);
    }
}

double get_imf_mass(double *m, double *alpha, int n, double *r, double *f)
{
    int i;
    double p = drand48();
    for(i = 0; f[i+1] <= p; i++);
    return pow((p - f[i]) * (alpha[i] + 1.) / r[i] + pow(m[i], alpha[i] + 1.), 1. / (alpha[i] + 1.));
}

//#define __UL_UTIL__GAUSS
#define __UL_UTIL__IMF
//#define __UL_UTIL__TEST
#ifdef __UL_UTIL__TEST

int main(int argc, char **argv)
{
    seedrand();
    #ifdef __UL_UTIL__GAUSS

    if(argc < 3)
    {
        fprintf(stderr, "Usage: gauss <a> <sigma> [n] \n");
        exit(-1);
    }

    double a = atof(argv[1]), sigma = atof(argv[2]);
    int n = 1, i;

    if(argc > 3)
    {
        n = atoi(argv[3]);
    }

    for(i = 0; i < n; i++)
    {
        printf("%3.8f\n", a + sigma * gauss());
    }

    #elif defined(__UL_UTIL__IMF)

    if(argc < 4 || argc % 2)
    {
        fprintf(stderr, "Usage: imf m_min m_max alpha0 [m_i alpha_i]* \n");
        exit(-1);
    }

    int i, j, n = argc / 2 - 1;
    double *m     = (double *)malloc((n+1) * sizeof(double));
    double *alpha = (double *)malloc(n     * sizeof(double));
    double *r     = (double *)malloc(n     * sizeof(double));
    double *f     = (double *)malloc((n+1) * sizeof(double));

    m[0] = atof(argv[1]);
    m[n] = atof(argv[2]);
    alpha[0] = -atof(argv[3]);

    for(i = 1; i < n; i++)
    {
        m[i]     =  atof(argv[2*i+2]);
        alpha[i] = -atof(argv[2*i+3]);
    }

    setup_imf(m, alpha, n, r, f);

    int d[] = {0, 0, 0, 0, 0, 0};
    /*
    for(i = 0; i < n; i++)
     printf("f_%d(m) = %1.8f * m ^ %1.1f, F(%1.3f) = %1.8f, F(%1.3f) = %1.8f\n",
     i, r[i], alpha[i], m[i], f[i], m[i+1], f[i+1]);
    */

    seedrand();
    for(j = 0; j < 1000000; j++)
    {
        double m_ = get_imf_mass(m, alpha, n, r, f);
        printf("%e\n", m_);
    }

    #endif
    exit(0);
}
#endif
