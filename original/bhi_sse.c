#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "bhi.h"

#ifdef USE_SSE
void evolv1_(int *kw, double *mass, double *mt, double *r, double *lum, double *mc, double *rc, double *menv, double *renv, double *ospin,
             double *epoch, double *tm, double *tphys, double *tphysf, double *dtp, double *z, double zpars[20], double *dtm);
void zcnsts_(double *z, double zpars[20]);
double _ul_sse__zpars[20];

void sse_update_z(double z)
{
    static double z_old = -1.;
    if(z != z_old)
    {
        z_old = z;
        zcnsts_(&z, _ul_sse__zpars);
    }
}

double evolv1(int *kw, double *mass, double *mt, double *epoch, double *tphys, double *dtm, double *z, double *r)
{
    double lum, mc, rc, menv, renv, ospin, tms, dtp, tphysf=*tphys + *dtm;
    // kw mass mt epoch tphys [z]
    sse_update_z(*z);
    ospin  = .0;
    dtp    = .0;
    evolv1_(kw, mass, mt, r, &lum, &mc, &rc, &menv, &renv, &ospin, epoch, &tms, tphys, &tphysf, &dtp, z, _ul_sse__zpars, dtm);

    return *r;
}
#endif
