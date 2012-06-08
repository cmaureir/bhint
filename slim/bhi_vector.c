#include "bhi.h"
/*
 * Return vector product of two vectors
 * in output vector _vp_.
 */

void vec_prod(double x[3], double y[3], double vp[3])
{
    vp[0] = x[1] * y[2] - x[2] * y[1];
    vp[1] = x[2] * y[0] - x[0] * y[2];
    vp[2] = x[0] * y[1] - x[1] * y[0];
    /*
        int i;
        for(i = 0; i < 3; i++)
            vp[i] = x[(i+1)%3] * y[(i+2)%3] - x[(i+2)%3] * y[(i+1)%3];
    */
}

/*
 * Return absolute distance between two vectors
 * to the power of _pot_.
 */

double v_dist(double x[3], double y[3], int pot)
{
    double dist;

    /*
    int i;
    double r[3];
    for(i = 0; i < DIMENSIONS; i++)
    r[i] = y[i] - x[i];
    dist = scal_prod(r, r);
    */
    dist = (y[0] - x[0]) * (y[0] - x[0]) + (y[1] - x[1]) * (y[1] - x[1]) + (y[2] - x[2]) * (y[2] - x[2]);

    switch(pot)
    {
        case 2:
            break;
        case 3:
            dist *= sqrt(dist);
            break;
        default:
            dist = sqrt(dist);
            if(pot != 1)
            {
                dist = pow(dist, pot);
            }
    }

    return dist;
}

/*
 * Return angle between three points in deg [0°,180°).
 */

double vec_angle(double o[3], double x[3], double y[3])
{
    int i;
    double sp=0.0, a[3], b[3];

    for(i = 0; i < 3; i++)
    {
        a[i] = x[i] - o[i];
        b[i] = y[i] - o[i];
        sp += a[i] * b[i];
    }

    return acos(sp / (v_abs(a) * v_abs(b))) * 180.0 /  M_PI;
}

/*
 * Return angle between three points in rad [0, 2pi)
 * with respect to normal vector j.
 */

double vec_angle_2pi(double o[3], double x[3], double y[3], double j[3])
{
    int i;
    double sp=0.0, a[3], b[3], c[3], ang;

    for(i = 0; i < 3; i++)
    {
        a[i] = x[i] - o[i];
        b[i] = y[i] - o[i];
        sp += a[i] * b[i];
    }

    if(sp == 0)
    {
        return ((a[0] * b[0] >= 0) && (a[1] * b[1] >= 0) && (a[2] * b[2] >= 0)) ? 0 : M_PI;
    }

    vec_prod(a, b, c);
    ang = acos(sp / (v_abs(a) * v_abs(b)));

    return ((c[0] * j[0] >= 0) && (c[1] * j[1] >= 0) && (c[2] * j[2] >= 0)) ? ang : 2*M_PI - ang;
}

/*
 * Return angle between between angular momentum of two bodies.
 */
double inclination(double r0[3], double v0[3], double r1[3], double v1[3], double r2[3], double v2[3])
{
    int i;
    double r1_[3], v1_[3], r2_[3], v2_[3], a[3], b[3];

    for(i = 0; i < 3; i++)
    {
        r1_[i] = r1[i] - r0[i];
        v1_[i] = v1[i] - v0[i];
        r2_[i] = r2[i] - r0[i];
        v2_[i] = v2[i] - v0[i];
    }

    vec_prod(r1_, v1_, a);
    vec_prod(r2_, v2_, b);

    return acos(scal_prod(a, b) / (v_abs(a) * v_abs(b))) * 180.0 /  M_PI;
}
