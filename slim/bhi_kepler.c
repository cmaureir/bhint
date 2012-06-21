#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bhi.h"

enum _function {_UL_KEPLER_GET_COMMON=0, _UL_KEPLER_GET_JACOBI, _UL_KEPLER_UNJACOBI, _UL_KEPLER_GET_CONSTANTS,
                _UL_KEPLER_GET_REDUCED, _UL_KEPLER_STEP_KEPLER_1, _UL_KEPLER_STEP_KEPLER, _UL_KEPLER_SOLVE_KEPLER,
                _UL_KEPLER_KEPLER};
#define DEL_E          9.e-16  /* maximum error in angle E for kepler equation            */
#define DEL_E_HYP      2.e-15  /* maximum error in angle E for hyperbolic kepler equation */
#define MAX_ITER       50      /* maximum number of iterations to solve kepler            */

static double
    _1_17_16 = 1./17./16.,
    _1_16_15 = 1./16./15.,
    _1_15_14 = 1./15./14.,
    _1_14_13 = 1./14./13.,
    _1_13_12 = 1./13./12.,
    _1_12_11 = 1./12./11.,
    _1_11_10 = 1./11./10.,
    _1_10_9 = 1./10./9.,
    _1_9_8 = 1./9./8.,
    _1_8_7 = 1./8./7.,
    _1_7_6 = 1./7./6.,
    _1_6_5 = 1./6./5.,
    _1_4_3 = 1./4./3.,
    _1_3_2 = 1./3./2.;

/*
 * Get center of mass / total momentum
 * of predicted values for first _pcount_ particles into _com_.
 * Possible values for _val_:
 *    'x' for center of mass
 *    'v' for total momentum
 */
void get_common(struct particle parts[], int pcount, double com[3], char val)
{
    _enter_function(_UL_KEPLER, _UL_KEPLER_GET_COMMON);
    int i;
    struct particle *p;

    for(i = 0; i < 3; i++)
    {
        com[i] = .0;
        for(p = parts; p < parts + pcount; p++)
        {
            switch(val)
            {
                case 'x' : com[i] += p->m * p->xp[i]; break;
                case 'v' : com[i] += p->m * p->vp[i]; break;
            }
        }
        com[i] /= parts[pcount-1].eta;
    }

    _exit_function();
}


/*
 * Convert predicted values from cartesian to Jacobi coordinates.
 */
void get_jacobi(struct particle parts[], int pcount)
{
    _enter_function(_UL_KEPLER, _UL_KEPLER_GET_JACOBI);
    int i, j;
    double x_[3], v_[3];

    for(j = pcount - 1; j > 0; j--)
    {
        if(parts[j].is_jacobi)
            continue;
        get_common(parts, j, x_, 'x');
        get_common(parts, j, v_, 'v');
        for(i = 0; i < DIMENSIONS; i++)
        {
            parts[j].xp[i] -= x_[i];
            parts[j].vp[i] -= v_[i];
        }
        parts[j].is_jacobi = 1;
    }

    for(i = 0; i < DIMENSIONS; i++)
    {
        parts[0].xp[i] = .0;
        parts[0].vp[i] = .0;
    }
    parts[0].is_jacobi = 1;

    _exit_function();
}

/*
 * Convert predicted values from Jacobi
 * back to cartesian coordinates.
 */
void unjacobi(struct particle parts[], int pcount, int keep_cartesian)
{
    _enter_function(_UL_KEPLER, _UL_KEPLER_UNJACOBI);
    int i;
    struct particle *p;
    double sumx[3] = {.0, .0, .0}, sumv[3] = {.0, .0, .0};

    for(p = parts + pcount-1; p >= parts; p--)
    {
        assert(p->is_jacobi);
        for(i = 0; i < DIMENSIONS; i++)
        {
            p->xpc[i] = p->m_ / p->m * p->xp[i] - sumx[i] / p->eta;
            p->vpc[i] = p->m_ / p->m * p->vp[i] - sumv[i] / p->eta;
            if(keep_cartesian)
            {
                p->xp[i] = p->xpc[i];
                p->vp[i] = p->vpc[i];
                p->is_jacobi = 0;
            }
            sumx[i] += p->m * p->xpc[i];
            sumv[i] += p->m * p->vpc[i];
        }
    }

    _exit_function();
}

/*
 * Retrieve constants of motion (as if) in 2-body problem.
 * Returns j, e, a, omega in vectors/fields provided.
 */
void get_constants( double x[3], double v[3], double m_central,
                    double j[3], double e[3], double *a, double *omega)
{
    _enter_function(_UL_KEPLER, _UL_KEPLER_GET_CONSTANTS);
    int i;

    vec_prod(x, v, j);
    vec_prod (v, j, e);
    double r = sqrt(scal_prod(x, x));
    for(i = 0; i < DIMENSIONS; i++)
        e[i] = e[i]/m_central - x[i]/r;
    *a = scal_prod(j, j) / m_central / fabs(1 - scal_prod(e, e));
    *omega = sqrt(m_central / pow(*a, 3.0));

    _exit_function();
}

/*
 * Retrieve constants of motion (as if) in 2-body problem
 * (of first two particles) reduced to 1-body-problem.
 * Returns |j|, |e|, a, T in fields provided.
 */
void get_reduced(struct particle parts[], int pos, double *red_e, double *red_a, double *red_t, double *red_j)
{
    _enter_function(_UL_KEPLER, _UL_KEPLER_GET_REDUCED);
    int i;
    double r_[3], v_[3], j[3], e[3], a, omega;

    for(i = 0; i < 3; i++)
    {
        r_[i] = parts[pos].xp[i] - parts[0].xp[i];
        v_[i] = parts[pos].vp[i] - parts[0].vp[i];
    }
    get_constants(r_, v_, parts[0].m//+parts[pos].m
                    , j, e, &a, &omega);
    *red_e = v_abs(e);
    *red_a = a;
    *red_t = 2*M_PI/omega;
    *red_j = v_abs(j);

    _exit_function();
}

/*
 * Calculate Kepler position and velocity for given timestep _dt_
 * for particle no. _pos_.
 * _xp_ and _vp_ will be updated.
 */
void step_kepler_1(struct particle parts[], int pcount, int pos, double dt,
                   double *out_a, double *out_a_, double *out_a_2, double *out_a_3,
                   double *curr_a, double *curr_e)
{
    _enter_function(_UL_KEPLER, _UL_KEPLER_STEP_KEPLER_1);
    int i;
    struct particle *p0 = parts, *p1 = parts + pos;
    double r_[3], v_[3], j_[3], ecc_[3], a_[3], b_[3], _1_r2, afact, v_r_, v_v_, r_a_, v_a_, r_a__;
    double ecc, a, r, v, b, omega, e, mean, cosp, sinp;
    double m_c=p0->m, _cosp_ecc, e2, _1_ecc, _cosp_1, de_dt;//+p1->m;

    // get relative position / motion
    for(i = 0; i < 3; i++)
    {
        r_[i] = p1->xp[i] - p0->xp[i];
        v_[i] = p1->vp[i] - p0->vp[i];
    }

    // calculate ellipse constants
    get_constants(r_, v_, m_c, j_, ecc_, &a, &omega);
    //printf("#  [%d]:\t%e\t%e\t%e\n", pos, v_abs(ecc_), a, omega);

    ecc = v_abs(ecc_);
    // b_ = a * sqrt(|1-e²|) * (j_ x e_) / |j_ x e_|
    vec_prod(j_, ecc_, b_);
    b = a * sqrt(fabs(1-ecc*ecc)) / v_abs(b_);
    for(i = 0; i < 3; i++)
    {
        a_[i]  = a*ecc_[i]/ecc;            // semi major vector
        b_[i] *= b;                        // semi minor vector
    }

    if(curr_a != NULL) *curr_a = a;
    if(curr_e != NULL) *curr_e = ecc;

    if(ecc < 1)
    // elliptical orbit
    {
        if(!p1->is_elliptical)
        {
            fprintf(get_file(FILE_WARNING),
                    "#### [t=%1.12e] Particle #%d captured onto elliptical orbit with e=%e ####\n",
                    t_total(p1->t), pos, ecc);
                    p1->is_elliptical = 1;
        }
        // calculate eccentric anomaly e at t+dt
        e = (a - v_abs(r_)) / (ecc * a);
        if(e >= 1.0) e = .0;
        else if(e <= -1.0) e = M_PI;
        else e = acos(e);
        if(scal_prod(r_, b_) < 0)
            e = 2*M_PI - e;
        mean = (e - ecc*sin(e)) + dt * omega;
        while(mean >= 2. * M_PI)
            mean -= 2. * M_PI;

        e = solve_kepler(mean, ecc);

        cosp = cos(e);
        sinp = sin(e);
        _cosp_ecc = cosp - ecc;
        de_dt = omega / (1. - ecc * cosp);
        if(ecc > .99)
        {
            e2 = (e > 2. * M_PI - 1e-3) ? e - 2. * M_PI : e;
            if(e2 < 1e-3)
            {
                e2 *= e2;
                _1_ecc    = scal_prod(j_, j_)/(p0->m*a*(1+ecc));
                _cosp_1   =  - .5 * e2 * (1 - e2 / 12. * (1 - e2 / 30.));
                _cosp_ecc = _1_ecc + _cosp_1;
                de_dt     = omega / (_1_ecc - ecc * _cosp_1);
            }
        }
        for(i = 0; i < DIMENSIONS; i++)
        {
            r_[i] =   a_[i] * _cosp_ecc + b_[i] * sinp ;  // new location
            v_[i] = (-a_[i] * sinp      + b_[i] * cosp) * de_dt;   // direction of v only
        }
    }
    else
    // hyperbolic orbit  // parabolic?
    {
        if(p1->is_elliptical)
        {
            fprintf(get_file(FILE_WARNING), "#### [t=%1.12e+%1.12e] Particle #%d thrown onto hyperbolic orbit with e=%e (E=%e, a=%e) ####\n",
                    t_total(p1->t), convert_time(dt, 0), pos, ecc, p1->energy, convert_length(a, 0));
            p1->is_elliptical = 0;
        }
        if(ecc == 1)
            fprintf(get_file(FILE_WARNING), "# # # %e\tParabolic orbit of m%d treated as hyperbolic: e=%e\t(x=%e)\n",
                    t_total(p1->t), pos, ecc, convert_length(v_abs(p1->xp), 0));

        // calculate eccentric anomaly e at t+dt
        e = (a + v_abs(r_)) / (ecc * a);
        if(e < 1.0) e = .0;
        else if(scal_prod(r_, v_) < 0) e = -acosh(e);
        else e = acosh(e);

        e = kepler(ecc, ecc * sinh(e) - e + dt * omega);
        cosp = cosh(e);
        sinp = sinh(e);
        de_dt = omega / (ecc * cosp - 1.);
        for(i = 0; i < DIMENSIONS; i++)
        {
            r_[i] =   a_[i] * (ecc - cosp)  + b_[i] * sinp;  // new location
            v_[i] = (-a_[i] * sinp          + b_[i] * cosp) * de_dt;  // direction of v only
        }
    }

    // get |v_| from j_ = r_ x v_
    v = v_abs(v_);
    r = v_abs(r_);
    v = v_abs(j_) / (r * v * sin(acos(scal_prod(r_, v_)/ (r * v))));

    for(i = 0; i < DIMENSIONS; i++)
    {
        //v_[i] *= v;
        // total motion relative to fix central mass
        p1->xp[i] = p0->xp[i] + r_[i];
        p1->vp[i] = p0->vp[i] + v_[i];
    }

    if(out_a != NULL)
    {
        _1_r2 = 1. / scal_prod(r_, r_);
        afact = - m_c * _1_r2 * sqrt(_1_r2);
        //printf("4  %e %e %e\n", *(out_a), *(out_a+1), *(out_a+2));
        for(i = 0; i < DIMENSIONS; i++)
            out_a[i] = afact * r_[i];
            if(out_a_ != NULL)
            {
                v_r_ = scal_prod(v_, r_);
                for(i = 0; i < DIMENSIONS; i++)
                    out_a_[i] = afact * (v_[i] - 3 * _1_r2 * v_r_ * r_[i]);
                    if(out_a_2 != NULL)
                    {
                        v_v_ = scal_prod(v_, v_);
                        r_a_ = scal_prod(r_, out_a);
                        for(i = 0; i < DIMENSIONS; i++)
                            out_a_2[i] = afact * (out_a[i] - 3. * _1_r2 * (v_r_ * (2. * v_[i] - 5. * v_r_ * r_[i] * _1_r2)
                                         + (v_v_ + r_a_) * r_[i]));
                        if(out_a_3 != NULL)
                        {
                            v_a_  = scal_prod(v_, out_a);
                            r_a__  = scal_prod(r_, out_a_);
                            for(i = 0; i < DIMENSIONS; i++)
                                out_a_3[i] = afact * (out_a_[i]
                                            - 3. * _1_r2 * (3. * v_r_ * out_a[i]
                                            + 3. * (v_v_ + r_a_)
                                            * (v_[i] - 5. * v_r_ * _1_r2 * r_[i])
                                            + (3. * v_a_ + r_a__) * r_[i]
                                            + v_r_ * v_r_ * _1_r2
                                            * (-15. * v_[i] + 35. * v_r_ * _1_r2 * r_[i])));
                        }
                    }
            }
    }

    _exit_function();
}

/*
 * Calculate Kepler step.
 */
int step_kepler(struct particle parts[], int pcount)
{
    _enter_function(_UL_KEPLER, _UL_KEPLER_STEP_KEPLER);
    int i, step=0;
    double t0=-1.0;
    double tmin = -1.0;
    struct particle *p;

    // find next particle to move
    for(p = parts + 1; p < parts + pcount; p++)
        if(tmin < 0 || p->t + p->dt < tmin)
            tmin = p->t + p->dt;
        for(p = parts + 1; p < parts + pcount; p++)
            if((p->t + p->dt <= tmin + DT_TOLERANCE) && (p->t < t0 || t0 < 0))
                t0 = p->t;

    // move all particles to be moved
    for(p = parts + 1; p < parts + pcount; p++)
    {
        if(p->t + p->dt > tmin + DT_TOLERANCE || p->t > t0 + DT_TOLERANCE)
            continue;
        step = 1;
        for(i = 0; i < 3; i++)
        {
            p->xp[i] = p->x[i];
            p->vp[i] = p->v[i];
        }
        step_kepler_1(parts, pcount, p - parts, tmin-p->t, NULL, NULL, NULL, NULL, &(p->curr_a), &(p->curr_e));
        p->t = tmin;
        // calculate new timestep
        // timestep for kepler orbit
        p->dt = get_timestep_simple(parts[0].xp, parts[0].vp, p->x, p->v);
    }

    // make sure has moved at least 1 particle
    assert(step == 1);

    for(p = parts + 1; p < parts + pcount; p++)
        for(i = 0; i < 3; i++)
        {
            p->x[i] = p->xp[i];
            p->v[i] = p->vp[i];
        }
    parts[0].t = tmin;

    _exit_function();
    return pcount-1;
}

/*
 * Solve kepler equation for given
 * mean anomaly and eccentricity.
 * Return eccentric anomaly E as solution.
 */
double solve_kepler(double mean, double ecc)
{
    _enter_function(_UL_KEPLER, _UL_KEPLER_SOLVE_KEPLER);
     double e = ecc > .8 ? M_PI : mean, d;//, e2=e;//, de, e1;

    int n_iter = 0;

    for(d = e - ecc * sin(e) - mean; fabs(d) > DEL_E; d = e - ecc * sin(e) - mean)
    {
        if(n_iter++ > MAX_ITER)
        {
            fprintf(get_file(FILE_WARNING),
                "### Aborting elliptical kepler solution after %d iterations, keeping error of %e (e=%e, M=%e, E_=%1.12e, sin(E_)=%1.12e)\n",
                MAX_ITER, d,
                ecc, mean, e, sin(e));
            break;
        }
        e -= d / (1.0 - ecc * cos(e));
    }

    _exit_function();

    return e;
}

/*
 * Following function taken from
 * http://www.projectpluto.com/kepler.htm
 */
double kepler(const double ecc, double mean_anom)
{
    double curr, err;
    int is_negative = 0, n_iter = 0;

    if(!mean_anom)
        return(0.);
    _enter_function(_UL_KEPLER, _UL_KEPLER_KEPLER);

    if(ecc < .3)     /* low-eccentricity formula from Meeus,  p. 195 */
    {
        curr = atan2(sin(mean_anom), cos(mean_anom) - ecc);
        /* one correction step,  and we're done */
        err = curr - ecc * sin(curr) - mean_anom;
        curr -= err / (1. - ecc * cos(curr));
        _exit_function();

        return(curr);
    }

    if(mean_anom < 0.)
    {
        mean_anom = -mean_anom;
        is_negative = 1;
    }

    curr = mean_anom;
    if((ecc > .8 && mean_anom < M_PI / 3.) || ecc > 1.)    /* up to 60 degrees */
    {
        double trial = mean_anom / fabs(1. - ecc);

        if(trial * trial > 6. * fabs(1. - ecc))   /* cubic term is dominant */
        {
            if(mean_anom < M_PI)
                trial = pow(6. * mean_anom, 1.0/3.0);
            else        /* hyperbolic w/ 5th & higher-order terms predominant */
                trial = asinh(mean_anom / ecc);
        }
        curr = trial;
    }

    double thresh = DEL_E_HYP * mean_anom;
    if(ecc < 1.)
    {
        err = (curr - ecc * sin(curr)) - mean_anom;
        while(fabs(err) > thresh)
        {
            n_iter++;
            curr -= err / (1. - ecc * cos(curr));
            err = (curr - ecc * sin(curr)) - mean_anom;

            if(n_iter > MAX_ITER) // amended
            {
                fprintf(get_file(FILE_WARNING),
                    "#### ! Aborting kepler solution after %d iterations, keeping error of %e (e=%e, M=%e, E_=%e) ####\n",
                    MAX_ITER, err,
                    ecc, mean_anom, curr);
                break;
            }
        }
    }
    else
    {
        curr = log(2.*mean_anom/ecc+1.8); // taken from Burkardt & Danby, CeMec 31 (1983), 317-328
        double curr_abs = fabs(curr);
        err = (ecc * sinh(curr) - curr) - mean_anom;
        while(fabs(err) > thresh)
        {
            n_iter++;
            if(curr_abs < .72 && ecc < 1.1)
            {
                // [e * sinh(E) - E] / E << 1, and/or e * cosh(E) - 1 << 1
                // so don't calculate it directly
                double curr2 = curr * curr;

                // relative error when omitting nth order term needs to be smaller than resolution 1.e-15:
                // .5 * E^2 > 1e15 * E^n/n!, i.e. E < (n!/2e15)^(1/(n-2))
                // n = 16: E < .72, n = 10: E < .08

                if(curr_abs > .08)
                    curr -= err / ((ecc - 1) * cosh(curr) +
                            (((((((_1_16_15 * curr2 + 1.)
                            * _1_14_13 * curr2 + 1.)
                            * _1_12_11 * curr2 + 1.)
                            * _1_10_9 * curr2 + 1.)
                            * _1_8_7 * curr2 + 1.)
                            * _1_6_5 * curr2 + 1.)
                            * _1_4_3 * curr2 + 1.)
                            * .5 *curr2);
                else
                    curr -= err / ((ecc - 1) * cosh(curr) +
                            ((((_1_10_9 * curr2 + 1.)
                            * _1_8_7 * curr2 + 1.)
                            * _1_6_5 * curr2 + 1.)
                            * _1_4_3 * curr2 + 1.)
                            * .5 *curr2);
                curr2 = curr * curr;
                curr_abs = fabs(curr);

                if(curr_abs > .08)
                    err = ((ecc - 1) * sinh(curr) +
                            (((((((_1_17_16 * curr2 + 1.)
                            * _1_15_14 * curr2 + 1.)
                            * _1_13_12 * curr2 + 1.)
                            * _1_11_10 * curr2 + 1.)
                            * _1_9_8 * curr2 + 1.)
                            * _1_7_6 * curr2 + 1.)
                            * .05 * curr2 + 1.)
                            * _1_3_2 * curr2 * curr) - mean_anom;
                else
                    err = ((ecc - 1) * sinh(curr) +
                            ((((_1_11_10 * curr2 + 1.)
                            * _1_9_8 * curr2 + 1.)
                            * _1_7_6 * curr2 + 1.)
                            * .05 * curr2 + 1.)
                            * _1_3_2 * curr2 * curr) - mean_anom;
            }
            else
            {
                curr -= err / (ecc * cosh(curr) - 1.);
                err = (ecc * sinh(curr) - curr) - mean_anom;
            }
            if(n_iter > MAX_ITER) // amended
            {
                fprintf(get_file(FILE_WARNING),
                    "### Aborting hyperbolic kepler solution after %d iterations, keeping error of %e (e=%e, M=%e, E_=%1.12e, sinh(E_)=%1.12e)\n",
                    MAX_ITER, err,
                    ecc, mean_anom, curr, sinh(curr));

                break;
            }
        }
     }

    _exit_function();
    return(is_negative ? -curr : curr);
}
