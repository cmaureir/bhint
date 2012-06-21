#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>
#include <time.h>
#include "bhi.h"

//#define SYNCHRONIZE_APPROACHING_TIMESTEPS // synchronize time-steps for approaching particles [TESTING]
#define MAX_APPROACH_FACTOR_2   8  // square of maximum approach factor per timestep
#define WARN_APPROACH_FACT      20 // multiple of impact parameter b
#define WARN_CLOSEENC              // warn for close encounters
//#define DEBUG_ALL

enum _function {_UL_TIMESTEP_CHECK_FAST_APPROACHES=0, _UL_TIMESTEP_NORMALIZE_DT};
unsigned long count_approach_reduce_t=0, count_approach_checks=0;
double count_approach_reduce_t_over=.0, count_approach_checks_over=0;
double _1_3 = 1. / 3., _1_27 = 1. / 27.;

/*
 * Return number of timestep reductions due to fast approaches.
 */
void get_approach_reduce_stats(double out[2])
{
    out[0] = count_approach_checks_over + count_approach_checks;
    out[1] = count_approach_reduce_t_over + count_approach_reduce_t;
}

/*
 * Check for fast approaches between p and pk and reduce p's timestep if necessary.
 * Return 1 if particles collided (merge!).
 */

int check_fast_approaches(  struct particle *parts,
                            struct particle *p, struct particle *pk/*,
                            double r_2*/)
{
    _enter_function(_UL_TIMESTEP, _UL_TIMESTEP_CHECK_FAST_APPROACHES);
    int i, collision=0;
    double temp, r_close_2=-1., r_temp_2, t_close, dt=.0, dt2=.0, dt3=.0, dt4=.0, dt5=.0;
    double *px=p->x, *pv=p->v, *pkx, *pkv, r_2;
    double x[3], v[3], a[3], a_[3];

    assert(p != pk);
    add_over(1, &count_approach_checks, &count_approach_checks_over);

    if(pk->active)
    {
        pkx = pk->x;
        pkv = pk->v;
    }

    else
    // have to use predicted x, v and corrected derivatives
    {
        pkx = pk->xp;
        pkv = pk->v;
        dt = p->t - pk->t;
        dt2 = .5 * dt * dt; dt3 = dt * dt2 * _1_3; dt4 = .25 * dt * dt3; dt5 = .2 * dt * dt4;
    }

  for(i = 0; i < 3; i++)
  {
      if(1)
      {
          x[i]  = px[i] - pkx[i];
      }
      else
      {
          x[i]  = px[i] - (pkx[i]
                  + dt * pk->v[i]
                  + dt2 * (pk->a[i] + pk->ha[i])
                  + dt3 * (pk->a_[i] + pk->ha_[i])
                  + dt4 * (pk->a_2[i] + pk->ha_2[i])
                  + dt5 * (pk->a_3[i] + pk->ha_3[i])
          );
      }

      if(pk->active)
      {
          v[i] = pv[i] - pkv[i];
          a[i]  = p->ha[i]  + p->a[i]  - pk->ha[i]  - pk->a[i];
          a_[i] = p->ha_[i] + p->a_[i] - pk->ha_[i] - pk->a_[i];
      }

      else
      {
          v[i] = pv[i] - (pkv[i]
                  + dt * (pk->a[i] + pk->ha[i])
                  + dt2 * (pk->a_[i] + pk->ha_[i])
                  + dt3 * (pk->a_2[i] + pk->ha_2[i])
                  + dt4 * (pk->a_3[i] + pk->ha_3[i])
                  );
          a[i]  = p->ha[i]  + p->a[i]  - (pk->a[i] + pk->ha[i]
                  + dt * (pk->a_[i] + pk->ha_[i])
                  + dt2 * (pk->a_2[i] + pk->ha_2[i])
                  + dt3 * (pk->a_3[i] + pk->ha_3[i])
                  );
          a_[i] = p->ha_[i] + p->a_[i] - (pk->a_[i] + pk->ha_[i]
                  + dt * (pk->a_2[i] + pk->ha_2[i])
                  + dt2* (pk->a_3[i] + pk->ha_3[i])
                  );
      }
    }

    r_2 = scal_prod(x, x);
    // linear approximation of time of closest encounter
    //t_close = -scal_prod(x, v) / scal_prod(v, v);
    // 2nd order approximation of time of closest encounter

    double xv, xa, v2, va, a2, _p, _q, _p3, _q2, _d, _u, _v, dy, t2, t3, _1_a2;
    xv = scal_prod(x, v);
    xa = scal_prod(x, a);
    v2 = scal_prod(v, v);
    va = scal_prod(v, a);
    a2 = scal_prod(a, a);
    _1_a2 = 1. / a2;
    dy = - va * _1_a2;
    _p = (a2 * 2.* (v2 + xa) - 3. * va * va) * (_1_a2 * _1_a2);
    _p3 = _p * _p * _p;
    _q = (2. * va * va * va - va * 2.*(v2 + xa) * a2 + 2. * xv * a2 * a2) * (_1_a2 * _1_a2 * _1_a2);
    _q2 = _q * _q;
    _d = 4. * _p3 + 27. * _q2;

    if(_d > 0)
    {
        _u = -.5 * _q;
        _v = sqrt(.25 * _q2 + _p3 * _1_27);
        t_close = cbrt(_u + _v) + cbrt(_u - _v) + dy;
    }
    else if(_d == 0)
    {
        t_close = cbrt(.5 * _q)  + dy;
        t3      = cbrt(-4. * _q) + dy;
        if(t3 > 0 && (t3 < t_close || t_close <= 0))
        {
            t_close = t3;
        }
    }

    else // _d < 0
    {
        _u = sqrt(-4. *_1_3 * _p);
        _v = acos(-.5 * _q * sqrt(-27. / _p3)) * _1_3;
        t_close =  _u *  _v               + dy;
        t2      = -_u * (_v + M_PI * _1_3) + dy;
        t3      = -_u * (_v - M_PI * _1_3) + dy;
        if(t2 > 0 && (t2 < t_close || t_close <= 0)) t_close = t2;
        if(t3 > 0 && (t3 < t_close || t_close <= 0)) t_close = t3;
    }
        while(1)
        {
            // check distance after next step
            r_temp_2 = .0;
            for(i = 0; i < DIMENSIONS; i++)
            {
                temp = x[i] + p->dt * (v[i] + p->dt * .5 * (a[i] /*+ p->dt / 3. * a_[i]*/));
                r_temp_2 += temp * temp;
            }

            if(r_2 > r_temp_2 * MAX_APPROACH_FACTOR_2 || r_2 * MAX_APPROACH_FACTOR_2 < r_temp_2)
            {
                add_over(1, &count_approach_reduce_t, &count_approach_reduce_t_over);
                #ifdef DEBUG_ALL
                fprintf(get_file(FILE_DEBUG),
                        "\t# halving dt: approach  m%d - m%d: \tr(t=%1.6e)=%1.2e\t\tr(t=%1.6e)=%1.2e\n",
                        pk->name, p->name,
                        t_total(p->t), convert_length(sqrt(r_2), 0),
                        t_total(p->t+p->dt), convert_length(sqrt(r_temp_2), 0));
                        fflush(get_file(FILE_DEBUG));
                #endif
                p->dt *= .5;

                #ifdef SYNCHRONIZE_APPROACHING_TIMESTEPS
                if(!pk->active)
                {
                    while(pk->htlast + .5 * pk->dt > p->t + p->dt + DT_TOLERANCE)
                    {
                        pk->dt *= .5;
                        #ifdef DEBUG_ALL
                        fprintf(get_file(FILE_DEBUG),
                                "#### [t=%1.12e] shrinking timestep for m%d as of close encounter with m%d to %e ####\n",
                                t_total(p->t),
                                pk->name,
                                p->name,
                                t_total(pk->dt));
                                fflush(get_file(FILE_DEBUG));
                        #endif
                    }
                }
                #endif

                continue;
                }

                if(r_2 < 9. * C_2G_C2 * C_2G_C2 * (pk->m + p->m) * (pk->m + p->m))
                {
                    // collision in 3 Schwarzschild-radii
                    collision = 1;
                    fprintf(get_file(FILE_WARNING), "#### [t=%1.12e] COLLISION of m%d and m%d at %1.12e: %e (r_S = %e) ####\n",
                            t_total(p->t),
                            p->name,
                            pk->name,
                            t_total(p->t + t_close),
                            convert_length(sqrt(r_2), 0),
                            convert_length(C_2G_C2 * (pk->m + p->m), 0));
                    fflush(get_file(FILE_WARNING));
                }


                if(t_close > .0 && t_close < p->dt)
                {
                    // close encounter will happen _during_ next step, now calculate distance
                    if(r_close_2 <.0)
                    {
                        r_close_2 = .0;
                        for(i = 0; i < DIMENSIONS; i++)
                        {
                            temp = (x[i] + t_close * (v[i] + .5 * t_close * a[i]));
                            r_close_2 += temp * temp;
                        }
                    }
                    if(r_close_2 < square(3. * C_2G_C2 * (pk->m + p->m)))
                    {
                        // collision in 3 Schwarzschild-radii
                        collision = 1;
                        fprintf(get_file(FILE_WARNING), "#### [t=%1.12e] COLLISION of m%d and m%d at %1.12e: %e (r_S = %e) ####\n",
                                t_total(p->t),
                                p->name,
                                pk->name,
                                t_total(p->t + t_close),
                                convert_length(sqrt(r_close_2), 0),
                                convert_length(C_2G_C2 * (pk->m + p->m), 0));
                        fflush(get_file(FILE_WARNING));
                    }

                    // approach to small multiple of impact parameter:
                    // r'_12 < warn_fact * b = warn_fact * 2 * r_1 * m / M

                    #ifdef WARN_CLOSEENC
                    if(r_close_2 * parts->m * parts->m < square(WARN_APPROACH_FACT * 2 * pk->m) * scal_prod(p->xp, p->xp)
                       && (N_MAX_DETAIL < -1 || p->name <= N_MAX_DETAIL || pk->name <= N_MAX_DETAIL)
                       )
                    {
                        fprintf(get_file(FILE_WARNING),
                                "\t# predicted close encounter m%d - m%d: \tr(t=%1.6e)=%1.2e\t\tr(t=%1.6e)=%1.2e=%1.2fb\tp.dt=%1.2e\tpk->dt=%1.2e (%1.2e el.) [%d:%d]\n",
                                pk->name, p->name,
                                t_total(p->t), convert_length(sqrt(r_2), 0),
                                t_total(p->t + t_close), convert_length(sqrt(r_close_2), 0),
                                sqrt(r_close_2) / (2. * v_abs(p->xp) * pk->m) * parts->m,
                                convert_time(p->dt, 0),
                                convert_time(pk->dt, 0),
                                convert_time(p->t - pk->t, 0),
                                pk->nearestneighbour, p->nearestneighbour);
                        fprintf(get_file(FILE_WARNING),
                                " PCE %1.12e\t%d\t%e\t%1.10e\t%1.10e\t%1.10e\t%1.10e\t%1.10e\t%1.10e\t%d\t%e\t%1.10e\t%1.10e\t%1.10e\t%1.10e\t%1.10e\t%1.10e\n",
                                t_total(p->t),
                                pk->name, convert_mass(pk->m, 0),
                                convert_length(pkx[0], 0), convert_length(pkx[1], 0), convert_length(pkx[2], 0),
                                convert_length(convert_time(pkv[0], 1), 0), convert_length(convert_time(pkv[1], 1), 0), convert_length(convert_time(pkv[2], 1), 0),
                                p->name, convert_mass(p->m, 0),
                                convert_length(px[0], 0), convert_length(px[1], 0), convert_length(px[2], 0),
                                convert_length(convert_time(pv[0], 1), 0), convert_length(convert_time(pv[1], 1), 0), convert_length(convert_time(pv[2], 1), 0)
                                );
                        fflush(get_file(FILE_WARNING));
                    }
                    #endif

                    if(r_2 > MAX_APPROACH_FACTOR_2 * r_close_2)
                    {
                        add_over(1, &count_approach_reduce_t, &count_approach_reduce_t_over);
                        #ifdef DEBUG_ALL
                        fprintf(get_file(FILE_DEBUG),
                                "\t# halving dt: encounter m%d - m%d: \tr(t=%1.6e)=%1.2e\t\tr(t=%1.6e)=%1.2e\tstep: t=%1.6e\n",
                                pk->name, p->name,
                                t_total(p->t), convert_length(sqrt(r_2), 0),
                                t_total(p->t + t_close), convert_length(sqrt(r_close_2), 0),
                                t_total(p->t + p->dt));
                                fflush(get_file(FILE_DEBUG));
                        #endif
                        p->dt *= .5;

                        #ifdef SYNCHRONIZE_APPROACHING_TIMESTEPS
                        if(!pk->active)
                            while(pk->htlast + .5 * pk->dt > p->t + p->dt + DT_TOLERANCE)
                            {
                                pk->dt *= .5;
                                #ifdef DEBUG_ALL
                                fprintf(get_file(FILE_DEBUG),
                                        "####  shrinking. timestep for m%d as of close encounter with m%d to %e ####\n",
                                        t_total(p->t),
                                        pk->name,
                                        p->name,
                                        t_total(pk->dt));
                                #endif
                            }
                        #endif
                        continue;
                        }
                }
            break;
        }

     _exit_function();
    return collision;
}

/*
 * Normalize timestep _dt_ such that _dt_ = 2^i
 * and _t_ / _dt_ = j for some i, j \elem Z.
 */

double normalize_dt(double t, double dt)
{
    int exp_;
    _enter_function(_UL_TIMESTEP, _UL_TIMESTEP_NORMALIZE_DT);
    //dt = pow(2., floor(log(dt) / log(2.) +.001));

    frexp(dt, &exp_);
    dt = ldexp(.5, exp_);

    while(fmod(t, dt) != .0)
    {
        //fprintf("### shrinking timestep to fit %e: %e -> %e\n", t, dt, dt/2.);
        dt *= .5;
    }
    assert(t + dt > t);
    _exit_function();

    return dt;
}
