/*                                        *
 * U. Loeckmann                           *
 * bhint integrator                       *
 * for n-body problem.                    *
 *                                        */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <signal.h>
#include "bhi.h"
#include <time.h>
#include <omp.h>

#define SOFTENING_PAR_2  .0        // softening length squared (eps^2) for perturbing forces [TESTING]
                                   // NOTE: IN INTERNAL UNITS! For 3e-4 pc, choose:
                                   // (convert_length(3.e-4, 1) * convert_length(3.e-4, 1))
#define INIT_TIME      convert_time(50, 1)   // time until usual timesteps are reached

static unsigned long   stepsize_count[2*MAX_STEPSIZE_POWER];
static double stepsize_count_over[2*MAX_STEPSIZE_POWER];
static double t_maxval=1., t_over=.0;
static double t_eval=.0;
static double _m_max=.0;
static double softening_par2=.0;
int ul_kill = 0;
int _global_module=-1, _global_function=-1;
double lost_energy = .0;
double count_steps_over=.0, count_blocks_over=.0;
unsigned long count_steps=0, count_blocks=0;

enum _function {_UL_NBODY_OUTPUT_STEPSIZE=0, _UL_NBODY_GET_ENERGY, _UL_NBODY_GET_EKIN, _UL_NBODY_GET_ACC,
                _UL_NBODY_GET_EPOT, _UL_NBODY_MOVE_CENTER, _UL_NBODY_PREDICT_PART, _UL_NBODY_PREDICT,
                _UL_NBODY_GET_TIMESTEP_SIMPLE, _UL_NBODY_RE_INIT, _UL_NBODY_INTEGRATE, _UL_NBODY_SIGPROC,
                _UL_NBODY_MAIN};
/*
 * Return physical time corresponding to specified system time.
 */
double t_total(double t)
{
    return(convert_time(t + t_over, 0));
}

/*
 * Returns maximum star mass among particles 1..N-1
 */
double m_max()
{
    return _m_max;
}

void add_over(unsigned long val, unsigned long *count, double *count_over)
{
    if(!(*count + val > *count))
    {
        *count_over += *count;
        *count = 0;
    }
    *count += val;
}

void inc_stepsize(int n)
{
    add_over(1, stepsize_count+n, stepsize_count_over+n);
}

void output_stepsize()
{
    _enter_function(_UL_NBODY, _UL_NBODY_OUTPUT_STEPSIZE);
    int n;

    for(n = 0; n < 2*MAX_STEPSIZE_POWER; n++)
    {
        if(stepsize_count[n] + stepsize_count_over[n] > 0.1)
        {
            printf( "###tc%d### \t%e\t%e\t%f\n", n,
                    pow(2., n-MAX_STEPSIZE_POWER),
                    convert_time(pow(2., n-MAX_STEPSIZE_POWER), 0),
                    stepsize_count[n] + stepsize_count_over[n]);
        }
    }

    _exit_function();
}

double get_energy(struct particle parts[], int pcount, int pos)
{
    _enter_function(_UL_NBODY, _UL_NBODY_GET_ENERGY);
    struct particle *p=parts+pos, *pk;
    double e;

    e = scal_prod(p->v,  p->v) * p->m * .5;
    for(pk = parts; pk < parts + pcount; pk++)
    {
        if(p != pk)
        {
            e -= p->m * pk->m / sqrt(v_dist(p->x,  pk->x,  2) + (pk > parts ? softening_par2 : .0)) * ((pk >= parts + 1) ? .5 : 1.);
        }
    }

    #ifdef EXT_POT
    double phi=.0;
    add_force_extpot(p->x, NULL, NULL, NULL, &phi);
    e += p->m * phi;
    //printf("%d\t%e\n", pos, e);
    #endif // EXT_POT


    _exit_function();

    return e;
}

/*
 * Return total kinetic energy T of the system.
 */
double get_ekin(struct particle parts[], int pcount, int pred)
{
    _enter_function(_UL_NBODY, _UL_NBODY_GET_EKIN);
    struct particle *p;
    double e_kin = 0.0;

    for(p = parts; p < parts + pcount; p++)
    {
        e_kin += (pred ? scal_prod(p->vp, p->vp) : scal_prod(p->v,  p->v)) * p->m * .5;
    }

    _exit_function();

    return e_kin;
}


/*
 * Return acceleration for specified particle _pos_ as F/m.
 */
void get_acc(struct particle parts[], int pcount, int pos, double acc[3])
{
    _enter_function(_UL_NBODY, _UL_NBODY_GET_ACC);
    int i;
    double dist3;
    struct particle *p=parts+pos, *p2;

    acc[0] = 0.0;
    acc[1] = 0.0;
    acc[2] = 0.0;

    for(p2 = parts; p2 < parts + pcount; p2++)
    {
        if(p2 == p)
        {
            continue;
        }
        dist3 = v_dist(p2->x, p->x, 3);
        for(i = 0; i < DIMENSIONS; i++)
        {
            acc[i] -= p2->m / dist3 * (p->x[i] - p2->x[i]);
        }
    }

    _exit_function();
}

/*
 * Return total potential energy U of the system.
 */
double get_epot(struct particle parts[], int pcount, int pred)
{
    _enter_function(_UL_NBODY, _UL_NBODY_GET_EPOT);
    double e_pot = 0.0;
    struct particle *p1;

    struct particle *p2;
    for(p1 = parts; p1 < parts + pcount - 1; p1++)
    {
        for(p2 = p1 + 1; p2 < parts + pcount; p2++)
        {
            e_pot -= p1->m * p2->m  / (pred
                     ? sqrt(v_dist(p1->xp, p2->xp, 2) + (p1 > parts ? softening_par2 : 0))
                     : sqrt(v_dist(p1->x,  p2->x,  2) + (p1 > parts ? softening_par2 : 0)));
        }
    }

    #ifdef EXT_POT
    for(p1 = parts + 1; p1 < parts + pcount; p1++)
    {
        double phi=.0;
        add_force_extpot(pred ? p1->xp : p1->x, NULL, NULL, NULL, &phi);
        e_pot += p1->m * phi;
    }
    #endif // EXT_POT

    _exit_function();

    return e_pot;
}

/*
 * Move center particle, if required, to present position
 * and velocity resulting in CoM and total momentum equal 0.
 * If t is specified (>= 0), the position, velocity and time
 * of particle are updated. If t<0, only prediction values
 * xp, vp are set.
 */
void move_center(struct particle parts[], int pcount, double t)
{
    _enter_function(_UL_NBODY, _UL_NBODY_MOVE_CENTER);
    if(t >= 0)
    {
        parts[0].t = t;
    }
    _exit_function();
}

/*
 * Predict x and v (as xp and vp) for single particle
 * for given timestep dt from taylor series.
 */
void predict_part(struct particle *p, double dt, int perturb_only)
{
    _enter_function(_UL_NBODY, _UL_NBODY_PREDICT_PART);
    int i;
    double dt2=.0, dt3=.0, dt4=.0, dt5=.0;

    if(fabs(dt) >= DT_TOLERANCE)
    {
        dt2 = dt * dt  / 2.0;
        dt3 = dt * dt2 / 3.0;
        dt4 = dt * dt3 / 4.0;
        dt5 = dt * dt4 / 5.0;
    }

    for(i = 0; i < DIMENSIONS; i++)
    {
        if(fabs(dt) < DT_TOLERANCE)
        // nothing to predict
        {
            if(perturb_only)
            {
                return;
            }

            p->xp[i] = p->x[i];
            p->vp[i] = p->v[i];
        }
        else if(perturb_only)
        {
            // add perturbations only
            p->xp[i] += dt2 * p->a[i] + dt3 * p->a_[i];
            p->vp[i] += dt  * p->a[i] + dt2 * p->a_[i];
        }
        else
        {
            // predict new position
            p->xp[i] = p->x[i] + dt * p->v[i] + dt2 * p->a[i] + dt3 * p->a_[i];
            p->vp[i] = p->v[i] + dt * p->a[i] + dt2 * p->a_[i];
        }
        if(!p->active && fabs(dt) >= DT_TOLERANCE)
        {
            // additional terms if not disturbing
            p->xp[i] += dt4 * p->a_2[i] + dt5 * p->a_3[i];
            p->vp[i] += dt3 * p->a_2[i] + dt4 * p->a_3[i];
        }
    }

    _exit_function();
}


/*
 * Predict x and v (as xp and vp) for all particles
 * at given time t.
 */
void predict(struct particle parts[], int pcount, double t)
{
    _enter_function(_UL_NBODY, _UL_NBODY_PREDICT);
    struct particle *p;

    for(p = parts+1; p < parts + pcount; p++)
    {
        predict_part_hermite2(p, t);
    }

    _exit_function();
}

/*
 * Calculate simple timestep for hermite integration.
 */

double get_timestep_simple(double x0[3], double v0[3], double x[3], double v[3])
{
    _enter_function(_UL_NBODY, _UL_NBODY_GET_TIMESTEP_SIMPLE);
    int i;
    double x_[3], v_[3], v_x[3], vx;

    for(i = 0; i < 3; i++)
    {
        x_[i] = x[i] - x0[i];
        v_[i] = v[i] - v0[i];
    }

    vx = scal_prod(v_, x_) * 3 / scal_prod(x_, x_);

    for(i = 0; i < 3; i++)
    {
        v_x[i] = v_[i] + vx * x_[i];
    }

    _exit_function();

    return sqrt( scal_prod(x_, x_) / scal_prod(v_x, v_x) );
}

void re_init(struct particle *parts, int pcount)
{
    _enter_function(_UL_NBODY, _UL_NBODY_RE_INIT);
    struct particle *p;
    int i;

    fprintf(get_file(FILE_DEBUG),
            "# re-init at t=%e\t[%18.0lf steps, %18.0lf blocks]\n",
            t_total(parts[1].t),
            count_steps + count_steps_over,
            count_blocks + count_blocks_over);
    fflush(get_file(FILE_DEBUG));

    count_steps  = 0;
    count_blocks = 0;
    count_steps_over  = .0;
    count_blocks_over = .0;

    for(i = 0; i < 2*MAX_STEPSIZE_POWER; i++)
    {
        stepsize_count[i] = 0;
        stepsize_count_over[i] = .0;
    }

    for(p = parts; p < parts + pcount; p++)
    {
        p->io_tlast   = p->t;
        p->io_steps_c = 0;
        p->io_steps_p = 0;
    }

    _exit_function();
}


/*
 * Integrate problem of _pcount_ particles.
 * If _orbits_ > 0, integrate for _orbits_ orbits of 1st particle.
 * If _orbits_ < 0, integrate for -_orbits_ units of time.
 * Output _t_steps_ steps, or all if zero.
 */

void integrate( struct particle parts[], int pcount,
                int method, int print, double orbits, double t_steps,
                FILE *dumpfile, int print_dump_only)
{
    _enter_function(_UL_NBODY, _UL_NBODY_INTEGRATE);
    double t_max, pt;
    double e_, a_, t_, j_, _help[4];
    double rv_=0.0, rv=0.0, del_t = 0.0;
    double orig_a, orig_e, orig_j, orig_t, m_tot=0.0, t_old=-1.0;
    int i, j, orb_count=0, order, init=1, printed=1;
    struct particle *p;
    int force_print=0, init_phase=1, steps=0;
    time_t lastdump=.0;
    int skip = 0; //false

    for(i = 0; i < 2 * MAX_STEPSIZE_POWER; i++)
    {
        stepsize_count[i] = 0;
        stepsize_count_over[i] = .0;
    }

    //
    // Read from dump file and continue if it's necessary restart the integration.
    //
    if(dumpfile != NULL)
    {
        // MAX-TIME
        fscanf(dumpfile,"%*s    %le", &t_max);
        t_max = convert_time(t_max, 1);
        // INTEGRATE-VARS
        fscanf(dumpfile,
                "%*s    %le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%d\t%d\t%d\t%ld\t%le\t%ld\t%le\t%d",
                &t_maxval, &t_eval, &pt, &e_, &a_, &t_, &j_, &rv_, &rv, &del_t,
                &orig_a, &orig_e, &orig_j, &orig_t, &m_tot, &t_old, &t_over,
                &orb_count, &order, &init,
                &count_steps, &count_steps_over, &count_blocks, &count_blocks_over,
                &force_print);
        // INTEGRATE-PARAMS
        fscanf(dumpfile,
                "%*s    %d\t%d\t%d\t%le\t%le",
                &pcount, &method, &print, &orbits, &t_steps);

        parts = (struct particle *)malloc(pcount * sizeof(struct particle));
        // ENERGY
        fscanf(dumpfile, "%*s    %le", &lost_energy);

        //PARTICLES
        fscanf(dumpfile, "%*s%*c");
        assert(pcount == fread(parts, sizeof(struct particle), pcount, dumpfile));

        for(p = parts; p < parts + pcount; p++)
            if(p > parts && p->m > _m_max
                #ifdef P_IMBH
                && p->name != P_IMBH
                #endif
                )
                 _m_max = p->m;

        if(print_dump_only)
        {
            for(p = parts; p < parts + pcount; p++)
                printf("%1.12e\t %d \t%1.5e\t%-1.6e\t%-1.6e\t%-1.6e\t%-1.6e\t%-1.6e\t%-1.6e\t%1.4e\t%1.4e\t%1.4e\t%1.6e\t%1.6e\t%1.6e\t%1.6e\t%1.9e\t%1.9e\t%d\n",
                        t_total(p->t), p->name, convert_mass(p->m, 0),
                        convert_length(p->x[0], 0), convert_length(p->x[1], 0), convert_length(p->x[2], 0),
                        convert_length(convert_time(p->v[0], 1), 0), convert_length(convert_time(p->v[1], 1), 0), convert_length(convert_time(p->v[2], 1), 0),
                        convert_length(v_abs(p->x), 0), convert_length(p->curr_a, 0), 1 - p->curr_e, convert_length(p->rmin, 0), convert_length(p->rmax, 0),
                        convert_length(p->r_apo, 0), convert_length(p->r_peri, 0), p->energy, p->energy + p->m * (.5 * p->phi_stars), 1
                        );
            return;
        }

        lastdump = time(NULL);
        fprintf(get_file(FILE_OUTPUT), "######## Loaded %d particles at t=%le.\n", pcount, t_total(parts[1].t));
        fflush(get_file(FILE_OUTPUT));
        lastdump = time(NULL);

        skip = 1; // true
    }
    //
    // END Read from dump
    //

    if(!skip){
        // calculate start parameters
        softening_par2 = SOFTENING_PAR_2;
        pt = -get_epot(parts, pcount, 0) - get_ekin(parts, pcount, 0);

        for(j = 0; j < pcount; j++)
            m_tot += parts[j].m;

        get_reduced(parts, 1, &orig_e, &orig_a, &orig_t, &orig_j);//red_, red_+1, red_+2, red_+3);
        t_max = (orbits >= 0) ? (double)orbits * orig_t : -orbits;

        // not need to be parallel
        for(j = 1; j < pcount; j++)
        {
            evaluate_1_2(parts, pcount, j, 0, 0, parts[j].ha, parts[j].ha_);
            evaluate_1_2(parts, pcount, j, 1, pcount - 1, parts[j].a, parts[j].a_);
            parts[j].dt = normalize_dt(parts[j].t, get_timestep_central(parts, parts + j, MIN_EVALS) * .01);
        }

        check_app(parts, pcount, .0);

        // print headlines
        fprintf(get_file(FILE_OUTPUT), "#t       \tN\t \tdE/E0/dt\tdE/E0   \tE_lost  \tE_emit  \t \tdt_out  \tavg-steps\tavg-timestep\tavg-blockstep\tyr/s   \t \te_max   \ta_min   \tincl(1,2)\tangle(1,2)\t \te_P     \ta_P     \tr_P     \tE_P     \tdt_P\n");
        print_header(FILE_DETAIL, pcount);

        if(t_steps >= 0)
            output(print, parts, pcount, t_old, pt, parts[1].t, t_over, orbits, t_steps, t_max,
                    count_steps + count_steps_over, count_blocks + count_blocks_over, orb_count);
    }

    // Main loop
    while(1)
    {
        // Write dump file
         if(ul_kill > 0 || (printed &&
                    ((time(NULL) - lastdump > DUMP_INTERVAL)
                    || (parts[1].t + t_over >= t_max && !force_print))))
                {
                    dumpfile = write_dump_io();
                    fprintf(dumpfile,
                            "MAX-TIME %1.32le\n",
                            convert_time(t_max, 0));
                    fprintf(dumpfile,
                            "INTEGRATE-VARS %1.32le\t%1.32le\t%1.32le\t%1.32le\t%1.32le\t%1.32le\t%1.32le\t%1.32le\t%1.32le\t%1.32le\t%1.32le\t%1.32le\t%1.32le\t%1.32le\t%1.32le\t%1.32le\t%1.32le\t%16d\t%16d\t%16d\t%32ld\t%1.32le\t%32ld\t%1.32le\t%16d\n",
                            t_maxval, t_eval, pt, e_, a_, t_, j_, rv_, rv, del_t, orig_a, orig_e, orig_j, orig_t, m_tot, t_old, t_over,
                            orb_count, order, init, count_steps, count_steps_over, count_blocks, count_blocks_over, force_print);
                    fprintf(dumpfile,
                            "INTEGRATE-PARAMS %32d\t%32d\t%32d\t%1.32le\t%1.32le\n",
                            pcount, method, print, orbits, t_steps);
                    fprintf(dumpfile, "ENERGY %1.32le\n", lost_energy);
                    fprintf(dumpfile, "PARTICLES ");
                    i = fwrite(parts, sizeof(struct particle), pcount, dumpfile);
                    if(i != pcount)
                    {
                        fprintf(get_file(FILE_WARNING), "### PROBLEM WRITING DUMP FILE: WROTE %d OF %d PARTICLES\n", i, pcount);
                        fflush(get_file(FILE_WARNING));
                    }
                    fclose(dumpfile);
                    lastdump = time(NULL);

                    if(ul_kill > 0) break;
                }

        // End condition
        if(parts[1].t + t_over >= t_max && !force_print) break;
        printed = 0;
        //printf("%d\tcount_blocks:%ld \tcount_blocks_over: %f\n",tmp, count_blocks, count_blocks_over);

        // struct particle parts
        // int pcount
        // int ETA
        // int MIN_EVALS
        // double t_eval
        steps = step_hermite_2(parts, &pcount, ETA, MIN_EVALS, &t_eval);

        add_over(steps, &count_steps, &count_steps_over);

        // Variable reset
        #ifdef INIT_TIME
        if( (init_phase)        &&
            (steps == pcount-1) &&
            (parts[1].t + t_over >= INIT_TIME))
        {
            re_init(parts, pcount);
            init_phase = 0;
        }
        #endif



        if(t_steps < 0)
        {
            // count apocenter passage
            predict(parts, pcount, parts[1].t);
            rv_ = (parts[1].xp[0]-parts[0].xp[0])*(parts[1].vp[0]-parts[0].vp[0])
                    + (parts[1].xp[1]-parts[0].xp[1])*(parts[1].vp[1]-parts[0].vp[1])
                       + (parts[1].xp[2]-parts[0].xp[2])*(parts[1].vp[2]-parts[0].vp[2]);
            if(rv_ != 0 && rv != 0)
            {
                if(rv_ < 0 && rv > 0) // 1st parti cle at apocenter
                //if(rv_ > 0 && rv < 0) // 1s  t particle at pericenter
                {
                    if(t_steps < 0)
                    {
                              force_print = output(force_print ? -1 : print, parts, pcount, t_old, pt, parts[1].t,
                                                   t_over, orbits, t_steps, t_max, count_steps + count_steps_over,
                                                count_blocks + count_blocks_over, orb_count);
                        if(!force_print) printed = 1;
                    }
                    orb_count++;
                    del_t = (parts[1].t +t_over) / orb_count / orig_t;
                }
             }
            if(rv_ != 0) rv = rv_;
        }
        if(parts[1].t > t_maxval)
        {
            i = 1;
            for(p = parts; p < parts + pcount; p++)
                if(p->t < t_maxval)
                {
                    i = 0;
                    break;
                }
                if(i)
                {
                    for(p = parts; p < parts + pcount; p++)
                    {
                        p->t -= t_maxval;
                        p->io_tlast -= t_maxval;
                        p->htlast -= t_maxval;
                    }
                    t_over += t_maxval;
                    fprintf(get_file(FILE_WARNING),"### [%1.12e] reset time by %e (%e system units) ###\n",
                        t_total(.0), convert_time(t_maxval, 0), t_maxval);
                }
            }

        if(t_old < t_max && parts[1].t + t_over >= t_max)
            force_print = 1;
        if(t_steps >= 0 )// && parts[1].t > 0)
        {
            force_print = output(force_print ? -1 : print, parts, pcount, t_old, pt, parts[1].t, t_over, orbits, t_steps,
                            t_max, count_steps + count_steps_over, count_blocks + count_blocks_over, orb_count);
            if(!force_print) printed = 1;
        }
        // Counter
        add_over(1, &count_blocks, &count_blocks_over);
    }
    // END Main loop

    // output initial & end values
    get_reduced(parts, 1, &e_, &a_, &t_, &j_);
    get_approach_reduce_stats(_help);

    double epot = get_epot(parts, pcount, 0);
    double ekin = get_ekin(parts, pcount, 0);
    fprintf(get_file(FILE_OUTPUT), "# 1-e\t\t a\t\t T\t\t de/e\t\t da/a\t\t dT");
    fprintf(get_file(FILE_OUTPUT), "\t\tdj\t\t dE\t\t T_obs/T_H\t n\t orb\tangle\n");
    fprintf(get_file(FILE_OUTPUT), "# %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %d %3.2f\n",
            1-orig_e, convert_length(orig_a, 0), convert_time(orig_t, 0), e_/orig_e-1.0, a_/orig_a-1.0,
            t_/orig_t-1.0, j_/orig_j-1.0, pt/(epot + ekin) -1., del_t, orbits, orb_count,
            pcount > 2 ? vec_angle(parts[1].xp, parts[0].xp, parts[2].xp) : 0);
    fprintf(get_file(FILE_OUTPUT), "# t1=%e\tE0=%e\t dE=%e\t dE/E0=%1.12e\tlost=%e\n",
            t_total(parts[1].t), -pt, pt + epot + ekin, (pt + epot + ekin) / -pt, lost_energy);
    fprintf(get_file(FILE_OUTPUT), "######## %1.4e steps\t%1.4e blocks\t%1.4e approach checks\t%1.4e approach reduces\n",
            count_steps + count_steps_over, count_blocks + count_blocks_over,
            _help[0], _help[1]);

    output_stepsize();
    _exit_function();
}



void sigproc(int signal)
{
    _enter_function(_UL_NBODY, _UL_NBODY_SIGPROC);
    switch(signal)
    {
        //case SIGSEGV:
        case SIGUSR1:
            fprintf(get_file(FILE_WARNING), "# Received signal %d in function %d.%d - continue.\n",
                    signal,
                    #ifdef TRACK_FUNCTIONS
                    __my_lastmodule, __my_lastfunction);
                    #else
                    -2, -2);
                    #endif
            fflush(get_file(FILE_WARNING));
            ul_kill = 4*SIGNAL_PRINT_POS;

            break;
    case SIGUSR2:
        fprintf(get_file(FILE_WARNING), "# Received signal %d in function %d.%d - trying to stop.\n",
                signal,
                #ifdef TRACK_FUNCTIONS
                __my_lastmodule, __my_lastfunction);
                #else
                -2, -2);
                #endif
        fflush(get_file(FILE_WARNING));

        fprintf(get_file(FILE_WARNING), "#### Falling asleep... ####\n"); fflush(get_file(FILE_WARNING));
        raise(SIGSTOP);
        fprintf(get_file(FILE_WARNING), "#### ...woken up again. ####\n"); fflush(get_file(FILE_WARNING));

        break;

    default:
        fprintf(get_file(FILE_WARNING), "# Received signal %d in function %d.%d - shutting down.\n",
                signal,
                #ifdef TRACK_FUNCTIONS
                __my_lastmodule, __my_lastfunction);
                #else
                -2, -2);
                #endif
        fflush(get_file(FILE_WARNING));
        ul_kill = signal;
    }

    _exit_function();
}

//
// Main
//
int main(int argc, char **argv)

     /* PARAMETERS:
      *   char   method    [m] create Model
      *                    [i] Integrate with bhint
      *                    [r] Restart from dump file
      *                    [d] print Details from dump file
      *                    [c] Comment data file
      *                    [n] convert to N-body units
      *                    [o] calculate Orbit parameters
      *                    [u] print Units
      *   char  *infile    input file name
      *   double orbits    >0: number of orbits of 1st particle
      *                    <0: negative total time
      *   double t_steps   stepcount for output
      *                       output all if 0
      *                       every -t_steps'th orbit if <0
      *   char  *outfile   output file path and/or name
      */

{
    _enter_function(_UL_NBODY, _UL_NBODY_MAIN);
    double t_steps=0.0;
    struct particle *parts[1];
    char mode='s', *infile_name="", *outfile_name=(char*)malloc(300), method='u', *c;
    double orbits=10;
    int pcount=0, i;
    FILE *infile=NULL;
    double t_=t_maxval+DT_TOLERANCE*.5;

    signal(SIGINT, sigproc);
    signal(SIGQUIT, sigproc);
    signal(SIGTERM, sigproc);
    signal(SIGUSR1, sigproc);
    signal(SIGUSR2, sigproc);
    //signal(SIGSEGV, sigproc);

    while(t_ > t_maxval)
    {
        t_maxval *= 8.;
        t_ = t_maxval + DT_TOLERANCE * .5;
    }
    t_maxval *= .5;

    if(argc > 1) method = argv[1][0];
    switch(method)
    {
        // calculate orbital values
        case 'o':
            calc(argv+2); exit(0);
        case 'c':
            comment_datfile(); exit(0);
        case 'u':
            print_conv(); exit(0);
        // generate new cluster
        case 'm':
            printf("# nbody output \n# PARAMETERS: ");
            for(i = 0; i < argc; i++) printf("%s ", argv[i]);
            printf("\n");
            create(argc-2, argv+2); exit(0);

        // continue from dump
        case 'r':
        case 'd':
            infile = fopen(argv[2], "r");
            read_dump_io(infile, method == 'r' ? 0 : 1);
            mode = method;
            break;
        // convert to nbody units
        case 'n':
            infile = fopen(argv[2], "r");
            if(get_params(parts, &orbits, infile, &pcount, &_m_max) != 0)
            {
                fprintf(stderr, "No input parameters read - exiting!\n");
                exit(-1);
            }
            convert_to_nbody_units(*parts, pcount);
            exit(0);

        default:
            if(argc > 2) infile_name = argv[2];
            if(argc > 3) orbits = atof(argv[3]);
            if(argc > 4) t_steps=atof(argv[4]);

            if(strlen(infile_name)) infile = fopen(infile_name, "r");

            c = strrchr(infile_name, '/');
            if(c == NULL) c = infile_name;
            else c++;

            sprintf(outfile_name, "%s%s_%1.2e_%1d_%1.1e.bhint", argc > 5 ? argv[5] : "",
                    c, ETA, MIN_EVALS, orbits);
            init_files(outfile_name, 0);

            if(get_params(parts, &orbits, infile, &pcount, &_m_max) != 0)
            {
                fprintf(stderr, "No input parameters read - exiting!\n");
                exit(-1);
            }
        }

        fprintf(get_file(FILE_OUTPUT), "# nbody output \n# PARAMETERS: ");
        for(i = 0; i < argc; i++)
        {
            fprintf(get_file(FILE_OUTPUT), "%s ", argv[i]);
        }

        fprintf(get_file(FILE_OUTPUT), "\n");
        fprintf(get_file(FILE_OUTPUT), "# ETA=%e\tMIN_EVALS=%d\n",ETA, MIN_EVALS);
        fprintf(get_file(FILE_OUTPUT), "# DT_TOLERANCE=%e\tT_MAX=%e\tP_IMBH=%d\tGRAPE=%d\tSSE=%d\tEXT_POT=%d\n",
                DT_TOLERANCE, t_maxval, P_IMBH
                , 0
                , 0
                #ifdef EXT_POT
                , 1
                #else
                , 0
                #endif
        );

        #ifdef EXT_POT
        fprintf(get_file(FILE_OUTPUT), "# EXTERNAL POTENTIAL: RHO0=%e\tGAMMA1=%e\tGAMMA2=%e\tRMIN=%e\tR0=%e\tRMAX=%e\n",
                EP_RHO0, EP_GAMMA1, EP_GAMMA2, EP_RMIN, EP_R0, EP_RMAX);
        #else
        fprintf(get_file(FILE_OUTPUT), "# NO EXTERNAL POTENTIAL\n");
        #endif
        #ifdef N_MAX_DETAIL
        fprintf(get_file(FILE_OUTPUT), "# N_MAX_DETAIL=%d\n", N_MAX_DETAIL);
        #endif
        fprintf(get_file(FILE_OUTPUT),
                "# enclosed mass for Kepler parameter: M(<.01pc)=%e\tM(<.04pc)=%e\tM(<.1pc)=%e\tM(<.4pc)=%e\tM(<1.pc)=%e\n",
                convert_mass(M_ENCL_IO(convert_length(.01, 1)), 0),
                convert_mass(M_ENCL_IO(convert_length(.04, 1)), 0),
                convert_mass(M_ENCL_IO(convert_length(.1, 1)), 0),
                convert_mass(M_ENCL_IO(convert_length(.4, 1)), 0),
                convert_mass(M_ENCL_IO(convert_length(1., 1)), 0));

        switch(mode)
        {
            case 's':
                // Default single run
                // with the option 'i' to integrate teh data
                integrate(*parts, pcount, method, 1, orbits, t_steps, NULL, 0);
                break;

            case 'r':
                /* continue from dump */
                integrate(NULL, 0, 0, 0, 0, 0, infile, 0);
                break;

            case 'd':
                /* print dump file state */
                integrate(NULL, 0, 0, 0, 0, 0, infile, 1);
                break;
        }

        fprintf(get_file(FILE_DEBUG), "### eval time: %8.2fs\n", t_eval);
        close_files();

        _exit_function();

    return 0;
}
