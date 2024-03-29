#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>
#include <time.h>
#include <signal.h>
#include "bhi.h"

#define CLOSE_MAX     100          // (initial) maximum number of neighbours; is increased dynamically
#define WARN_ENERGY_FACT .2
#define MAX_REDUCECOUNT 3

#define NO_USE_MASS_SUM

#define SWITCHON_PN  0
#define SWITCHOFF_PN 1

static struct timeval start, finish;
extern int ul_kill;

#define T_START \
gettimeofday(&start, NULL);

#define T_END \
gettimeofday(&finish, NULL);\
*t_eval += (double)(finish.tv_sec - start.tv_sec) + 1.e-6*(finish.tv_usec - start.tv_usec);

struct t_forceterm {
    double a[3], a_[3];
};
struct t_forceterm forceterm[N_MAX];

struct t_close {
    struct particle *p, *pk;
};

int close_n=0, close_max=CLOSE_MAX;
struct t_close *close_list = NULL;
int active[N_MAX], remove_part[N_MAX], movecount;

static int node_posmin=-1, node_posmax=-1;
static double _1_3 = 1. / 3., _1_6 = 1. / 6., _1_12 = 1. / 12., _1_14 = 1. / 14., _sqrt_mratio = .0;

static double step_size[2*MAX_STEPSIZE_POWER];
static int step_alloc[2*MAX_STEPSIZE_POWER], step_count[2*MAX_STEPSIZE_POWER],
           *step_part[2*MAX_STEPSIZE_POWER], step_min = MAX_STEPSIZE_POWER;
static double _1_LOG2 = -1.;
static int collisions, max_collisions=0;
static struct particle ***coll_vector=NULL;

#define PERTURBING_FORCE_RATIO .3  // at which ratio of central force is a particle perturber?

//
// add_close
//
void add_close(struct particle *p, struct particle *pk)
{
    //printf("%d\t%d\n", p->name, pk->name);
    if(close_list == NULL)
        close_list = (struct t_close *)malloc(CLOSE_MAX * sizeof(struct t_close));
    assert(p != pk);
    if(close_n >= close_max)
    {
        struct t_close *close_2 = (struct t_close *)realloc(close_list, 2*close_max * sizeof(struct t_close));
        if(close_2 == NULL)
        {
            fprintf(get_file(FILE_WARNING), "#### MEMORY ERROR (%d): Ignoring close encounter m%d and m%d ####\n",
                    close_max,
                    p->name, pk->name);
                    fflush(get_file(FILE_WARNING));
            return;
        }
        //printf("!!! %d\n", close_max);
        if(close_list != close_2)
            close_list = close_2;
       close_max *= 2;
    }
    close_list[close_n].p = p;
    close_list[close_n++].pk = pk;
}
// END
// add_close

void iteration_split(struct particle parts[],int pcount, int pos, int posmin, int posmax,
                    struct particle *p, double *px, double *pv,
                    double *a0, double *a1, double *a2, double *a_0,
                    double *a_1, double *a_2, double r_vic_2, double *phi,
                    double maxforce, double r_perturb_2, double r1_2,
                    double perturb, double rs_2, int *new_close_warn)
{
    double x_0, x_1, x_2, v_0, v_1, v_2;
    double min_r2 = 1.e99;
    double v_x_;
    double r_2 = .0;
    double _1_over_r2;
    int i;
    double afact[posmax-posmin];
    double afact_[posmax-posmin];
    double io_close[posmax-posmin];
    int close[posmax-posmin];
    int nearest[posmax-posmin];
    int energy[posmax-posmin];

    nearest[0] = p->nearestneighbour;
    for(i = posmin; i <= posmax; i++)
    {

        close[i-posmin] = 0;
        if(i == pos)
            continue;

        x_0 = parts[i].xp[0] - px[0];
        x_1 = parts[i].xp[1] - px[1];
        x_2 = parts[i].xp[2] - px[2];
        v_0 = parts[i].vp[0] - pv[0];
        v_1 = parts[i].vp[1] - pv[1];
        v_2 = parts[i].vp[2] - pv[2];

        r_2 = x_0*x_0 + x_1*x_1 + x_2*x_2;

        io_close[i-posmin] = p->io_close_warn;
        nearest[i-posmin] = p->nearestneighbour;

        if(r_2 < min_r2)
        {
            nearest[i-posmin] = i;
            min_r2 = r_2;
        }
        else if (i-posmin > 0)
        {
            if (i-1 == pos)
                nearest[i-posmin] = nearest[i-posmin-2];
            else
                nearest[i-posmin] = nearest[i-posmin-1];
        }

        _1_over_r2 = 1. / r_2;
        afact_[i-posmin] = parts[i].m * _1_over_r2;

        afact[i-posmin] = afact_[i-posmin] * sqrt(_1_over_r2);
        v_x_ = 3. * _1_over_r2 * (x_0*v_0 + x_1*v_1 + x_2*v_2);

        // TO DO: Divergence
        if(perturb)
            if(
                ((r_2 < r_perturb_2) &&
                (parts[i].m * r1_2 > PERTURBING_FORCE_RATIO * parts[0].m * r_2))
            )
            {
                if(io_close[i-posmin] <= 0 || io_close[i-posmin] > 16 * r_2)
                {
                    io_close[i-posmin] = r_2;
                }
                *new_close_warn = 1;
                energy[i-posmin] = get_energy(parts, pcount, pos);
            }

        if((r_2<rs_2)&&(r_2<9.*C_2G_C2*C_2G_C2*(p->m + parts[i].m)*(p->m + parts[i].m)))
        {
            // collision in 3 Schwarzschild-radii
            close[i-posmin] = 1;
        }
        else if(r_2 < r_vic_2 && v_x_ < 0)
        {
            // find approaching particles in vicinity
            close[i-posmin] = 1;
        }

        *a0  += afact[i-posmin] * x_0;
        *a_0 += afact[i-posmin] * (v_0 - v_x_ * x_0);
        *a1  += afact[i-posmin] * x_1;
        *a_1 += afact[i-posmin] * (v_1 - v_x_ * x_1);
        *a2  += afact[i-posmin] * x_2;
        *a_2 += afact[i-posmin] * (v_2 - v_x_ * x_2);
        *phi -= afact[i-posmin] * r_2;

    }

    for(i = posmin; i <= posmax; i++)
    {
        if(i == pos)
            continue;
        p->nearestneighbour = nearest[i-posmin];
        if(io_close[i-posmin] != 0)
            p->io_close_warn = io_close[i-posmin];
        if(energy[i-posmin] != 0 )
        {
            p->energy = energy[i-posmin];
        }
        if(close[i-posmin] == 1)
        {
            add_close(p, &parts[i]);
        }
        if(maxforce < .0 || afact_[i-posmin] > maxforce)
            maxforce = afact_[i-posmin];
    }

}



void iteration_normal(struct particle parts[],int pcount, int pos, int posmin, int posmax,
                    struct particle *p, double *px, double *pv,
                    double *a0, double *a1, double *a2, double *a_0,
                    double *a_1, double *a_2, double r_vic_2, double *phi,
                    double maxforce, double r_perturb_2, double r1_2,
                    double perturb, double rs_2, int *new_close_warn)
{
    double x_0, x_1, x_2, v_0, v_1, v_2;
    double min_r2 = 1.e99;
    double v_x_;
    double r_2, afact;
    double _1_over_r2;
    int i;

    for(i = posmin; i <= posmax; i++)
    {
        if(i == pos)
            continue;
        x_0 = parts[i].xp[0] - px[0]; v_0 = parts[i].vp[0] - pv[0];
        x_1 = parts[i].xp[1] - px[1]; v_1 = parts[i].vp[1] - pv[1];
        x_2 = parts[i].xp[2] - px[2]; v_2 = parts[i].vp[2] - pv[2];

        r_2 = x_0*x_0 + x_1*x_1 + x_2*x_2; //scal_prod(x_, x_);

        if(r_2 < min_r2)
        {
            p->nearestneighbour = i;
            min_r2 = r_2;
        }

        _1_over_r2 = 1. / r_2;
        afact = parts[i].m * _1_over_r2;
        if(maxforce < .0 || afact > maxforce)
            maxforce = afact;
        afact *= sqrt(_1_over_r2);
        v_x_ = 3. * _1_over_r2 * (x_0*v_0 + x_1*v_1 + x_2*v_2);

        if(perturb)
            if(
                ((r_2 < r_perturb_2) &&
                (parts[i].m * r1_2 > PERTURBING_FORCE_RATIO * parts[0].m * r_2))
            )
            {
                if(p->io_close_warn <= 0 || p->io_close_warn > 16 * r_2)
                {
                    //fprintf(get_file(FILE_DEBUG),
                    //        "#### [t=%1.12e] close encounter of m%d[%d] and m%d[%d]: %e (allowing %e, perturbing at %e) ####\n",
                    //        t_total(p->t),
                    //        p->name, pos,
                    //        parts[i].name, i,
                    //        convert_length(sqrt(r_2), 0),
                    //        convert_length(C_2G_C2 * (parts[i].m + p->m), 0),
                    //        convert_length(sqrt(parts[i].m / parts[0].m * scal_prod(px, px)), 0));
                    //fprintf(get_file(FILE_DEBUG),
                    //        " CE %1.12e\t%d\t%e\t%1.10e\t%1.10e\t%1.10e\t%1.10e\t%1.10e\t%1.10e\t%d\t%e\t%1.10e\t%1.10e\t%1.10e\t%1.10e\t%1.10e\t%1.10e\n",
                    //        t_total(p->t),
                    //        parts[i].name, convert_mass(parts[i].m, 0),
                    //        convert_length(parts[i].xp[0], 0), convert_length(parts[i].xp[1], 0), convert_length(parts[i].xp[2], 0),
                    //        convert_length(convert_time(parts[i].vp[0], 1), 0), convert_length(convert_time(parts[i].vp[1], 1), 0), convert_length(convert_time(parts[i].vp[2], 1), 0),
                    //        p->name, convert_mass(p->m, 0),
                    //        convert_length(px[0], 0), convert_length(px[1], 0), convert_length(px[2], 0),
                    //        convert_length(convert_time(pv[0], 1), 0), convert_length(convert_time(pv[1], 1), 0), convert_length(convert_time(pv[2], 1), 0)
                    //        );
                    //fflush(get_file(FILE_DEBUG));
                    p->io_close_warn = r_2;
                }
                *new_close_warn = 1;
                p->energy = get_energy(parts, pcount, p - parts);
            }
        if(
         (r_2 < rs_2) &&
         // calculate exact for star mass rather than maximum
         (r_2 < 9. * C_2G_C2 * C_2G_C2 * (p->m + parts[i].m) * (p->m + parts[i].m)))
        {
            // collision in 3 Schwarzschild-radii
            //fprintf(get_file(FILE_WARNING), "#### [t=%1.12e] COLLISION of m%d and m%d: %e (r_S = %e) ####\n",
            //        t_total(p->t),
            //        p->name,
            //        parts[i].name,
            //        convert_length(sqrt(r_2), 0),
            //        convert_length(C_2G_C2 * (parts[i].m + p->m), 0));
            //fflush(get_file(FILE_WARNING));
            add_close(p, &parts[i]);
        }

        // find approaching particles in vicinity
        else if(r_2 < r_vic_2 && v_x_ < 0)
            add_close(p, &parts[i]);

        *a0  += afact * x_0;
        *a_0 += afact * (v_0 - v_x_ * x_0);
        *a1  += afact * x_1;
        *a_1 += afact * (v_1 - v_x_ * x_1);
        *a2  += afact * x_2;
        *a_2 += afact * (v_2 - v_x_ * x_2);
        *phi -= afact * r_2;
    }

}




//
// register_dt
//
void register_dt(int pos, double dt)
{
    int n = floor(log(1.1 * dt) * _1_LOG2) + MAX_STEPSIZE_POWER;

    assert(n >= 0); assert(n < 2 * MAX_STEPSIZE_POWER);
    if(n < step_min)
        step_min = n;

    if(step_alloc[n] <= step_count[n])
    {
        step_alloc[n] = floor(1.4 * (float)step_alloc[n]);
        step_part[n] = (int*)realloc(step_part[n], step_alloc[n] * sizeof(int));
        assert(step_part[n] != NULL);
    }

    step_part[n][step_count[n]++] = pos;
    inc_stepsize(n);
}
// END
// register_dt

//
// init_dt
//
void init_dt(struct particle *parts, int pcount, int rebuild)
{
    static int init=0;
    if(init && !rebuild)
        return;

    int j, n;
    step_min = 2 * MAX_STEPSIZE_POWER;
    _1_LOG2 = 1. / log(2.);

    for(n = 0; n < 2 * MAX_STEPSIZE_POWER; n++)
    {
        if(!rebuild)
        {
            step_alloc[n] = 4;
            step_part[n]  = (int *)malloc(step_alloc[n] * sizeof(int));
            step_size[n] = pow(2., n - MAX_STEPSIZE_POWER);
        }
        step_count[n] = 0;
    }

    for(j = 1; j < pcount; j++)
        register_dt(j, parts[j].dt);

    fprintf(get_file(FILE_DEBUG), "### INITIALIZED TIMESTEPS:");
    for(n = 0; n < 2 * MAX_STEPSIZE_POWER; n++)
    if(step_count[n] > 0)
        fprintf(get_file(FILE_DEBUG), "\t%d:%d", n, step_count[n]);
    fprintf(get_file(FILE_DEBUG), "\n");
    fflush(get_file(FILE_DEBUG));

    init = 1;
}
// END
// init_dt

//
// predict_part_hermite2
//
// Hermite predict particles.
//
void predict_part_hermite2(struct particle *p, double t)
{
    int i;
    double dt, dt2, dt3, dt4, dt5;
    double dtn, dt2n, dt3n, dt4n, dt5n;

    dt = t - p->t;
    dt2 = dt * dt  * .5;
    dt3 = dt * dt2 * _1_3;
    dt4 = dt * dt3 * .25;
    dt5 = dt * dt4 * .2;
    dtn = t - p->htlast;
    dt2n = dtn * dtn  * .5; dt3n = dtn * dt2n * _1_3; dt4n = dtn * dt3n * .25; dt5n = dtn * dt4n * .2;

    //**/fprintf("## predicting #%d by %e for t=%1.12e\n", p->name, dt, tmin);      fflush(stdout);

   for(i = 0; i < DIMENSIONS; i++)
    {
        // orbital movement
        p->xp[i] = p->x[i] + (dt * p->v[i] ) + (dt2 * p->ha[i] ) + (dt3 * p->ha_[i]  )+ dt4 * p->ha_2[i] + (dt5 * p->ha_3[i]);
        p->vp[i] = p->v[i] + (dt * p->ha[i]) + (dt2 * p->ha_[i]) + (dt3 * p->ha_2[i] )+ dt4 * p->ha_3[i];
        // perturbing forces
        p->xp[i] += dt2n * p->a[i] + dt3n * p->a_[i] + dt4n * p->a_2[i] + dt5n * p->a_3[i];
        p->vp[i] += dtn  * p->a[i] + dt2n * p->a_[i] + dt3n * p->a_2[i] + dt4n * p->a_3[i];
    }

}
// END
// predict_part_hermite2


//
// add_collision
//
void add_collision(struct particle *p, struct particle *pk)
{
    int co;
    for(co = 0; co < collisions; co++)
        if(coll_vector[co][0] == p || coll_vector[co][0] == pk || coll_vector[co][1] == p || coll_vector[co][1] == pk)
        {
            return;
        }
    if(collisions + 1 > max_collisions)
    {
        max_collisions += 10;
        coll_vector = (struct particle ***)realloc(coll_vector, max_collisions*sizeof(struct particle **));
        for(co = max_collisions - 10; co < max_collisions; co++)
            coll_vector[co] = (struct particle **)malloc(2*sizeof(struct particle *));
    }
    if(p->m < pk->m)
    {
        coll_vector[collisions][0] = p;
        coll_vector[collisions][1] = pk;
    }
    else
    {
        coll_vector[collisions][0] = pk;
        coll_vector[collisions][1] = p;
    }
    collisions++;
}
// END
// add_collision

//
// check_app
//
void check_app(struct particle *parts, int pcount, double t)
{
    int j, collision;

    for(j = 0; j < close_n; j++)
    {
        if(close_list[j].p == NULL || close_list[j].p - parts >= pcount || close_list[j].pk - parts >= pcount )
            continue;
        collision = check_fast_approaches(parts, close_list[j].p, close_list[j].pk);
        if(collision && close_list[j].pk->active)
            add_collision(close_list[j].p, close_list[j].pk);

    }

}
// END
// check_app

//
// correct_timestep
//
void correct_timestep(struct particle *parts, int pcount, double tmin)
{
    int j;
    struct particle *p;

    for(j = 0; j < movecount; j++)
    {
        p = parts + active[j];
        if(p->dt < 2. * DT_TOLERANCE)
        {
            /* 	    if(dt >  2. * DT_TOLERANCE) */
            /* 	      fprintf(get_file(FILE_WARNING), "#### [t=%1.12e] [2] time step for particle #%d becoming too small: %e ####\n", */
            /* 		t_total(tmin), p->name, convert_time(p->dt, 0)); */
            p->dt = 2. * DT_TOLERANCE;
        }

        p->dt = normalize_dt(p->t, p->dt);

        register_dt(active[j], p->dt);
    }
}
// END
// correct_timestep

//
// sumforce
//
void sumforce(struct particle *parts)
{
    int i, j;
    struct particle *p;

    for(j = 0; j < movecount; j++)
    {
        p = parts + active[j];
        for(i = 0; i < DIMENSIONS; i++)
        {
            p->an[i]  = forceterm[j].a[i];
            p->a_n[i] = forceterm[j].a_[i];
        }
    }
}
// END
// sumforce

//
// find_move_particles
//
void find_move_particles(struct particle *parts, int pcount, double *tmin_out)
{
    struct particle *p;
    double tmin=-1.;
    //double help;
    //int pre_movecount=0,
    int j=0;

    p = parts + step_part[step_min][0];
    //printf("(%d-%d:%d)\n", step_min, step_count[step_min], step_part[step_min][0]);
    tmin = p->t + p->dt;
    for(movecount = 0; fmod(tmin, step_size[step_min]) == .0; step_min++)
    {
        for(j = 0; j < step_count[step_min]; j++)
        {
            active[movecount] = step_part[step_min][j];
            p = parts + active[movecount++];
            p->active = 1;
            p->dt = tmin - p->htlast;
            p->io_steps_p++;
            p->dtnext = 2.001 * p->dt;
        }
        step_count[step_min] = 0;
    }
    assert(movecount > 0);
    *tmin_out = tmin;
}
// END
// find_move_particles

// path_integral
//
void path_integral(int order, double dt, struct particle *p, double sign)
{

    double a[3]   = {p->a[0]   + p->ha[0]   + p->gr_a[0],   p->a[1]   + p->ha[1]   + p->gr_a[1],   p->a[2]   + p->ha[2]   + p->gr_a[2]};
    double a_[3]  = {p->a_[0]  + p->ha_[0]  + p->gr_a_[0],  p->a_[1]  + p->ha_[1]  + p->gr_a_[1],  p->a_[2]  + p->ha_[2]  + p->gr_a_[2]};
    double a_2[3] = {p->a_2[0] + p->ha_2[0] + p->gr_a_2[0], p->a_2[1] + p->ha_2[1] + p->gr_a_2[1], p->a_2[2] + p->ha_2[2] + p->gr_a_2[2]};

    switch(order)
    {
        case 0:
            break;
        case 1:
            //     trapezium rule              // F =     1/2 * (x1 - x0)   * (y1    + y0   )
            emit_energy(-.5 * dt * p->m * (scal_prod(p->v, p->gr_a)));
            break;
        case 3:
            //     PN_TRACK_3RD                // F =     1/2 * (x1 - x0)   * (y1    + y0   )
                                               //       -1/12 * (x1 - x0)^2 * (y1'   - y0'  )
            emit_energy(-.5 * dt * p->m * (scal_prod(p->v, p->gr_a) + sign * _1_6 * dt * (scal_prod(a, p->gr_a) + scal_prod(p->v, p->gr_a_))));
            break;
        case 5:
            //     PN_TRACK_5TH                // F =     1/2 * (x1 - x0)   * (y1    + y0   )
                                               //       -1/10 * (x1 - x0)^2 * (y1'   - y0'  )
                                               //      +1/120 * (x1 - x0)^3 * (y1(2) + y0(2))
            emit_energy(-.5 * dt * p->m * (scal_prod(p->v, p->gr_a) + sign * .2 * dt * (scal_prod(a, p->gr_a) + scal_prod(p->v, p->gr_a_)
                                 + sign * _1_12 * dt * (scal_prod(a_, p->gr_a) + 2. * scal_prod(a, p->gr_a_) + scal_prod(p->v, p->gr_a_2)))));
            break;
        case 7:
            //     PN_TRACK_7TH                // F =     1/2 * (x1 - x0)   * (y1    + y0   )
                                               //       -3/28 * (x1 - x0)^2 * (y1'   - y0'  )
                                               //       +1/84 * (x1 - x0)^3 * (y1(2) + y0(2))
                                               //     -1/1680 * (x1 - x0)^4 * (y1(3) - y0(3))
            emit_energy(-.5 * dt * p->m * (scal_prod(p->v, p->gr_a) + sign * _1_14 * dt * (3 * (scal_prod(a, p->gr_a) + scal_prod(p->v, p->gr_a_))
                            + sign * _1_3 * dt * (scal_prod(a_, p->gr_a) + 2. * scal_prod(a, p->gr_a_) + scal_prod(p->v, p->gr_a_2)
                            + sign * .05 * dt * (scal_prod(a_2, p->gr_a) + 3. * scal_prod(a_, p->gr_a_) + 3. * scal_prod(a, p->gr_a_2)
                            + scal_prod(p->v, p->gr_a_3))))));
    }
}
// END
// path_integral


//
// move_kepler
//
void move_kepler(struct particle *parts, int pcount, double tmin)
{
    struct particle *p;
    double *px, *pxp, *pv, *pvp;
    double dt;
    int j;
    printf("ITIME: %.10f\n",tmin);

    for(j = 0; j < movecount; j++)
        if(parts[active[j]].t <= tmin - DT_TOLERANCE)
        {
        //printf("Moving particle: %i\nMovecount: %d\nTime (t): %.10f\nTime(dt): %.10f\nTime (tmin): %.10f\n",
        //active[j],movecount, parts[active[j]].t, parts[active[j]].dt, tmin);
            p = parts + active[j];
            px = p->x;
            pv = p->v;

            pxp = p->xp;
            pvp = p->vp;

            pxp[0] = px[0];
            pxp[1] = px[1];
            pxp[2] = px[2];

            pvp[0] = pv[0];
            pvp[1] = pv[1];
            pvp[2] = pv[2];

            dt = tmin - p->t;

            std::cout << "Self (dt):"  << parts[active[j]].dt  << std::endl;
            std::cout << "Calculated (dt):"  << dt  << std::endl;

            for(;p->t < tmin; p->t += dt)
            {
                std::cout << "Loop time: " << p->t << std::endl;
                p->io_steps_c++;

                step_kepler_1(parts, pcount, p - parts, dt,
                              p->ha, p->ha_, p->ha_2, p->ha_3,
                              &(p->curr_a), &(p->curr_e));

            }
            std::cout << "End for" << std::endl;

            px[0] = pxp[0]; px[1] = pxp[1]; px[2] = pxp[2];
            pv[0] = pvp[0]; pv[1] = pvp[1]; pv[2] = pvp[2];
        }
}
// END
// move_kepler


//
// hermite_correct
//
void hermite_correct(struct particle *parts, int pcount)
{
    double dt=.0, dt2=.0, dt3=.0, dt4=.0, dt5=.0, _1_dt2=.0, _1_dt3=.0;
    int i, j;
    struct particle *p;

    for(j = 0; j < movecount; j++)
    {
        p = parts + active[j];

        // hermite corrector
        if(j == 0 || dt != p->dt)
        {
            dt     = p->dt;
            dt2    = dt * dt  * .5; dt3 = dt * dt2 * _1_3; dt4 = dt * dt3 * .25; dt5 = dt * dt4 * .2;
            _1_dt3 = 1. / dt3; _1_dt2 = _1_dt3 * dt * _1_3;
        }
        for(i = 0; i < DIMENSIONS; i++)
        {
            p->a_2[i] = -(3.*(p->a[i]-p->an[i]) + p->dt*(2.*p->a_[i]+p->a_n[i])) * _1_dt2;
            p->a_3[i] =  (2.*(p->a[i]-p->an[i]) + p->dt*(p->a_[i]+p->a_n[i])) * _1_dt3;
            p->x[i]  += dt4 * p->a_2[i] + dt5 * p->a_3[i];
            p->v[i]  += dt3 * p->a_2[i] + dt4 * p->a_3[i];
            p->a[i]   = p->an[i];
            p->a_[i]  = p->a_n[i];
        }
    }
}
// END
// hermite_correct


//
// get_timestep_aarseth
//
double get_timestep_aarseth(double a[3], double a_[3], double a_2[3], double a_3[3], double eta)
{
    double _a_ = v_abs(a_), _a_2 = v_abs(a_2);
    return sqrt(eta * ((v_abs(a)*_a_2 + _a_*_a_) / (_a_*v_abs(a_3) + _a_2*_a_2)));
}


//
// get_timestep_central
//
double get_timestep_central(struct particle *parts, struct particle *p, double min_evals)
{
    double r2 = scal_prod(p->x, p->x);
    return sqrt(r2 * sqrt(r2) / parts[0].m) * (2. * M_PI / min_evals);
}

//
// find_timesteps
//
void find_timesteps(struct particle *parts, int pcount, double tmin, double eta, double min_evals)
{
    double dt, dtn;
    struct particle *p;
    int j, exp_;

    for(j = 0; j < movecount; j++)
    {
        p = parts + active[j];

        dt = get_timestep_aarseth(p->a, p->a_, p->a_2, p->a_3, eta);
        // don't let perturbation timestep get bigger than ETA_FACT * central timestep

        //dtn = get_timestep_aarseth(p->ha, p->ha_, p->ha_2, p->ha_3, eta) * eta_fact;
        dtn = get_timestep_central(parts, p, min_evals);

        if(dtn < dt)
            dt = dtn;

        // allow timestep to double at maximum
        if(p->dtnext < dt)
            dt = p->dtnext;

        //p->t = tmin;
        p->htlast = tmin;

        // set to power of 2
        frexp(dt, &exp_);
        p->dt = ldexp(.5, exp_);

    }
}
// END
// find_timesteps


//
// herm_pred
//
void herm_pred(struct particle *parts, int pcount, double tmin, int pred_only)
{
    double dt, dt2, dt3, dt4, dt5;
    int i;
    struct particle *p;

    for(p = parts + 1; p < parts + pcount; p++)
    {
        if(!p->active && (p->t <= tmin - DT_TOLERANCE))
        {
            dt = tmin - p->t;
            dt2 = .5 * dt * dt;
            dt3 = dt * dt2 * _1_3;
            dt4 = .25 * dt * dt3;
            dt5 = .2 * dt * dt4;
            for(i = 0; i < DIMENSIONS; i++)
                p->xp[i] = p->x[i]+ dt * p->v[i] + dt2 * (p->a[i] + p->ha[i]) + dt3 * p->ha_[i] + dt4 * p->ha_2[i]+ dt5 * p->ha_3[i];
        }
        else
        {
            if(p->active)
            {
                dt = p->dt;
                dt2 = dt * dt  * .5;
                dt3 = dt * dt2 * _1_3;
                dt4 = dt * dt3 * .25;
                dt5 = dt * dt4 * .2;
                for(i = 0; i < DIMENSIONS; i++)
                {
                    p->x[i] += dt2 * p->a[i] + dt3 * p->a_[i];// + dt4 * p->a_2[i] + dt5 * p->a_3[i];
                    p->v[i] += dt  * p->a[i] + dt2 * p->a_[i];// + dt3 * p->a_2[i] + dt4 * p->a_3[i];
                }
            }
            for(i = 0; i < DIMENSIONS; i++)
            {
                p->xp[i] = p->x[i];
                p->vp[i] = p->v[i];
            }
        }
    }

}
// END
// herm_pred
//

//
// evaluate_1_2
//
// Calculate derivatives of v into a and a_
// for particle at position _pos_.
//
double evaluate_1_2(struct particle parts[], int pcount, int pos, int posmin, int posmax,
                    double a[3], double a_[3])
{
    int new_close_warn = 0;
    int perturb = (posmax>0?1:0);
    struct particle *p = parts+pos;
    double *px = p->xp, *pv = p->vp;
    double maxforce = -1.0, r_vic_2 = .0, r1_2 = scal_prod(px, px), rs_2, r_perturb_2;
    double a0 = .0, a1=.0, a2=.0, a_0=.0, a_1=.0, a_2=.0;
    double phi = .0;

    // determine sphere of vicinity (close encounters possible within 2dt)
    if(perturb)
    {
        if(_sqrt_mratio == .0)
            _sqrt_mratio = sqrt(m_max() / parts[0].m);

        r_vic_2 = _sqrt_mratio * v_abs(px) + (v_abs(p->vp) + sqrt(2. * parts[0].m / v_abs(px))) * 2. * p->dt;
        r_vic_2 *= r_vic_2;
    }
    else
        r_vic_2 = v_abs(p->vp) * 2 * p->dt;

    // calculate Schwarzschild sum of radii:  3 * rs = 3 * (2 * G * m1 / c� + 2 * G * m2 / c�)
    rs_2 =  9. * C_2G_C2 * C_2G_C2 * (p->m + m_max()) * (p->m + m_max());

    // calculate perturbing distance:
    //    F_perturb > fact * F_central, i. e.
    //    m2/r_2   > fact * m0/r1_2     or    r_2/m2 < r1_2 / fact / m0
    r_perturb_2 = r1_2 * m_max() / (PERTURBING_FORCE_RATIO * parts[0].m);

    // evaluate forces

    // approaching SMBH??
    if(px[0]*px[0] + px[1]*px[1] + px[2]*px[2] < 9. * C_2G_C2 * C_2G_C2 * (p->m + parts->m) * (p->m + parts->m))
    {
        // collision in 3 Schwarzschild-radii
        fprintf(get_file(FILE_WARNING), "#### [t=%1.12e] COLLISION of SMBH m0 and m%d: %e (r_S = %e) ####\n",
                t_total(p->t),
                p->name,
                convert_length(v_abs(px), 0),
                convert_length(C_2G_C2 * (parts->m + p->m), 0));
        fflush(get_file(FILE_WARNING));
        add_collision(p, parts);
    }

    //if (posmax - posmin > pcount*0.5)
    //{
    //    iteration_split(parts, pcount, pos, posmin,posmax, p,
    //                  px, pv,
    //                  &a0, &a1, &a2, &a_0, &a_1, &a_2, r_vic_2, &phi,
    //                  maxforce, r_perturb_2, r1_2, perturb, rs_2,
    //                  &new_close_warn);
    //}
    //else
    //{
        iteration_normal(parts, pcount, pos, posmin,posmax, p,
                      px, pv,
                      &a0, &a1, &a2, &a_0, &a_1, &a_2, r_vic_2, &phi,
                      maxforce, r_perturb_2, r1_2, perturb, rs_2,
                      &new_close_warn);
    //}

    a[0] = a0; a_[0] = a_0;
    a[1] = a1; a_[1] = a_1;
    a[2] = a2; a_[2] = a_2;

    //printf("%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n",a[0],a[1],a[2],a_[0],a_[1],a[2],phi);
    p->phi_stars = phi;
    p->phi_bgr = .0;

    if(perturb)
        if(!new_close_warn)
        {
            p->io_close_warn = -1.;
        }
        return maxforce;
}
// END
// evaluate_1_2



//
// step_hermite_2
//
// Calculate hermite step.
// _comp_ = 1 for composite, _comp_ = 0 otherwise.
int step_hermite_2(struct particle parts[], int *pcount, double eta, double min_evals, double *t_eval)
{
    int i, j, k, removecount=0; //, i;
    double tmin = -1.0, en;
    struct particle *p;
    collisions = 0;
    int re_evaluate = 0;
    struct particle *pk;

    close_n = 0;
    if(close_list == NULL)
        close_list = (struct t_close *)malloc(CLOSE_MAX * sizeof(struct t_close));
    if(node_posmin < 0)
    {
        // initialize integrator at first call
        _sqrt_mratio = sqrt(m_max() / parts[0].m);
        node_posmin = 1;
        node_posmax = *pcount - 1;
    }

    // find particles to move
    init_dt(parts, *pcount, 0);
    fprintf(stderr,"ITIME: %f\n",*t_eval);

    find_move_particles(parts, *pcount, &tmin);

    // make sure to move at least 1 particle
    assert(movecount > 0);

    // forward all particles to tmin along orbit
    move_kepler(parts, *pcount, tmin);

    // hermite predictor
    herm_pred(parts, *pcount, tmin, 0);

    T_START;

    for(j = 0; j < movecount; j++)
    {
        p = parts + active[j];
        // determine significant energy changes during close encounters
        en = p->m * (.5 * (scal_prod(p->v, p->v) + p->phi_stars) + p->phi_bgr - parts[0].m / v_abs(p->x));
        if(p->energy != .0)
        {
        //en = get_energy(parts, *pcount, p - parts);
        if(fabs(en - p->energy) > WARN_ENERGY_FACT * fabs(p->energy))
            fprintf(get_file(FILE_DEBUG),
                "#### [t=%1.12e+%1.4e] particle m%d gained %1.2f%% energy from %e to %e, at [%e\t%e\t%e\t%e\t%e\t%e\t] ####\n",
                t_total(tmin-p->dt), convert_time(p->dt, 0),
                p->name,
                100.*(en/p->energy-1.),
                p->energy, en,
                convert_length(p->x[0], 0), convert_length(p->x[1], 0), convert_length(p->x[2], 0),
                convert_length(convert_time(p->v[0], 1), 0),
                convert_length(convert_time(p->v[1], 1), 0),
                convert_length(convert_time(p->v[2], 1), 0));
        }
        p->energy = en;
    }

    // hermite evaluator
    for(j = 0; j < movecount; j++)
        evaluate_1_2(parts, *pcount, active[j],node_posmin, node_posmax, forceterm[j].a, forceterm[j].a_);

    sumforce(parts);

    hermite_correct(parts, *pcount);

    // track apocentre and pericentre passages
    for(j = 0; j < movecount; j++)
    {
        p = parts + active[j];
        double xv_new = scal_prod(p->x, p->v);
        if(1)//N_MAX_DETAIL < -1 || p->name <= N_MAX_DETAIL)
        {
            if(p->xv < .0 && xv_new >= .0)
            // pericentre passage
            {
                // linear approximation of time of closest encounter
                double t_close = -scal_prod(p->x, p->v) / (scal_prod(p->v, p->v) + scal_prod(p->x, p->a) + scal_prod(p->x, p->ha)), r=.0;
                int i;
                for(i = 0; i < 3; i++)
                    r += square(p->x[i] + t_close * (p->v[i] + .5 * t_close * (p->a[i] + p->ha[i] + _1_3 * t_close * (p->a_[i] + p->ha_[i]))));
                r = sqrt(r);

                if(p->r        < r) r = p->r;
                if(v_abs(p->x) < r) r = v_abs(p->x);
                p->r_peri = r;

                if(0)
                    if(N_MAX_DETAIL < -1 || p->name <= N_MAX_DETAIL)
                        fprintf(get_file(FILE_OTHER),
                            "PP  %e\t%d\t%e\t\t%e\t%e\t%e\t%e\n",
                            t_total(tmin), p->name, convert_length(r, 0),
                            convert_length(p->curr_a, 0), p->curr_e, convert_length(p->rmin, 0), convert_length(p->rmax, 0));
            }
            else if(p->xv > .0 && xv_new <= .0)
            // apocentre passage
            {
                double t_close = -scal_prod(p->x, p->v) / (scal_prod(p->v, p->v) + scal_prod(p->x, p->a) + scal_prod(p->x, p->ha)), r=.0;
                int i;
                for(i = 0; i < 3; i++)
                    r += square(p->x[i] + t_close * (p->v[i] + .5 * t_close * (p->a[i] + p->ha[i] + _1_3 * t_close * (p->a_[i] + p->ha_[i]))));
                r = sqrt(r);
                if(p->r        > r) r = p->r;
                if(v_abs(p->x) > r) r = v_abs(p->x);
                p->r_apo = r;
            }
        }
        p->xv = xv_new;
        p->r = v_abs(p->x);
        if(p->r < p->rmin)
            p->rmin = p->r;
        if(p->r > p->rmax)
            p->rmax = p->r;
    }
    fflush(get_file(FILE_OTHER));

    find_timesteps(parts, *pcount, tmin, eta, min_evals);

    // detect fast approaches
    check_app(parts, *pcount, tmin);

    // correct timesteps
    correct_timestep(parts, *pcount, tmin);

    for(j = 0; j < movecount; j++)
    {
        p = parts + active[j];
        p->active = 0;

        if(p->name == P_IMBH)
        {
            static double xv_imbh = 1.;
            double xv_new = scal_prod(p->x, p->v);
            if(xv_new <= .0 && xv_imbh > .0)
            {
                fprintf(get_file(FILE_OTHER),
                    "%1.9e\t%e\t%e\t%e\t%e\n",
                    t_total(p->t),
                    convert_length(p->curr_a, 0),
                    p->curr_e,
                    convert_length(v_abs(p->x), 0),
                    .0
                    );
                fflush(get_file(FILE_OTHER));
            }
            xv_imbh = xv_new;
        }

        double x = v_abs(p->x);

        if(p->nearestneighbour > 0 && p->nearestneighbour < *pcount && p != parts + p->nearestneighbour)
        {
            pk = parts + p->nearestneighbour;
            if(pk->active &&
                (.5 * (p->m * scal_prod(p->v, p->v) + (pk->m * scal_prod(pk->v, pk->v))) - p->m * pk->m / v_dist(p->x, pk->x, 1.) < 0
                // m1*v1_^2/2 + m2*v2_^2/2 = .5 * m1*m2 / (m1 + m2) * |v1-v2|^2 !< m1*m2 / |r1-r2|
                // (.5 * v_dist(p->v, pk->v, 2) * v_dist(p->x, pk->x, 1) < p->m + pk->m
                //(0
            ))
            {
                add_collision(p, pk);
                fprintf(get_file(FILE_WARNING),
                    "#### [t=%1.12e] particles m%d and m%d formed a binary - merge (r=%e\tE_kin=%e\tE_pot=%e)!\n",
                    t_total(tmin), p->name, pk->name,
                    convert_length(v_dist(p->x, pk->x, 1), 0) * 206265. / 4.7e-3,
                    .5 * (p->m * scal_prod(p->v, p->v) + (pk->m * scal_prod(pk->v, pk->v))),
                    //.5 * p->m * pk->m / (p->m + pk->m) * v_dist(p->v, pk->v, 2),
                    -p->m * pk->m / v_dist(p->x, pk->x, 1));
                fflush(get_file(FILE_WARNING));
            }
        }
        p->nearestneighbour = -1;
        int remove = 0;
        if(x > convert_length(MAX_X, 1))
            remove = 1;
        if(remove)
        {
            fprintf(get_file(FILE_WARNING), "### [%1.12e] m%d OUT OF RANGE - REMOVE: x=%e\tv=%e\te=%e\ta=%e\n",
                t_total(p->t),
                p->name,
                convert_length(x, 0),
                convert_length(convert_time(v_abs(p->v), 1), 0),
                p->curr_e,
                convert_length(p->curr_a, 0));
                fflush(get_file(FILE_WARNING));
        }
        //#endif

        for(i = 0; i < collisions; i++)
            if(coll_vector[i][0] == parts + active[j])
            {
                pk = coll_vector[i][1];
                if(pk != parts)
                {
                    fprintf(get_file(FILE_WARNING),
                        "# MERGE: m1 = %e\tr1 = %e\tm2 = %e\tr2 = %e\tr = %e pc = %e Rsun\tv = %e pc/Myr\n",
                        convert_mass(p->m,  0), .0, convert_mass(pk->m, 0), .0, convert_length(v_dist(p->x, pk->xp, 1), 0),
                        convert_length(v_dist(p->x, pk->xp, 1), 0) * 206265. / 4.7e-3, convert_length(convert_time(v_dist(p->v, pk->vp, 1), 1), 0) * 1.e6
                        );
                    p->energy = .5 * (p->m * scal_prod(p->v, p->v) + pk->m * (scal_prod(pk->vp, pk->vp)))
                      - p->m * pk->m / v_dist(p->x, pk->xp, 1);

                    pk->energy = .5 * p->m * pk->m / v_dist(p->x, pk->xp, 1);
                    pk->x[0] = (pk->m * pk->xp[0] + p->m * p->x[0]) / (p->m + pk->m);
                    pk->x[1] = (pk->m * pk->xp[1] + p->m * p->x[1]) / (p->m + pk->m);
                    pk->x[2] = (pk->m * pk->xp[2] + p->m * p->x[2]) / (p->m + pk->m);
                    pk->v[0] = (pk->m * pk->vp[0] + p->m * p->v[0]) / (p->m + pk->m);
                    pk->v[1] = (pk->m * pk->vp[1] + p->m * p->v[1]) / (p->m + pk->m);
                    pk->v[2] = (pk->m * pk->vp[2] + p->m * p->v[2]) / (p->m + pk->m);

                    if(p->dt < pk->dt)
                      pk->dt = p->dt;
                    pk->htlast = pk->t = tmin;
                    pk->m += p->m;

                    pk->energy += pk->m * (.5 * (scal_prod(pk->v, pk->v) + pk->phi_stars) + p->phi_bgr - parts[0].m / v_abs(pk->x));
                    // only valid if pk->phi is up to date
                    p->energy -= .5 * pk->m * scal_prod(pk->v, pk->v);
                }
                else // COLLISION with SMBH
                {
                    p->energy = p->m * (.5 * scal_prod(p->v, p->v) - pk->m / v_abs(p->x) + p->phi_bgr);
                    pk->m += p->m;
                }
                remove = 1;
                fprintf(get_file(FILE_WARNING),
                    "#### [t=%1.12e] particle m%d swallowed m%d. New mass: %e\n",
                    t_total(tmin), pk->name, p->name, convert_mass(pk->m, 0));
                    fflush(get_file(FILE_WARNING));
                    re_evaluate = 1;
                break; // only one collision per particle allowed
            }

        if(remove)
        {
                remove_part[removecount++] = active[j];
                active[j] = -1;
        }
    }

    if(removecount)
    {
        for(j = 0; j < removecount; j++)
        {
            lose_energy(parts[remove_part[j]].energy);
            fprintf(get_file(FILE_WARNING),
                   "#### [t=%1.12e] removing particle m%d at %d of %d loses energy %e\n",
                    t_total(tmin), parts[remove_part[j]].name, remove_part[j], *pcount, parts[remove_part[j]].energy);
            if(remove_part[j] < --(*pcount))
            {
                memcpy(parts + remove_part[j], parts + *pcount, sizeof(struct particle));
                fprintf(get_file(FILE_WARNING),
                "#### [t=%1.12e] particle at %d replaced by m%d\n",
                t_total(tmin), remove_part[j], parts[remove_part[j]].name);
            }
            node_posmax--;
            fflush(get_file(FILE_WARNING));
            for(k = 0; k < close_n; k++)
                if(close_list[k].p == parts + remove_part[j] || close_list[k].pk == parts + remove_part[j])
                    close_list[k].p = close_list[k].pk = NULL;
        }

        if(re_evaluate)
        {
            fprintf(get_file(FILE_WARNING),
                "#### [t=%1.12e] Need to re-evaluate %d moved particles.\n",
                    t_total(tmin), movecount);
            for(j = 0; j < movecount; j++)
                if(active[j] > 0 && active[j] < *pcount)
                {
                    evaluate_1_2(parts, *pcount, active[j], 0, 0, parts[active[j]].ha, parts[active[j]].ha_);
                    evaluate_1_2(parts, *pcount, active[j], 1, *pcount - 1, parts[active[j]].a, parts[active[j]].a_);
                }

            check_app(parts, *pcount, .0);
        }
        else
            fprintf(get_file(FILE_WARNING),
               "#### [t=%1.12e] No need to re-evaluate %d moved particles.\n",
               t_total(tmin), movecount);
        fflush(get_file(FILE_WARNING));

        init_dt(parts, *pcount, 1);
        fprintf(get_file(FILE_WARNING), "### NEW NUMBER OF PARTICLES: %d\n", *pcount);
        fflush(get_file(FILE_WARNING));
    }

    T_END;

    move_center(parts, *pcount, tmin);
    parts->v[0] = parts->v[1] = parts->v[2] = 0;
    parts->x[0] = parts->x[1] = parts->x[2] = 0;
    parts->vp[0] = parts->vp[1] = parts->vp[2] = 0;
    parts->xp[0] = parts->xp[1] = parts->xp[2] = 0;
    parts->t = tmin;

    // CUDA
    //freeCudaMem();

    return movecount;
}
// END
// step_hermite_2
