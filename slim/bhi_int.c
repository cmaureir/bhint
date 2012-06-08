/*                                        *
 * U. Loeckmann                           *
 * Kepler-Hermite integrator              *
 * for N-body problem.                    *
 *                                        */
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>
#include <time.h>
#include <signal.h>
#include "bhi.h"

#define CLOSE_MAX     100          // (initial) maximum number of neighbours; is increased dynamically
//#define DT_AARSETH_ALL
#define WARN_ENERGY_FACT .2
//#define WARN_ENERGYALL
#define MAX_REDUCECOUNT 3

#ifdef PN
#define PN_KEPLER
#define PN_TRACK_ORDER 3
double __1_c  = 1/C_C;
double __1_c2 = 1/(C_C*C_C);
double __1_c3 = 1/(C_C*C_C*C_C);
double __1_c4 = 1/(C_C*C_C*C_C*C_C);
double __1_c5 = 1/(C_C*C_C*C_C*C_C*C_C);
#define PN1
#define PN2
#define PN25

#ifdef PN25
double _17_3 = 17. / 3.;
#endif

#ifdef PN_KEPLER
#define PN_ETA_FACT 1.
#endif

#define NO_USE_MASS_SUM

#else // NOT PN
#define SWITCHON_PN  0
#define SWITCHOFF_PN 1
#endif // NOT PN

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
struct t_forceterm forcesum[N_MAX];

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

enum _function {_UL_HERMITE2_REGISTER_DT=0, _UL_HERMITE2_INIT_DT, _UL_HERMITE2_PREDICT_PART, _UL_HERMITE2_ADD_CLOSE,
                _UL_HERMITE2_ADD_COLLISION, _UL_HERMITE2_EVALUATE_1, _UL_HERMITE2_CHECK_APP, _UL_HERMITE2_CORRECT_TIMESTEP,
                _UL_HERMITE2_FIND_MOVE_PARTICLES, _UL_HERMITE2_GR_FORCE_COM, _UL_HERMITE2_GR_JERK_COM, _UL_HERMITE2_GR_FORCE,
                _UL_HERMITE2_GR_JERK, _UL_HERMITE2_PATH_INTEGRAL, _UL_HERMITE2_MOVE_KEPLER, _UL_HERMITE2_HERMITE_CORRECT,
                _UL_HERMITE2_FIND_TIMESTEPS, _UL_HERMITE2_HERM_PRED,
                _UL_HERMITE2_FIND_NEIGHBOURS, _UL_HERMITE2_STEP_HERMITE, _UL_HERMITE2_XXXXXXXXX,
                _G6_SET_TI, _G6_CALC_FIRSTHALF, _G6_CALC_LASTHALF, _G6_READ_NEIGHBOUR_LIST,
                _G6_GET_NEIGHBOUR_LIST, _UL_HERMITE2_ADD_FORCE_EXTPOT};

//
// add_force_extpot
//
void add_force_extpot(double x[3], double v[3], double a[3], double a_[3], double *phi)
{
    #ifdef EXT_POT
    _enter_function(_UL_HERMITE2, _UL_HERMITE2_ADD_FORCE_EXTPOT);
    double r_2 = scal_prod(x, x);
    double r = sqrt(r_2);

    // a = - G M / r^3 * r_
    double afact = -EP_M(r) / (r_2*r);
    if(a != NULL)
    {
        a[0] += afact * x[0];
        a[1] += afact * x[1];
        a[2] += afact * x[2];
    }

    // a_ = G [ -M / r^5 (v_ * r^2 - 3 * r_*v_ r_) - r_*v_ / r^2 * 4 pi rho r_ ]
    // a_ = -4 pi G rho0 / (3-alpha) * (r^(-alpha) * v_ - alpha * r^(-alpha-2) * r_ * (r_*v_))
    if(a_ != NULL)
    {
        assert(v != NULL);
        afact /= r_2;
        double xv = scal_prod(x, v);
        double xv4pirho_r2 = xv * 4. * M_PI * EP_RHO(r) / r_2;
        a_[0] += afact * (r_2 * v[0] - x[0] * 3. * xv) - xv4pirho_r2 * x[0];
        a_[1] += afact * (r_2 * v[1] - x[1] * 3. * xv) - xv4pirho_r2 * x[1];
        a_[2] += afact * (r_2 * v[2] - x[2] * 3. * xv) - xv4pirho_r2 * x[2];
        afact *= r_2;
    }

    // phi = 4 pi G rho0 / (3-alpha) / (2-alpha) * r^(2-alpha)
    if(phi != NULL)
        *phi += EP_PHI(r);

    _exit_function();
    #endif // EXT_POT
}
// END
// add_force_extpot

//
// register_dt
//
void register_dt(int pos, double dt)
{
    _enter_function(_UL_HERMITE2, _UL_HERMITE2_REGISTER_DT);
    int n = floor(log(1.1 * dt) * _1_LOG2) + MAX_STEPSIZE_POWER;

    assert(n >= 0); assert(n < 2 * MAX_STEPSIZE_POWER);
    if(n < step_min)
        step_min = n;

    //printf("%d[%d:%d]:%e\t", pos, n, step_min, dt);

    if(step_alloc[n] <= step_count[n])
    {
        step_alloc[n] = floor(1.4 * (float)step_alloc[n]);
        step_part[n] = realloc(step_part[n], step_alloc[n] * sizeof(int));
        assert(step_part[n] != NULL);
    }

    step_part[n][step_count[n]++] = pos;
    inc_stepsize(n);
    _exit_function();
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
    _enter_function(_UL_HERMITE2, _UL_HERMITE2_INIT_DT);

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
    _exit_function();
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
    _enter_function(_UL_HERMITE2, _UL_HERMITE2_PREDICT_PART);
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
        // relativistic forces
        if(p->use_pn)
        {
            p->xp[i] += dt2n * p->gr_a[i] + dt3n * p->gr_a_[i] + dt4n * p->gr_a_2[i] + dt5n * p->gr_a_3[i];
            p->vp[i] += dtn  * p->gr_a[i] + dt2n * p->gr_a_[i] + dt3n * p->gr_a_2[i] + dt4n * p->gr_a_3[i];
        }
    }

    _exit_function();
}
// END
// predict_part_hermite2

//
// add_close
//
void add_close(struct particle *parts, struct particle *p, struct particle *pk)
{
    _enter_function(_UL_HERMITE2, _UL_HERMITE2_ADD_CLOSE);
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
            _exit_function();
            return;
        }
        //printf("!!! %d\n", close_max);
        if(close_list != close_2)
            close_list = close_2;
       close_max *= 2;
    }
    close_list[close_n].p = p;
    close_list[close_n++].pk = pk;
    _exit_function();
}
// END
// add_close

//
// add_collision
//
void add_collision(struct particle *p, struct particle *pk)
{
    _enter_function(_UL_HERMITE2, _UL_HERMITE2_ADD_COLLISION);
    int co;
    for(co = 0; co < collisions; co++)
        if(coll_vector[co][0] == p || coll_vector[co][0] == pk || coll_vector[co][1] == p || coll_vector[co][1] == pk)
        {
            _exit_function();
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
        _exit_function();
}
// END
// add_collision

//
// check_app
//
void check_app(struct particle *parts, int pcount, double t)
{
    _enter_function(_UL_HERMITE2, _UL_HERMITE2_CHECK_APP);
    int j, collision;

    for(j = 0; j < close_n; j++)
    {
        if(close_list[j].p == NULL || close_list[j].p - parts >= pcount || close_list[j].pk - parts >= pcount )
            continue;
        #ifdef PRINT_2
        if(close_list[j].p == parts + PRINT_2 || close_list[j].pk == parts + PRINT_2)
            fprintf(get_file(FILE_DEBUG), "# [%1.8e] # # # CLOSE-CHECK : %d - %d\n",
                    t_total(t),
                    close_list[j].p->name,
                    close_list[j].p->name);
        #endif
        collision = check_fast_approaches(parts, close_list[j].p, close_list[j].pk);
        if(collision && close_list[j].pk->active)
            add_collision(close_list[j].p, close_list[j].pk);

    }

     _exit_function();
}
// END
// check_app

//
// correct_timestep
//
void correct_timestep(struct particle *parts, int pcount, double tmin)
{
    _enter_function(_UL_HERMITE2, _UL_HERMITE2_CORRECT_TIMESTEP);
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

        #ifdef PRINT_1
        if(p - parts == PRINT_1 || p - parts == PRINT_2)
        {
            for(i = 0; i < DIMENSIONS; i++)
                r_temp[i] =  p->x[i] - parts[PRINT_2 + PRINT_1 - (p - parts)].x[i]
                                + p->dt * (p->v[i] - parts[PRINT_2 + PRINT_1 - (p - parts)].v[i]
                                + p->dt * (p->ha[i] + p->a[i]
                                - parts[PRINT_2 + PRINT_1 - (p - parts)].ha[i]
                                - parts[PRINT_2 + PRINT_1 - (p - parts)].a[i]
                                + p->dt * (p->ha_[i] + p->a_[i]
                                - parts[PRINT_2 + PRINT_1 - (p - parts)].ha_[i]
                                - parts[PRINT_2 + PRINT_1 - (p - parts)].a_[i])));
            fprintf(get_file(FILE_WARNING),
                    "####  $$$$$   [t=%1.12e] time step %1.2e for approach of m%d and m%d (r0=%1.2e, r1=%1.2e)\n",
                    t_total(p->t),
                    convert_time(p->dt, 0),
                    p - parts, PRINT_1 + PRINT_2 - (p - parts),
                    convert_length(v_dist(p->x, parts[PRINT_1 + PRINT_2 - (p - parts)].x, 1), 0),
                    convert_length(v_abs(r_temp), 0));
            printf("#m%d\tt=%e\tx=( %1.4e %1.4e %1.4e )\t|v|=%1.2e\tr=%1.2e\tdt=%1.3e\thdt=%1.3e\tE1=%1.7e\tE=%1.7e\n",
                    p->name,
                    t_total(p->t),
                    convert_length(p->x[0], 0), convert_length(p->x[1], 0), convert_length(p->x[2], 0),
                    v_abs(p->v),
                    convert_length(sqrt((  parts[PRINT_1].x[0]-parts[PRINT_2].x[0])
                         * (parts[PRINT_1].x[0]-parts[PRINT_2].x[0])
                         + (parts[PRINT_1].x[1]-parts[PRINT_2].x[1])
                         * (parts[PRINT_1].x[1]-parts[PRINT_2].x[1])
                         + (parts[PRINT_1].x[2]-parts[PRINT_2].x[2])
                         * (parts[PRINT_1].x[2]-parts[PRINT_2].x[2])), 0),
                    convert_time(p->dt, 0),
                    convert_time(p->hdt, 0),
                    get_energy(parts, pcount, p - parts),
                    get_energy(parts, pcount, PRINT_1) + get_energy(parts, pcount, PRINT_2));
            printf("#\t|a|=%1.2e\t|a_|=%1.2e\t|ha|=%1.2e\t|ha_|=%1.2e\ta/a_=%1.2e\tha/ha_=%1.2e\tm=%1.2e\n",
                    v_abs(p->a), v_abs(p->a_),
                    v_abs(p->ha), v_abs(p->ha_),
                    sqrt(v_abs(p->a) / v_abs(p->a_)),
                    sqrt(v_abs(p->ha) / v_abs(p->ha_)),
                    p->m);
            printf("[[%e]]\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
                    t_total(p->t),
                    convert_length(parts[PRINT_1].x[0], 0),  convert_length(parts[PRINT_1].x[1], 0),  convert_length(parts[PRINT_1].x[2], 0),
                    convert_length(parts[PRINT_2].x[0], 0),  convert_length(parts[PRINT_2].x[1], 0),  convert_length(parts[PRINT_2].x[2], 0),
                    convert_length(parts[PRINT_1].xp[0], 0), convert_length(parts[PRINT_1].xp[1], 0), convert_length(parts[PRINT_1].xp[2], 0),
                    convert_length(parts[PRINT_2].xp[0], 0), convert_length(parts[PRINT_2].xp[1], 0), convert_length(parts[PRINT_2].xp[2], 0));
            printf("#         %e\n", 1. / v_dist(parts[PRINT_1].x, parts[PRINT_2].x, 2));
            printf("# CORRECT %5d:\t|x|=%e\t|v|=%e\n", p - parts, convert_length(v_abs(p->x), 0), v_abs(p->v));
        }
        fflush(stdout);
        #endif

        register_dt(active[j], p->dt);
    }
    _exit_function();
}
// END
// correct_timestep

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


// find_move_particles
//
void find_move_particles(struct particle *parts, int pcount, double *tmin_out)
{
    _enter_function(_UL_HERMITE2, _UL_HERMITE2_FIND_MOVE_PARTICLES);
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
    _exit_function();
}
// END
// find_move_particles

//
// gr_force_com
//
#ifdef PN
void gr_force_com(double m2, double m1, double _x[3], double v[3], double a[3])
{
    _enter_function(_UL_HERMITE2, _UL_HERMITE2_GR_FORCE_COM);
    //double __v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    double m, nu;
    int k;
    #ifdef NO_USE_MASS_SUM
    m = m2;
    nu = m1 / m2;
    #else
    m = m1 + m2;
    nu = m1 * m2 / (m * m);
    #endif
    double _1_r2 = 1. / (_x[0] * _x[0] + _x[1] * _x[1] + _x[2] * _x[2]);
    double _1_r  = sqrt(_1_r2), _1_r3 = _1_r * _1_r2;
    double n[3]  = {_x[0] * _1_r, _x[1] * _1_r, _x[2] * _1_r};
    double v_2   = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    double nv    = n[0] * v[0] + n[1] * v[1] + n[2] * v[2];
    #ifdef PN1
    double pn1_a = -v_2  * (1 + 3.*nu) + 1.5 * nv*nv * nu + m * _1_r * (4. + nu + nu);
    double pn1_b = nv * (4. - nu - nu);
    #endif
    #ifdef PN2
    double nv2 = nv * nv;
    double pn2_a = -2.*v_2*v_2 + 1.5*v_2*nv2*(3.-4.*nu) - 1.875*nv2*nv2*(1-3.*nu);
    double pn2_b = .5*v_2*nu*(11.+4.*nu) + 2.*nv2*(1. + nu*(12.+3.*nu));
    double pn2_c = (8.*v_2 - 1.5*nv2*(3.+nu+nu))*nv;
    #endif
    #ifdef PN25
    double pn25_a = v_2 + 3. * m * _1_r;
    double pn25_b = 3. * v_2 + _17_3 * m * _1_r;
    #endif
    for(k = 0; k < 3; k++)
    {
        a[k] = .0;
        #ifdef PN1
        a[k] += m * _1_r2 * (n[k] * pn1_a + v[k] * pn1_b) * __1_c2;
        #endif
        #ifdef PN2
        a[k] += m * _1_r2 * (n[k] * (nu * pn2_a + m * _1_r * pn2_b - m * m * _1_r2 * (9. + 21.75 * nu))
                   + v[k] * (nu * pn2_c - .5 * m * _1_r * nv * (4. + 43.*nu))) * __1_c4;
        #endif
        #ifdef PN25
        a[k] -= 1.6 * m * m * _1_r3 * nu * (v[k] * pn25_a - n[k] * nv * pn25_b) * __1_c5;
        #endif
    }

    _exit_function();
}
// END
// gr_force_com

//
// gr_jerk_com
//
void gr_jerk_com(double m2, double m1, double _x[3], double v[3], double an[3], double ha[3], double gr_a[3], double a_[3])
{
    _enter_function(_UL_HERMITE2, _UL_HERMITE2_GR_JERK_COM);
    //double __v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    double m, nu;
    // double a[3] = {an[0] + ha[0] + gr_a[0], an[1] + ha[1] + gr_a[1], an[2] + ha[2] + gr_a[2]};
    double a[3] = {ha[0] + gr_a[0], ha[1] + gr_a[1], ha[2] + gr_a[2]};
    int k;

    #ifdef NO_USE_MASS_SUM
    m = m2;
    nu = m1 / m2;
    #else
    m = m1 + m2;
    nu = m1 * m2 / (m * m);
    #endif
    double _1_r2 = 1. / (_x[0] * _x[0] + _x[1] * _x[1] + _x[2] * _x[2]);
    double _1_r  = sqrt(_1_r2), _1_r3 = _1_r * _1_r2;
    double n[3]  = {_x[0] * _1_r,  _x[1] *  _1_r,  _x[2] * _1_r};
    double v_2   =  v[0] *  v[0] +  v[1] *  v[1] +  v[2] * v[2];
    double v_2_  =    2. * (v[0] *  a[0] +  v[1] *  a[1] + v[2] * a[2]);
    double nv    =  n[0] *  v[0] +  n[1] *  v[1] +  n[2] * v[2];
    double na    =  n[0] *  a[0] +  n[1] *  a[1] +  n[2] * a[2];
    double xv    = _x[0] *  v[0] + _x[1] *  v[1] + _x[2] * v[2];
    double n_[3] = {v[0] *  _1_r -    xv * _x[0] * _1_r3,  v[1] * _1_r - xv * _x[1] * _1_r3, v[2] * _1_r - xv * _x[2] * _1_r3};
    double n_v   = n_[0] *  v[0] + n_[1] *  v[1] + n_[2] * v[2];
    double nv_   = na + n_v;
    #ifdef PN1
    double pn1_a  = -v_2  * (1 + 3.*nu) + 1.5 * nv *  nv * nu + m * _1_r * (4. + 2.*nu);
    double pn1_a_ = -v_2_ * (1 + 3.*nu) +  3. * nv * nv_ * nu - m *   xv * _1_r3 * (4. + 2.*nu);
    double pn1_b  =  nv  * (4. - 2.*nu);
    double pn1_b_ =  nv_ * (4. - 2.*nu);
    #endif
    #ifdef PN2
    double nv2 = nv * nv;
    double pn2_a1  = -2. *  v_2 *  v_2 + 1.5*v_2*nv2*(3.-4.*nu) - 1.875*nv2*nv2*(1-3.*nu);
    double pn2_a1_ = -4. *  v_2 * v_2_ + 1.5*v_2_*nv2*(3.-4.*nu) + 3.*v_2*nv*nv_*(3.-4.*nu) - 7.5*nv2*nv*nv_*(1-3.*nu);
    double pn2_a2  = .5  *  v_2 *   nu * (11.+4.*nu)  + 2.*nv2*(1. + nu*(12.+3.*nu));
    double pn2_a2_ = .5  * v_2_ *   nu * (11.+4.*nu) + 4.*nv*nv_*(1. + nu*(12.+3.*nu));
    double pn2_a   = nu  * pn2_a1  + m * _1_r * pn2_a2 - m * m * _1_r2 * (9. + 21.75 * nu);
    double pn2_a_  = nu  * pn2_a1_ + m * _1_r * pn2_a2_ - m * _1_r3 * xv * pn2_a2 + 2. * m * m * _1_r2 * _1_r2 * xv * (9. + 21.75 * nu);
    double pn2_b1  =  8. * v_2  * nv - 1.5*nv2*nv*(3.+2.*nu);
    double pn2_b1_ =  8. * v_2_ * nv + 8.*v_2*nv_ - 4.5*nv2*nv_*(3.+nu+nu);
    double pn2_b   = nu  * pn2_b1  - .5 * m * _1_r *  nv * (4. + 43.*nu);
    double pn2_b_  = nu  * pn2_b1_ - .5 * m * _1_r * nv_ * (4. + 43.*nu) + .5 * m * _1_r3 * xv * nv * (4. + 43.*nu);
    #endif
    #ifdef PN25
    double pn25_a  = v_2 + 3. * m * _1_r;
    double pn25_a_ = v_2_ - 3. * m * xv * _1_r3;
    double pn25_b  = 3. * v_2 + _17_3 * m * _1_r;
    double pn25_b_ = 3. * v_2_ - _17_3 * m * xv * _1_r3;
    #endif
    for(k = 0; k < 3; k++)
    {
        a_[k] = .0;
        #ifdef PN1
        a_[k] += m * _1_r2 * (n_[k] * pn1_a + n[k] * pn1_a_
                    + a[k] * pn1_b + v[k] * pn1_b_ - 2. * _1_r2 * xv * (n[k] * pn1_a + v[k] * pn1_b)) * __1_c2;
        #endif
        #ifdef PN2
        a_[k] += m * _1_r2 * (n_[k] * pn2_a + n[k] * pn2_a_ + a[k] * pn2_b + v[k] * pn2_b_
                    - 2. * xv * _1_r2 * (n[k] * pn2_a + v[k] * pn2_b)) * __1_c4;
        #endif
        #ifdef PN25
        a_[k] -= 1.6*m*m*_1_r3*nu * (a[k] * pn25_a + v[k] * pn25_a_ - n_[k] * nv * pn25_b - n[k] * nv_ * pn25_b - n[k] * nv * pn25_b_
                       - 3. * xv * _1_r2 * (v[k] * pn25_a - n[k] * nv * pn25_b)) * __1_c5;
        #endif
    }
    _exit_function();
}
#endif // PN
// END
// gr_jerk_com

//
// path_integral
//
void path_integral(int order, double dt, struct particle *p, double sign)
{
    _enter_function(_UL_HERMITE2, _UL_HERMITE2_PATH_INTEGRAL);

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
    _exit_function();
}
// END
// path_integral


//
// move_kepler
//
void move_kepler(struct particle *parts, int pcount, double tmin)
{
    _enter_function(_UL_HERMITE2, _UL_HERMITE2_MOVE_KEPLER);
    struct particle *p;
    double *px, *pxp, *pv, *pvp;
    double dt;
    int j;
    #ifdef PN_KEPLER
    double gr_an[3], gr_a_n[3], dt2, dt3, dt4, dt5, _1_dt2, _1_dt3;
    int i;
    #endif

    for(j = 0; j < movecount; j++)
        if(parts[active[j]].t <= tmin - DT_TOLERANCE)
        {
            p = parts + active[j];
            px = p->x; pv = p->v; pxp = p->xp; pvp = p->vp;
            //p->io_steps_c++;
            pxp[0] = px[0]; pxp[1] = px[1]; pxp[2] = px[2];
            pvp[0] = pv[0]; pvp[1] = pv[1]; pvp[2] = pv[2];

            dt = tmin - p->t;

            #ifdef PN_KEPLER
            if(p->use_pn)
            {
                dt *= PN_KEPLER_FACT;
                if(p->switch_pn)
                {
                    gr_force_com(parts->m, p->m, pxp, pvp, p->gr_a);
                    gr_jerk_com (parts->m, p->m, pxp, pvp, p->a, p->ha, p->gr_a, p->gr_a_);
                }
            }
            #endif

            for(;p->t < tmin; p->t += dt)
            {
                p->io_steps_c++;
                #ifdef PN_KEPLER
                if(p->use_pn)
                {
                    path_integral((PN_TRACK_ORDER > 3 && p->switch_pn) ? 3 : PN_TRACK_ORDER, dt, p, +1.);
                }
                #endif

                step_kepler_1(parts, pcount, p - parts, dt, p->ha, p->ha_, p->ha_2, p->ha_3, &(p->curr_a), &(p->curr_e));

                #ifdef PN_KEPLER
                if(p->use_pn)// && !p->switch_pn)
                {
                    dt2    = dt * dt  * .5; dt3 = dt * dt2 * _1_3; dt4 = dt * dt3 * .25; dt5 = dt * dt4 * .2;
                    _1_dt3 = 1. / dt3; _1_dt2 = _1_dt3 * dt * _1_3;

                    // hermite predictor
                    for(i = 0; i < DIMENSIONS; i++)
                    {
                        pxp[i] += dt2 * p->gr_a[i] + dt3 * p->gr_a_[i];// + dt4 * p->gr_a_2[i] + dt5 * p->gr_a_3[i];
                        pvp[i] += dt  * p->gr_a[i] + dt2 * p->gr_a_[i];// + dt3 * p->gr_a_2[i] + dt4 * p->gr_a_3[i];
                    }

                    // hermite evaluator
                    gr_force_com(parts->m, p->m, pxp, pvp, gr_an);

                    //gr_jerk(parts->m, p->m, pxp, pvp, p->a, p->ha, gr_an, gr_a_n);
                    gr_jerk_com(parts->m, p->m, pxp, pvp, p->a, p->ha, gr_an, gr_a_n);

                    // hermite corrector
                    for(i = 0; i < DIMENSIONS; i++)
                    {
                        p->gr_a_2[i] = -(3.*(p->gr_a[i]-gr_an[i]) + dt*(2.*p->gr_a_[i]+gr_a_n[i])) * _1_dt2;
                        p->gr_a_3[i] =  (2.*(p->gr_a[i]-gr_an[i]) + dt*(p->gr_a_[i]+gr_a_n[i])) * _1_dt3;
                        p->gr_a[i]   = gr_an[i];
                        p->gr_a_[i]  = gr_a_n[i];
                        pxp[i]      += dt4 * p->gr_a_2[i] + dt5 * p->gr_a_3[i];
                        pvp[i]      += dt3 * p->gr_a_2[i] + dt4 * p->gr_a_3[i];
                        px[i] = pxp[i];
                        pv[i] = pvp[i];
                    }

                    // GR track energy
                    // add only y1 here, y0 already added above
                    path_integral((PN_TRACK_ORDER > 3 && p->switch_pn) ? 3 : PN_TRACK_ORDER, dt, p, -1.);
                }
                #endif // PN_KEPLER
            }

            px[0] = pxp[0]; px[1] = pxp[1]; px[2] = pxp[2];
            pv[0] = pvp[0]; pv[1] = pvp[1]; pv[2] = pvp[2];
        }
        _exit_function();
}
// END
// move_kepler


//
// hermite_correct
//
void hermite_correct(struct particle *parts, int pcount)
{
    _enter_function(_UL_HERMITE2, _UL_HERMITE2_HERMITE_CORRECT);
    double dt=.0, dt2=.0, dt3=.0, dt4=.0, dt5=.0, _1_dt2=.0, _1_dt3=.0;
    int i, j;
    struct particle *p;

    for(j = 0; j < movecount; j++)
    {
        p = parts + active[j];

        #ifdef PN
        if(p->use_pn && SWITCHOFF_PN)
        {
            p->switch_pn = 1;
            p->use_pn = 0;
            //if(p->name == P_IMBH)
            fprintf(get_file(FILE_WARNING), "# [%1.8e] switching off PN for m%d\t(a=%e\te=%e\tT_inspiral=%e)\n",
                    t_total(p->t + p->dt), p->name,
                    convert_length(p->curr_a, 0), p->curr_e, convert_time(T_INSPIRAL, 0) );
        }
        else if (!p->use_pn && SWITCHON_PN && p->v_thresh_2 <= .0)
        // ensure PN is switched on at increasing speed only, i.e. speed has fallen under threshold since last switchoff
        {
            p->switch_pn = 1;
            p->use_pn = 1;
            p->v_thresh_2 = scal_prod(p->vp, p->vp);
            p->gr_a[0]  = p->gr_a[1]  = p->gr_a[2]  = .0;
            p->gr_a_[0] = p->gr_a_[1] = p->gr_a_[2] = .0;
            //if(p->name == P_IMBH)
            fprintf(get_file(FILE_WARNING), "# [%1.8e] switching on PN for m%d\t(a=%e\te=%e\tT_inspiral=%e)\n",
                t_total(p->t + p->dt), p->name,
                convert_length(p->curr_a, 0), p->curr_e, convert_time(T_INSPIRAL, 0));
        }
        else p->switch_pn = 0;
        if(p->use_pn)
        {
            // ???
        }
        else if(p->v_thresh_2 > 0 && !SWITCHON_PN)
        // reset PN entry speed once speed has fallen below threshold
        p->v_thresh_2 = .0;

        #endif // PN

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
    _exit_function();
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
    _enter_function(_UL_HERMITE2, _UL_HERMITE2_FIND_TIMESTEPS);
    double dt, dtn;
    #ifdef DT_AARSETH_ALL
    double a, a_, a_2, a_3;
    #endif
    struct particle *p;
    int j, exp_;

    for(j = 0; j < movecount; j++)
    {
        p = parts + active[j];

        #ifdef DT_AARSETH_ALL
        a = sqrt((p->a[0] + p->ha[0]) * (p->a[0] + p->ha[0])
               + (p->a[1] + p->ha[1]) * (p->a[1] + p->ha[1])
               + (p->a[2] + p->ha[2]) * (p->a[2] + p->ha[2]));
        a_ = sqrt((p->a_[0] + p->ha_[0]) * (p->a_[0] + p->ha_[0])
               + (p->a_[1] + p->ha_[1]) * (p->a_[1] + p->ha_[1])
               + (p->a_[2] + p->ha_[2]) * (p->a_[2] + p->ha_[2]));
        a_2 = sqrt((p->a_2[0] + p->ha_2[0]) * (p->a_2[0] + p->ha_2[0])
               + (p->a_2[1] + p->ha_2[1]) * (p->a_2[1] + p->ha_2[1])
               + (p->a_2[2] + p->ha_2[2]) * (p->a_2[2] + p->ha_2[2]));
        a_3 = sqrt((p->a_3[0] + p->ha_3[0]) * (p->a_3[0] + p->ha_3[0])
               + (p->a_3[1] + p->ha_3[1]) * (p->a_3[1] + p->ha_3[1])
               + (p->a_3[2] + p->ha_3[2]) * (p->a_3[2] + p->ha_3[2]));
        dt = sqrt(eta * (a * a_2 + a_ * a_) / (a_ * a_3 + a_2 * a_2));

        #else // NOT DT_AARSETH_ALL

        dt = get_timestep_aarseth(p->a, p->a_, p->a_2, p->a_3, eta);
        // don't let perturbation timestep get bigger than ETA_FACT * central timestep

        //dtn = get_timestep_aarseth(p->ha, p->ha_, p->ha_2, p->ha_3, eta) * eta_fact;
        dtn = get_timestep_central(parts, p, min_evals);

        if(dtn < dt)
            dt = dtn;
        #endif // NOT DT_AARSETH_ALL

        // allow timestep to double at maximum
        // GRAPE:    also limit to max. number of neighbours
        if(p->dtnext < dt)
            dt = p->dtnext;

        //p->t = tmin;
        p->htlast = tmin;

        // set to power of 2
        frexp(dt, &exp_);
        p->dt = ldexp(.5, exp_);

    }
    _exit_function();
}
// END
// find_timesteps


//
// herm_pred
//
void herm_pred(struct particle *parts, int pcount, double tmin, int pred_only)
{
    _enter_function(_UL_HERMITE2, _UL_HERMITE2_HERM_PRED);
    double dt, dt2, dt3, dt4, dt5;
    int i;
    struct particle *p;

    for(p = parts + 1; p < parts + pcount; p++)
    {
        if(!p->active && (p->t <= tmin - DT_TOLERANCE))
        {
            dt = tmin - p->t;
            dt2 = .5 * dt * dt; dt3 = dt * dt2 * _1_3; dt4 = .25 * dt * dt3; dt5 = .2 * dt * dt4;
            for(i = 0; i < DIMENSIONS; i++)
                p->xp[i] = p->x[i]+ dt * p->v[i] + dt2 * (p->a[i] + p->ha[i]) + dt3 * p->ha_[i]
                                       + dt4 * p->ha_2[i]+ dt5 * p->ha_3[i];
        }
        else
        {
            if(p->active)
            {
                dt = p->dt;
                dt2 = dt * dt  * .5; dt3 = dt * dt2 * _1_3; dt4 = dt * dt3 * .25; dt5 = dt * dt4 * .2;
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
    _exit_function();
}
// END
// herm_pred

//
// evaluate_1_2
//
// Calculate derivatives of v into a and a_
// for particle at position _pos_.
//
double evaluate_1_2(struct particle parts[], int pcount, int pos, int posmin, int posmax,
                    double a[3], double a_[3])
{
    _enter_function(_UL_HERMITE2, _UL_HERMITE2_EVALUATE_1);
    int new_close_warn = 0;
    int perturb = (posmax>0?1:0);
    double r_2, afact, v_x_;
    struct particle *p = parts+pos, *pk;
    double *px = p->xp, *pv = p->vp, *pkx, *pkv;
    double maxforce = -1.0, r_vic_2 = .0, r1_2 = scal_prod(px, px), _1_over_r2, rs_2, r_perturb_2;
    double a0 = .0, a1=.0, a2=.0, a_0=.0, a_1=.0, a_2=.0;
    double min_r2 = 1.e99;
    double phi = .0;
    double x_0, x_1, x_2, v_0, v_1, v_2;
    double px0=px[0], px1=px[1], px2=px[2], pv0=pv[0], pv1=pv[1], pv2=pv[2];

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

    // calculate Schwarzschild sum of radii:  3 * rs = 3 * (2 * G * m1 / c² + 2 * G * m2 / c²)
    rs_2 =  9. * C_2G_C2 * C_2G_C2 * (p->m + m_max()) * (p->m + m_max());

    // calculate perturbing distance:
    //    F_perturb > fact * F_central, i. e.
    //    m2/r_2   > fact * m0/r1_2     or    r_2/m2 < r1_2 / fact / m0
    r_perturb_2 = r1_2 * m_max() / (PERTURBING_FORCE_RATIO * parts[0].m);

    // evaluate forces

    // approaching SMBH??
    if(px0*px0 + px1*px1 + px2*px2 < 9. * C_2G_C2 * C_2G_C2 * (p->m + parts->m) * (p->m + parts->m))
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

    // FOR
    for(pk = parts + posmin; pk <= parts + posmax; pk++)
    {
        if(pk == p)
            continue;
        pkx = pk->xp; pkv = pk->vp;
        x_0 = pkx[0] - px0; v_0 = pkv[0] - pv0;
        x_1 = pkx[1] - px1; v_1 = pkv[1] - pv1;
        x_2 = pkx[2] - px2; v_2 = pkv[2] - pv2;

        // calculate factors needed
        r_2 = x_0*x_0 + x_1*x_1 + x_2*x_2; //scal_prod(x_, x_);

        if(r_2 < min_r2)
        {
            p->nearestneighbour = pk - parts;
            min_r2 = r_2;
        }

        _1_over_r2 = 1. / r_2;
        afact = pk->m * _1_over_r2;
        if(maxforce < .0 || afact > maxforce)
            maxforce = afact;
        afact *= sqrt(_1_over_r2);
        v_x_ = 3. * _1_over_r2 * (x_0*v_0 + x_1*v_1 + x_2*v_2); //scal_prod(v_, x_);

        if(perturb)
            if(
                ((r_2 < r_perturb_2) &&
                // calculate exact for star mass rather than maximum
                (pk->m * r1_2 > PERTURBING_FORCE_RATIO * parts[0].m * r_2))
            )
            {
                if(p->io_close_warn <= 0 || p->io_close_warn > 16 * r_2)
                {
                    fprintf(get_file(FILE_DEBUG),
                            "#### [t=%1.12e] close encounter of m%d[%d] and m%d[%d]: %e (allowing %e, perturbing at %e) ####\n",
                            t_total(p->t),
                            p->name, (int)(p - parts),
                            pk->name, (int)(pk - parts),
                            convert_length(sqrt(r_2), 0),
                            convert_length(C_2G_C2 * (pk->m + p->m), 0),
                            convert_length(sqrt(pk->m / parts[0].m * scal_prod(px, px)), 0));
                    fprintf(get_file(FILE_DEBUG),
                            " CE %1.12e\t%d\t%e\t%1.10e\t%1.10e\t%1.10e\t%1.10e\t%1.10e\t%1.10e\t%d\t%e\t%1.10e\t%1.10e\t%1.10e\t%1.10e\t%1.10e\t%1.10e\n",
                            t_total(p->t),
                            pk->name, convert_mass(pk->m, 0),
                            convert_length(pkx[0], 0), convert_length(pkx[1], 0), convert_length(pkx[2], 0),
                            convert_length(convert_time(pkv[0], 1), 0), convert_length(convert_time(pkv[1], 1), 0), convert_length(convert_time(pkv[2], 1), 0),
                            p->name, convert_mass(p->m, 0),
                            convert_length(px[0], 0), convert_length(px[1], 0), convert_length(px[2], 0),
                            convert_length(convert_time(pv[0], 1), 0), convert_length(convert_time(pv[1], 1), 0), convert_length(convert_time(pv[2], 1), 0)
                            );
                    fflush(get_file(FILE_DEBUG));
                    p->io_close_warn = r_2;
                }
                new_close_warn = 1;
                p->energy = get_energy(parts, pcount, p - parts);
            }
        if(
         (r_2 < rs_2) &&
         // calculate exact for star mass rather than maximum
         (r_2 < 9. * C_2G_C2 * C_2G_C2 * (p->m + pk->m) * (p->m + pk->m)))
        {
            // collision in 3 Schwarzschild-radii
            fprintf(get_file(FILE_WARNING), "#### [t=%1.12e] COLLISION of m%d and m%d: %e (r_S = %e) ####\n",
                    t_total(p->t),
                    p->name,
                    pk->name,
                    convert_length(sqrt(r_2), 0),
                    convert_length(C_2G_C2 * (pk->m + p->m), 0));
            fflush(get_file(FILE_WARNING));
            add_close(parts, p, pk);
        }

        // find approaching particles in vicinity
        else if(r_2 < r_vic_2 && v_x_ < 0)
            add_close(parts, p, pk);

        #ifdef PRINT_1
        if(perturb)
            if((p->name == PRINT_1 && pk->name == PRINT_2)
                || (pk->name == PRINT_1 && p->name == PRINT_2))
                    printf("# [t=%1.8e][%1.2e] PARTICLES  m%d / m%d :\ta=%e\t1/r²=%e\tr=%e\tr_vic=%e\n",
                        t_total(p->t),
                        convert_time(p->dt, 0),
                        p->name, pk->name,
                        afact * sqrt(r_2),
                        convert_length(convert_length(afact / pk->m * sqrt(r_2), 1), 1),
                        convert_length(sqrt(r_2), 0),
                        convert_length(sqrt(r_vic_2), 0)
                        );
        #endif
        a0  += afact * x_0;
        a_0 += afact * (v_0 - v_x_ * x_0);
        a1  += afact * x_1;
        a_1 += afact * (v_1 - v_x_ * x_1);
        a2  += afact * x_2;
        a_2 += afact * (v_2 - v_x_ * x_2);
        phi -= afact * r_2;
    }
    //
    // END FOR
    //

    a[0] = a0; a_[0] = a_0;
    a[1] = a1; a_[1] = a_1;
    a[2] = a2; a_[2] = a_2;

    p->phi_stars = phi;
    p->phi_bgr = .0;
    add_force_extpot(px, pv, a, a_, &(p->phi_bgr));

    if(perturb)
        if(!new_close_warn)
        {
            p->io_close_warn = -1.;
        }
        _exit_function();
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
    _enter_function(_UL_HERMITE2, _UL_HERMITE2_STEP_HERMITE);
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
                    #ifdef PN
                    convert_time(T_INSPIRAL, 0)
                    #else
                    .0
                    #endif
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
                        convert_mass(p->m,  0), p->sse_r, convert_mass(pk->m, 0), pk->sse_r, convert_length(v_dist(p->x, pk->xp, 1), 0),
                        convert_length(v_dist(p->x, pk->xp, 1), 0) * 206265. / 4.7e-3, convert_length(convert_time(v_dist(p->v, pk->vp, 1), 1), 0) * 1.e6
                        );
                    /*
                      double d_e = .5 * (p->m * scal_prod(p->v, p->v) + pk->m * scal_prod(pk->vp, pk->vp)) - p->m * pk->m / v_dist(p->x, pk->xp, 1);
                      struct particle *d_p1;

                      for(d_p1 = parts; d_p1 < parts + (*pcount); d_p1++)
                      if(d_p1 != p && d_p1 != pk)
                      d_e -= d_p1->m * (p->m / v_dist(p->x, d_p1->xp, 1) + pk->m / v_dist(pk->xp, d_p1->xp, 1));

                      fprintf(get_file(FILE_WARNING), "************** ENERGY: %1.12e\t", d_e);

                    */

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

                    // stellar collision
                    pk->sse_mass += p->sse_mass / ((double) p->sse_multiple);
                    pk->sse_mt += p->sse_mt / ((double) p->sse_multiple);

                    pk->energy += pk->m * (.5 * (scal_prod(pk->v, pk->v) + pk->phi_stars) + p->phi_bgr - parts[0].m / v_abs(pk->x));
                    // only valid if pk->phi is up to date (i.e. pk is active and USE_GRAPE ?)
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

        _exit_function();
        return movecount;
}
// END
// step_hermite_2
