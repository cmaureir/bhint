#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <omp.h>
#include "bhi_config.h"

#ifdef  USE_M_ENCL
#define M_ENCL_IO(r) M_ENCL(r)
#else
#define M_ENCL_IO(r) .0
#endif

#define min(x, y) ((x) <= (y) ? (x) : (y))

#define DIMENSIONS     3
#define DT_TOLERANCE   1e-10

#define __C_AU      1.495978707e+11                // m
#define __C_GRAVIT  6.6742e-11                     // m³/kg/s²
#define __C_SOLMASS 1.9891e+30                     // kg
#define __C_SOLMASS_AU3DAY2 2.9591220828559115e-04 // au³/day²
#define __C_2G_C2_DAY2AU2 6.671322399e-5           // day²/au²
#define __C_C_AUDAY   173.1446327                  // au/day

//#define PLAIN_UNITS

#ifdef  PLAIN_UNITS // [TESTING]
#define CONV_M_SOLMASS 1.
#define CONV_X_PARSEC  1.
#define CONV_X_KM      1.
#define CONV_T_YEAR    1.
#define CONV_T_SECOND  1.

#define CONVINT_X_AU_SQRT      1.
#define CONVINT_M_AU3DAY2_SQRT 1.

#else // NOT PLAIN_UNITS

/* mass constants - base unit: au³/day² */
#define CONV_M_SOLMASS (1./__C_SOLMASS_AU3DAY2)

/* length constants - base unit au */
#define CONV_X_PARSEC  4.848136813e-6
#define CONV_X_KM      149597870.7

/* time constants - base unit day */
#define CONV_T_YEAR    0.2737909330e-2
#define CONV_T_SECOND  86400.

#define CONVINT_X_AU_SQRT      32.  //sqrt(1000)
#define CONVINT_M_AU3DAY2_SQRT 100.

#endif //PLAIN_UNITS

#define CONVINT_X_AU           (CONVINT_X_AU_SQRT * CONVINT_X_AU_SQRT)
#define CONVINT_M_AU3DAY2      (CONVINT_M_AU3DAY2_SQRT * CONVINT_M_AU3DAY2_SQRT)
#define CONVINT_T_DAY          (CONVINT_X_AU_SQRT * CONVINT_X_AU_SQRT * CONVINT_X_AU_SQRT / CONVINT_M_AU3DAY2_SQRT)
#define C_2G_C2                (__C_2G_C2_DAY2AU2 * CONVINT_X_AU * CONVINT_X_AU / (CONVINT_T_DAY * CONVINT_T_DAY))
#define C_C                    (__C_C_AUDAY / CONVINT_X_AU * CONVINT_T_DAY)

#define FILE_OUTPUT   1
#define FILE_DETAIL   2
#define FILE_DEBUG    3
#define FILE_WARNING  4
#define FILE_INIT     5
#define FILE_RESTART  6
#define FILE_DUMP     7
#define FILE_OTHER    8

#define EXT_OUTPUT   "out"
#define EXT_DETAIL   "detail"
#define EXT_DEBUG    "debug"
#define EXT_WARNING  "warn"
#define EXT_INIT     "init"
#define EXT_RESTART  "restart"
#define EXT_DUMP     "dump"
#define EXT_OTHER    "other"

#define MAX_STEPSIZE_POWER 100

enum _module {_UL_NBODY=0, _UL_HERMITE2, _UL_IO, _UL_TIMESTEP, _UL_KEPLER};

#define SIGNAL_PRINT_POS  -10

extern int _global_module, _global_function;

#define TRACK_FUNCTIONS

#ifdef TRACK_FUNCTIONS

#define _enter_function(module, function)   \
  int __my_lastmodule   = _global_module;   \
  int __my_lastfunction = _global_function; \
  _global_module        = module;           \
  _global_function      = function;

#define _exit_function()		\
  _global_module   = __my_lastmodule;   \
  _global_function = __my_lastfunction;

#define _print_position(v_kill, v_text, v_1, v_2, v_3, v_f)                     \
  if(v_kill <= SIGNAL_PRINT_POS) {                                              \
    fprintf(get_file(FILE_WARNING), "### WORKING AT %s ### %d\t%d\t%d\t%e\t\n", \
    v_text, v_1, v_2, v_3, v_f);                                                \
    fflush(get_file(FILE_WARNING));                                             \
    v_kill -= SIGNAL_PRINT_POS;}


#else  // NOT TRACK_FUNCTIONS

#define _enter_function(module, function)
#define _exit_function()
#define _print_position()

#endif // TRACK_FUNCTIONS

struct particle {
  double x[3], xc[3], xp[3], xpc[3], dx[3], v[3], vc[3], vp[3], vpc[3], dv[3];
  double a[3], a_[3], a_2[3], a_3[3], an[3], a_n[3];
  double ha[3], ha_[3], ha_2[3], ha_3[3], hdt, htlast;

  double xv, r, rmin, rmax, r_apo, r_peri;

  double phi_stars, phi_bgr;

  double m, m_, eta, dt, dtnext;
  double t, energy;
  double io_tlast;
  double orig_a, orig_e, orig_t, orig_j;
  double curr_a, curr_e, curr_t, curr_j;
  double io_close_warn;
  long io_steps_c, io_steps_p;
  int active, is_jacobi, is_elliptical;

  int nearestneighbour, name;

  // GR
  int use_pn;
  double gr_a[3], gr_a_[3];

  double gr_a_2[3], gr_a_3[3], v_thresh_2;
  int switch_pn;

  int sse_on, sse_kw;
  double sse_mass, sse_mt, sse_epoch, sse_tphys, sse_z, sse_dtm, sse_r;
  int sse_multiple;
};

/*
 * Vector functions.
 */

#define square(x) ((x)*(x))
#define scal_prod(x, y)  (x[0]*y[0] + x[1]*y[1] + x[2]*y[2])
#define v_abs(x) sqrt((x)[0]*(x)[0] + (x)[1]*(x)[1] + (x)[2]*(x)[2])
void vec_prod(double x[3], double y[3], double vp[3]);
double v_dist(double x[3], double y[3], int pot);
double vec_angle(double o[3], double x[3], double y[3]);
double vec_angle_2pi(double o[3], double x[3], double y[3], double j[3]);
double inclination(double r0[3], double v0[3], double r1[3], double v1[3], double r2[3], double v2[3]);

/*
 * Kepler integration functions.
 */

int step_kepler(struct particle parts[], int pcount);
void get_constants(double x[3], double v[3], double m_central, double j[3], double e[3], double *a, double *omega);
void get_reduced(struct particle parts[], int pos, double *red_e, double *red_a, double *red_t, double *red_j);
void step_kepler_1(struct particle parts[], int pcount, int pos, double dt, double *out_a, double *out_a_, double *out_a_2, double *out_a_3, double *curr_a, double *curr_e); ;
void move_center(struct particle parts[], int pcount, double t);
void predict(struct particle parts[], int pcount, double t);
double solve_kepler(double mean, double ecc);
double kepler(const double ecc, double mean_anom);
double get_timestep_simple(double x0[3], double v0[3], double x[3], double v[3]);

/*
 * I/O functions.
 */

double convert_time(double time, int forward);
double convert_length(double length, int forward);
double convert_mass(double mass, int forward);
double convert_energy(double energy, int forward);
double t_total(double t);
double m_max(void);
int get_params(struct particle *parts[], double *time, FILE *infile, int *pcount, double *m_max);
int output(int print, struct particle parts[], int pcount, double t_old,
    double pt, double t, double t_over,
    double orbits, double t_steps, double t_max, double steps, double blocks,
    int orb_count);
int local_io();
void print_header(int filetype, int pcount);
void calc(char **argv);
void init_files(const char *filename, int append);
void close_files();
void create(int argc, char **arg);
void add_over(unsigned long val, unsigned long *count, double *count_over);
void convert_to_nbody_units(struct particle parts[], int pcount);
void lose_energy(double de);
void emit_energy(double de);
void comment_datfile();
void print_conv();
FILE *get_file(int filetype);

/*
 *  Force & energy functions.
 */

double get_energy(struct particle parts[], int pcount, int pos);
double get_ekin(struct particle parts[], int pcount, int pred);
double get_epot(struct particle parts[], int pcount, int pred);
void get_acc(struct particle parts[], int pcount, int pos, double acc[3]);

/*
 *  Functions for Hermite integrator.
 */

double evaluate_1_2(struct particle parts[], int pcount, int pos, int posmin, int posmax,double a[3], double a_[3]);
double get_timestep_aarseth(double a[3], double a_[3], double a_2[3], double a_3[3], double eta);
double get_timestep_central(struct particle *parts, struct particle *p, double min_evals);
void predict_part_hermite2(struct particle *p, double t);
void check_app(struct particle *parts, int pcount, double t);
void add_force_extpot(double x[3], double v[3], double a[3], double a_[3], double *phi);
int step_hermite_2(struct particle parts[], int *pcount, double eta, double min_evals, double *t_eval);

/*
 * DUMP functions.
 */
FILE *write_dump_io();
void read_dump_io(FILE *dumpfile, int tmp);
void inc_stepsize(int n);

/*
 *  Timestep functions.
 */
double normalize_dt(double t, double dt);
double drand(double a, double b);
double gauss();
double get_imf_mass(double *m, double *alpha, int n, double *r, double *f);
double evolv1(int *kw, double *mass, double *mt, double *epoch, double *tphys, double *dtm, double *z, double *r);
int check_fast_approaches(struct particle *parts,
struct particle *p, struct particle *pk);
void get_approach_reduce_stats(double out[2]);
void seedrand();
void setup_imf(double *m, double *alpha, int n, double *r, double *f);
