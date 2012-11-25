#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bhi.h"
#include <time.h>

/*                                        *
 * U. Loeckmann                           *
 * Input/output functions                 *
 * for N-body problem.                    *
 *                                        */

//#define DUMP_BEFORE_ERR [TESTING]
#define OUTPUT_INITDATA       0    // maximum particle name to print initial data for; -1 for all
#define MAX_FILES 20
//#define _CONV_M 1. //CONV_M_SOLMASS
//#define _CONV_X 1. //CONV_X_PARSEC
//#define _CONV_T 1. //CONV_T_YEAR
#define _CONV_M (CONV_M_SOLMASS * CONVINT_M_AU3DAY2)
#define _CONV_X (CONV_X_PARSEC  * CONVINT_X_AU)
#define _CONV_T (CONV_T_YEAR    * CONVINT_T_DAY)
#define _CONV_V (_CONV_X / _CONV_T)
#define _CONV_E (_CONV_M * _CONV_M / _CONV_X)
//#define PRINT_DETAIL_AB            // in addition to _XV, print also semi-major and -minor vectors a, b

static char *filename;
static double t_last=.0, e_last=.0;
static double blocks_last=.0, steps_last=.0;
static FILE *file[MAX_FILES];
static int dump_no=0;
static double _conv_m=_CONV_M, _conv_x=_CONV_X, _conv_t=_CONV_T, _conv_v=_CONV_V, _conv_e=_CONV_E;
extern double lost_energy;
static double emitted_energy=.0;
static int detailout = 0;

/*
 *  Track energy balance for energy removed from the system (i.e., mergers, out-of-range particles, or emission due to stellar evolution)
 */
void lose_energy(double de)
{
    lost_energy += de;
}

/*
 *  Track energy balance for energy emitted from the system as gravitational waves
 */
void emit_energy(double de)
{
    lost_energy += de;
    emitted_energy += de;
}

/*
 *  Return subsequent file name from continuation
 */
char *get_next_filename(char *fname)
{
    char *suffix=strrchr(fname, '.'), *filename=(char *)malloc(strlen(fname)+10);
    int no;

    strcpy(filename, fname);
    if(suffix != NULL && (no = strtol(suffix + 1, NULL, 10)))
    {
        sprintf(filename + (suffix - fname), ".%.4d", no + 1);
    }
    else
    {
        strcat(filename, ".0001");
    }

    return filename;
}

/*
 *  Return output file for specified type
 */
FILE *get_file(int filetype)
{
    return file[filetype];
}

/*
 *  Open specified project file for writing
 */
FILE *open_file(const char *ext, int append, int no)
{
    char *newname = (char *)malloc(strlen(filename)+50);

    strcpy(newname, filename);
    if(no > 0)
    {
        char n[10];
        sprintf(n, "_%d", no);
        strcat(newname, n);
    }
    strcat(newname, ".");
    strcat(newname, ext);

    return fopen(newname, append ? "a" : "w");
}

/*
 *  Open all project files for writing
 */
void init_files(const char *fname, int append)
{
    int i;

    filename = (char *)malloc(strlen(fname)+1);
    strcpy(filename, fname);

    for(i = 0; i < MAX_FILES; i++)
    {
        file[i] = NULL;
    }

    file[FILE_OUTPUT]  = open_file(EXT_OUTPUT,  append, 0);
    file[FILE_DETAIL]  = open_file(EXT_DETAIL,  append, 0);
    file[FILE_DEBUG]   = open_file(EXT_DEBUG,   append, 0);
    file[FILE_WARNING] = open_file(EXT_WARNING, append, 0);
    file[FILE_OTHER]   = open_file(EXT_OTHER,   append, 0);
}

/*
 *  Close all open project files
 */
void close_files()
{
    int i;

    for(i = 0; i < MAX_FILES; i++)
        if(file[i] != NULL)
            fclose(file[i]);
}


/*
 *  Write I/O configuration to dump file
 */
FILE* write_dump_io()
{
    FILE *dumpfile = open_file(EXT_DUMP, 0, ++dump_no);
    fprintf(dumpfile,
            "IO-VALUES %s\n%1.32le\t%1.32le\t%1.32le\t%1.32le\t%1.32le\t%1.32le\t%1.32le\t%1.32le\t%1.32le\n",
            filename,
            _conv_m, _conv_x, _conv_t, _conv_v, _conv_e,
            t_last, e_last,
            blocks_last, steps_last);

    return dumpfile;
}

/*
 *  Read I/O configuration from dump file
 */
void read_dump_io(FILE *dumpfile, int tmp)
{
    char fname[3000];

    fscanf(dumpfile,
            "%*s %s %le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le",
            fname,
            &_conv_m, &_conv_x, &_conv_t, &_conv_v, &_conv_e,
            &t_last, &e_last,
            &blocks_last, &steps_last);

    init_files(tmp ? "/tmp/uln_t.0001" : get_next_filename(fname), 0);
}

/*
 *  Convert mass from internal value to external output value (forward=0) or vice versa (forward=1)
 */
double convert_mass(double mass, int forward)
{
    return forward ? mass / _conv_m : mass * _conv_m;
}

/*
 *  Convert time from internal value to external output value (forward=0) or vice versa (forward=1)
 */
double convert_time(double time, int forward)
{
    return forward ? time / _conv_t : time * _conv_t;
}

/*
 *  Convert energy from internal value to external output value (forward=0) or vice versa (forward=1)
 */
double convert_energy(double energy, int forward)
{
    return forward ? energy / _conv_e : energy * _conv_e;
}

/*
 *  Convert distance from internal value to external output value (forward=0) or vice versa (forward=1)
 */
double convert_length(double length, int forward)
{
    return forward ? length / _conv_x : length * _conv_x;
}

/*
 *  Read problem's initial parameters.
 */
int get_params(struct particle *parts[], double *time, FILE *infile, int *pcount, double *m_max)
{
    int i, j, k;

    *m_max = .0;
    *pcount = 2;
    if(*time < 0)
        *time = convert_time(*time, 1);

    double m_tot=.0, eta=.0, _1_r;
    char s[300];

    if(infile == NULL)
    {
        return -1;
    }

    *pcount = 0;
    for(i=0; !feof(infile);)
    {
        if(fgets(s, 299, infile) == NULL) continue;
        if(s[0] != '#' && s[0] != '\n') (*pcount)++;
    }

    // initialize arrays
    *parts = (struct particle *)malloc(*pcount * sizeof(struct particle));
    //read data from file and return

    rewind(infile);
    for(j = 0; !feof(infile); )
    {
        if(fgets(s, 299, infile) == NULL) continue;
        if(s[0] != '#' && s[0] != '\n')
        {
            sscanf(s, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
                    &((*parts)[j].m),
                    (*parts)[j].x, (*parts)[j].x + 1, (*parts)[j].x + 2,
                    (*parts)[j].v, (*parts)[j].v + 1, (*parts)[j].v + 2);

                    if((*parts)[j].m < 0)
                    {
                        (*parts)[j].m = -(*parts)[j].m;
                    }

                    //#endif
                    (*parts)[j].m /= _conv_m;
                    (*parts)[j].x[0] /= _conv_x;
                    (*parts)[j].x[1] /= _conv_x;
                    (*parts)[j].x[2] /= _conv_x;
                    (*parts)[j].v[0] /= _conv_v;
                    (*parts)[j].v[1] /= _conv_v;
                    (*parts)[j].v[2] /= _conv_v;

                    if(j > 0 && (*parts)[j].m > *m_max
                        #ifdef P_IMBH
                        && j != P_IMBH
                        #endif
                    )

                        *m_max = (*parts)[j].m;
                        j++;
        }
    }

    // output data found
    if(OUTPUT_INITDATA < 0 || *pcount <= OUTPUT_INITDATA)
        for(j = 0; j < *pcount; j++)
            fprintf(get_file(FILE_OUTPUT),
                    "# m%d:\t%e\t%-1.5e\t%-1.5e\t%-1.5e\t%-1.5e\t%-1.5e\t%-1.5e\n",
                    j,
                    (*parts)[j].m * _conv_m,
                    (*parts)[j].x[0] * _conv_x, (*parts)[j].x[1] * _conv_x, (*parts)[j].x[2] * _conv_x,
                    (*parts)[j].v[0] * _conv_v, (*parts)[j].v[1] * _conv_v, (*parts)[j].v[2] * _conv_v);

            // reset other values

            for(j = 0; j < *pcount; j++)
            {
                for(i = 0; i < DIMENSIONS; i++)
                {
                    (*parts)[j].name     = j;
                    (*parts)[j].xp[i]    = (*parts)[j].x[i];
                    (*parts)[j].vp[i]    = (*parts)[j].v[i];
                    (*parts)[j].a[i]     = .0;
                    (*parts)[j].a_[i]    = .0;
                    (*parts)[j].an[i]    = .0;
                    (*parts)[j].a_n[i]   = .0;
                    (*parts)[j].ha[i]    = .0;
                    (*parts)[j].ha_[i]   = .0;
                    (*parts)[j].ha_2[i]  = .0;
                    (*parts)[j].ha_3[i]  = .0;
                    (*parts)[j].gr_a[i]  = .0;
                    (*parts)[j].gr_a_[i] = .0;
                    (*parts)[j].r      = v_abs((*parts)[j].x);
                    (*parts)[j].xv     = scal_prod((*parts)[j].x, (*parts)[j].v);
                    (*parts)[j].rmin   = 1.e99;
                    (*parts)[j].rmax   = .0;
                    (*parts)[j].r_apo  = -1.;
                    (*parts)[j].r_peri = -1.;
                }

                (*parts)[j].v_thresh_2 = .0;
                (*parts)[j].phi_stars    = .0;
                (*parts)[j].phi_bgr    = .0;
                (*parts)[j].t      = .0;
                (*parts)[j].dt     = .0;
                (*parts)[j].hdt    = .0;
                (*parts)[j].htlast = .0;
                (*parts)[j].dtnext = .0;
                (*parts)[j].m_ = (j == 0 ? m_tot : eta * (*parts)[j].m / (eta + (*parts)[j].m));
                eta += (*parts)[j].m;
                (*parts)[j].eta = eta;
                (*parts)[j].io_tlast   = .0;
                (*parts)[j].io_steps_c = 0;
                (*parts)[j].io_steps_p = 0;
                (*parts)[j].is_elliptical = 1;
                (*parts)[j].io_close_warn = -1.;

                if(j > 0)
                {
                    get_reduced(*parts, j, &((*parts)[j].orig_e), &((*parts)[j].orig_a), &((*parts)[j].orig_t), &((*parts)[j].orig_j));
                }
            }

            for(j = 1; j < *pcount; j++)
            {
                for(k = j + 1; k < *pcount; k++)
                {
                    _1_r = 1 / v_dist((*parts)[j].x, (*parts)[k].x, 1);
                    (*parts)[j].phi_stars -= (*parts)[k].m * _1_r;
                    (*parts)[k].phi_stars -= (*parts)[j].m * _1_r;
                }

                (*parts)[j].energy = (*parts)[j].m * (.5 * (scal_prod((*parts)[j].v, (*parts)[j].v) + (*parts)[j].phi_stars) + (*parts)[j].phi_bgr
                                     - (*parts)[0].m / v_abs((*parts)[j].x));
            }

            if(OUTPUT_INITDATA < 0 || *pcount <= OUTPUT_INITDATA)
                fprintf(get_file(FILE_OUTPUT), "# e_pot=%e\te_kin=%e\n",
                        get_epot(*parts, *pcount, 0) * _conv_e,
                        get_ekin(*parts, *pcount, 0) * _conv_e);

            if(OUTPUT_INITDATA < 0 || *pcount <= OUTPUT_INITDATA)
                fprintf(get_file(FILE_OUTPUT), "#\n## canonical coordinates ##\n");

            // output result
            if(OUTPUT_INITDATA < 0 || *pcount <= OUTPUT_INITDATA)
                for(j = 0; j < *pcount; j++)
                    fprintf(get_file(FILE_OUTPUT), "###### m%d=%e\tm_=%e\teta=%e\t%-1.5e\t%-1.5e\t%-1.5e\t%1.5e\t%1.5e\t%1.5e\n",
                            j,
                            (*parts)[j].m * _conv_m, (*parts)[j].m_ * _conv_m, (*parts)[j].eta * _conv_m,
                            (*parts)[j].x[0] * _conv_x, (*parts)[j].x[1] * _conv_x, (*parts)[j].x[2] * _conv_x,
                            (*parts)[j].v[0] * _conv_v, (*parts)[j].v[1] * _conv_v, (*parts)[j].v[2] * _conv_v);



            return 0;
}

void print_header(int filetype, int pcount)
{
    if(filetype == FILE_DETAIL)
    {
        #ifdef PRINT_DETAIL
        fprintf(get_file(FILE_DETAIL), "#t      \t\tname\tm       \tx       \ty       \tz       \tvx      \tvy      \tvz      \tr       \ta       \t1-e     \tr_min   \tr_max   \tr_apo   \tr_peri   \t");
        #ifdef PRINT_DETAIL_AB
        fprintf(get_file(FILE_DETAIL), "ax      \tay      \taz      \tbx      \tby      \tbz      \t");
        #endif
        fprintf(get_file(FILE_DETAIL), "energy_share\tenergy  \tmultiple\n");
        #endif
    }
}

/*
 *  Output data.
 *
 *  print == -1 : print always
 *            1 : print at next step (see below)
 *            2 : print only at beginning and end
 *  t_old       : last output time in internal units
 *  pt          : negative system energy at beginning (for comparison)
 *  t           : current time in internal units
 *  t_over      : carry-over for total internal time
 *  orbits  > 0 : number of orbits of 1st particle (assuming constant orbital period)
 *          < 0 : negative total time in external units
 *  t_steps > 0 : number of times to output
 *          = 0 : output at all synchronized times
 *          < 0 : output after every -t_steps'th orbit of 1st particle (assuming constant orbital period)
 *  t_max       : internat time until end of calculation
 *  steps       : total number of steps performed (i.e., total number of particles moved)
 *  blocks      : total number of blocks performed (i.e., total number of times at which particles were moved)
 *  orb_count   : number of completed orbits of 1st particle
 */
int output( int print, struct particle parts[], int pcount,
            double t_old,
            double pt,
            double t, double t_over, double orbits, double t_steps, double t_max, double steps, double blocks,
            int orb_count)
{
    double e, t_tot=t+t_over, stepsc=.0, stepsp=.0, sigma_a=.0, sigma_e=.0, sigma_j=.0, e_max=.0, a_min=.0, sigma_v=.0;
    double vr_2=.0, v_2=.0, vr[3];
    struct particle *p;
    int ellipticals=0, stepscount=0;
    double kep_a, kep_e[3], kep_j[3], kep_omega;
    static int t_real_last = -1;

    // output if necessary
    if((print == 1
        && (t_steps == 0
        || (t_steps > 0 && floor(t_tot*t_steps/t_max) > floor(t_old*t_steps/t_max))
        || (t_steps < 0 && (orb_count % (int)(-t_steps) == 0)))
        )
        || (print == 2
        && floor(t_tot*10.0/t_max) > floor((t_tot-parts[1].dt)*10.0/t_max)
        && (t_tot==0 || t_tot>parts[1].dt))
        || print == -1)
    {
        for(p = parts + 1; p < parts + pcount; p++)
            if(fabs(t - p->t) > DT_TOLERANCE)
            {
                //printf("%e\t\t", p->t);

                return -1;
            }

        if(t_tot <= .0)
        {
            return 0;
        }

        //printf("#[t=%e]\n", convert_time(t_tot, 0));
        predict(parts, pcount, t);
        assert(P_OUTPUT > 0 && P_OUTPUT < pcount);

        for(p = parts + 1; p < parts + pcount; p++)
        {
            get_constants(p->x, p->v, parts->m + M_ENCL_IO(v_abs(p->x)), kep_j, kep_e, &kep_a, &kep_omega);
            p->curr_t = 2. * M_PI / kep_omega;
            p->curr_j = v_abs(kep_j);
            p->curr_e = v_abs(kep_e);
            p->curr_a = kep_a;

            if(p->t > p->io_tlast)
            {
                stepsc += p->curr_t * p->io_steps_c / (p->t - p->io_tlast);
                stepsp += p->curr_t * p->io_steps_p / (p->t - p->io_tlast);
                stepscount++;
                p->io_steps_c = 0;
                p->io_steps_p = 0;
                p->io_tlast = p->t;
            }
            if(p->curr_e < 1.)
            // particle still on elliptical orbit
            {
                ellipticals++;
                sigma_a += (p->curr_a / p->orig_a - 1) * (p->curr_a / p->orig_a - 1);
                sigma_e += (p->curr_e - p->orig_e)     * (p->curr_e - p->orig_e);
                sigma_j += (p->curr_j / p->orig_j - 1) * (p->curr_j / p->orig_j - 1);
                sigma_v += scal_prod(p->v, p->v);
                vr[0] = vr[1] = vr[2] = scal_prod(p->xp, p->vp) / scal_prod(p->xp, p->xp);
                vr[0] *= p->xp[0]; vr[1] *= p->xp[1]; vr[2] *= p->xp[2];
                vr_2 += scal_prod(vr, vr);
                v_2  += (p->vp[0] - vr[0]) * (p->vp[0] - vr[0])
                + (p->vp[1] - vr[1]) * (p->vp[1] - vr[1])
                + (p->vp[2] - vr[2]) * (p->vp[2] - vr[2]);

                if(p->curr_e > e_max)
                    e_max = p->curr_e;
                if(p == parts + 1 || p->curr_a < a_min)
                    a_min = p->curr_a;
            }
        }

        sigma_a = sqrt(sigma_a / (ellipticals + 1));
        sigma_e = sqrt(sigma_e / (ellipticals + 1));
        sigma_j = sqrt(sigma_j / (ellipticals + 1));
        sigma_v = sqrt(sigma_v / ellipticals);

        e = .0;
        //e = get_epot(parts, pcount, 1) + get_ekin(parts, pcount, 1);
        for(p = parts + 1; p < parts + pcount; p++)
        {
            p->energy = get_energy(parts, pcount, p - parts);
            e += p->energy;
        }

        e += lost_energy;
        stepsc /= stepscount;
        stepsp /= stepscount;

        fprintf(get_file(FILE_OUTPUT),
                " %e\t%d\t \t%e\t%e\t%e\t%e\t \t%e\t%e\t%e\t%e\t%e\t \t%e\t%e\t%3.6f\t%3.6f\t \t%e\t%e\t%e\t%e\t%e\n",
                convert_time(t_tot, 0),                                  // current time
                pcount,                                                  // number of particles

                (e - e_last)/fabs(pt) / (t_tot - t_last) / _conv_t,      // current relative energy error per time unit
                (e + pt)/fabs(pt),                                       // relative energy error
                -lost_energy / pt,                                       // fraction of energy removed from system
                emitted_energy / C_C / C_C * _conv_m,                    // emitted energy in c² M_sun

                convert_time(t_tot - t_last, 0),                         // model time since last output
                stepsc,                                                  // average number of steps since last output
                convert_time((t_tot - t_last)/(steps-steps_last)*(double)pcount, 0), // average time-step
                convert_time((t_tot - t_last)/(blocks-blocks_last), 0),  // average block-step
                t_real_last > 0 ? convert_time((t_tot - t_last), 0)
                / (double)(time(NULL) - t_real_last) : .0/.0,            // progress in model years per wall clock second

                e_max,                                                   // maximum eccentricity among all particles
                a_min * _conv_x,                                         // minimum semi-major axis among all particles
                inclination(parts[0].xp, parts[0].vp,
                            parts[1].xp, parts[1].vp,
                            parts[2].xp, parts[2].vp),                   // inclination between first two orbits in list
                vec_angle(parts[1].xp, parts[0].xp, parts[2].xp),        // angle between first two particles' position vectors

                                                                         // properties of output particle:
                parts[P_OUTPUT].curr_e,                                  // eccentricity
                convert_length(parts[P_OUTPUT].curr_a, 0),               // semi-major axis
                v_dist(parts[0].x, parts[P_OUTPUT].x, 1)*_conv_x,        // central distance
                parts[P_OUTPUT].energy
                + .5 * parts[P_OUTPUT].m * parts[P_OUTPUT].phi_stars,    // energy
                convert_time(t_tot - t_old, 0)                           // time-step
                );

        #ifdef PRINT_DETAIL
        if(!(detailout--))
        {
            #ifdef PRINT_DETAIL_AB
            double j_[3], ecc_[3], ecc, a, b, omega, a_[3] = {.0,.0,.0}, b_[3] = {.0,.0,.0};
            int i;
            //fprintf(get_file(FILE_DETAIL), "\t\t%1.12e", convert_time(t_tot, 0));
            #endif

            for(p = parts; p < parts + pcount; p++)
            if(N_MAX_DETAIL < -1 || p->name <= N_MAX_DETAIL)
            {
                #ifdef PRINT_DETAIL_AB
                if(p > parts)
                {
                  get_constants(p->x, p->v, parts->m, j_, ecc_, &a, &omega);
                  ecc = v_abs(ecc_);
                  vec_prod(j_, ecc_, b_);
                  b = a * sqrt(fabs(1-ecc*ecc)) / v_abs(b_);
                  for(i = 0; i < 3; i++)
                    {
                      a_[i]  = a*ecc_[i]/ecc;            // semi major vector
                      b_[i] *= b;                        // semi minor vector
                    }
                }
                #endif
                fprintf(get_file(FILE_DETAIL),
                        #ifdef PRINT_DETAIL_AB
                        "%1.12e\t %d \t%1.5e\t%-1.6e\t%-1.6e\t%-1.6e\t%-1.6e\t%-1.6e\t%-1.6e\t%1.4e\t%1.4e\t%1.4e\t%1.6e\t%1.6e\t%1.6e\t%1.6e\t%-1.6e\t%-1.6e\t%-1.6e\t%-1.6e\t%-1.6e\t%-1.6e\t%1.6e\t%1.6e\t%d\n",
                        #else
                        "%1.12e\t %d \t%1.5e\t%-1.6e\t%-1.6e\t%-1.6e\t%-1.6e\t%-1.6e\t%-1.6e\t%1.4e\t%1.4e\t%1.4e\t%1.6e\t%1.6e\t%1.6e\t%1.6e\t%1.6e\t%1.6e\t%d\n",
                        #endif
                        convert_time(t_tot, 0),
                        p->name,
                        p->m * _conv_m,
                        p->x[0]*_conv_x, p->x[1]*_conv_x, p->x[2]*_conv_x,
                        p->v[0]*_conv_v, p->v[1]*_conv_v, p->v[2]*_conv_v,
                        v_abs(p->x)*_conv_x,
                        p->curr_a*_conv_x,
                        1 - p->curr_e,
                        p->rmin*_conv_x,
                        p->rmax*_conv_x,
                        p->r_apo*_conv_x,
                        p->r_peri*_conv_x
                        #ifdef PRINT_DETAIL_AB
                        , a_[0]*_conv_x, a_[1]*_conv_x, a_[2]*_conv_x
                        , b_[0]*_conv_x, b_[1]*_conv_x, b_[2]*_conv_x
                        #endif
                        , p->energy
                        , p->energy + p->m * (.5 * p->phi_stars /*+ p->phi_bgr*/)
                        , 1
                        );

                p->rmin = 1.e99;
                p->rmax = .0;
            }
            //fprintf(get_file(FILE_DETAIL), "\n");
            fflush(get_file(FILE_DETAIL));
            detailout += PRINT_DETAIL_INTERVAL;
        }
        #endif
        fflush(get_file(FILE_OUTPUT));
        fflush(get_file(FILE_WARNING));
        t_last = t_tot;
        e_last = e;
        blocks_last = blocks;
        steps_last  = steps;
        t_real_last = time(NULL);

        #ifdef DUMP_BEFORE_ERR
        dump_no = 0;
        static int new = 1;
        if(!new && (e + pt) /fabs(pt)> 1.e-4)
        {
            exit(0);
        }
        new = 0;
        #endif
    }

    return 0;
}

/*
 *  Calculate orbital parameters given central mass as well as particle position and velocity.
 *  Arguments passed from command line as {m0 rx ry rz vx vy vz}
 */
void calc(char **argv)
{
    int i;
    double r_[3], v_[3], j_[3], e_[3], omega, e, j, a_[3], b_[3], a, b, m=atof(argv[0])/_conv_m;

    for(i = 0; i < 3; i++)
    {
        r_[i] = atof(argv[i+1]) / _conv_x;
        v_[i] = atof(argv[i+4]) / _conv_v;
    }
    get_constants(r_, v_, m, j_, e_, &a, &omega);
    j = v_abs(j_); e = v_abs(e_);
    printf("x=(%e, %e, %e), |x|=%e\n", r_[0] * _conv_x, r_[1] * _conv_x, r_[2] * _conv_x, v_abs(r_) * _conv_x);
    printf("v=(%e, %e, %e), |v|=%e\n", v_[0] * _conv_v, v_[1] * _conv_v, v_[2] * _conv_v, v_abs(v_) * _conv_v);
    //printf("j=(%e, %e, %e), |j|=%e\n", j_[0], j_[1], j_[2], j);
    printf("e=(%e, %e, %e), |e|=%1.22e\n", e_[0], e_[1], e_[2], e);
    //printf("p=(%e, %e, %e)\n", e_[0]/e*j*j/(1-e),  e_[1]/e*j*j/(1-e),  e_[2]/e*j*j/(1-e));
    printf("peri=(%e, %e, %e)\n", a * (1-e) * e_[0]/e * _conv_x, a * (1-e) * e_[1]/e * _conv_x, a * (1-e) * e_[2]/e * _conv_x);
    printf("apo=(%e, %e, %e)\n", a * (-e-1) * e_[0]/e * _conv_x, a * (-e-1) * e_[1]/e * _conv_x, a * (-e-1) * e_[2]/e * _conv_x);
    //printf("center=(%e, %e, %e)\n", -a * e_[0] * _conv_x, -a * e_[1] * _conv_x, -a * e_[2] * _conv_x);
    //-e[0]/e_*j_*j_/(1+e_) * _conv_x,  -e[1]/e_*j_*j_/(1+e_) * _conv_x,  -e[2]/e_*j_*j_/(1+e_) * _conv_x);

    vec_prod(j_, e_, b_);
    b = a * sqrt(fabs(1-e*e)) / v_abs(b_);
    for(i = 0; i < 3; i++)
    {
        a_[i]  = a*e_[i]/e;            // semi major vector
        b_[i] *= b;                        // semi minor vector
    }

    printf("T=%e\na=%2.18f\t\n", 2*M_PI/omega * _conv_t, a*_conv_x);

  /*
    printf("a=(%e, %e, %e), |a|=%e\n", a_[0] * _conv_x, a_[1] * _conv_x, a_[2] * _conv_x, v_abs(a_) * _conv_x);
    printf("b=(%e, %e, %e), |b|=%e\n", b_[0] * _conv_x, b_[1] * _conv_x, b_[2] * _conv_x, v_abs(b_) * _conv_x);
    printf("#	, %e*cos(t)+%e*sin(t)+%e, %e*cos(t)+%e*sin(t)+%e notitle with lines 2\\\n",
    a_[0] * _conv_x, b_[0] * _conv_x, -a * e_[0] * _conv_x,
    a_[1] * _conv_x, b_[1] * _conv_x, -a * e_[1] * _conv_x);
    printf("#	, %e*cos(t)+%e*sin(t)+%e, %e*cos(t)+%e*sin(t)+%e notitle with lines 2\\\n",
    a_[0] * _conv_x, b_[0] * _conv_x, -a * e_[0] * _conv_x,
    a_[2] * _conv_x, b_[2] * _conv_x, -a * e_[2] * _conv_x);
    printf("#	, .001*cos(t)+%e, .001*sin(t)+%e notitle with lines 2\\\n",
    r_[0] * _conv_x, r_[1] * _conv_x);
    printf("#	, .001*cos(t)+%e, .001*sin(t)+%e notitle with lines 2\\\n",
    r_[0] * _conv_x, r_[2] * _conv_x);
  */
}

int smaller(const void *a, const void *b)
{
    double x = *((double *)a), y = *((double *)b);
    return  x == y ? 0 : (x < y ? -1 : 1);
}

/*
 * Create random distribution of particles.
 * Parameters:
 *              n         number of particles
 *              m_0       mass of central particle
 *              p         3D power index of mass distribution
 *              a_min     minimum semimajor axis
 *              a_max     maximum semimajor axis
 *              e_min     minimum eccentricity
 *              e_max     maximum eccentricity
 *              m_min     minimum particle mass
 *              m_max     maximum particle mass
 *              imf_slope
 *              sigma     deviation angle from disc (if given and not equal to "_", disc in xy-plane; spherically symmetric otherwise)
 *              delta     final rotation of system around x-axis (i.e. inclination to xy-plane)
 *              m1        starting interval of 2nd IMF section
 *              alpha1    slope of 2nd IMF section
 *              m2/alpha2...
 */

void create(int argc, char **arg)
{

    if(argc < 9 || (argc > 12 && argc % 1))
    {
        fprintf(stderr, "too few arguments! Usage: nbody ! <N> <m_0> <alpha_3D> <a_min> <a_max> <e_min> <e_max> <m_min> <m_max> [imf_slope] [sigma] [delta] {m_i alpha_i}*\n");
        exit(-1);
    }

    int i, k, n=atoi(arg[0]);
    double m0=atof(arg[1]) / _conv_m, p=atof(arg[2])+2., mi;
    double a_min=atof(arg[3]) / _conv_x, a_max=atof(arg[4]) / _conv_x, e_min=atof(arg[5]), e_max=atof(arg[6]);
    double m_min=atof(arg[7]), m_max=atof(arg[8]);
    double *a, e, phi, int_=pow(a_max,p+1)-pow(a_min,p+1), x[3], v[3], x_, v_, g, h;
    double j[3], e_[3], a_, om, _a[3], _n[3], _b[3], _j[3], e1[3], e2[3];

    // for power-law mass spectrum
    int imf_n = 1 + (argc > 12 ? (argc - 12) / 2 : 0); // number of power-law sections
    double *imf_m     = (double *)malloc((imf_n+1) * sizeof(double));
    double *imf_alpha = (double *)malloc(imf_n     * sizeof(double));
    double *imf_r     = (double *)malloc(imf_n     * sizeof(double));
    double *imf_f     = (double *)malloc((imf_n+1) * sizeof(double));

    imf_m[0] = m_min;
    imf_m[imf_n] = m_max;
    imf_alpha[0] = (argc > 9 ? -atof(arg[9]) : .0);
    for(i = 1; i < imf_n; i++)
    {
        imf_m[i]     =  atof(arg[2*i+10]);
        imf_alpha[i] = -atof(arg[2*i+11]);
    }

    setup_imf(imf_m, imf_alpha, imf_n, imf_r, imf_f);

    for(i = 0; i < imf_n; i++)
    {
        printf("#  f_%d(m) = %1.8f * m ^ %1.1f, F(%1.3f) = %1.8f, F(%1.3f) = %1.8f\n",
                i, imf_r[i], imf_alpha[i], imf_m[i], imf_f[i], imf_m[i+1], imf_f[i+1]);
    }

    printf("#\tm\t\tx\t\ty\t\tz\t\tvx\t\tvy\t\tvz\n");
    printf("# central particle\n\t%e\t\t.0\t\t.0\t\t.0\t\t.0\t\t.0\t\t.0\n", m0 * _conv_m);
    seedrand();

    // create list of semi-major axes and sort for output
    a = (double *)malloc(n*sizeof(double));
    for(i = 1; i < n; i++)
    {
        a[i] = pow(drand(.0, int_) + pow(a_min, p+1), 1 / (p + 1));
    }
    qsort(a+1, n-1, sizeof(double), smaller);

    // choose remaining parameters
    for(i = 1; i < n; i++)
    {
        e  = sqrt(drand(e_min*e_min, e_max*e_max));
        //e = square(drand(e_min, e_max));
        //do { e = .1 + gauss() * .1; } while (e < .0 || e >= 1.);
        mi = get_imf_mass(imf_m, imf_alpha, imf_n, imf_r, imf_f) / _conv_m;

        if(drand(.0,1.) > (1. - BINARY_FRACTION))
        //if(i < .5 * n)
        {
            printf("#binary");
            mi = -2.*mi;
        }
        //*/

        /*
            # calculate mean mass
            alpha:=1.35:
            a:=0.08:b:=120:
            p:=m->m^(-alpha);
            m_:=int(m*p(m),m=a..b)/int(p(m),m=a..b);
        */

        double r = a[i];
        //a[i] /= 1. + .5*e*e;
        r = a[i] * (1. + .5*e*e);
        x_ = (1 + e) * a[i];
        v_ = sqrt((m0 + M_ENCL_IO(r) + fabs(mi)) * (1 - e) / a[i] / (1 + e));

        if(argc > 10 && arg[10][0] != '_')     // for discs: gaussian distribution of inclination vector
        {
            double sigma = atof(arg[10]);       // standard variation of angle of normal vector to z-axis in degree
            if(sigma >= .0)
            {
                do
                {
                    h = gauss() * sigma;
                } while(fabs(h) > 90.);
                h = cos(h/180.*M_PI);
            }

            else // uniform distribution within angle
            {
                h = drand(cos(sigma/180*M_PI), 1.);
                if(drand(-1., 1.) < 0)
                    h = -h;
            }
        }
        else                             // spherically symmetric distribution
        {
            h = drand(-1., 1.);
        }
        g   = sqrt(1. - h * h);
        phi = drand(.0, 2. * M_PI);
        _n[0] = g * cos(phi);
        _n[1] = g * sin(phi);
        _n[2] = h;

        if(fabs(h) >= 1.)
        {
            e1[0] = e2[1] = 1.;
            e1[1] = e1[2] = e2[0] = e2[2] = .0;
        }
        else
        {
            // assuming |z| < 1.
            // orthogonal unit vectors: e1 = (-y, x, 0) / g, e2 = ea x e1 = (-xz, -yz, g^2) / g
            // g^2 = x^2 + y^2
            e1[0] = -_n[1] / g;         e1[1] = _n[0] / g;          e1[2] = .0;
            e2[0] = -_n[0] * _n[2] / g; e2[1] = -_n[1] * _n[2] / g; e2[2] = g;
        }

        // set ea = e1 * sin(phi) + e2 * cos(phi), eb = e1 * cos(phi) - e2 * sin(phi)
        phi = drand(-M_PI, M_PI);

        _a[0] = a[i] * (e1[0] * sin(phi) + e2[0]* cos(phi));
        _a[1] = a[i] * (e1[1] * sin(phi) + e2[1]* cos(phi));
        _a[2] = a[i] * (e1[2] * sin(phi) + e2[2]* cos(phi));
        _b[0] = a[i] * sqrt(fabs(1 - e*e)) * (e1[0] * cos(phi) - e2[0] * sin(phi));
        _b[1] = a[i] * sqrt(fabs(1 - e*e)) * (e1[1] * cos(phi) - e2[1] * sin(phi));
        _b[2] = a[i] * sqrt(fabs(1 - e*e)) * (e1[2] * cos(phi) - e2[2] * sin(phi));

        phi = solve_kepler(drand(.0, 2 * M_PI), e);
        vec_prod(_a, _b, _j);
        for(k = 0; k < DIMENSIONS; k++)
        {
            x[k] =  _a[k] * cos(phi) + _b[k] * sin(phi) - _a[k]*e;  // location
            v[k] = -_a[k] * sin(phi) + _b[k] * cos(phi);            // direction of v only
            _j[k] *= x_ * v_ / a[i] / a[i] / sqrt(fabs(1 - e*e));
        }
        v_ = v_abs(v);
        v_ = v_abs(_j) / v_abs(x) / sin(acos(scal_prod(x, v)/ (v_abs(x)*v_))) / v_;

        for(k = 0; k < DIMENSIONS; k++)
        {
            v[k] *= v_;
        }
        //printf("#\t%e\t%e\t%e\n", a, e, scal_prod(x, v));
        get_constants(x, v, m0 + M_ENCL_IO(r) + fabs(mi), j, e_, &a_, &om);
        printf("# particle no. %d, a=%1.8e, e=%1.8f\t<r>=%1.8e M_encl=%1.8e\tr=%1.8e=%1.8e\n",
                i, a_ * _conv_x, v_abs(e_), r*_conv_x, M_ENCL_IO(r)*_conv_m,
                _conv_x * sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]),
                _conv_x * a[i] * (1. - e*cos(phi)));
        // printf("#\t%e\t%e\t%e\n", a, v_abs(e_), v_abs(j));

        double delta = .0;
        if(argc > 11 && arg[10][0] != '_' && arg[11][0] != '_')
            delta = atof(arg[11]) * M_PI / 180.;
        printf("\t%1.12e \t%1.12e \t%1.12e \t%1.12e \t%1.12e \t%1.12e \t%1.12e\n",
        mi * _conv_m,
        x[0] * _conv_x, (x[1]*cos(delta) - x[2]*sin(delta)) * _conv_x, (x[1]*sin(delta) + x[2]*cos(delta)) * _conv_x,
        v[0] * _conv_v, (v[1]*cos(delta) - v[2]*sin(delta)) * _conv_v, (v[1]*sin(delta) + v[2]*cos(delta)) * _conv_v);
    }
    free(a);

}

/*
 *  Convert input file from external units to N-body units, where
 *     G = 1
 *     M_tot = 1
 *     center of mass = total linear momentum = 0
 *     potential energy = -0.5
 *
 *  conversion factors are printed as comments
 */
void convert_to_nbody_units(struct particle parts[], int pcount)
{

    int i;
    double c_m=.0, c_r=.0, c_v=.0, e_pot=.0, c_v_[3]={0,0,0};
    struct particle *p, *p2;

    for(p = parts; p < parts + pcount; p++)
        c_m += p->m;
    for(p = parts; p < parts + pcount; p++)
        p->m /= c_m;
    for(p = parts; p < parts + pcount - 1; p++)
        for(p2 = p + 1; p2 < parts + pcount; p2++)
            e_pot -= p->m * p2->m / v_dist(p->x, p2->x, 1);

    c_r = -.5 / e_pot;
    c_v = sqrt(c_m / c_r);

    for(p = parts; p < parts + pcount; p++)
        for(i = 0; i < 3; i++)
        {
            p->x[i] /= c_r;
            p->v[i] /= c_v;
            c_v_[i] += p->v[i]*p->m;
        }
        for(p = parts; p < parts + pcount; p++)
            for(i = 0; i < 3; i++)
                p->v[i] -= c_v_[i];

        printf("###\tE_kin=%e\tE_pot=%e\tE=%e\n",
                get_ekin(parts, pcount, 0),
                get_epot(parts, pcount, 0),
                get_ekin(parts, pcount, 0) + get_epot(parts, pcount, 0));
        printf("### 1 nbody unit = %e pc\n", _conv_x * c_r);
        printf("### 1 nbody unit = %e M_sun\n", _conv_m * c_m);
        printf("### 1 nbody unit = %e yr\n", _conv_t * sqrt(c_r*c_r*c_r/c_m));
        printf("### 1 yr = %e nbody units\n", 1. / _conv_t / sqrt(c_r*c_r*c_r/c_m));

        for(p = parts; p < parts + pcount; p++)
            printf("%18.10f\t%18.10f\t%18.10f\t%18.10f\t%18.10f\t%18.10f\t%18.10f\n",
                    p->m,
                    p->x[0], p->x[1], p->x[2],
                    p->v[0], p->v[1], p->v[2]);

}

/*
 *  Read input file from standard input and add orbital parameters as comments
 */
void comment_datfile()
{
    double m_cent = .0, _1_m_cent;
    double m, x[3], v[3], j[3], e[3], r, a, ecc2;
    char s[300];
    int i;

    for(i=0; !feof(stdin);)
    {
        if(fgets(s, 299, stdin) == NULL) continue;
        if(s[0] != '#' && s[0] != '\n')
        {
            sscanf(s, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
                    &m, x, x+1, x+2, v, v+1, v+2);
                    m /= _conv_m;
                    x[0] /= _conv_x;
                    x[1] /= _conv_x;
                    x[2] /= _conv_x;
                    v[0] /= _conv_v;
                    v[1] /= _conv_v;
                    v[2] /= _conv_v;
                    if(i == 0)
                         m_cent = m;
                    else
                    {
                        vec_prod(x, v, j);
                        vec_prod (v, j, e);
                        r = sqrt(scal_prod(x, x));
                        _1_m_cent = 1. / (m_cent + m + M_ENCL_IO(r));
                        e[0] = e[0] * _1_m_cent - x[0]/r;
                        e[1] = e[1] * _1_m_cent - x[1]/r;
                        e[2] = e[2] * _1_m_cent - x[2]/r;
                        ecc2 = scal_prod(e, e);
                        a = scal_prod(j, j) * _1_m_cent / fabs(1 - ecc2);
                        printf("# T= %2.12f\ta= %2.15f\te= %2.12f \tM_encl=%e ##### \n",
                           convert_time(2 * M_PI * sqrt(pow(a, 3.0) * _1_m_cent), 0),
                           convert_length(a, 0),
                           sqrt(ecc2),
                           M_ENCL_IO(r) * _conv_m);
                    }
                    i++;
        }

        printf("%s", s);
    }

}



/*
 *  Print conversion factors
 */
void print_conv()
{
  printf("1 mass unit \t= %e M_sun\n1 length unit \t= %e pc\n1 time unit \t= %e yr\n1 vel. unit \t= %e pc/yr \t\t= %e km/s\n1 energy unit \t= %e G*M_sun^2/pc \t= %e M_sun*pc^2/yr^2\t= %e erg\n",
        _CONV_M, _CONV_X, _CONV_T, _CONV_V, 977813.952*_CONV_V, _CONV_E, _CONV_M*_CONV_V*_CONV_V, 8.55470098e40 * _CONV_E);
}
