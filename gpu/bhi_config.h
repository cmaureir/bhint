/******************************************************************************/
/*                                                                            */
/*                     CONFIGURATION PARAMETERS FOR BHINT                     */
/*                                                                            */
/******************************************************************************/

/*          ****** PRECISION PARAMETERS ******          */
#define MIN_EVALS     80           // precision parameter: minimum number of evaluations per orbit
#define ETA           0.1          // precision parameter: \eta from Aarseth (1985)


/*          ************ OUTPUT **************          */
#define PRINT_DETAIL               // print particle details for all particles (up to N_MAX_DETAIL), one line each
#define PRINT_DETAIL_INTERVAL 100  // print details of all particles for every Nth output line

#define P_IMBH        0            // Name of special particle (IMBH) for tracking; preferably 1, if present

#define N_MAX_DETAIL  -2           // maximum particle name to print details for, -1 for none, -2 for all

#if (defined(P_IMBH) && P_IMBH>0)  // IMBH name, to be defined in ul_nbody.h
#define P_OUTPUT P_IMBH            // name of particle to print details (a, e) for in main output file
#else                              // usually IMBH, 1st or any other special particle
#define P_OUTPUT 1
#endif

/*          ************ VARIOUS *************          */
#define DUMP_INTERVAL    43200     // interval for creating dump files (in seconds of real time)

#define N_MAX            150000    // maximum number of particles

#define MAX_X            10.       // Maximum central distance in pc (particle will be removed if further out)

#define BINARY_FRACTION  .0        // fraction of stars marked as binaries (i.e., with double initial mass)
                                   // for model creation only


/******************************************************************************
 *
 * Defines the density profile of the model (used for optimisation and enclosed mass).
 *
 * Following parameters represent the findings of Schoedel et al. (A&A 469, 125 (2007), Eq. 7),
 * the model includes stars up to ME_RMAX = 0.22 pc.
 *
 ******************************************************************************/

#define ME_R0            .22       // break radius in external units [pc]; assumes R_MIN = 0
#define ME_RMAX          .22       // outer radius of model in external units [pc]; may be 0 for no cusp
#define ME_GAMMA1        -1.2      // power-law exponent inside break radius
#define ME_GAMMA2        -1.75     // power-law exponent outside break radius

                                   // enclosed mass in stars (from model), to estimate orbital parameters
#define M_ENCL_(r)       convert_mass(((r < convert_length(ME_R0, 1) || ME_RMAX < ME_R0) \
				       ? 2.1e5 * pow(min(convert_length(r, 0), ME_RMAX) / ME_R0, 3. + ME_GAMMA1) \
				       : 3.e5 * (pow(min(convert_length(r, 0), ME_RMAX) / ME_R0, 3. + ME_GAMMA2) - .3)), 1)


#define M_ENCL(r)        M_ENCL_(r)

/******************************************************************************/
