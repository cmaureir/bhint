/******************************************************************************/
/*                                                                            */
/*                     CONFIGURATION PARAMETERS FOR BHINT                     */
/*                                                                            */
/******************************************************************************/

/*          ****** PRECISION PARAMETERS ******          */
#define MIN_EVALS     80           // precision parameter: minimum number of evaluations per orbit
#define ETA           0.1          // precision parameter: \eta from Aarseth (1985)


/*          ******* USE PN TREATMENT *********          */
//#define PN                         // Use post-Newtonian treatment
                                   // ********** CONFIGURE BELOW **********


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

#ifdef PN

/******************************************************************************
 *
 * Parameters for use of post-Newtonian terms.
 *
 * The SWITCH-parameters specify when (and for which particles) to switch PN treatment on and off.
 * Choose or combine nay of the examples.
 * Default: PN for all particles, but only if estimated inspiral time (Peters 1964) is below 100 Myr; switch off again if it exceeds 1 Gyr.
 * Make sure PN is not constantly switched for any particle (e.g., every revolution or so).
 *
 ******************************************************************************/

#define PN_KEPLER_FACT .00390625   // Requires a reciprocal integer!

/*
  Peters 1964, Eq. 5.14: T = 12/19 * c0^4 / beta * INT(...) = T_c * 48/19 * (c0/a0)^4 * INT(...)
  T / 4T_c * (1-e^2)^(-7/2) = 12/19 * (c0/a0)^4 * (1-e^2)^(-7/2) *INT(...) =: g(e) \in {0.2485..0.44} => g(e) ~= 0.3
  => T = 4 T_c * (1-e^2)^(7/2) * g(e)

  e_:=0.2: c0_a0:=(1-e_^2)/e_^(12/19)/(1+121/304*e_^2)^(870/2299):
  T_Tc:=12/19*c0_a0^4*int(e^(29/19)*(1+121/304*e^2)^(1181/2299)/(1-e^2)^(3/2),e=0..e_)*4; # Fig. 2
  evalf(T_Tc/4/(1-e_^2)^(7/2)); # =g(e), 0.24 < g(e) < 0.46
*/

#define T_INSPIRAL (square(square(p->curr_a)) / (__1_c5*.2*64.*parts[0].m*(parts[0].m+p->m)*p->m) * exp(3.5 * log(1. - square(p->curr_e))) *.3)


/***** PN ONLY IF INSPIRAL TIME FALLS BELOW 100 MYR (switch off again if it exceeds 1 Gyr) *****/
#define SWITCHON_PN  (T_INSPIRAL < convert_time(1.e8, 1))
#define SWITCHOFF_PN (T_INSPIRAL > convert_time(1.e10, 1))


/***** PN FOR ECCENTRIC PARTICLES ONLY: switch on for e>0.9, switch off when e falls below 0.8 *****/
//#define SWITCHON_PN  (p->curr_e > .9 && (N_MAX_DETAIL <  -1 || p->name <= N_MAX_DETAIL))
//#define SWITCHOFF_PN (p->curr_e < .8 || (N_MAX_DETAIL >= -1 && p->name >  N_MAX_DETAIL))


/***** PN FOR DISPLAY PARTICLES ONLY *****/
//#define SWITCHON_PN  (N_MAX_DETAIL <  -1 || p->name <= N_MAX_DETAIL)
//#define SWITCHOFF_PN (N_MAX_DETAIL >= -1 && p->name >  N_MAX_DETAIL)


/***** PN FOR PARTICLE NAME 1 ONLY *****/
//#define SWITCHON_PN  (p->name == 1)
//#define SWITCHOFF_PN (p->name != 1)


/***** PN ALWAYS ON *****/
//#define SWITCHON_PN  1
//#define SWITCHOFF_PN 0

/******************************************************************************/

#endif
