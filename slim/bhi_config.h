/******************************************************************************/
/*                                                                            */
/*                     CONFIGURATION PARAMETERS FOR BHINT                     */
/*                                                                            */
/******************************************************************************/

/*          ****** PRECISION PARAMETERS ******          */
#define MIN_EVALS     80           // precision parameter: minimum number of evaluations per orbit
#define ETA           0.1          // precision parameter: \eta from Aarseth (1985)


/*          ******* EXTERNAL POTENTIAL *******          */
//#define EXT_POT                    // Use external spherical potential defined below for integration.
                                   // Requires density, enclosed mass, and potential at radius r
                                   // (EP_RHO(r), EP_M(r), EP_PHI(r), all in internal units)
                                   // Use of USE_M_ENCL recommended.
                                   // ********** CONFIGURE BELOW **********


/*          ********* MASS PROFILE ***********          */
//#define USE_M_ENCL                 // Use mass profile defined below (e.g., for stellar cusp or external potential)
                                   // to calculate orbital parameters.
                                   // (for creating models and output only; does not affect simulation)
                                   // ********** CONFIGURE BELOW **********

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


/******************************************************************************/

#ifdef EXT_POT

/******************************************************************************
 *
 * Defines external spherical potential for integration.
 * Requires density, enclosed mass, and potential at radius r (EP_RHO(r), EP_M(r), EP_PHI(r), all in internal units).
 * For convenience, we provide properties for a broken power-law density profile in external units:
 *   rho(r) = EP_RHO0 * (r / EP_R0)^EP_GAMMA1  ,  r <= EP_R0
 *   rho(r) = EP_RHO0 * (r / EP_R0)^EP_GAMMA2  ,  r >= EP_R0
 * restricted to EP_RMIN <= r <= EP_RMAX.
 *
 * Following parameters represent the findings of Schoedel et al. (A&A 469, 125 (2007), Eq. 7),
 * the external potential only includes the region outside EP_RMIN=0.22 pc
 * (assuming the model includes the stars inside 0.22 pc, see also below).
 *
 ******************************************************************************/

#define EP_RHO0   2.8e6            // density at break radius in external units [M_sun pc^-3]
#define EP_RMIN   .22              // inner radius of external profile in external units [pc]
#define EP_R0     .22              // break radius in external units [pc]
#define EP_RMAX   10.              // outer radius of external profile in external units [pc]; requires EP_RMAX >= EP_RMIN
#define EP_GAMMA1 -1.2             // power-law exponent inside break radius
#define EP_GAMMA2 -1.75            // power-law exponent outside break radius

#define EP_MIN(a, b)             (a < b ? a : b)
#define EP_MAX(a, b)             (a > b ? a : b)
#define EP_MIN3(a, b, c)         EP_MIN(EP_MIN(a, b), c)
#define EP_MAX3(a, b, c)         EP_MAX(EP_MAX(a, b), c)

                                   // density in external units
#define EP_RHO_PHYS(r)           ((r < EP_RMIN || r > EP_RMAX) ? .0 : EP_RHO0 * pow(r/EP_R0, r < EP_R0 ? EP_GAMMA1 : EP_GAMMA2))


                                    // mass integral contribution of interval r1..r2, power-law exponent g
#define EP_M_INT(r1, r2, g)      (4.*M_PI*EP_RHO0*pow(EP_R0, -g)/(3.+g) * (pow(r2, 3. + g) - pow(r1, 3. + g)))

                                  // potential integral contribution of interval a..b, with r>=b, power-law exponent g
#define EP_PHI_INT1(r, g, a, b)  (EP_RHO0 / ((3. + g) * r * pow(EP_R0, g)) * (pow(b, 3. + g) - pow(a, 3. + g)))
                                   // potential integral contribution of interval a..b, with r<=a, power-law exponent g
#define EP_PHI_INT2(g, a, b)     (EP_RHO0 / ((2. + g)     * pow(EP_R0, g)) * (pow(b, 2. + g) - pow(a, 2. + g)))

                                   // enclosed mass at r in external units
                                   // first interval:  r_min           ... max(r_min, min(r_0, r, r_max))
                                   // second interval: max(r_0, r_min) ... max(r_min, r_0, min(r, r_max))

#define EP_M_PHYS(r)             (EP_M_INT(EP_RMIN,                EP_MAX (EP_RMIN, EP_MIN3(EP_R0, r, EP_RMAX)), EP_GAMMA1) + \
				  EP_M_INT(EP_MAX(EP_R0, EP_RMIN), EP_MAX3(EP_RMIN, EP_R0, EP_MIN (r, EP_RMAX)), EP_GAMMA2))

                                   // potential at r in external units
                                   // first interval inside r:   r_min              ... max(r_min, min(r, r_0, r_max))
                                   // first interval outside r:  max(r_min, r)      ... max(r_min, r, min(r_0, r_max))
                                   // second interval inside r:  max(r_0, r_min)    ... max(r_0, r_min, min(r, r_max))
                                   // second interval outside r: max(r_0, r_min, r) ... max(r_0,            r, r_max )

#define EP_PHI_PHYS(r)           (-4. * M_PI *				\
				  (EP_PHI_INT1(r, EP_GAMMA1, EP_RMIN,                    EP_MAX (EP_RMIN, EP_MIN3(r, EP_R0, EP_RMAX))) + \
				   EP_PHI_INT2(   EP_GAMMA1, EP_MAX(EP_RMIN, r),         EP_MAX3(EP_RMIN, r, EP_MIN (EP_R0, EP_RMAX))) + \
				   EP_PHI_INT1(r, EP_GAMMA2, EP_MAX (EP_R0, EP_RMIN),    EP_MAX3(EP_R0, EP_RMIN, EP_MIN (r, EP_RMAX))) + \
				   EP_PHI_INT2(   EP_GAMMA2, EP_MAX3(EP_R0, EP_RMIN, r), EP_MAX3(EP_R0,                  r, EP_RMAX )) ) )

#define EP_RHO(r)                convert_mass(convert_length(convert_length(convert_length(EP_RHO_PHYS(convert_length(r, 0)), 0), 0), 0), 1)
#define EP_M(r)                  convert_mass(EP_M_PHYS(convert_length(r, 0)), 1)
#define EP_PHI(r)                convert_mass(convert_length(EP_PHI_PHYS(convert_length(r, 0)), 0), 1)

/******************************************************************************/

#endif // EXT_POT

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


#ifdef EXT_POT                     // in case of an external potential, add it to the enclosed mass for orbit calculation
#define M_ENCL(r)        (EP_M(r) + M_ENCL_(r))
#else
#define M_ENCL(r)        M_ENCL_(r)

#endif

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
