#define  PHYSICS                        RMHD
#define  DIMENSIONS                     2
#define  COMPONENTS                     3
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     NO
#define  FORCED_TURB                    NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  DIMENSIONAL_SPLITTING          NO
#define  NTRACER                        0
#define  USER_DEF_PARAMETERS            8

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   CONSTRAINED_TRANSPORT

/* -- user-defined parameters (labels) -- */

#define  PRS_IN                         0
#define  PRS_OUT                        1
#define  RHO_IN                         2
#define  RHO_OUT                        3
#define  BMAG                           4
#define  THETA                          5
#define  PHI                            6
#define  RADIUS                         7

/* [Beg] user-defined constants (do not change this line) */

#define  LIMITER                        VANLEER_LIM
#define  CT_EMF_AVERAGE                 CT_CONTACT
#define  RECONSTRUCT_4VEL               NO

/* [End] user-defined constants (do not change this line) */
