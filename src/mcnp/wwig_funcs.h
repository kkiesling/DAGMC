#ifndef WWIG_MCNP_IFACE_H
#define WWIG_MCNP_IFACE_H

#include <map>
#include <string>

#include "moab/Interface.hpp"
#include "moab/CartVect.hpp"

#include <fstream>
#include <sstream>

#include "pyne.h"
#include <unistd.h>

#ifdef __cplusplus
extern "C" {
#endif

/* initialize WWIG from FORTRAN main
 * @param max_pbl - The maximum index of the pblcm (temporary particle state) array
 *                  This is the largest n that will arrive in calls to savpar and getpar
 */

void dagmcinitww_(char* cdir, int* clen, int* max_pbl);

moab::ErrorCode wwiginit(std::string filename, int egrp);

/* Add the current particle state to the bank */
void wwig_bank_push_(int *nbnk, int *icl, int *jsu);

/* Revert to the most recently banked particle state */
void wwig_bank_usetop_(int *icl, int *jsu);

/* Remove and forget the most recently banked particle state.
 * Called after wwig_bank_usetop_() */
void wwig_bank_pop_(int* nbnk);

/* Remove all entries in the particle bank */
void wwig_bank_clear_() ;

/* Save the current particle state to temporary index n */
void wwig_savpar_(int* n);

/* Reset the current particle state using temporary index n*/
void wwig_getpar_(int* n);

/* write facet file after initialization and OBBTree generation */
void wwigwritefacets_(char* ffile, int* flen);

/* Get normal of surface with id *jsu at location (*xxx,*yyy,*zzz) and store
   in three doubles at ang (an arry of length 3) */
void wwigangl_(int* jsu, double* xxx, double* yyy, double* zzz, double* ang);


/* Given point and direction, determine if the particle is going into
 * cell i1.  Assume the point is on the surface jsu.  Return j=0 if
 * the particle is directed into the cell, j=1 if it is leaving the cell.
 * This is used as a point-in-volume query for points that are known
 * to be on a surface.
 */
void wwigchkcel_by_angle_(double* uuu, double* vvv, double* www,
                           double* xxx, double* yyy, double* zzz,
                           int* jsu, int* i1, int* j);

/* Point-in-volume query.  Determine if the particle at given coordinates
 * is inside or outside of cell i1.  Return j=1 if outside or on boundary,
 * and j=0 if inside.
 */
void wwigchkcel_(double* uuu, double* vvv, double* www, double* xxx,
                  double* yyy, double* zzz, int* i1, int* j);

/* Determine distance to nearest surface
 * *ih - current RefVolume ID
 * *xxx, *yyy, *zzz - current point
 * *huge - passed definition of a large number
 * *dbmin - Output, distance to nearest surface
 */
void wwigdbmin_(int* ih,
                 double* xxx, double* yyy, double* zzz,
                 double* huge, double* dbmin);


/* Get next Cell
 * *jsu - Surface to cross
 * *icl - Previous cell ID
 * *iap - Next cell ID (output parameter)
 */
void wwignewcel_(int* jsu, int* icl, int* iap);

/**
 * Tell wwig that a particle has changed direction at the most recently reached surface boundary,
 * and reset ray history as necessary.
 * The parameters are the new particle direction and a flag.
 * If *verify_dir_change==0, assume without checking that a direction change has actually happened.
 * If *verify_dir_change==1, compare the new uvw to the last known one, and hold on to
 * all ray history if not.  (This mode is used during electron transport,
 * which occasionally calls this function without having changed the direction;
 * it's less invasive to do the check in C++ instead of Fortran.)
 */
void wwig_surf_reflection_(double* uuu, double* vvv, double* www, int* verify_dir_change);

void wwig_particle_terminate_();

/* Do ray fire
 * *ih  - Volume ID to do ray fire against
 * *jsu - ? (RefFace ID)
 * *nps - ?
 * (*uuu,*vvv,*www) - Ray direction vector
 * (*xxx,*yyy,*zzz) - Ray point
 * *huge - ?
 * *dls - output distnace of intersection
 * *jap - Next intersected surface, or zero if none
 */
void wwigtrack_(int* ih, double* uuu, double* vvv, double* www, double* xxx,
                 double* yyy, double* zzz, double* huge, double* dls, int* jap, int* jsu,
                 int* nps);

/* Measure entities
 * vols - 2xN array where first column contains, as output, measure of every volume.
 * aras - 2xN array where first column contains, as output, measure of every surface
 */
void wwigvolume_(int* mxa, double* vols, int* mxj, double* aras);

/* Set distance limit */
void wwig_setdis_(double* d);

void wwig_set_settings_(int* use_dist_limits);

void wwig_init_settings_(int* use_dist_limit);

void wwig_teardown_();

void wwig_set_energy_group_(double* erg, int* icl);

void wwig_find_cell_(double *x, double *y, double *z,
                double *u, double *v, double *w, int* icl);

void wwig_lookup_(int *jap, double *wwval);

#ifdef __cplusplus
} // extern "C"
#endif


#endif /* WWIG_MCNP_IFACE_H */
