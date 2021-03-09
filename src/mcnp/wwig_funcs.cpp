#include "wwig_funcs.h"

#include "DagMC.hpp"
using moab::DagMC;

#include <limits>
#include <iostream>
#include <algorithm>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>

#ifdef CUBIT_LIBS_PRESENT
#include <fenv.h>
#endif

// globals

std::map<int, moab::DagMC*> WWIG;
moab::DagMC* CURRENT_WWIG;
std::map<int, std::pair<double, double>> ww_bounds; // group: <lower, upper>

#define DGFM_SEQ   0
#define DGFM_READ  1
#define DGFM_BCAST 2

#ifdef ENABLE_WW_RAYSTAT_DUMPS

#include <fstream>
#include <numeric>

static std::ostream* ww_raystat_dump = NULL;
#endif


/* Static values used by wwigtrack_ */

static DagMC::RayHistory historyww;
static int last_nps_ww = 0;
static double last_uvw_ww[3] = {0, 0, 0};
static std::vector< std::pair<moab::DagMC*, DagMC::RayHistory > > historyww_bank;
static std::vector< DagMC::RayHistory > pblcm_historyww_stack;
static std::vector< moab::DagMC * > pblcm_wwig_stack;
static std::vector< std::pair<int, int>> pblww_bank;

static bool visited_surface_ww = false;

static bool use_dist_limit_ww = false;
static double dist_limit_ww; // needs to be thread-local


void dagmcinitww_(char* cdir, int* clen, int* max_pbl) {
  /* Load each WWIG geometry as a separate DAGMC instance.
   *
   * cdir : path to directory containing each energy group geometry.
   *
   * Files must be named in this format: <VAR>_<GROUP>.h5m, where <VAR>
   * is some arbitrary variable and <GROUP> is the integer corresponding
   * to the energy group number. Example: ww_p_003.h5m.
   * Directory should not contain any other files besides the WWIG files.
   */

    moab::ErrorCode rval;

    //char directory = *cdir;
   // For every file in cfile (dir), load into instance pointer
    struct dirent *entry = nullptr;
    DIR *dp = nullptr;
    char* dir_term = cdir;
    dir_term[*clen] = '\0';
    dp = opendir(dir_term);

    // add check here if cdir couldn't be opened
    if (dp != nullptr) {
        while (entry = readdir(dp)){

            // get file name
            char* fname = entry->d_name;
            std::string wfile = std::string(fname);

            if (wfile != "." && wfile != "..") {

                // get energy group number
                std::size_t start = wfile.rfind('_')+1;
                std::size_t len = wfile.rfind('.') - start;
                std::string grp_str = wfile.substr(start, len);
                int egrp = std::stoi(grp_str);

                // concatenate directory and filename to load as char*
                std::string full_path_str = cdir;
                full_path_str += '/';
                full_path_str += fname;
                full_path_str += '\0';
                char full_path_char[full_path_str.length() + 1];
                strcpy(full_path_char, full_path_str.c_str());

                // TO-DO: Add a check here that file ends w/ '.h5m'
                rval = wwiginit(full_path_str, egrp);
            }
        }
    }

    closedir(dp);

    pblcm_historyww_stack.resize(*max_pbl + 1);
    pblcm_wwig_stack.resize(*max_pbl + 1);
    historyww.reset();
}

moab::ErrorCode wwiginit(std::string filename, int egrp) {

  moab::ErrorCode rval;

  // make new DagMC

  WWIG[egrp] = new moab::DagMC();
  CURRENT_WWIG = WWIG[egrp];

#ifdef ENABLE_WW_RAYSTAT_DUMPS
  // the file to which ray statistics dumps will be written
  ww_raystat_dump = new std::ofstream("wwig_ww_raystat_dump.csv");
#endif

  // read geometry
  rval = CURRENT_WWIG->load_file(filename.c_str());
  if (moab::MB_SUCCESS != rval) {
    std::cerr << "WWIG failed to read input file: " << filename << std::endl;
    exit(EXIT_FAILURE);
  }

  // initialize geometry
  rval = CURRENT_WWIG->init_OBBTree();
  if (moab::MB_SUCCESS != rval) {
    std::cerr << "WWIG failed to initialize geometry and create OBB tree" <<  std::endl;
    exit(EXIT_FAILURE);
  }

  // get group energy bounds and store into ww_bounds
  moab::EntityHandle rs = CURRENT_WWIG->moab_instance()->get_root_set();
  moab::Tag el_tag;
  moab::Tag eu_tag;
  double el;
  double eu;
  rval = CURRENT_WWIG->moab_instance()->tag_get_handle("E_LOW_BOUND", 1, moab::MB_TYPE_DOUBLE, el_tag);
  rval = CURRENT_WWIG->moab_instance()->tag_get_handle("E_UP_BOUND", 1, moab::MB_TYPE_DOUBLE, eu_tag);
  rval = CURRENT_WWIG->moab_instance()->tag_get_data(el_tag, &rs, 1, &el);
  rval = CURRENT_WWIG->moab_instance()->tag_get_data(eu_tag, &rs, 1, &eu);
  ww_bounds[egrp] = std::make_pair(el,eu);

  return moab::MB_SUCCESS;

}

void wwigwritefacets_(char* ffile, int* flen) { // facet file
  // terminate all filenames with null char
  ffile[*flen]  = '\0';

  moab::ErrorCode rval = CURRENT_WWIG->write_mesh(ffile, *flen);
  if (moab::MB_SUCCESS != rval) {
    std::cerr << "WWIG failed to write mesh file: " << ffile <<  std::endl;
    exit(EXIT_FAILURE);
  }

  return;

}


void wwigangl_(int* jsu, double* xxx, double* yyy, double* zzz, double* ang) {
  moab::EntityHandle surf = CURRENT_WWIG->entity_by_index(2, *jsu);
  double xyz[3] = {*xxx, *yyy, *zzz};
  moab::ErrorCode rval = CURRENT_WWIG->get_angle(surf, xyz, ang, &historyww);
  if (moab::MB_SUCCESS != rval) {
    std::cerr << "WWIG: failed in calling get_angle" <<  std::endl;
    exit(EXIT_FAILURE);
  }

#ifdef TRACE_WWIG_CALLS
  std::cout << "angl: " << *xxx << ", " << *yyy << ", " << *zzz << " --> "
            << ang[0] << ", " << ang[1] << ", " << ang[2] << std::endl;
  moab::CartVect uvw(last_uvw_ww);
  moab::CartVect norm(ang);
  double aa = angle(uvw, norm) * (180.0 / M_PI);
  std::cout << "    : " << aa << " deg to uvw" << (aa > 90.0 ? " (!)" : "")  << std::endl;
#endif

}

void wwigchkcel_by_angle_(double* uuu, double* vvv, double* www,
                           double* xxx, double* yyy, double* zzz,
                           int* jsu, int* i1, int* j) {


#ifdef TRACE_WWIG_CALLS
  std::cout << " " << std::endl;
  std::cout << "chkcel_by_angle: vol=" << CURRENT_WWIG->id_by_index(3, *i1) << " surf=" << CURRENT_WWIG->id_by_index(2, *jsu)
            << " xyz=" << *xxx  << " " << *yyy << " " << *zzz << std::endl;
  std::cout << "               : uvw = " << *uuu << " " << *vvv << " " << *www << std::endl;
#endif

  double xyz[3] = {*xxx, *yyy, *zzz};
  double uvw[3] = {*uuu, *vvv, *www};

  moab::EntityHandle surf = CURRENT_WWIG->entity_by_index(2, *jsu);
  moab::EntityHandle vol  = CURRENT_WWIG->entity_by_index(3, *i1);

  int result;
  moab::ErrorCode rval = CURRENT_WWIG->test_volume_boundary(vol, surf, xyz, uvw, result, &historyww);
  if (moab::MB_SUCCESS != rval) {
    std::cerr << "WWIG: failed calling test_volume_boundary" << std::endl;
    exit(EXIT_FAILURE);
  }

  switch (result) {
    case 1:
      *j = 0; // inside==  1 -> inside volume -> j=0
      break;
    case 0:
      *j = 1; // outside== 0  -> outside volume -> j=1
      break;
    default:
      std::cerr << "Impossible result in wwigchkcel_by_angle" << std::endl;
      exit(EXIT_FAILURE);
  }

#ifdef TRACE_WWIG_CALLS
  std::cout << "chkcel_by_angle: j=" << *j << std::endl;
#endif

}

void wwigchkcel_(double* uuu, double* vvv, double* www, double* xxx,
                  double* yyy, double* zzz, int* i1, int* j) {


#ifdef TRACE_WWIG_CALLS
  std::cout << " " << std::endl;
  std::cout << "chkcel: vol=" << CURRENT_WWIG->id_by_index(3, *i1) << " xyz=" << *xxx
            << " " << *yyy << " " << *zzz << std::endl;
  std::cout << "      : uvw = " << *uuu << " " << *vvv << " " << *www << std::endl;
#endif

  int inside;
  moab::EntityHandle vol = CURRENT_WWIG->entity_by_index(3, *i1);
  double xyz[3] = {*xxx, *yyy, *zzz};
  double uvw[3] = {*uuu, *vvv, *www};
  moab::ErrorCode rval = CURRENT_WWIG->point_in_volume(vol, xyz, inside, uvw);

  if (moab::MB_SUCCESS != rval) {
    std::cerr << "WWIG: failed in point_in_volume" <<  std::endl;
    exit(EXIT_FAILURE);
  }

  if (moab::MB_SUCCESS != rval)
    *j = -2;
  else
    switch (inside) {
      case 1:
        *j = 0; // inside==  1 -> inside volume -> j=0
        break;
      case 0:
        *j = 1; // outside== 0  -> outside volume -> j=1
        break;
      case -1:
        *j = 1; // onboundary== -1 -> on boundary -> j=1 (assume leaving volume)
        break;
      default:
        std::cerr << "Impossible result in wwigchkcel" << std::endl;
        exit(EXIT_FAILURE);
    }

#ifdef TRACE_WWIG_CALLS
  std::cout << "chkcel: j=" << *j << std::endl;
#endif

}


void wwigdbmin_(int* ih, double* xxx, double* yyy, double* zzz, double* huge, double* dbmin) {
  double point[3] = {*xxx, *yyy, *zzz};

  // get handle for this volume (*ih)
  moab::EntityHandle vol  = CURRENT_WWIG->entity_by_index(3, *ih);

  // get distance to closest surface
  moab::ErrorCode rval = CURRENT_WWIG->closest_to_location(vol, point, *dbmin);

  // if failed, return 'huge'
  if (moab::MB_SUCCESS != rval) {
    *dbmin = *huge;
    std::cerr << "WWIG: error in closest_to_location, returning huge value from dbmin_" <<  std::endl;
  }

#ifdef TRACE_WWIG_CALLS
  std::cout << "dbmin " << CURRENT_WWIG->id_by_index(3, *ih) << " dist = " << *dbmin << std::endl;
#endif

}

void wwignewcel_(int* jsu, int* icl, int* iap) {

  moab::EntityHandle surf = CURRENT_WWIG->entity_by_index(2, *jsu);
  moab::EntityHandle vol  = CURRENT_WWIG->entity_by_index(3, *icl);
  moab::EntityHandle newvol = 0;

  moab::ErrorCode rval = CURRENT_WWIG->next_vol(surf, vol, newvol);
  if (moab::MB_SUCCESS != rval) {
    *iap = -1;
    std::cerr << "WWIG: error calling next_vol, newcel_ returning -1" << std::endl;
  }

  *iap = CURRENT_WWIG->index_by_handle(newvol);

  visited_surface_ww = true;

#ifdef TRACE_WWIG_CALLS
  std::cout << "newcel: prev_vol=" << CURRENT_WWIG->id_by_index(3, *icl) << " surf= "
            << CURRENT_WWIG->id_by_index(2, *jsu) << " next_vol= " << CURRENT_WWIG->id_by_index(3, *iap) << std::endl;

#endif
}

void wwig_surf_reflection_(double* uuu, double* vvv, double* www, int* verify_dir_change) {


#ifdef TRACE_WWIG_CALLS
  // compute and report the angle between old and new
  moab::CartVect oldv(last_uvw_ww);
  moab::CartVect newv(*uuu, *vvv, *www);

  std::cout << "surf_reflection: " << angle(oldv, newv)*(180.0 / M_PI) << std::endl;;
#endif

  // a surface was visited
  visited_surface_ww = true;

  bool update = true;
  if (*verify_dir_change) {
    if (last_uvw_ww[0] == *uuu && last_uvw_ww[1] == *vvv && last_uvw_ww[2] == *www)
      update = false;
  }

  if (update) {
    last_uvw_ww[0] = *uuu;
    last_uvw_ww[1] = *vvv;
    last_uvw_ww[2] = *www;
    historyww.reset_to_last_intersection();
  }

#ifdef TRACE_WWIG_CALLS
  else {
    // mark it in the log if nothing happened
    std::cout << "(noop)";
  }

  std::cout << std::endl;
#endif

}

void wwig_particle_terminate_() {
  historyww.reset();

#ifdef TRACE_WWIG_CALLS
  std::cout << "particle_terminate:" << std::endl;
#endif
}

// *ih              - volue index
// *uuu, *vvv, *www - ray direction
// *xxx, *yyy, *zzz - ray point
// *huge            - passed to ray_fire as 'huge'
// *dls             - output from ray_fire as 'dist_traveled'
// *jap             - intersected surface index, or zero if none
// *jsu             - previous surface index
void wwigtrack_(int* ih, double* uuu, double* vvv, double* www, double* xxx,
                 double* yyy, double* zzz, double* huge, double* dls, int* jap, int* jsu,
                 int* nps) {
  if (CURRENT_WWIG == NULL) {
    *dls = -(*huge);
    *jap = 0;
    return;
  }
  // Get data from IDs
  moab::EntityHandle vol = CURRENT_WWIG->entity_by_index(3, *ih);
  moab::EntityHandle prev = CURRENT_WWIG->entity_by_index(2, *jsu);
  moab::EntityHandle next_surf = 0;
  double next_surf_dist;

#ifdef ENABLE_WW_RAYSTAT_DUMPS
  moab::OrientedBoxTreeTool::TrvStats trv;
#endif

  double point[3] = {*xxx, *yyy, *zzz};
  double dir[3]   = {*uuu, *vvv, *www};

  /* detect streaming or reflecting situations */
  if (last_nps_ww != *nps || prev == 0) {
    // not streaming or reflecting: reset historyww
    historyww.reset();
#ifdef TRACE_WWIG_CALLS
    std::cout << "track: new historyww" << std::endl;
#endif

  } else if (last_uvw_ww[0] == *uuu && last_uvw_ww[1] == *vvv && last_uvw_ww[2] == *www) {
    // streaming -- use historyww without change
    // unless a surface was not visited
    if (!visited_surface_ww) {
    // When is this relevant?
    //  * when the previous flight ended at a transport geometry boundary,
    //    the particle will continue in the same direction from a different point
      historyww.rollback_last_intersection();
#ifdef TRACE_WWIG_CALLS
      std::cout << "     : (rbl)" << std::endl;
#endif
    }
#ifdef TRACE_WWIG_CALLS
    std::cout << "track: streaming " << historyww.size() << std::endl;
#endif
  } else {
    // not streaming or reflecting
    historyww.reset();

#ifdef TRACE_WWIG_CALLS
    std::cout << "track: reset" << std::endl;
#endif

  }

  moab::ErrorCode result = CURRENT_WWIG->ray_fire(vol, point, dir,
                                         next_surf, next_surf_dist, &historyww,
                                         (use_dist_limit_ww ? dist_limit_ww : 0)
#ifdef ENABLE_WW_RAYSTAT_DUMPS
                                         , ww_raystat_dump ? &trv : NULL
#endif
                                        );


  if (moab::MB_SUCCESS != result) {
    std::cerr << "WWIG: failed in ray_fire" << std::endl;
    exit(EXIT_FAILURE);
  }


  for (int i = 0; i < 3; ++i) {
    last_uvw_ww[i] = dir[i];
  }
  last_nps_ww = *nps;

  // Return results: if next_surf exists, then next_surf_dist will be nearer than dist_limit_ww (if any)
  if (next_surf != 0) {
    *jap = CURRENT_WWIG->index_by_handle(next_surf);
    *dls = next_surf_dist;
  } else {
    // no next surface
    *jap = 0;
    if (use_dist_limit_ww) {
      // Dist limit on: return a number bigger than dist_limit_ww
      *dls = dist_limit_ww * 2.0;
    } else {
      // Dist limit off: return huge value, triggering lost particle
      *dls = *huge;
    }
  }

  visited_surface_ww = false;

#ifdef ENABLE_WW_RAYSTAT_DUMPS
  if (ww_raystat_dump) {

    *ww_raystat_dump << *ih << ",";
    *ww_raystat_dump << trv.ray_tri_tests() << ",";
    *ww_raystat_dump << std::accumulate(trv.nodes_visited().begin(), trv.nodes_visited().end(), 0) << ",";
    *ww_raystat_dump << std::accumulate(trv.leaves_visited().begin(), trv.leaves_visited().end(), 0) << std::endl;

  }
#endif

#ifdef TRACE_WWIG_CALLS

  std::cout << "track: vol=" << CURRENT_WWIG->id_by_index(3, *ih) << " prev_surf=" << CURRENT_WWIG->id_by_index(2, *jsu)
            << " next_surf=" << CURRENT_WWIG->id_by_index(2, *jap) << " nps=" << *nps << std::endl;
  std::cout << "     : xyz=" << *xxx << " " << *yyy << " " << *zzz << " dist = " << *dls << std::flush;
  if (use_dist_limit_ww && *jap == 0)
    std::cout << " > distlimit" << std::flush;
  std::cout << std::endl;
  std::cout << "     : uvw=" << *uuu << " " << *vvv << " " << *www << std::endl;
#endif

}

void wwig_bank_push_(int* nbnk, int* icl, int* jsu) {
  if (((unsigned)*nbnk) != historyww_bank.size()) {
    std::cerr << "bank push size mismatch: F" << *nbnk << " C" << historyww_bank.size() << std::endl;
  }
  historyww_bank.push_back(std::make_pair(CURRENT_WWIG,historyww));
  pblww_bank.push_back(std::make_pair(*icl, *jsu));

#ifdef TRACE_WWIG_CALLS
  std::cout << "bank_push (" << *nbnk + 1 << ")" << std::endl;
#endif
}

void wwig_bank_usetop_(int* icl, int* jsu) {

#ifdef TRACE_WWIG_CALLS
  std::cout << "bank_usetop" << std::endl;
#endif
std::pair<moab::DagMC*, DagMC::RayHistory > banked_history;
std::pair<int, int> banked_pblww;

  if (historyww_bank.size()) {
    banked_history = historyww_bank.back();
    CURRENT_WWIG = banked_history.first;
    historyww = banked_history.second;
  } else {
    std::cerr << "wwig_bank_usetop_() called without bank historyww!" << std::endl;
  }
  if (pblww_bank.size())
  {
    banked_pblww = pblww_bank.back();
    *icl = banked_pblww.first;
    *jsu = banked_pblww.second;
  }
  else
  {
    std::cerr << "wwig_bank_usetop_() called without bank historyww!" << std::endl;
  }
}

void wwig_bank_pop_(int* nbnk) {

  if (((unsigned)*nbnk) != historyww_bank.size()) {
    std::cerr << "bank pop size mismatch: F" << *nbnk << " C" << historyww_bank.size() << std::endl;
  }

  if (historyww_bank.size()) {
    historyww_bank.pop_back();
  }
  if (pblww_bank.size()) {
    pblww_bank.pop_back();
  }


#ifdef TRACE_WWIG_CALLS
  std::cout << "bank_pop (" << *nbnk - 1 << ")" << std::endl;
#endif

}

void wwig_bank_clear_() {
  historyww_bank.clear();
#ifdef TRACE_WWIG_CALLS
  std::cout << "bank_clear" << std::endl;
#endif
}

void wwig_savpar_(int* n) {
#ifdef TRACE_WWIG_CALLS
  std::cout << "savpar: " << *n << " (" << historyww.size() << ")" << std::endl;
#endif
  pblcm_historyww_stack[*n] = historyww;
  pblcm_wwig_stack[*n] = CURRENT_WWIG;
}

void wwig_getpar_(int* n) {
#ifdef TRACE_WWIG_CALLS
  std::cout << "getpar: " << *n << " (" << pblcm_historyww_stack[*n].size() << ")" << std::endl;
  std::cout << "getpar: " << *n << " (" << pblcm_wwig_stack[*n].size() << ")" << std::endl;
#endif

  CURRENT_WWIG = pblcm_wwig_stack[*n];
  historyww = pblcm_historyww_stack[*n];
}


void wwigvolume_(int* mxa, double* vols, int* mxj, double* aras) {
  moab::ErrorCode rval;

  // get size of each volume
  int num_vols = CURRENT_WWIG->num_entities(3);
  for (int i = 0; i < num_vols; ++i) {
    rval = CURRENT_WWIG->measure_volume(CURRENT_WWIG->entity_by_index(3, i + 1), vols[i * 2]);
    if (moab::MB_SUCCESS != rval) {
      std::cerr << "WWIG: could not measure volume " << i + 1 << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  // get size of each surface
  int num_surfs = CURRENT_WWIG->num_entities(2);
  for (int i = 0; i < num_surfs; ++i) {
    rval = CURRENT_WWIG->measure_area(CURRENT_WWIG->entity_by_index(2, i + 1), aras[i * 2]);
    if (moab::MB_SUCCESS != rval) {
      std::cerr << "WWIG: could not measure surface " << i + 1 << std::endl;
      exit(EXIT_FAILURE);
    }
  }

}

void wwig_setdis_(double* d) {
  dist_limit_ww = *d;
#ifdef TRACE_WWIG_CALLS
  std::cout << "setdis: " << *d << std::endl;
#endif
}

void wwig_set_settings_(int* fort_use_dist_limit_ww) {

  if (*fort_use_dist_limit_ww) {
    std::cout << "WWIG distance limit optimization is ENABLED" << std::endl;
    use_dist_limit_ww = true;
  }

}

void wwig_init_settings_(int* fort_use_dist_limit_ww) {

  *fort_use_dist_limit_ww = use_dist_limit_ww ? 1 : 0;

}

// delete the stored data
void wwig_teardown_() {

  for (std::map<int, moab::DagMC*>::iterator ww_iter = WWIG.begin(); ww_iter != WWIG.end(); ww_iter++){
    delete ww_iter->second;
  }

}



void wwig_set_energy_group_(double* erg, int* icl) {

  std::map<int, std::pair<double, double>>::iterator it = ww_bounds.begin();
  bool found = false;
  int group;
  while (it != ww_bounds.end() && found == false)
  {
    double el = it->second.first;
    double eu = it->second.second;
    int grp = it->first;

    // if within the bounds of the group, then update group number
    if (*erg > el && *erg <= eu)
    {
      group = grp;
      found = true;
    }
    it++;
  }
  if (it == ww_bounds.end() && found == false)
  {
    // std::cerr << "WWIG: WW group look up failed E= " << *erg << std::endl;
    group = -1;
  }

  /* check if the group changed*/
  moab::DagMC* next_wwig = NULL;
  if (group > -1) {
    next_wwig = WWIG[group];
  }

  if (next_wwig != CURRENT_WWIG) {
    /* changing group */
    // reset/update rayhistory
    historyww.reset();
    // find current cell - can't do this here
    *icl = -1;
  } else {

    /* same group */

  }

  CURRENT_WWIG = next_wwig;

}

void wwig_find_cell_(double *x, double *y, double *z,
                double *u, double *v, double *w, int* icl){

  if (CURRENT_WWIG == NULL) {
    return;
  }

  double xyz[3] = { *x, *y, *z};
  double uvw[3] = { *u, *v, *w};
  moab::EntityHandle vol;
  bool found_vol = false;

  int num_cells = CURRENT_WWIG->num_entities(3);

  for (int i = 1; i <= num_cells; ++i)
  {

    vol = CURRENT_WWIG->entity_by_index(3, i);

    // check point_in_volume for each until location is found
    int inside = 0;
    moab::ErrorCode result = CURRENT_WWIG->point_in_volume(vol, xyz,
                                                           inside, uvw);

    if (moab::MB_SUCCESS != result)
    {
      std::cerr << "WWIG failed in point_in_volume" << std::endl;
      exit(EXIT_FAILURE);
    }

    if (inside == 1)
    {
      // inside volume, can stop searching
      // update current cell
      *icl = CURRENT_WWIG->index_by_handle(vol);
      found_vol = true;
      break;
    }
  }

  if (found_vol == false)
  {
    // no volume found
    std::cerr << "WWIG failed to find current cell" << std::endl;
    exit(EXIT_FAILURE);
  }

  // MOAB call does not work, use above iteration
  // CURRENT_WWIG->find_volume(xyz, volume, uvw);
  // *icl = CURRENT_WWIG->index_by_handle(volume);
}

void wwig_lookup_(int *jap, double *wwval)
{
  // look up the WW value on the geometry surface

  moab::EntityHandle surf = CURRENT_WWIG->entity_by_index(2, *jap);
  std::vector<moab::Tag> tag_handles;
  moab::ErrorCode result = CURRENT_WWIG->moab_instance()->tag_get_tags_on_entity(surf, tag_handles);

  if (moab::MB_SUCCESS != result)
  {
    std::cerr << "DAGMC: WW tag list lookup failed" << std::endl;
    exit(EXIT_FAILURE);
  }

  /* TODO: move identification of relevant tag to file opening */
  std::string wwn = "ww_n";
  std::string wwp = "ww_p";
  std::string wwv = "ww_val";
  double data;
  for (int i = 0; i < tag_handles.size(); i++)
  {
    std::string name;
    CURRENT_WWIG->moab_instance()->tag_get_name(tag_handles[i], name);
    if (name == wwn || name == wwp || name == wwv)
    {
      moab::ErrorCode rval = CURRENT_WWIG->moab_instance()->tag_get_data(tag_handles[i], &surf, 1, (void *)&data);
      if (moab::MB_SUCCESS != rval) {
        std::cerr << "WWIG failed to look up WW surface value for tag " << name << std::endl;
      }
      *wwval = data;
      return;
    }
  }
}

void ww_surf_check_(double *erg, double *x, double *y, double *z, double *wgt, int *jap)
{
  /* check for each ww look up. (Call from MCNP)
   Print:
    - position: x, y, z
    - energy: group id (get from mcnp), energy bounds (from ww_bounds), CURRENT_WWIG energy bounds
    - ww info: CURRENT_WWIG ID, ww value, particle weight
  */

  std::cout << "" << std::endl;
  std::cout << "****WWIG SURFACE CHECK: ****" << std::endl;

  // Position
  std::cout << "Position x y z (mcnp): " << *x << " " << *y << " " << *z << std::endl;

  // Energy bounds:
  std::cout << "Energy (mcnp):      " << *erg << std::endl;
  // look up bounds on current wwig, make sure they match the ones we are using
  moab::EntityHandle rs = CURRENT_WWIG->moab_instance()->get_root_set();
  moab::Tag el_tag;
  moab::Tag eu_tag;
  double el;
  double eu;
  CURRENT_WWIG->moab_instance()->tag_get_handle("E_LOW_BOUND", 1, moab::MB_TYPE_DOUBLE, el_tag);
  CURRENT_WWIG->moab_instance()->tag_get_handle("E_UP_BOUND", 1, moab::MB_TYPE_DOUBLE, eu_tag);
  CURRENT_WWIG->moab_instance()->tag_get_data(el_tag, &rs, 1, &el);
  CURRENT_WWIG->moab_instance()->tag_get_data(eu_tag, &rs, 1, &eu);
  std::cout << "E_low (tags):      " << el << std::endl;
  std::cout << "E_upper (tags):    " << eu << std::endl;

  // weight values
  double wwval;
  wwig_lookup_(jap, &wwval);
  std::cout << "Pre-WW Weight (mcnp):  " << *wgt << std::endl;
  std::cout << "WW Low bound (wwig):   " << wwval << std::endl;
  std::cout << "WW Up bound (calc x5):   " << wwval*5.0 << std::endl;
}

void ww_wgt_check_(double *wgt){
  /* Check that weight is as expected after each WW lookup event */
  std::cout << "Post-WW Weight (mcnp): " << *wgt << std::endl;
}