#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "allvars.h"

/* some global parameters */
int ahf_file_num;                            /* the number of AHF files produced by multi CPU */
int gadget_file_num;                         /* the file number of gadget files */
PPT_FILE_INFO p_file_info;                   /* point to the file info array */
unsigned total_halo_num;                     /* total halo number in all files */
long int total_part_num;                     /* total particle numbe in the gadget file */
PPT_HALO_INFO p_halo;                        /* point to the halo array */

GADGET_HEAD g_head;                          /* the header gadget file */
double boxsize;                              /* boxsize in unit of kpc/h */
PPT_PART_INFO p_part;                        /* point to the array of all particles */

/* related to the files */
FILE **fp_ahf_halo;                          /* point to multi AHF.halos files */
#ifndef NO_AHF_PROFILE
FILE **fp_ahf_profile;                       /* point to multi AHF.profiles files */
#endif
FILE **fp_ahf_part;                          /* point to multi AHF.particles files */
FILE **fp_ahf_subhalo;                       /* point to multi AHF.substructures files */
FILE **fp_gadget;                            /* point to multi gadget files, if any */

#ifdef DECAY_DARK_MATTER
double total_mass;                           /* total mass in the file */
int    decay_freq;                           /* the decay frequency */ 
double global_decayed_mass;                  /* decayed mass in the file */
double global_frac;                          /* mass fraction of decayed dark matter */
#endif

double total_halo_mass;                      /* total mass in the halos, no matter it's good or not by our criteria */
double total_envir_mass;                     /* total mass in the environment */

#ifdef DECAY_DARK_MATTER
double decayed_halo_mass;                    /* DM mass in all halos */
double decayed_envir_mass;                   /* DM mass in the environment */
#endif

