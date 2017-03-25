#ifndef ALLVAR_H
#define ALLVAR_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "define.h"

typedef struct gadget_head GADGET_HEAD;
typedef struct gadget_head *PPT_GADGET_HEAD;
typedef struct part_info PART_INFO;
typedef struct part_info *PPT_PART_INFO;
typedef struct file_info FILE_INFO;
typedef struct file_info *PPT_FILE_INFO;
typedef struct ahf_halo_info AHF_HALO_INFO;
typedef struct ahf_halo_info *PPT_AHF_HALO_INFO;

#ifndef NO_AHF_PROFILE
typedef struct ahf_profile_info AHF_PROFILE_INFO;
typedef struct ahf_profile_info *PPT_AHF_PROFILE_INFO;
#endif

typedef struct halo_profile_info HALO_PROFILE_INFO;
typedef struct halo_profile_info *PPT_HALO_PROFILE_INFO;

typedef struct ahf_subhalo_info AHF_SUBHALO_INFO;
typedef struct ahf_subhalo_info *PPT_AHF_SUBHALO_INFO;
typedef struct decay_info DECAY_INFO;
typedef struct decay_info *PPT_DECAY_INFO;
typedef struct halo_info HALO_INFO;
typedef struct halo_info *PPT_HALO_INFO;

struct gadget_head
{  
  unsigned int Npart[6];
  double Massarr[6];
  double Time;
  double Redshift;
  int FlagSfr;
  int FlagFeedback;
  int Nall[6];
  int  FlagCooling;
  int NumFiles;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  int FlagAge;
  int FlagMetals;
  int NallHW[6];
  int Flag_entr_ics;
#ifndef DECAY_DARK_MATTER
  char unused[60];         
#else
  unsigned int Original_num;                    /* the original particle number before decay */
  double Original_soften_length;                /* original particle's softening length  */
  char fill[48];                                /* new unused */
#endif
}; 

struct part_info
{
  unsigned id;                                  /* ID in int */
  float pos[3];
  float vel[3];
  float mass;
  int envir_flag;                               /* whether this particle belongs to the environment, 1 for yes, 0 for no */  
};

struct file_info
{
  unsigned file_tag;                            /* tag of the file, start from 0 */
  unsigned halo_num_in_file;        			      /* total number of halos in the AHF file */
};	

struct ahf_halo_info
{
  unsigned halo_id;
  int host_flag;                                /* -1 if it is host, otherwise, ID of host halo */
  unsigned sub_num;                             /* number of subhalos inside */
  double Mvir;                                  /* mass of the halo, in unit of M_sun/h */
  unsigned npart;                               /* number of particles of the halo */
  double pos[3];                                /* position of the center, in unit of kpc/h */
  double vel[3];                                /* peculiar velocity of halo, in unit of km/s */
  double Rvir;                                  /* virial radius in comoving scale */
  double Rvmax;                                 /* position of rotation curve maximum */
  double r2;                                    /* position where %rho r^2 peaks */
  double mbp_offset;                            /* offset between most bound particle and halo center */
  double com_offset;                            /* offset between cnter-of-mass and halo center */
  double Vmax;                                  /* maximum of rotation curve */
  double v_esc;                                 /* escape velocity at Rvir */
  double sigV;                                  /* 3D velocity dispersion */
  double lambda;                                /* spin parameter (Bullock et al. 2001) */
  double lambdaE;                               /* classical spin parameter (Peebles) */
  double vec_L[3];                              /* orientation of angular momentum */
  double b_c;                                   /* b over c, shape parameter, a<b<c */
  double a_c;                                   /* a over c, shape parameter*/
  double vec_c[3];                              /* direction of a, shape parameter */
  double vec_b[3];                              /* direction of b, shape parameter */
  double vec_a[3];                              /* direction of c, shape parameter */
  double ovdens;                                /* overdensity at virial radius */
  int nbins;                                    /* number of bins used for ahf_profile_info, and decay_info */
  double fMhires;                               /* mass fraction in high resolution particles for zoom simulations */
  double Ekin;                                  /* kinetic energy, in uit of M_sun/h*(km/s)^2 */
  double Epot;                                  /* potential energy, same unit as Ekin */
  double SurfP;                                 /* surface pressure(Shaw et al. 2006), same unit as Ekin */
  double Phi0;                                  /* %Phi_0 used in unbinding procedure, in unit of (km/s)^2 */ 
  double cNFW;                                  /* NFW concentration (Prada et al. 2012) */     
};

#ifndef NO_AHF_PROFILE
struct ahf_profile_info
{
  double r_right;                               /* right edge of the radial bin */
  unsigned npart;                               /* number of particles inside sphere of radius r */
  double M_in_r;                                /* mass inside sphere of radius r */
  double ovdens;                                /* over density, something not clear */
  double dens;                                  /* dnesity, don't understand what it is */
  double vcirc;                                 /* rotation curve */
  double v_esc;                                 /* escape velocity from material inside r-sphere */
  double sigV;                                  /* velocity dispersion of material inside r-sphere */
  double L[3];                                  /* angular momentum of material inside r-sphere */
  double b_c;                                   /* b/c, shape parameter of the sphere, a<b<c */
  double a_c;                                   /* a/c, shape parameter of the sphere, a<b<c */
  double vec_c[3];                              /* direction of c */
  double vec_b[3];                              /* direction of b */
  double vec_a[3];                              /* direction of a */
  double Ekin;                                  /* kinetic energy inside r-sphere, in unit of M_sun/h * (km/s)^2 */
  double Epot;                                  /* potential energy inside r-sphere, the same unit as Ekin */		
};
#endif

struct halo_profile_info
{
  double r_right;                               /* right boundary of the bin, physical distance, no log */	
  double r_density;                             /* the log mid of the bin, for density plot physical distance, no log */
  double shell_dens;                            /* total density of the shell */
  double sphere_dens;                           /* the mean density inside the sphere */
  double shell_mass;                            /* total mass in the shell */
  double sphere_mass;                           /* total mass in the sphere */
  double cir_vel;                               /* circular velocity at r_right, peculiar velocity in unit km/s */
#ifdef DECAY_DARK_MATTER
  double shell_frac;                            /* dm fraction within the shell */
  double sphere_frac;                           /* dm fraction within the sphere of r_right */
#endif  
};

struct ahf_subhalo_info
{
  unsigned  sub_id;                             /* subhalo position in the halo array */			
};

struct decay_info
{
  int ddm_num;                                   /* the particle number of undecayed dm, mother particle */
  double over_all_frac[2];                       /* mass fraction of DM inside r_s[0] and r_vir[1] */
};	

struct halo_info
{ 
  int trust_flag;                                /* controlled only by particle number, good 1, bad 0 */     
  unsigned int trust_sub_num;                    /* the number of trusted subhaloes, may be smaller than ahf_halo.sub_num */ 
  double cNFW_Rmax;                              /* concentration from the radius, not used anymore */
  double cNFW_Vmax;                              /* concentration from the max circular velocity */
                                            
  AHF_HALO_INFO ahf_halo;
#ifdef DECAY_DARK_MATTER  
  DECAY_INFO dl_decay;
#endif

/* pointers */  
  unsigned int * p_ahf_part;                     /* point to particle ID */
#ifndef NO_AHF_PROFILE  
  PPT_AHF_PROFILE_INFO  p_ahf_profile;           /* point to ahf_profile array */
#endif  
  PPT_HALO_PROFILE_INFO p_halo_profile;
  PPT_AHF_SUBHALO_INFO  p_ahf_subhalo;           /* point to the subhalo array */    	
};

/* some global parameters */
extern int ahf_file_num;                         /* the number of AHF files produced by multi CPU */
extern int gadget_file_num;                      /* the file number of gadget files */
extern PPT_FILE_INFO p_file_info;                /* point to the file info array */
extern unsigned total_halo_num;                  /* total halo number in all files */
extern long int total_part_num;                  /* total particle numbe in the gadget file */
extern PPT_HALO_INFO p_halo;                     /* point to the halo array */

extern GADGET_HEAD g_head;                       /* the header gadget file */
extern double boxsize;                           /* boxsize in unit of kpc/h */
extern PPT_PART_INFO p_part;                     /* point to the array of all particles */

/* related to the files */
extern FILE **fp_ahf_halo;                       /* point to multi AHF.halos files, if any */
#ifndef NO_AHF_PROFILE
extern FILE **fp_ahf_profile;                    /* point to multi AHF.profiles files, if any */
#endif
extern FILE **fp_ahf_part;                       /* point to multi AHF.particles files, if any */
extern FILE **fp_ahf_subhalo;                    /* point to multi AHF.substructures files, if any */
extern FILE **fp_gadget;                         /* point to multi gadget files, if any */

#ifdef DECAY_DARK_MATTER
extern double total_mass;                        /* total mass in the file */
extern int    decay_freq;                        /* the decay frequency */ 
extern double global_decayed_mass;               /* decayed mass in the file */
extern double global_frac;                       /* mass fraction of decayed dark matter */
#endif

extern double total_halo_mass;                   /* total mass in the halos, no matter it's good or not by our criteria */
extern double total_envir_mass;                  /* total mass in the environment */

#ifdef DECAY_DARK_MATTER
extern double decayed_halo_mass;                 /* DM mass in all halos */
extern double decayed_envir_mass;                /* DM mass in the environment */
#endif

#endif
