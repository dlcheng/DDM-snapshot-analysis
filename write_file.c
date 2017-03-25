#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

void out_put()
{
  FILE *fp_halo;
  FILE *fp_profile;

  char file_name_1[500];
  char file_name_2[500];

  unsigned j;           /* run over halo number */
  int n;                /* run over subhalos */
  unsigned sub_id;      /* temperary subhalo ID */
  int m;                /* run over profile bins */

  sprintf(file_name_1, "%s%s", OUTPUT_PATH, "halo_catalog.txt");
  sprintf(file_name_2, "%s%s", OUTPUT_PATH, "halo_profile.txt");

  fp_halo = fopen(file_name_1, "w+");
  fp_profile = fopen(file_name_2, "w+");

#ifndef DECAY_DARK_MATTER
  fprintf(fp_halo,    "%s\n", "#(1)HaloID  (2)HostFlag (3)Trusted_Nsub (4)Npart (5)Mvir (6)Rvir (7)cNFW_ahf (8)cNFW_Vmax (9)cNFW_Rmax (10)X (11)Y (12)Z (13)Nbin");
  fprintf(fp_profile, "%s\n", "#(1)r_right (2)shell_density (3)sphere_density (4)Vcir");
  for(j=0; j<total_halo_num; j++)
    {
    if(p_halo[j].trust_flag == 1 && p_halo[j].ahf_halo.host_flag == -1)     /* for trusted host halos */
      {
      fprintf(fp_halo, "%u %d %d %u %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %d\n", 
                        p_halo[j].ahf_halo.halo_id, p_halo[j].ahf_halo.host_flag, 
                        p_halo[j].trust_sub_num, p_halo[j].ahf_halo.npart,
                        p_halo[j].ahf_halo.Mvir, p_halo[j].ahf_halo.Rvir, 
                        p_halo[j].ahf_halo.cNFW, p_halo[j].cNFW_Vmax, p_halo[j].cNFW_Rmax,
                        p_halo[j].ahf_halo.pos[0], p_halo[j].ahf_halo.pos[1], p_halo[j].ahf_halo.pos[2],
                        PROFILE_BIN_NUM);
      for(m=0; m<PROFILE_BIN_NUM; m++)
        {
        fprintf(fp_profile, "%.6e %.6e %.6e %.6e\n", 
               p_halo[j].p_halo_profile[m].r_right, p_halo[j].p_halo_profile[m].shell_dens, 
               p_halo[j].p_halo_profile[m].sphere_dens, p_halo[j].p_halo_profile[m].cir_vel);
        }
      for(n=0; n<p_halo[j].ahf_halo.sub_num; n++)
        {
         sub_id = p_halo[j].p_ahf_subhalo[n].sub_id;
         if(p_halo[sub_id].trust_flag)
           {
           fprintf(fp_halo, "%u %d %d %u %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %d\n", 
                            p_halo[sub_id].ahf_halo.halo_id, p_halo[sub_id].ahf_halo.host_flag, 
                            p_halo[sub_id].trust_sub_num, p_halo[sub_id].ahf_halo.npart,
                            p_halo[sub_id].ahf_halo.Mvir, p_halo[sub_id].ahf_halo.Rvir, 
                            p_halo[sub_id].ahf_halo.cNFW, p_halo[sub_id].cNFW_Vmax, p_halo[sub_id].cNFW_Rmax,
                            p_halo[sub_id].ahf_halo.pos[0], p_halo[sub_id].ahf_halo.pos[1], p_halo[sub_id].ahf_halo.pos[2],
                            PROFILE_BIN_NUM);
           for(m=0; m<PROFILE_BIN_NUM; m++)
              {
               fprintf(fp_profile, "%.6e %.6e %.6e %.6e\n", 
                                   p_halo[sub_id].p_halo_profile[m].r_right, p_halo[sub_id].p_halo_profile[m].shell_dens, 
                                   p_halo[sub_id].p_halo_profile[m].sphere_dens, p_halo[sub_id].p_halo_profile[m].cir_vel);
              }
           }
        }
      }
    }
#else
  fprintf(fp_halo,    "%s\n", "#(1)HaloID  (2)HostFlag (3)Trusted_Nsub (4)Npart (5)Mvir (6)Rvir  (7)cNFW_ahf (8)cNFW_Vmax (9)cNFW_Rmax (10)X (11)Y (12)Z (13)Nbin (14)fs (15)fvir (16)fglobal (17)Nddm");
  fprintf(fp_profile, "%s\n", "#(1)r_right (2)shell_density (3)sphere_density (4)Vcir (5)dm_shell_density (6)dm_sphere_density (7)dm_shell_frac (8)dm_sphere_frac (9)fglobal");
  for(j=0; j<total_halo_num; j++)
    {
    if(p_halo[j].trust_flag == 1 && p_halo[j].ahf_halo.host_flag == -1)     /* for trusted host halos */
      {
      fprintf(fp_halo, "%u %d %d %u %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %d %.6e %.6e %.6e %u\n", 
                        p_halo[j].ahf_halo.halo_id, p_halo[j].ahf_halo.host_flag, 
                        p_halo[j].trust_sub_num, p_halo[j].ahf_halo.npart,
                        p_halo[j].ahf_halo.Mvir, p_halo[j].ahf_halo.Rvir, 
                        p_halo[j].ahf_halo.cNFW, p_halo[j].cNFW_Vmax, p_halo[j].cNFW_Rmax,
                        p_halo[j].ahf_halo.pos[0], p_halo[j].ahf_halo.pos[1], p_halo[j].ahf_halo.pos[2],
                        PROFILE_BIN_NUM, 
                        p_halo[j].dl_decay.over_all_frac[0], p_halo[j].dl_decay.over_all_frac[1],
                        global_frac, p_halo[j].dl_decay.ddm_num);
      for(m=0; m<PROFILE_BIN_NUM; m++)
        {
        fprintf(fp_profile, "%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n", 
               p_halo[j].p_halo_profile[m].r_right, p_halo[j].p_halo_profile[m].shell_dens, 
               p_halo[j].p_halo_profile[m].sphere_dens, p_halo[j].p_halo_profile[m].cir_vel,
               p_halo[j].p_halo_profile[m].shell_dens * p_halo[j].p_halo_profile[m].shell_frac,
               p_halo[j].p_halo_profile[m].sphere_dens * p_halo[j].p_halo_profile[m].sphere_frac,
               p_halo[j].p_halo_profile[m].shell_frac, p_halo[j].p_halo_profile[m].sphere_frac,               
               global_frac);
        }
      for(n=0; n<p_halo[j].ahf_halo.sub_num; n++)
        {
         sub_id = p_halo[j].p_ahf_subhalo[n].sub_id;
         if(p_halo[sub_id].trust_flag)
           {
           fprintf(fp_halo, "%u %d %d %u %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %d %.6e %.6e %.6e %u\n", 
                            p_halo[sub_id].ahf_halo.halo_id, p_halo[sub_id].ahf_halo.host_flag, 
                            p_halo[sub_id].trust_sub_num, p_halo[sub_id].ahf_halo.npart,
                            p_halo[sub_id].ahf_halo.Mvir, p_halo[sub_id].ahf_halo.Rvir, 
                            p_halo[sub_id].ahf_halo.cNFW, p_halo[sub_id].cNFW_Vmax, p_halo[sub_id].cNFW_Rmax,
                            p_halo[sub_id].ahf_halo.pos[0], p_halo[sub_id].ahf_halo.pos[1], p_halo[sub_id].ahf_halo.pos[2],
                            PROFILE_BIN_NUM,
                            p_halo[sub_id].dl_decay.over_all_frac[0], p_halo[sub_id].dl_decay.over_all_frac[1],
                            global_frac, p_halo[sub_id].dl_decay.ddm_num);
           for(m=0; m<PROFILE_BIN_NUM; m++)
             {
              fprintf(fp_profile, "%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n", 
                                  p_halo[sub_id].p_halo_profile[m].r_right, p_halo[sub_id].p_halo_profile[m].shell_dens, 
                                  p_halo[sub_id].p_halo_profile[m].sphere_dens, p_halo[sub_id].p_halo_profile[m].cir_vel,
                                  p_halo[sub_id].p_halo_profile[m].shell_dens * p_halo[sub_id].p_halo_profile[m].shell_frac,
                                  p_halo[sub_id].p_halo_profile[m].sphere_dens * p_halo[sub_id].p_halo_profile[m].sphere_frac,
                                  p_halo[sub_id].p_halo_profile[m].shell_frac, p_halo[sub_id].p_halo_profile[m].sphere_frac,               
                                  global_frac);
             }           
           }
        }
      }
    }
#endif  

  fclose(fp_halo);
  fclose(fp_profile);
  state("Finish output.");
}            /*end out_put */