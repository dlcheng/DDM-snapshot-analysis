#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>


void init_halo_profile_array(unsigned j)
{
  int m;  
  /* 
   * The most inner region is assumed to be 1kpc/h,
   * Halo with 100kpc/h is about 2.3e11 Msun/h.
   * Alternatively, we can use an adaptive distance, where the smallest radius to define log scale r_right is 1% of
   * the virial radius.
   * In practice, particles are collected from radius 0.
   */
  double log_co_vr = log10(p_halo[j].ahf_halo.Rvir);                       /* log10 virial radius, which have been converted to physical radius */
  double log_co_ir = log_co_vr - 2.0;                                      /* 0.01 of the virial radius in log */
  double log_radius_dis = 2.0 / (double) PROFILE_BIN_NUM;
  
  for(m = 0; m < PROFILE_BIN_NUM; m++)
    {
     p_halo[j].p_halo_profile[m].r_right = pow(10, (m + 1) * log_radius_dis + log_co_ir);
     if(m !=0 )
       p_halo[j].p_halo_profile[m].r_density = pow(10, (m + 0.5) * log_radius_dis + log_co_ir);
     else
       p_halo[j].p_halo_profile[m].r_density = p_halo[j].p_halo_profile[m].r_right / 2.0;
     p_halo[j].p_halo_profile[m].shell_dens = 0;
     p_halo[j].p_halo_profile[m].sphere_dens = 0;
     p_halo[j].p_halo_profile[m].shell_mass = 0;
     p_halo[j].p_halo_profile[m].sphere_mass = 0;
     p_halo[j].p_halo_profile[m].cir_vel = 0;
#ifdef DECAY_DARK_MATTER
     p_halo[j].p_halo_profile[m].shell_frac = 0;
     p_halo[j].p_halo_profile[m].sphere_frac = 0;
#endif     
    }

#ifdef DECAY_DARK_MATTER
  p_halo[j].dl_decay.ddm_num = 0;
  p_halo[j].dl_decay.over_all_frac[0] = 0;              /* dm fraction within rs */
  p_halo[j].dl_decay.over_all_frac[1] = 0;              /* dm fraction within rvir */
#endif
}                                  /* end init_halo_profile_array */

void cal_halo_profile(unsigned j)
{
  unsigned int p;                 /* particle index */	
  int b;                          /* bin index */
  int m;
  double dis;                     /* distance to the halo center */
  double r_left;
  double pi = 3.14159265359;
  unsigned local_id;
  char s;
#ifdef DECAY_DARK_MATTER
  int flag_dm = 0;
  double mass_in_rs = 0;
  double mass_in_rvir = 0;
  double rs;
#endif  
	
  for(p = 0; p < p_halo[j].ahf_halo.npart; p++)
    {
    local_id = p_halo[j].p_ahf_part[p];                                     /* the local_id records the posion of a particle in the particle array */
	  dis = distance_to_center(p_part[local_id].pos, p_halo[j].ahf_halo.pos); /* physical distance of a particle to the halo center */
	
	  if(dis > p_halo[j].ahf_halo.Rvir && (fabs(dis - p_halo[j].ahf_halo.Rvir) > 0.01 * p_halo[j].ahf_halo.Rvir))
	    {
	    printf("Halo[%d], particle[%d], dis = %f, Rvir = %f\n", j, p, dis, p_halo[j].ahf_halo.Rvir);
	    printf("Press any key to continue.");
	    scanf("%c", &s);
      }

	  for(b = 0; b < PROFILE_BIN_NUM; b++)
	    {
	    if(p_halo[j].p_halo_profile[b].r_right >= dis)
	      break;
      }
      
    if(dis >  p_halo[j].p_halo_profile[PROFILE_BIN_NUM - 1].r_right)   
      b = PROFILE_BIN_NUM - 1;  
      
	  p_halo[j].p_halo_profile[b].shell_mass += p_part[local_id].mass;

#ifdef DECAY_DARK_MATTER
    mass_in_rvir += p_part[local_id].mass;   
    flag_dm = 0;                                    /* before judgement, it is consider as ddm particle */  
    if(local_id > (g_head.Original_num - 1))
      {
       flag_dm = 1;                                 /* this is a dm particle */              
       p_halo[j].p_halo_profile[b].shell_frac += p_part[local_id].mass;
       p_halo[j].dl_decay.over_all_frac[1] += p_part[local_id].mass;
      }
    else
      p_halo[j].dl_decay.ddm_num++;                  /* this is a ddm particle */

    if(p_halo[j].ahf_halo.cNFW > 0 )
       {
        rs = p_halo[j].ahf_halo.Rvir / p_halo[j].ahf_halo.cNFW;
        if(dis <= rs)
          {
          mass_in_rs += p_part[local_id].mass;
          if(flag_dm)
            p_halo[j].dl_decay.over_all_frac[0] += p_part[local_id].mass;
          }
       }      
#endif      
    }

#ifdef DECAY_DARK_MATTER
  p_halo[j].dl_decay.over_all_frac[0] /= mass_in_rs;
  p_halo[j].dl_decay.over_all_frac[1] /= mass_in_rvir;
#endif     

  for(b = 0; b < PROFILE_BIN_NUM; b++)
    {
      for(m = 0; m <= b; m++)
       {
        p_halo[j].p_halo_profile[b].sphere_mass += p_halo[j].p_halo_profile[m].shell_mass;
#ifdef DECAY_DARK_MATTER
        p_halo[j].p_halo_profile[b].sphere_frac += p_halo[j].p_halo_profile[m].shell_frac;
#endif        
       }
        
      if(b == 0)
        r_left = 0;                            /* the left boundary of the first shell is just 0 */
      else
        r_left = p_halo[j].p_halo_profile[b-1].r_right;
        
      p_halo[j].p_halo_profile[b].shell_dens = 3 * p_halo[j].p_halo_profile[b].shell_mass / 4.0 / pi /(pow(p_halo[j].p_halo_profile[b].r_right,3) - pow(r_left,3));
      p_halo[j].p_halo_profile[b].sphere_dens = 3 * p_halo[j].p_halo_profile[b].sphere_mass / 4.0 / pi / pow(p_halo[j].p_halo_profile[b].r_right,3);
      p_halo[j].p_halo_profile[b].cir_vel = sqrt(43007.1 * p_halo[j].p_halo_profile[b].sphere_mass / 1e10 / p_halo[j].p_halo_profile[b].r_right);
    }
#ifdef DECAY_DARK_MATTER   
  for(b = 0; b < PROFILE_BIN_NUM; b++)                    /* must be out of the previous for */
    {
      p_halo[j].p_halo_profile[b].shell_frac /= p_halo[j].p_halo_profile[b].shell_mass;
      p_halo[j].p_halo_profile[b].sphere_frac /= p_halo[j].p_halo_profile[b].sphere_mass;
    }
#endif
}                                  /* end cal_halo_profile */


void get_cNFW(unsigned j)
{
  double Rmax = p_halo[j].p_halo_profile[0].r_right;
  double Vmax = p_halo[j].p_halo_profile[0].cir_vel;
  double Vvir = p_halo[j].p_halo_profile[PROFILE_BIN_NUM - 1].cir_vel;	
  double Rvir = p_halo[j].ahf_halo.Rvir;
  
  int m;
  
  for(m = 1; m < PROFILE_BIN_NUM; m++)      /* starting from the second shell */
    {
	  if(p_halo[j].p_halo_profile[m].cir_vel > Vmax)
	    {
		    Vmax = p_halo[j].p_halo_profile[m].cir_vel;
		    Rmax = p_halo[j].p_halo_profile[m].r_right;  
	    }   
    }
  
  if(Rmax/Rvir < 1)
    {
     p_halo[j].cNFW_Rmax = 2.15 * Rvir / Rmax;  
    }
  else
    p_halo[j].cNFW_Rmax = -1;    

#ifdef DEBUG
  printf("Vmax/Vvir = %.6e\n", Vmax/Vvir);
#endif
  if(Vmax/Vvir > 1.0 && Vmax/Vvir < 45)   
    {
     p_halo[j].cNFW_Vmax = cNFW_from_vel(Vmax/Vvir);
    }
   else
    p_halo[j].cNFW_Vmax = -1;          
}                               /* end cal_cNFW */	


double distance_to_center(float *a, double *c)
{
  int i,j,k;
  double d = 0;
  double dmin = 0;

  for(i=0; i<3 ; i++)
    dmin += (a[i] - c[i]) * (a[i] - c[i]);
    
  for(i=-1; i<2; i++)
    for(j=-1; j<2; j++)
      for(k=-1; k<2; k++)
        {
         d =  (a[0] - c[0] + i * g_head.BoxSize) * (a[0] - c[0] + i * g_head.BoxSize);
         d += (a[1] - c[1] + j * g_head.BoxSize) * (a[1] - c[1] + j * g_head.BoxSize);
         d += (a[2] - c[2] + k * g_head.BoxSize) * (a[2] - c[2] + k * g_head.BoxSize);

         if(d < dmin)
           dmin = d;		
	      }
    
  return sqrt(dmin) / (1.0 + g_head.Redshift); /* this is the phyiscal distance, corrected with periodic boundary condition */
}                   /* end distance_to_center */
