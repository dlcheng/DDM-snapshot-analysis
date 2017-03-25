#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

/* This function fills the p_halo array */

void init_halo_array()
{
  unsigned i, j;
  unsigned filled_halo_num = 0;
  
  alloc_halo_array();
   
  for(i=0; i<ahf_file_num; i++)
     {
     for(j=0; j<p_file_info[i].halo_num_in_file; j++)
       {
       p_halo[filled_halo_num + j].trust_sub_num = 0;       /* set to 0 at first */
#ifdef DEBUG       
       printf("Construct AHF halo INFO for halo %d.\n", filled_halo_num + j);    
#endif  	   
       construct_ahf_halo_info(i, filled_halo_num + j);

#ifdef DEBUG       
       printf("Construct AHF halo PART for halo %d.\n", filled_halo_num + j);   
#endif     
       construct_ahf_halo_part(i, filled_halo_num + j);  
             
#ifndef NO_AHF_PROFILE       
#ifdef DEBUG       
       printf("Construct AHF halo PROF for halo %d.\n", filled_halo_num + j); 
#endif 
       construct_ahf_halo_profile(i, filled_halo_num + j);      
#endif

#ifdef DEBUG       
       printf("Construct OUR halo PROF for halo %d.\n", filled_halo_num + j);
#endif       
       construct_halo_profile(filled_halo_num + j);           /* need no AHF files here */ 

       init_halo_trust_flag(filled_halo_num + j);             /* this is done for both hosts and subs */
	     }
	   filled_halo_num += p_file_info[i].halo_num_in_file; 
     }
     
  construct_ahf_subhalo();                                    /* lots of limitations now */
    
  state("Finish initializing halo array.");   
}                    /* end init_halo_array */

void construct_ahf_halo_info(unsigned i, unsigned j)
{
  fscanf(fp_ahf_halo[i], 
            "%u %d \
             %u %lf \
             %u \
             %lf %lf %lf \
             %lf %lf %lf \
             %lf %lf %lf \
             %lf %lf \
             %lf %lf %lf \
             %lf %lf \
             %lf %lf %lf \
             %lf %lf \
             %lf %lf %lf \
             %lf %lf %lf \
             %lf %lf %lf \
             %lf \
             %d \
             %lf \
             %lf %lf %lf \
             %lf \
             %lf\n", 
             &p_halo[j].ahf_halo.halo_id, &p_halo[j].ahf_halo.host_flag,
             &p_halo[j].ahf_halo.sub_num, &p_halo[j].ahf_halo.Mvir,
             &p_halo[j].ahf_halo.npart,
             &p_halo[j].ahf_halo.pos[0], &p_halo[j].ahf_halo.pos[1], &p_halo[j].ahf_halo.pos[2], 
             &p_halo[j].ahf_halo.vel[0], &p_halo[j].ahf_halo.vel[1], &p_halo[j].ahf_halo.vel[2],
             &p_halo[j].ahf_halo.Rvir, &p_halo[j].ahf_halo.Rvmax, &p_halo[j].ahf_halo.r2, 
             &p_halo[j].ahf_halo.mbp_offset, &p_halo[j].ahf_halo.com_offset, 
             &p_halo[j].ahf_halo.Vmax, &p_halo[j].ahf_halo.v_esc, &p_halo[j].ahf_halo.sigV, 
             &p_halo[j].ahf_halo.lambda, &p_halo[j].ahf_halo.lambdaE, 
             &p_halo[j].ahf_halo.vec_L[0], &p_halo[j].ahf_halo.vec_L[1], &p_halo[j].ahf_halo.vec_L[2],
             &p_halo[j].ahf_halo.b_c, &p_halo[j].ahf_halo.a_c, 
             &p_halo[j].ahf_halo.vec_c[0], &p_halo[j].ahf_halo.vec_c[1], &p_halo[j].ahf_halo.vec_c[2],
             &p_halo[j].ahf_halo.vec_b[0], &p_halo[j].ahf_halo.vec_b[1], &p_halo[j].ahf_halo.vec_b[2],
             &p_halo[j].ahf_halo.vec_a[0], &p_halo[j].ahf_halo.vec_a[1], &p_halo[j].ahf_halo.vec_a[2],
             &p_halo[j].ahf_halo.ovdens,
             &p_halo[j].ahf_halo.nbins,
             &p_halo[j].ahf_halo.fMhires,
             &p_halo[j].ahf_halo.Ekin, &p_halo[j].ahf_halo.Epot, &p_halo[j].ahf_halo.SurfP,
             &p_halo[j].ahf_halo.Phi0,
             &p_halo[j].ahf_halo.cNFW);

  p_halo[j].ahf_halo.Rvir  /= (1.0 + g_head.Redshift);     /* the native AHF halo radius is in comoving corrdinates
                                                            * change to physical one.
                                                            */
  p_halo[j].ahf_halo.Rvmax /= (1.0 + g_head.Redshift);
  p_halo[j].ahf_halo.r2    /= (1.0 + g_head.Redshift);  
  p_halo[j].ahf_halo.mbp_offset /= (1.0 + g_head.Redshift);        
  p_halo[j].ahf_halo.com_offset /= (1.0 + g_head.Redshift);                                                                      		
}                    /* end construct_ahf_halo_info */	

void construct_ahf_halo_part(unsigned i, unsigned j)
{
  int halo_part_num;
  unsigned local_id, n;
  int temp_type;
  
  fscanf(fp_ahf_part[i], "%d\n", &halo_part_num);	
  
  if( halo_part_num != p_halo[j].ahf_halo.npart)
     end_run("MATCHING HALO PARTICLE NUMBER");
	
  alloc_halo_part_array(j);
  
  for(n=0; n<halo_part_num; n++)
    {
    fscanf(fp_ahf_part[i],"%u %d\n", &local_id, &temp_type);
   
    if(local_id != p_part[local_id].id)
      end_run("MATCHING THE PARTICLE ID 1");
    else	
      {	  
      p_halo[j].p_ahf_part[n] = local_id;            /* only record the ID, starting from 0 */
      if(p_halo[j].ahf_halo.Mvir > ENVIR_MASS_CRIT)  /* subhalo with smaller mass inside the host will not change the flag */
        p_part[local_id].envir_flag = 0;
      }       
    }
}                    /* end construct_ahf_halo_part */

#ifndef NO_AHF_PROFILE
void construct_ahf_halo_profile(unsigned i, unsigned j)
{
  int n;
  
  alloc_ahf_halo_profile_array(j, p_halo[j].ahf_halo.nbins);	
	
  for(n=0; n<p_halo[j].ahf_halo.nbins; n++)
    {
   	fscanf(fp_ahf_profile[i], 
	          "%lf \
	           %u \
	           %lf \
	           %lf %lf \
	           %lf %lf %lf \
	           %lf %lf %lf \
	           %lf %lf \
	           %lf %lf %lf \
	           %lf %lf %lf \
	           %lf %lf %lf \
	           %lf %lf\n", 
	           &p_halo[j].p_ahf_profile[n].r_right, 
	           &p_halo[j].p_ahf_profile[n].npart, 
	           &p_halo[j].p_ahf_profile[n].M_in_r, 
	           &p_halo[j].p_ahf_profile[n].ovdens, &p_halo[j].p_ahf_profile[n].dens, 
	           &p_halo[j].p_ahf_profile[n].vcirc, &p_halo[j].p_ahf_profile[n].v_esc, &p_halo[j].p_ahf_profile[n].sigV, 
	           &p_halo[j].p_ahf_profile[n].L[0], &p_halo[j].p_ahf_profile[n].L[1], &p_halo[j].p_ahf_profile[n].L[2], 
	           &p_halo[j].p_ahf_profile[n].b_c, &p_halo[j].p_ahf_profile[n].a_c,
	           &p_halo[j].p_ahf_profile[n].vec_c[0], &p_halo[j].p_ahf_profile[n].vec_c[1], &p_halo[j].p_ahf_profile[n].vec_c[2],
	           &p_halo[j].p_ahf_profile[n].vec_b[0], &p_halo[j].p_ahf_profile[n].vec_b[1], &p_halo[j].p_ahf_profile[n].vec_b[2], 
	           &p_halo[j].p_ahf_profile[n].vec_a[0], &p_halo[j].p_ahf_profile[n].vec_a[1], &p_halo[j].p_ahf_profile[n].vec_a[2],
	           &p_halo[j].p_ahf_profile[n].Ekin, &p_halo[j].p_ahf_profile[n].Epot);			
    }
}                   /* end construct_ahf_halo_part */	
#endif

void construct_halo_profile(unsigned j)
{	
  alloc_halo_profile_array(j, PROFILE_BIN_NUM);
  init_halo_profile_array(j);
  cal_halo_profile(j);
  get_cNFW(j);	
}                  	/* end construct_halo_profile */


void construct_ahf_subhalo()
{
  unsigned total_host_num_here, i, j;           
  unsigned host_id, sub_num;
  char temp_a;	
  
  if(ahf_file_num == 1)                          /* lots of limitation these days */
    {	
	  fscanf(fp_ahf_subhalo[0], "%u\n", &total_host_num_here);
	  for(i=0; i<total_host_num_here; i++)
	    {
	    fscanf(fp_ahf_subhalo[0], "%u %u\n", &host_id, &sub_num);
	  
	    if(host_id < total_halo_num)                /* host ID here may be larger than total_halo_num -.-! */
	      {
	      if(p_halo[host_id].ahf_halo.sub_num != sub_num)
	        {
		      printf("Host halo %u, subnum from AHF halo file is %u, but from subhalo file is %u.\n", 
		             host_id, p_halo[host_id].ahf_halo.sub_num, sub_num); 
	        end_run("MATCHING SUBHALO DATA, ERROR 1");
	        }
	      alloc_halo_subhalo_array(host_id, sub_num);
	  
	      for(j=0; j<sub_num; j++)
	        {
	        fscanf(fp_ahf_subhalo[0], "%d", &p_halo[host_id].p_ahf_subhalo[j].sub_id);	      
	        if(p_halo[p_halo[host_id].p_ahf_subhalo[j].sub_id].ahf_halo.host_flag != host_id)
	          end_run("MATCHING SUBHALO DATA, ERROR 2");
          if(p_halo[p_halo[host_id].p_ahf_subhalo[j].sub_id].trust_flag)
            p_halo[host_id].trust_sub_num++;
          }
	      }    	
	    do{
         temp_a = fgetc(fp_ahf_subhalo[0]);
        }while(temp_a != 10);     /* read remains of the line */          
      }
#ifdef DEBUG       
  printf("Complete initialize sub-halo INFO.\n");
#endif       
    }
}                  /* end construct_ahf_subhalo */	

/* 
 * This function will determine whether the halo is suitable for analysis.
 * For normal simulation, it requires the total particle number larger than the RESOLUTED_HALO_PART_NUMBER,
 * while for ddm, the requirement is for the number of mother particles. 
 */ 

void init_halo_trust_flag(unsigned j)   
{
#ifndef DECAY_DARK_MATTER
  if(p_halo[j].ahf_halo.npart >= RESOLUTED_HALO_PART_NUM)
    p_halo[j].trust_flag = 1;
  else
    p_halo[j].trust_flag = 0;
#else
  if(p_halo[j].dl_decay.ddm_num >= RESOLUTED_HALO_PART_NUM)
    p_halo[j].trust_flag = 1;
  else
    p_halo[j].trust_flag = 0;
#endif		
}                   /* end init_halo_trust_flag */
