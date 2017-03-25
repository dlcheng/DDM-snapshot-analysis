#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"


void alloc_gadget_file_array()
{
  fp_gadget = (FILE **)calloc(gadget_file_num, sizeof(FILE *));		
}                                    /* end alloc_gadget_file_array */

void alloc_ahf_file_arrays()
{
  fp_ahf_halo = (FILE **)calloc(ahf_file_num, sizeof(FILE *));
#ifndef NO_AHF_PROFILE  
  fp_ahf_profile = (FILE **)calloc(ahf_file_num, sizeof(FILE *));
#endif  
  fp_ahf_part = (FILE **)calloc(ahf_file_num, sizeof(FILE *));
  fp_ahf_subhalo = (FILE **)calloc(ahf_file_num, sizeof(FILE *));

  p_file_info= (PPT_FILE_INFO)calloc(ahf_file_num, sizeof(FILE_INFO));	
		
}                                    /* end alloc_ahf_file_arrays */	

void alloc_halo_array()
{
  if((p_halo = (PPT_HALO_INFO)calloc(total_halo_num, sizeof(HALO_INFO))) == NULL)	
	  end_run("ALLOCATING HALO ARRAY");	
}                                    /* end alloc_haloc_arrays */	

void alloc_part_array()
{
  if((p_part = (PPT_PART_INFO)calloc(total_part_num, sizeof(PART_INFO))) == NULL)	
	  end_run("ALLOCATING GADGET PARTICLE ARRAY");
}	                                 /* alloc_part_array */

void alloc_halo_part_array(unsigned j)
{
  if((p_halo[j].p_ahf_part = (unsigned int*)calloc(p_halo[j].ahf_halo.npart, sizeof(unsigned int))) == NULL)
    end_run("ALLOCATING MEMORY FOR HALO PARTICLE ARRAY");
}                                    /* end alloc_halo_part_array */

#ifndef NO_AHF_PROFILE
void alloc_ahf_halo_profile_array(unsigned j, int nbins)
{
  if((p_halo[j].p_ahf_profile = (PPT_AHF_PROFILE_INFO)calloc(nbins, sizeof(AHF_PROFILE_INFO))) == NULL)
	  end_run("ALLOCATING MEMORY FOR AHF HALO PROFILE ARRAY");
}                                    /* end alloc_ahf_halo_profile_array */	
#endif

void alloc_halo_profile_array(unsigned j, int nbins)
{
  if((p_halo[j].p_halo_profile = (PPT_HALO_PROFILE_INFO)calloc(nbins, sizeof(HALO_PROFILE_INFO))) == NULL)
	  end_run("ALLOCATING MEMORY FOR HALO PROFILE ARRAY");
}                                    /* end alloc_ahf_halo_profile_array */	

void alloc_halo_subhalo_array(unsigned host_id, unsigned sub_num)
{
  if((p_halo[host_id].p_ahf_subhalo = (PPT_AHF_SUBHALO_INFO)calloc(sub_num, sizeof(AHF_SUBHALO_INFO))) == NULL)	
	  end_run("ALLOCATING MEMORY FOR SUBHALO");	
}                                    /* end alloc_halo_subhalo_array */	

void free_all_arrays()
{
  unsigned j;
  
  free(fp_ahf_halo);
  free(fp_ahf_part);
#ifndef NO_AHF_PROFILE  
  free(fp_ahf_profile); 
#endif  
  free(fp_ahf_subhalo);    
  free(p_file_info);
  free(p_part);
  
  for(j=0; j<total_halo_num; j++)
    {
	  free(p_halo[j].p_ahf_part);	
#ifndef NO_AHF_PROFILE
	  free(p_halo[j].p_ahf_profile);
#endif	
    free(p_halo[j].p_halo_profile);	    
    free(p_halo[j].p_ahf_subhalo);	    
    }
   
  free(p_halo);   
  state("Free all arrays.");
}                                    /* end free_all_arrays */	
