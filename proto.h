#ifndef ALLVAR_H
 #include "allvars.h"
#endif

void detect_and_link_gadget_file();
void load_gedget_part();
void transfer_to_peculiar_vel(float *vel);
void check_sorted_part_array();
void sort_gadget_part(long int l, long int r);
void swap_gadget_part(PPT_PART_INFO a , PPT_PART_INFO b);
void copy_float_array(float *a, float *b, int n);
void check_gadget_file(char *file_name);

void detect_and_link_ahf_file();
void effective_rewind(FILE **p);

void init_halo_array();
void construct_ahf_halo_info(unsigned i, unsigned j);
void construct_ahf_halo_part(unsigned i, unsigned j);
#ifndef NO_AHF_PROFILE
void construct_ahf_halo_profile(unsigned i, unsigned j);
#endif
void construct_halo_profile(unsigned j);
void construct_ahf_subhalo();
void init_halo_trust_flag(unsigned j);

double cNFW_from_vel(double v_ratio);
double function_c(double c, void *param);
double d_function_c(double c, void *param);
void fd_function_c(double c, void *param, double *f, double *df);

void init_halo_profile_array(unsigned j);
void cal_halo_profile(unsigned j);
void get_cNFW(unsigned j);
double distance_to_center(float *a, double *c);

void init_all();

void alloc_gadget_file_array();
void alloc_ahf_file_arrays();
void alloc_halo_array();
void alloc_part_array();
void alloc_halo_part_array(unsigned j);
#ifndef NO_AHF_PROFILE
void alloc_ahf_halo_profile_array(unsigned j, int nbins);
#endif
void alloc_halo_profile_array(unsigned j, int nbins);
void alloc_halo_subhalo_array(unsigned host_id, unsigned sub_num);
void free_all_arrays();

void out_put();

void close_all_files();

void state(char *s);
void end_run(char *s);