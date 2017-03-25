/* Only valid for AHF version >=1.0 
 * and dark matter only simulations 
 * with total number of particles less than 2^32.
 * Units:
 *  1 length = 1 kpc/h
 *  1 mass   = 1 M_sun/h
 *  1 velocity = 1 km/s
 */

#define AHF_INPUT_PATH               "/users/s0902248/Install/ahf-v1.0-049/Result/Old-cosmo-param/Decay-test-21/z-0/" 
#define GADGET_INPUT_PATH            "/users/s0902248/Install/Gadget2/Decay-dark-matter-cosmology/Result/Old-cosmo-param/Decay-test-21/"
#define AHF_FILE_BASE                "LTic_decay_test_21"
#define GADGET_FILE_BASE             "snapshot_009"    
#define OUTPUT_PATH                  "/users/s0902248/Install/Decay-halo-info/v-0.2.5/Result/Old-cosmo-param/Decay-test-21/z-0/" 

#define RESOLUTED_HALO_PART_NUM       50         /* number of undecayed particle in a halo to be considered as trustworthy */
#define PROFILE_BIN_NUM               30         /* the profile bin number in log scale, physical length */

#define ENVIR_MASS_CRIT               1e8        /* halo with mass less than this value
                                                  * will be classified as enviroment, in unit of M_sun/h */


