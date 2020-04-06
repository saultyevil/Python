/* ************************************************************************** */
/**
 * @file  import.h
 * @author EJP
 * @date   Feb 2019
 *
 * @brief    Global structures for containing imported models
 *
 * @details
 *
 * The constants and structures which require global scope throughout the wind
 * import functions are defined here. This includes the structure which will
 * temporarily hold the imported model data.
 *
 * ************************************************************************** */

#define NDIM_MAX2D NDIM_MAX * NDIM_MAX      // Maximum dimensions for 2D importing

/*
 * The following definitions are used to try and improve the readability for
 * some of the code used to read in 1D and 2D grids from an input file. They
 * refer to how many columns are in the data file and what will be read in.
 */

#define READ_NO_TEMP_1D          4
#define READ_ELECTRON_TEMP_1D    READ_NO_TEMP_1D + 1
#define READ_BOTH_TEMP_1D        READ_NO_TEMP_1D + 2

#define READ_NO_TEMP_2D          9
#define READ_ELECTRON_TEMP_2D    READ_NO_TEMP_2D + 1
#define READ_BOTH_TEMP_2D        READ_NO_TEMP_2D + 2

/*
 * The following structure will contain all of the required information for
 * creating a wind grid using an imported model.
 */

struct
{
  int init_temperature;               // initialise to t.wind.init if TRUE
  int ncell;                          // the total number of cells read in
  int ndim, mdim;                     // the number of coordinates in the n and m dimensions
  int *i, *j;                         // the i (row) and j (column) elements
  int *inwind;                        // flag for the cell being inwind or not inwind
  double *x, *z, *r, *theta;          // the x/r or z/theta coordinates of the grid in cgs units
  double *v_x, *v_y, *v_z;            // the velocity in Cartesian coordinates in cgs units
  double *v_r;                        // the radial velocity in cgs units
  double *mass_rho;                   // the mass density in cgs units
  double *t_e, *t_r;                  // the electron and radiation temperature in Kelvin
  double *wind_x, *wind_z;            // the wind grid coordinates
  double *wind_midx, *wind_midz;      // the wind grid mid points
} *imported_model;             // MaxDom is defined in python.h and as such import.h has to be included after
