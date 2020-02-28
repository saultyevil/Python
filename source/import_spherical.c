
/***********************************************************/
/** @file  import_spherical.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief
 * General purpose routines for reading in an arbitrary wind model
 * in spherical coordinates.
 *
 * The basic data we need to read in are,
 *
 *     i r v_r mass_rho (and optionally T_r)
 *
 *  where,
 *
 *  * i are the element numbers (increasing outwards)
 *  * r is the radial coordinates
 *  * v_r is the velocity in the radial direction
 *  * mass_rho is the density in cgs units
 *  * T_r is the radiation temperature in Kelvin
 *
 * We assume that all of the physical quantities are centered, that is
 * we are assuming that we are giving mass_rho/T_r at the center of
 * a cell. However, r and v_r should be given at the edges of a cell.
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"
#include "import.h"


/**********************************************************/
/**
 * @brief      Read an arbitrary wind model in spherical symmetry.
 *
 * @param [in out] int     ndom       The domain number for the imported model
 * @param [in out] char *  filename   The file containing the model to import
 * @return                            Always returns 0
 *
 * @details
 *
 * This routine just reads in the data and stores it in arrays
 *
 * ### Notes ###
 *
 * The basic data we need to read in are,
 *
 *     i r v_r mass_rho (and optionally T_r)
 *
 *  where,
 *
 *  * i are the element numbers (increasing outwards)
 *  * r is the radial coordinates
 *  * v_r is the velocity in the radial direction
 *  * mass_rho is the density in cgs units
 *  * T_r is the radiation temperature in Kelvin
 *
 * We assume that all of the physical quantities are centered, that is
 * we are assuming that we are giving mass_rho/T_r at the center of
 * a cell. However, r and v_r should be given at the edges of a cell.
 *
 **********************************************************/

int
import_1d (ndom, filename)
     int ndom;
     char *filename;
{
  FILE *fptr;
  char line[LINELENGTH];
  int n, icell, ncell;
  double r, v_r, mass_rho, t_r;

  Log ("Reading a 1d model %s\n", filename);

  if ((fptr = fopen (filename, "r")) == NULL)
  {
    Error ("import_1d: No such file\n");
    Exit (0);
  }

  ncell = 0;
  while (fgets (line, LINELENGTH, fptr) != NULL)
  {
    n = sscanf (line, " %d %le %le %le %le", &icell, &r, &v_r, &mass_rho, &t_r);
    if (n < 4)
    {
      continue;
    }
    else
    {
      import_model_1d.element[ncell] = icell;
      import_model_1d.r[ncell] = r;
      import_model_1d.v_r[ncell] = v_r;
      import_model_1d.mass_rho[ncell] = mass_rho;
      if (n > 4)
      {
        import_model_1d.t_r[ncell] = t_r;
      }
      else
      {
        import_model_1d.t_r[ncell] = DEFAULT_IMPORT_TEMPERATURE;
      }

      ncell++;

      if (ncell > NDIM_MAX)
      {
        Error ("%s : %i : trying to read in more grid points than allowed (%i). Try changing NDIM_MAX and recompiling.\n", __FILE__,
               __LINE__, NDIM_MAX);
        Exit (1);
      }

    }
  }

  import_model_1d.ncell = ncell;

  /* Although much of the initialization of zdom can be postponed
   * one has to define mdim and ndim of zdom here, so that the correct
   * number of wind cells will be allocated */

  zdom[ndom].ndim2 = zdom[ndom].ndim = import_model_1d.ncell;
  zdom[ndom].mdim = 1;

  return (0);
}


/* The next section contains routines to make the grids for imported models.

   We make some assumptions here.  We assume that every cell is in the wind.
   and we assume that r refers to the inside edge of the cell.

 * */


/**********************************************************/
/**
 * @brief      Use the imported data to initialize various
 * portions of the Wind and Domain structures
 *
 *
 * @param [out] WindPtr  w   The entire wind
 * @param [in] int  ndom   The domain for the imported model
 * @return     Always returns 0
 *
 * @details
 * This routine initializes the portions of the wind structure
 * using the imported model, specifically those portions having
 * to do with positions.
 *
 * ### Notes ###
 *
 **********************************************************/

int
spherical_make_grid_import (w, ndom)
     WindPtr w;
     int ndom;
{

  int j, n;

  zdom[ndom].wind_rho_min = zdom[ndom].rho_min = 0;
  zdom[ndom].rmin = import_model_1d.r[0];
  zdom[ndom].wind_rho_max = zdom[ndom].zmax = zdom[ndom].rho_max = zdom[ndom].rmax = import_model_1d.r[import_model_1d.ncell - 1];
  zdom[ndom].wind_thetamin = zdom[ndom].wind_thetamax = 0.;

  for (j = 0; j < import_model_1d.ncell; j++)
  {
    n = j + zdom[ndom].nstart;
    w[n].r = import_model_1d.r[j];
    /* Put the radial velocity in v[0] */
    w[n].v[0] = import_model_1d.v_r[j];
  }

  /* Need to define the midpoints of the grid */

  for (j = 0; j < zdom[ndom].ndim; j++)
  {
    n = j + zdom[ndom].nstart;
    if (j < zdom[ndom].ndim - 1)
    {
      w[n].rcen = 0.5 * (w[n].r + w[n + 1].r);
    }
    else
    {
      w[n].rcen = w[n].r * 1.005;
    }
    w[n].x[1] = w[n].xcen[1] = 0.0;
    w[n].x[0] = w[n].x[2] = w[n].r * sin (PI / 4.);
    w[n].xcen[0] = w[n].xcen[2] = w[n].rcen * sin (PI / 4.);
  }

  /* Since we assume all of the cells are in the wind in a spherical wind
   * we can use the standard routine to finish everything off
   */

  spherical_wind_complete (ndom, w);

  return (0);
}


/* The next section calculates velocities.
 *
 * One could follow the zeus_hydro approach of getting those velocities from the original grid.
 * but for consistency with the 2d case we get it by interpolating on values in the cells
 *
 *
 */


/**********************************************************/
/**
 * @brief      The velocity at any position in an imported spherical
 * model
 *
 *
 * @param [in] int  ndom   The domain of the imported model
 * @param [in] double *  x   A position (3d)
 * @param [out] double *  v   The velocity at x
 * @return     The speed at x
 *
 * @details
 * This routine interpolates on the values read in for the
 * imported model to give one a velocity
 *
 *
 * ### Notes ###
 * Note that v_r is stored in v_0
 *
 * Not also that In practice this routine is only used to initallize v in
 * wind structure.  This is consistent with the way velocities
 * are treated throughout Python
 *
 **********************************************************/

double
velocity_1d (ndom, x, v)
     int ndom;
     double *x, *v;
{
  double speed;
  double r;
  int nelem, nn, nnn[4];
  double frac[4];

  r = length (x);
  coord_fraction (ndom, 0, x, nnn, frac, &nelem);

  speed = 0;
  for (nn = 0; nn < nelem; nn++)
  {
    speed += wmain[zdom[ndom].nstart + nnn[nn]].v[0] * frac[nn];
  }

  v[0] = x[0] / r * speed;
  v[1] = x[1] / r * speed;
  v[2] = x[2] / r * speed;

  return (speed);
}




/**********************************************************/
/**
 * @brief      Get the density for an imported spherical model at x
 *
 * @param [in] int  ndom   The domain for the imported model
 * @param [in] double *  x   A position (3d)
 * @return     The density in cgs units is returned
 *
 * @details
 * This routine finds rho from the imported model
 * at a position x.  The routine does not interpolate rho, but
 * simply locates the cell associated with x
 *
 *
 * ### Notes ###
 * This routine is really only used to intialize rho in the
 * Plasma structure.  In reality, once the Plasma structure is
 * initialized we always interpolate within the plasma structure
 * and do not access the original data.
 *
 * This routine is peculiar because it depends on the assumption
 * that x is the center of the cell.
 *
 **********************************************************/

double
rho_1d (ndom, x)
     int ndom;
     double *x;
{
  double rho = 0;
  double r;
  int n;

  r = length (x);

  n = 0;
  while (r >= import_model_1d.r[n] && n < import_model_1d.ncell)
  {
    n++;
  }
  n--;

  if (n < import_model_1d.ncell)
  {
    rho = import_model_1d.mass_rho[n];
  }
  else
  {
    rho = import_model_1d.mass_rho[import_model_1d.ncell - 1];
  }

  return (rho);
}
