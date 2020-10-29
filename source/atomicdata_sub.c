
/***********************************************************/
/** @file  atomicdata_sub.c
 * @author ksl
 * @date   March, 2018
 *
 * @brief  Subroutines related directly to the atomic data 
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"
#include "log.h"
// If routines are added cproto > atomic_proto.h should be run
#include "atomic_proto.h"

#define LINELENGTH 400
#define MAXWORDS    20

/**********************************************************/
/**
 * @brief      Write out the atomic data to a file   
 *
 * @return     Always returns 0
 *
 * @details
 *
 * ### Notes ###
 *
 * Note that there are a number of "variables" such as 
 * nelements, and nions, that are global variables which
 * are part of atomic.h
 **********************************************************/

int
atomicdata2file ()
{
  FILE *fptr;
  int nelem;
  int n;

  if ((fptr = fopen ("atomic_data.out.txt", "w")) == NULL)
  {
    Error ("get_atomic data:  Could not open atomic_data.out.txt\n");
    exit (0);
  }

  fprintf (fptr, "# This file contains data which was read in by get_atomicdata\n");

  fprintf (fptr, "# Data of %d elements, %d ions, and %d levels %d lines\n", nelements, nions, nlevels, nlines);


  /* Write the element array */
  fprintf (fptr, "# Element Data:\n");
  for (nelem = 0; nelem < nelements; nelem++)
  {
    fprintf (fptr, "Element %2d %5s firstion %4d nions %2d\n", nelem, ele[nelem].name, ele[nelem].firstion, ele[nelem].nions);
  }

  /* Write the ion array */
  fprintf (fptr, "# Ion data:\n");
  for (n = 0; n < nions; n++)
  {
    fprintf (fptr,
             "Ion %3d z %3d istate %3d firstlevel %5d nlevels %3d potential %8.3g  phot_info %3d\n",
             n, ion[n].z, ion[n].istate, ion[n].firstlevel, ion[n].nlevels, ion[n].ip / EV2ERGS, ion[n].phot_info);
  }

  /* Write the excitation level data */
  fprintf (fptr, "# Excitation levels: There are %d levels\n", nlevels);
  for (n = 0; n < nlevels; n++)
    fprintf (fptr, "Level n %5d z %2d istate %2d q %.1f g %4.0f ex %8.3g  bb %2d %2d bf %2d %2d\n", n,
             config[n].z, config[n].istate, config[n].q_num, config[n].g, config[n].ex,
             config[n].n_bbu_jump, config[n].n_bbd_jump, config[n].n_bfu_jump, config[n].n_bfd_jump);

  /* Write the resonance line data to the file */

  fprintf (fptr, "# Line data: There are %d lines\n", nlines);

  for (n = 0; n < nlines; n++)
  {
    fprintf (fptr, "Line n %5d ion_no %3d z %2d ion %2d freq %8.1e f %6.3f macro %2d coll %4d up_down %2d %2d\n",
             n, line[n].nion, line[n].z, line[n].istate, line[n].freq, line[n].f, line[n].macro_info, line[n].coll_index,
             line[n].down_index, line[n].up_index);
  }

  /* Write the photoionization data  */
  fprintf (fptr, "# Photoionization data: There are %d edges\n", ntop_phot + nxphot);
  for (n = 0; n < ntop_phot + nxphot; n++)
  {
    fprintf (fptr, "Phot n %5d z %2d istate %3d sigma %8.2e freq[0] %8.2e nlev %5d uplev %2d macro %2d  %2d %2d use %2d\n",
             n, phot_top[n].z, phot_top[n].istate, phot_top[n].x[0], phot_top[n].freq[0],
             phot_top[n].nlev, phot_top[n].uplev, phot_top[n].macro_info, phot_top[n].down_index, phot_top[n].up_index, phot_top[n].use);
  }

  /* Write the ground fraction data to the file */
  fprintf (fptr, "# Ground frac data (just first and last fracs here as a check):\n");

  for (n = 0; n < NIONS; n++)
  {
    fprintf (fptr, "Ground %3d %3d %6.3f %6.3f\n", ground_frac[n].z, ground_frac[n].istate, ground_frac[n].frac[0],
             ground_frac[n].frac[19]);
  }

  fclose (fptr);

  return (0);
}


/**********************************************************/
/**
 * @brief      sort the lines into frequency order
 *
 * @return     Always returns 0
 *
 * @details
 *
 * ### Notes ###
 * All use a Numerical recipes routine indexx for this
 *
 **********************************************************/

int
index_lines ()
{
  float *freqs, foo;
  int *index, ioo;
  int n;
  void indexx ();

  /* Allocate memory for some modestly large arrays */
  freqs = calloc (sizeof (foo), NLINES + 2);
  index = calloc (sizeof (ioo), NLINES + 2);

  freqs[0] = 0;
  for (n = 0; n < nlines; n++)
    freqs[n + 1] = line[n].freq;        /* So filled matrix
                                           elements run from 1 to nlines */

  indexx (nlines, freqs, index);        /* Note that this math recipes routine
                                           expects arrays to run from 1 to nlines inclusive */

  /* The for loop indices are complicated by the numerical recipes routine,
     which is a simple translation of a fortran routine.
     Specifically, index array elements 1 to nlines are now filled,
     and the numbers run from 1 to nlines, but the
     pointer array is only filled from elements 0 to nlines -1 */

  /* SS - adding quantity "where_in_list" to line structure so that it is easy to from emission
     in recombination line to correct place in line list. */


  for (n = 0; n < nlines; n++)
  {
    lin_ptr[n] = &line[index[n + 1] - 1];
    line[index[n + 1] - 1].where_in_list = n;
  }

  /* Free the memory for the arrays */
  free (freqs);
  free (index);

  return (0);
}


/**********************************************************/
/**
 * @brief      Index the topbase photoionzation crossections by frequency
 *
 *
 * @return     Always returns 0
 *
 * @details
 *
 * The results are stored in phot_top_ptr
 *
 * ### Notes ###
 * Adapted from index_lines as part to topbase
 *             addition to python
 *
 **********************************************************/

int
index_phot_top ()
{
  float *freqs, foo;
  int *index, ioo;
  int n;
  void indexx ();

  /* Allocate memory for some modestly large arrays */
  freqs = calloc (sizeof (foo), ntop_phot + nxphot + 2);
  index = calloc (sizeof (ioo), ntop_phot + nxphot + 2);

  freqs[0] = 0;
  for (n = 0; n < ntop_phot + nxphot; n++)
    freqs[n + 1] = phot_top[n].freq[0]; /* So filled matrix
                                           elements run from 1 to ntop_phot */

  indexx (ntop_phot + nxphot, freqs, index);    /* Note that this math recipes routine
                                                   expects arrays to run from 1 to nlines inclusive */

  /* The for loop indices are complicated by the numerical recipes routine,
     which is a simple translation of a fortran routine.
     Specifically, index array elements 1 to nlines are now filled,
     and the numbers run from 1 to nlines, but the
     pointer array is only filled from elements 0 to nlines -1 */

  for (n = 0; n < ntop_phot + nxphot; n++)
  {
    phot_top_ptr[n] = &phot_top[index[n + 1] - 1];
  }

  /* Free the memory for the arrays */
  free (freqs);
  free (index);

  return (0);

}


/**********************************************************/
/**
 * @brief      Index inner shell xsections in frequency order
 *
 * @return     Alwasy returns 0
 *
 * @details
 * The rusults are stored in inner_cross_ptr
 *
 * ### Notes ###
 * ??? NOTES ???
 *
 **********************************************************/

int
index_inner_cross ()
{
  float *freqs, foo;
  int *index, ioo;
  int n;
  void indexx ();

  /* Allocate memory for some modestly large arrays */
  freqs = calloc (sizeof (foo), n_inner_tot + 2);
  index = calloc (sizeof (ioo), n_inner_tot + 2);

  freqs[0] = 0;
  for (n = 0; n < n_inner_tot; n++)
    freqs[n + 1] = inner_cross[n].freq[0];      /* So filled matrix
                                                   elements run from 1 to ntop_phot */

  indexx (n_inner_tot, freqs, index);   /* Note that this math recipes routine
                                           expects arrays to run from 1 to nlines inclusive */

  /* The for loop indices are complicated by the numerical recipes routine,
     which is a simple translation of a fortran routine.
     Specifically, index array elements 1 to nlines are now filled,
     and the numbers run from 1 to nlines, but the
     pointer array is only filled from elements 0 to nlines -1 */

  for (n = 0; n < n_inner_tot; n++)
  {
    inner_cross_ptr[n] = &inner_cross[index[n + 1] - 1];
  }

  /* Free the memory for the arrays */
  free (freqs);
  free (index);

  return (0);

}


/* Numerical recipes routine used by index_lines which in turn is used by get_atomic_data */


/**********************************************************/
/**
 * @brief      Numerical recipes routine used by various routines to sort the atomic data
 *
 * @param [in] int  n   The dimension of the array that needs sorting
 * @param [in] float  arrin[]   The array to sort
 * @param [out] int  indx[]   The array that contains poionters to the sorted array
 * @return     N/A
 *
 * @details
 *
 * ### Notes ###
 * See nuerical recipes
 *
 * @bug Should probably replaces this with the equivalent gsl routine
 *
 **********************************************************/

void
indexx (n, arrin, indx)
     int n, indx[];
     float arrin[];
{
  int l, j, ir, indxt, i;
  float q;
/* NSH 1408 - This routine fails in the very odd circumstance that n=1 so we now do a little test here */
  if (n < 2)
  {
    Log_silent ("Nothing for indexx to do! Only one element\n");
    indx[0] = 0;
    indx[1] = 1;                /* NSH 1707 - We still need to populate the array */
    return;
  }

  for (j = 1; j <= n; j++)
    indx[j] = j;
  l = (n >> 1) + 1;
  ir = n;
  for (;;)
  {
    if (l > 1)
      q = arrin[(indxt = indx[--l])];
    else
    {
      q = arrin[(indxt = indx[ir])];
      indx[ir] = indx[1];
      if (--ir == 1)
      {
        indx[1] = indxt;
        return;
      }
    }
    i = l;
    j = l << 1;
    while (j <= ir)
    {
      if (j < ir && arrin[indx[j]] < arrin[indx[j + 1]])
        j++;
      if (q < arrin[indx[j]])
      {
        indx[i] = indx[j];
        j += (i = j);
      }
      else
        j = ir + 1;
    }
    indx[i] = indxt;
  }
}




/**********************************************************/
/**
 * @brief      uses freqmin and freqmax to set limits on which lines in the
 * frequency ordered list of lines need to be
 * considered when, for example,  calculating which lines have Sobolev surfaces
 * along a path. 
 *
 *
 * @param [in] double  freqmin   The minimum frequency we are interested in
 * @param [in] double  freqmax   The maximum frequency we are interested in
 * @return     the number of lines that are potentially in resonance.
 *
 * limit_lines sets the external variables nline_min and nline_max for use with
 * other routines
 *
 *
 * If limit_lines
 * 	returns 0 there are no lines of interest and one does not need to worry about any
 * 	resonances at this frequency.  If limit_lines returns a number greater than 0, then
 * 	the lines of interest are defined by nline_min and nline_max (inclusive) in atomic.h
 * 	nline_delt is also set which is the number of lines that are in the specified range
 *
 * @details
 * limit_lines  define the lines that are close to a given frequency.  The degree of closeness
 * 	is defined by v. The routine must be used in conjuction with get_atomic_data which
 * 	will have created an ordered list of the lines.
 *
 * ### Notes ###
 * Limit_lines needs to be used somewhat carefully.  Carefully means checking the
 * 	return value of limit_lines or equivalently nline_delt=0.  If nline_delt=0 then there
 * 	were no lines in the region of interest.  Assuming there were lines in thte range,
 * 	one must sum over lines from nline_min to nline_max inclusive.
 *
 * 	One might wonder why nline_max is not set to one larger than the last line which
 * 	is in range.  This is because depending on how the velocity is trending you may
 * 	want to sum from the highest frequency line to the lowest.
 *
 *
 **********************************************************/

int
limit_lines (freqmin, freqmax)
     double freqmin, freqmax;
{

  int nmin, nmax, n;
  double f;


  if (freqmin > lin_ptr[nlines - 1]->freq || freqmax < lin_ptr[0]->freq)
  {
    nline_min = 0;
    nline_max = 0;
    nline_delt = 0;
    return (0);
  }

  f = freqmin;

  nmin = 0;
  nmax = nlines - 1;
  n = (nmin + nmax) >> 1;       // Compute a midpoint >> is a bitwise right shift

  while (n != nmin)
  {
    if (lin_ptr[n]->freq < f)
      nmin = n;
    if (lin_ptr[n]->freq >= f)
      nmax = n;
    n = (nmin + nmax) >> 1;     // Compute a midpoint >> is a bitwise right shift
  }

  nline_min = nmin;

  f = freqmax;
  nmin = 0;
  nmax = nlines - 1;
  n = (nmin + nmax) >> 1;       // Compute a midpoint >> is a bitwise right shift

  while (n != nmin)
  {
    if (lin_ptr[n]->freq <= f)
      nmin = n;
    if (lin_ptr[n]->freq > f)
      nmax = n;
    n = (nmin + nmax) >> 1;     // Compute a midpoint >> is a bitwise right shift
  }

  nline_max = nmax;


  return (nline_delt = nline_max - nline_min + 1);
}




/**********************************************************/
/**
 * @brief      Perform sanity checks on xsection data
 *
 * @return     Alway returns 0; errors are written to the log file
 *
 * @details
 * Results are written to the log file
 *
 * ### Notes ###
 *
 * @bug  A careful look at this routine is warranted as get_atomic_data has changed over the years.
 *
 **********************************************************/

int
check_xsections ()
{
  int nion, n;

  for (n = 0; n < nphot_total; n++)
  {
    nion = phot_top[n].nion;
    if (ion[nion].phot_info == 1)
      Debug
        ("Topbase Ion %i Z %i istate %i nground %i ilv %i ntop %i f0 %8.4e IP %8.4e\n",
         nion, ion[nion].z, ion[nion].istate, ion[nion].ntop_ground, phot_top[n].nlev, ion[nion].ntop, phot_top[n].freq[0], ion[nion].ip);
    else if (ion[nion].phot_info == 0)
      Debug ("Vfky Ion %i Z %i istate %i nground %i f0 %8.4e IP %8.4e\n",
             nion, ion[nion].z, ion[nion].istate, ion[nion].nxphot, phot_top[n].freq[0], ion[nion].ip);

    /* some simple checks -- could be made more robust */
    if (ion[nion].n_lte_max == 0 && ion[nion].phot_info == 1)
    {
      Error
        ("get_atomicdata: not tracking levels for ion %i z %i istate %i, yet marked as topbase xsection!\n",
         nion, ion[nion].z, ion[nion].istate);
    }
    if (ion[nion].phot_info != 1 && ion[nion].macro_info)
    {
      Error
        ("get_atomicdata: macro atom but no topbase xsection! ion %i z %i istate %i, yet marked as topbase xsection!\n",
         nion, ion[nion].z, ion[nion].istate);
    }

    if (ion[nion].macro_info != phot_top[n].macro_info)
    {
      Error ("Macro info for ion doesn't match with xsection!! Could be serious.\n");
    }
  }

  return 0;
}




/*

   q21 calculates and returns the collisional de-excitation coefficient q21

   Notes:  This uses a simple approximation for the effective_collision_strength omega.
   ?? I am not sure where this came from at present and am indeed no sure that
   omega is omega(2->1) as it should be.  This should be carefully checked.??
   This has no been improved, and wherever possible we use tabulated collision
   stength data from Chianti (after Burgess and Tully 1992)

   c21=n_e * q21
   q12 = g_2/g_1 q21 exp(-h nu / kT )

   History:
   98aug        ksl     Recoded from other routines so that one always calculates q21 in
			the same way.
	01nov	ksl	Add tracking mechanism so if called with same conditions it
			returns without recalculation
	12oct	nsh	Added, then commented out approximate gaunt factor given in
			hazy 2.
	17jan	nsh Added code to use the actual collision strength date from chianti

 */


/// (8*PI)/(sqrt(3) *nu_1Rydberg
#define ECS_CONSTANT 4.773691e16

struct lines *q21_line_ptr;
double q21_a, q21_t_old;


/**********************************************************/
/**
 * @brief      calculates and returns the collisional de-excitation coefficient q21
 *
 * @param [in] struct lines *  line_ptr   A single line
 * @param [in] double  t   The temperature at wich to calculate the coefficient
 * @return    The collisional de-excitiation coefficient
 *
 * @details
 * This routine uses stored coefficients if we have them, or the Van Regemorter
 * approximation if we do not.
 *
 * ### Notes ###
 * the relevant paper to consult here is Van Regemorter 1962 (ApJ 136 906). We use an effective gaunt
 * factor to calculate collision strengths. There is one regime in which kt < hnu. 
 * 
************************************************************/

double
q21 (line_ptr, t)
     struct lines *line_ptr;
     double t;
{
  double gaunt;
  double omega;
  double u0;
  double upsilon ();


  if (q21_line_ptr != line_ptr || t != q21_t_old)
  {


    u0 = (BOLTZMANN * t) / (PLANCK * line_ptr->freq);

    if (line_ptr->istate == 1 && u0 < 2)        // neutrals at low energy. Used 2 to give continuous function.
      gaunt = u0 / 10.0;
    else                        // low energy electrons, positive ions
      gaunt = 0.2;



    if (line_ptr->coll_index < 0)       //if we do not have a collision strength for this line use the g-bar formulation
    {
      omega = ECS_CONSTANT * line_ptr->gl * gaunt * line_ptr->f / line_ptr->freq;
    }
    else                        //otherwise use the collision strength directly. NB what we call omega, most people including hazy call upsilon.
    {
      omega = upsilon (line_ptr->coll_index, u0);
    }


    q21_a = 8.629e-6 / (sqrt (t) * line_ptr->gu) * omega;
    q21_t_old = t;
  }

  return (q21_a);
}


/**********************************************************/
/**
 * @brief      Calculate the collisional excitation coefficient for a line
 *
 * @param [in] struct lines *  line_ptr   The line of interest
 * @param [in] double  t   The temperature of interest
 * @return     The collisional excitation coeffient
 *
 * @details
 * The routine calls q21 and uses detailed balance to derive q12
 *
 * ### Notes ###
 *
 **********************************************************/

double
q12 (line_ptr, t)
     struct lines *line_ptr;
     double t;
{
  double x;
  double q21 ();
  double exp ();

  x = line_ptr->gu / line_ptr->gl * q21 (line_ptr, t) * exp (-H_OVER_K * line_ptr->freq / t);

  return (x);
}


/*
   a21 alculates and returns the Einstein A coefficient
   History:
   98aug        ksl     Coded and debugged
   99jan        ksl Modified so would shortcircuit calculation if
   called multiple times for same a
 */
#define A21_CONSTANT 7.429297e-22       // 8 * PI * PI * E * E / (MELEC * C * C * C)

struct lines *a21_line_ptr;
double a21_a;


/**********************************************************/
/**
 * @brief      Calculate the Einstein A coefficient for a line
 *
 * @param [in out] struct lines *  line_ptr   The structure that describes a single line
 * @return     The Einstein A coefficient
 *
 * @details
 * Using the oscillator strecnth and the multiplicity (which are
 * read in) calculate A
 *
 * ### Notes ###
 *
 **********************************************************/

double
a21 (line_ptr)
     struct lines *line_ptr;
{
  double freq;

  if (a21_line_ptr != line_ptr)
  {
    freq = line_ptr->freq;
    a21_a = A21_CONSTANT * line_ptr->gl / line_ptr->gu * freq * freq * line_ptr->f;
    a21_line_ptr = line_ptr;
  }

  return (a21_a);
}


/**********************************************************/
/**
 * @brief      calculates the thermally averaged collision strength thermally excited line emission.
 *
 * @param [in] int  n_coll   the index of the collision strength record we are working with
 * @param [in] double  u0  - kT_e/hnu - where nu is the transition frequency for the line of interest and T_e is the electron temperature
 * @return     upsilon - the thermally averaged collision strength for a given line at a given temp
 *
 * @details
 *
 * ### Notes ###
 * It uses data extracted from Chianti stored in coll_stren.
 * The paper to consult is Burgess and Tully A&A 254,436 (1992).
 * u0 is the ratio of Boltzmans constant times the electron temperature in the cell
 * divided by Plancks constant times the frequency of the line under analysis.
 *
 *
 **********************************************************/

double
upsilon (n_coll, u0)
     int n_coll;
     double u0;
{
  double x;                     //The scaled temperature
  double y;                     //The scaled collision strength
  double upsilon;               //The actual collision strength
  int linterp ();

  /* first we compute x. This is the "reduced temperature" from
     Burgess & Tully 1992. */
  x = 0.0;
  if (coll_stren[n_coll].type == 1 || coll_stren[n_coll].type == 4)
  {
    x = 1. - (log (coll_stren[n_coll].scaling_param) / log (u0 + coll_stren[n_coll].scaling_param));
  }
  else if (coll_stren[n_coll].type == 2 || coll_stren[n_coll].type == 3)
  {
    x = u0 / (u0 + coll_stren[n_coll].scaling_param);
  }
  else
  {
    Error ("upsilon - coll_stren %i has no type %g\n", coll_stren[n_coll].type);
    exit (0);
  }


  /* we now compute y from the interpolation formulae
     y is the reduced upsilon from Burgess & Tully 1992. */
  linterp (x, coll_stren[n_coll].sct, coll_stren[n_coll].scups, coll_stren[n_coll].n_points, &y, 0);

  /*  now we extract upsilon from y  - there are four different parametrisations */

  upsilon = 0.0;
  if (coll_stren[n_coll].type == 1)
  {
    upsilon = y * (log (u0 + exp (1)));
  }
  else if (coll_stren[n_coll].type == 2)
  {
    upsilon = y;
  }
  else if (coll_stren[n_coll].type == 3)
  {
    upsilon = y / (u0 + 1);
  }
  else if (coll_stren[n_coll].type == 4)
  {
    upsilon = y * (log (u0 + coll_stren[n_coll].scaling_param));
  }
  else
  {
    Error ("upsilon - coll_stren %i has no type %g\n", coll_stren[n_coll].type);
    exit (0);
  }
  return (upsilon);
}






/**********************************************************/
/**
 * @brief Perform linear/logarithmic interpolation of an array
 *
 * @param [in] double  value   The value used for interpoaltion
 * @param [in] double  array[]   An array containing a set of asceding values
 * @param [in] int  npts   The size of the array
 * @param [out] int *  ival  The lower index to use in the interpolation
 * of the array which need to be used to interpolate on
 * @param [out] double *  f   The fraction of the upper point to use in the interpolation
 * @param [in] int  mode  A switch to choose linear(0)  or lograrithmic(1) interpolation
 * @return     Usually returns 0, but returns -1 if the input value is less
 * than the first elememnt in the array, and 1 if it is greater than the
 * last element int he array.  In either of these cases, the fractions are
 * set up only to access one array element
 *
 * @details
 *
 *
 * Typically one has two parallel arrays, one containing a set of values,
 * in the original case frequencies, and another containing some function
 * of those values.  This routine finds the fraction of the function at
 * two point in the data array to interpolate.
 *
 * The values for the input array can be interpolated linear or logarithmically.
 * that is the fraction that is returned are based on the values in the
 * array or the logarithm of them.
 *
 *
 * The routine uses bisection
 *
 *
 *
 * ### Notes ###
 *
 * fraction is a utility written to speed up the search for the bracketing
 * topbase photionization x-sections, but in principle it should work in
 * other situations within python (which at one time at least contributed
 * significantly to the runtime for Python).
 *
 * Today, fraction is called directly in the routines that are used to set up coordinate
 * systems, and indirectly through linterp (below) which should be inspected
 * to see how the logarithmic interpolation is supposed to work.
 *
 * The routine is similar to the numerical * recipes routine locate.
 * There may be a gsl routine as well.  The routine
 * should probably be replaced.
 *
 *
 **********************************************************/

int
fraction (value, array, npts, ival, f, mode)
     double array[];            // The array in we want to search
     int npts, *ival;           // ival is the lower point
     double value;              // The value we want to index
     double *f;                 // The fractional "distance" to the next point in the array
     int mode;                  // 0 = compute in linear space, 1=compute in log space
{
  int imin, imax, ihalf;

  if (value < array[0])
  {
    *ival = 0;
    *f = 0.0;
    return (-1);
  }

  imax = npts - 1;
  if (value > array[imax])
  {
    *ival = npts - 2;
    *f = 1.0;
    return (1);
  }



  imin = 0;

/* In what follows, there is a specific issue as to how to treat
the situation where the value is exactly on an array element. Right
now this is set to try to identify the array element below the one
on which the value sits, and to set the fraction to 1.  This was
to reflect the behavior of the search routine in where_in_grid. */

  while (imax - imin > 1)
  {
    ihalf = (imin + imax) >> 1; // Compute a midpoint >> is a bitwise right shift
    if (value > array[ihalf])
    {
      imin = ihalf;
    }
    else
      imax = ihalf;
  }

// So array[imin] just <= value

  if (mode == 0)
    *f = (value - array[imin]) / (array[imax] - array[imin]);   //linear interpolation
  else if (mode == 1)
    *f = (log (value) - log (array[imin])) / (log (array[imax]) - log (array[imin]));   //log interpolation
  else
  {
    Error ("Fraction - unknown mode %i\n", mode);
    exit (0);
    return (0);
  }

  *ival = imin;

  return (0);
}







/**********************************************************/
/**
 * @brief      Perform a linear interpolation on two parallel arrays, the first
 * of which contains a set of values to be interpolated and the second of
 * which has the function at those values
 *
 * @param [in] double  x   A value
 * @param [in] double  xarray[]   The array that is interpolated
 * @param [in] double  yarray[]   The array that contains a function of the values in xarray
 * @param [in] int  xdim   The length of the two arrays
 * @param [out] double *  y   The resulting interpolated value
 * @param [in] int  mode   A switch to choose linear(0) or "logarithmic" (1)
 * interpolation
 *
 * @return     The number of the array element that is used for the lower of the two
 * elements that are interpolated on.
 *
 * @details
 * Given a number x, and an array of x's in xarray, and functional
 * values y = f(x) in yarray and the dimension of xarray and yarray,
 * linterp calculates y, the linearly interpolated value of y(x). It
 * also returns the element in the xarray so that one can consider
 * not doing a search to find the right array element in certain
 * circumstances

 *
 * ### Notes ###
 *
 * For mode 0, the value that is returned is
 *
 * (1-f)*y[nelem]+f*[nelem+1)
 *
 * For mode 1, the value returned is
 * exp ((1. - f) * log (y[nelem]) + f * log (y[nelem + 1]))
 *
 *
 **********************************************************/

int
linterp (x, xarray, yarray, xdim, y, mode)
     double x;                  // The value that we wish to index i
     double xarray[], yarray[];
     int xdim;
     double *y;
     int mode;                  //0 = linear, 1 = log
{
  int nelem = 0;
  double frac;


  fraction (x, xarray, xdim, &nelem, &frac, mode);

  if (mode == 0)
    *y = (1. - frac) * yarray[nelem] + frac * yarray[nelem + 1];
  else if (mode == 1)
    *y = exp ((1. - frac) * log (yarray[nelem]) + frac * log (yarray[nelem + 1]));
  else
  {
    Error ("linterp - unknown mode %i\n", mode);
    exit (0);
  }

  return (nelem);

}
