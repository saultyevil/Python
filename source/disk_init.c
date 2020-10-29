
/***********************************************************/
/** @file  disk_init.c 
 * @author ksl
 * @date   April, 2020
 *
 * @brief  Primary routines for initializing the disk
 * as on a luminoisity weighted basis
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


#define STEPS 100000



/**********************************************************/
/**
 * @brief      calculates the total luminosity and the luminosity between freqqmin and freqmax
 * 	of the disk.  More importantly  divides the disk into annulae such that each
 * 	annulus contributes and equal amount to the lumionosity of the disk (within the frequency
 * 	limits).  Thus  initializes the structure "disk".
 *
 * @param [in] double  rmin   The minimum radius of the disk
 * @param [in] double  rmax   The maximum radius of the disk
 * @param [in] double  m   mass of central object
 * @param [in] double  mdot   mass accretion rate
 * @param [in] double  freqmin   The minimum frequency
 * @param [in] double  freqmax   The maximum frequency
 * @param [in] int  ioniz_or_final   A flag indicating whether this is an ionization or
 * a detailed spectral cycle (used to determine the spectral type to use)
 * @param [out] double *  ftot   The band limited luminosity in the freqency interval
 * @return     the total luminosity of the disk
 *
 * @details
 * This routine assumes the temperature distribution for the disk is
 * that of a simple Shakura-Sunyaev disk, and uses this to determine
 * the band limited luminosity of the disk.  Additionally, it divides
 * the disk in the rings of the same band-limited luminosity, so that
 * equal numbers of photons can be generated from each ring.  (The
 * reason the disk has to be initilaized mulitple times is because
 * the rings are different for different freqency intervals.)
 *
 * ### Notes ###
 * The information needed to generate photons from the disk is stored
 * in the disk structure.
 * The positional parameters x and v are at the edge of the ring,
 * but many of the other parameters (like temperature) are at the mid point.
 *
 *
 **********************************************************/

double
disk_init (rmin, rmax, m, mdot, freqmin, freqmax, ioniz_or_final, ftot)
     double rmin, rmax, m, mdot, freqmin, freqmax, *ftot;
     int ioniz_or_final;
{
  double t, tref;
  double log_g, gref;
  double v, dr, r;
  double logdr, logrmin, logrmax, logr;
  double f, ltot;
  double q1;
  int nrings, i, icheck;
  int spectype;
  double emit, emittance_bb (), emittance_continuum ();
  double factor;

  /* Calculate the reference temperature and luminosity of the disk */
  tref = tdisk (m, mdot, rmin);
  gref = gdisk (m, mdot, rmin);

  /* Now compute the apparent luminosity of the disk.  This is not actually used
     to determine how annulae are set up.  It is just used to populate geo.ltot.
     It can change if photons hitting the disk are allowed to raise the temperature
   */

  logrmax = log (rmax);
  logrmin = log (rmin);
  logdr = (logrmax - logrmin) / STEPS;

  for (nrings = 0; nrings < NRINGS; nrings++)
  {
    disk.nphot[nrings] = 0;
    disk.nphot[nrings] = 0;
    disk.r[nrings] = 0;
    disk.t[nrings] = 0;
    disk.nhit[nrings] = 0;
    disk.heat[nrings] = 0;
    disk.ave_freq[nrings] = 0;
    disk.w[nrings] = 0;
    disk.t_hit[nrings] = 0;
  }




  ltot = 0;

  for (logr = logrmin; logr < logrmax; logr += logdr)
  {
    r = exp (logr);
    dr = exp (logr + logdr) - r;
    t = teff (tref, (r + 0.5 * dr) / rmin);
    ltot += t * t * t * t * (2. * r + dr) * dr;
  }
  geo.lum_disk_init = ltot *= 2. * STEFAN_BOLTZMANN * PI;


  /* Now establish the type of spectrum to create */

  if (ioniz_or_final == 1)
    spectype = geo.disk_spectype;       /* type for final spectrum */
  else
    spectype = geo.disk_ion_spectype;   /*type for ionization calculation */

/* Next compute the band limited luminosity ftot */

/* The area of an annulus is  PI*((r+dr)**2-r**2) = PI * (2. * r +dr) * dr.
   The extra factor of two arises because the disk radiates from both of its sides.
   */

  q1 = 2. * PI;

  (*ftot) = 0;
  icheck = 0;


  for (logr = logrmin; logr < logrmax; logr += logdr)
  {
    r = exp (logr);
    dr = exp (logr + logdr) - r;
    t = teff (tref, (r + 0.5 * dr) / rmin);
    log_g = log10 (geff (gref, (r + 0.5 * dr) / rmin));
    v = sqrt (GRAV * geo.mstar / r);
    v /= VLIGHT;
    if (rel_mode == REL_MODE_FULL)
    {
      factor = sqrt (1. - v * v);
    }
    else
    {
      factor = 1.0;
    }

    if (spectype > -1)
    {
      emit = emittance_continuum (spectype, freqmin, freqmax, t, log_g);
    }
    else
    {
      emit = emittance_bb (freqmin, freqmax, t);

    }
    (*ftot) += emit * (2. * r + dr) * dr * factor;
  }

  (*ftot) *= q1;



  /* If *ftot is 0 in this energy range then all the photons come elsewhere, e. g. the star or BL  */

  if ((*ftot) < EPSILON)
  {
    Log_silent ("disk_init: Warning! Disk does not radiate enough to matter in this wavelength range\n");
    return (ltot);
  }

  /* Now find the boundaries of the each annulus, which depends on the band limited flux.
     Note that disk.v is calculated at the boundaries, because vdisk() interporlates on
     the actual radius. */

  disk.r[0] = rmin;
  disk.v[0] = sqrt (GRAV * geo.mstar / rmin);
  nrings = 1;
  f = 0;

  i = 0;
  for (logr = logrmin; logr < logrmax; logr += logdr)
  {
    r = exp (logr);
    dr = exp (logr + logdr) - r;
    t = teff (tref, (r + 0.5 * dr) / rmin);
    log_g = log10 (geff (gref, (r + 0.5 * dr) / rmin));
    v = sqrt (GRAV * geo.mstar / r);
    v /= VLIGHT;
    if (rel_mode == REL_MODE_FULL)
    {
      factor = sqrt (1. - v * v);
    }
    else
    {
      factor = 1.0;
    }

    if (spectype > -1)
    {
      emit = emittance_continuum (spectype, freqmin, freqmax, t, log_g);
    }
    else
    {
      emit = emittance_bb (freqmin, freqmax, t);
    }

    f += q1 * emit * (2. * r + dr) * dr * factor;
    i++;
    /* EPSILON to assure that roundoffs don't affect result of if statement */
    if (f / (*ftot) * (NRINGS - 1) >= nrings)
    {
      if (r <= disk.r[nrings - 1])      //If the radius we have reached is smaller than or equal to the last assigned radius - we make a tiny annulus
      {
        r = disk.r[nrings - 1] * (1. + 1.e-10);
      }
      disk.r[nrings] = r;
      disk.v[nrings] = sqrt (GRAV * geo.mstar / r);
      nrings++;
      if (nrings >= NRINGS)
      {
        break;
      }
    }
  }
  if (nrings < NRINGS - 1)
  {
    Error ("error: disk_init: Integration on setting r boundaries got %d nrings instead of %d\n", nrings, NRINGS - 1);
    Exit (0);
  }


  disk.r[NRINGS - 1] = exp (logrmax);
  disk.v[NRINGS - 1] = sqrt (GRAV * geo.mstar / disk.r[NRINGS - 1]);


  /* Now calculate the temperature and gravity of the annulae */

  for (nrings = 0; nrings < NRINGS - 1; nrings++)
  {
    r = 0.5 * (disk.r[nrings + 1] + disk.r[nrings]);
    disk.t[nrings] = teff (tref, r / rmin);
    disk.g[nrings] = geff (gref, r / rmin);
  }

  /* Wrap up by zeroing other parameters */
  for (nrings = 0; nrings < NRINGS; nrings++)
  {
    disk.nphot[nrings] = 0;
    disk.nhit[nrings] = 0;
    disk.heat[nrings] = 0;
    disk.ave_freq[nrings] = 0;
    disk.w[nrings] = 0;
    disk.t_hit[nrings] = 0;
  }
  geo.lum_disk = ltot;
  return (ltot);
}




/**********************************************************/
/** 
 * @brief      Initialize a structure (qdisk) for recording information about photons/energy impinging
 * 	on the disk, which is stored in a disk structure called qdisk.
 *
 * @return     Always return zero
 *
 *
 * ###Notes###
 *
 * The information stored in qdisk can be used to modify the effective temperature
 * of the disk
 *
 * The reason the qdisk is needed as well as disk structure is that whenever the
 * 	wavelength bands are changed the radii in the disk structure are recalibratied.
 * 	We want qdisk to have fixed boundaries when this happens.
 *
 *
 * The annular are defined differently than in disk_init.  Here they are simply
 * logarithmically spaced; there they are spaced so that equal amounts of emission
 * are emitted from each annulus.
 *
 **********************************************************/

int
qdisk_init (rmin, rmax, m, mdot)
     double rmin, rmax, m, mdot;
{
  int nrings;
  double log_rmin, log_rmax, dlog_r, log_r;
  double r;
  double tref, gref;

  log_rmin = log10 (disk.r[0]);
  log_rmax = log10 (disk.r[NRINGS - 1]);
  dlog_r = (log_rmax - log_rmin) / (NRINGS - 1);


  for (nrings = 0; nrings < NRINGS; nrings++)
  {
    log_r = log_rmin + dlog_r * nrings;
    qdisk.r[nrings] = pow (10, log_r);
  }

  /* Calculate the reference temperature and luminosity of the disk */
  tref = tdisk (m, mdot, rmin);
  gref = gdisk (m, mdot, rmin);

  for (nrings = 0; nrings < NRINGS; nrings++)
  {
    if (nrings < NRINGS - 1)
    {
      r = 0.5 * (qdisk.r[nrings + 1] + qdisk.r[nrings]);
    }
    else
    {
      r = qdisk.r[nrings];
    }
    qdisk.t[nrings] = teff (tref, r / rmin);
    qdisk.g[nrings] = geff (gref, r / rmin);
    qdisk.v[nrings] = sqrt (GRAV * geo.mstar / r);
    qdisk.heat[nrings] = 0.0;
    qdisk.nphot[nrings] = 0;
    qdisk.nhit[nrings] = 0;
    qdisk.w[nrings] = 0;
    qdisk.ave_freq[nrings] = 0;
    qdisk.t_hit[nrings] = 0;
  }

  return (0);
}


/**********************************************************/
/** 
 * @brief      Save the information in the qdisk structure to a file
 *
 * @param [in] char *  diskfile   Name of the file whihc is writteen
 * @param [in] double  ztot   The total luminosity of the disk as
 *                      calculate over multiple cycles
 * @return     Always returns 0
 *
 * The routine reformats the data about disk heating which has
 * been accumulated and writes it to a file
 *
 * ###Notes###
 *
 * The data concerning heating by the disk is built up during
 * a the ionization cycles
 *
 * The file that is produced should be readable as an astropy
 * table
 *
 **********************************************************/

int
qdisk_save (diskfile, ztot)
     char *diskfile;
     double ztot;
{
  FILE *qptr;
  int n;
  double area, theat, ttot;
  qptr = fopen (diskfile, "w");
  fprintf (qptr, "r          zdisk      t_disk    heat       nhit nhit/nemit  t_heat    t_irrad  W_irrad  t_tot\n");

  for (n = 0; n < NRINGS; n++)
  {
    area = (2. * PI * (qdisk.r[n + 1] * qdisk.r[n + 1] - qdisk.r[n] * qdisk.r[n]));
    theat = qdisk.heat[n] / area;
    theat = pow (theat / STEFAN_BOLTZMANN, 0.25);       // theat is temperature if no internal energy production
    if (qdisk.nhit[n] > 0)
    {

      qdisk.ave_freq[n] /= qdisk.heat[n];
      qdisk.t_hit[n] = PLANCK * qdisk.ave_freq[n] / (BOLTZMANN * 3.832);        // Basic conversion from freq to T
      qdisk.w[n] = qdisk.heat[n] / (4. * PI * STEFAN_BOLTZMANN * area * qdisk.t_hit[n] * qdisk.t_hit[n] * qdisk.t_hit[n] * qdisk.t_hit[n]);
    }

    ttot = pow (qdisk.t[n], 4) + pow (theat, 4);
    ttot = pow (ttot, 0.25);
    fprintf (qptr,
             "%9.4e %9.4e %8.3e %8.3e %5d %8.3e %8.3e %8.3e %8.3e %8.3e\n",
             qdisk.r[n], zdisk (qdisk.r[n]), qdisk.t[n],
             qdisk.heat[n], qdisk.nhit[n], qdisk.heat[n] * NRINGS / ztot, theat, qdisk.t_hit[n], qdisk.w[n], ttot);
  }

  fclose (qptr);
  return (0);
}





/**********************************************************/
/** 
 * @brief      Read the temperature profile from a file
 *
 * @param [in, out] char *  tprofile   Name of the input file
 * @return     Always returns 0
 *
 * The input format for the file to be read in is quite spefic
 * The first line should contina the number of points n
 * the profile
 *
 * Subsequent lines should contian a radius and a temperarue
 *
 * The radius values shoule be in units of 1e11 cm
 * The temperature should be in units of 1000K
 *
 * ###Notes###
 *
 * @bug This routine which was written for the YSO study
 *  needs to be made less YSO centric. It should also be
 *  retested.
 *
 **********************************************************/

int
read_non_standard_disk_profile (tprofile)
     char *tprofile;
{

  FILE *fopen (), *fptr;
  int n;
  float dumflt1, dumflt2;

  char *line;
  size_t buffsize = LINELENGTH;

  if ((fptr = fopen (tprofile, "r")) == NULL)
  {
    Error ("Could not open filename %s\n", tprofile);
    Exit (0);
  }

  line = (char *) malloc (buffsize * sizeof (char));
  blmod.n_blpts = 0;


  while (getline (&line, &buffsize, fptr) > 0)
  {
    n = sscanf (line, "%g %g", &dumflt1, &dumflt2);
    if (n == 2)
    {
      blmod.r[blmod.n_blpts] = dumflt1;
      blmod.t[blmod.n_blpts] = dumflt2;
      blmod.n_blpts += 1;
    }
    else
    {
      Error ("read_non_standard_disk_file: could not convert a line in %s, OK if comment\n", tprofile);
    }

    if (blmod.n_blpts == NBLMODEL)
    {
      Error ("read_non_standard_disk_file: More than %d points in %s; increase NBLMODEL\n", NBLMODEL, tprofile);
      Exit (1);

    }
  }


  if (geo.diskrad > blmod.r[blmod.n_blpts - 1])
  {
    Error ("read_non_standard_disk_profile: The disk radius (%.2e) exceeds rmax (%.2e) in the temperature profile\n", geo.diskrad,
           blmod.r[blmod.n_blpts - 1]);
    Log ("read_non_standard_disk_profile: Portions of the disk outside are treated as part of a steady state disk\n");
  }

  fclose (fptr);

  return (0);
}
