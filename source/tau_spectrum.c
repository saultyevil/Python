/** ************************************************************************* */
/**
 * @file     tau_diag.c
 * @author   Edward Parkinson
 * @date     April 2019
 *
 * @brief    The main resting place for most functions related to providing
 *           optical depth diagnostics.
 *
 * @details
 *
 * To add more PI optical depths to evalulate, one has to add another entry
 * to the PIEDGES array. The rest of the code should be able to cope with any
 * new additions. Remember to put frequency to be slightly blue of the PI edge.
 *
 * ************************************************************************** */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#include "atomic.h"
#include "python.h"

struct sightlines
{
  char name[50];
  double lmn[3];
} *SIGHTLINES;

static int NSIGHTLINES;         // The number of inclination angles

struct edges
{
  char name[50];
  double freq;
} PIEDGES[] = {
  {"HLymanEdge", 3.387485e+15},
  {"HBalmerEdge", 8.293014e+14},
  {"HeIIEdge", 1.394384e+16}
};

static const int NPIEDGES = sizeof PIEDGES / sizeof *PIEDGES;   // The number of optical depths for the simple calculation

#define LAUNCH_DOMAIN 0         // For now, assume we only care about photons starting in domain 0
#define MAXDIFF VCHECK / VLIGHT // For linear velocity requirement for photon transport

/* ************************************************************************** */
/**
 * @brief           Print the various optical depths calculated using this
 *                  routine
 *
 * @param[in]       double   optical_depths    The 2d array containing the optical
 *                                             depth for each observer and tau
 * @param[in]       double   column_densities  The 2d array containing the column
 *                                             densities for each observer
 *
 * @return          void
 *
 * @details
 *
 * Prints the different optical depths for each angle and optical depth.
 * Historically, this used to also print to file. But I found that it was mostly
 * useless to do this, considering the information is already in the diag files
 * by using Log().
 *
 * ************************************************************************** */

static void
print_optical_depths (const double *pi_optical_depths, const double *column_densities)
{
  int i, j;
  int linelen;
  char str[LINELENGTH];
  const int MAX_COL = 120;

  if (rank_global != 0)
    return;

  Log ("Optical depths along the defined line of sights for domain %i:\n\n", LAUNCH_DOMAIN);

  for (i = 0; i < NSIGHTLINES; i++)
  {
    Log ("%s: Hydrogen column density: %3.2e cm^-2\n", SIGHTLINES[i].name, column_densities[i] / MPROT);

    linelen = 0;
    for (j = 0; j < NPIEDGES; j++)
    {
      linelen += sprintf (str, "tau_%-9s: %3.2e  ", PIEDGES[j].name, pi_optical_depths[i * NPIEDGES + j]);
      if (linelen > MAX_COL)
      {
        linelen = 0;
        Log ("\n");
      }

      Log ("%s", str);
    }

    Log ("\n\n");
  }
}

/* ************************************************************************* */
/**
 * @brief           Write the various optical depth spectra to file
 *
 * @param[in]       double tau_spectrum     The various optical depth spectra
 * @param[in]       double wave_min         The smallest wavelength in the
 *                                          spectra
 * @param[in]       double dwave            The size of the wavelength bins
 *
 * @return          void
 *
 * @details
 *
 * Simply write the optical depth spectra to the file named root.tau_spec.diag.
 * This file will be located in the diag folder.
 *
 * ************************************************************************** */

static void
write_optical_depth_spectrum (const double *tau_spectrum, const double freq_min, const double dfreq)
{
  int i, j;
  double current_wavelength, current_freq;
  char filename[1024];
  FILE *fp;

  if (rank_global != 0)
    return;

  sprintf (filename, "diag_%s/%s.tau_spec.diag", files.root, files.root);

  if ((fp = fopen (filename, "w")) == NULL)
  {
    Error ("write_optical_depth_spectrum: unable to open tau spectrum diag file\n");
    Exit (EXIT_FAILURE);
  }

  /*
   * Write out the tau spectrum header
   */


  fprintf (fp, "%-12s %-12s ", "Freq.", "Lambda");
  for (j = 0; j < NSIGHTLINES; j++)
    fprintf (fp, "%-12s ", SIGHTLINES[j].name);
  fprintf (fp, "\n");

  /*
   * Write out the tau spectrum for each inclination angle
   */

  current_freq = freq_min;

  for (i = 0; i < NWAVE; i++)
  {
    current_wavelength = VLIGHT / current_freq / ANGSTROM;
    fprintf (fp, "%-12e %-12e ", current_freq, current_wavelength);

    for (j = 0; j < NSIGHTLINES; j++)
      fprintf (fp, "%-12e ", tau_spectrum[j * NWAVE + i]);
    fprintf (fp, "\n");

    current_freq += dfreq;
  }

  fflush (fp);

  if (fclose (fp))
  {
    Error ("write_optical_depth_spectrum: could not close tau spectrum diag file\n");
    Exit (EXIT_FAILURE);
  }
}

/* ************************************************************************** */
/**
 * @brief  Gather the separate optical depth spectra back to the root process
 *
 * @param[in,out] double *spec  The optical depth spectrum
 * @param[in]     int    nspec  The number of inclination angles
 *
 * @details
 *
 * Uses MPI_Reduce to add all the arrays together. This should be fine as all
 * arrays are initialised to zero. Hence, any empty array entries from the
 * parallelisation design, shouldn't interfere with the fact that we use
 * MPI_SUM.
 *
 * This is most likely poor design. It should be changed in the future to use
 * MPI_Send + MPI_Recv.
 *
 * spec is a constant pointer. The contents of spec can be modified, but the
 * pointer to spec can't be changed.
 *
 * ************************************************************************** */

static void
mpi_gather_spectra (double *const spec, const int nspec)
{
#ifdef MPI_ON
  if (rank_global == 0)
    MPI_Reduce (MPI_IN_PLACE, spec, NWAVE * nspec, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  else
    MPI_Reduce (spec, spec, NWAVE * nspec, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
}

/* ************************************************************************* */
/**
 * @brief           Calculate the total optical depth a photon experiences
 *                  across the cell of distance smax.
 *
 * @param[in]       WindPtr   w           A pointer to the entire wind
 * @param[in]       PhotPtr   pextract    The photon packet to extract
 * @param[in,out]   double    *col_den    The column density the photon has
 *                                        translated through
 * @param[in,out]   double    *tau        The optical depth experienced by the
 *                                        photon
 *
 * @return          int       istat       The current photon status or
 *                                        EXIT_FAILURE on failure.
 *
 * @details
 *
 * This function is concerned with finding the opacity of the photon's current
 * cell, the distance the photon can move in the cell and hence it increments
 * the optical depth tau a photon has experienced as it moves through the wind.
 *
 * The photon status istat is updated and returned for some reason.
 *
 * ************************************************************************** */

static int
calculate_tau_across_cell (WindPtr w, PhotPtr photon, double *current_col_den, double *current_tau)
{
  int istat;
  int ndom, nplasma;
  double kappa_tot, density;
  double smax, diff;
  double freq_inner, freq_outer;
  double mean_freq;
  WindPtr wind_cell;
  PlasmaPtr plasma_cell;
  struct photon p_start, p_stop, p_now;

  if ((photon->grid = where_in_grid (w[photon->grid].ndom, photon->x)) < 0)
  {
    Error ("calculate_tau_across_cell: pextract is not in grid\n");
    return EXIT_FAILURE;
  }

  /*
   * Determine where in the plasma grid the cell is. This is required so a
   * column density can be calculated. Also calculate the maximum distance the
   * photon can traverse across the cell
   */

  wind_cell = &w[photon->grid];
  ndom = wind_cell->ndom;
  nplasma = w[photon->grid].nplasma;
  plasma_cell = &plasmamain[nplasma];
  density = plasma_cell->rho;
  smax = smax_in_cell (photon);
  if (smax < 0)
  {
    Error ("calculate_tau_across_cell: abnormal value of smax for photon\n");
    return EXIT_FAILURE;
  }

  // Transform photon at starting location to local frame
  observer_to_local_frame (photon, &p_start);

  // Move photon smax and transform to local frame
  stuff_phot (photon, &p_stop);
  move_phot (&p_stop, smax);
  observer_to_local_frame (&p_stop, &p_stop);

  /* At this point p_start and p_stop are in the local frame
   * at the and p_stop is at the maximum distance it can
   * travel. We want to check that the frequency shift is
   * not too great along the path that a linear approximation
   * to the change in frequency is not reasonable
   */

  while (smax > DFUDGE)
  {
    stuff_phot (photon, &p_now);
    move_phot (&p_now, smax * 0.5);
    observer_to_local_frame (&p_now, &p_now);
    diff = fabs (p_now.freq - 0.5 * (p_start.freq + p_stop.freq)) / p_start.freq;
    if (diff < MAXDIFF)
    {
      break;
    }
    stuff_phot (&p_now, &p_stop);
    smax *= 0.5;
  }

  freq_inner = p_start.freq;
  freq_outer = p_stop.freq;
  mean_freq = 0.5 * (freq_inner + freq_outer);

  /*
   * Now we can finally calculate the opacity due to all the continuum
   * processes. In macro-atom mode, we need to calculate the continuum opacity
   * using kappa_bf and kappa_ff using the macro treatment. For simple mode, we
   * can simply use radiation which **SHOULD** return the continuum opacity
   * as well, plus something from induced Compton heating. In either cases,
   * we still then need to add the optical depth from electron scattering at
   * the end.
   */

  kappa_tot = 0;

  if (geo.rt_mode == RT_MODE_2LEVEL)
  {
    kappa_tot += radiation (photon, smax);
  }
  else if (geo.rt_mode == RT_MODE_MACRO)
  {
    if (wind_cell->vol > 0)
    {
      kappa_tot += kappa_bf (plasma_cell, freq_inner, 0);
      kappa_tot += kappa_ff (plasma_cell, freq_inner);
    }
  }

  kappa_tot += klein_nishina (mean_freq) * plasma_cell->ne * zdom[ndom].fill;

  /*
   * Increment the optical depth and column density variables and move the
   * photon to the edge of the cell
   */

  *current_col_den += smax * density;
  *current_tau += smax * kappa_tot;
  move_phot (photon, smax);
  istat = photon->istat;

  return istat;
}

/* ************************************************************************* */
/**
 * @brief           Extract the optical depth the photon packet porig must
 *                  travel through to reach the observer.
 *
 * @param[in]       WindPtr   w          A pointer to the entire wind
 * @param[in]       PhotPtr   porig      The photon packet to extract
 * @param[out]      double    *col_den   The column depth of the extracted
 *                                       photon angle
 * @param[out]      double    *tau       The optical depth from photon origin to
 *                                       the observer
 *
 * @return          int                  EXIT_SUCCESS or EXIT_FAILURE
 *
 * @details
 *
 * The direction of the observer is set in porig.lmn and should be set prior
 * to passing a photon to this function. If any values of lmn are negative, the
 * photon will translate in the negative direction. However, have no fear as this
 * is normal and is fine due to the assumed symmetry of models in Python.
 *
 * ************************************************************************** */

static int
extract_tau (WindPtr w, PhotPtr photon, double *current_col_den, double *current_tau)
{
  int err, ndom;
  enum istat_enum istat;
  int ninspace;
  const int max_space = 10;
  double norm[3];
  struct photon pextract;

  istat = P_INWIND;             // assume photon is in wind for initialisation reasons
  stuff_phot (photon, &pextract);

  ninspace = 0;

  while (istat == P_INWIND)
  {
    // Move the photon in the grid, but not in the wind
    if (where_in_wind (pextract.x, &ndom) < 0)
    {
      translate_in_space (&pextract);
      if (++ninspace > max_space)
      {
        Error ("extract_tau: tau_extract photon transport ended due to too many translate_in_space\n");
        return EXIT_FAILURE;
      }
    }
    // Move the photon in the wind
    else if ((pextract.grid = where_in_grid (ndom, pextract.x)) >= 0)
    {
      err = calculate_tau_across_cell (w, &pextract, current_col_den, current_tau);
      if (err)
      {
        pextract.istat = -1;
        return EXIT_FAILURE;
      }
    }
    else
    {
      Error ("extract_tau: photon in unknown location, grid stat %i\n", pextract.grid);
      pextract.istat = P_ERROR;
      return EXIT_FAILURE;
    }

    istat = walls (&pextract, photon, norm);
  }

  if (istat == P_HIT_STAR || istat == P_HIT_DISK)
  {
    Error ("extract_tau: photon hit central source or disk incorrectly istat = %i\n", istat);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

/* ************************************************************************* */
/**
 * @brief           Reposition a photon depending on various criteria
 *
 * @param[in, out]  PhotPtr pout    The photon to reposition
 *
 * @return          void
 *
 * @details
 *
 * This function assumes that we only care about DEFAULT_DOMAIN, whatever
 * that has been set as.
 *
 * For photons which are pointing in the polar direction, these are launched
 * from most dense region of the lower wind.
 *
 * ************************************************************************** */

static void
reposition_tau_photon (PhotPtr pout)
{
  int i;
  int wind_index, plasma_index;
  double xloc;
  static double upward_xloc;
  double current_density, max_density;
  static bool xloc_init = false;

  /*
   * Hacky code for when photon is pointing upwards - this photon will not
   * be launched from the origin, but will be launched from the most dense
   * cell at the bottom of the disk.
   */

  if (xloc_init == false)
  {
    max_density = 0;
    for (i = 0; i < zdom[LAUNCH_DOMAIN].ndim - 1; i++)
    {
      wind_index = i * zdom[LAUNCH_DOMAIN].mdim + zdom[LAUNCH_DOMAIN].mdim;
      xloc = wmain[wind_index].x[0];
      plasma_index = wmain[wind_index].nplasma;
      current_density = plasmamain[plasma_index].rho;
      if (current_density > max_density)
      {
        upward_xloc = xloc + wmain[wind_index].dfudge;
        max_density = current_density;
      }
    }

    xloc_init = true;
  }

  if (pout->lmn[0] == 0 && pout->lmn[1] == 0 && pout->lmn[2] == 1)
  {
    pout->x[0] = upward_xloc;
    pout->x[2] = DFUDGE;
  }
}

/* ************************************************************************* */
/**
 * @brief           Generate a photon packet with a given frequency nu for use
 *                  with the optical depth diagnostic routines.
 *
 * @param[in,out]   PhotPtr  pout     The photon packet to initialise
 *
 * @param[in]       double   nu       The frequency of the photon packet
 *
 * @return          int               EXIT_SUCCESS or EXIT_FAILURE
 *
 * @details
 *
 * This routine assumes that the user wishes for the photon to be generated from
 * the surface of the central source, taking into account the current direction
 * of the sight-line being extracted. It's currently not possible, without a
 * tedious re-write, to place the photon at the very origin of the grid when
 * there is a central source because of how we check boundaries.
 *
 * Note that photons are initialised with a weight of f_tot as photons are
 * required to have weight, but since functions do not care about the weight of
 * the photon, it is set to something large to make sure it does not get
 * destroyed by accident somehow.
 *
 * ************************************************************************** */

static int
create_tau_photon (PhotPtr pout, double freq, double *lmn)
{
  if (freq < 0)
  {
    Error ("create_tau_photon: photon can't be created with negative frequency\n");
    return EXIT_FAILURE;
  }

  stuff_v (lmn, pout->lmn);
  pout->freq = pout->freq_orig = freq;
  pout->origin = pout->origin_orig = PTYPE_DISK;
  pout->istat = P_INWIND;
  pout->w = pout->w_orig = geo.f_tot;
  pout->tau = 0.0;
  pout->x[0] = pout->x[1] = pout->x[2] = 0;
  pout->frame = F_OBSERVER;
  move_phot (pout, geo.rstar + DFUDGE);
  reposition_tau_photon (pout);

  return EXIT_SUCCESS;
}

/* ************************************************************************* */
/**
 * @brief           Initialise the viewing angles for the optical depth
 *                  diagnostic routines.
 *
 * @param           void
 *
 * @return          void
 *
 * @details
 *
 * This purpose of this function is to initialise the angles of which to extract
 * optical depth information from. If xxspec is initialised, then the optical
 * depth information will be extracted from the observer angles provided in the
 * parameter file for spectrum generation. Otherwise, the routines will instead
 * use a set of pre-defined angles.
 *
 * ************************************************************************** */

static void
init_sightlines (void)
{
  int iangle, mem_req;
  const double default_phase = 0.5;
  const double default_angles[] = { 0.0, 10.0, 30.0, 45.0, 60.0, 75.0, 85.0, 90.0 };
  const int n_default_angles = sizeof default_angles / sizeof default_angles[0];

  /*
   * Use the angles specified for by the user for spectrum generation, this
   * requires for xxspec to be initialised.
   */

  if (geo.nangles > 0 && xxspec != NULL)
  {
    NSIGHTLINES = geo.nangles;
    SIGHTLINES = calloc (geo.nangles, sizeof *SIGHTLINES);
    if (SIGHTLINES == NULL)
    {
      mem_req = geo.nangles * (int) sizeof *SIGHTLINES;
      Error ("init_sightlines: cannot allocate %d bytes for observers array\n", mem_req);
      Exit (EXIT_FAILURE);
    }
    else                        // Use else to avoid compiler warning
    {
      for (iangle = MSPEC; iangle < MSPEC + geo.nangles; iangle++)
      {
        strcpy (SIGHTLINES[iangle - MSPEC].name, xxspec[iangle].name);
        stuff_v (xxspec[iangle].lmn, SIGHTLINES[iangle - MSPEC].lmn);
      }
    }
  }

  /*
   * If no spectrum cycles are being run, or if the model is being restarted
   * without initialising xxspec, then a default set of angles is used instead.
   */

  else
  {
    Log_silent ("tau_spectrum: as there are no spectrum cycles or observers defined, a set of default angles will be used instead\n");
    NSIGHTLINES = n_default_angles;
    SIGHTLINES = calloc (n_default_angles, sizeof *SIGHTLINES);
    if (SIGHTLINES == NULL)
    {
      mem_req = n_default_angles * (int) sizeof *SIGHTLINES;
      Error ("init_sightlines: cannot allocate %d bytes for observers array\n", mem_req);
      Exit (EXIT_FAILURE);
    }
    else
    {
      for (iangle = 0; iangle < n_default_angles; iangle++)
      {
        sprintf (SIGHTLINES[iangle].name, "A%02.0fP%04.2f", default_angles[iangle], default_phase);
        SIGHTLINES[iangle].lmn[0] = sin (default_angles[iangle] / RADIAN) * cos (-default_phase * 360.0 / RADIAN);
        SIGHTLINES[iangle].lmn[1] = sin (default_angles[iangle] / RADIAN) * sin (-default_phase * 360.0 / RADIAN);
        SIGHTLINES[iangle].lmn[2] = cos (default_angles[iangle] / RADIAN);
      }
    }
  }
}

/* ************************************************************************* */
/**
 * @brief           Create spectra of tau vs lambda for each observer angle
 *
 * @param[in]       WindPtr    w         A pointer to the wind
 *
 * @return          void
 *
 * @details
 *
 * This is the main function which will generate the optical depth spectra for
 * each observer angle in xxspec. The algorithm is similar to extract and the
 * tau_diag algorithm which this function is called in.
 *
 * A photon is generated at the central source of the model and is extracted
 * from this location towards the observer where it escapes, where extract_tau
 * returns the integrated optical depth along its path to escape. This is done
 * for a range of photon frequencies to find how optical depth changes with
 * frequency.
 *
 * This processes can take some time compared to tau_evalulate_photo_edges. But,
 * since NBINS photons are being generated for each spectrum and the fact that
 * these photons do not interact, the spectra does not actually take that long
 * to complete. If NBINS is large enough, it can take a bit longer than one
 * would like, hence the reason why this routine is parallelised.
 *
 * NOTES
 * MPI parallelisation has been disabled for now. It doesn't work when the
 * number of processes > 24.
 *
 * ************************************************************************** */

static void
create_optical_depth_spectrum (WindPtr w)
{
  int i, j;
  int err;
  int nbins, mpi_lower, mpi_upper;
  double tau, column;
  double *tau_spectrum;
  double current_freq;
  double freq_min, freq_max, dfreq;
  struct photon photon;

  Log ("Creating optical depth spectra:\n");

  tau_spectrum = calloc (NSIGHTLINES * NWAVE, sizeof *tau_spectrum);
  if (tau_spectrum == NULL)
  {
    Error ("create_optical_depth_spectrum: cannot allocate %d bytes for tau_spectrum\n", NSIGHTLINES * NWAVE * sizeof *tau_spectrum);
    return;
  }

  /*
   * Define the limits of the spectra in frequency space. If xxpsec is NULL,
   * then the frequency range will be over a default 100 - 10,000 Angstrom
   * band.
   */

  if ((geo.nangles == 0 && xxspec == NULL) || (geo.swavemax == 0 && geo.swavemin == 0))
  {
    Error ("create_optical_depth_spectrum: xxspec is uninitialized, defaulting spectral wavelength range to 100 - 10,000 Angstrom\n");
    freq_min = VLIGHT / (10000 * ANGSTROM);
    freq_max = VLIGHT / (100 * ANGSTROM);
  }
  else
  {
    freq_min = VLIGHT / (geo.swavemax * ANGSTROM);
    freq_max = VLIGHT / (geo.swavemin * ANGSTROM);
    if (sane_check (freq_min))
    {
      freq_min = VLIGHT / (10000 * ANGSTROM);
      Error ("create_optical_depth_spectrum: freq_min has an invalid value setting to %e\n", freq_min);
    }
    if (sane_check (freq_max))
    {
      freq_max = VLIGHT / (100 * ANGSTROM);
      Error ("create_optical_depth_spectrum: freq_min has an invalid value setting to %e\n", freq_max);
    }
  }

  dfreq = (freq_max - freq_min) / NWAVE;
  kbf_need (freq_min, freq_max);

  /*
   * Initialise the variables to split the calculation across multiple
   * processors
   */

  // MPI has been temporarily disabled
  // nbins = ceil ((double) NWAVE / np_mpi_global);
  // mpi_lower = nbins * rank_global;
  // mpi_upper = nbins * (rank_global + 1);
  // if (mpi_upper > NWAVE)
  //   mpi_upper = NWAVE;

  mpi_lower = 0;
  mpi_upper = NWAVE;

  /*
   * Now create the optical depth spectra for each sightline
   */

  for (i = 0; i < NSIGHTLINES; i++)
  {
    Log ("  - Creating spectrum: %s\n", SIGHTLINES[i].name);
    current_freq = mpi_lower * dfreq;

    for (j = mpi_lower; j < mpi_upper; j++)
    {
      tau = 0.0;
      column = 0.0;
      err = create_tau_photon (&photon, current_freq, SIGHTLINES[i].lmn);
      if (err == EXIT_FAILURE)
      {
        Log ("create_optical_depth_spectrum: skipping photon of frequency %e when creating spectrum\n", current_freq);
        continue;
      }

      err = extract_tau (w, &photon, &column, &tau);
      if (err == EXIT_FAILURE)
        continue;

      tau_spectrum[i * NWAVE + j] = tau;
      current_freq += dfreq;
    }
  }

  // mpi_gather_spectra (tau_spectrum, NANGLES);
  write_optical_depth_spectrum (tau_spectrum, freq_min, dfreq);
  free (tau_spectrum);
}

/* ************************************************************************* */
/**
 * @brief           The main controlling function for the optical depth
 *                  diagnostics.
 *
 * @param[in]       WindPtr    w       A pointer to the wind
 *
 * @return          void
 *
 * @details
 *
 * This is the main function which will control the procedure for calculating
 * various diagnostic numbers for the optical depth's experienced in the current
 * model. Namely, this function aims to show the total integrated optical depth
 * to each observer angle using the following optical depths:
 *
 *  - Lymann edge
 *  - Balmer edge
 *  - Helium II edge
 *
 * Once these integrated optical depths have been calculated for each angle, a
 * spectrum of optical depth vs wavelength is created for each angle.
 *
 * The aim of these diagnostic numbers it to provide some sort of metric on the
 * optical thickness of the current model.
 *
 * TODO: method for calling this for individual cells instead of the entire wind
 *
 * ************************************************************************** */

static void
calculate_pi_optical_depth (WindPtr w)
{
  int err;
  int j, i;
  double current_freq, current_tau, current_col_density;
  double *pi_optical_depths, *column_densities;
  struct photon photon;

  pi_optical_depths = calloc (NSIGHTLINES * NPIEDGES, sizeof *pi_optical_depths);
  if (pi_optical_depths == NULL)
  {
    Error ("calculate_photoion_edge_optical_depth: cannot allocate %d bytes for optical_depths\n",
           NSIGHTLINES * NPIEDGES * sizeof *pi_optical_depths);
    return;
  }

  /*
   * Allocate storage for the column densities for each edge and sight line
   */

  column_densities = calloc (NSIGHTLINES, sizeof *pi_optical_depths);
  if (column_densities == NULL)
  {
    Error ("calculate_photoion_edge_optical_depth: cannot allocate %d bytes for column_densities\n",
           NSIGHTLINES * sizeof *pi_optical_depths);
    free (pi_optical_depths);
    return;
  }

  /*
   * Now extract the optical depths and mass column densities
   */

  for (i = 0; i < NSIGHTLINES; i++)
  {
    for (j = 0; j < NPIEDGES; j++)
    {
      current_tau = 0.0;
      current_col_density = 0.0;
      current_freq = PIEDGES[j].freq;

      err = create_tau_photon (&photon, current_freq, SIGHTLINES[i].lmn);
      if (err == EXIT_FAILURE)
      {
        Error ("calculate_photoion_edge_optical_depth: skipping photon of frequency %e\n", current_freq);
        continue;
      }

      err = extract_tau (w, &photon, &current_col_density, &current_tau);
      if (err == EXIT_FAILURE)  // Do not throw another "extra" error
        continue;

      pi_optical_depths[i * NPIEDGES + j] = current_tau;
      column_densities[i] = current_col_density;
    }
  }

  print_optical_depths (pi_optical_depths, column_densities);
  free (pi_optical_depths);
  free (column_densities);
}

/* ************************************************************************* */
/**
 * @brief  Main control function for create optical depth diagnostics.
 *
 * @param[in]  WindPtr w  A pointer to the wind
 *
 * @details
 *
 * This function is the main steering function for generating the optical depth
 * diagnostics. In MPI mode, all MPI processes should be doing an identical
 * amount of work, however only the master process (0) will write any diagnostics
 * out to file. To avoid any complications due to the processes being out of
 * sync with one another, a MPI_Barrier is used at the end to allow the processes
 * to resynchronise.
 *
 * ************************************************************************** */

void
optical_depth_diagnostics (WindPtr w)
{
  /*
   * We have to set this flag, so radiation returns without updating the
   * radiation field and etc estimators. Don't worry, we switch it back to
   * what it was later...
   */

  geo.ioniz_or_extract = 0;
  xsignal (files.root, "%-20s Optical depth diagnostics beginning\n", "TAU");

  init_sightlines ();
  calculate_pi_optical_depth (w);
  create_optical_depth_spectrum (w);

#ifdef MPI_ON
  MPI_Barrier (MPI_COMM_WORLD);
#endif

  free (SIGHTLINES);
  geo.ioniz_or_extract = 1;
  xsignal (files.root, "%-20s Optical depth diagnostics completed\n", "TAU");
  Log (" Completed optical depth diagnostics. The elapsed TIME was %f\n", timer ());
}
