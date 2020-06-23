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
 * ### HOW TO ADD MORE OPTICAL DEPTHS ###
 *
 * In theory it should be relatively straight forwards to add more optical
 * depths to the photoionization edge diagnostics if desired. One should only
 * need to update the OPACITY_EDGES[] array with the relevant edge frequency
 * and name of that edge. Note that the frequency should actually be slightly
 * blueward of the actual edge frequency - this is mostly to take into account
 * any blueshifting.
 *
 * If none of this makes sense, just try to follow how everything else is done!
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
};

struct edges
{
  double nu;
  char name[50];
};

static struct sightlines *SIGHTLINES;
static const struct edges EDGES[] = {
  {3.387485e+15, "HLymanEdge"},
  {8.293014e+14, "HBalmerEdge"},
  {1.394384e+16, "HeIIEdge"}
};

static int NANGLES;             // The number of inclination angles
static const int NTAU = sizeof EDGES / sizeof *EDGES;   // The number of optical depths for the simple calculation
static int NBINS_SPEC = NWAVE;  // The number of bins for the spectra - modified for large frequency range

#define LAUNCH_DOMAIN 0         // For now, assume we only care about photons starting in domain 0
#define MAXDIFF VCHECK/VLIGHT   // The same as our old velocity requirement

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
print_optical_depths (const double *optical_depths, const double *column_densities)
{
  int itau, ispec;
  int line_len;
  char tmp_str[LINELENGTH];
  char observer_name[LINELENGTH];

  const int MAX_COL = 120;

  if (rank_global != 0)
    return;

  Log ("Optical depths along the defined line of sights for domain %i:\n\n", LAUNCH_DOMAIN);

  for (ispec = 0; ispec < NANGLES; ispec++)
  {
    strcpy (observer_name, SIGHTLINES[ispec].name);
    Log ("%s: Hydrogen column density: %3.2e cm^-2\n", observer_name, column_densities[ispec] / MPROT);

    line_len = 0;
    for (itau = 0; itau < NTAU; itau++)
    {
      line_len += sprintf (tmp_str, "tau_%-9s: %3.2e  ", EDGES[itau].name, optical_depths[ispec * NTAU + itau]);
      if (line_len > MAX_COL)
      {
        line_len = 0;
        Log ("\n");
      }

      Log ("%s", tmp_str);
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
  int ispec, ifreq;
  double wavelength, frequency;
  char filename[1024];
  FILE *filep;

  if (rank_global != 0)
    return;

  sprintf (filename, "diag_%s/%s.tau_spec.diag", files.root, files.root);

  if ((filep = fopen (filename, "w")) == NULL)
  {
    Error ("%s : %i : unable to open tau spectrum diag file\n", __FILE__, __LINE__);
    Exit (1);
  }

  /*
   * Write out the tau spectrum header
   */


  fprintf (filep, "%-12s %-12s ", "Freq.", "Lambda");
  for (ispec = 0; ispec < NANGLES; ispec++)
    fprintf (filep, "%-12s ", SIGHTLINES[ispec].name);
  fprintf (filep, "\n");

  /*
   * Write out the tau spectrum for each inclination angle
   */

  frequency = freq_min;

  for (ifreq = 0; ifreq < NBINS_SPEC; ifreq++)
  {
    wavelength = VLIGHT / frequency / ANGSTROM; // TODO check correct conversion
    fprintf (filep, "%-12e %-12e ", frequency, wavelength);
    for (ispec = 0; ispec < NANGLES; ispec++)
      fprintf (filep, "%-12e ", tau_spectrum[ispec * NBINS_SPEC + ifreq]);
    fprintf (filep, "\n");
    frequency += dfreq;
  }

  fflush (filep);

  if (fclose (filep))
  {
    Error ("%s : %i : could not close tau spectrum diag file\n", __FILE__, __LINE__);
    Exit (1);
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
    MPI_Reduce (MPI_IN_PLACE, spec, NBINS_SPEC * nspec, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  else
    MPI_Reduce (spec, spec, NBINS_SPEC * nspec, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
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
calculate_tau_across_cell (WindPtr w, PhotPtr p, double *col_den, double *tau)
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

  if ((p->grid = where_in_grid (w[p->grid].ndom, p->x)) < 0)
  {
    Error ("%s : %i : pextract is not in grid\n", __FILE__, __LINE__);
    return EXIT_FAILURE;
  }

  /*
   * Determine where in the plasma grid the cell is. This is required so a
   * column density can be calculated. Also calculate the maximum distance the
   * photon can traverse across the cell
   */

  wind_cell = &w[p->grid];
  ndom = wind_cell->ndom;
  nplasma = w[p->grid].nplasma;
  plasma_cell = &plasmamain[nplasma];
  density = plasma_cell->rho;
  smax = smax_in_cell (p);
  if (smax < 0)
  {
    Error ("%s : %i : abnormal value of smax for photon\n", __FILE__, __LINE__);
    return EXIT_FAILURE;
  }

  // Transform photon at starting location to local frame
  observer_to_local_frame (p, &p_start);

  // Move photon smax and transform to local frame
  stuff_phot (p, &p_stop);
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
    stuff_phot (p, &p_now);
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
    kappa_tot += radiation (p, smax);
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

  *col_den += smax * density;
  *tau += smax * kappa_tot;
  move_phot (p, smax);
  istat = p->istat;

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
extract_tau (WindPtr w, PhotPtr porig, double *col_den, double *tau)
{
  int ierr, ndom, istat;
  int nspace;
  const int max_space = 10;
  double norm[3];
  struct photon pextract;

  istat = P_INWIND;             // assume photon is in wind for initialisation reasons
  stuff_phot (porig, &pextract);

  nspace = 0;

  while (istat == P_INWIND)
  {
    // Move the photon in the grid, but not in the wind
    if (where_in_wind (pextract.x, &ndom) < 0)
    {
      translate_in_space (&pextract);
      if (++nspace > max_space)
      {
        Error ("%s : %i : tau_extract photon transport ended due to too many translate_in_space\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
      }
    }
    // Move the photon in the wind
    else if ((pextract.grid = where_in_grid (ndom, pextract.x)) >= 0)
    {
      ierr = calculate_tau_across_cell (w, &pextract, col_den, tau);
      if (ierr)
      {
        pextract.istat = -1;
        return EXIT_FAILURE;
      }
    }
    else
    {
      Error ("%s : %i : photon in unknown location, grid stat %i\n", __FILE__, __LINE__, pextract.grid);
      pextract.istat = P_ERROR;
      return EXIT_FAILURE;
    }

    istat = walls (&pextract, porig, norm);
  }

  if (istat == P_HIT_STAR || istat == P_HIT_DISK)
  {
    Error ("%s : %i : photon hit central source or disk incorrectly istat = %i\n", __FILE__, __LINE__, istat);
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
  int icell;
  int wind_index, plasma_index;
  static bool xloc_init = false;
  double x_loc;
  static double upward_x_loc;
  double plasma_density, max_plasma_density;

  /*
   * Hacky code for when photon is pointing upwards - this photon will not
   * be launched from the origin, but will be launched from the most dense
   * cell at the bottom of the disk.
   */

  if (xloc_init == false)
  {
    max_plasma_density = 0;
    for (icell = 0; icell < zdom[LAUNCH_DOMAIN].ndim - 1; icell++)
    {
      wind_index = icell * zdom[LAUNCH_DOMAIN].mdim + zdom[LAUNCH_DOMAIN].mdim;
      x_loc = wmain[wind_index].x[0];
      plasma_index = wmain[wind_index].nplasma;
      plasma_density = plasmamain[plasma_index].rho;
      if (plasma_density > max_plasma_density)
      {
        upward_x_loc = x_loc + wmain[wind_index].dfudge;
        max_plasma_density = plasma_density;
      }
    }

    xloc_init = true;
  }

  if (pout->lmn[0] == 0 && pout->lmn[1] == 0 && pout->lmn[2] == 1)
  {
    pout->x[0] = upward_x_loc;
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
create_tau_photon (PhotPtr pout, double nu, double *lmn)
{
  if (nu < 0)
  {
    Error ("%s : %i : photon can't be created with negative nu\n", __FILE__, __LINE__);
    return EXIT_FAILURE;
  }

  stuff_v (lmn, pout->lmn);
  pout->freq = pout->freq_orig = nu;
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
    NANGLES = geo.nangles;
    SIGHTLINES = calloc (geo.nangles, sizeof *SIGHTLINES);
    if (SIGHTLINES == NULL)
    {
      mem_req = geo.nangles * (int) sizeof *SIGHTLINES;
      Error ("%s : %i : cannot allocate %d bytes for observers array\n", __FILE__, __LINE__, mem_req);
      Exit (1);
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
    NANGLES = n_default_angles;
    SIGHTLINES = calloc (n_default_angles, sizeof *SIGHTLINES);
    if (SIGHTLINES == NULL)
    {
      mem_req = n_default_angles * (int) sizeof *SIGHTLINES;
      Error ("%s : %i : cannot allocate %d bytes for observers array\n", __FILE__, __LINE__, mem_req);
      Exit (1);
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
 * ************************************************************************** */

static void
create_optical_depth_spectrum (WindPtr w)
{
  int ierr;
  int ispec, ifreq;
  int nbins, mpi_lower, mpi_upper;
  double tau, column;
  double freq;
  double freq_min, freq_max, dfreq;
  double *current_observer;
  double *tau_spectrum;

  struct photon ptau;

  Log ("Creating optical depth spectra:\n");

  tau_spectrum = calloc (NANGLES * NBINS_SPEC, sizeof *tau_spectrum);
  if (tau_spectrum == NULL)
  {
    Error ("%s : %i : cannot allocate %d bytes for tau_spectrum\n", __FILE__, __LINE__, NANGLES * NBINS_SPEC * sizeof *tau_spectrum);
    return;
  }

  /*
   * Define the limits of the spectra in frequency space. If xxpsec is NULL,
   * then the frequency range will be over the generated photon bands. Otherwise,
   * the frequency range of the extracted spectrum is used.
   * TODO used logarithmically spaced frequency bins
   */

  if (xxspec == NULL)
  {
    freq_min = xband.f1[0];
    freq_max = xband.f2[xband.nbands - 1];
    NBINS_SPEC = (int) pow (10, log10 (freq_max / freq_min));
  }
  else
  {
    freq_min = VLIGHT / (geo.swavemax * ANGSTROM);
    freq_max = VLIGHT / (geo.swavemin * ANGSTROM);
  }

  if (NBINS_SPEC < NWAVE)
    NBINS_SPEC = NWAVE;

  dfreq = (freq_max - freq_min) / NBINS_SPEC;
  kbf_need (freq_min, freq_max);

  /*
   * Initialise the variables to split the calculation across multiple
   * processors
   */

  nbins = ceil ((double) NBINS_SPEC / np_mpi_global);
  mpi_lower = nbins * rank_global;
  mpi_upper = nbins * (rank_global + 1);
  if (mpi_upper > NBINS_SPEC)
    mpi_upper = NBINS_SPEC;

  /*
   * Now create the optical depth spectra for each sightline
   */

  for (ispec = 0; ispec < NANGLES; ispec++)
  {
    Log ("  - Creating spectrum: %s\n", SIGHTLINES[ispec].name);
    current_observer = SIGHTLINES[ispec].lmn;
    freq = mpi_lower * dfreq;
    if (freq < freq_min)
      freq = freq_min;

    for (ifreq = mpi_lower; ifreq < mpi_upper; ifreq++)
    {
      tau = 0;
      column = 0;
      ierr = create_tau_photon (&ptau, freq, current_observer);
      if (ierr == EXIT_FAILURE)
      {
        Log ("%s : %i : skipping photon of frequency %e when creating spectrum\n", __FILE__, __LINE__, freq);
        continue;
      }

      ierr = extract_tau (w, &ptau, &column, &tau);
      if (ierr == EXIT_FAILURE)
        continue;

      tau_spectrum[ispec * NBINS_SPEC + ifreq] = tau;
      freq += dfreq;
    }
  }

  mpi_gather_spectra (tau_spectrum, NANGLES);
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
calculate_photoion_edge_optical_depth (WindPtr w)
{
  int ierr;
  int itau, ispec;
  double nu, tau, column;
  double *current_observer;
  double *optical_depths, *column_densities;
  struct photon ptau;

  optical_depths = calloc (NANGLES * NTAU, sizeof *optical_depths);
  if (optical_depths == NULL)
  {
    Error ("%s : %i : cannot allocate %d bytes for optical_depths\n", __FILE__, __LINE__, NANGLES * NTAU * sizeof *optical_depths);
    return;
  }

  /*
   * Allocate storage for the column densities for each edge and sight line
   */

  column_densities = calloc (NANGLES, sizeof *optical_depths);
  if (column_densities == NULL)
  {
    Error ("%s : %i : cannot allocate %d bytes for column_densities\n", __FILE__, __LINE__, NANGLES * sizeof *optical_depths);
    free (optical_depths);
    return;
  }

  /*
   * Now extract the optical depths and mass column densities
   */

  for (ispec = 0; ispec < NANGLES; ispec++)
  {
    current_observer = SIGHTLINES[ispec].lmn;
    for (itau = 0; itau < NTAU; itau++)
    {
      tau = 0.0;
      column = 0.0;
      nu = EDGES[itau].nu;
      ierr = create_tau_photon (&ptau, nu, current_observer);
      if (ierr == EXIT_FAILURE)
      {
        Error ("%s : %i : skipping photon of frequency %e\n", __FILE__, __LINE__, nu);
        continue;
      }

      ierr = extract_tau (w, &ptau, &column, &tau);
      if (ierr == EXIT_FAILURE) // Do not throw another "extra" error
        continue;

      optical_depths[ispec * NTAU + itau] = tau;
      column_densities[ispec] = column;
    }
  }

  print_optical_depths (optical_depths, column_densities);
  free (optical_depths);
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
  calculate_photoion_edge_optical_depth (w);
  create_optical_depth_spectrum (w);

#ifdef MPI_ON
  MPI_Barrier (MPI_COMM_WORLD);
#endif

  free (SIGHTLINES);
  geo.ioniz_or_extract = 1;
  xsignal (files.root, "%-20s Optical depth diagnostics completed\n", "TAU");
  Log (" Completed optical depth diagnostics. The elapsed TIME was %f\n", timer ());
}
