/* ************************************************************************** */
/**
 *  @file    vr.c
 *
 *  @author  Edward Parkinson
 *
 *  @date    July 2019
 *
 *  @brief
 *
 *  Contains the various routines dedicated for variance reduction operations.
 *  The purpose of these routines are to improve the program's ability to deal
 *  with high optical depths as well as to constrain photon weights.
 *
 * ************************************************************************** */


#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "atomic.h"
#include "python.h"


/* ************************************************************************** */
/**
 *  @brief      Initialise the variance reduction variables and structure.
 *
 *  @details
 *
 *  Simply takes in parameters from the user, as well as allocates memory for
 *  the photon store structure.
 *
 * ************************************************************************** */

void
init_vr (void)
{
  char default_usage[LINELENGTH];
  size_t memory_requirement;

  strcpy (default_usage, "no");
  vr_configuration.on = rdchoice ("VR.use_packet_split_russian_roulette(yes,no)", "1,0", default_usage);

  /*
   * Set the default values for the number of photons which are created with
   * packet splitting, as well as the kill probability for Russian Roulette.
   * These are parameters which are also set by the user.
   */

  vr_configuration.ps_nsplit = 10;
  vr_configuration.rr_pkill = 0.25;
  rdint ("VR.packet_splitting_nsplit", &vr_configuration.ps_nsplit);
  rddoub ("VR.russian_roulette_pkill", &vr_configuration.rr_pkill);

  /*
   * If the photon store structure is NULL, then we need to allocate memory for
   * it. We will be re-using this structure over and over again, so we shouldn't
   * need to re-allocate as we limit the number of low weight photons which
   * can be produced by ps_nsplit.
   */

  if (vr_configuration.ps_photstore == NULL)
  {
    memory_requirement = vr_configuration.ps_nsplit * sizeof (p_dummy);
    if (!(vr_configuration.ps_photstore = calloc (sizeof (p_dummy), vr_configuration.ps_nsplit)))
    {
      Error ("%s : %i : could not allocate %d bytes for packet splitting photon store\n", __FILE__, __LINE__, memory_requirement);
      Exit (1);
    }

    Log_silent ("%d bytes allocated for %i low weight photons\n", memory_requirement, vr_configuration.ps_nsplit);
  }
}


/* ************************************************************************** */
/**
 *  @brief      Cleans up any mess left behind by variance reduction functions.
 *
 *  @details
 *
 *  Simply frees the memory used by the packet splitting photon store.
 *
 * ************************************************************************** */

void
clean_up_vr (void)
{
  free (vr_configuration.ps_photstore);
}

/* ************************************************************************** */
/**
 *  @brief     Split a high weight photon packet into multiple lower weight
 *             photon packets
 *
 *  @param[in,out]    struct photon *pin    The high weight photon to be split
 *  @param[in]        double tau_scat       The optical depth to the next
 *                                          interaction to be given to the
 *                                          child photons
 *
 *  @return       void
 *
 *  @details
 *
 *  Any high weight photon given to this function will be split into multiple
 *  low weight photons pointing in different directions, with the same
 *  frequency, etc. as the high weight photon. The only difference is in the
 *  fact that each low weight photon has an identical lower weight than the high
 *  weight photon.
 *
 *  Photons are stored in a PhotPtr array located in the vr_configuration
 *  structure. To free this memory, the function clean_up_vr () should be called
 *  as this is basically a wrapper for cleaning up any mess left behind by the
 *  VR techniques.
 *
 *  TODO: ability to track child photons
 *
 * ************************************************************************** */

void
split_photon_packet (struct photon *pin)
{
  int iphot;
  double new_weight;

  new_weight = pin->w / vr_configuration.ps_nsplit;

  for (iphot = 0; iphot < vr_configuration.ps_nsplit; iphot++)
  {
    stuff_phot (pin, &vr_configuration.ps_photstore[iphot]);
    vr_configuration.ps_photstore[iphot].w = vr_configuration.ps_photstore[iphot].w_orig = new_weight;
    vr_configuration.ps_photstore[iphot].daughter += 1;
    scatter (&vr_configuration.ps_photstore[iphot], &vr_configuration.ps_photstore[iphot].nres,
             &vr_configuration.ps_photstore[iphot].nnscat);
  }

  pin->istat = P_PS_SPLIT;
}

/* ************************************************************************** */
/**
 *  @brief    Play Russian Roulette with a photon packet which will either
 *            "kill" it or increase its weight
 *
 *  @param[in,out]   struct photon *pin   The photon to play Russian Roulette
 *                                        with
 *  @param[in]       double p_kill        The kill probability of Russian
 *                                        Roulette
 *
 *  @return          istat                The status of the photon
 *
 *  @details
 *
 *  Generates a random number, if this is less than or equal to p_kill then
 *  the weight of the photon is reduced to zero which (hopefully) should cause
 *  Python to cease tracking the photon. If, however, the photon "survives",
 *  then the weight of this photon is increased. This function should be called
 *  for low weight photons, but there is nothing to stop it from being used
 *  on high weight photons.
 *
 * ************************************************************************** */

int
play_russian_roulette (struct photon *pin, double p_kill)
{
  double xi;

  xi = random_number (0.0, 1.0);

  if (xi < p_kill)
  {
    pin->w = 0;
    pin->istat = P_RR_KILLED;
  }
  else
  {
    pin->w *= 1.0 / (1.0 - p_kill);
  }

  return pin->istat;
}

/* ************************************************************************** */
/**
 *  @brief      Assign the importance value for plasma grid cells.
 *
 *  @details
 *
 *
 * ************************************************************************** */

void
update_importance_map (void)
{
  int iplasma;
  PlasmaPtr pcell;

  for (iplasma = 0; iplasma < NPLASMA; iplasma++)
  {
    pcell = &plasmamain[iplasma];
    pcell->importance = 1;
  }
}


/* ************************************************************************** */
/**
 *  @brief
 *
 *  @details
 *
 *  TODO: debug code to be removed
 *
 * ************************************************************************** */

void
vr_debug_print_weights (void)
{
  int i;
  char fname1[LINELENGTH];
  FILE *f1;
  struct photon p;
  double tot_weight = 0;
  double tot_orig_weight = 0;

  if (rank_global != 0)
    return;

  sprintf (fname1, "%s_%i_photon_w.txt", files.root, rank_global);
  f1 = fopen (fname1, "w");
  if (!f1)
  {
    Error ("Unable to open file to print photon weights\n");
    Exit (1);
  }

  fprintf (f1, "n np w w0 istat\n");

  int icount = 0;
  for (i = 0; i < NPHOT; i++)
  {
    p = photmain[i];
    fprintf (f1, "%i %i %e %e ", icount, p.np, p.w, p.w_orig);
    if (p.istat == P_ESCAPE)
    {
      fprintf (f1, "P_ESCAPE\n");
      icount++;
      tot_weight += p.w;
      tot_orig_weight += p.w_orig;
    }
    else if (p.istat == P_ABSORB)
    {
      fprintf (f1, "P_ABSORB\n");
      icount++;
    }
    else if (p.istat == P_HIT_STAR)
    {
      fprintf (f1, "P_HIT_STAR\n");
      icount++;
    }
    else
    {
      fprintf (f1, "OTHER\n");
      icount++;
    }
  }

  fclose (f1);

  Log ("Total weight escaped %e\n", tot_weight);
  Log ("Total original weight of escaped photons %e\n", tot_orig_weight);
  Log ("tot_orig_weight - tot_weight =  %e\n", tot_orig_weight - tot_weight);
}

/* ************************************************************************** */
/**
 *  @brief      Generate a value for the optical depth to the next interaction
 *              event
 *
 *  @param[in]  struct photon *pin    The photon packet to generate tau_scat for
 *
 *  @return     double                The generated optical depth, tau_scat
 *
 *  @details
 *
 *  Generates a value for tau_scat, i.e. the optical depth to the next
 *  interaction site in the wind.
 *
 *  ### Notes ###
 *
 *  This function is sampling from the exponential distribution,
 *
 *      p(tau) = e^(-tau)
 *
 *  using the inversion method.
 *
 *  The photon taken in as an argument is unsed. The function is defined this
 *  way such that a function pointer can be used to point to the unbiased
 *  and biased distribution used in path stretching.
 *
 * ************************************************************************** */

/*
double
generate_tau_scat (struct photon *pdummy)
{
  double tau_scat;

  tau_scat = -log (1.0 - random_number (0.0, 1.0));
  return tau_scat;
}
*/

/* ************************************************************************** */
/**
 *  @brief      Re-weight a photon packet which has taken a stretched path
 *
 *  @params[in,out]   struct photon *pin      The photon which took a stretched
 *                                            path
 *  @params[in]       double tau_scat         The stretched optical depth sample
 *  @params[in]       double alpha            The stretching parameter alpha
 *
 *  @return           void
 *
 *  @details
 *
 *  Re-weights a photon which has taken a stretched path via sampling an optical
 *  depth to the next interaction using a stretched exponential distribution.
 *
 *  The change in weight is given by,
 *
 *      W_new = W_old * p(tau) / q(tau)
 *
 *  where p(tau) = e^(-tau) and q(tau) = alpha * e^(-alpha * tau). This
 *  expression simplifies into the code below,
 *
 *      W_new = W_old * e^(tau * (alpha - 1)) / alpha
 *
 * ************************************************************************** */

/*
void
reweight_biased_photon (struct photon *pin, double tau_scat, double alpha)
{
  pin->w *= exp (tau_scat * (alpha - 1.0)) / alpha;
}
*/

/* ************************************************************************** */
/**
 *  @brief      Generate a value for the optical depth to the next interaction
 *              event from a biased distribution
 *
 *  @param[in,out]  struct photon *pin    The photon packet to generate tau_scat
 *
 *  @return         double                The generated optical depth, tau_scat
 *
 *  @details
 *
 *  Generates a value for tau_scat, i.e. the optical depth to the next
 *  interaction site in the wind. The photon packet is then re-weighted to
 *  account for sampling from the biased distribution.
 *
 *  ### Notes ###
 *
 *  This function is sampling from the exponential distribution,
 *
 *      p(tau) = alpha * e^(-alpha * tau)
 *
 *  using the inversion method. Alpha here is the stretching parameter and takes
 *  the value 0 < alpha < 1. Alpha can either be set as a fixed value throughout
 *  the grid, or can be calculated in the fly. Either way, the current value
 *  of alpha to use is found by querying the plasma grid cell the photon packet
 *  is currently within.
 *
 * ************************************************************************** */

/*
double
generate_biased_tau_scat (struct photon *pin)
{
  double alpha;
  double tau_scat;
  double in_weight;

  in_weight = pin->w;
  alpha = plasmamain[wmain[pin->grid].nplasma].path_stretch_alpha;
  tau_scat = -log (1.0 - random_number (0.0, 1.0)) / alpha;
  reweight_biased_photon (pin, tau_scat, alpha);

  if (pin->w < EPSILON * in_weight)
    play_russian_roulette (pin, vr_configuration.rr_pkill);

  if (pin->w > in_weight / EPSILON)
    split_photon_packet (pin, vr_configuration.ps_nsplit);

  return tau_scat;
}
*/
