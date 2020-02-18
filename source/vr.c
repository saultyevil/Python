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
init_variance_reduction_optimisations (void)
{
  char default_usage[LINELENGTH];

  /*
   * Initialise the Russian Roulette variables
   */

  strcpy (default_usage, "no");
  RussianRoulette.enabled = rdchoice ("RR.enable(yes,no)", "1,0", default_usage);

  if (RussianRoulette.enabled)
  {
    RussianRoulette.kill_probability = 0.1;
    rddoub ("RR.kill_probability", &RussianRoulette.kill_probability);
    if (RussianRoulette.kill_probability < 0 && RussianRoulette.kill_probability > 1)
    {
      Error ("%s : %i: the kill probability has to be 0 < pkill < 1\n", __FILE__, __LINE__);
      Exit (1);
    }

    RussianRoulette.weight_limit = 1e3;
    rddoub ("RR.weight_limit(original_photon_weight)", &RussianRoulette.weight_limit);
    if (RussianRoulette.weight_limit <= 1)
    {
      Error ("%s : %i: the upper weight limit of photon packets has to be more than 1\n", __FILE__, __LINE__);
      Exit (1);
    }

    RussianRoulette.critcal_tau = 20;
    rddoub ("RR.critical_optical_depth", &RussianRoulette.critcal_tau);
    if (RussianRoulette.critcal_tau < 0)
    {
      Error ("%s : %i: the critical optical depth to play Russian Roulette has to be positive and non-zero\n", __FILE__, __LINE__);
      Exit (1);
    }

    /*
     * TODO
     * Don't really need this when fully implemented
     */

    strcpy (default_usage, "no");
    RussianRoulette.debug_messages = rdchoice ("RR.debug_messages(yes,no)", "1,0", default_usage);
  }

  /*
   * Initialise the Packet Splitting variables -- reinitialise default_usage
   * just in case it has been modified
   */

  strcpy (default_usage, "no");
  PacketSplitting.enabled = rdchoice ("PS.enable(yes,no)", "1,0", default_usage);

  if (PacketSplitting.enabled)
  {
    PacketSplitting.nsplit = 5;
    rdint ("PS.max_low_weight_photons_to_create", &PacketSplitting.nsplit);
    if (PacketSplitting.nsplit <= 1)
    {
      Error ("%s : %i : need to create at least 2 low weight photons when splitting\n", __FILE__, __LINE__);
      Exit (1);
    }

    PacketSplitting.critical_tau = 1;
    rddoub ("PS.critical_tau", &PacketSplitting.critical_tau);
    if (PacketSplitting.critical_tau <= 0)
    {
      Error ("%s : %i: the critical optical depth to split photons at has to be positive and non-zero\n", __FILE__, __LINE__);
      Exit (1);
    }

    PacketSplitting.weight_limit = 5e2;
    if (PacketSplitting.weight_limit <= 0)
    {
      Error ("%s : %i: the weight threshold to split photons has to be positive and non-zero\n", __FILE__, __LINE__);
      Exit (1);
    }

    /*
     * Allocate memory for the storage of low weight photons. We can't create
     * more photons than the value nsplit, as this is the maximum we have
     * allowed to be created.
     */

    PacketSplitting.photons = calloc (PacketSplitting.nsplit, sizeof *PacketSplitting.photons);
    if (!PacketSplitting.photons)
    {
      unsigned long mem = PacketSplitting.nsplit * sizeof *PacketSplitting.photons;
      Error ("%s : %i : unable to allocate storage (%d bytes) for low weight photons\n", __FILE__, __LINE__, mem);
      Exit (1);
    }

    /*
     * TODO to be remove at later date
     */

    strcpy (default_usage, "no");
    PacketSplitting.debug_messages = rdchoice ("PS.debug_messages(yes,no)", "1,0", default_usage);
  }
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
 *  @brief      Split an input, high weight, photon packet into multiple low
 *              weight photon packets.
 *
 *  @param[in]  PhotPtr p    The photon which is to be split into multiple
 *                           low weight photon packets
 *
 *  @param[in]  int nsplit   The number of low weight photon packets to create
 *
 *  @returns    EXIT_FAIL or EXIT_SUCESS
 *
 *  @details
 *
 *  The purpose of this function is to create nsplit low weight photon packets
 *  which are child packets of the original p photon packet. The created photon
 *  packets are identical, apart from they all have different directions.
 *
 * ************************************************************************** */

int
split_photon_packet (PhotPtr p, int nsplit)
{
  int i;
  double new_weight;
  struct photon new_photon;

  if (nsplit > PacketSplitting.nsplit)
  {
    Error ("%s : %i : nsplit %i exceeds maximum allowed split photons %i. Not splitting photon\n", __FILE__, __LINE__, nsplit,
           PacketSplitting.nsplit);
    return EXIT_FAILURE;
  }

  /*
   * Check if this is a split photon, because we do not want to split this again
   * and create situations where we are continuously creating low weight photon
   * packets.
   */

  if (p->split_child)
    return EXIT_SUCCESS;

  /*
   * Now create the low weight photons and mark them as being children
   */

  new_weight = p->w / nsplit;

  for (i = 0; i < nsplit; ++i)
  {
    stuff_phot (p, &new_photon);
    new_photon.w = new_photon.w_orig = new_weight;
    scatter (p, &p->nres, NULL);
    new_photon.split_child = TRUE;
    PacketSplitting.photons[i] = new_photon;
  }

  p->istat = P_PS_SPLIT;

  return EXIT_SUCCESS;
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
record_photon (PhotPtr p)
{
  char filename[LINELENGTH];
  static int init = 0;
  static FILE *f;

  if (rank_global != 0)
    return;

  if (!init)
  {
    sprintf (filename, "%s_rr_diag.log.txt", files.root);
    f = fopen (filename, "w");
    if (!f)
    {
      Error ("Couldn't open RR diag log file???\n");
      Exit (1);
    }
    fprintf (f, "np w_i w istat\n");
    init++;
  }

  fprintf (f, "%i %e %e %i\n", p->np, p->w_rr_orig, p->w, p->istat);
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

  for (i = 0; i < NPHOT; i++)
  {
    p = photmain[i];
    fprintf (f1, "%i %i %e %e ", i, p.np, p.w, p.w_orig);
    if (p.istat == P_ESCAPE)
    {
      fprintf (f1, "P_ESCAPE\n");
      tot_weight += p.w;
      tot_orig_weight += p.w_orig;
    }
    else if (p.istat == P_ABSORB)
    {
      fprintf (f1, "P_ABSORB\n");
    }
    else if (p.istat == P_HIT_STAR)
    {
      fprintf (f1, "P_HIT_STAR\n");
    }
    else if (p.istat == P_RR_KILLED)
    {
      fprintf (f1, "P_RR_KILLED\n");
    }
    else
    {
      fprintf (f1, "OTHER\n");
    }
  }

  fclose (f1);

  Log ("Total weight of escaped photons             %e\n", tot_weight);
  Log ("Total original weight of escaped photons    %e\n", tot_orig_weight);
  Log ("tot_orig_weight - tot_weight =  %e\n", tot_orig_weight - tot_weight);
}
