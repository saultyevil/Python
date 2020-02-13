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
init_russian_roulette (void)
{
  char default_usage[LINELENGTH];

  strcpy (default_usage, "no");
  vr_configuration.use_russian_roulette = rdchoice ("VR.use_russian_roulette(yes,no)", "1,0", default_usage);

  if (vr_configuration.use_russian_roulette)
  {
    /*
     * TODO
     * kill probability for russian roulette
     */

    vr_configuration.rr_pkill = 0.1;
    rddoub ("VR.russian_roulette_pkill", &vr_configuration.rr_pkill);

    /*
     * TODO
     * the optical depth in the cell to consider using RR with, this should be
     * an advanced diag probably
     */

    vr_configuration.rr_tau_crit = 20;
    rddoub ("VR.russian_roulette_tau_crit", &vr_configuration.rr_tau_crit);

    /*
     * TODO
     * Don't really need this
     */

    strcpy (default_usage, "no");
    vr_configuration.debug_messages = rdchoice ("VR.debug_messages(yes,no)", "1,0", default_usage);
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
    fprintf (f, "np w_i w\n");
    init++;
  }

  fprintf (f, "%i %e %e\n", p->np, p->w, p->w_rr_orig);
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
