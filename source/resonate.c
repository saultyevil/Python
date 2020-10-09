
/***********************************************************/
/** @file  resonate.c
 * @author ksl
 * @date   May, 2018
 *
 * @brief
 * These routines have to do with calculating where resonances are and
 * how strong they are.
 *
 * ###Notes###
 * These routines should be made geometry independent and therefore
 * should work with
 * very minor (header only) modifications in other configurations.
 *
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

struct photon cds_phot_old_observer, cds_phot_old_loc;
//OLD double cds_v2_old, cds_dvds2_old;


/**********************************************************/
/**
 * @brief     calculate the distance a photon can travel
 * within a single shell without scattering
 *
 * @param [in] WindPtr  w   the entire wind structure
 * @param [in] PhotPtr  p   A photon (bundle)
 * @param [in] double  tau_scat   the optical depth at which the photon
 * will scatter
 * @param [in,out] double *  tau   Initially the current optical depth for
 * the photon; finally the optical depth at the distance the photon can be
 * moved.
 * @param [out] int *  nres   the number of the resonance, -1 if it was electron 
 * scatttering, -2 if it was ff, -99 if the
 * distance is not limited by some kind of scatter.
 * @param [in] double  smax   the maximum distance the photon can
 * travel in the cell
 * @param [out] int *  istat   A flag indicating whether the
 * photon should scatter if it travels the distance estimated, 0 if no, TAU_SCAT if yes.
 * @return     The distance the photon can travel
 *
 * calculate_ds finds the distance the photon can travel subject to a
 * number of conditions, which include reach a distance where tau is the
 * distance where a scatter can occur, or a maximum distance set to ensure
 * we do not cross into another cell, or indeed to go so far that one is
 * not sure that the velocity can be approximated as a linear function of
 * distance.
 *
 * The routine returns the distance that the photon can travel subject to
 * these conditions and information about why the photon was halted at the
 * distance it was halted.
 *
 * @details
 *
 * ### Notes ###
 * Calculate_ds does not modify the p in any way!!
 *
 * Any paramaters that depend explicitly on the
 * coordinate gridding, such as the maximum distance the photon
 * can travel before hitting the edge of the
 * shell should be calculated outside of this routine.
 *
 * 180519 - ksl - I modified the behavior ot this routine so that nres
 * indicates that a scattering event is not limiting the maximum
 * distance.  nres is now -99 if there was no scatter, instead of
 * -1
 *
 *
 **********************************************************/
double
calculate_ds (w, p, tau_scat, tau, nres, smax, istat)
     WindPtr w;
     PhotPtr p;
     double tau_scat, *tau;
     int *nres;
     double smax;
     int *istat;
{
  int kkk;
  double kap_es;
  double freq_inner, freq_outer, dfreq, ttau, freq_av;
  double mean_freq;             //A mean freq used in compton calculations.
  int n, nn, nstart, ndelt;
  double x;
  double ds_current, ds;
  double dvds, dd;
  double dvds1, dvds2;
  struct photon p_start, p_stop, p_now;
  int init_dvds;
  double kap_bf_tot, kap_ff, kap_cont;
  double tau_sobolev;
  WindPtr one, two;
  int check_in_grid;
  int nplasma;
  PlasmaPtr xplasma, xplasma2;
  int ndom;
  double normal[3];
  double diff;

  one = &w[p->grid];            //pointer to the cell where the photon bundle is located.

  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];
  ndom = one->ndom;

  ttau = *tau;
  ds_current = 0;
  init_dvds = 0;
  dvds1 = dvds2 = 0.0;
  *nres = -1;
  *istat = P_INWIND;

  if (ttau < 0.0)
  {
    Error ("calculate_ds:  odd tau  %8.2e at %g entering calculate_ds\n", ttau, p->freq);
  }


  if (comp_phot (&cds_phot_old_observer, p))
  {
    observer_to_local_frame (p, &p_start);
  }
  else
  {
    stuff_phot (&cds_phot_old_loc, &p_start);
//OLD  Log ("This helped  %d\n", p_start.np);
  }

  stuff_phot (p, &p_stop);
  move_phot (&p_stop, smax);
  stuff_phot (&p_stop, &cds_phot_old_observer);
  observer_to_local_frame (&p_stop, &p_stop);
  stuff_phot (&p_stop, &cds_phot_old_loc);



  /* At this point p_start and pstop are in the local frame 
   * at the and p_stop is at the maximum distance it can 
   * travel.  We want to check that the frequncy shift is 
   * not too great along the path that a linear approximation
   * to the change in frequency is not reasonable
   */


  diff = 1;
#define  MAXDIFF  VCHECK/VLIGHT /* The same as our old velocity requirement */

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
//OLD    Log ("Shorten %d %d %e %e\n", p_now.np, p_now.grid, smax, diff);
  }

  freq_inner = p_start.freq;
  freq_outer = p_stop.freq;



/* We use the doppler shifted frequency to compute the Klein-Nishina cross
 * section, if the frequency is high enough, otherwise we just use the
 * Thompson cross section.  For the time being, use the average frequency.
 * If we want true fidelity, perhaps we could compute the cross section
 * for every little path section between resonances */

  mean_freq = 0.5 * (freq_inner + freq_outer);
  dfreq = freq_outer - freq_inner;



/* The next section checks to see if the frequency difference on
 * the two sides is very small and if not limits the resonances
 * one has to worry about
 */

  if (fabs (dfreq) < EPSILON)
  {
    Error ("calculate_ds: freq same at both sides of cell %d\n", one->nwind);
    x = -1;
    return (smax);              // This is not really the best thing to do, but it avoids disaster below

  }
  else if (dfreq > 0)
  {
    limit_lines (freq_inner, freq_outer);
    nstart = nline_min;
    ndelt = 1;
  }
  else
  {
    limit_lines (freq_outer, freq_inner);
    nstart = nline_max;
    ndelt = (-1);
  }

/*nline_min, nline_max, and nline_delt are found in atomic.h and
 * are set by limit_lines()
 */

  /*Compute the angle averaged electron scattering cross section.  Note electron scattering  is always
     treated as a scattering event. */

  //XFRAME In principle, lines like this need a transformation: kappa transforms
  kap_es = klein_nishina (mean_freq) * xplasma->ne * zdom[ndom].fill;

/* If in macro-atom mode, calculate the bf and ff opacities, becuase in macro-atom mode
 * everthing including bf is calculated as a scattering process.  The routine
 * kappa_bound stores the individual opacities as well as the total, because when
 * there is more than one opacity contributin to the total, these are needed to choose
 * which particular bound-free transition to activate.  
 * For the two level approximation, none of this needed. 
 */

  kap_bf_tot = 0;
  kap_ff = 0;


  if (geo.rt_mode == RT_MODE_MACRO)
  {

    freq_av = freq_inner;

    // XFRAME need to do better than this perhaps but okay for star - comoving frequency (SS)
    //(freq_inner + freq_outer) * 0.5;  


    kap_bf_tot = kappa_bf (xplasma, freq_av, 0);
    kap_ff = kappa_ff (xplasma, freq_av);

    /* Okay the bound free contribution to the opacity is now sorted out (SS) */
  }

  if (one->vol == 0)
  {
    kap_bf_tot = kap_ff = 0.0;
    Error_silent ("ds_calculate vol = 0: cell %d position %g %g %g\n", p->grid, p->x[0], p->x[1], p->x[2]);
  }




  kap_cont = kap_es + kap_bf_tot + kap_ff;      //total continuum opacity
//XFRAME Probably here to transformation on all opacities for continuum, once summed up? ds (below) is observer frame

/* Finally begin the loop over the resonances that can interact
 * with the photon in the cell
 */

  for (n = 0; n < nline_delt; n++)
  {
    nn = nstart + n * ndelt;    /* So if the frequency of resonance increases as we travel through
                                   the grid cell, we go up in the array, otherwise down */
    x = (lin_ptr[nn]->freq - freq_inner) / dfreq;

    if (0. < x && x < 1.)
    {                           /* this particular line is in resonance */
      ds = x * smax;


/* Before checking for a resonant scatter, need to check for scattering
 * due to a continuum
 * process.
 */
      if (ttau + (kap_cont) * (ds - ds_current) > tau_scat)
      {
/* then the photon was scattered by the continuum before reaching the
 * resonance.  Need to randomly select the continumm process which caused
 * the photon to scatter.  The variable threshold is used for this. */

        *nres = select_continuum_scattering_process (kap_cont, kap_es, kap_ff, xplasma);
        //XFRAME This call will depend on kap_cont, kap_es being consistently in the same frame
        *istat = P_SCAT;        //flag as scattering
        ds_current += (tau_scat - ttau) / (kap_cont);   //distance travelled
        ttau = tau_scat;
        *tau = ttau;
        return (ds_current);
      }
      else
      {

/* increment tau by the continuum optical depth to this point */
        ttau += kap_cont * (ds - ds_current);
        //XFRAME could also transform kap_cont here instead of above


        ds_current = ds;        /* At this point ds_current is exactly the position of the resonance */
        kkk = lin_ptr[nn]->nion;


/* The density is calculated in the wind array at the center of a cell.
 * We use that as the first estimate of the density.  */
        stuff_phot (p, &p_now);
        move_phot (&p_now, ds_current); // So p_now contains the current position of the photon



        dd = get_ion_density (ndom, p_now.x, kkk);

        if (dd > LDEN_MIN)
        {
/* If we have reached this point then we have to initalize dvds1 and dvds2.
 * Otherwise there is no need to do this, especially as dvwind_ds is an
 * expensive calculation time wise */

          if (init_dvds == 0)
          {
            dvds1 = dvwind_ds (p);
            dvds2 = dvwind_ds (&p_stop);
            init_dvds = 1;
          }

          dvds = (1. - x) * dvds1 + x * dvds2;



          tau_sobolev = sobolev (one, p->x, dd, lin_ptr[nn], dvds);
          //XFRAME I think tau_sobolev is invariant, but will depend on all the ingredients being in the same frame

/* tau_sobolev now stores the optical depth. This is fed into the next statement for the bb estimator calculation. SS March 2004 */

          ttau += tau_sobolev;


          if (geo.rt_mode == RT_MODE_MACRO)     //Macro Atom case (SS)
          {

/* Because push through distance may take us out of the cell we want,
 * need to make sure that the cell is correct before incrementing the
 * heating rate/estimators. So 1st check if it's still in the wind and
 * second get a pointer to the grid cell where the resonance really happens.
 */

            check_in_grid = walls (&p_now, p, normal);

            if (check_in_grid != P_HIT_STAR && check_in_grid != P_HIT_DISK && check_in_grid != P_ESCAPE)
            {
              /* The next line may be redundant.  */
              two = &w[where_in_grid (wmain[p_now.grid].ndom, p_now.x)];

              if (two->vol == 0)
              {
                /* See issue #389 - Sometimes DFUDGE pushes a photon into a cell with no volume.  Note that this
                 * should be very rare, so if this error occurs in significant numbers the problem should be
                 * investigated further.  ksl -180626
                 */
                Error ("calculate_ds: Macro atom problem when photon moved into cell with no volume\n");
              }
              else if (geo.ioniz_or_extract == 1)
              {
                /* we only want to increment the BB estimators if we are in the ionization cycles (e.g. #730) */

                if (lin_ptr[nn]->macro_info == 1 && geo.macro_simple == 0)
                {
                  /* The line is part of a macro atom so increment the estimator if desired */
                  bb_estimators_increment (two, p, tau_sobolev, dvds, nn);
                  /* XFRAME ultimately need to consider transformaions in estimators
                     either transform quantities before passing (p) or else in the increment 
                     routines (choose a principle and apply to all?) */
                }
                else
                {
                  /* The line is from a simple ion. Record the heating contribution and move on. */
                  xplasma2 = &plasmamain[two->nplasma];

                  bb_simple_heat (xplasma2, p, tau_sobolev, dvds, nn);
                  /* XFRAME ultimately need to consider transformaions in estimators */

                }
              }
            }
          }
          /* Completed special calculations for the Macro Atom case */


        }



        /* Check to see whether the photon should scatter at this point */
        if (ttau > tau_scat)
        {
          *istat = P_SCAT;
          *nres = nn;
          *tau = ttau;




          return (ds_current);
        }

        /* End of loop to process an individual resonance */
      }
      *tau = ttau;
    }
  }



/* If the photon reaches this point it was not scattered by resonances.
 * ds_current is either 0 if there were no resonances or the position of the
 * "last" resonance if there were resonances.  But we need to check one
 * last time to see if it was scattered by continuum process.
 * Note: ksl -- It is generally a bad policy to have a bunch of repeated code
 * like this.  We should probably rewrite some of this in terms of subroutines
 * for clarity, especially the bit that calculates where the continuum
 * scattering event occurred.  04 apr
 */

  if (ttau + kap_cont * (smax - ds_current) > tau_scat)
  {
    *nres = select_continuum_scattering_process (kap_cont, kap_es, kap_ff, xplasma);
    //XFRAME as above -- need consistency of frames for kappas

    /* A scattering event has occurred in the shell  and we
     * remain in the same shell */
    ds_current += (tau_scat - ttau) / (kap_cont);
    *istat = P_SCAT;            /* Flag for scattering (SS) */
    ttau = tau_scat;
  }
  else
  {                             /* Then we did hit the other side of the shell
                                   (or possibly the another wall of the same shell) */
    *istat = P_INWIND;
    ttau += kap_cont * (smax - ds_current);     /* kap_es replaced with kap_cont (SS) */
    ds_current = smax;

  }

  *tau = ttau;

//OLD   stuff_phot (&phot, &cds_phot_old);    // Store the final photon position

//OLD   cds_v2_old = v2;              // and the velocity along the line of sight

  return (ds_current);
}



/**********************************************************/
/**
 * @brief      determine what the continuum
 * process was that caused the scatter and returns this to
 * the main routine.
 *
 * @param [in] double  kap_cont   The continuum opacity
 * @param [in] double  kap_es   The electron scattering opacity
 * @param [in] double  kap_ff   The free free opacity
 * @param [in] PlasmaPtr  xplasma   The plasma cell where everything is being
 * calculated
 * @return     The process that cause the photon to stop/scatter at a particular
 * point
 *
 *  * -1 implies electron scattering
 *  * -2 implies free-free
 *  * 0 or greater implies a specific photionization process was responsible (and 
 *  also that the program was operating in macro-atom mode.
 *
 * @details
 * This routine is called to determine which of several continuum proceseses
 * cause a photon to be scattered or absorbed.  In addition to electron scattering
 * and free-free absorption, the routine can identify which photoionization process
 * is implicated.
 *
 * ### Notes ###
 *
 * select_continuum_scattering_proces is called when it has already been
 * determined that a photon will stop or be scattered due to some continuum
 * process.  This determination is based on the total continuum opacity.  In this
 * routine a random number is generated and this is used to determine
 * which of the processes was responsible.  Data for the opacity due to
 * photonionization is passed remotely via the PlasmaPtr.
 *
 * In a program running in the two level approximation, only electron scattering
 * and ff and bf are treated as absorption processes.  In macro atom, ff and
 * photoionization are treated as a scattering 
 * processes.
 *
 **********************************************************/

int
select_continuum_scattering_process (kap_cont, kap_es, kap_ff, xplasma)
     double kap_cont, kap_es, kap_ff;
     PlasmaPtr xplasma;
{
  int nres;
  double threshold;
  double run_tot;
  int ncont;

  //XFRAME think this works provided that all the incoming kappas are in the same frame. shouldn't matter which frame though really - i.e. just consistency

  threshold = random_number (0.0, 1.0) * (kap_cont);


  /* First check for electron scattering. */
  if (kap_es > threshold)
  {                             /* electron scattering event occurred (SS) */
    nres = -1;                  // flag electron scatterin (SS) */
  }
  /* Now chech for ff. */
  else if ((kap_es + kap_ff) > threshold)
  {
    nres = -2;
  }
  /* Now check for bf. */
  else
  {
/* use a running sum to find which photoionisation process it was */
/*
If a non-macro-atom run is being done this part should never be reached.
Just do a check that all is well - this can be removed eventually (SS)
*/
    if (geo.rt_mode == RT_MODE_2LEVEL)
    {
      Error ("calculate_ds: Not using macro atoms but trying to excite one? Aboort.\n");
      Exit (0);
    }
    run_tot = kap_es + kap_ff;
    ncont = 0;
    while (run_tot < threshold)
    {
      run_tot += kap_bf[ncont];
      ncont++;
    }
    /* When it gets here know that excitation is in photoionisation labelled by ncont */
    nres = NLINES + 1 + xplasma->kbf_use[ncont - 1];    //modified SS Nov 04
  }
  return (nres);
}



/**********************************************************/
/**
 * @brief      calculates the bf opacity in a specific
 * 	cell.
 *
 * @param [in] PlasmaPtr  xplasma   The plasma cell of interest
 * @param [in] double  freq   The frequency at which the opacity is calculated
 * @param [in] int  macro_all   1--> macro_atoms only, 0 all topbase ions
 * @return     The bf opacity
 *
 * @details
 *
 * The routine calculates the bf 
 *
 * ### Notes ###
 * The routine allows for clumping, reducing kappa_bf by the filling
 * factor.
 *
 **********************************************************/

double
kappa_bf (xplasma, freq, macro_all)
     PlasmaPtr xplasma;
     double freq;
     int macro_all;


{
  double kap_bf_tot;
  double ft;
  double x;
  int nconf;
  double density;
  int n;
  int nn;
  int ndom;

  //XFRAME this will be giving a cmf kappa, I think - provided that level populations are densities in cmf (should be)

  kap_bf_tot = 0;

  macro_all--;                  // Subtract one from macro_all to avoid >= in for loop below.


  ndom = wmain[xplasma->nwind].ndom;

  for (nn = 0; nn < xplasma->kbf_nuse; nn++)    // Loop over photoionisation processes.
  {
    n = xplasma->kbf_use[nn];
    ft = phot_top[n].freq[0];   //This is the edge frequency (SS)

    kap_bf[nn] = 0.0;

    if (freq > ft && freq < phot_top[n].freq[phot_top[n].np - 1] && phot_top[n].macro_info > macro_all)
    {
      /* Need the appropriate density at this point. */

      nconf = phot_top[n].nlev; //Returning lower level = correct (SS)

      //XXX For Debugging
      if (nconf < 0)
      {
        Error ("kappa_bf: nconf %d for phot_top %d for ion %d of z %d and istate %d\n", nconf, n, phot_top[n].nion, phot_top[n].z,
               phot_top[n].istate);
        continue;
      }
      //XXX For Debugging

      density = den_config (xplasma, nconf);    //Need to check what this does (SS)


      if (density > DENSITY_PHOT_MIN || phot_top[n].macro_info == 1)
      {

        /* JM1411 -- added filling factor - density enhancement cancels with zdom[ndom].fill */
        kap_bf[nn] = x = sigma_phot (&phot_top[n], freq) * density * zdom[ndom].fill;   //stimulated recombination? (SS)
        kap_bf_tot += x;
      }
    }
  }


  return (kap_bf_tot);
}


/**********************************************************/
/**
 * @brief      computes and stores the set of bound-free processes which
 *         make significiant contributions to the opacity in a grid cell.
 *
 * @param [in] double  fmin   the lowest frequency of interest to the calculation
 * @param [in] double  fmax   the highest frequency of interest to the calculations
 * @return     Always returns 0
 *
 * The purpose of this routine is to speed up calculations by idenfifying which
 * bound-free x-sections are important enough to be included when calculationg the
 * bound-fee opacity, and which an be ignored because the density of the particular
 * ion is so low it will not contribute.
 *
 * For each cell, the routine determines what bf transitons are important
 * and stores them in one->kbf_use[n].
 *
 * The total number of such transitions
 * is given in one->kbf_nuse.
 *
 * @details
 *
 * This routine is now called before various cylces of the Monte Carlo calculation.
 * (in run.c) * It determines which
 * bf processes are worth considering during the calculation that follows.
 *
 * To do this
 * is uses the edge position (must be inside, or close to, the spectral region of
 * interest) and the edge opacity (must be greater than a threshold value, currently
 * 10^-6.
 *
 * The need for the routine is just to prevent wasting time on unimportant bf
 * transitions.  It is called (currently from resonate).
 *
 * ### Notes ###
 * The dimensionalty of kbf_use is currently set to NTOP_PHOT, which is large enough
 * to store every transition if necessary.
 *
 **********************************************************/

int
kbf_need (fmin, fmax)
     double fmin, fmax;


{
  int nconf;
  double density;
  double tau_test, ft;
  int n;
  int nuse;

  PlasmaPtr xplasma;
  WindPtr one;
  int nplasma, nion;
  //XFRAME this will be giving a cmf kappa, I think - provided that level populations are densities in cmf (should be)


  for (nplasma = 0; nplasma < NPLASMA; nplasma++)       // Loop over all the cells in the wind
  {
    xplasma = &plasmamain[nplasma];
    one = &wmain[xplasma->nwind];
    nuse = 0;

    for (n = 0; n < nphot_total; n++)   // Loop over photoionisation processes.
    {

      ft = phot_top[n].freq[0]; //This is the edge frequency (SS)

      if ((ft > (fmin / 3.)) && (ft < fmax))
      {
        nion = phot_top[n].nion;

        if (ion[nion].phot_info == 0)   // vfky
        {
          density = xplasma->density[nion];
        }
        // else if (ion[nion].phot_info > 0)  // topbase or hybrid
        else
        {
          nconf = phot_top[n].nlev;     //Returning lower level = correct (SS)

          // XXX This if for debugging 
          if (nconf < 0)
          {
            Error ("kbf_need: nconf %d for phot_top %d for ion %d of z %d and istate %d\n", nconf, n, nion, phot_top[n].z,
                   phot_top[n].istate);
          }
          // XXX This if for debugging 

          density = den_config (xplasma, nconf);
        }

        tau_test = phot_top[n].x[0] * density * SMAX_FRAC * length (one->xcen);


        if (tau_test > 1.e-6 || phot_top[n].macro_info == 1)
        {
          /* Store the bf transition and increment nuse */
          xplasma->kbf_use[nuse] = n;
          nuse += 1;
        }
      }

    }
    xplasma->kbf_nuse = nuse;
  }


  return (0);
}

int sobolev_error_counter = 0;
/**********************************************************/
/**
 * @brief      calculates tau in the sobolev approxmation for a resonance, given the
 * conditions in the wind and the direction of the photon.
 *
 * It does not modify any of the variables that are passed to it, including for example
 * the photon.
 *
 * @param [in] WindPtr  one   A single wind cell
 * @param [in] double  x[]   A position
 * @param [in] double  den_ion   The density of the ion.  If less than 0, the routine calculates
 * the density at x
 * @param [in] struct lines *  lptr   A pointer to a particular ion
 * @param [in] double  dvds   the velocity gradient in the direction of travel of the photon
 * @return     The optical depth associated with a transition
 *
 * @details
 * In the sobolev approx  tau = PI e^2 / m_e c  * NL * f * lambda * 1 / (|dv/dx|) where dv/ds
 * is the velocity gradient in the direction the photon is travelling.  This can be made slightly
 * simpler by noting that lambda = freq / c and so  tau = PI e^w / m_e
 *
 * ### Notes ###
 *
 * The routine includes a correcton for the filling factor which
 * reduces tau
 *
 **********************************************************/
double
sobolev (one, x, den_ion, lptr, dvds)
     WindPtr one;               // This is a single cell in the wind
     double x[];
     double den_ion;
     struct lines *lptr;
     double dvds;
{
  double tau, xden_ion, tau_x_dvds, levden_upper;
  double two_level_atom (), d1, d2;
  int nion;
  double d_hold;
  int nplasma;
  int ndom;
  PlasmaPtr xplasma;

  nplasma = one->nplasma;
  xplasma = &plasmamain[nplasma];
  ndom = wmain[plasmamain->nwind].ndom;
  nion = lptr->nion;

  if ((dvds = fabs (dvds)) == 0.0)      // This forces dvds to be positive -- a good thing!
  {
    d1 = d2 = 0.;               // Eliminate warning when complied with clang
    tau = VERY_BIG;
    Error ("Sobolev: Surprise tau = VERY_BIG as dvds is 0\n");
  }

  else if (lptr->macro_info == 1 && geo.rt_mode == RT_MODE_MACRO && geo.macro_simple == 0)
  {
    // macro atom case SS
    d1 = den_config (xplasma, lptr->nconfigl);
    d2 = den_config (xplasma, lptr->nconfigu);
    levden_upper = xplasma->levden[config[lptr->nconfigu].nden];
  }

  else
  {
/* Next few steps to allow used of better calculation of density of this particular
ion which was done above in calculate ds.  It was made necessary by a change in the
calls to two_level atom
*/
    d_hold = xplasma->density[nion];    // Store the density of this ion in the cell

    if (den_ion < 0)
    {
      xplasma->density[nion] = get_ion_density (ndom, x, lptr->nion);   // Forced calculation of density
    }
    else
    {
      xplasma->density[nion] = den_ion; // Put den_ion into the density array
    }
    two_level_atom (lptr, xplasma, &d1, &d2);   // Calculate d1 & d2
    xplasma->density[nion] = d_hold;    // Restore w
    levden_upper = d2 / xplasma->density[nion];
  }

/* At this point d1 and d2 are known for all of the various ways sobolev can be called, and whether
 * macro atoms are involved or not
 */

  /* Check whether both d1 and d2 are below a minium value where we expect tau to be zero and where 
   * we can be subject to the effects of roundoff errors in terms of the determination of densities.
   * If densities are this low we expect the sobolev optical depth to be extremely small in any event
   * JM -- I've added something that checks if the fractional population for the upper level is below 
   * or equal to the minimum too.
   */

  if ((d1 < DENSITY_PHOT_MIN && d2 < DENSITY_PHOT_MIN) || (levden_upper <= DENSITY_MIN))
  {
    return (0);
  }


  xden_ion = (d1 - lptr->gl / lptr->gu * d2);

  if (xden_ion < 0)
  {
    sobolev_error_counter++;
    if (sobolev_error_counter < 100)
    {
      Error ("sobolev: VERY BAD population inversion in cell %d: d1 %g d2 %g g1 %g g2  %g freq %g f %g frac_upper %g\n",
             xplasma->nplasma, d1, d2, lptr->gl, lptr->gu, lptr->freq, lptr->f, levden_upper);
    }
    else if (sobolev_error_counter == 100)
    {
      Error ("sobolev: suppressing population inversion errors\n");
    }

    /* With the changes above to limit the densities the above error should not be happening, and if this does occur then 
     * we should determine why.  When we become convinced this problem has been dealt with effectively we can simplify this
     * code and just quit if the error occurs
     * ksl 181127
     */

    /* SS July 08: With macro atoms, the population solver can default to d2 = gu/gl * d1 which should
       give exactly zero here but can be negative, numerically.
       So I'm modyfying this to set tau to zero in such cases, when the populations are vanishingly small anyway. */
    tau_x_dvds = PI_E2_OVER_M * d1 * lptr->f / (lptr->freq);
    tau = tau_x_dvds / dvds;

    tau *= zdom[ndom].fill;     // filling factor is on a domain basis

    if (tau > 1.e-3)
    {
      /* JM -- I'm not sure why this particular value of tau is chosen, but I've added 
         an error message to give more information and make sure no ambiguity for code exit */
      Error ("sobolev: tau is >1e-3 and nu < gu/gl * nl. Exiting.\n");
      Exit (0);
    }

    else
    {
      //  Error ("sobolev: continuing by setting tau to zero\n");
      return (0.0);
    }
  }


  tau_x_dvds = PI_E2_OVER_M * xden_ion * lptr->f / (lptr->freq);
  tau = tau_x_dvds / dvds;

  /* JM 1411 -- multiply the optical depth by the filling factor */
  tau *= zdom[ndom].fill;

  return (tau);
}





/**********************************************************/
/**
 * @brief      determine a new direction and frequency for a photon
 * that is in the wind
 *
 * @param [in,out] PhotPtr  p   the  photon of interest
 * @param [in] int *  nres   either the number of the scatter
 * or a nonresonant scatter if nres < 0
 * @param [out] int *  nnscat   Returned from anisotropic thermal scattering model
 * @return  Always returns 0
 *
 * The results are stored in the PhotPtr p which contains the direction and frequency
 * of the scattered photon.
 *
 * @details
 * The routine calculates a new direction and frequency for a photon in both the
 * resonant and non-resonant cases.  In the frame of the wind, scattering is assumed
 * to be isotropic.
 *
 * ### Notes ###
 * This is the routine that is called when a resonant scatter does occur.  It is
 * relevant for both simple and macro atoms
 *
 * The equations for the frequency shifts are accurate only to first order in beta
 *
 * This routine should not move the photon at all, because other routines need to
 * take this photon in differing directions, and if one moves it here they may
 * encounter this resonance again
 *
 **********************************************************/

int
scatter (p, nres, nnscat)
     PhotPtr p;
     int *nres;
     int *nnscat;
{
  double z_prime[3];
  int which_out;
  struct photon p_orig;
  int i, n;
  double p_init[3], p_final[3], dp[3], dp_cyl[3];
  WindPtr one;
  double prob_kpkt, kpkt_choice, freq_comoving;
  double gamma_twiddle, gamma_twiddle_e, stim_fact;
  int m, llvl, ulvl;
  PlasmaPtr xplasma;
  MacroPtr mplasma;
  int ndom;


  /* get wind+plasma ptrs and domain number */
  one = &wmain[p->grid];
  xplasma = &plasmamain[one->nplasma];
  ndom = wmain[p->grid].ndom;

  /* Scatter assumes it has been passed a photon in the observer frame */

  if (check_frame (p, F_OBSERVER, "BeginScatter(err)\n") && modes.save_photons)
  {
    save_photons (p, "BeginScatter(err)");
  }

  stuff_phot (p, &p_orig);

  n = where_in_grid (ndom, p->x);       // Find out where we are

  if (n < 0)
  {
    Error ("scatter: Trying to scatter a photon in grid cell %d\n", n);
    return (-1);
  }

//OLD  if (modes.save_photons)
//OLD  {
//OLD    save_photons (p, "ScatIn");
//OLD  }

  if (observer_to_local_frame (p, p))
  {
    Error ("scatter: observer to local frame error (begin)\n");
  }
  freq_comoving = p->freq;



  /* So p is now in the local, or co-moving frame */


  /* On entering this subroutine we know that a photon packet has been
     absorbed. nres tells us which process absorbed it. There are currently
     four "absorption" processes: electron scattering (flagged -1) line
     absorption (flagged by +ve integer < NLINES), bf absorption
     (flagged by +ve integer > NLINES) and ff absorption (flagged -2). (SS) */

  /* BEGINNING OF SECTION FOR HANDLING MACRO-ATOMS 
     If the macro atom method is being used then the following section must be
     performed to select a macro atom deactivation process. If not then the
     deactivation process is always the same as the activation process and so
     nothing needs to be done. */

  if (geo.rt_mode == RT_MODE_MACRO)     //check if macro atom method in use
  {

    mplasma = &macromain[xplasma->nplasma];

    /* Electron scattering is the simplest to deal with. The co-moving
       frequency is unchanged so it's just a randomisation of the direction.
       For b-b and b-f processes it is first necessary to determine the
       process by which re-emission occurs. (SS). */

    if (*nres > (-1) && *nres < NLINES)
    {
      /* It's a bb line - we can go straight to macro_gov since we know that
         we don't want a k-packet immediately. macro_gov now makes the decision
         regarding the treament (simple or full macro). */
      /*XFRAME may need to check, but I believe all the macro atom (including macro gov) should already be 
         formulated in the cmf: i.e. provided that the photon was transformed to the cmf before this call 
         (see above) then nothing should be needed here before going in - need to consider reverse 
         transfer when coming back out though (likely in macro_gov?)
       */

      macro_gov (p, nres, 1, &which_out);
    }

    else if (*nres > NLINES)
    {
      /* It was a photoionisation process.
         For this case we need to decide first whether to excite
         a macro atom directly or to create a k-packet. */

      /*
         The probability of creating a k-packet is given by the
         mc estimators gamma, gamma_e, alpha_st, alpha_st_e.
         Start by identifying which estimators we want and then
         by computing gamma_twiddle (in Leon's notation -
         Lucy 2003 A&A 403 261 */

      /* Now, (Apr04) I'm adding the possibility that the continuum
         is not from a macro ion but from one that we don't have/want
         a macro atom treatment. If it's non macro atom all that happens
         is an on-the-spot decision about whether to excite a fake
         bf macro atom or create a k-packet. Since we've not recorded
         monte carlo estimators for simple ions the decision about creating
         a k-packet is based only on the frequency of the scattered packet.
         (SS) */

      if (phot_top[*nres - NLINES - 1].macro_info == 1 && geo.macro_simple == 0)
      {
        /* Macro ion case (SS) This is the case bound free interaction for a 
           xsection that is part of a macro atom (and we have not defaulted 
           to the simplifed approach) */

        /* Note:  NLINES-1 in the lines below is correct.  This is becasue
           the 1st bf is identified by nres = NLINES+1 and this is
           the zeroth element of phot_top: hence the -1.  SS
         */

        llvl = phot_top[*nres - NLINES - 1].nlev;       //lower level
        ulvl = phot_top[*nres - NLINES - 1].uplev;      //upper level

        for (m = 0; m < config[llvl].n_bfu_jump; m++)
        {
          if (config[llvl].bfu_jump[m] == *nres - NLINES - 1)
          {
            break;
          }
        }

        // m should now be the label to identify which of the bf processes from llvl
        // this is. Check that it is reasonable

        if (m > config[llvl].n_bfu_jump - 1)
        {
          Error ("scatter (resonate.c): could not identify bf transition. Abort. \n");
          Exit (0);
        }

        /* Need to compute the factor needed for the stimulated term. */

        stim_fact = den_config (xplasma, ulvl) / den_config (xplasma, llvl) / xplasma->ne;

        gamma_twiddle =
          mplasma->gamma_old[config[llvl].bfu_indx_first + m] - (mplasma->alpha_st_old[config[llvl].bfu_indx_first + m] * stim_fact);
        gamma_twiddle_e =
          mplasma->gamma_e_old[config[llvl].bfu_indx_first + m] - (mplasma->alpha_st_e_old[config[llvl].bfu_indx_first + m] * stim_fact);

        /* Both gamma_twiddles must be greater that zero if this is going to work. If they
           are zero then it's probably because this is the first iteration and so the've not
           been computed yet. For that first iteration k-packets will be ignored. If the
           gamma_twiddles are negative then something has gone wrong.
         */

        prob_kpkt = 0.0;        // initialise value
        if (gamma_twiddle > 0 && gamma_twiddle_e > 0)
        {
          prob_kpkt = 1. - gamma_twiddle / gamma_twiddle_e;
        }
        else if (gamma_twiddle == 0 && gamma_twiddle_e == 0)
        {
          prob_kpkt = 0.;
        }
        else
        {
          Error ("scatter (resonate.c): a gamma_twiddle is negative. Abort.\n");
          Exit (0);
        }

        /* Having got here we have calculated the probability of a k-packet
           being created. Now either make a k-packet or excite a macro atom. */

        kpkt_choice = random_number (0.0, 1.0); //random number for kpkt choice

        if (prob_kpkt > kpkt_choice)
        {
          macro_gov (p, nres, 2, &which_out);   //routine to deal with kpkt
        }
        else
        {
          macro_gov (p, nres, 1, &which_out);   //routine to deal with macro atom excitation
        }

        /* This ends the calculation for dealing with free-bound absorption of a macro atom */
      }
      else if (phot_top[*nres - NLINES - 1].macro_info == 0 || geo.macro_simple == 1)
      {
        /* Simple ion case.  It's a bf interaction in a calculation involving macro-atoms, 
           but this photoionization x-section is is associated with one of the levels
           a simple atom, not one that is part of a full blown macro atom.

           (Alternatively we are treating all atoms in a simplfied mode)

           Need to make decision about making a k-packet. Get the fraction of the energy
           that goes into the electron rather than being stored as ionisation energy: this
           fraction gives the selection probability for creating a k packet. It's given by 
           (1 - ionization edge frequecy / photon frequency)  */

        prob_kpkt = 1. - (phot_top[*nres - NLINES - 1].freq[0] / freq_comoving);

        if (prob_kpkt < 0)
        {
          /* only report an error for a negative prob_kpkt if it's large-ish in magnitude. see #436 discussion */
          if (prob_kpkt < -1e-3)
          {
            Error ("scatter: kpkt probability (%8.4e) < 0, zeroing\n", prob_kpkt);
            Log ("scatter: photon edge frequency: %8.4e, comoving frequency %8.4e\n", phot_top[*nres - NLINES - 1].freq[0], freq_comoving);
          }
          prob_kpkt = 0.0;
        }

        /* Now choose whether or not to make a k-packet. */

        kpkt_choice = random_number (0.0, 1.0); //random number for kpkt choice

        /* what happens next depends on whether we want to use the 
           altered mode for bound-free in "simple-macro mode". If we do,
           then the photon weight gets multiplied down by a factor (nu-nu_0)/nu
           and we force a kpkt to be created */

#if BF_SIMPLE_EMISSIVITY_APPROACH
        /* This is the new approach which does not explicityly conserve energy.
           Re record the amount of energy going into the simple ion ionization pool.  This
           version does not produce nan if prob_kpt is 0 and then reduce the photon weight
           to allow for the portion of the energy that went into the ionization pool before
           generating a kpkt.  In this approach we always generate a kpkt */

        xplasma->bf_simple_ionpool_in += p->w * (1 - prob_kpkt);
        p->w *= prob_kpkt;

        macro_gov (p, nres, 2, &which_out);     //routine to deal with kpkt
#else
        /* This is the old apporach.  Process the BF photon for a simple atom.  In this
           approach we generate a kpkt or an r-packet depending on whether the probility
           of creating a kpkt, namely prob_kpkt
         */

        if (prob_kpkt > kpkt_choice)
        {
          macro_gov (p, nres, 2, &which_out);   //routine to deal with kpkt
        }
        else
        {
          macro_gov (p, nres, 1, &which_out);   //routine to deal with fake macro atom bf excitation
        }
#endif

      }
      else
      {
        /* Our best-laid schemes have gang agley. It should never get here unless the input has been
           messed up in some way. (SS) */
        Error ("scatter (resonate.c): continuum scatter - seems to be neither macro nor simple. Abort.\n");
        Exit (0);
      }
    }
    else if (*nres == -2)
    {                           /* This is a ff event (SS). */
      macro_gov (p, nres, 2, &which_out);       //ff always make a k-packet
    }



  }
  /*XFRAME need to make sure that macro_gov leaves the photon in the Local frame. We will do all
     of the transformations back to the observer frame here.
   */


  /* END OF SECTION FOR HANDLING ASPECTS OF SCATTERING PROCESSES THAT ARE SPECIFIC TO MACRO-ATOMS. */

  /* Set nres  correctly and make sure the frequency is correct
     for a resonanant scatter. Note that nres may have changed especially for macro-atoms */

  p->nres = *nres;
  if (*nres > -1 && *nres < nlines)
  {
    p->freq = lin_ptr[*nres]->freq;
  }


  /* Now determine the direction of the scattered photon, for electrons scatering (-1), ff emision (-2), or
     bound free emission (>NLINES), allowing depending on the scattering mode for thermal trapping. 
     Note that this portion of the code is identical for both simple and macro atoms, except for the fact
     that ff and bf are only treated as scattering processin in macro-atom mode.
   */

  if (*nres == -1)
  {

    compton_dir (p);            // uses the KN formula

  }
  else if (*nres == -2 || *nres > NLINES || geo.scatter_mode == SCATTER_MODE_ISOTROPIC)
  {
    randvec (z_prime, 1.0);
    stuff_v (z_prime, p->lmn);
  }
  else
  {
    randwind_thermal_trapping (p, nnscat);      // the thermal trapping case
  }



  /* Finally put everything back in the observer frame */

  local_to_observer_frame (p, p);


  /* If we are in macro-atom mode, add the photon to the created wind spectrum.  For simple
     atoms the wind spectrum is constructed in sectra.c */

  if (geo.rt_mode == RT_MODE_MACRO && *nres != -1)
  {
    p->nmacro++;
    spec_add_one (p, SPEC_CWIND);
  }

/*Now calculate the momentum transfer.  What follows appears to be
correct only if there was no energy absorbed at the scattering site.
The rest of this is only needed during ionization cycles, before the wind itself
if fixed.  
*/

//OLD  if (modes.save_photons)
//OLD  {
//OLD    save_photons (p, "ExScatter");
//OLD  }

//XFRAME - Is this  the correct way to calculate the momentum transfer allowing for CMF calculation ? 
// ioniz_or_extract is True for ionization cycles, false for spectral cycles

  if (geo.ioniz_or_extract)
  {
    stuff_v (p_orig.lmn, p_init);
    renorm (p_init, p_orig.w / VLIGHT);
    stuff_v (p->lmn, p_final);
    renorm (p_final, p->w / VLIGHT);
    vsub (p_final, p_init, dp);

    project_from_xyz_cyl (p_orig.x, dp, dp_cyl);



    if (p_orig.x[2] < 0)
      dp_cyl[2] *= (-1);
    for (i = 0; i < 3; i++)
    {
      xplasma->dmo_dt[i] += dp_cyl[i];
    }

  }


  return (0);
}
