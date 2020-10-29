
/***********************************************************/
/** @file  trans_phot.c
 * @author ksl
 * @date   April, 2018
 *
 * @brief  This file contains high level routines which carry
 * out the propagation of photons through the wind/wind domains
 * both for ionization and spectrum (extract) cycles
 *
 * ### Notes ###
 *
 * The routines contained here are central to Python, and anyone
 * who wants to understand Python in general should spend time
 * understanding how they work
 *
 * There are two basic options associated with the routines.
 *
 * * live or die is the case in which one is only concerned
 * with how a photon goes through the wind and random numbers
 * are used to generate the scattering directions.  The photon is
 * followed until it leaves the system.  This is the only option
 * used in ionization cycles (currently).
 *
 * * extract is the mode where whenever a photon scatters we
 * also calculate what the photon weight would be if it scattered
 * in a set of certain directions.  When this happens the weight
 * of the photon is changed to reflect the fact that the direction
 * is not random.
 *
 * The extract option is used
 * normally during the spectral extraction cycles.
 * However, as an advanced option one can use the live or die
 * to construct the detailed spectrum.  One would not normally
 * want to do this, as many photons are "wasted" since they
 * don't scatter at the desired angle.  With sufficient numbers
 * of photons however the results of the two methods should
 * (by construction) be identical (or at least very very
 * similar).
 *
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "atomic.h"
#include "python.h"

long n_lost_to_dfudge = 0;


/**********************************************************/
/**
 * @brief      int (w,p,iextract) oversees the propagation of a "flight" of photons
 *
 * @param [in] WindPtr  w   The entire wind domain
 * @param [in, out] PhotPtr  p   A pointer to a "fligh" of photons
 * @param [in] int  iextract   An integer controlling whether we are to process the
 * flight in the live or die option (0) or whether we also need to extract photons in
 * specific directions (which is usually the case in constructing spectra
 * @return   Normally returns 0
 *
 * @details
 * This routine oversees the propagation of  individual photons.  The main loop
 * covers an entire "flight" of photons.   The routine generates the random
 * optical depth a photon can travel before scattering and monitors the
 * progress of the photon through the grid.
 *
 * The real physics is done elsewhere, in lower level routines.
 *
 * ### Notes ###
 *
 * At the end of the routine the position for each of the photons in p is the
 * last point where the photon was in the wind, * not the outer boundary of
 * the radiative transfer
 *
 **********************************************************/

int
trans_phot (WindPtr w, PhotPtr p, int iextract)
{
  int nphot;
  struct photon pp, pextract;
  int absorb_reflect;           /* this is a variable used to store geo.absorb_reflect during exxtract */
  int nreport;
  struct timeval timer_t0;
  double rho;

  nreport = 100000;
  if (nreport < NPHOT / 100)
  {
    nreport = NPHOT / 100;
  }


  Log ("\n");

  timer_t0 = init_timer_t0 ();

  /* Beginning of loop over photons */

  for (nphot = 0; nphot < NPHOT; nphot++)
  {

    check_frame (&p[nphot], F_OBSERVER, "trans_phot_start\n");
    if (modes.save_photons)
    {
      save_photons (&p[nphot], "trans_phot_start");
      rho = sqrt (p[nphot].x[0] * p[nphot].x[0] + p[nphot].x[1] * p[nphot].x[1]);
      if (fabs (p[nphot].x[2]) <= zdisk (rho))
      {
        Diag ("ZDISK %d  %e < = %e delta= %e\n", nphot, fabs (p[nphot].x[2]), rho, rho - fabs (p[nphot].x[2]));
      }
    }



    /* This is just a watchdog method to tell the user the program is still running */
    if (nphot % nreport == 0)
    {
      if (geo.ioniz_or_extract == CYCLE_IONIZ)
        Log (" Ion. Cycle %d/%d of %s : Photon %10d of %10d or %6.1f per cent \n", geo.wcycle + 1, geo.wcycles, basename, nphot, NPHOT,
             nphot * 100. / NPHOT);
      else
        Log ("Spec. Cycle %d/%d of %s : Photon %10d of %10d or %6.1f per cent \n", geo.pcycle + 1, geo.pcycles, basename, nphot, NPHOT,
             nphot * 100. / NPHOT);
    }
    Log_flush ();

    stuff_phot (&p[nphot], &pp);
    absorb_reflect = geo.absorb_reflect;


    /* The next if statement is executed if we are calculating the detailed spectrum and makes sure we always run extract on
       the original photon no matter where it was generated */
    if (iextract)
    {

      stuff_phot (&p[nphot], &pextract);
      extract (w, &pextract, pextract.origin);

    }

    p[nphot].np = nphot;

    /* Transport a single photon */
    trans_phot_single (w, &p[nphot], iextract);

  }

  /* This is the end of the loop over all of the photons; after this the routine returns */

  /* Line to complete watchdog timer */
  Log ("\n");

  print_timer_duration ("!!python: photon transport completed in", timer_t0);
  /* sometimes photons scatter near the edge of the wind and get pushed out by DFUDGE. We record these */
  if (n_lost_to_dfudge > 0)
    Error
      ("trans_phot: %ld photons were lost due to DFUDGE (%8.4e) pushing them outside of the wind after scatter\n",
       n_lost_to_dfudge, DFUDGE);

  n_lost_to_dfudge = 0;         // reset the counter

  return (0);
}






/**********************************************************/
/**
 * @brief      Transport a single photon photon through the wind.
 *
 * @param [in] WindPtr  w   The entire wind
 * @param [in, out] PhotPtr  p   A single photon
 * @param [in] int  iextract   If 0, then process this photon in the live or die option, without
 * calling extract
 *
 * @return     Always returns 0
 *
 * @details
 * This routine oversees the propagation of an individual photon through
 * the wind.  As the photon moves through the wind, its position, direction
 * and weight are updated.  In reverberation mode, the flight time is also
 * tracked.
 *
 * Basically what the routine does is generate a random number which is used to
 * determine the optical depth to a scatter, and then it calles translate
 * multiple times.   translate involves moving the photon only a single cell
 * (or alternatively a single tranfer in the windless region), and returns
 * a status.  Depending on what this status is, trans_phot_single calls
 * trans_phot again doing nothing, but if the scattering depth has been
 * reached, then trans_phot_single causes the photon to scatter,
 * which changes its direction.  This process continues until the photon
 * exits the system or hits a barrier.  
 
 * If the photon hits a radiating
 * surface, the disk or star, then the photon may either be absorbed,
 * or scattered depending on the reflection/absorption mode. If the
 * reflection/absorption mode is set to reflection, then the 
 * program is redirected in the main loop, but if the surfaces are
 * set to aborb, hitting a surface will exist the routine.
 *
 *
 *
 * ### Notes ###
 *
 * This routine is called by trans_phot once for each photon in a flight of photons
 * Internally, there are two main PhotPtrs p, and pp.  pp is a place that
 * the photon has reached, and p is the location where it is going.  At
 * the end of the main loop before a new cycle, pp is updated.
 *
 *
 **********************************************************/

int
trans_phot_single (WindPtr w, PhotPtr p, int iextract)
{
  double tau_scat, tau;
  int istat;
  int ierr;
  int icell;
  int current_nres;
  int kkk, n;
  double weight_min;
  struct photon pp, pextract;
  int nnscat;
  double x_dfudge_check[3];
  int ndom;
  double normal[3];
  double rho, dz;


  check_frame (p, F_OBSERVER, "trans_phot_single: Starting Error\n");

  /* Initialize parameters that are needed for the flight of the photon through the wind */
  stuff_phot (p, &pp);
  tau_scat = -log (1. - random_number (0.0, 1.0));

  weight_min = EPSILON * pp.w;
  istat = P_INWIND;
  tau = 0;
  icell = 0;
  n = 0;


  /* Next section is really to allow one to follow a specific photon for diagnositic purposes
     This may be useful more generally but it should at least be commented out one diagnositcs
     are complete */
//OLD if (geo.wcycle == 1 && p->np == 5474)
//OLD  {
//OLD     Log ("Hello: Starting photon  %d in cycle %d we want to follow\n", p->np, geo.wcycle);
//OLD   }


  if (modes.save_photons)
  {
    save_photons (p, "trans_phot_single:Begin");
  }
  /* This is the beginning of the loop for a single photon and executes until the photon leaves the wind */

  while (istat == P_INWIND)
  {

    /* The call to translate below involves only a single cell (or alternatively a single transfer 
       in the windless region). 

       istat as returned by should either be P_INWIND in which case the photon hit the other side 
       of the cell without scattering or P_SCAT in which case there was a scattering event in the shell, 
       P_ESCAPE in which case the photon reached the outside edge of the grid and escaped, P_STAR in
       which case it reach the inner central object, etc. If the photon escapes then we leave the 
       photon at the position of it's last scatter.  In most other cases though we store the final 
       position of the photon. */


    if (modes.save_photons)
    {
      save_photons (&pp, "BeforeTranslate");
    }

    istat = translate (w, &pp, tau_scat, &tau, &current_nres);

    if (modes.save_photons)
    {
      save_photons (&pp, "AfterTranslate");
    }


    icell++;
    istat = walls (&pp, p, normal);
    /* pp is where the photon is going, p is where it was  */

    if (istat == P_ERROR)
    {
      Error_silent ("trans_phot: Abnormal return from translate on photon %d\n", p->np);
      break;
    }

    if (pp.w < weight_min)
    {
      istat = pp.istat = P_ABSORB;      /* This photon was absorbed by continuum opacity within the wind */
      pp.tau = VERY_BIG;
      stuff_phot (&pp, p);
      break;
    }

    if (istat == P_HIT_STAR)
    {                           /* It hit the star */
      geo.lum_star_back += pp.w;
      spec_add_one (&pp, SPEC_HITSURF);
      if (geo.absorb_reflect == BACK_RAD_SCATTER)
      {
        /* If we got here, the a new photon direction needs to be defined that will cause the photon
         * to continue in the wind.  Since this is effectively a scattering event we also have to
         * extract a photon to construct the detailed spectrum
         */


        randvcos (pp.lmn, normal);
        move_phot (&pp, DFUDGE);
        stuff_phot (&pp, p);
        p->ds = 0;
        tau_scat = -log (1. - random_number (0.0, 1.0));
        istat = pp.istat = P_INWIND;    /* Set the status back to P_INWIND so the photon will continue */
        tau = 0;
        if (iextract)
        {
          stuff_phot (&pp, &pextract);
          extract (w, &pextract, PTYPE_STAR);   // Treat as stellar photon for purpose of extraction
        }
      }
      else
      {                         /*In this case, photons that hit the star are simply absorbed 
                                   so this is the end of the line. */
        stuff_phot (&pp, p);
        break;
      }
    }

    if (istat == P_HIT_DISK)
    {
      /* It hit the disk */
      /* ZFRAME - this next section assumes that disk heating is supposed to be carried out
         in the local frame of the disk.  That this is the correct thing to do needs to
         be confirmed.
       */



      /* Store the energy of the photon bundle into a disk structure so that one 
         can determine later how much and where the disk was heated by photons.
         Note that the disk is defined from 0 to NRINGS-2. NRINGS-1 contains the position 
         of the outer radius of the disk. */


      rho = sqrt (pp.x[0] * pp.x[0] + pp.x[1] * pp.x[1]);
      kkk = 0;
      while (rho > qdisk.r[kkk] && kkk < NRINGS - 1)
        kkk++;
      kkk--;                    /* So that the heating refers to the heating between kkk and kkk+1 */
      qdisk.nhit[kkk]++;
      geo.lum_disk_back = qdisk.heat[kkk] += pp.w;
      qdisk.ave_freq[kkk] += pp.w * pp.freq;

      if (modes.save_photons && rho > geo.diskrad)
      {
        Diag ("trans_phot: Photon %d hit disk at %.3e \n", pp.np, rho);
        save_photons (p, "BeforeHitDisk");
        save_photons (&pp, "AfterHitDisk");

      }

      if (geo.absorb_reflect == BACK_RAD_SCATTER)
      {

        if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
        {
          dz = (zdisk (rho) - fabs (pp.x[2]));
          if (dz > 0)
          {
            if (modes.save_photons)
            {
              Diag ("trans_phot: Photon %d is still in the disk, zdisk %e diff %e\n", pp.np, zdisk (rho), dz);
            }
            if (pp.x[2] > 0)
            {
              pp.x[2] += (dz + 1000.);
            }

            else
            {
              pp.x[2] -= (dz + 1000.);
            }

            if (zdisk (rho) > fabs (pp.x[2]))
            {
              if (modes.save_photons)
              {
                Diag ("trans_phot: Photon %d is still in the disk, zdisk %e diff %e\n", pp.np, zdisk (rho), zdisk (rho) - fabs (pp.x[2]));
              }
            }
          }
        }


        if (modes.save_photons)
        {
          save_photons (&pp, "HitDisk");
        }


        spec_add_one (&pp, SPEC_HITSURF);



        /* If we got here, the a new photon direction needs to be defined that will cause the photon
         * to continue in the wind.  Since this is effectively a scattering event we also have to
         * extract a photon to construct the detailed spectrum
         */
        randvcos (pp.lmn, normal);

        stuff_phot (&pp, p);

        p->ds = 0;
        tau_scat = -log (1. - random_number (0.0, 1.0));
        istat = pp.istat = P_INWIND;
        tau = 0;
        stuff_phot (&pp, p);
        if (iextract)
        {
          stuff_phot (&pp, &pextract);
          extract (w, &pextract, PTYPE_DISK);
        }
      }
      else
      {                         /* In this case, photons that hit the disk are to be absorbed 
                                   so this is the end of the line. */
        stuff_phot (&pp, p);
        break;
      }

      if (modes.save_photons)
      {
        save_photons (&pp, "HitDisk");
      }
    }

    if (istat == P_SCAT)
    {                           /* Cause the photon to scatter and reinitialize */


      pp.grid = n = where_in_grid (wmain[pp.grid].ndom, pp.x);

      if (n < 0)
      {
        Error ("trans_phot: Trying to scatter a photon which is not in the wind\n");
        Error ("trans_phot: %d grid %3d x %8.2e %8.2e %8.2e (%8.2e)\n", pp.np, pp.grid, pp.x[0], pp.x[1], pp.x[2],
               sqrt (pp.x[0] * pp.x[0] + pp.x[1] * pp.x[1] + pp.x[2] * pp.x[2]));
        Error ("trans_phot: This photon is effectively lost!\n");
        istat = pp.istat = p->istat = P_ERROR;
        stuff_phot (&pp, p);
        break;
      }

      /* 1506 JM -- If the next error reoccurs, see Issue #154 for discussion */

      if (wmain[n].nplasma == NPLASMA)
      {
        Error ("trans_phot: Trying to scatter a photon which is not in a cell in the plasma structure\n");
        Error ("trans_phot: %d grid %3d x %8.2e %8.2e %8.2e\n", pp.np, pp.grid, pp.x[0], pp.x[1], pp.x[2]);
        Error ("trans_phot: This photon is effectively lost!\n");
        istat = pp.istat = p->istat = P_ERROR;
        stuff_phot (&pp, p);
        break;
      }


      if (wmain[n].inwind < 0)
      {
        Error ("trans_phot: Trying to scatter a photon in a cell with no wind volume\n");
        Error ("trans_phot: %d grid %3d x %8.2e %8.2e %8.2e\n", pp.np, pp.grid, pp.x[0], pp.x[1], pp.x[2]);
        Error ("trans_phot: istat %d\n", pp.istat);
        Error ("trans_phot: This photon is effectively lost!\n");
        istat = pp.istat = p->istat = P_ERROR;
        stuff_phot (&pp, p);
        break;

      }

      /* Add path lengths for reverberation mapping */
      if ((geo.reverb == REV_WIND || geo.reverb == REV_MATOM) && geo.ioniz_or_extract == CYCLE_IONIZ && geo.wcycle == geo.wcycles - 1)
      {
        wind_paths_add_phot (&wmain[n], &pp);
      }


      /* SS July 04 - next lines modified so that the "thermal trapping" model of anisotropic scattering is included in the
         macro atom method. What happens now is all in scatter - within that routine the "thermal trapping" model is used to
         decide what the direction of emission is before returning here.
       */

      nnscat = 1;
      pp.nscat++;


      if (modes.save_photons)
      {
        save_photons (&pp, "BeforeScat");
      }

      ierr = scatter (&pp, &current_nres, &nnscat);

      if (modes.save_photons)
      {
        save_photons (&pp, "AfterScat");
      }
      if (ierr)
      {
        Error ("trans_phot: bad return from scatter %d at point 2\n", ierr);
      }


      /* SS June 04: During the spectrum calculation cycles, photons are thrown away when they interact with macro atoms or
         become k-packets. This is done by setting their weight to zero (effectively they have been absorbed into either
         excitation energy or thermal energy). Since they now have no weight there is no need to follow them further. */

      if (geo.matom_radiation == 1 && geo.rt_mode == RT_MODE_MACRO && pp.w < weight_min)
        /* Flag for the spectrum calculations in a macro atom calculation SS */
      {
        istat = pp.istat = P_ABSORB;
        pp.tau = VERY_BIG;
        stuff_phot (&pp, p);
        break;
      }

      // Calculate the line heating and if the photon was absorbed break finish up
      // XXXX ??? Need to modify line_heat for multiple scattering but not yet
      // Condition that nres < nlines added (SS)

      if (current_nres > -1 && current_nres < nlines)
      {
        pp.nrscat++;

        /* This next statement writes out the position of every resonant scattering event to a file */
        if (modes.track_resonant_scatters)
          track_scatters (&pp, wmain[n].nplasma, "Resonant");


        plasmamain[wmain[n].nplasma].scatters[line[current_nres].nion] += 1;

        if (geo.rt_mode == RT_MODE_2LEVEL)      // only do next line for non-macro atom case
        {
          line_heat (&plasmamain[wmain[n].nplasma], &pp, current_nres);
        }

        if (pp.w < weight_min)
        {
          istat = pp.istat = P_ABSORB;  /* This photon was absorbed by continuum opacity within the wind */
          pp.tau = VERY_BIG;
          stuff_phot (&pp, p);
          break;
        }
      }



      /* Now extract photons if we are in detailed the detailed spectrum portion of the program */

      /* N.B. To use the anisotropic scattering option, extract needs to follow scatter.  
       * This is because the reweighting which occurs in extract needs the pdf for scattering 
       * to have been initialized.
       */

      if (iextract)
      {
        stuff_phot (&pp, &pextract);
        pextract.nnscat = nnscat;
        extract (w, &pextract, PTYPE_WIND);     // Treat as wind photon for purpose of extraction
      }


      /* OK we have completed extract, if that had to be done, 
       * need to reinitialize parameters for the scattered photon so it can 
       * can continue throug the wind 
       */

      tau_scat = -log (1. - random_number (0.0, 1.0));
      istat = pp.istat = P_INWIND;
      tau = 0;

      stuff_v (pp.x, x_dfudge_check);   // this is a vector we use to see if dfudge moved the photon outside the wind cone

      pp.ds = 0;
      stuff_phot (&pp, p);

      reposition (&pp);

      /* JM 1506 -- call walls again to account for instance where DFUDGE
         can take photon outside of the wind and into the disk or star
         after scattering. Note that walls updates the istat in pp as well.
         This may not be necessary but I think to account for every eventuality
         it should be done. This *does not* update istat if the photon scatters
         outside of the wind- I guess P_IN_WIND is really in wind or empty space
         but not escaped. translate_in_space will take care of this next time
         round. All a bit convoluted but should work. */

      istat = walls (&pp, p, normal);

      if (istat != p->istat)
      {
        Error ("trans_phot:Status of %9d changed from %d to %d after reposition\n", p->np, p->istat, istat);
      }


      /*ksl - eliminated reposition_lost_photon from code, as this should not happen anymore.  If it
         does then it needs to be investigated. */

      if (istat == P_REPOSITION_ERROR)
      {
        Error ("Got reposition error for %d.  INVESTIGATE\n", p->np);
      }

      /* JM 1506 -- we don't throw errors here now, but we do keep a track
         of how many 4 photons were lost due to DFUDGE pushing them
         outside of the wind after scatter */

      if (where_in_wind (pp.x, &ndom) != W_ALL_INWIND && where_in_wind (x_dfudge_check, &ndom) == W_ALL_INWIND)
      {
        n_lost_to_dfudge++;     // increment the counter (checked at end of trans_phot)
      }

      stuff_phot (&pp, p);
      icell = 0;
    }

    /* This completes the portion of the code that handles the scattering of a photon.
     * What follows is a simple check to see if
     * this particular photon has gotten stuck in the wind */

    if (pp.nscat == MAXSCAT)
    {
      istat = pp.istat = P_TOO_MANY_SCATTERS;
      stuff_phot (&pp, p);
      break;
    }

    if (pp.istat == P_ERROR_MATOM || pp.istat == P_LOFREQ_FF || pp.istat == P_ADIABATIC)
    {
      istat = p->istat = pp.istat;
      stuff_phot (&pp, p);
      break;
    }

    /* This appears partly to be an insurance policy. It is not obvious that for example nscat
     * and nrscat need to be updated */

    p->istat = istat;
    p->nscat = pp.nscat;
    p->nrscat = pp.nrscat;
    p->w = pp.w;                // Assure that final weight of photon is returned.
  }
  /* This is the end of the loop over a photon */

  /* The next section is for diagnostic purposes.  There are two possibilities.  If you wish to know where
   * the photon was last while in the wind, you want to track p; if you wish to know where it hits the
   * outer boundary of the calculation you would want pp.  So one should keep both lines below, and comment
   * out the one you do not want. */

  if (modes.save_photons)
  {
    save_photons (&pp, "Final");        //The position of the photon where it exits the calculation
  }

  return (0);
}
