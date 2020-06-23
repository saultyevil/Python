
/***********************************************************/
/** @file  extract.c
 * @author ksl
 * @date   April, 2018
 *
 * @brief  Maing routines for extracting photons during the 
 * spectral generation phase.
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/**********************************************************/
/** 
 * @brief      A supervisory routine called to 
 * 	builds detailed spectra in the normal (extract) mode.
 *
 * @param [in] WindPtr  w   The entire wind
 * @param [in] PhotPtr  p   The photon to extract
 * @param [in] int  itype   An integer representing the type of photon 
 * for the purpose of being extracted.
 * @return     Always returns 0
 *
 * @details
 * extract is called when a photon begins its flight and every time that photon
 * scatters, unless the user has exercised the "live or die" option, in
 * which case it is not called.  
 *
 * extract carries out several preperatory steps for the extraction
 * and for each spectrum one wants to build, it calls extract_one,
 * where the actual incrementing of the spectrum is done.
 *
 * itype takes on the following values:
 * * PTYPE_STAR->the photon came for the star 
 * * PTYPE_BL->the photon came from the boundary layer
 * * PTYPE_DISK->the photon being redirected arose in the disk,
 * * PTYPE_WIND->the photon being redirected arose in the wind,
 * 
 * extract uses the types to prepare the photon for extraction, including
 * doppler shifting the photon if is of PTYPE_WIND or PTYPE_DISK.
 * 
 * Usually, Python constructs a spectrum of all photons, but there are
 * advanced options which allone to restict the spectrum created to
 * those produced with a certain number of scatters or from photons 
 * that arise from the above or below the disk.  extract enforces
 * those choices before calling extract_one.
 * The parameters for this option all come in through
 * python.h, and are contained in the spectrum structure. 
 *
 * The basic choices, some of which can be used in tandom are as follows:
 *
 * *	If s[n].nscat>999; then one obtains the entire spectrum
 * *	if s[n].nscat is a positive number, then the spectrum is composed just 
 * 		of photons with that number of scatters
 * *	if s[n].nscat is a negative number, then the spectrum contains all photons
 * 		with >=  the absolute value of s[n].nscat
 * * 	s[n].top_bot=0 -> accept all photons
 * * 	s[n].top_bot    >0  -> accept photons last above the disk
 * * 	s[n].top_bot    <0  -> accept photons below the disk
 * *    One can also select photons that are to be extracted from a particular
 * spherical region.
 * ### Notes ###
 * 
 *
 * @bug This is also commented in the text, but there is a rather bizarre separation
 * for where the photon frequency is updated and where the weight is update. The former is
 * done in extract, while the latter is done in extract_one.  This is not an error precisely
 * but makes the code more confusing than it needs to be.
 *
 **********************************************************/

int
extract (w, p, itype)
     WindPtr w;
     PhotPtr p;
     int itype;
{
  int n, mscat, mtopbot;
  struct photon pp, p_in;
  int yep;
  double xdiff[3];
  double p_norm, tau_norm;



  /* Make sure the input photon is not modified */

  stuff_phot (p, &p_in);


/* The next section was moved from trans_phot 200518 */

  /* We increase weight to account for the number of scatters. This is done because in extract we multiply by the escape
     probability along a given direction, but we also need to divide the weight by the mean escape probability, which is
     equal to 1/nnscat.  See issue #710 for a more extended explanation of how the weight is renormalized stocahstically. */




  if (itype == PTYPE_WIND)
  {
    if (geo.scatter_mode == SCATTER_MODE_THERMAL && p_in.nres <= NLINES && p_in.nres > -1)
    {
      /* we normalised our rejection method by the escape probability along the vector of maximum velocity gradient.
         First find the sobolev optical depth along that vector. The -1 enforces calculation of the ion density */

      tau_norm = sobolev (&wmain[p_in.grid], p_in.x, -1.0, lin_ptr[p_in.nres], wmain[p_in.grid].dvds_max);

      /* then turn into a probability */
      p_norm = p_escape_from_tau (tau_norm);

    }
    else
    {
      p_norm = 1.0;

      /* throw an error if nnscat does not equal 1 */
      if (p_in.nnscat != 1)
        Error
          ("trans_phot: nnscat is %i for photon %i in scatter mode %i! nres %i NLINES %i\n",
           p_in.nnscat, p_in.np, geo.scatter_mode, p_in.nres, NLINES);
    }

    p_in.w *= p_in.nnscat / p_norm;

  }


  if (itype == PTYPE_WIND)
  {
    observer_to_local_frame (&p_in, &p_in);
  }
  if (itype == PTYPE_DISK)
  {
    observer_to_local_frame_disk (&p_in, &p_in);
  }



  for (n = MSPEC; n < nspectra; n++)
  {
    /* If statement allows one to choose whether to construct the spectrum
       from all photons or just from photons that have scattered a specific number
       of times or in specific regions of the wind. A region is specified by a position
       and a radius. */

    yep = 1;                    // Start by assuming it is a good photon for extraction

    if ((mscat = xxspec[n].nscat) > 999 || p_in.nscat == mscat || (mscat < 0 && p_in.nscat >= (-mscat)))
      yep = 1;
    else
      yep = 0;

    if (yep)
    {
      if ((mtopbot = xxspec[n].top_bot) == 0)
        yep = 1;                // Then there are no positional parameters and we are done
      else if (mtopbot == -1 && p_in.x[2] < 0)
        yep = 1;
      else if (mtopbot == 1 && p_in.x[2] > 0)
        yep = 1;
      else if (mtopbot == 2)    // Then to count, the photom must originate within sn.r of sn.x
      {
        vsub (p_in.x, xxspec[n].x, xdiff);
        if (length (xdiff) > xxspec[n].r)
          yep = 0;

      }
      else
        yep = 0;
    }

    if (yep)                    //Then we want to extract this photon
    {


/* Create a photon pp to use here and in extract_one.  This assures we
 * have not modified p as part of extract.
 *
 * If it is a wind photon, it will have be in the observer frame in 
 * a different direction so we need to put it into the local frame
 * This needs to be done before we stuff the new direction in
 */

      stuff_phot (&p_in, &pp);
      stuff_v (xxspec[n].lmn, pp.lmn);  /* Stuff new photon direction into pp */

/* 

Need to frequency shift the disk photons as well as the wind 
photons.    

 */

      if (itype == PTYPE_DISK)
      {
//OLD        pp.freq = pp.freq_orig;
//OLD        pp.frame = F_LOCAL;
        local_to_observer_frame_disk (&pp, &pp);

      }
      if (itype == PTYPE_WIND)
      {                         /* If the photon was scattered in the wind, 
                                   the frequency also must be shifted */

/* XFRAME  Doppler shift the photon  to new direction.  In what follows
   we make the assumption which seems explicit in the old doppler routine that we 
   are in the observe frame
 */
//OLD Lines below look like belt and suspenders, but we should already be in the local
//OLD and so we should trust.
//OLD       if (pp.nres > -1 && pp.nres < nlines)
//OLD        {
//OLD          pp.freq = lin_ptr[pp.nres]->freq;
//OLD        }
        local_to_observer_frame (&pp, &pp);

      }

      if (modes.save_extract_photons && 1545.0 < 2.997925e18 / pp.freq && 2.997925e18 / pp.freq < 1565.0)
      {
        save_extract_photons (n, p, &pp);
      }


      /* Now extract the photon */
//OLD      if (modes.save_photons)
//OLD      {
//OLD        Diag ("BeforeExtract freq  %10.3e itype %d  nres %d\n", pp.freq, itype, pp.nres);
//OLD        save_photons (&pp, "BeforeExtract");
//OLD      }


      extract_one (w, &pp, itype, n);

//OLD      if (modes.save_photons)
//OLD      {
//OLD        save_photons (&pp, "AfterExtract");
//OLD      }


    }

  }
  return (0);
}





/**********************************************************/
/** 
 * @brief      Extract a single photon along a single line of sight.
 *
 * @param [in] WindPtr  w   The entire wind
 * @param [in] PhotPtr  pp  The photon to be extracted
 * @param [in] int  itype   The type of photon (star, disk, wind, etc)
 * @param [in] int  nspec   the spectrum which will be incremented
 * @return     The photon status after translation
 *
 * @details
 * extract_one is analogous to the detailed portion of transphot except here the
 * basic point is to calculate the optical depth through the plasma in a certain
 * direction, and to increment the appropriate spectrum.  
 *
 * Unlike trans_phot, this routine also checks to see if 
 * whether the photon hits the secondary star, if one exists.
 *
 * ### Notes ###
 * The logic behind the weighting of the photons is described in Christian Knigge's thesis in
 * section 2.3.4.  According to equation 2.19
 * 	Pc/Pw=12 cos(theta)*(1+b cos(theta)/(3+2b) where b=1.5 corresponds to the
 * Eddington approximation.
 * 
 * In Python, and in extract and transphot in particular, tau generally refers to the tau associated
 * with scattering processes, and the weight contains the effect of dimunition of the energy of
 * the photon bundle due to pure absorption processes.  So, in extract, we add pp->w * exp(-tau)
 * to the spectrum.
 *
 * Note that both linearly and logarithmically spaced spectra are produced.
 *
 **********************************************************/

int
extract_one (w, pp, itype, nspec)
     WindPtr w;
     PhotPtr pp;
     int itype, nspec;

{
  int istat, nres;
  struct photon pstart;
  double weight_min;
  int icell;
  int k, k1;
  double x[3];
  double tau;
  double zz;
  double dvds;
  double lfreqmin, lfreqmax, ldfreq;
  int ishell;
  double normal[3];


  weight_min = EPSILON * pp->w;
  istat = P_INWIND;
  tau = 0;
  icell = 0;

/* Preserve the starting position of the photon so one can use this to determine whether the
 * photon encountered the disk or star as it tried to exist the wind.
 */

  stuff_phot (pp, &pstart);

/* Reweight the photons. Note that photons have already been frequency shifted prior 
to entering extract */

  if (itype == PTYPE_STAR || itype == PTYPE_BL)
  {                             /* It was an unscattered photon from the star */
    stuff_v (pp->x, x);
    renorm (x, 1.);
    zz = fabs (dot (x, xxspec[nspec].lmn));
    pp->w *= zz * (2.0 + 3.0 * zz);     /* Eqn 2.19 Knigge's thesis */
  }
  else if (itype == PTYPE_DISK)
  {                             /* It was an unscattered photon from the disk */
    zz = fabs (xxspec[nspec].lmn[2]);
    pp->w *= zz * (2.0 + 3.0 * zz);     /* Eqn 2.19 Knigge's thesis */
  }
  else if (pp->nres > -1 && pp->nres < NLINES)  // added < NLINES condition for macro atoms (SS)
  {

/* It was a wind photon.  In this case, what we do depends
on whether it is a photon which arose via line radiation or 
some other process.

If geo.scatter_mode==SCATTER_MODE_ISOTROPIC then there is no need 
to reweight.  This is the isotropic assumption.  Otherwise, one
needs to reweight
*/

    if (geo.scatter_mode == SCATTER_MODE_THERMAL)
    {

      dvds = dvwind_ds (pp);
      ishell = pp->grid;
      tau = sobolev (&w[ishell], pp->x, -1.0, lin_ptr[pp->nres], dvds);
      if (tau > 0.0)
        pp->w *= (1. - exp (-tau)) / tau;
      tau = 0.0;
    }

/* But in any event we have to reposition wind photons so that they don't go through
the same resonance again */

    reposition (pp);            // Only reposition the photon if it was a wind photon
  }

  if (tau > TAU_MAX)
    istat = P_ABSORB;           /* Check to see if tau already too large */
  else if (geo.binary == TRUE)
    istat = hit_secondary (pp); /* Check to see if it hit secondary */


/* Now we can actually extract the reweighted photon */

  while (istat == P_INWIND)
  {
    istat = translate (w, pp, 20., &tau, &nres);
    icell++;

    istat = walls (pp, &pstart, normal);
    if (istat == -1)
    {
      Error ("Extract_one: Abnormal return from translate\n");
      break;
    }

    if (pp->w < weight_min)
    {
      istat = P_ABSORB;         /*This photon was absorbed within the wind */
      break;
    }

    if (istat == P_HIT_STAR)
    {                           /* It was absorbed in the photosphere */
      break;
    }
    if (istat == P_HIT_DISK)
    {                           /* It was absorbed in the disk */
      break;
    }
    if (istat == P_SCAT)
    {                           /* Cause the photon to scatter and reinitilize */
      break;
    }
  }

  if (istat == P_ESCAPE)
  {

    if (!(0 <= tau && tau < 1.e4))
      Error_silent ("Warning: extract_one: ignoring very high tau  %8.2e at %g\n", tau, pp->freq);
    else
    {
      k = (pp->freq - xxspec[nspec].freqmin) / xxspec[nspec].dfreq;

      /* Force the frequency to be in range of that recorded in the spectrum */

      if (k < 0)
        k = 0;
      else if (k > NWAVE - 1)
        k = NWAVE - 1;


      lfreqmin = log10 (xxspec[nspec].freqmin);
      lfreqmax = log10 (xxspec[nspec].freqmax);
      ldfreq = (lfreqmax - lfreqmin) / NWAVE;

      /* find out where we are in log space */
      k1 = (log10 (pp->freq) - log10 (xxspec[nspec].freqmin)) / ldfreq;
      if (k1 < 0)
      {
        k1 = 0;
      }
      if (k1 > NWAVE - 1)
      {
        k1 = NWAVE - 1;
      }

      /* Increment the spectrum.  Note that the photon weight has not been diminished
       * by its passage through th wind, even though it may have encounterd a number
       * of resonance, and so the weight must be reduced by tau
       */

      xxspec[nspec].f[k] += pp->w * exp (-(tau));       //OK increment the spectrum in question
      xxspec[nspec].lf[k1] += pp->w * exp (-(tau));     //And increment the log spectrum



      /* If this photon was a wind photon, then also increment the "reflected" spectrum */
      if (pp->origin == PTYPE_WIND || pp->origin == PTYPE_WIND_MATOM || pp->nscat > 0)
      {

        xxspec[nspec].f_wind[k] += pp->w * exp (-(tau));        //OK increment the spectrum in question
        xxspec[nspec].lf_wind[k1] += pp->w * exp (-(tau));      //OK increment the spectrum in question

      }


      /* Records the total distance travelled by extracted photon if in reverberation mode */
      if (geo.reverb != REV_NONE)
      {
        if (pstart.nscat > 0 || pstart.origin > 9 || (pstart.nres > -1 && pstart.nres < nlines))
        {                       //If this photon has scattered, been reprocessed, or originated in the wind it's important
          pstart.w = pp->w * exp (-(tau));      //Adjust weight to weight reduced by extraction
          stuff_v (xxspec[nspec].lmn, pstart.lmn);
          delay_dump_single (&pstart, nspec);   //Dump photon now weight has been modified by extraction
        }
      }


    }

  }


  if (istat > -1 && istat < 9)
    xxspec[nspec].nphot[istat]++;
  else
    Error
      ("Extract: Abnormal photon %d %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e\n",
       istat, pp->x[0], pp->x[1], pp->x[2], pp->lmn[0], pp->lmn[1], pp->lmn[2]);

  return (istat);
}
