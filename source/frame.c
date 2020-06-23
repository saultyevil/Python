/***********************************************************/
/** @file  frame.c
 * @author ksl
 * @date   May, 2020
 *
 * @brief
 * These routines have to do with switching from the observer 
 * to the local frame and back
 *
 * ###Notes###
 *
 * These routines are written so that if p and p_obs point 
 * to the same variable then everything occurs in place
 *
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/**********************************************************/
/**
 * @brief  Check that a photon is in the desired frame.
 *
 * @param [in] PhotPtr  p      The photon whose frame is to be checked
 * @param [in] frame    frame  The desired frame of the photon
 * @param [in] char    *msg    A comment to print as an error message if the
 *                             photon is in the incorrect frame
 * @return 0 (False) if the photons is in the desired frame, 1 (True) if the
 *  photon is NOT in the desired frame.
 *
 * @details
 *
 * The purpose of this function is to avoid situations where photons are being
 * transformed incorrectly between the local and observer frame.
 *
 * If the photon is in the incorrect frame and save_photons mode is enabled,
 * the photon will be dumped to an external file.
 *
 **********************************************************/

int ncheck_frame = 0;

int
check_frame (p, desired_frame, msg)
     PhotPtr p;
     enum frame desired_frame;
     char *msg;
{
  if (p->frame == desired_frame)
  {
    return (0);
  }
  else if (ncheck_frame < 100)
  {
    Error ("check_frame: %s :Photon (%d) of type (%d) not in frame %d\n", msg, p->np, p->istat, desired_frame);
    ncheck_frame++;

    if (modes.save_photons)
      save_photons (p, "PhotonInIncorrectFrame");

    return (1);
  }
  else
  {
    return (0);
  }
}



/**********************************************************/
/**
 * @brief      carries out the transformation of a all the quantities
 *      in a photon structure
 *      from the observer (or global) frame to the local (or co-moving)
 *      frame
 *
 * @param [in] PhotPtr  p_in   The photon in the observer frame
 * @param [out] PhotPtr  p_out   The photon in the local frame
 *
 * @return    The routine returns 0, unless there was an error
 *
 *
 * @details
 * The routine copies all of the quantities contain in the p_in to
 * p_out, obtains the velocity at that point in the wind, and then
 * performs an in place transformation from the * global to the local
 * frame of the photon allowing for special
 * relativity.
 *
 *
 * ### Notes ###
 *
 * It p_in and p_out are the same then the transformation will
 * be performed in place.
 *
 **********************************************************/

int
observer_to_local_frame (p_in, p_out)
     PhotPtr p_in, p_out;
{
  WindPtr one;
  int ndom;
  double f_out, f_in;
  double x;
  double v[3], vel;
  double gamma;
  int i, ierr;



  ierr = check_frame (p_in, F_OBSERVER, "Observer_to_local_frame");

  /* Initialize the output photon */
  stuff_phot (p_in, p_out);
  p_out->frame = F_LOCAL;

  /* Calculate the local velocity of the wind at this position */
  one = &wmain[p_in->grid];
  ndom = one->ndom;
  vwind_xyz (ndom, p_in, v);

  vel = dot (p_in->lmn, v);


  f_in = p_in->freq;

  if (rel_mode == REL_MODE_LINEAR)
  {
    f_out = p_out->freq = f_in * (1. - vel / VLIGHT);
    return (ierr);
  }



  // beta=length(v)/VLIGHT;

  gamma = sqrt (1 - (dot (v, v) / (VLIGHT * VLIGHT)));

  f_out = p_out->freq = f_in * gamma * (1. - vel / VLIGHT);

  x = gamma / VLIGHT * (1.0 - (gamma * vel / ((gamma + 1) * VLIGHT)));

  for (i = 0; i < 3; i++)
  {
    p_out->lmn[i] = f_out / f_in * (p_in->lmn[i] - x * v[i]);
  }

  p_out->w *= (f_out / f_in);

  return (ierr);
}



/**********************************************************/
/**
 * @brief      carries out the transformation of all the quantities
 *      in a photon structure 
 *      from the local (or co-moving) frame to the observer (or global)
 *      frame
 *
 * @param [in] PhotPtr  p_in   The photon in the local frame                   
 * @param [out] PhotPtr  p_out   The photon in the global frame                  
 *
 * @return    The routine returns 0 unless there was an error       
 *
 *
 * @details
 * The routine copies all of the quantities contain in the p_in to 
 * p_out, obtains the velocity at that point in the wind, and then 
 * performs an in place transformation from the local to the observer
 * frame of the photon allowing for special 
 * relativity.
 *
 *
 * ### Notes ###
 *
 * It p_in and p_out are the same then the transformation will
 * be performed in place.
 *
 **********************************************************/


int
local_to_observer_frame (p_in, p_out)
     PhotPtr p_in, p_out;
{
  WindPtr one;
  int ndom;
  double f_out, f_in;
  double x;
  double v[3], vel;
  double gamma;
  int i;
  int ierr;




  ierr = check_frame (p_in, F_LOCAL, "local_to_observer_frame");

  /* Initialize the output photon */
  stuff_phot (p_in, p_out);
  p_out->frame = F_OBSERVER;

  /* Calculate the local velocity of the wind at this position */
  one = &wmain[p_in->grid];
  ndom = one->ndom;
  vwind_xyz (ndom, p_in, v);

  f_in = p_in->freq;

  vel = dot (p_in->lmn, v);
  if (rel_mode == REL_MODE_LINEAR)
  {
    f_out = p_out->freq = f_in / (1. - vel / VLIGHT);
    return (ierr);
  }
  gamma = sqrt (1 - (dot (v, v) / (VLIGHT * VLIGHT)));
  f_out = p_out->freq = f_in * gamma * (1. + vel / VLIGHT);


/* Need to worry about sign changes, etc. here */
  x = gamma / VLIGHT * (1.0 + (gamma * vel / ((gamma + 1) * VLIGHT)));

  for (i = 0; i < 3; i++)
  {
    p_out->lmn[i] = f_out / f_in * (p_in->lmn[i] + x * v[i]);
  }

  p_out->w *= (f_out / f_in);

  return (ierr);
}



/**********************************************************/
/**
 * @brief      carries out the transformation of a all the quantities
 *      in a photon structure
 *      from the observer (or global) frame to the local (or co-moving)
 *      frame
 *
 * @param [in] PhotPtr  p_in   The photon in the observer frame
 * @param [out] PhotPtr  p_out   The photon in the local frame
 *
 * @return    The routine returns 0, unless there was an error
 *
 *
 * @details
 * The routine copies all of the quantities contain in the p_in to
 * p_out, obtains the velocity at that point in the wind, and then
 * performs an in place transformation from the * global to the local
 * frame of the photon allowing for special
 * relativity.
 *
 *
 * ### Notes ###
 *
 * It p_in and p_out are the same then the transformation will
 * be performed in place.
 *
 **********************************************************/

int
observer_to_local_frame_disk (p_in, p_out)
     PhotPtr p_in, p_out;
{
  double f_out, f_in;
  double x;
  double v[3], vel;
  double gamma;
  int i, ierr;



  ierr = check_frame (p_in, F_OBSERVER, "Observer_to_local_frame_disk");

  /* Initialize the output photon */
  stuff_phot (p_in, p_out);
  p_out->frame = F_LOCAL;

  /* Calculate the local velocity of the disk at this position */
  vdisk (p_in->x, v);

  vel = dot (p_in->lmn, v);


  f_in = p_in->freq;

  if (rel_mode == REL_MODE_LINEAR)
  {
    f_out = p_out->freq = f_in * (1. - vel / VLIGHT);
    return (ierr);
  }



  // beta=length(v)/VLIGHT;

  gamma = sqrt (1 - (dot (v, v) / (VLIGHT * VLIGHT)));

  f_out = p_out->freq = f_in * gamma * (1. - vel / VLIGHT);

  x = gamma / VLIGHT * (1.0 - (gamma * vel / ((gamma + 1) * VLIGHT)));

  for (i = 0; i < 3; i++)
  {
    p_out->lmn[i] = f_out / f_in * (p_in->lmn[i] - x * v[i]);
  }

  p_out->w *= (f_out / f_in);

  return (ierr);
}



/**********************************************************/
/**
 * @brief      carries out the transformation of a all the quantities
 *      in a photon emitted from the surface of the disk
 *      from the local (or co-moving) frame to the observer (or global)
 *      frame
 *
 * @param [in] PhotPtr  p_in   The photon in the local frame                   
 * @param [out] PhotPtr  p_out   The photon in the global frame                  
 *
 * @return    The routine returns 0 unless the was an error         
 *
 *
 * @details
 * The routine copies all of the quantities contain in the p_in to 
 * p_out, obtains the velocity at that point in the wind, and then 
 * performs an in place transformation from the local to the observer
 * frame of the photon allowing for special 
 * relativity.
 *
 *
 * ### Notes ###
 *
 * It p_in and p_out are the same then the transformation will
 * be performed in place.  This is normally the way the 
 * routine is used.   One could have omitted p_out in this case
 * but to make the interface coisistent, both photon structures
 * are included.  
 *
 **********************************************************/

int
local_to_observer_frame_disk (p_in, p_out)
     PhotPtr p_in, p_out;
{
  double f_out, f_in;
  double x;
  double v[3], vel;
  double gamma;
  int i;
  int ierr;

  ierr = check_frame (p_in, F_LOCAL, "local_to_observer_frame_disk");

  /* Initialize the output photon */
  stuff_phot (p_in, p_out);
  p_out->frame = F_OBSERVER;


  /* Calculate the local velocity of the disk at this position */
  vdisk (p_in->x, v);
  vel = dot (p_in->lmn, v);

//OLD  f_in = p_in->freq_orig_loc;
  f_in = p_in->freq;



  if (rel_mode == REL_MODE_LINEAR)
  {
    f_out = p_out->freq = f_in / (1. - vel / VLIGHT);
    return (ierr);
  }

  gamma = sqrt (1 - (dot (v, v) / (VLIGHT * VLIGHT)));
  f_out = p_out->freq = f_in * gamma * (1. + vel / VLIGHT);

/* Need to worry about sign changes, etc. here */
  x = gamma / VLIGHT * (1.0 + (gamma * vel / ((gamma + 1) * VLIGHT)));

  for (i = 0; i < 3; i++)
  {
    p_out->lmn[i] = f_out / f_in * (p_in->lmn[i] + x * v[i]);
  }

  p_out->w *= (f_out / f_in);

  return (ierr);
}
