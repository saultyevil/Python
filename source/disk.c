

/***********************************************************/
/** @file  disk.c
 * @author ksl
 * @date   January, 2018
 *
 * @brief  Various routines related setting up the disk and calculating
 * when a photon encounters the disk surface
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"


/**********************************************************/
/** 
 * @brief      Calculate the reference temperarure for a standard Shakura-Sunyaeev disk
 *
 * @param [in] double  m   mass of the central object
 * @param [in] double  mdot   mass accretion rate of the disk
 * @param [in] double  r   radiis of the central object
 * @return     The reference temperature for the disk
 *
 * @details
 * All of the inputs are in cgs units
 *
 * ### Notes ###
 * For an SS disk,  The maximum temperature
 *    actuallly seen in the disk is 0.488 tdisk
 *
 *
 **********************************************************/

double
tdisk (m, mdot, r)
     double m, mdot, r;
{
  double t;
  t = 3. * GRAV / (8. * PI * STEFAN_BOLTZMANN) * m * mdot / (r * r * r);
  t = pow (t, 0.25);
  return (t);
}




/**********************************************************/
/** 
 * @brief      Calculate the temperature of the disk at a normalised distance x
 *
 * @param [in] double  t   reference temperature of the disk in degrees
 * @param [in] double  x   distance from center of disk in units of r/rstar
 * @return     The temperature
 *
 * The routine returns the effective temperature for the disk as a distance
 * which is measured in terms of the ratio of the actual distance to the inner
 * edge of the disk (rmin).
 *
 * Modifications to a simple steady state disk are possible, depending
 * on the variable geo.disk_tprofile:
 *
 *  - DISK_TPROFILE_READIN  the returned temperature is interpolated from values that
 *    been read in
 *
 *  In addition for a standard profile, this is where disk heading can optionally
 *  be taken into account
 *
 *
 * ###Notes###
 *
 * A reference for the standard steady state disk is Wade, 1984 MNRAS 208, 381
 *
 * An analytic profile used for Sim+05 to model certain YSO winds was removed in
 * 1907
 *
 **********************************************************/

double
teff (t, x)
     double t, x;
{
  double q = 0;
  double theat, r;
  double temp;
  double pow ();
  int kkk;



  if (x < 1)
  {
    Error ("teff: x %f less than 1.0\n", x);
    return (0.0);
  }


  if ((geo.disk_tprofile == DISK_TPROFILE_READIN) && ((x * geo.rstar) < blmod.r[blmod.n_blpts - 1]))
  {
    /* This is the case where the temperature profile is read in as an array, and so we
       simply find the array elements that bracked the requested radius and do a linear
       interpolation to calcualte the temperature at the requested radius. */
    if ((r = (x * geo.rstar)) < blmod.r[0])
    {
      return (blmod.t[0]);
    }
    else
    {
      linterp (r, &blmod.r[0], &blmod.t[0], blmod.n_blpts, &temp, 0);
      return (temp);
    }
  }
  else
  {
    /* This is a standard accretion disk */

    q = (1.e0 - pow (x, -0.5e0)) / (x * x * x);
    q = t * pow (q, 0.25e0);

    if (geo.absorb_reflect == BACK_RAD_ABSORB_AND_HEAT && geo.wcycle > 0)       /* Absorb photons and increase t so that heat is radiated
                                                                                   but only do this if there has been at least one
                                                                                   ionization cycle */
    {
      /* qdisk is initialized only once (in init_qdisk) and does not generally have the same
       * values for r as does the disk structure, whose annulae vary as the 
       * frequency limits are set. Here we just search for a radius that is just above
       * the requested r
       * 
       */
      r = x * geo.rstar;        // 04aug -- Requires fix if disk does not extend to rstar
      kkk = 1;                  // photon cannot hit the disk at r<qdisk.r[0]
      while (r > qdisk.r[kkk] && kkk < NRINGS - 1)
        kkk++;

      /* Note that disk has 2 sides */
      theat = qdisk.heat[kkk - 1] / (2. * PI * (qdisk.r[kkk] * qdisk.r[kkk] - qdisk.r[kkk - 1] * qdisk.r[kkk - 1]));

      /* T_eff is given by T_eff**4= T_disk**4+Heating/area/STEFAN_BOLTZMANN */
      q = pow (q * q * q * q + (theat / STEFAN_BOLTZMANN), 0.25);

    }
  }
  return (q);
}


/**********************************************************/
/** 
 * @brief      Calculate a reference gravity for the disk
 *
 * @param [in] double  mass   Mass of the central object
 * @param [in] double  mdot   Mass accretion rate of the disk
 * @param [in] double  rmin   Minimum radius for the disk
 * @return     A fiducial gravity for the disk
 *
 * The gravity is needed when one constructs a disk spectrum
 * from spectra from a grid of stellar atmospheres
 *
 * The input units are all cgs
 *
 * ###Notes###
 *
 * See Long & Knigge for details
 *
 *
 **********************************************************/

double
gdisk (mass, mdot, rmin)
     double mass, rmin, mdot;
{
  double g0;
  g0 = 0.625 * log10 (mass / MSOL) - 1.875 * log10 (rmin / 1.e9) + 0.125 * log10 (mdot / 1.e16);
  g0 = 5.96e5 * pow (10., g0);
  return (g0);
}


/**********************************************************/
/** 
 * @brief      Calculate the effective gravity at a specfic position in the disk
 *
 * @param [in] double  g0  reference gravity in cm s**-2
 * @param [in] double  x   distance from center in units of r/rmin
 * @return     gravity at x in cm s**-2
 *
 * The gravity is needed when one constructs a disk spectrum
 * from spectra from a grid of stellar atmospheres
 *
 *
 * ###Notes###
 *
 * See Long & Knigge for details
 *
 **********************************************************/

double
geff (g0, x)
     double g0, x;
{
  double q;
  q = (1.0e0 - pow (x, -0.5e0));
  q = pow (x, -1.875e0) * pow (q, 0.125);
  q = g0 * q;
  return (q);
}


double north[] = { 0.0, 0.0, 1.0 };


/**********************************************************/
/** 
 * @brief      Calculate the speed and velocity v at which material given
 *    specific position in the disk
 *
 * @param [in] double  x[]   A position in the disk
 * @param [out] double  v[]   The velocity at that position
 * @return     The speed
 *
 * This routine interpolates the velocities from a grid which gives
 * the run of velocity with radius. Although the we generally
 * assume keplerian velocities for the disk those velocities
 * are not calculated here.
 *
 * ###Notes###
 *
 * The routine projects the input variable x on to the xy plane
 * before it calculates velocities
 *
 * It then takes the xproduct of a unit vector pointing north and
 * the position, to determine the direction of motion at this point
 * in the disk
 *
 *
 *  Finally it rescales this to get the actual velocity.
 *
 *
 *
 *
 **********************************************************/

double
vdisk (x, v)
     double x[];
     double v[];
{
  double xhold[3];
  double r, speed;

  stuff_v (x, xhold);
  xhold[2] = 0.0;
  r = length (xhold);
  linterp (r, disk.r, disk.v, NRINGS, &speed, 0);       //interpolate in linear space
  cross (north, xhold, v);      /* The velocity vector direction is given by north x r */
  renorm (v, speed);
  return (speed);
}




/**********************************************************/
/** 
 * @brief      Calculate the height of the disk at a certain radius
 *
 * @param [in] double  r   a radial position in the disk.
 * @return     The vertical height of the disk at that point in the xy
 * 	plane.
 *
 *
 * If the disk is vertically extended then the disk surface is
 * defined in terms of the h the fraction of the disk radius
 * at the edge of the disk, and a power law exponent.  So
 * a disk thst us a sunoke wedgge wiykd gave an exponent
 * of 1. A flat but thick disk would have an exponent of 0.
 * A flared disk wouuld generally have an exponent greater than
 * 1.
 *
 * ###Notes###
 *
 * 	zdisk returns a number that will positive or zero, so one
 * 	often needs to take this account if for example one has
 * 	a photon that hits the disk below the z=0 plane.
 *
 *      zdisk does not take the maximum radius of the disk
 *      into account
 *
 **********************************************************/

double
zdisk (r)
     double r;
{
  double z;
  z = geo.disk_z0 * pow (r / geo.diskrad, geo.disk_z1) * geo.diskrad;
  return (z);
}


int ds_to_disk_init = 0;
struct photon ds_to_disk_photon;
struct plane diskplane, disktop, diskbottom;


/**********************************************************/
/**
 * @brief      Calculates the distance that a photon
 *  	would need to travel from its current position to hit disk.
 *
 * @param [in] struct photon *  p   a photon pointer.
 * @param [in] int  allow_negative   if nonzero, permits a
            negative distance to
 *          be returned
 * @return     The distance to the disk.
 *
 *
 * Usually, ds_to_disk returns the distance along the line of
 * sight to the disk in the direction a photon is currently travelling.
 *
 * Usually, if the photon misses the disk going in the positive
 * direction, a very large number is returned.
 *
 * If allow_negative is non-zero (true), then ds_to_disk returns
 * a negative distance if it has not hit the disk going in the
 * positive direction.
 *
 * The routine allows both for a flat disk, and a vertically
 * extended disk
 *
 *
 *
 * ###Notes###
 *
 * The z-height of a vertically extended disk is defined by
 * zdisk.  The outside edge of the disk is assumed to be
 * a cylinder at geo.diskrad.  The routine does not say which
 * surface of the disk has been hit (although for vertically
 * extended disks this might be a good idea).
 *
 * The need to allow for negative distances arises
 * because several of the parameterization for the wind (SV, KWD) depend
 * on the distance between the current position and footpoint
 * along a streamline (and these streamlines are defined in the
 * outgoing direction).
 *
 * If the position of the photon is on the disk the routine returns
 * 0 and it is up to the calling routine to interpret what to do
 * in this situation.  For very large distances one may need to worry
 * about round off errors if moving a photon to the disk position.
 **********************************************************/

double
ds_to_disk (p, allow_negative, hit)
     struct photon *p;
     int allow_negative;
     int *hit;
{

  /*
   * Respectively, the distance of the photon to the disk plane, the
   * distance to the plane that defines the top of a pill box
   * surrounding the disk, the distance to the bottom, and the distance
   * to the cylinder that forms the outer edge of a vertically extended
   * disk. s will be the distance the photon must be moved to hit a
   * boundary.
   */

  double s_diskplane, s_top, s_bot, s_cyl, s_cyl2, s;
  double z_cyl, z_cyl2;

  /* A variable used to describe whether the photon is inside a vertically extended disk */
  double delta_z;

  /* The next set of variables represent rho at various positions */
  double r_diskplane, r_top, r_bot, r_hit;
  double r_phot;

  double smin, smax;

  struct photon phit;

  int location = -1;
  smin = smax = 0.0;


  /* Simply return is there is no disk */
  if (geo.disk_type == DISK_NONE)
  {

    *hit = DISK_MISSED;
    return (VERY_BIG);
  }


  /*
   * Initialize 3 structures that define the plane of the disk, and two
   * other planes that encompass the disk.  
   */

  if (ds_to_disk_init == 0)
  {
    diskplane.x[0] = diskplane.x[1] = diskplane.x[2] = 0.0;
    diskplane.lmn[0] = diskplane.lmn[1] = 0.0;
    diskplane.lmn[2] = 1.0;


    disktop.x[0] = disktop.x[1] = 0.0;
    disktop.x[2] = geo.diskrad * geo.disk_z0;
    disktop.lmn[0] = disktop.lmn[1] = 0.0;
    disktop.lmn[2] = 1.0;

    diskbottom.x[0] = diskbottom.x[1] = 0.0;
    diskbottom.x[2] = (-geo.diskrad * geo.disk_z0);
    diskbottom.lmn[0] = diskbottom.lmn[1] = 0.0;
    diskbottom.lmn[2] = 1.0;

    ds_to_disk_init++;

  }

  /* Find where the photon would hit the disk for the case of a FLAT disk */

  s_diskplane = ds_to_plane (&diskplane, p);
  stuff_phot (p, &phit);
  move_phot (&phit, s_diskplane);
  r_diskplane = sqrt (phit.x[0] * phit.x[0] + phit.x[1] * phit.x[1]);


  if (geo.disk_type == DISK_FLAT)
  {
    if (r_diskplane > geo.diskrad)
    {
      *hit = DISK_MISSED;
      return (VERY_BIG);
    }
    else if (s_diskplane > 0 || allow_negative)
    {
      if (p->x[2] > 0)
      {
        *hit = DISK_HIT_TOP;
      }
      else
      {
        *hit = DISK_HIT_BOT;
      }
      return (s_diskplane);
    }
    else
    {
      *hit = DISK_MISSED;
      return (VERY_BIG);
    }
  }


  /*
   * Begin the calculation for a VERTICALLY EXTENDED DISK. 
   * 
   * For the vertically extended disk we have to keep track of the
   * smallest positive value and the smallest (in absolute terms
   * negative value).
   * 
   * We would like to avoid actually having to calculate the intercept to
   * the disk if we can because is time consuming.  So we first try
   * various checks to see if the photons misses the disk.  
   * The first check we make is to see whether the photon misses 
   * the pill box that surrounds the disk.
   */

  r_phot = sqrt (p->x[0] * p->x[0] + p->x[1] * p->x[1]);

  s_top = ds_to_plane (&disktop, p);
  stuff_phot (p, &phit);
  move_phot (&phit, s_top);
  r_top = sqrt (phit.x[0] * phit.x[0] + phit.x[1] * phit.x[1]);

  s_bot = ds_to_plane (&diskbottom, p);
  stuff_phot (p, &phit);
  move_phot (&phit, s_bot);
  r_bot = sqrt (phit.x[0] * phit.x[0] + phit.x[1] * phit.x[1]);

  /* When a photon is outside of the cylinder of the disk
     there are instances where we need to know 
     both of the cylinder boundaries.   

     In this case s_cyl will be the distance to the closest
     intersection and s_cyl2 will be the distance to the further
     intersection.  

     z_cyl and z_cyl2 are the zheights of the intersections with
     the cylinders
   */

  s_cyl = ds_to_cylinder (geo.diskrad, p);
  stuff_phot (p, &phit);
  move_phot (&phit, s_cyl);
  z_cyl = phit.x[2];

  z_cyl2 = VERY_BIG;
  s_cyl2 = VERY_BIG;

  if (r_phot > geo.diskrad)
  {
    move_phot (&phit, DFUDGE);
    s_cyl2 = ds_to_cylinder (geo.diskrad, &phit);
    move_phot (&phit, s_cyl2);
    z_cyl2 = phit.x[2];
    s_cyl2 += s_cyl;
  }



  /*
   * So now we know the distances to the top and bottom planes, as well
   * as to the disk plane, and we also know the distance to the
   * cylinder.  There are going to be two possibilites that we want to
   * deal with, the normal one where the photon is outside of the
   * disk, and the slightly perverse one where the photon is already
   * inside the disk
   *
   * Note that if delta_z below is positive the photon is inside the 
   * disk
   */

  delta_z = zdisk (r_phot) - fabs (p->x[2]);

  /* We assume if we are within a cm of the disk surface we don't need 
   * do anything.
   */

  if ((r_phot < geo.diskrad) && (fabs (delta_z) < 1.0))
  {
    return (0);
  }

  /*
   * Handle the case where the photon is already inside the
   * disk.  The main issue here is we want to push the
   * photon back where it came from now allow it to go
   * forward.
   *  
   * We will make the assumption that if you are moving away from 
   * the disk plane (z=0) you want to go
   * forward, but if you are moving towards the plane you want to go backwards
   *
   * Note that for all of these cases we are ultimately going to solve
   * for where we hit the disk, and so what is done in this section
   * is to resolve what the minimum and maximum distance possiblities are.
   */


  if ((r_phot < geo.diskrad) && delta_z > 0)
  {

    if (p->x[2] * p->lmn[2] > 0)
    {
      /* The photon is moving away from the disk plane and we want to continue in that direction */
      location = 1;

      smin = 0;

      smax = s_top;
      if (s_bot > 0)
      {
        smax = s_bot;
      }
      if (s_cyl < smax)
      {
        smax = s_cyl;
      }

    }
    else
    {
      /* We are moving toward the disk plane, and we want to go back to a position on the disk */

      location = 2;

      smax = 0;

      /* One or the other of s_top or s_bot should be negative */

      if (s_top < 0)
      {
        smin = s_top;
      }
      else
      {
        smin = s_bot;
      }

      /* This ignores the possibility that a photon inside the disk goes out 
         the edge of the disk, so we have to deal with this possibility */
      stuff_phot (p, &phit);
      phit.lmn[0] *= (-1);
      phit.lmn[1] *= (-1);
      phit.lmn[2] *= (-1);
      s_cyl = ds_to_cylinder (geo.diskrad, &phit);
      s_cyl *= (-1);

      /* s_cyl should be negative and but it may be closer to 0 than 
         where the photons hit the disk. If that is the case we want to
         use s_cyl */

      if (s_cyl > smin)
      {
        smin = s_cyl;
      }
    }

  }

  /* At this point we haave settled on the limits for solving
     for the intercept for the case where the photon was inside 
     the disk.  We now address the case where the 
     This is the case where the photon location is outside 
     the disk.  

     Begin work on case where the photon is OUTSIDE the 
     disk

     The strategy here is to first try to eliminate the 
     possibility that the photon hits the disk and 
     to return a large value.  Only if we cannot do this
     will we actually set limits for where the photon might hit
     and calculate the intersection.
   */

  else if (r_phot > geo.diskrad || fabs (p->x[2]) > zdisk (r_phot))
  {

    /* Identify photons moving away from the disk that do not hit 
       the edge of the disk */

    if (s_diskplane < 0 && fabs (z_cyl) > disktop.x[2])
    {
      *hit = DISK_MISSED;
      return (VERY_BIG);
    }

    /* Identify photons that are outside the pill box, and 
       never even hit the cylinder that encloses the disk
     */
    else if (s_cyl == VERY_BIG)
    {
      *hit = DISK_MISSED;
      return (VERY_BIG);
    }

    /* Identify photons that are soutside the disk radius and the photon hits 
       the outer cylindrical edge of the disk */
    else if (r_diskplane > geo.diskrad && fabs (z_cyl) < disktop.x[2])
    {
      *hit = DISK_HIT_EDGE;
      return (s_cyl);
    }

    /* Identify photons outside the pill box, and do not  hit either
       the disk or the edges of the disk
     */
    else if (r_diskplane > geo.diskrad && fabs (z_cyl) > disktop.x[2] && fabs (z_cyl2) > disktop.x[2])
    {
      *hit = DISK_MISSED;
      return (VERY_BIG);
    }

    /* Finally set limits for photons that we could not eliminate
     * simply.
     *
     * Set limits smin and smax on what the possibilities are. The
     * maximum is the distance to the disk or cylinder
     * whichever is smaller.  The minimum is the current
     * position if we are inside the pillbox.
     */
    else
    {

      location = 0;


      /* We certainly do not need to go further than
         the distance to the disk and (if we are inside the cylinder) 
         the distance to the cylinder wall.
       */

      smax = VERY_BIG;
      if (s_diskplane > 0)
      {
        smax = s_diskplane;
      }
      else if (r_phot < geo.diskrad && s_cyl < smax)
      {
        smax = s_cyl;
      }
      else
      {
        smax = 2 * geo.diskrad; /* this is a cheat */
      }

      /* For the minimum, we need to either to be inside the
         pill box (in which case smin is 0), or we can take correct boundary of the
         pill box
       */

      if (fabs (p->x[2]) < disktop.x[2] && r_phot < geo.diskrad)
      {
        smin = 0;
      }
      else
      {

        smin = smax;
        if (s_top > 0 && s_top < smin)
        {
          smin = s_top;
        }
        if (s_bot > 0 && s_bot < smin)
        {
          smin = s_bot;
        }
        if (s_cyl < smin)
        {
          smin = s_cyl;
        }
      }
    }

  }

  /* This is a failure bexause it implies we have a case we have not accounted for */
  else
  {
    Error ("ds_to_disk: Unknnown situation\n");

  }

  /* If we have reached this point in the program we have to SOLVE FOR 
     THE INTERCEPT.  Note that we can have arrived here if the photon
     was inside the disk, or if the photon was outside the disk, but
     we could not eliminate the possibility that it actually hit the
     disk.
   */


  /*
   * Now we need to find exactly where we have hit the disk.  This is
   * not trivial.  Basically we need to bracket the possibilities
   * 
   * At this point we know the photon is on a path that passes through (or
   * passed through the disk) and we must locate it. There are two
   * possibilities.  It hit the disk edge, or it hit the face of the
   * disk.  Most of the time it will hit the disk face so we will
   * calculate this first, and then check if s is less than s_diskplanerad.
   * 
   * To setup rtsafe, we have to find distances which bracket what we
   * want.  smax should be no more than the distance to the disk_plane
   * (assuming we have a photon headed toward the plane)
   * 
   * Note that the next little section is intended to select between the
   * top and bottom face of the disk
   */


  stuff_phot (p, &ds_to_disk_photon);
  move_phot (&ds_to_disk_photon, smin);


  if ((smax - smin) > 0.)
  {
    s = zero_find (disk_height, 0.0, smax - smin, 1e-8);
  }
  else
  {
    s = zero_find (disk_height, smax - smin, 0.0, 1e-8);
  }

  /* Check if our solution is the maxium or the minimum distance.
     This could mean that we have not selected smin or smax
     correctly, but (see #763) there are instances where we
     do have this condition where we are intersecting the outer
     edge of the disk
   */

  if (s == smin || s == smax)
  {
    stuff_phot (p, &phit);
    move_phot (&phit, s);
    r_hit = sqrt (phit.x[0] * phit.x[0] + phit.x[1] * phit.x[1]);


    Error ("ds_to_disk: Phot %d loc %d smin %e smax %e s at limit %e r_phot %e z_phot %e zdisk %e r_hit %e z_hit %e zdisk_hit %e\n",
           p->np, location, smin, smax, s, r_phot, fabs (p->x[2]), zdisk (r_phot), r_hit, fabs (phit.x[2]), zdisk (r_hit));
    if (modes.save_photons)
    {
      save_photons (p, "ds_to_disk");
      Diag ("ds_to_disk: Phot %d loc %d smin %e smax %e s at limit %e r_phot %e z_phot %e zdisk %e r_hit %e z_hit %e zdisk_hit %e\n",
            p->np, location, smin, smax, s, r_phot, fabs (p->x[2]), zdisk (r_phot), r_hit, fabs (phit.x[2]), zdisk (r_hit));
    }
  }

  if (p->x[2] > 0)
  {
    *hit = DISK_HIT_TOP;
  }
  else
  {
    *hit = DISK_HIT_BOT;
  }


  s += smin;



  return (s);
}



/**********************************************************/
/**
 * @brief      Calculate the change in disk height as s function of s
 *
 * @param [in] double  s  A distance along the path of a photon
 *           void params  unused variable required to present the correct function to gsl

 * @param [out] double *  value   the disk height at s
 * @return     N/A
 *
 * Used in ds_to_disk by zero_find
 * the place a photon hits a vertically extended disk
 *
 * ###Notes###
 *
 *
 *  The function we are trying to zero is the
 *  difference between the height of the disk and the z height of the photon.
 *  Modified 2019 to allow use of gsl rather than old numerical recipies
 *
 *  Note that the disk height is defined even outside or the region
 *  actually occupied by the disk, so one needs to make sure
 *  this is not called unless there is a valid intersection with
 *  the disk, or alternatively one needs to check afterwards
 *
 **********************************************************/


double
disk_height (double s, void *params)
{
  struct photon phit;
  double z1, r1;

  (void) params;

  stuff_phot (&ds_to_disk_photon, &phit);
  move_phot (&phit, s);
  r1 = sqrt (phit.x[0] * phit.x[0] + phit.x[1] * phit.x[1]);
  z1 = zdisk (r1) - fabs (phit.x[2]);
  //this is the function

  return (z1);

}
