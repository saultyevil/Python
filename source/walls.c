
/***********************************************************/
/** @file  walls.c
 * @author ksl
 * @date   May, 2020  
 *
 * @brief  This file contains the routine which determines
 * whether a there is a boundary between two photon 
 * positions 
 *
 *
 * ### Notes ###
 *
 * This was split off from photon2d.c because the logic is 
 * quite tricky for vertically extended disks.  
 *
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atomic.h"
#include "python.h"

double xnorth[] = {
  0., 0., 1.
};

double xsouth[] = {
  0., 0., -1.
};


/**********************************************************/
/**
 * @brief      determines whether the photon has reached a boundary, the central object, the disk or
 * reached the edges of the wind and, if so, modify the position of the photon, and return the appropriate
 * status.
 *
 * @param [in, out] PhotPtr  p   A photon at its new proposed location.  On exiting the routine
 * this will contain the position of the photon after taking boundaries (e.g wind cones)
 * into account.
 * @param [in] PhotPtr  pold   the previous position and direction of the photon. .
 * @param [out] double *  normal   A vector when the star or disk has been hit, which contains
 * the normal for the reflecting surface at the point the photon encountered the the boundary
 * @return   A status
 *
 * * The photon has hit the star				P_HIT_STAR
 * * The photon has hit the disk				P_HIT_DISK
 * * The photon has reached the edge of grid 		P_ESCAPE
 * * The status is undeterminable                   P_ERROR
 *
 * If a photon does not fall into one of these categories, walls returns the old status, which is stored in
 * p->istat
 *
 * @details
 *
 * walls is generally called after one has calculated how far a photon can travel in a cell
 * in order to see whether it has exited the wind region.  This is necessary because some
 * grid cells are not actually part of the wind, even though they are needed so that the
 * grid is regular.
 *
 * pold is the place where the photon was before the last attempt to move the photon forward.
 * p on input is a proposed location for photon before considering whether one has hit a boundary. The
 * location of p is either at the edge of a cell, or at the position of a resonance.  So pold should
 * be a valid position for the photon, but p may need to be adjusted.
 *
 * If one of the walls has been hit, the routine moves the photon back to that wall, but does not
 * othewise changed it.
 *
 * The routine also calculates the normal to the surface that was hit, which is intended to
 * be used by trans_phot to redirect the photon (if it has hit the disk or the star and is to
 * be scattered at that surface).
 *
 * ### Notes ###
 *
 * Note that p and pold must be travelling in the same directon to work properly.  There is an
 * error check at the beginning of the routine to make sure this is the case.  
 *
 *
 **********************************************************/
int
walls (p, pold, normal)
     PhotPtr p, pold;
     double *normal;
{
  double r, rho, rho_sq;
  double xxx[3];
  double s_disk, s_star, z;
  double theta, phi;
  double xpath;
  int hit_disk;

  /* Check to see if the photon has hit the star. If so
   * put the photon at the star surface and use that position
   * to determine the normal to the surface, the assumption
   * being that the star is located at the center of the
   * coordinate grid.
   * 
   * Note that we check both whether the new position is inside
   * the star and whether the path the photon has travelled 
   * passed thourht the star. 
   */

  r = dot (p->x, p->x);
  s_star = ds_to_sphere (geo.rstar, pold);
  vsub (p->x, pold->x, xxx);
  xpath = length (xxx);

  if (r < geo.rstar_sq || (s_star < VERY_BIG && xpath > s_star))
  {


    stuff_phot (pold, p);
    move_phot (p, s_star);


    stuff_v (p->x, normal);
    return (p->istat = P_HIT_STAR);
  }

  /* Check to see if it has hit the disk.
   * 
   * Deal with the simpler case of a flat disk first
   */

  if (geo.disk_type == DISK_FLAT && p->x[2] * pold->x[2] < 0.0)
  {                             /* Then the photon crossed the xy plane and probably hit the disk */
    s_disk = (-(pold->x[2])) / (pold->lmn[2]);

    if (s_disk < 0 && fabs (pold->x[2]) < wmain[pold->grid].dfudge && pold->lmn[2] * p->lmn[2] < 0.0)
    {
      Error ("walls: Reposition error\n");
      return (p->istat = P_REPOSITION_ERROR);
    }

    if (s_disk < 0)
    {
      Error ("walls: distance %g<0. Position %g %g %g \n", s_disk, p->x[0], p->x[1], p->x[2]);
      return (-1);
    }

    /* Check whether it hit the disk plane beyond the geo.diskrad**2 */
    vmove (pold->x, pold->lmn, s_disk, xxx);

    if (dot (xxx, xxx) < geo.diskrad_sq)
    {                           /* The photon has hit the disk */
      stuff_phot (pold, p);     /* Move the photon to the point where it hits the disk */
      move_phot (p, s_disk - DFUDGE);

      /* Now fill in the direction for the normal to the surface */
      if (pold->x[2] > 0)
      {
        stuff_v (xnorth, normal);
      }
      else
      {
        stuff_v (xsouth, normal);
      }
      return (p->istat = P_HIT_DISK);
    }

  }

  else if (geo.disk_type == DISK_VERTICALLY_EXTENDED)
  {
    /* For a vertically extended disk these means checking whether
     * we are inside the maximum radius of the disk and then lookng
     * comparing the z position of z to the height of the disk at
     * that point
     *
     * What we want to know is whether one hits the disk along
     * the line of sight between pold and p.  
     */

    rho = sqrt (pold->x[0] * pold->x[0] + pold->x[1] * pold->x[1]);
    z = zdisk (rho);
    if (z - fabs (pold->x[2]) > 1 && rho < geo.diskrad)
    {
      Error ("walls: %d The previous position %11.4e %11.4e %11.4e is inside the disk by %e\n", pold->np, pold->x[0], pold->x[1],
             pold->x[2], z - fabs (pold->x[2]));

      if (modes.save_photons)
      {
        save_photons (pold, "walls:old");
        save_photons (p, "walls:new");
        Diag ("walls: %d At rho %11.4e, zdisk is %e and zpos %e is inside the disk by %e\n", pold->np, rho, z, fabs (pold->x[2]),
              z - fabs (pold->x[2]));
      }

    }

    if (pold->lmn[0] != p->lmn[0])
    {
      Error ("The assumption is that p and pold are travelling along the same line of sight\n");
    }


    /* Check if the new position (p) is inside the disk */


    rho = sqrt (p->x[0] * p->x[0] + p->x[1] * p->x[1]);
    z = zdisk (rho);

    if (rho < geo.diskrad && fabs (p->x[2]) <= z)
    {
      /* This is the case where the proposed photon is within the disk */

      p->istat = P_HIT_DISK;
    }

    s_disk = ds_to_disk (pold, 0, &hit_disk);   /* The 0 imples that s cannot be negative */

    if (s_disk > 0 && p->ds > s_disk)
    {
      /* This is the case where the photon hits the disk along 
         the line of sight between pold and p
       */
      p->istat = P_HIT_DISK;
    }


    if (p->istat == P_HIT_DISK)
    {


      stuff_phot (pold, p);
      move_phot (p, s_disk - DFUDGE);

      rho = sqrt (p->x[0] * p->x[0] + p->x[1] * p->x[1]);
      /* This leaves the photon just outside the disk */

      /* Finally, we must calculate the normal to the disk at this point to be able to calculate the scattering direction */

      phi = atan2 (p->x[1], p->x[0]);

      if (hit_disk == DISK_HIT_EDGE)
      {
        normal[0] = cos (phi);
        normal[1] = sin (phi);
        normal[2] = 0;
      }
      else
      {
        theta = atan ((zdisk (rho * (1. + EPSILON)) - zdisk (rho)) / (EPSILON * rho));

        normal[0] = (-cos (phi) * sin (theta));
        normal[1] = (-sin (phi) * sin (theta));
        normal[2] = cos (theta);

        if (p->x[2] < 0)
        {
          normal[2] *= -1;
        }
      }

      if (modes.save_photons)
      {
        Diag ("NORMAL   %d   %e %e %e (%e %e %e %e %e)\n", p->np, normal[0], normal[1], normal[2], theta,
              (zdisk (rho * (1. + EPSILON)) - z), (EPSILON * rho), rho, EPSILON);

      }


      return (p->istat = P_HIT_DISK);
    }
  }
  /* At this point we know the photon has not hit the disk or the star, so we now
   * need to check if it has escaped the grid.  See note above regarding whether
   * we ought to be checking this differently.  This definition is clearly coordinate
   * system dependent.
   */

  rho_sq = (p->x[0] * p->x[0] + p->x[1] * p->x[1]);
  if (rho_sq > geo.rmax_sq)
    return (p->istat = P_ESCAPE);       /* The photon is coursing through the universe */

  if (fabs (p->x[2]) > geo.rmax)
    return (p->istat = P_ESCAPE);

  return (p->istat);            /* The photon is still in the wind */
}
