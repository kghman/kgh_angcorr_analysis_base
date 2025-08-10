#include "TargetLadder.h"

//for use in an MC to include shadowing effects due to target ladder geometry
//assuming an azimuthally symmetric target mount with the following geometry:
/*
                   fsw (far side width)
      ___         ___
     |   \       /   |
     |    \     /    | th (thickness)
     |_____=====_____|
            tw   tsw (target side width)
       (target width)

There is further some rotation in theta (around the vertical axis to the rxn plane).

The assumption is also made that the beam strikes the center of the target with
dimensions negligible compared to those of the mount itself.

In addition there are magnets on either side, in the following configuration:


   |     |
   |     |
   |  M  |
   |     |
   |  T  |
   |     |
   |  M  |
   |     |
   |     |

where the magnets are rectangular, and are defined by the following quantities:

         mw (magnet width)
     ___________
    |           | mth (magnet thickness) as well as md (magnet depth)
    |___________|
          ^
          |  medft (magnet edge distance from target)
          |
          x  TARGET


*/

using namespace std;

TargetLadder::TargetLadder(double tw, double th, double tsw, double fsw, double var,
			   double medft, double mth, double md, double mw) {

  target_width = tw;
  thickness = th;
  target_side_width = tsw;
  far_side_width = fsw;
  vert_axis_rotation = var;

  magnet_edge_distance_from_target = medft;
  magnet_thickness = mth;
  magnet_depth = md;
  magnet_width = mw;

  if (medft*mth*md*mw==0) include_magnets=false;
  else                    include_magnets=true;

  shadow_min_theta = atan2(th,(tw/2 + tsw - fsw)); shadow_min_theta = M_PI/2 - shadow_min_theta;
  shadow_max_theta = M_PI/2; //no matter what

  shadow_min_phi = 0;
  shadow_max_phi = 2*M_PI;

  if (include_magnets) {

    magnet_shadow_max_theta = atan2(md,medft);
    magnet_shadow_min_theta = -1*magnet_shadow_max_theta;
    magnet_shadow_max_theta = M_PI/2+magnet_shadow_max_theta;
    magnet_shadow_min_theta = M_PI/2+magnet_shadow_min_theta;

    top_magnet_shadow_max_phi = atan2(mw/2,medft);  
    top_magnet_shadow_min_phi = -1*top_magnet_shadow_max_phi;
    top_magnet_shadow_max_phi += M_PI/2;
    top_magnet_shadow_min_phi += M_PI/2;

    bot_magnet_shadow_min_phi = top_magnet_shadow_min_phi + M_PI;
    bot_magnet_shadow_max_phi = top_magnet_shadow_max_phi + M_PI;

  }

}

TargetLadder::~TargetLadder() {};

bool TargetLadder::DoesEscape(double theta, double phi) {

  //first we have to rotate the input event as though facing forward...
  double theta_rot = acos(cos(theta)*cos(vert_axis_rotation) + cos(phi)*sin(theta)*sin(vert_axis_rotation));
  double phi_rot   = atan2(sin(phi)*sin(theta),cos(phi)*sin(theta)*cos(vert_axis_rotation)-cos(theta)*sin(vert_axis_rotation));

  if (phi_rot < 0) phi_rot += 2*M_PI;

  if (!include_magnets &&
      theta_rot > shadow_min_theta && theta_rot < shadow_max_theta &&
      phi_rot   > shadow_min_phi   && phi_rot   < shadow_max_phi)
    return false;
  else if (include_magnets &&
	   (theta_rot > shadow_min_theta && theta_rot < shadow_max_theta &&
	    phi_rot   > shadow_min_phi   && phi_rot   < shadow_max_phi) ||
	   (theta_rot > magnet_shadow_min_theta   && theta_rot < magnet_shadow_max_theta &&
	    phi_rot   > top_magnet_shadow_min_phi && phi_rot   < top_magnet_shadow_max_phi) ||
	   (theta_rot > magnet_shadow_min_theta   && theta_rot < magnet_shadow_max_theta &&
	    phi_rot   > bot_magnet_shadow_min_phi && phi_rot   < bot_magnet_shadow_max_phi))
    return false;

  return true;

}
