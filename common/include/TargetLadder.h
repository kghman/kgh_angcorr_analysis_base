#include <cmath>

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

class TargetLadder {

 public:

  TargetLadder(double tw, double th, double tsw, double fsw, double var,
	       double medft=0, double mth=0, double md=0, double mw=0);
  ~TargetLadder();

  bool DoesEscape(double theta, double phi);

 private:

  double
    target_width,
    thickness,
    target_side_width,
    far_side_width,
    vert_axis_rotation;

  double
    magnet_edge_distance_from_target,
    magnet_thickness,
    magnet_depth,
    magnet_width;

  double shadow_min_theta, shadow_max_theta;
  double shadow_min_phi,   shadow_max_phi;

  double magnet_shadow_min_theta, magnet_shadow_max_theta;
  double top_magnet_shadow_min_phi, top_magnet_shadow_max_phi;
  double bot_magnet_shadow_min_phi, bot_magnet_shadow_max_phi;

  bool include_magnets;

};
