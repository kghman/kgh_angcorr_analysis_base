/*

  Class which represents a single MMM detector in the SABRE array at FSU. Origial code by KGH, re-written by
  GWM.

  Distances in meters, angles in radians.

  The channel arrays have four points, one for each corner. The corners are
  as follows, as if looking BACK along beam (i.e. from the target's pov):

  0---------------------1      y
  |                     |      ^
  |                     |      |
  |                     |      |
  |                     |      |
  |                     |      |---------> x
  3---------------------2      
                               (z is hence positive along beam direction) 

  The channel numbers, also as looking back from target pov, are:

  >> rings are 0 -- 15 from inner to outer:

    15 -------------------
    14 -------------------
    13 -------------------
       .
       .
       .
     2 -------------------
     1 -------------------
     0 -------------------

  >> wedges are 0 -- 7 moving counterclockwise (i.e. pos phi)

      7 6 ... 1 0
     | | |   | | |
     | | |   | | |
     | | |   | | |
     | | |   | | |
     | | |   | | |
     | | |   | | |


  >> Note that the detector starts centered on the x-axis (central phi = 0) untilted, and then is rotated to wherever the frick
  	 it is supposed to go; phi = 90 is centered on y axis, pointing down towards the bottom of the scattering chamber

  -- GWM, Dec 2020; based on the og code from kgh

*/

#ifndef SABREDETECTOR_H
#define SABREDETECTOR_H

#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include "Vec3.h"
#include "Rotation.h"

class SabreDetector {
public:

	SabreDetector();
	SabreDetector(double Rin, double Rout, double deltaPhi_flat, double phiCentral, double tiltFromVert, double zdist, double xdist=0, double ydist=0);
	~SabreDetector();

	/*Return coordinates of the corners of each ring/wedge in SABRE*/
	inline Mask::Vec3 GetRingFlatCoords(int ch, int corner) { return CheckRingLocation(ch, corner) ? m_ringCoords_flat[ch][corner] : Mask::Vec3(); };
	inline Mask::Vec3 GetWedgeFlatCoords(int ch, int corner) { return CheckWedgeLocation(ch, corner) ? m_wedgeCoords_flat[ch][corner] : Mask::Vec3(); };
	inline Mask::Vec3 GetRingTiltCoords(int ch, int corner) { return CheckRingLocation(ch, corner) ? m_ringCoords_tilt[ch][corner] : Mask::Vec3(); };
	inline Mask::Vec3 GetWedgeTiltCoords(int ch, int corner) { return CheckWedgeLocation(ch, corner) ? m_wedgeCoords_tilt[ch][corner] : Mask::Vec3(); };

	inline double GetRingTiltedTheta(int ch, int corner) { return CheckRingLocation(ch, corner) ? m_ringCoords_tilt[ch][corner].GetTheta() : 0; };
	inline double GetWedgeTiltedPhi(int ch, int corner) { return CheckWedgeLocation(ch, corner) ? m_wedgeCoords_tilt[ch][corner].GetPhi() : 0; };

	Mask::Vec3 GetTrajectoryCoordinates(double theta, double phi);
	std::pair<int, int> GetTrajectoryRingWedge(double theta, double phi);
	Mask::Vec3 GetHitCoordinates(int ringch, int wedgech);

	//kgh add--20201214 for simpler access (basic types)
	double GetHitTheta(int ringch, int wedgech);
	double GetHitPhi(int ringch, int wedgech);

	/*Basic getters*/
	inline int GetNumberOfWedges() { return m_nWedges; };
	inline int GetNumberOfRings() { return m_nRings; };
	inline double GetInnerRadius() { return m_Rinner; };
	inline double GetOuterRadius() { return m_Router; };
	inline double GetPhiCentral() { return m_phiCentral; };
	inline double GetTiltAngle() { return m_tilt; };
	inline Mask::Vec3 GetTranslation() { return m_translation; };


private:

	/*Class constants*/
	static constexpr int m_nRings = 16;
	static constexpr int m_nWedges = 8;
	static constexpr double deg2rad = M_PI/180.0;
	/*These are implicitly the width of the spacing between detector active strips*/
	static constexpr double POSITION_TOL = 0.0001; //0.1 mm position tolerance
	static constexpr double ANGULAR_TOL = 0.1*M_PI/180.0; // 0.1 degree angular tolerance

	void CalculateCorners();

	/*Performs the transformation to the tilted,rotated,translated frame of the SABRE detector*/
	inline Mask::Vec3 TransformToTiltedFrame(Mask::Vec3& vector) { return (m_ZRot*(m_YRot*vector)) + m_translation; };

	/*Determine if a given channel/corner combo is valid*/
	inline bool CheckRingChannel(int ch) { return (ch<m_nRings && ch>=0) ? true : false; };
	inline bool CheckWedgeChannel(int ch) { return (ch<m_nWedges && ch >=0) ? true : false; };
	inline bool CheckCorner(int corner) { return (corner < 4 && corner >=0) ? true : false; };
	inline bool CheckRingLocation(int ch, int corner) { return CheckRingChannel(ch) && CheckCorner(corner); };
	inline bool CheckWedgeLocation(int ch, int corner) { return CheckWedgeChannel(ch) && CheckCorner(corner); };

	/*
		For all of the calculations, need a limit precision to determine if values are actually equal or not
		Here the approx. size of the strip spacing is used as the precision.
	*/
	inline bool CheckPositionEqual(double val1,double val2) { return fabs(val1-val2) > POSITION_TOL ? false : true; };
	inline bool CheckAngleEqual(double val1,double val2) { return fabs(val1-val2) > ANGULAR_TOL ? false : true; };

	/*Determine if a hit is within the bulk detector*/
	inline bool IsInside(double r, double phi) { 
		double phi_1 = m_deltaPhi_flat/2.0;
		double phi_2 = M_PI*2.0 - m_deltaPhi_flat/2.0;
		return (((r > m_Rinner && r < m_Router) || CheckPositionEqual(r, m_Rinner) || CheckPositionEqual(r, m_Router)) && (phi > phi_2 || phi < phi_1 || CheckAngleEqual(phi, phi_1) || CheckAngleEqual(phi, phi_2)));
	};

	/*
		For a given radius/phi are you inside of a given ring/wedge channel,
		or are you on the spacing between these channels
	*/
	inline bool IsRing(double r, int ringch) {
		double ringtop = m_Rinner + m_deltaR_flat_ring*(ringch + 1);
		double ringbottom = m_Rinner + m_deltaR_flat_ring*(ringch);
		return (r>ringbottom && r<ringtop); 
	};

	inline bool IsRingTopEdge(double r, int ringch) {
		double ringtop = m_Rinner + m_deltaR_flat_ring*(ringch + 1);
		return CheckPositionEqual(r, ringtop); 
	};

	inline bool IsRingBottomEdge(double r, int ringch) {
		double ringbottom = m_Rinner + m_deltaR_flat_ring*(ringch);
		return CheckPositionEqual(r, ringbottom); 
	};

	inline bool IsWedge(double phi, int wedgech) {
		double wedgetop = -m_deltaPhi_flat/2.0 + m_deltaPhi_flat_wedge*(wedgech+1);
		double wedgebottom = -m_deltaPhi_flat/2.0 + m_deltaPhi_flat_wedge*(wedgech);
		return ((phi>wedgebottom && phi<wedgetop));
	};

	inline bool IsWedgeTopEdge(double phi, int wedgech) {
		double wedgetop = -m_deltaPhi_flat/2.0 + m_deltaPhi_flat_wedge*(wedgech+1);
		return CheckAngleEqual(phi, wedgetop);
	}

	inline bool IsWedgeBottomEdge(double phi, int wedgech) {
		double wedgebottom = -m_deltaPhi_flat/2.0 + m_deltaPhi_flat_wedge*(wedgech);
		return CheckAngleEqual(phi, wedgebottom);
	}

	/*Class data*/
	double m_Router, m_Rinner, m_deltaPhi_flat, m_phiCentral, m_tilt;
	Mask::Vec3 m_translation;
	Mask::YRotation m_YRot;
	Mask::ZRotation m_ZRot;
	double m_deltaR_flat, m_deltaR_flat_ring, m_deltaPhi_flat_wedge;

	std::vector<std::vector<Mask::Vec3>> m_ringCoords_flat, m_wedgeCoords_flat;
	std::vector<std::vector<Mask::Vec3>> m_ringCoords_tilt, m_wedgeCoords_tilt;

};

#endif
