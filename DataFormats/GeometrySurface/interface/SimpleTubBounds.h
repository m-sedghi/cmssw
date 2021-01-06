#ifndef Geom_SimpleTubBounds_H
#define Geom_SimpleTubBounds_H

#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometrySurface/interface/Bounds.h"

/** \class SimpleTubBounds
 * Plane bounds that define a Tub with a concentric hole in the middle. 
 * Starts anf angle phi1 and ends with phi2
 */
 

class SimpleTubBounds final : public Bounds {
public:
  /// Construct the bounds from min and max R and Z in LOCAL coordinates.
  // dz equals to HalfThickness
  SimpleTubBounds(float rmin, float rmax, float dz, float startPhi, float deltaPhi);

  float length() const override { return 2 * halfThickness; }
  float width() const override { return 2 * theRmax; }
  float thickness() const override { return theZmax - theZmin; } // = 2 * halfThickness

  bool inside(const Local3DPoint& p) const override {
    return (    (p.z() > theZmin) && (p.z() < theZmax) &&
				(p.perp2() > theRmin * theRmin) && (p.perp2() < theRmax * theRmax) &&
				(p.barePhi() > theStartPhi) && (p.barePhi() < theEndPhi)
		   );
  }

  using Bounds::inside;

  bool inside(const Local3DPoint& p, const LocalError& err, float scale) const override;

  virtual bool inside(const Local2DPoint& p, const LocalError& err) const;

  Bounds* clone() const override;

  /// Extension of the Bounds interface
  float innerRadius() const { return theRmin; }
  float outerRadius() const { return theRmax; }

  float minZ() const { return theZmin; }
  float maxZ() const { return theZmax; }
  
  float StartPhi() const { return theStartPhi; }
  float EndPhi() const { return theEndPhi; }

private:
  float theRmin;
  float theRmax;
  float theZmin;
  float theZmax;
  float halfThickness;
  float theStartPhi;
  float theDeltaPhi;
  float theEndPhi;
  
};

#endif  // Geom_SimpleTubBounds_H
