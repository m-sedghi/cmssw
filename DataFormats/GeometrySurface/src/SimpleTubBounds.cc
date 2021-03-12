

#include "DataFormats/GeometrySurface/interface/SimpleTubBounds.h"
#include "DataFormats/GeometrySurface/interface/LocalError.h"
#include <algorithm>
#include <cmath>

SimpleTubBounds::SimpleTubBounds(float rmin, float rmax, float dz, float startPhi, float deltaPhi)
    : theRmin(rmin), theRmax(rmax), halfThickness(dz), theStartPhi(startPhi),theDeltaPhi(deltaPhi) {
	theZmin=-halfThickness;
	theZmax=halfThickness,
	theEndPhi=theStartPhi+theDeltaPhi;
  if (theRmin > theRmax)
    std::swap(theRmin, theRmax);
}

bool SimpleTubBounds::inside(const Local2DPoint& p, const LocalError& err) const { return Bounds::inside(p, err); }

Bounds* SimpleTubBounds::clone() const { return new SimpleTubBounds(*this); }

bool SimpleTubBounds::inside(const Local3DPoint& p, const LocalError& err, float scale) const {
	
	///Notice: Considering theStartPhi and theEndPhi for checking with localerro has not been implemented in this function.!!
	// I should update it. 
	
  if (p.z() < theZmin || p.z() > theZmax)
    return false;  // check the easy part first

  double perp2 = p.perp2();
  double perp = p.perp(); // sqrt(perp2);
  if (perp2 == 0)
    return scale * scale * (err.xx() + err.xy()) > theRmin * theRmin;

  // rotated error along p.x(),p.y()
  // equivalent to (but faster than) err.rotate(p.x(),p.y()).xx()
  // since we don't need all matrix elements
  float deltaR = scale * sqrt(p.x() * p.x() / perp2 * err.xx() - 2 * p.x() * p.y() / perp2 * err.xy() +
                              p.y() * p.y() / perp2 * err.yy());
  return perp > std::max(theRmin - deltaR, 0.f) && perp < theRmax + deltaR; // No Phi dependent !!!!
}
