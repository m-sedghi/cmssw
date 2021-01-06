///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------

#include "Geometry/FbcmGeometryBuilder/src/FbcmGeometryBuilderFromDDD.h"
#include "Geometry/FbcmGeometry/interface/FbcmGeometry.h"
#include "Geometry/FbcmGeometry/interface/FbcmSiPadSpecs.h"

#include "DetectorDescription/Core/interface/DDFilter.h"
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "DetectorDescription/Core/interface/DDSolid.h"

#include "DataFormats/GeometrySurface/interface/RectangularPlaneBounds.h"
#include "DataFormats/GeometrySurface/interface/SimpleTubBounds.h"
#include "DataFormats/GeometryVector/interface/Basic3DVector.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include <algorithm>
#include <iostream>
#include <string>


#define DONTCARE 0

FbcmGeometryBuilderFromDDD::FbcmGeometryBuilderFromDDD() {}

FbcmGeometryBuilderFromDDD::~FbcmGeometryBuilderFromDDD() {}

void FbcmGeometryBuilderFromDDD::build(FbcmGeometry& theGeometry, const DDCompactView* cview) {
  std::string attribute = "FbcmDDDGeom" ; 
  std::string value = "FbcmTopLevelGeometry"; 
  
  //std::string attribute = "MuStructure" ; 
  //std::string value = "MuonEndCapME0"; 
  
//  std::cout << "attribute: " << attribute << ", Tag: " << value << "\n";
  
  DDSpecificsMatchesValueFilter filter{DDValue(attribute, value, 0.0)};
  // I should make sure that this filter leads to a valid fview
  // Otherwise and error shoould be rised
  
  
  //std::cout << "Hello from FbcmGeometryBuilderFromDDD_build\n";
  
  DDFilteredView fview(*cview, filter);
// before utilization of fview, it refers to the Main RootNode, i.e. OCMS, 
// after one move down (child), it refers to the Filter
  fview.firstChild(); // this is essential to apply the filter to the Node
    
  buildGeometry(theGeometry, fview);
}

void FbcmGeometryBuilderFromDDD::buildGeometry(FbcmGeometry& theFbcmGeometry, DDFilteredView& fv) {
  //FbcmGeometry* geometry = new FbcmGeometry();

  
   	DDValue NumOfStations("nStations");
	const std::vector<const DDsvalues_type*>& specs = fv.specifics();
	//std::cout << "SpecSizeSattion: " << specs.size() << "\n" ;
	//std::cout << specs << "\n" ;
	int nStations=0;
	for (auto const& is : specs) {
		if (DDfetch(is, NumOfStations)) 
			nStations = (int) (NumOfStations.doubles()[0]);
	}
	theFbcmGeometry.SetNumOfStations(nStations); 
	//std::cout << "SetNumOfStations read as: " << theFbcmGeometry.NumOfStations() << "\n";
	
	
	
//std::cout << "Hello from FbcmGeometryBuilderFromDDD_buildGeometry\n";

//std::cout << "fviewTopName logical Name: " <<  fv.logicalPart().name().name() << "\n" ;

/* 
 	std::cout << "fviewTopName: " <<  fv.name() << "\n" ;
	std::cout << "fviewTopName logical Name: " <<  fv.logicalPart().name().name() << "\n" ;
	fv.firstChild();
	std::cout << "fviewfirstChild: " <<  fv.logicalPart().name().name() << "\n" ;
		fv.firstChild();
	std::cout << "fviewfirstChild: " <<  fv.logicalPart().name().name() << "\n" ;
		fv.firstChild();
	std::cout << "fviewfirstChild: " <<  fv.logicalPart().name().name() << "\n" ;
		fv.firstChild();
	std::cout << "fviewfirstChild: " <<  fv.logicalPart().name().name() << "\n" ;
		fv.firstChild();
	std::cout << "fviewfirstChild: " <<  fv.logicalPart().name().name() << "\n" ;
			fv.nextSibling();
	std::cout << "fviewa after next: " <<  fv.logicalPart().name().name() << "\n" ;
				fv.firstChild();
	std::cout << "fviewfirstChild: " <<  fv.logicalPart().name().name() << "\n" ;
	  */
	

  LogTrace("FbcmGeometryBuilderFromDDD") << "Building the FBCM geometry:";
  LogTrace("FbcmGeometryBuilderFromDDD") << "Top level logical part: " << fv.logicalPart().name().name();

 /// --------- FBCM Geometry Builder --------------
   /// Caution: there is a aditinal level in fbcm.xml which is ignored in the the FBCMGeometry Class. 
  /// (i.e. SiliconCols). 
  // For the moment, in the fbxm.xml we have the following Parent-Child relationship: 
  /// FBCM --> Station --> SiliconDie --> SensorRow ---> SiPad
  // but the FbcmGeometry class is something like this: 
  /// FBCM --> Station --> SiliconDie -----------------> SiPad

  
  bool doFbcmSides = true;
  int SideID =1; // This is a temporary solving the code!!. It should be read out from the FBCM Side CopyNo. 
while (doFbcmSides) {
   //int FbcmCopyNo = fv.copyno(); // FbcmCopyNo is always 1, we should get the CopyNo of its parent.
   //fv.parent(); // its parent is the volum at endcap
    //= fv.copyno(); // the CopyNo of the Sides, i.e. the parent of FBCM
   //fv.firstChild(); // returing to FBCM
   
   //std::cout << "SideID: " << SideID << "\n"; 
   
   
  bool doStations = fv.firstChild(); // move to Stations level
  while (doStations) {
	int StationCopyNo = fv.copyno(); // Stations
	
	//std::cout << "Station: " <<  fv.logicalPart().name().name() << ",  ";
	//std::cout << StationCopyNo << "\n";
	
	FbcmDetId StationDetId_ = FbcmDetId(SideID,StationCopyNo,DONTCARE,DONTCARE);
	FbcmStationGeom* NewStation = buildStation(fv, StationDetId_.StationDetId());
    theFbcmGeometry.add(NewStation);
	
	//std::cout << "nRings: " << NewStation->NumOfRings() << ",nDiePerRing: " << NewStation->NumOfDiesPerRing() << " \n";


    // in FbcmGeometry: loop over SiliconDies of the Station
    bool doSiDies = fv.firstChild(); // SiliconDie
    while (doSiDies) {
      int SilconDieCopyNo = fv.copyno(); // SiliconDie
      	  
		//std::cout << "SiliconDie: " <<  fv.logicalPart().name().name() << ",  ";
		//std::cout << SilconDieCopyNo << "\n";
		
		FbcmDetId SiDieDetId = FbcmDetId(SideID,StationCopyNo,SilconDieCopyNo,DONTCARE);				
		FbcmSiliconDieGeom* NewSiDie = buildSiliconDie(fv, SiDieDetId.SiliconDieDetId());
		NewStation->add(NewSiDie);
		theFbcmGeometry.add(NewSiDie); 
		//unsigned int nCols=NewSiDie->NumOfCols();
		unsigned int nRows=NewSiDie->NumOfRows();
		
		//std::cout << SiDieDetId;
		//std::cout << "nCols: " << NewSiDie->NumOfCols() << ",nRows: " << NewSiDie->NumOfRows() << " \n";
	
	// loop over SiliconCols of each SiliconDie
	// SiliconCol has not a geometry in the FBCMGeometry !
    bool doSiCols = fv.firstChild(); // SiliconCols
    while (doSiCols) {
		int SilconColCopyNo = fv.copyno(); // SiliconCol
     
	  
	  // loop over SiPads of the SiliconCol
      bool doSiPads = fv.firstChild(); // SiPads
      while (doSiPads) {
        
       int SiPadCopyNo = fv.copyno(); // SiPad
	   	
		int SensorPadID=SilconColCopyNo*nRows+SiPadCopyNo; 
				
		//std::cout << "SiPad: " <<  fv.logicalPart().name().name() << ",  ";
		//std::cout << SensorPadID << "\n";
				
        FbcmDetId SiPadDetId = FbcmDetId(SideID,StationCopyNo,SilconDieCopyNo,SensorPadID);
		//std::cout << SiPadDetId << "\n"; 
        
		FbcmSiPadGeom* NewSiPad = buildSiPad(fv, SiPadDetId);
        NewSiDie->add(NewSiPad);
        theFbcmGeometry.add(NewSiPad); 

        doSiPads = fv.nextSibling();
      }
	  
	  
	  fv.parent(); // get back to SiliconCol
      doSiCols = fv.nextSibling(); // Next SiliconCol
	  
	}
      fv.parent(); // get back to SiliconDie
      doSiDies = fv.nextSibling();
    }
    fv.parent(); // get back to Station
    doStations = fv.nextSibling(); // Next Station
  }
    fv.parent(); // get back to FbcmSides
    doFbcmSides = fv.nextSibling(); // Other Side
	
	SideID++; // Notice: this line is just for temporary fix!, it should be read out from FilterView
}
  //return theFbcmGeometry;
}

FbcmStationGeom* FbcmGeometryBuilderFromDDD::buildStation(DDFilteredView& fv, FbcmDetId StationDetId_) const {
  LogTrace("FbcmGeometryBuilderFromDDD") << "buildStation " << fv.logicalPart().name().name() << ", StationDetId_: " << StationDetId_ << std::endl;
  
  //std:: cout << "buildStation " << fv.logicalPart().name().name() << ", StationDetId_: " << StationDetId_ << std::endl;

	//std::cout << " Hi buildStation: "<< StationDetId_;

  DDTubs  solid = (DDTubs)(fv.logicalPart().solid());
  
  std::vector<double> dpar = solid.parameters();

	double halfZ  = solid.zhalf() / cm;  // halfThickness 
	double rMin = solid.rIn() / cm; 
    double rMax = solid.rOut() / cm ;  
    double StartPhi = solid.startPhi();
    double DeltaPhi = solid.deltaPhi();


    // std::cout << " name of logical part = " << fv.logicalPart().name().name() << std::endl;
    // std::cout << " dpar is vector with size = " << dpar.size() << std::endl;
    // for (unsigned int i = 0; i < dpar.size(); ++i) {
      // std::cout << " dpar [" << i << "] = " << dpar[i] << " mm or rad" << std::endl;
	// }
	// std::cout << "halfZ: " << halfZ << "cm, rMin: " << rMin << "cm,  rMax: " << rMax << "cm, StartPhi: " << StartPhi << ", DeltaPhi: " << DeltaPhi << std::endl;
	// std::cout << "\n";

// #ifdef EDM_ML_DEBUG
  // LogTrace("FbcmGeometryBuilderFromDDD") << " name of logical part = " << fv.logicalPart().name().name() << std::endl;
  // LogTrace("FbcmGeometryBuilderFromDDD") << " dpar is vector with size = " << dpar.size() << std::endl;
  // for (unsigned int i = 0; i < dpar.size(); ++i) {
    // LogTrace("FbcmGeometryBuilderFromDDD") << " dpar [" << i << "] = " << dpar[i] / cm << " cm " << std::endl;
  // }
  // LogTrace("FbcmGeometryBuilderFromDDD") << "size  halfWidth: " << w << "cm, halfLength: " << h << "cm,  halfThickness: " << t << "cm" << std::endl;
// #endif

	DDValue NumOfDiesPerRing("nSiDiesPerStationPerRing");
	DDValue NumOfRings("nRingsPerStation");
	const std::vector<const DDsvalues_type*>& specs = fv.specifics();
	//std::cout << specs << "\n" ;
	int nDiesPerRing = 0, nRings = 0;
	for (auto const& is : specs) {
		if (DDfetch(is, NumOfDiesPerRing)) 
			nDiesPerRing = (int) (NumOfDiesPerRing.doubles()[0]);
		if (DDfetch(is, NumOfRings))
			nRings = (int)(NumOfRings.doubles()[0]);
	}



  // SimpleTubBounds(float rmin, float rmax, float dz, float startPhi, float deltaPhi);  
  FbcmBoundPlane surf(boundPlane(fv, new SimpleTubBounds(rMin, rMax , halfZ, StartPhi, DeltaPhi)));
  FbcmStationGeom* Station = new FbcmStationGeom(StationDetId_, surf, nDiesPerRing, nRings);
  
  return Station;
}

FbcmSiliconDieGeom* FbcmGeometryBuilderFromDDD::buildSiliconDie(DDFilteredView& fv, FbcmDetId SiDieDetId) const {
  LogTrace("FbcmGeometryBuilderFromDDD") << "buildSiliconDie " << fv.logicalPart().name().name() << ", SiliconDieId: " << SiDieDetId << std::endl;

	//std::cout << "buildSiliconDie: "<< SiDieDetId;

  DDBox solid = (DDBox)(fv.logicalPart().solid());
  std::vector<double> dpar = solid.parameters();
  
    // double w = dpar[0] / cm;  //  halfWidth
  // double h = dpar[1] / cm;  // halfLength
  // double t = dpar[2] / cm;  // halfThickness 
  
    double w = solid.halfX() / cm;  //  halfWidth
    double h = solid.halfY() / cm;  // halfLength
    double t = solid.halfZ() / cm;  // halfThickness 
  
  
  	DDValue NumOfCols("nCols");
	DDValue NumOfRows("nRows");
	const std::vector<const DDsvalues_type*>& specs = fv.specifics();
	//std::cout << specs << "\n" ;
	unsigned int nCols = 0, nRows = 0;
	for (auto const& is : specs) {
		if (DDfetch(is, NumOfCols))
			nCols = (unsigned int)(NumOfCols.doubles()[0]);
		if (DDfetch(is, NumOfRows)) 
			nRows = (unsigned int)(NumOfRows.doubles()[0]);
	}
  
  
  // std::cout << " name of logical part = " << fv.logicalPart().name().name() << std::endl;
  // std::cout << " dpar is vector with size = " << dpar.size() << std::endl;
  // for (unsigned int i = 0; i < dpar.size(); ++i) {
    // std::cout << " dpar [" << i << "] = " << dpar[i] / cm << " cm " << std::endl;
  // }
  // std::cout << "size  halfWidth: " << w << "cm, halfLength: " << h << "cm,  halfThickness: " << t << "cm" << std::endl;
  // std::cout << "\n";
  
  
#ifdef EDM_ML_DEBUG
  LogTrace("FbcmGeometryBuilderFromDDD") << " name of logical part = " << fv.logicalPart().name().name() << std::endl;
  LogTrace("FbcmGeometryBuilderFromDDD") << " dpar is vector with size = " << dpar.size() << std::endl;
  for (unsigned int i = 0; i < dpar.size(); ++i) {
    LogTrace("FbcmGeometryBuilderFromDDD") << " dpar [" << i << "] = " << dpar[i] / cm << " cm " << std::endl;
  }
  LogTrace("FbcmGeometryBuilderFromDDD") << "size  halfWidth: " << w << "cm, halfLength: " << h << "cm,  halfThickness: " << t << "cm" << std::endl;
#endif

  FbcmBoundPlane surf(boundPlane(fv, new RectangularPlaneBounds(w, h, t)));
  FbcmSiliconDieGeom* SiliconDie = new FbcmSiliconDieGeom(SiDieDetId, surf, nCols, nRows);
  return SiliconDie;
}

FbcmSiPadGeom* FbcmGeometryBuilderFromDDD::buildSiPad(DDFilteredView& fv, FbcmDetId detId) const {
  LogTrace("FbcmGeometryBuilderFromDDD") << "buildSiPad " << fv.logicalPart().name().name() << ", FbcmDetId: " << detId << std::endl;
  
  //std::cout << "buildSiPad: "<< detId;
  
  DDBox solid = (DDBox)(fv.logicalPart().solid());
  std::vector<double> dpar = solid.parameters();
  
    // double w = dpar[0] / cm;  //  halfWidth
  // double h = dpar[1] / cm;  // halfLength
  // double t = dpar[2] / cm;  // halfThickness 
  
    double w = solid.halfX() / cm;  //  halfWidth
    double h = solid.halfY() / cm;  // halfLength
    double t = solid.halfZ() / cm;  // halfThickness 
  
    // std::cout << " name of logical part = " << fv.logicalPart().name().name() << std::endl;
  // std::cout << " dpar is vector with size = " << dpar.size() << std::endl;
   // for (unsigned int i = 0; i < dpar.size(); ++i) {
   // std::cout << " dpar [" << i << "] = " << dpar[i] / cm << " cm " << std::endl;
   // }
   // std::cout << "size  halfWidth: " << w << "cm, halfLength: " << h << "cm,  halfThickness: " << t << "cm" << std::endl;
   // std::cout << "\n";
  
#ifdef EDM_ML_DEBUG
  LogTrace("FbcmGeometryBuilderFromDDD") << " name of logical part = " << fv.logicalPart().name().name() << std::endl;
  LogTrace("FbcmGeometryBuilderFromDDD") << " dpar is vector with size = " << dpar.size() << std::endl;
  for (unsigned int i = 0; i < dpar.size(); ++i) {
    LogTrace("FbcmGeometryBuilderFromDDD") << " dpar [" << i << "] = " << dpar[i] / cm << " cm " << std::endl;
  }
  LogTrace("FbcmGeometryBuilderFromDDD") << "size  halfWidth: " << w << "cm, halfLength: " << h << "cm,  halfThickness: " << t << "cm" << std::endl;
#endif

  std::vector<float> pars;
  pars.emplace_back(2*w); // pitchX = pars[0];
  pars.emplace_back(2*h); // pitchY = pars[1];

  //RectangularPlaneBounds::RectangularPlaneBounds(float w, float h, float t), halfWidth(w), halfLength(h), halfThickness(t) 
  FbcmBoundPlane surf(boundPlane(fv, new RectangularPlaneBounds(w, h, t)));
  
  std::string name = fv.logicalPart().name().name();
  FbcmSiPadSpecs* SiPad_Specs = new FbcmSiPadSpecs(GeomDetEnumerators::FBCM, name, pars);

  FbcmSiPadGeom* SiPad = new FbcmSiPadGeom(detId, surf, SiPad_Specs);
  return SiPad;
}

FbcmGeometryBuilderFromDDD::FbcmBoundPlane FbcmGeometryBuilderFromDDD::boundPlane(const DDFilteredView& fv, Bounds* bounds) const {
  // extract the position
  const DDTranslation& trans(fv.translation());
  const Surface::PositionType posResult(float(trans.x() / cm), float(trans.y() / cm), float(trans.z() / cm));
  
  //std::cout << " name of logical part = " << fv.logicalPart().name().name() << std::endl;
//  std::cout << "x: " << posResult.x() << ", y: " << posResult.y() << ", z: " << posResult.z() << "\n";



  // now the rotation
  //  DDRotationMatrix tmp = fv.rotation();
  // === DDD uses 'active' rotations - see CLHEP user guide ===
  //     ORCA uses 'passive' rotation.
  //     'active' and 'passive' rotations are inverse to each other
  //  DDRotationMatrix tmp = fv.rotation();
  
  const DDRotationMatrix& rotation = fv.rotation();  
  DD3Vector x, y, z;
  rotation.GetComponents(x, y, z);
  
  
    
  // std::cout << "translation: "<< fv.translation() << std::endl;
  // std::cout << "rotation   : "<< fv.rotation() << std::endl;
  
  // std::cout << "INVERSE rotation manually: \n"
         // << x.X() << ", " << x.Y() << ", " << x.Z() << std::endl
         // << y.X() << ", " << y.Y() << ", " << y.Z() << std::endl
         // << z.X() << ", " << z.Y() << ", " << z.Z() << std::endl;
		 
  
  // LogTrace("GEMGeometryBuilderFromDDD") << "translation: "<< fv.translation() << std::endl;
  // LogTrace("GEMGeometryBuilderFromDDD") << "rotation   : "<< fv.rotation() << std::endl;
  // LogTrace("GEMGeometryBuilderFromDDD") << "INVERSE rotation manually: \n"
  //        << x.X() << ", " << x.Y() << ", " << x.Z() << std::endl
  //        << y.X() << ", " << y.Y() << ", " << y.Z() << std::endl
  //        << z.X() << ", " << z.Y() << ", " << z.Z() << std::endl;

  Surface::RotationType rotResult(float(x.X()),
                                  float(x.Y()),
                                  float(x.Z()),
                                  float(y.X()),
                                  float(y.Y()),
                                  float(y.Z()),
                                  float(z.X()),
                                  float(z.Y()),
                                  float(z.Z()));

//for shape other than BOXs, the conversion of local GEANT position
// to global needs some modifications to change axis. 

// Im not sure about the rotation Axis conversion. I should double check it later. !!

//--------------------------
// it seems there was also a need to chagnge XZ to XY rotation axis
// for release earlier than 11 for Dbox as well. something like this:
  // if (!theCompatiblity11Flag) {
        // if (tran[2] > -1500.) {
          // Basic3DVector<float> newX(-1., 0., 0.);
          // Basic3DVector<float> newY(0., -1., 0.);
          // Basic3DVector<float> newZ(0., 0., 1.);
          // rot.rotateAxes(newX, newY, newZ);
        // }
      // }

//--------------------------
   // // Change of axes for the forward for other shapse other than DBox
   // Basic3DVector<float> newX(1., 0., 0.);
   // Basic3DVector<float> newY(0., 0., 1.);
   // Basic3DVector<float> newZ(0., 1., 0.);
   // newY *= -1;
   // rotResult.rotateAxes(newX, newY, newZ);
//--------------------------


  //std::cout << "rotation final :\n"<< rotResult << std::endl;

  return FbcmBoundPlane(new BoundPlane(posResult, rotResult, bounds));
}
