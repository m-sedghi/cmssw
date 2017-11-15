#include "Geometry/GEMGeometryBuilder/src/ME0GeometryParsFromDD.h"
#include "DataFormats/MuonDetId/interface/ME0DetId.h"
#include "CondFormats/GeometryObjects/interface/RecoIdealGeometry.h"
#include "DetectorDescription/Core/interface/DDFilter.h"
#include "DetectorDescription/Core/interface/DDFilteredView.h"
#include "DetectorDescription/Core/interface/DDSolid.h"
#include "Geometry/MuonNumbering/interface/MuonDDDNumbering.h"
#include "Geometry/MuonNumbering/interface/MuonBaseNumber.h"
#include "Geometry/MuonNumbering/interface/ME0NumberingScheme.h"
#include "DataFormats/GeometryVector/interface/Basic3DVector.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"

void
ME0GeometryParsFromDD::build( const DDCompactView* cview, 
			      const MuonDDDConstants& muonConstants,
			      RecoIdealGeometry& rgeo )
{
  std::string attribute = "MuStructure";
  std::string value     = "MuonEndCapME0";

  // Asking only for the MuonME0's
  DDSpecificsMatchesValueFilter filter{DDValue(attribute, value, 0.0)};
  DDFilteredView fview(*cview,filter);

  this->buildGeometry( fview, muonConstants, rgeo );
}

void
ME0GeometryParsFromDD::buildGeometry(DDFilteredView& fv,
				     const MuonDDDConstants& muonConstants,
				     RecoIdealGeometry& rgeo)
{  
  LogDebug("ME0GeometryParsFromDD") <<"Building the geometry service";
  LogDebug("ME0GeometryParsFromDD") << "About to run through the ME0 structure\n" 
				    <<" First logical part "
				    <<fv.logicalPart().name().name(); 
  
  bool doChambers = fv.firstChild();
  LogDebug("ME0GeometryParsFromDD") << "doSuperChamber = " << doChambers;
  // loop over superchambers
  while (doChambers){

    // getting chamber id from eta partitions
    fv.firstChild();fv.firstChild();
    MuonDDDNumbering mdddnumCh(muonConstants);
    ME0NumberingScheme me0NumCh(muonConstants);
    int rawidCh = me0NumCh.baseNumberToUnitNumber(mdddnumCh.geoHistoryToBaseNumber(fv.geoHistory()));
    ME0DetId detIdCh = ME0DetId(rawidCh);
    // back to chambers
    fv.parent();fv.parent();

    buildChamber(fv, detIdCh, rgeo);

    // loop over chambers
    // only 1 chamber
    bool doLayers = fv.firstChild();
    while (doLayers){
      

      fv.firstChild();
      MuonDDDNumbering mdddnum(muonConstants);
      ME0NumberingScheme me0Num(muonConstants);
      int rawId = me0Num.baseNumberToUnitNumber(mdddnum.geoHistoryToBaseNumber(fv.geoHistory()));
      ME0DetId detId = ME0DetId(rawId);
      ME0DetId detIdLa = detId.layerId();
      fv.parent();

      // build layer
      buildLayer(fv, detIdLa, rgeo);
      
      // loop over ME0EtaPartitions
      bool doEtaPart = fv.firstChild();
      while (doEtaPart){

	MuonDDDNumbering mdddnum(muonConstants);
	ME0NumberingScheme me0Num(muonConstants);
	int rawid = me0Num.baseNumberToUnitNumber(mdddnum.geoHistoryToBaseNumber(fv.geoHistory()));
	ME0DetId detId = ME0DetId(rawid);

	buildEtaPartition(fv, detId, rgeo);
	
	doEtaPart = fv.nextSibling();
      }
      fv.parent();
      doLayers = fv.nextSibling();
    }
    fv.parent();
    doChambers = fv.nextSibling();
  }  
}

void
ME0GeometryParsFromDD::buildChamber(DDFilteredView& fv, ME0DetId detId, RecoIdealGeometry& rgeo)
{
  LogDebug("ME0GeometryParsFromDD") << "buildChamber "<<fv.logicalPart().name().name()
				    <<" "<< detId <<std::endl;
  
  DDBooleanSolid solid = (DDBooleanSolid)(fv.logicalPart().solid());
  std::vector<double> dpar = solid.parameters(); 
  
  double dy = dpar[0]/cm;//length is along local Y
  double dz = dpar[3]/cm;// thickness is long local Z
  double dx1= dpar[4]/cm;// bottom width is along local X
  double dx2= dpar[8]/cm;// top width is along local X
  //dpar = solid.solidB().parameters();
  //dz += dpar[3]/cm;// chamber thickness

  ME0DetId me0id = detId.chamberId();
  std::vector<std::string> strpars;
  std::string name = fv.logicalPart().name().name();
  strpars.emplace_back(name);
  
  std::vector<double> pars{dx1, dx2, dy, dz};
  std::vector<double> vtra = getTranslation(fv);
  std::vector<double> vrot = getRotation(fv);

  rgeo.insert(me0id.rawId(), vtra, vrot, pars, strpars);
}

void
ME0GeometryParsFromDD::buildLayer(DDFilteredView& fv, ME0DetId detId, RecoIdealGeometry& rgeo)
{
  LogDebug("ME0GeometryParsFromDD") << "buildLayer "<<fv.logicalPart().name().name()
				    <<" "<< detId <<std::endl;
  
  DDBooleanSolid solid = (DDBooleanSolid)(fv.logicalPart().solid());
  std::vector<double> dpar = solid.parameters(); 
  
  double dy = dpar[0]/cm;//length is along local Y
  double dz = dpar[3]/cm;// thickness is long local Z
  double dx1= dpar[4]/cm;// bottom width is along local X
  double dx2= dpar[8]/cm;// top width is along local X
  //dpar = solid.solidB().parameters();
  //dz += dpar[3]/cm;// chamber thickness

  ME0DetId me0id = detId.chamberId();
  std::vector<std::string> strpars;
  std::string name = fv.logicalPart().name().name();
  strpars.emplace_back(name);
  
  std::vector<double> pars{dx1, dx2, dy, dz};
  std::vector<double> vtra = getTranslation(fv);
  std::vector<double> vrot = getRotation(fv);

  rgeo.insert(me0id.rawId(), vtra, vrot, pars, strpars);
}

void
ME0GeometryParsFromDD::buildEtaPartition(DDFilteredView& fv, ME0DetId detId, RecoIdealGeometry& rgeo)
{
  LogDebug("ME0GeometryParsFromDD") << "buildEtaPartition "<<fv.logicalPart().name().name()
				    <<" "<< detId <<std::endl;
  
  // EtaPartition specific parameter (nstrips and npads) 
  DDValue numbOfStrips("nStrips");
  DDValue numbOfPads("nPads");
  std::vector<const DDsvalues_type* > specs(fv.specifics());
  std::vector<const DDsvalues_type* >::iterator is = specs.begin();
  double nStrips = 0., nPads = 0.;
  for (;is != specs.end(); is++){
    if (DDfetch( *is, numbOfStrips)) nStrips = numbOfStrips.doubles()[0];
    if (DDfetch( *is, numbOfPads))   nPads = numbOfPads.doubles()[0];
  }
  LogDebug("ME0GeometryParsFromDD") 
    << ((nStrips == 0. ) ? ("No nStrips found!!") : ("Number of strips: " + std::to_string(nStrips))); 
  LogDebug("ME0GeometryParsFromDD") 
    << ((nPads == 0. ) ? ("No nPads found!!") : ("Number of pads: " + std::to_string(nPads)));
  
  // EtaPartition specific parameter (size) 
  std::vector<double> dpar = fv.logicalPart().solid().parameters();

  double dy = dpar[0]/cm;//length is along local Y
  double dz = dpar[3]/cm;// thickness is long local Z
  double dx1= dpar[4]/cm;// bottom width is along local X
  double dx2= dpar[8]/cm;// top width is along local X

  std::vector<std::string> strpars;
  std::string name = fv.logicalPart().name().name();
  strpars.emplace_back(name);
  
  std::vector<double> pars{dx1, dx2, dy, dz, nStrips, nPads};
  std::vector<double> vtra = getTranslation(fv);
  std::vector<double> vrot = getRotation(fv);

  rgeo.insert(detId.rawId(), vtra, vrot, pars, strpars);
}

std::vector<double> ME0GeometryParsFromDD::getTranslation(DDFilteredView& fv)
{
  const DDTranslation& tran = fv.translation();
  return {tran.x(), tran.y(), tran.z()};
}

std::vector<double> ME0GeometryParsFromDD::getRotation(DDFilteredView& fv)
{
  const DDRotationMatrix& rota = fv.rotation();//.Inverse();
  DD3Vector x, y, z;
  rota.GetComponents(x,y,z);  
  return
    { x.X(), x.Y(), x.Z(),
      y.X(), y.Y(), y.Z(),
      z.X(), z.Y(), z.Z() };
}
