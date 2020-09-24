#ifndef Tracker_FbcmSD_H
#define Tracker_FbcmSD_H
#include "SimG4CMS/Forward/interface/TimingSD.h"
#include "DataFormats/FbcmDetId/interface/FbcmDetId.h"

#include <string>

class SimTrackManager;
class G4Step;

class FbcmSD : public TimingSD {
public:
  FbcmSD(const std::string &,
          const edm::EventSetup &,
          const SensitiveDetectorCatalog &,
          edm::ParameterSet const &,
          const SimTrackManager *);
  ~FbcmSD() override;

  uint32_t setDetUnitId(const G4Step *) override;

protected:
  bool checkHit(const G4Step *, BscG4Hit *) override;

private:
  float energyCut;
  float energyHistoryCut;
};

#endif
