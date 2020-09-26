#ifndef DataFormats_FbcmDigi_SiPadAmplitude_H
#define DataFormats_FbcmDigi_SiPadAmplitude_H
///-------------------------------------------
//  Author: Mohammad Sedghi, msedghi@cern.ch
//  Isfahan University of Technology
//  Date created: September 2020
///-------------------------------------------

#include <cstdint>
#include <utility>
#include <cassert>

class SiPadAmplitude {
public:
 explicit SiPadAmplitude(
    unsigned int Side,
  unsigned int Station,
  unsigned int SiliconDie,
  unsigned int SiPad,
  float Ampl_A,
  float Ampl_SP,
  float x,
  float y,
  float sigma_x,
  float sigma_y,
  float sig_ToT,
  float Time,
  float ToF,
  float energyLoss,
  unsigned int rawId,
  unsigned int trackId) {      
			  Side_=Side;
			  Station_=Station;
			  SiliconDie_=SiliconDie;
			  SiPad_=SiPad;
			  Ampl_A_=Ampl_A;
			  Ampl_SP_=Ampl_SP;
			  x_=x;
			  y_=y;
			  sigma_x_=sigma_x;
			  sigma_y_=sigma_y;
			  sig_ToT_=sig_ToT;
			  time_=Time;
			  ToF_=ToF;
			  energyLoss_=energyLoss;
			  rawId_=rawId;
			  trackId_=trackId;
  };

  SiPadAmplitude() {      
			  Side_=0;
			  Station_=0;
			  SiliconDie_=0;
			  SiPad_=0;
			  Ampl_A_=0;
			  Ampl_SP_=0;
			  x_=0;
			  y_=0;
			  sigma_x_=0;
			  sigma_y_=0;
			  sig_ToT_=0;
			  time_=0;
			  ToF_=0;
			  energyLoss_=0;
			  rawId_=0;
			  trackId_=0;
  };


  ~SiPadAmplitude(){};
  
  unsigned int SideId() const { return   Side_  ; };
  unsigned int StationId()  const { return   Station_  ; };
  unsigned int SiliconDieId()  const { return  SiliconDie_   ; };
  unsigned int SiPadId()  const { return    SiPad_ ; };
  float AmplitudeA()  const { return    Ampl_A_ ; };
  float AmplitudeSP()  const { return    Ampl_SP_ ; };
  float xPos()  const { return    x_ ; };
  float yPos()  const { return   y_  ; };
  float Sigma_x()  const { return   sigma_x_  ; };
  float Sigma_y()  const { return   sigma_y_  ; };
  float Sig_ToT()  const { return  sig_ToT_   ; };
  float time()  const { return   time_  ; };
  float TimeOfFlight()  const { return  ToF_   ; };
  float EnergyLoss()  const { return   energyLoss_  ; };
  unsigned int RawId()  const { return    rawId_ ; };
  unsigned int TrackId()  const { return    trackId_ ; };
    

  inline bool operator<(const SiPadAmplitude& other) const { return AmplitudeA() < other.AmplitudeA(); }

private:
  unsigned int Side_;
  unsigned int Station_;
  unsigned int SiliconDie_;
  unsigned int SiPad_;
  float Ampl_A_;
  float Ampl_SP_;
  float x_;
  float y_;
  float sigma_x_;
  float sigma_y_;
  float sig_ToT_;
  float time_;
  float ToF_;
  float energyLoss_;
  unsigned int rawId_;
  unsigned int trackId_;
};

#include <iostream>
inline std::ostream& operator<<(std::ostream& o, const SiPadAmplitude& SiPadAmpl) {
  return o 
  << "SideId:" << SiPadAmpl.SideId() << ", " 
  << "StationId:" << SiPadAmpl.StationId() << ", " 
  << "SiliconDieId:" << SiPadAmpl.SiliconDieId() << ", " 
  << "SiPadId:" << SiPadAmpl.SiPadId() << "\n" 
  << "AmplA:" << SiPadAmpl.AmplitudeA() << ", " 
  << "AmplSP:" << SiPadAmpl.AmplitudeSP() << ", " 
  << "Sig_ToT:" << SiPadAmpl.Sig_ToT() << ", " 
  << "time:" << SiPadAmpl.time() << ", " 
  << "ToF:" << SiPadAmpl.TimeOfFlight() ;
}

#endif