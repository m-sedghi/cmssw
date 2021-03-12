import FWCore.ParameterSet.Config as cms

# def getSiPadDigitizer(process):
    # if hasattr(process,'mix') and hasattr(process.mix,'digitizers') and hasattr(process.mix.digitizers,'SiPad'):
        # return process.mix.digitizers.SiPad
    # return None
def no_aging(process):
	return process

def _1000invfb(process):
	#SiPadDigi=getSiPadDigitizer(process)
	if hasattr(process.mix.digitizers,'SiPad'):
		process.mix.digitizers.SiPad.SiPadSimParam.ReadoutNoiseInElec = cms.double(808.0)
		process.mix.digitizers.SiPad.SiPadSimParam.chargeCollectionEfficiency = cms.double(0.5534)
		for SiPadFE in process.mix.digitizers.SiPad.SiPadFrontEndParam :
			SiPadFE.SensorCapPerCm2  = cms.double(91.61)
			SiPadFE.SensorCouplingCapacitance = cms.double(314.10)
	return process
	
def _3000invfb(process):
	#SiPadDigi=getSiPadDigitizer(process)
	if hasattr(process.mix.digitizers,'SiPad'):
		process.mix.digitizers.SiPad.SiPadSimParam.ReadoutNoiseInElec = cms.double(1409.0)
		process.mix.digitizers.SiPad.SiPadSimParam.chargeCollectionEfficiency = cms.double(0.4060)
		for SiPadFE in process.mix.digitizers.SiPad.SiPadFrontEndParam :
			SiPadFE.SensorCapPerCm2  = cms.double(115.69)
			SiPadFE.SensorCouplingCapacitance = cms.double(312.29)
	return process

def _4000invfb(process):
	#SiPadDigi=getSiPadDigitizer(process)
	if hasattr(process.mix.digitizers,'SiPad'):
		process.mix.digitizers.SiPad.SiPadSimParam.ReadoutNoiseInElec = cms.double(1628.0)
		process.mix.digitizers.SiPad.SiPadSimParam.chargeCollectionEfficiency = cms.double(0.3744)
		for SiPadFE in process.mix.digitizers.SiPad.SiPadFrontEndParam :
			SiPadFE.SensorCapPerCm2  = cms.double(126.02)
			SiPadFE.SensorCouplingCapacitance = cms.double(311.38)
	return process
