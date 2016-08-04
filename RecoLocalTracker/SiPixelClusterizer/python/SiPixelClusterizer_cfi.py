
import FWCore.ParameterSet.Config as cms

#
from CondTools.SiPixel.SiPixelGainCalibrationService_cfi import *
siPixelClusters = cms.EDProducer("SiPixelClusterProducer",
    SiPixelGainCalibrationServiceParameters,
    src = cms.InputTag("siPixelDigis"),
    ChannelThreshold = cms.int32(1000),
    MissCalibrate = cms.untracked.bool(True),
    SplitClusters = cms.bool(False),
    VCaltoElectronGain = cms.int32(65),
    VCaltoElectronOffset = cms.int32(-414),                          
    # **************************************
    # ****  payLoadType Options         ****
    # ****  HLT - column granularity    ****
    # ****  Offline - gain:col/ped:pix  ****
    # **************************************
    payloadType = cms.string('Offline'),
    SeedThreshold = cms.int32(1000),
    ClusterThreshold = cms.double(4000.0),
    # **************************************
    maxNumberOfClusters = cms.int32(-1), # -1 means no limit.
)


# This customization will be removed once we have phase1 pixel digis
from Configuration.StandardSequences.Eras import eras
eras.phase1Pixel.toModify(siPixelClusters, # FIXME
    src = 'simSiPixelDigis',
    MissCalibrate = False
)

# Need these until phase2 pixel templates are used
eras.phase2_tracker.toModify(siPixelClusters, # FIXME
    src = cms.InputTag('simSiPixelDigis', "Pixel")
)
