import FWCore.ParameterSet.Config as cms
from RecoCTPPS.TotemRPLocal.totemRPLocalReconstruction_cff import *
# from RecoCTPPS.TotemRPLocal.ctppsDiamondLocalReconstruction_cff import *
# from RecoCTPPS.TotemRPLocal.ctppsLocalTrackLiteProducer_cfi import ctppsLocalTrackLiteProducer

recoCTPPS = cms.Sequence(
    totemRPLocalReconstruction
    # totemRPLocalReconstruction *
    # ctppsDiamondLocalReconstruction removed for now, as the simulation does not support it yet
    # ctppsDiamondLocalReconstruction *
    # ctppsLocalTrackLiteProducer
)
