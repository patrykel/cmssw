from SimG4Core.Application.hectorParameter_cfi import *
import FWCore.ParameterSet.Config as cms
import copy
process = cms.Process("TestFlatGun")
process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring('warnings',
        'errors',
        'infos',
        'debugs'),
    categories = cms.untracked.vstring('ForwardSim',
        'TotemRP'),
    debugModules = cms.untracked.vstring('*'),
    errors = cms.untracked.PSet(
        threshold = cms.untracked.string('ERROR')
    ),
    warnings = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING')
    ),
    infos = cms.untracked.PSet(
        threshold = cms.untracked.string('INFO')
    ),
    debugs = cms.untracked.PSet(
        threshold = cms.untracked.string('DEBUG'),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        DEBUG = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        TotemRP = cms.untracked.PSet(
            limit = cms.untracked.int32(1000000)
        ),
        ForwardSim = cms.untracked.PSet(
            limit = cms.untracked.int32(1000000)
        )
    )
)

# Specify the maximum events to simulate
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(20)
)

# Configure the output module (save the result in a file)
process.o1 = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
    fileName = cms.untracked.string('file:test.root')
)
process.outpath = cms.EndPath(process.o1)

################## STEP 1 - process.generator
process.source = cms.Source("EmptySource")

# Use random number generator service
#process.load("Configuration.TotemCommon.RandomNumbers_cfi")

# particle generator paramteres
process.load("Configuration.StandardSequences.Services_cff")

process.load("Configuration.Generator.FlatLogKsiLogT_pp_13TeV_cfi")
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    g4SimHits = cms.PSet(initialSeed = cms.untracked.uint32(9876)),
    SimG4Object = cms.PSet(initialSeed =cms.untracked.uint32(9876)),
    RPSiDetDigitizer = cms.PSet(initialSeed =cms.untracked.uint32(137137)),
    sourceSeed = cms.PSet(initialSeed =cms.untracked.uint32(98765)),
    generator = cms.PSet(initialSeed = cms.untracked.uint32(98766)),
    SmearingGenerator = cms.PSet(initialSeed =cms.untracked.uint32(3849)),
    T2Digis = cms.PSet(initialSeed =cms.untracked.uint32(98765)),
    T2MCl = cms.PSet(initialSeed =cms.untracked.uint32(24141)),
    RPFastStationSimulation = cms.PSet(initialSeed =cms.untracked.uint32(12)),
    RPFastFullSimulation = cms.PSet(initialSeed =cms.untracked.uint32(13)),
    mix = cms.PSet(initialSeed = cms.untracked.uint32(24141)),
    LHCTransport = cms.PSet(initialSeed = cms.untracked.uint32(24143), engineName = cms.untracked.string('TRandom3')
  )

)
################# STEP 2 process.SmearingGenerator
process.SmearingGenerator = cms.EDProducer("GaussEvtVtxEnergyGenerator",
    src   = cms.InputTag("generator"),
    MeanX = cms.double(0.0),
    MeanY = cms.double(0.0),
    MeanZ = cms.double(0.0),
    SigmaX = cms.double(0.000001),
    SigmaY = cms.double(0.000001),
    SigmaZ = cms.double(0.000001),
    TimeOffset = cms.double(0.0)
)
# declare optics parameters
#process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_6500GeV_0p8_145urad_cfi")
process.BeamOpticsParamsESSource = cms.ESSource("BeamOpticsParamsESSource",
    BeamEnergy = cms.double(6500.0), # Gev
    ProtonMass = cms.double(0.938272029), # Gev
    LightSpeed = cms.double(300000000.0),
    NormalizedEmittanceX = cms.double(3.75e-06),
    NormalizedEmittanceY = cms.double(3.75e-06),
    BetaStarX = cms.double(0.8), # m
    BetaStarY = cms.double(0.8), # m
    CrossingAngleX = cms.double(145e-6),
    CrossingAngleY = cms.double(0.0),
    BeamDisplacementX = cms.double(0.0), # m
    BeamDisplacementY = cms.double(0.0), # m
    BeamDisplacementZ = cms.double(0.0), # m
    BunchSizeZ = cms.double(0.07), # m
    MeanXi = cms.double(0.0), # energy smearing
    SigmaXi = cms.double(0.0001)
)

process.ProtonTransportFunctionsESSource = cms.ESProducer("ProtonTransportFunctionsESSource",
    opticsFile = cms.string(''), # automatic
    maySymmetrize = cms.bool(True), # this optic is assymmetric
    verbosity = cms.untracked.uint32(10)
)

BeamProtTransportSetup = cms.PSet(
    Verbosity = cms.bool(True),
    ModelRootFile = cms.string('Geometry/VeryForwardProtonTransport/data/parametrization_6500GeV_90_transp_75.root'),
    Model_IP_150_R_Name = cms.string('ip5_to_beg_150_station_lhcb1'),
    Model_IP_150_L_Name = cms.string('ip5_to_beg_150_station_lhcb1'),

    # in m, should be consistent with geometry xml definitions
    Model_IP_150_R_Zmin = cms.double(0.0),
    Model_IP_150_R_Zmax = cms.double(202.769),
    Model_IP_150_L_Zmax = cms.double(-202.769),
    Model_IP_150_L_Zmin = cms.double(0.0),
)


# ################## STEP 3 process.g4SimHits
#
# # Geometry - beta* specific
# process.load("Geometry.VeryForwardGeometry.geometryRP_cfi")

# DDL geometry (ideal)
totemGeomXMLFiles = cms.vstring(
        'Geometry/CMSCommonData/data/materials.xml',
        'Geometry/CMSCommonData/data/rotations.xml',
        'Geometry/CMSCommonData/data/extend/cmsextent.xml',
        'Geometry/CMSCommonData/data/cms.xml',
        'Geometry/CMSCommonData/data/PhaseI/beampipe.xml',
        'Geometry/CMSCommonData/data/cmsBeam.xml',
        'Geometry/CMSCommonData/data/cmsMother.xml',
        'Geometry/CMSCommonData/data/mgnt.xml',
        'Geometry/ForwardCommonData/data/forward.xml',
        'Geometry/ForwardCommonData/data/totemRotations.xml',
        'Geometry/ForwardCommonData/data/totemMaterials.xml',
        'Geometry/ForwardCommonData/data/totemt1.xml',
        'Geometry/ForwardCommonData/data/totemt2.xml',
        'Geometry/ForwardCommonData/data/ionpump.xml',
        'Geometry/VeryForwardData/data/RP_Box.xml',
        'Geometry/VeryForwardData/data/RP_Box/RP_Box_000.xml',
        'Geometry/VeryForwardData/data/RP_Box/RP_Box_001.xml',
        'Geometry/VeryForwardData/data/RP_Box/RP_Box_002.xml',
        'Geometry/VeryForwardData/data/RP_Box/RP_Box_003.xml',
        'Geometry/VeryForwardData/data/RP_Box/RP_Box_004.xml',
        'Geometry/VeryForwardData/data/RP_Box/RP_Box_005.xml',
        'Geometry/VeryForwardData/data/RP_Box/RP_Box_020.xml',
        'Geometry/VeryForwardData/data/RP_Box/RP_Box_021.xml',
        'Geometry/VeryForwardData/data/RP_Box/RP_Box_022.xml',
        'Geometry/VeryForwardData/data/RP_Box/RP_Box_023.xml',
        'Geometry/VeryForwardData/data/RP_Box/RP_Box_024.xml',
        'Geometry/VeryForwardData/data/RP_Box/RP_Box_025.xml',
        'Geometry/VeryForwardData/data/RP_Box/RP_Box_100.xml',
        'Geometry/VeryForwardData/data/RP_Box/RP_Box_101.xml',
        'Geometry/VeryForwardData/data/RP_Box/RP_Box_102.xml',
        'Geometry/VeryForwardData/data/RP_Box/RP_Box_103.xml',
        'Geometry/VeryForwardData/data/RP_Box/RP_Box_104.xml',
        'Geometry/VeryForwardData/data/RP_Box/RP_Box_105.xml',
        'Geometry/VeryForwardData/data/RP_Box/RP_Box_120.xml',
        'Geometry/VeryForwardData/data/RP_Box/RP_Box_121.xml',
        'Geometry/VeryForwardData/data/RP_Box/RP_Box_122.xml',
        'Geometry/VeryForwardData/data/RP_Box/RP_Box_123.xml',
        'Geometry/VeryForwardData/data/RP_Box/RP_Box_124.xml',
        'Geometry/VeryForwardData/data/RP_Box/RP_Box_125.xml',
        'Geometry/VeryForwardData/data/RP_Hybrid.xml',
        'Geometry/VeryForwardData/data/RP_Materials.xml',
        'Geometry/VeryForwardData/data/RP_Transformations.xml',
        'Geometry/VeryForwardData/data/RP_Detectors_Assembly.xml',
        'Geometry/VeryForwardData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_000.xml',
        'Geometry/VeryForwardData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_001.xml',
        'Geometry/VeryForwardData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_002.xml',
        'Geometry/VeryForwardData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_003.xml',
        'Geometry/VeryForwardData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_004.xml',
        'Geometry/VeryForwardData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_005.xml',
        'Geometry/VeryForwardData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_020.xml',
        'Geometry/VeryForwardData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_021.xml',
        'Geometry/VeryForwardData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_022.xml',
        'Geometry/VeryForwardData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_023.xml',
        'Geometry/VeryForwardData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_024.xml',
        'Geometry/VeryForwardData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_025.xml',
        'Geometry/VeryForwardData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_100.xml',
        'Geometry/VeryForwardData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_101.xml',
        'Geometry/VeryForwardData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_102.xml',
        'Geometry/VeryForwardData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_103.xml',
        'Geometry/VeryForwardData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_104.xml',
        'Geometry/VeryForwardData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_105.xml',
        'Geometry/VeryForwardData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_120.xml',
        'Geometry/VeryForwardData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_121.xml',
        'Geometry/VeryForwardData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_122.xml',
        'Geometry/VeryForwardData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_123.xml',
        'Geometry/VeryForwardData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_124.xml',
        'Geometry/VeryForwardData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_125.xml',
        'Geometry/VeryForwardData/data/RP_Device.xml',
        'Geometry/VeryForwardData/data/RP_Vertical_Device.xml',
        'Geometry/VeryForwardData/data/RP_Horizontal_Device.xml',
        'Geometry/VeryForwardData/data/RP_220_Right_Station.xml',
        'Geometry/VeryForwardData/data/RP_220_Left_Station.xml',
        'Geometry/VeryForwardData/data/RP_147_Right_Station.xml',
        'Geometry/VeryForwardData/data/RP_147_Left_Station.xml',
        'Geometry/VeryForwardData/data/RP_Stations_Assembly.xml',
        'Geometry/VeryForwardData/data/RP_Sensitive_Dets.xml',
        'Geometry/VeryForwardData/data/RP_Cuts_Per_Region.xml',
        'Geometry/VeryForwardData/data/RP_Param_Beam_Region.xml')

process.XMLIdealGeometryESSource = cms.ESSource("XMLIdealGeometryESSource",
    geomXMLFiles = totemGeomXMLFiles,
    rootNodeName = cms.string('cms:CMSE')
)

# position of RPs
process.XMLIdealGeometryESSource.geomXMLFiles.append("Geometry/VeryForwardData/data/2016_ctpps_15sigma_margin0/RP_Dist_Beam_Cent.xml")

# extended geometries
process.TotemRPGeometryESModule = cms.ESProducer("TotemRPGeometryESModule",
    verbosity = cms.untracked.uint32(10)
)


#
# # TODO Change to the LowBetaSettings
# process.XMLIdealGeometryESSource_CTPPS.geomXMLFiles.append('Geometry/VeryForwardData/data/RP_Dist_Beam_Low_Betha/RP_Dist_Beam_Cent.xml')
#
# # misalignments
# # https://github.com/cms-sw/cmssw/blob/d40b17c22838d14956c0399f357148f05ecdb369/Geometry/VeryForwardGeometryBuilder/python/TotemRPIncludeAlignments_cfi.py
#process.load("Geometry.VeryForwardGeometryBuilder.TotemRPIncludeAlignments_cfi")
#process.TotemRPIncludeAlignments.MisalignedFiles = cms.vstring()
#process.TotemRPIncludeAlignments.RealFiles = cms.vstring(
#'Alignment/RPData/LHC/2015_10_18_fill4511/version2/sr/45.xml',
#'Alignment/RPData/LHC/2015_10_18_fill4511/version2/sr/56.xml'
#)
# # Magnetic Field, by default we have 3.8T
process.load("Configuration.StandardSequences.MagneticField_cff")
#
# # G4 simulation & proton transport
# # SimG4Core/Application/python/g4SimHits_cfi.py

process.load("SimG4Core.Application.g4SimHits_cfi")
process.g4SimHits.Physics.BeamProtTransportSetup = BeamProtTransportSetup
#process.g4SimHits.Generator.HepMCProductLabel = 'generator'    # The input source for G4 module is connected to "process.source".
process.g4SimHits.G4TrackingManagerVerbosity = cms.untracked.int32(10)
process.g4SimHits.OverrideUserStackingAction = cms.bool(True)   # HINT: TOTEM specific
process.g4SimHits.TransportParticlesThroughWholeBeampipe = cms.bool(True)
# TOTEM uses MeasuredGeometryRecord instead of IdealGeometryRecord
process.g4SimHits.UseMeasuredGeometryRecord = cms.untracked.bool(False)  # HINT: TOTEM specific
process.g4SimHits.SteppingVerbosity = cms.int32(5)
process.g4SimHits.G4EventManagerVerbosity = cms.untracked.int32(10)
process.g4SimHits.G4StackManagerVerbosity = cms.untracked.int32(0)
process.g4SimHits.Watchers = cms.VPSet(
#        cms.PSet( # HINT: TOTEM specific
#            type = cms.string('SimTracer'),
#            SimTracer = cms.PSet(verbose = cms.bool(True)),
#        ),
#       cms.PSet( # HINT: TOTEM specific
#           type = cms.string('TotemRP'),
#           TotemRP = cms.PSet(
#               Names = cms.vstring('TotemHitsRP'),
#               FileName = cms.string('TotemTestRP_Hits.root'),
#               RPDebugFileName = cms.string('TotemDebugRP.root'),
#               FileNameOLD = cms.string('TotemTestRP_Hits_Old.root'),
#               Verbosity = cms.bool(True)
#           )
#       )
    )
process.g4SimHits.HepMCProductLabel = cms.InputTag("generator")
process.g4SimHits.Physics.DefaultCutValue = cms.double(100.0)
process.g4SimHits.Generator = cms.PSet(
    HectorEtaCut,
    HepMCProductLabel=cms.string('SmearingGenerator'),
    ApplyPCuts=cms.bool(False),
    ApplyPtransCut=cms.bool(False),
    MinPCut=cms.double(0.04),  ## the cut is in GeV
    MaxPCut=cms.double(99999.0),  ## the pmax=99.TeV
    ApplyEtaCuts=cms.bool(False),
    MinEtaCut=cms.double(-5.5),
    MaxEtaCut=cms.double(5.5),
    RDecLenCut=cms.double(2.9),  ## (cm) the cut on vertex radius
    LDecLenCut=cms.double(30.0),  ## (cm) decay volume length
    ApplyPhiCuts=cms.bool(False),
    MinPhiCut=cms.double(-3.14159265359),  ## (radians)
    MaxPhiCut=cms.double(3.14159265359),  ## according to CMS conventions
    ApplyLumiMonitorCuts=cms.bool(False),  ## primary for lumi monitors
    Verbosity=cms.untracked.int32(0),
    LeaveScatteredProtons=cms.untracked.bool(True),
    ## HINT: TOTEM specific - Leave intact protons after scattering for further near beam transport
    LeaveOnlyScatteredProtons=cms.untracked.bool(False)
    ## HINT: TOTEM specific - Leave only intact protons and reject all the other particles
)

process.g4SimHits.StackingAction.SaveFirstLevelSecondary = cms.untracked.bool(True)
process.g4SimHits.Totem_RP_SD = cms.PSet(  # HINT: TOTEM specific
    Verbosity=cms.int32(0)
)
process.g4SimHits.CastorSD.nonCompensationFactor = cms.double(0.85)
process.g4SimHits.PPS_Timing_SD = cms.PSet(
        Verbosity = cms.int32(0)
    )
#
# # Use particle table
process.load("SimGeneral.HepPDTESSource.pdt_cfi")
#
process.g4SimHits.PPSSD = cms.PSet(
  Verbosity = cms.untracked.int32(10)
)
#
# ################## Step 3 - Magnetic field configuration
# # todo declare in standard way (not as hardcoded raw config)
#
process.magfield = cms.ESSource("XMLIdealGeometryESSource",
    geomXMLFiles = cms.vstring('Geometry/CMSCommonData/data/normal/cmsextent.xml',
        'Geometry/CMSCommonData/data/cms.xml',
        'Geometry/CMSCommonData/data/cmsMagneticField.xml',
        'MagneticField/GeomBuilder/data/MagneticFieldVolumes_1103l.xml',
        'Geometry/CMSCommonData/data/materials.xml'),
    rootNodeName = cms.string('cmsMagneticField:MAGF')
)
process.prefer("magfield")

process.ParametrizedMagneticFieldProducer = cms.ESProducer("ParametrizedMagneticFieldProducer",
    version = cms.string('OAE_1103l_071212'),
    parameters = cms.PSet(
        BValue = cms.string('3_8T')
    ),
    label = cms.untracked.string('parametrizedField')
)

process.VolumeBasedMagneticFieldESProducer = cms.ESProducer("VolumeBasedMagneticFieldESProducer",
    scalingVolumes = cms.vint32(14100, 14200, 17600, 17800, 17900,
        18100, 18300, 18400, 18600, 23100,
        23300, 23400, 23600, 23800, 23900,
        24100, 28600, 28800, 28900, 29100,
        29300, 29400, 29600, 28609, 28809,
        28909, 29109, 29309, 29409, 29609,
        28610, 28810, 28910, 29110, 29310,
        29410, 29610, 28611, 28811, 28911,
        29111, 29311, 29411, 29611),
    scalingFactors = cms.vdouble(1, 1, 0.994, 1.004, 1.004,
        1.005, 1.004, 1.004, 0.994, 0.965,
        0.958, 0.958, 0.953, 0.958, 0.958,
        0.965, 0.918, 0.924, 0.924, 0.906,
        0.924, 0.924, 0.918, 0.991, 0.998,
        0.998, 0.978, 0.998, 0.998, 0.991,
        0.991, 0.998, 0.998, 0.978, 0.998,
        0.998, 0.991, 0.991, 0.998, 0.998,
        0.978, 0.998, 0.998, 0.991),
    useParametrizedTrackerField = cms.bool(True),
    label = cms.untracked.string(''),
    version = cms.string('grid_1103l_090322_3_8t'),
    debugBuilder = cms.untracked.bool(False),
    paramLabel = cms.string('parametrizedField'),
    geometryVersion = cms.int32(90322),
    gridFiles = cms.VPSet(cms.PSet(
        path = cms.string('grid.[v].bin'),
        master = cms.int32(1),
        volumes = cms.string('1-312'),
        sectors = cms.string('0')
    ),
        cms.PSet(
            path = cms.string('S3/grid.[v].bin'),
            master = cms.int32(3),
            volumes = cms.string('176-186,231-241,286-296'),
            sectors = cms.string('3')
        ),
        cms.PSet(
            path = cms.string('S4/grid.[v].bin'),
            master = cms.int32(4),
            volumes = cms.string('176-186,231-241,286-296'),
            sectors = cms.string('4')
        ),
        cms.PSet(
            path = cms.string('S9/grid.[v].bin'),
            master = cms.int32(9),
            volumes = cms.string('14,15,20,21,24-27,32,33,40,41,48,49,56,57,62,63,70,71,286-296'),
            sectors = cms.string('9')
        ),
        cms.PSet(
            path = cms.string('S10/grid.[v].bin'),
            master = cms.int32(10),
            volumes = cms.string('14,15,20,21,24-27,32,33,40,41,48,49,56,57,62,63,70,71,286-296'),
            sectors = cms.string('10')
        ),
        cms.PSet(
            path = cms.string('S11/grid.[v].bin'),
            master = cms.int32(11),
            volumes = cms.string('14,15,20,21,24-27,32,33,40,41,48,49,56,57,62,63,70,71,286-296'),
            sectors = cms.string('11')
        )),
    cacheLastVolume = cms.untracked.bool(True)
)

# ################## STEP 4 mix pdt_cfi
#
from SimGeneral.MixingModule.mixObjects_cfi import *
from SimGeneral.MixingModule.pixelDigitizer_cfi import *
from SimGeneral.MixingModule.stripDigitizer_cfi import *
from SimGeneral.MixingModule.ecalDigitizer_cfi import *
from SimGeneral.MixingModule.hcalDigitizer_cfi import *
from SimGeneral.MixingModule.castorDigitizer_cfi import *
from SimGeneral.MixingModule.trackingTruthProducer_cfi import *

process.mix = cms.EDProducer("MixingModule",
   moduleLabel = cms.untracked.string("mix"),
   digitizers = cms.PSet(
#     pixel = cms.PSet(
#       pixelDigitizer
#     )
#                         ,
#       strip = cms.PSet(
#     stripDigitizer
#       )
#                           ,
#     ecal = cms.PSet(
#       ecalDigitizer
#     ),
# #     hcal = cms.PSet(
# #       hcalDigitizer
# #     ),
#     castor  = cms.PSet(
#       castorDigitizer
#     ),
#     mergedtruth = cms.PSet(
#         trackingParticles
#     )
    ),
    LabelPlayback = cms.string(''),
    maxBunch = cms.int32(3),
    minBunch = cms.int32(-5), ## in terms of 25 ns

    bunchspace = cms.int32(25),
    mixProdStep1 = cms.bool(False),
    mixProdStep2 = cms.bool(False),

    playback = cms.untracked.bool(False),
    useCurrentProcessOnly = cms.bool(False),
    mixObjects = cms.PSet(
        mixCH = cms.PSet(
            mixCaloHits
        ),
        mixTracks = cms.PSet(
            mixSimTracks
        )
                          ,
        mixVertices = cms.PSet(
            mixSimVertices
        ),
        mixSH = cms.PSet(
            mixSimHits
        ),
        mixHepMC = cms.PSet(
            mixHepMCProducts
        )
    )
)


#from SimGeneral/MixingModule/python/mix_Objects_cfi.py
process.mix.mixObjects.mixSH.input =  cms.VInputTag(  # note that this list needs to be in the same order as the subdets
        #cms.InputTag("g4SimHits","BSCHits"), cms.InputTag("g4SimHits","BCM1FHits"), cms.InputTag("g4SimHits","PLTHits"), cms.InputTag("g4SimHits","FP420SI"),
        cms.InputTag("g4SimHits","MuonCSCHits"), cms.InputTag("g4SimHits","MuonDTHits"), cms.InputTag("g4SimHits","MuonRPCHits"),
        cms.InputTag("g4SimHits","TotemHitsRP"), cms.InputTag("g4SimHits","PPSTrackerHits"),
        cms.InputTag("g4SimHits","TrackerHitsPixelBarrelHighTof"), cms.InputTag("g4SimHits","TrackerHitsPixelBarrelLowTof"),
        cms.InputTag("g4SimHits","TrackerHitsPixelEndcapHighTof"), cms.InputTag("g4SimHits","TrackerHitsPixelEndcapLowTof"),
	cms.InputTag("g4SimHits","TrackerHitsTECHighTof"), cms.InputTag("g4SimHits","TrackerHitsTECLowTof"), cms.InputTag("g4SimHits","TrackerHitsTIBHighTof"),
        cms.InputTag("g4SimHits","TrackerHitsTIBLowTof"), cms.InputTag("g4SimHits","TrackerHitsTIDHighTof"),
	cms.InputTag("g4SimHits","TrackerHitsTIDLowTof"), cms.InputTag("g4SimHits","TrackerHitsTOBHighTof"), cms.InputTag("g4SimHits","TrackerHitsTOBLowTof"))

process.mix.mixObjects.mixSH.subdets = cms.vstring(
       # 'BSCHits',
       # 'BCM1FHits',
       # 'PLTHits',
       # 'FP420SI',
        'MuonCSCHits',
        'MuonDTHits',
        'MuonRPCHits',
        'TotemHitsRP',
        'PPSTrackerHits',
        'TrackerHitsPixelBarrelHighTof',
        'TrackerHitsPixelBarrelLowTof',
        'TrackerHitsPixelEndcapHighTof',
        'TrackerHitsPixelEndcapLowTof',
        'TrackerHitsTECHighTof',
        'TrackerHitsTECLowTof',
        'TrackerHitsTIBHighTof',
        'TrackerHitsTIBLowTof',
        'TrackerHitsTIDHighTof',
        'TrackerHitsTIDLowTof',
        'TrackerHitsTOBHighTof',
        'TrackerHitsTOBLowTof')

process.mix.mixObjects.mixSH.crossingFrames = cms.untracked.vstring('MuonCSCHits',
'MuonDTHits',
'MuonRPCHits',
'TotemHitsRP',
'PPSTrackerHits')


# Use particle table
process.load("SimGeneral.HepPDTESSource.pdt_cfi")


 ################## STEP 5 RPDigiProducer

process.load("SimTotem.RPDigiProducer.RPSiDetConf_cfi")

################### STEP 6 reco
#
# process.load("Configuration.TotemStandardSequences.RP_Digi_and_TrackReconstruction_cfi")
#
# ################## STEP 7 TotemNtuplizer
#
# process.load("TotemAnalysis.TotemNtuplizer.TotemNtuplizer_cfi")
# process.TotemNtuplizer.outputFileName = "test.ntuple.root"
# process.TotemNtuplizer.RawEventLabel = 'source'
# process.TotemNtuplizer.RPReconstructedProtonCollectionLabel = cms.InputTag('RP220Reconst')
# process.TotemNtuplizer.RPReconstructedProtonPairCollectionLabel = cms.InputTag('RP220Reconst')
# process.TotemNtuplizer.RPMulFittedTrackCollectionLabel = cms.InputTag("RPMulTrackNonParallelCandCollFit")
# process.TotemNtuplizer.includeDigi = cms.bool(True)
# process.TotemNtuplizer.includePatterns = cms.bool(True)
#
# #######
process.load("RecoCTPPS.Configuration.recoCTPPS_cff")
process.totemRPClusterProducer.tagDigi = cms.InputTag("RPSiDetDigitizer")
process.dump = cms.EDAnalyzer("EventContentAnalyzer")
#
process.p1 = cms.Path(
	process.generator
#*process.VtxSmeared
	*process.SmearingGenerator
	*process.g4SimHits
       	*process.mix
	*process.RPSiDetDigitizer
        #*process.RPClustProd
        #*process.RPHecoHitProd
	# *process.RPSinglTrackCandFind
	# *process.RPSingleTrackCandCollFit
#	*process.RP220Reconst
	*process.recoCTPPS
)
