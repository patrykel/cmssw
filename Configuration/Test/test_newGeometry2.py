import FWCore.ParameterSet.Config as cms


# https://twiki.cern.ch/twiki/pub/TOTEM/CMSShow81/config.py.txt

process = cms.Process("DUMP")

import Geometry.VeryForwardGeometry.geometryRP_cfi as geometry_rp
process.XMLIdealGeometryESSource = geometry_rp.XMLIdealGeometryESSource_CTPPS
process.TotemRPGeometryESModule = geometry_rp.TotemRPGeometryESModule

process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/VeryForwardData/data/RP_Dist_Beam_Low_Betha/RP_Dist_Beam_Cent.xml')
process.XMLIdealGeometryESSource.rootNodeName = cms.string('cms:OCMS')

if 'Geometry/CMSCommonData/data/normal/cmsextent.xml' in process.XMLIdealGeometryESSource.geomXMLFiles:
    process.XMLIdealGeometryESSource.geomXMLFiles.remove('Geometry/CMSCommonData/data/normal/cmsextent.xml')
    process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/CMSCommonData/data/extend/cmsextent.xml')


process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger = cms.Service("MessageLogger",
                                    debugModules = cms.untracked.vstring('*'),
                                    destinations = cms.untracked.vstring('cout'),
                                    cout = cms.untracked.PSet
                                        (
                                        threshold = cms.untracked.string('DEBUG'),
                                        noLineBreaks = cms.untracked.bool(True)
                                        )
                                    )


process.source = cms.Source("EmptySource")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

process.load("SimG4Core.Application.g4SimHits_cfi")
process.g4SimHits.FileNameGDML = cms.untracked.string('totem_geom_evc2.gdml')


process.add_(cms.ESProducer("TGeoMgrFromDdd",
                            verbose = cms.untracked.bool(False),
                            level   = cms.untracked.int32(25)
                            ))


process.dump = cms.EDAnalyzer("DumpSimGeometry")


process.p = cms.Path(process.dump)
