/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors:
*   Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#ifndef TotemCondFormatsDataRecordProtonTransportRcd_h
#define TotemCondFormatsDataRecordProtonTransportRcd_h

#include "boost/mpl/vector.hpp"
#include "FWCore/Framework/interface/DependentRecordImplementation.h"
#include "TotemCondFormats/DataRecord/interface/BeamOpticsParamsRcd.h"
#include "Geometry/Records/interface/VeryForwardRealGeometryRecord.h"
 
/**
 *\brief EventSetup record for proton transport data.
 **/

class ProtonTransportRcd : public edm::eventsetup::DependentRecordImplementation<ProtonTransportRcd,
  boost::mpl::vector<BeamOpticsParamsRcd, VeryForwardRealGeometryRecord> >
{
};
 
#endif

