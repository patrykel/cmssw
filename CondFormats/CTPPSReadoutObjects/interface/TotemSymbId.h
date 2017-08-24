/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*   Jan Kašpar (jan.kaspar@cern.ch)
*
****************************************************************************/

#ifndef CondFormats_CTPPSReadoutObjects_TotemSymbId
#define CondFormats_CTPPSReadoutObjects_TotemSymbId

#include <iostream>

/**
 *\brief Symbolic ID describing an entity of a TOTEM subdetector.
 **/
class TotemSymbID
{
  public:
    /// identifies the TOTEM subsystem
    enum {RP, T1, T2} subSystem;

    /// chip ID, raw integer representation of DetId class
    unsigned int symbolicID;

    bool operator < (const TotemSymbID &sid) const
    {
	  return (symbolicID < sid.symbolicID);
    }

    bool operator == (const TotemSymbID &sid) const
    {
      return (symbolicID == sid.symbolicID);
    }
    
    friend std::ostream& operator << (std::ostream& s, const TotemSymbID &sid);
};

#endif
