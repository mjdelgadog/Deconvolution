// -*- mode: c++; c-basic-offset: 2; -*-
#ifndef OPHITALG_DECO_H
#define OPHITALG_DECO_H
/*!
 * Title:   OpHit Algorithims
 * Author:  Ben Jones, MIT ( Edited by wketchum@lanl.gov, gleb.sinev@duke.edu
 *                           and kevin.wood@stonybrook.edu )
 *
 * Description:
 * These are the algorithms used by OpHit to produce optical hits.
 */

#include "larana/OpticalDetector/OpHitFinder/PMTPulseRecoBase.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"

#include <vector>

namespace calib {
  class IPhotonCalibrator;
}
namespace detinfo {
  class DetectorClocksData;
}
namespace geo {
  class GeometryCore;
}
namespace pmtana {
  class PulseRecoManager;
}

namespace opdet {

 void RunHitFinder_deco(std::vector<recob::OpWaveform> const&,                 
                    std::vector<recob::OpHit>&,
                    pmtana::PulseRecoManager const&,
                    pmtana::PMTPulseRecoBase const&,
                    geo::GeometryCore const&,
                    float,
                    detinfo::DetectorClocksData const&,
                    calib::IPhotonCalibrator const&,
                    bool use_start_time = false);

  void ConstructHit(float,
                    int,
                    double,
                    pmtana::pulse_param const&,
                    std::vector<recob::OpHit>&,
                    detinfo::DetectorClocksData const&,
                    calib::IPhotonCalibrator const&,
                    bool use_start_time = false);

} // End opdet namespace

#endif
