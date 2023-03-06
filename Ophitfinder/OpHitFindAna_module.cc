// This analyzer writes out a TTree containing the properties of
// each reconstructed ophit
//

#ifndef OpHitFindAna_h
#define OpHitFindAna_h

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/OpDetPulse.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpWaveform.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// ART includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes
#include "TH1.h"
#include "THStack.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"

// C++ Includes
#include <map>
#include <vector>
#include <iostream>
#include <cstring>
#include <sstream>
#include "math.h"
#include <climits>


namespace opdet {

  class OpHitFindAna : public art::EDAnalyzer {
  public:
    // Standard constructor and destructor for an ART module.
    OpHitFindAna(const fhicl::ParameterSet&);
    virtual ~OpHitFindAna();
    // The analyzer routine, called once per event.
    void analyze(const art::Event&);

  private:
    // The stuff below is the part you'll most likely have to change to
    // go from this custom example to your own task.

    // The parameters we'll read from the .fcl file.
    std::string fInputModule; // Input tag for OpDet collection
    std::string fInstanceName;             // Input tag for OpDet collection
    double fSampleFreq;               // in MHz
    float fTimeBegin;                // in us
    float fTimeEnd;                  // in us

    // Flags to enable or disable output of debugging TH1 / TH2s
    //bool fMakeHistPerChannel;
    //bool fMakeHistAllChannels;
    //bool fMakeHitTree;

    // Output TTree and its branch variables
    TTree* fHitTree;
    Int_t fEventID;
    Int_t fOpChannel;
    Float_t fPeakTime;
    Float_t fNPe;
    
    double GetDriftWindow(detinfo::DetectorPropertiesData const& detProp) const;
    
  };

}

#endif

namespace opdet {

  DEFINE_ART_MODULE(OpHitFindAna)

}


namespace opdet {

  //-----------------------------------------------------------------------
  // Constructor
  OpHitFindAna::OpHitFindAna(fhicl::ParameterSet const& pset) : EDAnalyzer(pset)
  {

    // Indicate that the Input Module comes from .fcl
    fInputModule  =   pset.get< std::string >("InputModule");
    fInstanceName =   pset.get< std::string >("InstanceName");
    // Obtain parameters from DetectorClocksService
    
    // Obtaining parameters from the DetectorClocksService
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clockData);
    fSampleFreq = clockData.OpticalClock().Frequency();
    fTimeEnd   = detProp.ReadOutWindowSize() / clockData.TPCClock().Frequency();
    // Assume the readout is symmetrical around 0
    fTimeBegin = -1.*fTimeEnd;
    
    //fMakeHistPerChannel = pset.get<bool>("MakeHistPerChannel");
    //fMakeHistAllChannels = pset.get<bool>("MakeHistAllChannels");
    //fMakeHitTree = pset.get<bool>("MakeHitTree");

    // If required, make TTree for output

    //if (fMakeHitTree) {
      art::ServiceHandle<art::TFileService const> tfs;
      fHitTree = tfs->make<TTree>("HitTree", "HitTree");
      fHitTree->Branch("EventID", &fEventID, "EventID/I");
      fHitTree->Branch("OpChannel", &fOpChannel, "OpChannel/I");
      fHitTree->Branch("PeakTime", &fPeakTime, "PeakTime/F");
      fHitTree->Branch("NPe", &fNPe, "NPe/F");
     //}
   //  art::ServiceHandle< art::TFileService > tfs;
  }
 
  //-----------------------------------------------------------------------
  // Destructor
  OpHitFindAna::~OpHitFindAna() {  
  }    
  
  //-----------------------------------------------------------------------
  void OpHitFindAna::analyze(const art::Event& evt)
  {

    // Create a handle for our vector of pulses
    art::Handle<std::vector<recob::OpHit>> HitHandle;

    
    // Read in HitHandle
    evt.getByLabel(fInputModule, HitHandle);

    // Access ART's TFileService, which will handle creating and writing
    // histograms for us.
    art::ServiceHandle<art::TFileService const> tfs;

    art::ServiceHandle<geo::Geometry const> geom;
    int NOpChannels = geom->NOpChannels();
    
    
    for (int i = 0; i != NOpChannels; ++i) {
      //space for create histogram
    }

    fEventID = evt.id().event();

    // For every OpHit in the vector
    for (unsigned int i = 0; i < HitHandle->size(); ++i) {
      // Get OpHit
      art::Ptr<recob::OpHit> TheHitPtr(HitHandle, i);
      recob::OpHit TheHit = *TheHitPtr;

      fOpChannel = TheHit.OpChannel();
      fNPe = TheHit.PE();
      fPeakTime = TheHit.PeakTime();

      fHitTree->Fill();

     }
  }

} // namespace opdet

