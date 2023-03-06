// Gleb Sinev, Duke, 2016
//
// This module finds periods of time-localized activity
// on each channel of the optical system and creates OpHits.
//
// Split from OpFlashFinder_module.cc
// by Ben Jones, MIT, 2013
//

#ifndef OpHitFinderSPE_h
#define OpHitFinderSPE_h

// Framework includes
 
// LArSoft includes
#include "larana/OpticalDetector/OpHitFinder/AlgoCFD.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoFixedWindow.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoSiPM.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoSlidingWindow.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoThreshold.h"
#include "OpHitAlg_deco.h"
#include "larana/OpticalDetector/OpHitFinder/PMTPulseRecoBase.h"
#include "larana/OpticalDetector/OpHitFinder/PedAlgoEdges.h"
#include "larana/OpticalDetector/OpHitFinder/PedAlgoRollingMean.h"
#include "larana/OpticalDetector/OpHitFinder/PedAlgoUB.h"
#include "larana/OpticalDetector/OpHitFinder/PulseRecoManager.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/Simulation/BeamGateInfo.h"
#include "larreco/Calibrator/IPhotonCalibrator.h"
#include "larreco/Calibrator/IPhotonCalibratorService.h"
#include "larreco/Calibrator/PhotonCalibratorStandard.h"
 
// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/Exception.h"
#include "fhiclcpp/ParameterSet.h"
 
// ROOT includes
 
// C++ Includes
#include <map>
#include <memory>
#include <string>
#include "TTree.h"
#include <TSpectrum.h>
#include <algorithm>

 namespace opdet {
  
    class OpHitFinderSPE : public art::EDProducer {
    public:
      // Standard constructor and destructor for an ART module.
      explicit OpHitFinderSPE(const fhicl::ParameterSet&);
      virtual ~OpHitFinderSPE();
    /*  art::InputTag fRawOpWaveformLabel;
      art::InputTag fWienerOpWaveformLabel;
      art::InputTag fImpulseOpWaveformLabel;  */ 

      // The producer routine, called once per event.
      void produce(art::Event&);
   
    private:
      std::map<int, int> GetChannelMap();
      std::vector<double> GetSPEScales();
      std::vector<double> GetSPEShifts();
    
      // The parameters we'll read from the .fcl file.
      std::string fInputModule; // Input tag for OpDetWaveform collection
      std::string fGenModule;
      std::vector<std::string> fInputLabels;
      std::set<unsigned int> fChannelMasks;
    
      pmtana::PulseRecoManager fPulseRecoMgr;
      pmtana::PMTPulseRecoBase* fThreshAlg;
      pmtana::PMTPedestalBase* fPedAlg;
    
      Float_t fHitThreshold;
      unsigned int fMaxOpChannel;
      bool fUseStartTime;
    
      calib::IPhotonCalibrator const* fCalib = nullptr;
      //fRawOpWaveformLabel(p.get<art::InputTag>("RawOpWaveformLabel")),
      //fWienerOpWaveformLabel(p.get<art::InputTag>("WienerOpWaveformLabel")),
      //fImpulseOpWaveformLabel(p.get<art::InputTag>("ImpulseOpWaveformLabel"))

      };
    
    }
    
  #endif
    
  namespace opdet {
   DEFINE_ART_MODULE(OpHitFinderSPE)
  }
 
  namespace opdet {
 
   //----------------------------------------------------------------------------
   // Constructor
   OpHitFinderSPE::OpHitFinderSPE(const fhicl::ParameterSet& pset) : EDProducer{pset}, fPulseRecoMgr()
   {
    // Indicate that the Input Module comes from .fcl
       fInputModule = pset.get<std::string>("InputModule");
       fGenModule = pset.get<std::string>("GenModule");
       fInputLabels = pset.get<std::vector<std::string>>("InputLabels");
       fUseStartTime = pset.get<bool>("UseStartTime", false);
    
       for (auto const& ch :
             pset.get<std::vector<unsigned int>>("ChannelMasks", std::vector<unsigned int>()))
          fChannelMasks.insert(ch);
   
       fHitThreshold = pset.get<float>("HitThreshold");
       bool useCalibrator = pset.get<bool>("UseCalibrator", false);
   
       auto const& geometry(*lar::providerFrom<geo::Geometry>());
       fMaxOpChannel = geometry.MaxOpChannel();
   
       if (useCalibrator) {
         // If useCalibrator, get it from ART
         fCalib = lar::providerFrom<calib::IPhotonCalibratorService>();
       }
       else {
         // If not useCalibrator, make an internal one based
         // on fhicl settings to hit finder.
         bool areaToPE = pset.get<bool>("AreaToPE");
         float SPEArea = pset.get<float>("SPEArea");
         float SPEShift = pset.get<float>("SPEShift", 0.);
   
         // Reproduce behavior from GetSPEScales()
         if (!areaToPE) SPEArea = 20;
   
         // Delete and replace if we are reconfiguring
         if (fCalib) { delete fCalib; }
   
         fCalib = new calib::PhotonCalibratorStandard(SPEArea, SPEShift, areaToPE);
       }

       // Initialize the rise time calculator tool
       auto const rise_alg_pset = pset.get_if_present<fhicl::ParameterSet>("RiseTimeCalculator");
   
       // Initialize the hit finder algorithm
       auto const hit_alg_pset = pset.get<fhicl::ParameterSet>("HitAlgoPset");
       std::string threshAlgName = hit_alg_pset.get<std::string>("Name");
       if (threshAlgName == "Threshold")
         fThreshAlg = new pmtana::AlgoThreshold(hit_alg_pset);
       else if (threshAlgName == "SiPM")
         fThreshAlg = new pmtana::AlgoSiPM(hit_alg_pset);
       else if (threshAlgName == "SlidingWindow")
         fThreshAlg = new pmtana::AlgoSlidingWindow(hit_alg_pset);
       else if (threshAlgName == "FixedWindow")
         fThreshAlg = new pmtana::AlgoFixedWindow(hit_alg_pset);
       else if (threshAlgName == "CFD")
         fThreshAlg = new pmtana::AlgoCFD(hit_alg_pset);
       else
         throw art::Exception(art::errors::UnimplementedFeature)
           << "Cannot find implementation for " << threshAlgName << " algorithm.\n";
       
       // Initialize the pedestal estimation algorithm
       auto const ped_alg_pset = pset.get<fhicl::ParameterSet>("PedAlgoPset");
       std::string pedAlgName = ped_alg_pset.get<std::string>("Name");
       if (pedAlgName == "Edges")
         fPedAlg = new pmtana::PedAlgoEdges(ped_alg_pset);
       else if (pedAlgName == "RollingMean")
         fPedAlg = new pmtana::PedAlgoRollingMean(ped_alg_pset);
       else if (pedAlgName == "UB")
         fPedAlg = new pmtana::PedAlgoUB(ped_alg_pset);
       else
         throw art::Exception(art::errors::UnimplementedFeature)
           << "Cannot find implementation for " << pedAlgName << " algorithm.\n";
   
       produces<std::vector<recob::OpHit>>();
       fPulseRecoMgr.AddRecoAlgo(fThreshAlg);
       fPulseRecoMgr.SetDefaultPedAlgo(fPedAlg);
       
     }
   
     //----------------------------------------------------------------------------
     // Destructor
     OpHitFinderSPE::~OpHitFinderSPE()
     {
   
       delete fThreshAlg;
       delete fPedAlg;
     }
   
     //----------------------------------------------------------------------------
     void OpHitFinderSPE::produce(art::Event& evt)
     {
       std:: cout << "OPHITFINDER RECO: " << std::endl; 
       
       art::Handle<std::vector<recob::OpWaveform>> wienerHandle;
       evt.getByLabel(fInputModule, fInputLabels.front(),wienerHandle);
       
      std:: cout << "Number of waveforms reco: " << wienerHandle->size() << std::endl;

       // These is the storage pointer we will put in the event
       std::unique_ptr<std::vector<recob::OpHit>> HitPtr(new std::vector<recob::OpHit>);

       //New OpHit vector
      // std::unique_ptr<std::vector<recob::OpHit>> HitPtr_D(new std::vector<recob::OpHit>);
      std:: cout << "OPHITFINDER RECO1: " << std::endl;
   
       std::vector<const sim::BeamGateInfo*> beamGateArray;
       try {
         evt.getView(fGenModule, beamGateArray);
       }
       catch (art::Exception const& err) {
         if (err.categoryCode() != art::errors::ProductNotFound) throw;
       }
     
       auto const& geometry(*lar::providerFrom<geo::Geometry>());
       auto const clock_data =
         art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
       auto const& calibrator(*fCalib);
       
        //Get the pulses from the event
        // Load pulses into DecoVector
        std:: cout << "OPHITFINDER RECO2.1: " << std::endl;
       if (fChannelMasks.empty() && fInputLabels.size() < 2) {
         art::Handle<std::vector<recob::OpWaveform>> wienerHandle;
         if (fInputLabels.empty())
           evt.getByLabel(fInputModule, wienerHandle);
         else
           evt.getByLabel(fInputModule, fInputLabels.front(),wienerHandle);
         assert(wienerHandle.isValid());
         RunHitFinder_deco(//*wfHandle,
                      *wienerHandle,
                      *HitPtr,
                      fPulseRecoMgr,
                      *fThreshAlg,
                      geometry,
                      fHitThreshold,
                      clock_data,
                      calibrator,
                      fUseStartTime);
       }
      else {
         // Reserve a large enough array/reco
         int totalsize_dec = 0;
         for (auto label : fInputLabels) {
           art::Handle<std::vector<recob::OpWaveform>> wienerHandle;
           evt.getByLabel(fInputModule, label, wienerHandle);
           if (!wienerHandle.isValid()) continue; // Skip non-existent collections
           totalsize_dec += wienerHandle->size();
          }

         std::vector<recob::OpWaveform> DecoWaveformVector;  //recob
         DecoWaveformVector.reserve(totalsize_dec);
                  
         for (auto label : fInputLabels) {    
           art::Handle<std::vector<recob::OpWaveform>> wienerHandle;
           evt.getByLabel(fInputModule, label, wienerHandle);
           if (!wienerHandle.isValid()) continue; // Skip non-existent collections
           
           for (auto const& dwf : *wienerHandle) {   //reco
             if (fChannelMasks.find(dwf.Channel()) != fChannelMasks.end()) continue;
             DecoWaveformVector.push_back(dwf);
           }
         }
   
         RunHitFinder_deco(//WaveformVector,
                      DecoWaveformVector,      
                      *HitPtr,
                      fPulseRecoMgr,
                      *fThreshAlg,
                      geometry,
                      fHitThreshold,
                      clock_data,
                      calibrator,
                      fUseStartTime);  //incluir risetimecalculator
          }
       
       evt.put(std::move(HitPtr));
    }
       
 } // namespace opdet
