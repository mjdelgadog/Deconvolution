// =========================================================
// OpHitDeconvolution Module
// @authors     : Daniele Guffanti, Maritza Delgado
// @created     : Jan 26, 2022 
// Filter wiener-FFT
//========================================================= 
 
#ifndef OpHitDeconvolution_h
#define OpHitDeconvolution_h
 
// Framework includes

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ART extensions
#include "nurandom/RandomUtils/NuRandomService.h"

// LArSoft includes
#include "larana/OpticalDetector/OpHitFinder/AlgoCFD.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoFixedWindow.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoSiPM.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoSlidingWindow.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoThreshold.h"
#include "larana/OpticalDetector/OpHitFinder/OpHitAlg.h"
#include "larana/OpticalDetector/OpHitFinder/PMTPulseRecoBase.h"
#include "larana/OpticalDetector/OpHitFinder/PedAlgoEdges.h"
#include "larana/OpticalDetector/OpHitFinder/PedAlgoRollingMean.h"
#include "larana/OpticalDetector/OpHitFinder/PedAlgoUB.h"
#include "larana/OpticalDetector/OpHitFinder/PulseRecoManager.h"
#include "larcore/CoreUtils/ServiceUtil.h" 
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpWaveform.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/Simulation/BeamGateInfo.h"
#include "larreco/Calibrator/IPhotonCalibrator.h"
#include "larreco/Calibrator/IPhotonCalibratorService.h"
#include "larreco/Calibrator/PhotonCalibratorStandard.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"

 // CLHEP includes
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandFlat.h"

// C++ Includes
#include <map>
#include <memory>
#include <string>
#include <stdio.h>
#include <iostream>
#include <vector>

// ROOT includes
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "RooGaussian.h"
#include "TVirtualFFT.h"
#include "TStyle.h"
#include "TLine.h"
#include "TF1.h"
#include "TComplex.h"
#include "TFFTComplexReal.h"
#include <Rtypes.h> 

using std::string;
using std::vector;


namespace opdet {
    class OpHitDeconvolution : public art::EDProducer{
        public:
            explicit OpHitDeconvolution(const fhicl::ParameterSet&);
            virtual ~OpHitDeconvolution();
            void produce(art::Event& evt);
        
            // Parameters we'll read from the fcl-file
            std::string fInputModule;        // Module used to create OpDetWaveforms
            std::string fInstanceName;       // Input tag for OpDetWaveforms collection
            double fSampleFreq;              // Sampling frequency in MHz 
            double  fTimeBegin;              // Beginning of waveform in us
            double  fTimeEnd;                // End of waveform in us
            size_t fReadoutWindow;           // In ticks
            int    fBaselineSubtract;        // Baseline to subtract from each waveform
            short  fPedestal;                // In ADC counts
            double  fLineNoiseRMS;           // Pedestal RMS in ADC counts
            size_t fPreTrigger;              // In ticks
            std::vector<double> fSinglePEWaveform;   // Template for a single PE in ADC
            unsigned int WfDeco;
            std::string fDigiDataFile;
            double ScintYieldRatio;          // liquid argon scintillation yield ratio
            double ScintFastTimeConst;
            double ScintSlowTimeConst;
            double ScintResolutionScale;  
            double fScale;
            int fSamples;
            bool WienerFilter;
            bool GaussFilter;
            double fFrequencyCutOff;
            double fTickCutOff;
            //----------------------------------------------------
            // Declare member functions
            std::vector<raw::OpDetWaveform> RunOpHitDeconvolution(std::vector<raw::OpDetWaveform> const& wfHandle);
            
            //------------------------------------------------------
            //Load TFileService service
            art::ServiceHandle<art::TFileService> tfs;
      };
    }
#endif
  
namespace opdet{
    DEFINE_ART_MODULE(OpHitDeconvolution)
}

namespace opdet {
    //---------------------------------------------------------------------------
    // Constructor
    OpHitDeconvolution::OpHitDeconvolution(const fhicl::ParameterSet& pset)
    : EDProducer{pset}
    {  
        //read fhicl paramters

        fLineNoiseRMS      = pset.get< double  >("LineNoiseRMS" ); //noise for FFT
        fInputModule       = pset.get< std::string >("InputModule");
        fPreTrigger        = pset.get< size_t >("PreTrigger" );
        fTimeBegin         = pset.get< double >("TimeBegin");
        fTimeEnd           = pset.get< double >("TimeEnd"  );
        fPedestal          = pset.get< short  >("Pedestal"   );
        fReadoutWindow     = pset.get< size_t >("ReadoutWindow" );
        fBaselineSubtract  = pset.get< int >("fBaselineSubtract", 0 );
        fInputModule       = pset.get< std::string >("InputModule");
        fInstanceName      = pset.get< std::string >("InstanceName");
        fDigiDataFile      = pset.get< std::string >("DigiDataFile");
        fSamples           = pset.get< int >("fSamples");
        fScale             = pset.get< double >("Scale");
        WienerFilter       = pset.get< bool   >("WienerFilter");
        GaussFilter        = pset.get< bool   >("GaussFilter");
        fFrequencyCutOff   = pset.get< double >("GaussFilterCutOff");
        WfDeco=0;  
      
        // auto const *LarProp = lar::providerFrom<detinfo::LArPropertiesService>();
        art::ServiceHandle< art::TFileService > tfs; 

        // This module produces
        produces< std::vector< raw::OpDetWaveform > >();  
        produces< std::vector< recob::OpWaveform> > ();   


        // Obtain parameters from TimeService
        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
        fSampleFreq = clockData.OpticalClock().Frequency();
                 
        auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clockData);
       
        fTimeBegin = 0; 
        fTimeEnd   = detProp.ReadOutWindowSize() / clockData.TPCClock().Frequency();  
    }
    
    //---------------------------------------------------------------------------
    // Destructor
    OpHitDeconvolution::~OpHitDeconvolution(){
    }    
    
    //-------------------------------------------------------------------------
    void OpHitDeconvolution::produce(art::Event& evt){
    
        art::Handle< std::vector< raw::OpDetWaveform > > wfHandle;
        evt.getByLabel(fInputModule, fInstanceName, wfHandle);

        // Access ART's TFileService, which will handle creating and writing
        art::ServiceHandle< art::TFileService > tfs;

        //******************************
        //--- Read SPE ----
        //******************************
        std::ifstream SPEData;
        SPEData.open(fDigiDataFile);
        SPEData.is_open();   
        std::vector <double> SinglePEVec_x;   
        Double_t x = 0.0;         
        while (SPEData >> x ) { SinglePEVec_x.push_back(x); }
        fSinglePEWaveform = SinglePEVec_x;
        fSinglePEWaveform.resize(fSamples, 0);

        double*xh = new double[fSamples];
        
        for (int i=0; i< fSamples; i++) {
            xh[i] = fSinglePEWaveform[i];
        }
   
        //******************************
        //-- Read Waveform----
        //****************************** 
          
        std::vector<raw::OpDetWaveform> digi_wave = *wfHandle;
        int NOpDetWaveform = digi_wave.size();
        std::cout << NOpDetWaveform << std::endl;

        //pointer that will store produced Waveform
        auto out_wave = std::make_unique< std::vector< raw::OpDetWaveform > >();
        auto out_decowave = std::make_unique< std::vector< recob::OpWaveform > >();

        std::vector<short unsigned int > out_digiwave(NOpDetWaveform); //vector in which the waveform will be saved
        std::vector<float> out_recowave(NOpDetWaveform);               //vector in which the decowaveform will be saved, using float
        double*xv = new double[fSamples]; 
     
        for (auto const& wf: digi_wave) {
            
            // Resize wvfs
            if (static_cast<int>(wf.Waveform().size()) <= fSamples) { 
                out_recowave.resize(fSamples,0); 
            }

            else {
                printf("\nWARNING: waveform size is %lu, which is larger than fSamples (%i)\n", 
                wf.Waveform().size(), fSamples); 
                out_recowave.resize(fSamples); 
            }

            //Loop through the waveforms
            for (Int_t i= 0; i < fSamples; i++){
                if (i < static_cast<int>(wf.Waveform().size())) xv[i] = (wf[i]-fPedestal);
                else xv[i] = CLHEP::RandGauss::shoot(0, fLineNoiseRMS); 
            }

            
            //Found maximum peak in the Waveform
            Double_t SPE_Max = 0;
            double maxADC=*max_element(wf.begin(),wf.end());
            double maxAmplit= maxADC-fPedestal;
            SPE_Max = maxAmplit/8.0 ;
                       
            //******************************
            //--Original Signal---
            //******************************  
            std::vector<double>xs(fSamples,0.);
            //ST profile 
            auto const *Larprop = lar::providerFrom<detinfo::LArPropertiesService>();
            std::vector<double>SignalTime={(Larprop->ScintFastTimeConst()* 0.001),(Larprop->ScintSlowTimeConst()* 0.001)};
            std::vector<double>SignalScint={Larprop->ScintYieldRatio(),1.-Larprop->ScintYieldRatio()};

            //******************************
            //--Noise---
            //******************************
            double*xn = new double[fSamples]; 

            //******************************
            //--deconvolved signal---
            //******************************     
            double* xy;   
        
            //******************************
            //TIME
            //******************************
            double dt = 1/fSampleFreq; // in us       
            double t0 = fTimeBegin;      // in us
            double t = t0; 
            double* xt = new double[fSamples];      

            for (int i=0; i < fSamples; i++) {     
                xs[i] = SPE_Max*TMath::Gaus(i, 8, 0.1,true);//0.09//gauss function  
                //NEW Definition 
                double lightsignal=0;
                for (size_t j=0; j<SignalTime.size();j++){
                    lightsignal+=SignalScint[j]*exp(-t/SignalTime[j]);
                }
                xs[i] = lightsignal; 
                xt[i] = t;
                xn[i] = CLHEP::RandGauss::shoot(fPedestal, fLineNoiseRMS);//gRandom->Gaus(0, fLineNoiseRMS);//white noise  (0, hNoiseAmpl->GetRMS())
                t+=dt;
            }
            
            //******************************
            // Get FFT service.
            //******************************
            TVirtualFFT* fft = TVirtualFFT::FFT(1, &fSamples, "M R2C");

            // TVirtualFFT::SetTransform(0);
            
            // xH: FFT of spe response
            std::vector<TComplex> xH;   xH.resize(fSamples, TComplex(0, 0)); 
            double*xH_re = new double[fSamples];  double*xH_im =new double[fSamples];
            // xV: FFT of waveform
            std::vector<TComplex> xV;   xV.resize(fSamples, TComplex(0, 0));  
            double*xV_re = new double[fSamples];  double*xV_im =new double[fSamples];
            // xS: FFT of original signal 
            std::vector<TComplex> xS;   xS.resize(fSamples, TComplex(0, 0));  
            double*xS_re = new double[fSamples];  double*xS_im =new double[fSamples];  
            // xN: FFT of original noise 
            std::vector<TComplex> xN;   xN.resize(fSamples, TComplex(0, 0));  
            double*xN_re = new double[fSamples];  double*xN_im =new double[fSamples];
            // xY: FFT of the filtered signal
            std::vector<TComplex> xY;   xY.resize(fSamples, TComplex(0, 0));  
            double*xY_re = new double[fSamples];  double*xY_im =new double[fSamples] ;
            // G: FFT of the filtered signal
            std::vector<TComplex> G ;   G.resize(fSamples, TComplex(0, 0));
   
            //*******************************Despues de TComplex los vectores=0*************************
            // spe template FFT
            fft->SetPoints(xh);  
            fft->Transform();
            fft->GetPointsComplex(xH_re, xH_im);
            
            // waveform template FFT
            fft->SetPoints(xv);
            fft->Transform();
            fft->GetPointsComplex(xV_re, xV_im);

            // Original signal FFT 
            fft->SetPoints(&xs[0]);
            //fft->SetPoints(xs);
            fft->Transform();
            fft->GetPointsComplex(xS_re, xS_im);

            // Noise FFT 
            fft->SetPoints(xn);
            fft->Transform();
            fft->GetPointsComplex(xN_re, xN_im);

            //Spectral density
      
            double*H2 =new double[fSamples];   //spe response spectral density
            double*S2 =new double[fSamples];   //original signal spectral density
            double*N2 =new double[fSamples];   // noise spectral density 
      
            //******************************
            // Compute filters.
            //******************************                             
            for (int i=0; i<fSamples*0.5+1; i++) {
                // fill FFT arrays
                xH[i] = TComplex(xH_re[i], xH_im[i]);
                xV[i] = TComplex(xV_re[i], xV_im[i]);
                xS[i] = TComplex(xS_re[i], xS_im[i]);
                xN[i] = TComplex(xN_re[i], xN_im[i]);
                
                // Compute spectral density
                H2[i] = xH[i].Rho2();
                S2[i] = xS[i].Rho2();
                N2[i] = fLineNoiseRMS * fLineNoiseRMS * fSamples ;
                
                if (WienerFilter == true){
                    // Compute Wiener filter
                    G[i] = TComplex::Conjugate(xH[i])*S2[i] / (H2[i]*S2[i] + N2[i]);
                }

                else if (GaussFilter == true){
                    // Compute gauss filter
                    G[0] = TComplex(0,0);
                    fTickCutOff = fSamples*fFrequencyCutOff*1e6*fSampleFreq;
                    G[i] = TComplex::Exp(-0.5*TMath::Power(i*1e6*fSampleFreq/(fSamples*fTickCutOff),2))/xH[i];
                }

                else{
                    // Compute dec signal
                    G[i] = TComplex::Power(xH[i],-1);// Standard dec is just the division of signal and SPE template in Fourier space
                }

                // Correct template pretrigger
                TComplex phase = TComplex(0., -TMath::TwoPi()*i*fPreTrigger/(fSamples)); 
                G[i] = G[i]*TComplex::Exp(phase);
                
                // Compute dec signal
                xY[i] = (G[i]*xV[i]);
                xY_re[i] = xY[i].Re(); xY_im[i] = xY[i].Im();
            }

            //Transform of the filtered signal
            fft = TVirtualFFT::FFT(1, &fSamples, "M C2R");
            fft->SetPointsComplex(xY_re, xY_im);
            fft->Transform();
            xy = fft->GetPointsReal();
            
            // Correct baseline after OpHitDeconvolution
            double DecPedestal = 0;
            for (size_t i=0; i<fPreTrigger; i++){
                DecPedestal = DecPedestal + xy[i];
            }
            DecPedestal = DecPedestal/int(fPreTrigger);

            // Write dec wvf to output vector
            for (int i=0; i<fSamples; i++){                 
                out_recowave[i] = (xy[i]-DecPedestal)*fScale;
                // std::cout << xy[i] << std::endl;
            }
            
            raw::OpDetWaveform dgwave( wf.TimeStamp(), wf.ChannelNumber(), out_digiwave );
            out_wave->emplace_back(std::move(dgwave));
        
            recob::OpWaveform decwav(wf.TimeStamp(), wf.ChannelNumber(), out_recowave );
            out_decowave->emplace_back(std::move(decwav));   
        }//for

    // Push the OpDetWaveforms and OpWaveform into the event
    evt.put(std::move(out_wave));
    evt.put(std::move(out_decowave));
    WfDeco++;
    }//void
} // namespace opdet