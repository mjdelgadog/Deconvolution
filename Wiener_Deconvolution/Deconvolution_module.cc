// =========================================================
// Deconvolution Module
// @authors     : Daniele Guffanti, Maritza Delgado
// @created     : Jan 26, 2022 
// Filter wiener-FFT
//========================================================= 
 
#ifndef Deconvolution_h
#define Deconvolution_h
 
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
  class Deconvolution : public art::EDProducer{
    public:
      struct Config {
        public:
          fhicl::Atom<std::string> module_type{ fhicl::Name("module_type") }; 
          fhicl::Atom<std::string> InputModule{ fhicl::Name("InputModule") }; 
          fhicl::Atom<std::string> InstanceName{ fhicl::Name("InstanceName") }; 
          fhicl::Atom<Double_t>    LineNoiseRMS{ fhicl::Name("LineNoiseRMS"), 1.0 };
          fhicl::Atom<Double_t>    TimeBegin{ fhicl::Name("TimeBegin"), 0.}; 
          fhicl::Atom<Double_t>    TimeEnd{ fhicl::Name("TimeEnd"), 16.0}; 
          fhicl::Atom<size_t>      PreTrigger{ fhicl::Name("PreTrigger"), 0};
          fhicl::Atom<size_t>      ReadoutWindow{ fhicl::Name("ReadoutWindow"), 1000}; 
          fhicl::Atom<short>       Pedestal{ fhicl::Name("Pedestal"), 1500}; 
          fhicl::Atom<std::string> DigiDataFile{ fhicl::Name("DigiDataFile") }; 
          fhicl::Atom<size_t>      DigiDataColumn{ fhicl::Name("DigiDataColumn"), 1 }; 
          fhicl::Atom<Int_t>       Samples{ fhicl::Name("Samples"), 699 }; 
          fhicl::Atom<Double_t>    Scale{ fhicl::Name("Scale") }; 
          fhicl::Atom<bool>        ApplyPrefilter{ fhicl::Name("ApplyPrefilter"), false };
          struct Filter {
            fhicl::Atom<std::string> Name{fhicl::Name("Name"), "Gauss"};
            fhicl::Atom<Double_t>    Cutoff{fhicl::Name("Cutoff"), 32};  
          };
          fhicl::Table<Config::Filter> Filter{ fhicl::Name("WfmFilter") };

          struct Prefilter {
            fhicl::Atom<std::string> Name{fhicl::Name("Name"), "Gauss"}; 
            fhicl::Atom<Double_t>    Cutoff{fhicl::Name("Cutoff")}; 
          }; 
          fhicl::Table<Config::Prefilter> Prefilter{ fhicl::Name("WfmPrefilter") };
      };

      using Parameters = art::EDProducer::Table<Config>;

      //! Input signal shape model 
      enum EInputShape {kDelta = 0, kScint = 1};
      //! Waveform filter type
      enum EFilterType {kOther = 0, kWiener = 1, kGauss = 2}; 

      struct WfmPrefilter_t {
        TString fName; 
        Double_t fCutoff; 

        WfmPrefilter_t() {fName = ""; fCutoff = 1.0;}

        WfmPrefilter_t(const char* name, Double_t cutoff) {
          fName = name; fCutoff = cutoff;
        }

        WfmPrefilter_t(const struct Config::Prefilter& config) {
          fName = config.Name();
          fCutoff = config.Cutoff(); 
        }
      }; 

      struct WfmFilter_t {
        TString fName; 
        EFilterType fType; 
        Double_t fCutoff; 

        WfmFilter_t() : fName("Unknown_filter"), fType(kOther), fCutoff(32) {}
        WfmFilter_t(const char* name) : fName(name), fType(kOther), fCutoff(32) {}
        WfmFilter_t(const struct Config::Filter& config) {
          fName = config.Name(); 
          fCutoff = config.Cutoff(); 
          if (fName == "Wiener") fType = kWiener; 
          else if (fName == "Gauss") fType = kGauss; 
          else fType = kOther; 
        }
      };

      /**
       * @brief Helper struct to handle waveforms' FT
       *
       * Helper data struct to store and handle waveforms' Fourier Transform. 
       */
      struct CmplxWaveform_t {
        //! vector of Fourier coefficients
        std::vector<TComplex> fCmplx; 
        //! vector of the real part of Fourier coefficients
        std::vector<Double_t> fRe; 
        //! vector of the imaginary part of Fourier coefficients
        std::vector<Double_t> fIm; 

        //! Empty constructor
        inline CmplxWaveform_t() {}
        //! Constructor
        inline CmplxWaveform_t(int size) {
          fCmplx = std::vector(size, TComplex(0., 0)); 
          fRe = std::vector(size, 0.); 
          fIm = std::vector(size, 0.); 
        }

        //! Copy constructor
        inline CmplxWaveform_t(const CmplxWaveform_t& cwf) {
          fCmplx = std::vector( cwf.fCmplx ); 
          fRe = std::vector( cwf.fRe ); 
          fIm = std::vector( cwf.fIm ); 
        }

        //! Destructor
        inline ~CmplxWaveform_t() {
          fCmplx.clear(); fRe.clear(); fIm.clear(); 
        }

        /**
         * @brief Point-side product of complex waveforms
         *
         * Compute the product of two complex waveforms point 
         * by point
         */
        inline CmplxWaveform_t operator*(const CmplxWaveform_t& cwf) const {
          CmplxWaveform_t result(fCmplx.size()); 
          if (cwf.fCmplx.size() != fCmplx.size()) {
            printf("opdet::Deconvolution_module::CmplxWaveform_t::operator* ERROR"); 
            printf(" waveforms size does not match.\n"); 
            return result;
          }

          for (size_t i=0; i<fCmplx.size(); i++) {
            result.fCmplx.at(i) = fCmplx.at(i) * cwf.fCmplx.at(i); 
            result.fRe.at(i) = result.fCmplx.at(i).Re(); 
            result.fRe.at(i) = result.fCmplx.at(i).Im(); 
          }

          return result;
        }

        //! Set the complex coefficient `i` given its real and imaginary part
        inline void MakeCmplx(size_t i) { 
          fCmplx.at(i) = TComplex(fRe.at(i), fIm.at(i)); 
        }

        //! Set the complex coefficients given real and imaginary parts
        //!
        //! Note that the second half of the waveform is set to zero
        inline void MakeCmplx() {
          for (size_t i=0; i<fCmplx.size()*0.5+1; i++) {
            MakeCmplx(i);
          }
          for (size_t i=fCmplx.size()*0.5+1; i<fCmplx.size(); i++) {
            fCmplx.at(i) = TComplex(0, 0); 
            fRe.at(i) = 0.; 
            fIm.at(i) = 0.;
          }
        }

        //! Set the real and imaginary part of the coefficient `i` 
        inline void MakeReAndIm(size_t i) {
          fRe.at(i) = fCmplx.at(i).Re(); 
          fIm.at(i) = fCmplx.at(i).Im(); 
        }

        //! Set real and imaginary parts from the complex coefficients
        //!
        //! Note that the second half of the waveform is set to zero
        inline void MakeReAndIm() {
          for (size_t i=0; i<fCmplx.size()*0.5+1; i++) {
            MakeReAndIm(i); 
          } 
          for (size_t i=fCmplx.size()*0.5+1; i<fCmplx.size(); i++) {
            fCmplx.at(i) = TComplex(0., 0.); 
            fRe.at(i) = 0.;
            fIm.at(i) = 0.; 
          }
        }
      };


      explicit Deconvolution(Parameters const&);
      virtual ~Deconvolution();
      void produce(art::Event& evt);

      // Parameters we'll read from the fcl-file
      std::string fInputModule;                 //!< Module used to create OpDetWaveforms
      std::string fInstanceName;                //!< Input tag for OpDetWaveforms collection
      double fSampleFreq;                       //!< Sampling frequency in MHz 
      double  fTimeBegin;                       //!< Beginning of waveform in us
      double  fTimeEnd;                         //!< End of waveform in us
      short  fPedestal;                         //!< In ADC counts
      double  fLineNoiseRMS;                    //!< Pedestal RMS in ADC counts
      size_t fPreTrigger;                       //!< In ticks
      std::vector<double> fSinglePEWaveform;    //!< Template for a single PE in ADC
      double fSinglePEAmplitude;                //!< single PE amplitude
      unsigned int WfDeco;                      //!< nr of waveform processed
      std::string fDigiDataFile;                //!< single p.e. template source file
      size_t fDigiDataColumn;                      //!< single p.e. template source file column    
      double fScale;                            //!< ???
      size_t fReadoutWindow;                    //!< In ticks
      int fSamples;                             //!< (Same as ReadoutWindow?)
      bool fApplyPrefilter;
      WfmPrefilter_t fPrefilterConfig;
      WfmFilter_t fFilterConfig; 
      EInputShape fInputShape = kDelta; 

      //----------------------------------------------------
      // Declare member functions
      //std::vector<raw::OpDetWaveform> RunDeconvolution(std::vector<raw::OpDetWaveform> const& wfHandle);

      //------------------------------------------------------
      //Load TFileService service
      art::ServiceHandle<art::TFileService> tfs;

    private:
      int  CountFileColumns(const char* file_path);
      void SourceSPEDigiDataFile(); 
      void BuildPrefilter(CmplxWaveform_t& xF0); 
      void ComputeExpectedInput(std::vector<double>& s, double nmax);
      Double_t ComputeNormalization(CmplxWaveform_t& xGH); 
      TVirtualFFT* fft_r2c; 
      TVirtualFFT* fft_c2r; 
  };
}
#endif

namespace opdet{
  DEFINE_ART_MODULE(Deconvolution)
}

namespace opdet {
  //---------------------------------------------------------------------------
  // Constructor
  Deconvolution::Deconvolution(const Parameters& pars)
    : EDProducer{pars}, 
    fInputModule{ pars().InputModule()}, 
    fInstanceName{ pars().InstanceName()},
    
    fPedestal{ pars().Pedestal()},
    fLineNoiseRMS{ pars().LineNoiseRMS() },
    fPreTrigger{ pars().PreTrigger()},

    fDigiDataFile{ pars().DigiDataFile()},
    fDigiDataColumn{ pars().DigiDataColumn()},
    fScale{ pars().Scale()},
    fReadoutWindow{ pars().ReadoutWindow()},
    fSamples{ pars().Samples()},
    fPrefilterConfig{ WfmPrefilter_t( pars().Prefilter()) }, 
    fFilterConfig{ WfmFilter_t( pars().Filter() ) }
  {  

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

    fft_r2c = TVirtualFFT::FFT(1, &fSamples, "M R2C K");
    fft_c2r = TVirtualFFT::FFT(1, &fSamples, "M C2R K");

    SourceSPEDigiDataFile(); 
  }

  //---------------------------------------------------------------------------
  // Destructor
  Deconvolution::~Deconvolution(){
    delete fft_r2c; 
    delete fft_c2r;
  }    


  //-------------------------------------------------------------------------
  void Deconvolution::produce(art::Event& evt){

    art::Handle< std::vector< raw::OpDetWaveform > > wfHandle;
    evt.getByLabel(fInputModule, fInstanceName, wfHandle);

    // Access ART's TFileService, which will handle creating and writing
    art::ServiceHandle< art::TFileService > tfs;

    //******************************
    //-- Read Waveform----
    //****************************** 

    std::vector<raw::OpDetWaveform> digi_wave = *wfHandle;
    int NOpDetWaveform = digi_wave.size();
    std::cout << NOpDetWaveform << std::endl;

    //pointer that will store produced Waveform
    auto out_wave = std::make_unique< std::vector< raw::OpDetWaveform > >();
    auto out_decowave = std::make_unique< std::vector< recob::OpWaveform > >();
    auto out_prefilter = std::make_unique< std::vector< recob::OpWaveform >>(); 

    std::vector<short unsigned int > out_digiwave(NOpDetWaveform); //vector in which the waveform will be saved
    std::vector<float> out_recowave(NOpDetWaveform);               //vector in which the decowaveform will be saved, using float
    std::vector<float> out_recopref(NOpDetWaveform); 
    std::vector<double> xv(fSamples, 0.); 

    //******************************
    //--- Setup filter's components 
    //******************************
    CmplxWaveform_t xG0(fSamples); 
    if (fApplyPrefilter) BuildPrefilter(xG0);

    // xH: FFT of spe response
    std::vector<double> xh(fSinglePEWaveform);
    CmplxWaveform_t xH(fSamples);  
    fft_r2c->SetPoints(&xh[0]);  
    fft_r2c->Transform();
    fft_r2c->GetPointsComplex(&xH.fRe[0], &xH.fIm[0]);
    xH.MakeCmplx(); 
    if (fApplyPrefilter) { xH = xG0 * xH; }

    //******************************
    //--- Process waveforms
    //******************************
    for (auto const& wf: digi_wave) {
      CmplxWaveform_t xV(fSamples); 
      CmplxWaveform_t xS(fSamples); 
      CmplxWaveform_t xG(fSamples); 
      CmplxWaveform_t xY(fSamples); 
      CmplxWaveform_t xGH(fSamples); 
      //----------------------------------------------------------- Resize wvfs
      if (static_cast<int>(wf.Waveform().size()) <= fSamples) { 
        out_recowave.resize(fSamples,0); 
        out_recopref.resize(fSamples,0); 
      }

      else {
        printf("\nWARNING: waveform size is %lu, which is larger than fSamples (%i)\n", 
            wf.Waveform().size(), fSamples); 
        out_recowave.resize(fSamples);
        out_recopref.resize(fSamples); 
      }

      for (Int_t i= 0; i < fSamples; i++){
        // Remove baseline 
        if (i < static_cast<int>(wf.Waveform().size())) xv[i] = (wf[i]-fPedestal);
        // if waveform is shorter than fSamples fill the rest with noise
        else xv[i] = CLHEP::RandGauss::shoot(0, fLineNoiseRMS); 
      }


      //---------------------------------------------------- Guess input signal
      // Found maximum peak in the Waveform
      Double_t SPE_Max = 0;
      double maxADC=*max_element(wf.begin(),wf.end());
      double maxAmplit= maxADC-fPedestal;
      SPE_Max = maxAmplit/fSinglePEAmplitude;
      //printf("SPE_Max = (%g - %i) / %g = %g\n", 
          //maxADC, fPedestal, fSinglePEAmplitude, SPE_Max); 

      std::vector<double>xs(fSamples,0.);
      // Compute expected input (using a delta or the scint tile profile
      // as a model)
      ComputeExpectedInput(xs, SPE_Max); 

      //-------------------------------------------------- Compute waveform FFT
      fft_r2c->SetPoints(&xv[0]);
      fft_r2c->Transform();
      fft_r2c->GetPointsComplex(&xV.fRe[0], &xV.fIm[0]);
      xV.MakeCmplx(); 

      //----------------------------------------------------- Compute input FFT
      fft_r2c->SetPoints(&xs[0]);
      fft_r2c->Transform();
      fft_r2c->GetPointsComplex(&xS.fRe[0], &xS.fIm[0]);
      xS.MakeCmplx(); 

      //******************************
      // Compute filters.
      //******************************                             
      for (int i=0; i<fSamples*0.5+1; i++) {
        // Compute spectral density
        double H2 = xH.fCmplx.at(i).Rho2();
        double S2 = xS.fCmplx.at(i).Rho2();
        double N2 = fLineNoiseRMS * fLineNoiseRMS * fSamples ;

        if (fApplyPrefilter) {
          Double_t prefilter_PSD = xG0.fCmplx.at(i).Rho2(); 
          N2 *= prefilter_PSD;
        }

        if (fFilterConfig.fType == Deconvolution::kWiener){
          // Compute Wiener filter
          xG.fCmplx.at(i) = TComplex::Conjugate(xH.fCmplx.at(i))*S2 / (H2*S2 + N2);
        }
        else if (fFilterConfig.fType == Deconvolution::kGauss){
          // Compute gauss filter
          Double_t gauss_cutoff = fFilterConfig.fCutoff; 
          xG.fCmplx[0] = TComplex(0,0);
          xG.fCmplx.at(i) = TComplex::Exp(
              -0.5*TMath::Power(i*1e-6*fSampleFreq/(fSamples*gauss_cutoff),2))
            /xH.fCmplx.at(i);
        }
        else{
          // Compute dec signal
          xG.fCmplx.at(i) = TComplex::Power(xH.fCmplx.at(i),-1);// Standard dec is just the division of signal and SPE template in Fourier space
        }

        // Correct template pretrigger (= phase shift in the freq domain)
        TComplex phase = TComplex(0., -TMath::TwoPi()*i*fPreTrigger/(fSamples)); 
        xG.fCmplx.at(i) = xG.fCmplx.at(i)*TComplex::Exp(phase);
        xG.MakeReAndIm(i); 
        //printf("[%i] H2 = %g, N2 = %g, S2 = %g -> W2 = %g \n", 
            //i, H2, N2, S2, xG.fCmplx.at(i).Rho2()); 
      }

      if (fApplyPrefilter) {xV = (xG0 * xV); xV.MakeReAndIm();}

      // Apply filter to the waveform
      xY  = xG * xV; 
      // Apply filter to the detector response (for normalization)
      xGH = xG * xH; 
      xGH.MakeReAndIm(); 
      xY.MakeReAndIm();

      fft_c2r->SetPointsComplex(&xV.fRe.at(0), &xV.fIm.at(0)); 
      fft_c2r->Transform(); 
      Double_t* xPrefilter = fft_c2r->GetPointsReal(); 
      Double_t scale = 1.0 / fSamples; 
      for (int i=0; i<fSamples; i++) {
        out_recopref.at(i) = xPrefilter[i]*scale; 
      }

      Double_t filter_norm = ComputeNormalization(xGH); 
      scale = filter_norm / (Double_t)fSamples;  
      //printf("scale = %g/%i =  %g\n", filter_norm, fSamples, scale);  
      //getchar(); 

      //Transform of the filtered signal
      fft_c2r->SetPointsComplex(&xY.fRe[0], &xY.fIm[0]);
      fft_c2r->Transform();
      double *xy = fft_c2r->GetPointsReal();
      for (int i=0; i<fSamples; i++){ 
        out_recowave[i] = xy[i]*scale;
        // out_recowave[i] = xy[i];
        // std::cout << xy[i] << std::endl;
      }

      raw::OpDetWaveform dgwave( wf.TimeStamp(), wf.ChannelNumber(), out_digiwave );
      out_wave->emplace_back(std::move(dgwave));

      recob::OpWaveform decpref(wf.TimeStamp(), wf.ChannelNumber(), out_recopref );
      out_prefilter->emplace_back(std::move(decpref));   

      recob::OpWaveform decwav(wf.TimeStamp(), wf.ChannelNumber(), out_recowave );
      out_decowave->emplace_back(std::move(decwav));   
    }//waveforms loop

    // Push the OpDetWaveforms and OpWaveform into the event
    evt.put(std::move(out_wave));
    //evt.put(std::move(out_prefilter)); 
    evt.put(std::move(out_decowave));
    WfDeco++;
  }


  /**
   * @brief Build a filter to be applied prior the deconvolution
   *
   * Construct a filter to be applied before the deconvolution process.
   * Different filters can be implemented by switching the flag `fPrefilterConfig.fName`
   * via the Config::Prefilter::name parameter. 
   *
   * @param xF
   */
  void Deconvolution::BuildPrefilter(CmplxWaveform_t& xF) {
    if (fPrefilterConfig.fName != "Gauss") {
      printf("Deconvolution::BuildPrefilter WARNING: Unknown filter model %s. Skip.\n", 
          fPrefilterConfig.fName.Data()); 
      return;
    } 

    // Compute sigma corresponding to the given cutoff frequency
    const Double_t df       = fSampleFreq / (Double_t)fSamples; 
    const Double_t cutoff   = fPrefilterConfig.fCutoff / df; 
    const Double_t k_cutoff = sqrt(log(2)); 
    const double sigma = fSamples * k_cutoff / (TMath::TwoPi() * cutoff);
    const int    mu    = 4*sigma; 

    printf("Deconvolution::BuildPrefilter sigma is %g\n", sigma); 

    std::vector<Double_t> xf(fSamples, 0.);
    for (int i=0; i<fSamples; i++) xf.at(i) = TMath::Gaus(i, mu, sigma, kTRUE); 

    std::vector<Double_t> re_(fSamples, 0.); 
    std::vector<Double_t> im_(fSamples, 0.);
  
    fft_r2c->SetPoints(&xf[0]); 
    fft_r2c->Transform(); 
    fft_r2c->GetPointsComplex(&re_[0], &im_[0]);

    for (int i=0; i<0.5*fSamples+1; i++) {
      TComplex F(re_.at(i), im_.at(i));  
      TComplex phase = TComplex(0., -TMath::TwoPi()*i*mu/(fSamples));  
      xF.fCmplx.at(i) = F*TComplex::Exp(phase); 
      xF.MakeReAndIm(i); 
    }



    //std::ofstream dump_filter; 
    //dump_filter.open("prefilter_dump.txt"); 
    //for (int i=0; i<0.5*fSamples + 1; i++) dump_filter << xF.fCmplx.at(i).Rho2() <<"\n";
    //dump_filter.close(); 

    return; 
  }

  
  /**
   * @brief Compute expected input signal
   *
   * Produce a waveform representing a guess of the input signal based on the
   * estimated nr of p.e. in the waveform peak. Based on the `fInputShape`
   * variable, the input shape can be assumed as a δ-function scaled for 
   * the nr of p.e. or as the LAr scintillation time profile
   *
   * @param s input signal
   * @param nmax estimated nr of p.e. at the waveform max
   */
  void Deconvolution::ComputeExpectedInput(std::vector<double>& s, double nmax) {
    if (fInputShape == kScint) {
      //ST profile 
      auto const *Larprop = lar::providerFrom<detinfo::LArPropertiesService>();
      std::vector<double>SignalTime={(Larprop->ScintFastTimeConst()* 0.001),(Larprop->ScintSlowTimeConst()* 0.001)};
      std::vector<double>SignalScint={Larprop->ScintYieldRatio(),1.-Larprop->ScintYieldRatio()};

      double dt = 1/fSampleFreq; 
      double t = fTimeBegin;

      for (size_t i=0; i<s.size(); i++) {
        double lightsignal=0;
        for (size_t j=0; j<SignalTime.size();j++){
          lightsignal+=SignalScint[j]*exp(-t/SignalTime[j]);
        }
        s.at(i) = nmax*lightsignal; 
        t+=dt; 
      }
    }
    else if (fInputShape == kDelta) {
      s.at(1) = nmax; 
    }
    else {
      printf("Deconvolution::ComputeExpectedInput WARNING\n"); 
      printf("Unknown input shape: assuming Deconvolution::kDelta\n"); 
      s.at(1) = nmax;
    }
    return; 
  }

  /**
   * @brief Source the single p.e. response from file
   *
   * Source the single p.e. response template from the txt file set by 
   * `fDigiDataFile` and set the variable `fSinglePEAmplitude` with the 
   * amplitude of the single p.e. response. In case of a multi-column
   * template file, the relevant column can be selected by setting the
   * varible `fDigiDataColumn`.
   */
  void Deconvolution::SourceSPEDigiDataFile() {
    size_t n_columns = CountFileColumns(fDigiDataFile.c_str()); 
    if (fDigiDataColumn >= n_columns) {
      printf("Deconvolution::SourceSPETemplate ERROR: "); 
      printf("The module is supposed to select column %lu, but only %lu columns are present.\n", 
          fDigiDataColumn, n_columns); 
      throw art::Exception(art::errors::InvalidNumber); 
    }

    std::ifstream SPEData; 
    SPEData.open(fDigiDataFile); 
    Double_t buff[100] = {0};  

    std::string temp_str; 
    if (SPEData.is_open()) {
      while (std::getline(SPEData, temp_str)) {
        std::stringstream ss; ss << temp_str; 
        int  icol = 0; 
        while (ss) {ss >> buff[icol]; ++icol;}

        fSinglePEWaveform.push_back(buff[fDigiDataColumn]);
      } 
    } else {
      printf("Deconvolution::produce ERROR "); 
      printf("Cannot open SPE template file.\n"); 

      throw art::Exception(art::errors::FileOpenError); 
    }

    fSinglePEWaveform.resize(fSamples, 0);

    SPEData.close(); 

    // Set single p.e. maximum value
    fSinglePEAmplitude = TMath::Max(1.0, 
      *(std::max_element(fSinglePEWaveform.begin(), fSinglePEWaveform.end()))); 

    return;
  }

  /**
   * @brief Count the nr of column in a txt file
   *
   * @param file_path
   *
   * @return nr of columns
   */
  int Deconvolution::CountFileColumns(const char* file_path) {
    std::ifstream file_; 
    file_.open(file_path); 
    
    if (!file_.is_open()) {
      printf("Deconvolution::CountFileColumns(%s) ERROR:\n", 
          file_path); 
      printf("Unable to open file."); 
      throw art::Exception(art::errors::FileOpenError); 
    }

    int N_COLUMNS = 0; 
    std::string line; 
    int iline = 0; 
    while ( std::getline(file_, line) ) {
      std::stringstream sstream; 
      sstream << line;
      std::string sub;
      int n_columns = 0; 

      while (sstream) {
        sstream >> sub; 
        if (sub.length()) ++n_columns; 
      }

      if (iline == 1) {N_COLUMNS = n_columns;}
      else if (iline > 1) {
        if (n_columns != N_COLUMNS) {
          printf("Deconvolution::CountFileColumns(%s): WARNING ", 
              file_path); 
          printf("Nr of columns change along the file!\n"); 
          N_COLUMNS = n_columns; 
        }
      }
      iline++; 
    }
    file_.close();  
    return N_COLUMNS; 
  }

  /**
   * @brief Compute normalization factor for a given filter
   *
   * The filter normalization factor is obtained by applying the 
   * filter to the single p.e. response template (i.e., to a 1 p.e. noiseless signal).
   * the product of the filter should be a the best approximation of the 
   * input function achievable with the signal to noise ratio in such waveform, 
   * thus, the normalization consists in the integral of the product 
   * around the peak in a region defined as the one where the signal is positive. 
   * This factor is supposed to tend to 1 for high SNR (>= 10).
   *
   * @param xGH
   *
   * @return filter normalization
   */
  Double_t Deconvolution::ComputeNormalization(CmplxWaveform_t& xGH) {
    double norm = 0; 

    fft_c2r->SetPointsComplex(&xGH.fRe[0], &xGH.fIm[0]); 
    fft_c2r->Transform(); 
    Double_t* x_ =  fft_c2r->GetPointsReal(); 
    std::vector<Double_t> x(x_, x_ + fSamples);

    Int_t imax = std::distance(x.begin(), std::max_element(x.begin(), x.end())); 
    Int_t ileft = imax; 
    Int_t iright = imax; 

    while (x[ileft] >= 0.0 && ileft > 0) ileft--; 
    while (x[iright] >=0.0 && iright < fSamples) iright++; 

    for (Int_t k=ileft; k<=iright; k++) norm += x[k]; 

    norm /= (Double_t)fSamples; 

    if (norm > 1.0) norm = 1.0; 
    else if (norm <= 0.){
      printf("Deconvolution::ComputeNormalization() WARNING: "); 
      printf(" bad normalization (%g), force to 1.0\n", norm);  
      norm = 1.0; 
    }
    //printf("imax = %i, ileft = %i, iright = %i -> integral = %g\n", 
        //imax, ileft, iright, norm*(fSamples)); 


    return 1.0 / norm; 
  }

} // namespace opdet
