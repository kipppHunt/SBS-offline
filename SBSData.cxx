#include "SBSData.h"
#include "TMath.h"
#include <iostream>
#define SBS_ADC_MODE_SINGLE 0 //< Simple ADC with only integral
#define SBS_ADC_MODE_MULTI  1 //< FADC 250 mode 7

namespace SBSData {

  /////////////////////////////////////////////////////////////////////////////
  // ADC data functions
  ADC::ADC(Double_t ped, Double_t gain, Double_t tcal) :
    fHasData(false), fMode(SBS_ADC_MODE_SINGLE)
  {
    SetPed(ped);
    SetGain(gain);
    SetTimeCal(tcal);
  }

  void ADC::Process(Double_t val)
  {
    SingleData zero = { 0.0, 0.0 };
    SingleData integral = { val, (val-fADC.ped)*fADC.cal };
    fADC.hits.push_back({integral,zero,zero});
    fHasData = true;
    fMode = SBS_ADC_MODE_SINGLE; //< Mode gets set to simple if this function is called
  }

  void ADC::Process(Double_t integral, Double_t time, Double_t amp, Double_t ped) {
    //fADC.push_back({ped,fGlobalCal,val,val-ped,(val-ped)*fGlobalCal});
    // convert to pC, assume tcal is in ns, and 50ohm resistance
    Double_t pC_Conv = fADC.tcal/50.;
    Double_t PedVal = ped*GetChanTomV();
    Double_t TimeVal= time*fADC.tcal/64. + fADC.timeoffset;
    Double_t IntRaw=  integral*GetChanTomV()*pC_Conv;
    Double_t IntVal=  (IntRaw-PedVal*(GetNSA()+GetNSB()+1)*pC_Conv)*GetGain();
    Double_t AmpRaw=  amp*GetChanTomV();
    Double_t AmpVal=  (AmpRaw-PedVal)*GetAmpCal();
    SingleData t_integral = { IntRaw, IntVal   };
    SingleData t_time     = { time, TimeVal };
    SingleData t_amp     = { AmpRaw, AmpVal };
    fADC.hits.push_back({t_integral,t_time,t_amp } );
    SetPed(PedVal);
    fHasData = true;
    fMode = SBS_ADC_MODE_MULTI; //< Mode gets set to multi if this function is called
  }

  void ADC::Clear()
  {
    fADC.good_hit = 0;
    fADC.hits.clear();
    fHasData = false;
  }

  /////////////////////////////////////////////////////////////////////////////
  // TDC data functions
  TDC::TDC(Double_t offset, Double_t cal, Double_t GoodTimeCut) : fHasData(false)
  {
    fEdgeIdx[0] = fEdgeIdx[1]=0;
    SetOffset(offset);
    SetCal(cal);
    SetGoodTimeCut(GoodTimeCut);
    fTrigPhaseCorr = 0.0; 
  }

  void TDC::ProcessSimple(Int_t elemID, Double_t val, Int_t nhit,UInt_t TrigTime)
  {
    fTDC.hits.push_back(TDCHit());
    Int_t size = fTDC.hits.size();
    TDCHit *hit = &fTDC.hits[size-1];
    hit->elemID = elemID;
    hit->TrigTime = TrigTime;
    hit->le.raw = val;
    hit->le.val = (val-fTDC.offset)*fTDC.cal;
    hit->te.raw = 0;
    hit->te.val = 0;
    hit->ToT.raw=0;
    hit->ToT.val=0;
    fHasData = true;
  }

  void TDC::Process(Int_t elemID, Double_t val, Double_t fedge)
  {
    Int_t edge = int(fedge);
    // std::cout << " tdc process " << val << " " << edge  << " ftdc hits size = " <<fTDC.hits.size() << " hits in edge "  << fEdgeIdx[edge]<< std::endl;
    if(edge < 0 || edge>1) {
      std::cerr << "Edge specified is not valid!" << std::endl;
      edge = 0;
    }
    size_t idx = fEdgeIdx[edge]++;
    if(idx >= fTDC.hits.size()) {
      // Must grow the hits array to accomodate the new hit
      // if ( edge ==1)   std::cout << " First edge is TE , this is not right: " << "  idx = " << idx  << " ftdc hits size = " <<fTDC.hits.size() << " hits with LE ="  << fEdgeIdx[0] << " hits with TE ="  << fEdgeIdx[1] << std::endl;
      fTDC.hits.push_back(TDCHit());
    }
    TDCHit *hit = &fTDC.hits[idx];
    hit->elemID = elemID;
    hit->TrigTime = 0;
    if( edge == 0 ) { // Leading edge
      hit->le.raw = val;
      hit->le.val = (val-fTDC.offset)*fTDC.cal;
    } else {
      hit->te.raw = val;
      hit->te.val = (val-fTDC.offset)*fTDC.cal;
    }
    if(fEdgeIdx[0] == fEdgeIdx[1]) { // Both leading and trailing edges now found
      hit->ToT.raw = hit->te.raw - hit->le.raw;
      hit->ToT.val = hit->te.val - hit->le.val;
    }
    if(fEdgeIdx[1] > fEdgeIdx[0]) fEdgeIdx[0] = fEdgeIdx[1]; // if TE found first force LE count to increase
    fHasData = true;
  }

  void TDC::Clear()
  {
    fEdgeIdx[0] = fEdgeIdx[1] = 0;
    fTDC.hits.clear();
    fHasData = false;
    fTDC.good_hit = 0;
  }

  /////////////////////////////////////////////////////////////////////////////
  // Waveform data functions
  Waveform::Waveform(Double_t ped, Double_t gain, Double_t ChanTomV,Double_t GoodTimeCut, Double_t tcal) :
    fHasData(false)
  {
    SetPed(ped);
    SetGain(gain);
    SetChanTomV(ChanTomV);
    SetTimeCal(tcal);
    SetGoodTimeCut(GoodTimeCut);
  }

  void Waveform::Process(std::vector<Double_t> &vals)
  {
    //printf("vals size %d, samples raw size %d\n", vals.size(), fSamples.samples_raw.size());
    if( vals.size() != fSamples.samples_raw.size()) {
      // Resize our data vector
      fSamples.samples_raw.resize(vals.size());
      fSamples.samples.resize(vals.size());
      Clear();
    }
    //
    for(size_t i = 0; i < vals.size(); i++ ) {
      fSamples.samples_raw[i] = vals[i]*fSamples.ChanTomV;
    }
    // Determining pedestal in first and last four samples. Will choose the
    // minimum of the two.
    Double_t pedsum_fsamps=0, pedsum_lsamps=0;
    Int_t NPedsum=GetNPedBin(), totTimeSamps=vals.size();
    
    NPedsum= TMath::Min(NPedsum,int(vals.size()));
    // std::cout << " Npedsum = " << NPedsum << " " << GetNPedBin() << std::endl;
    
    for(Int_t i = 0; i < NPedsum ; i++ ) {
      // Computing pedestal using first N samples
      pedsum_fsamps += fSamples.samples_raw[i];
      // Computing pedestal using last N samples
      pedsum_lsamps += fSamples.samples_raw[(totTimeSamps-1)-i];
    }
   
    SetPed( TMath::Min(pedsum_fsamps,pedsum_lsamps)/NPedsum );
    for(size_t i = 0; i < vals.size(); i++ ) {
      fSamples.samples[i] = (fSamples.samples_raw[i]-fSamples.ped)*fSamples.cal;
    }
    
    // Try and fixd sample threshold crossing above the pedestal
    Double_t ThresVal=GetThres(); // mV
    UInt_t ThresCrossBin=TMath::Max(NPedsum-1,0);
    //    std::cout << " ped = " << fSamples.ped << " thres = " << ThresVal << std::endl ;
    while ( ThresCrossBin < vals.size() && fSamples.samples_raw[ThresCrossBin] < fSamples.ped+ThresVal ) {
      ThresCrossBin++;
    }
    //
    // if threshold crossing found
    UInt_t NSB = GetNSB();
    UInt_t NSA = GetNSA();
    UInt_t FixedThresCrossBin=GetFixThresBin();
    Double_t FineTime = 0;
    Double_t max  = 0;
    Double_t sum = 0;
    Double_t sped = 0;
    
    Bool_t PeakFound= kFALSE;
    UInt_t PeakBin= 0;
    UInt_t IntMinBin= 0;
    UInt_t IntMaxBin= vals.size();
    if (ThresCrossBin < vals.size()) {
      IntMinBin= TMath::Max(ThresCrossBin-NSB,IntMinBin);
      IntMaxBin= TMath::Min(ThresCrossBin+NSA-1,IntMaxBin);
    } else {
      IntMinBin = TMath::Max(FixedThresCrossBin-NSB,IntMinBin);
      IntMaxBin = TMath::Min(FixedThresCrossBin+NSA-1,IntMaxBin);
    }
    // convert to pC, assume tcal is in ns, and 50ohm resistance
    Double_t pC_Conv = fSamples.tcal/50.;
    
    for(size_t i =IntMinBin ; i <IntMaxBin ; i++ ) {
      sped+=fSamples.ped*pC_Conv;
      sum+=fSamples.samples_raw[i]*pC_Conv;
      if ( i >= ThresCrossBin && !PeakFound) {
	if (fSamples.samples_raw[i] > max) {
	  max = fSamples.samples_raw[i];
	} else {
	  PeakFound= kTRUE;
	  PeakBin = i-1;
	}
      }
    }

    
    //
    //    std::cout << " Int = " << IntMinBin << " " << IntMaxBin<< " ThresCrossBin =   " << ThresCrossBin << " peak-found " << PeakFound << std::endl ;
    //
    Double_t VMid = (max+fSamples.ped)/2.;
    if (PeakFound) {
      for(size_t i =IntMinBin ; i <PeakBin+1 ; i++ ) {
	if (VMid >= fSamples.samples_raw[i]  && VMid < fSamples.samples_raw[i+1]) {
	  FineTime = i+(VMid-fSamples.samples_raw[i])/(fSamples.samples_raw[i+1]-fSamples.samples_raw[i]);
	}
      }
    }
    
    /*
      if (ThresCrossBin==vals.size()) {
      std::cout << " Threshold = " << fThresVal << " ped = " << fSamples.ped << " element = " << ThresCrossBin << " " << " Vmid = " << VMid << "  sum = " << sum << " max = " << max << " time = " << FineTime << " tcal = " << fSamples.tcal <<  std::endl;
      }
    */
    /*
    if(FineTime == 0.0){
      std::cout << "time 0 at finetime and peakfound is: " << PeakFound << std::endl;
    }
    */
    fSamples.pulse.integral.raw = sum;
    fSamples.pulse.integral.val = (sum-sped)*fSamples.cal;
    fSamples.pulse.time.raw = FineTime;
    fSamples.pulse.time.val = (FineTime)*fSamples.tcal + fSamples.timeoffset;
    fSamples.pulse.amplitude.raw = max;
    fSamples.pulse.amplitude.val = (max-fSamples.ped)*fSamples.acal;
    if (max==0) fSamples.pulse.amplitude.val=max;
    //
    fHasData = (vals.size() > 0);
  }

  void Waveform::ProcessMulti(std::vector<Double_t> &vals)
  {
    //printf("vals size %d, samples raw size %d\n", vals.size(), fSamples.samples_raw.size());
    if( vals.size() != fSamples.samples_raw.size()) {
      // Resize our data vector
      fSamples.samples_raw.resize(vals.size());
      fSamples.samples.resize(vals.size());
      Clear();
    }
    //
    for(size_t i = 0; i < vals.size(); i++ ) {
      fSamples.samples_raw[i] = vals[i]*fSamples.ChanTomV;
    }
    // Determining pedestal in first and last four samples. Will choose the
    // minimum of the two.
    Double_t pedsum_fsamps=0, pedsum_lsamps=0;
    Int_t NPedsum=GetNPedBin(), totTimeSamps=vals.size();
    
    NPedsum= TMath::Min(NPedsum,int(vals.size()));
    // std::cout << " Npedsum = " << NPedsum << " " << GetNPedBin() << std::endl;
    
    for(Int_t i = 0; i < NPedsum ; i++ ) {
      // Computing pedestal using first N samples
      pedsum_fsamps += fSamples.samples_raw[i];
      // Computing pedestal using last N samples
      pedsum_lsamps += fSamples.samples_raw[(totTimeSamps-1)-i];
    }
   
    SetPed( TMath::Min(pedsum_fsamps,pedsum_lsamps)/NPedsum );
    for(size_t i = 0; i < vals.size(); i++ ) {
      fSamples.samples[i] = (fSamples.samples_raw[i]-fSamples.ped)*fSamples.cal;
    }
    
    // Look for local maxima instead of just threshold crossings
    Double_t ThresVal=GetThres(); // mV
    //max betwen NPedsum-1 and 1 since we want to compare to a previous sample as well
    //should it be just NPedsum instead of subtracting 1?
    UInt_t LocalMaxBin=TMath::Max(NPedsum-1,1);
    std::vector<UInt_t> LocalMaxBinList;
    std::vector<UInt_t> ThresCrossBinList;
    while ( LocalMaxBin < vals.size() ) {
      Double_t currSamp = fSamples.samples_raw[LocalMaxBin], prevSamp = fSamples.samples_raw[LocalMaxBin-1], nextSamp = fSamples.samples_raw[LocalMaxBin+1];
      //find threshold crossings
      if( currSamp >= fSamples.ped+ThresVal && prevSamp < fSamples.ped+ThresVal ){
	ThresCrossBinList.push_back(LocalMaxBin);
      }
      //look for local max above thres
      if( currSamp >= fSamples.ped+ThresVal && currSamp > prevSamp && currSamp > nextSamp){
	LocalMaxBinList.push_back(LocalMaxBin);
      }
      LocalMaxBin++;
    }
    //

    //loop over maxima and do peak erasing, taken from SBSGEMModule::find_clusters_1D_experimental

    vector<int> peaks_to_erase;

    int n_peaks_erased = 0;
    
    for( auto maxima_bin_iter : LocalMaxBinList ){

      double ADCmax = fSamples.samples_raw[maxima_bin_iter];
      double prominence = ADCmax;
      int samplo = maxima_bin_iter, samphi = maxima_bin_iter;
      double ADCminright = ADCmax, ADCminleft = ADCmax;

      bool higherpeakleft = false, higherpeakright = false;
      int peakright = -1, peakleft = -1;

      //sigma_max should be based on the individual ch threshold which is defined as something like 3x the pedrms for that ch (atleast 3x for ecal, need to check hcal)
      double sigma_max = ThresVal/3.0;

      //search for potentially overlapping pulses to the right of current local max
      while( fSamples.samples_raw.find( samphi + 1 ) != fSamples.samples_raw.end() ){
	
	samphi++;
	double ADCtest = fSamples.samples_raw[samphi];

	//want to find lowest sample before next peak
	if( ADCtest < ADCminright && !higherpeakright ){
	  ADCminright = ADCtest;
	}

	//if this eval true then we have found a higher peak to the right of this local maxima so we will want to erase it if it doesnt pass prominence threshold
	if( LocalMaxBinList.find(samphi) != LocalMaxBinList.end() && ADCtest > ADCmax ){
	  higherpeakright = true;
	  peakright = samphi;
	  break;
	}

      }

      //search for potentially overlapping pulses to the left of current local max
      while( fSamples.samples_raw.find( samplo - 1 ) != fSamples.samples_raw.end() ){
	
	samplo++;
	double ADCtest = fSamples.samples_raw[samplo];

	//want to find lowest sample before next peak
	if( ADCtest < ADCminleft && !higherpeakleft ){
	  ADCminleft = ADCtest;
	}

	//if this eval true then we have found a higher peak to the left of this local maxima so we will want to erase it if it doesnt pass prominence threshold
	if( LocalMaxBinList.find(samplo) != LocalMaxBinList.end() && ADCtest > ADCmax ){
	  higherpeakleft = true;
	  peakleft = samplo;
	  break;
	}

      }

      bool peak_close = false;
      if( !higherpeakleft ) ADCminleft = 0.0;
      if( !higherpeakright ) ADCminright = 0.0;

      double sigma = ThresVal/3.0;
      double prom_thres_sigma_factor = 1.0;

      int max_peak_time_sep = 6;
 
      //prominence is vertical diff btwn peaks, adc val diff
      //max separation here will be time difference btwn peaks, horizontal
      if( higherpeakright || higherpeakleft ){
	
	prominence = ADCmax - std::max ( ADCminleft, ADCminright );// define prom as diff btwn current max and highest valley to left or right

	if( higherpeakleft && std::abs( peakleft - maxima_bin_iter ) <= 2*max_peak_time_sep ) peak_close = true;
	if( higherpeakright && std:: abs( peakright - maxima_bin_iter ) <= 2*max_peak_time_sep ) peak_close = true;

	//for now just defining prom threshold, include fraction of current peak later
	if( peak_close && prominence < prom_thres_sigma_factor*sigma ) n_peaks_erased++;

	peaks_to_erase.push_back( maxima_bin_iter )
	
      }

    }//end of loop over local max

    for( auto peak_ind : peaks_to_erase ){
      LocalMaxBinList.erase(peak_ind);
    }

    fSamples.multipulse.resize(LocalMaxBinList.size());

    //now, similar to the SBSGEMModule method, we want to split remaining significant pulses and split any overlapping pulses contributions
   
    //first idea: scale example elastic pulse per ch to match their first samp above threshold and their peak samp, then use that scaled version of the example elastic pulse to estimate the actual decay of the first pulse to be subtracted from the samples convolved with following pulse
    
    //need some provision in DB or somewhere to put example elastic pulses

    for(int local_max_iter = 0; local_max_iter < LocalMaxBinList.size(); local_max_iter++){

      //starting assumption: if we can find a clean start to the rise of the first pulse in the waveform and its peak sample, then we can scale the "reference" pulse linearly to fit the shape determined by the starting point and peak of the pulse
      //find ratio of number of samps between start and peak of "reference" pulse and number of samps between start and peak of the actual pulse here to find horizontal scale factor
      //find ratio of peak amplitudes between "reference" pulse and actual pulse here to find vertical scale factor
      //then use these factors to scale the "reference" pulse to match the actual pulse here and then find the amplitude at the time-values (x-axis) coming after the peak sample
      //this is currently giving the sample values of the template pulses
      std::vector<Double_t> template_pulse = GetExamplePulseParams();
      TSpline3 *cubic_spline = new TSpline3("cubic_spline",0,template_pulse.size(),template_pulse.data(),template_pulse.size());
      Int_t template_dist_to_peak = std::max_element(template_pulse.begin(), template_pulse.end());
      Int_t SampsFromStartToPeak = LocalMaxBinList[local_max_iter] - ThresCrossBinList[local_max_iter];
      Double_t template_peak = template_pulse[template_dist_to_peak];
      Double_t real_pulse_peak = fSamples.samples_raw[LocalMaxBinList[local_max_iter]];
      Double_t h_scale = SampsFromStartToPeak / template_dist_to_peak;
      Double_t v_scale = real_pulse_peak / template_peak;
      Int_t max_samps_after_peak = 5;
      for(int samp_after_peak = SampsFromStartToPeak; samp_after_peak < max_samps_after_peak; samp_after_peak++){
	Double_t scaled_samp = v_scale*cubic_spline->Eval(samp_after_peak*h_scale);
      }
      
    }
    
    //Loop over all threshold crossings
    for(int ncross = 0; ncross < ThresCrossBinList.size(); ncross++){
      
      // if threshold crossing found
      UInt_t NSB = GetNSB();
      UInt_t NSA = GetNSA();
      UInt_t FixedThresCrossBin=GetFixThresBin();
      Double_t FineTime = 0;
      Double_t max  = 0;
      Double_t sum = 0;
      Double_t sped = 0;
      
      Bool_t PeakFound= kFALSE;
      UInt_t PeakBin= 0;
      UInt_t IntMinBin= 0;
      UInt_t IntMaxBin= vals.size();
      
      IntMinBin= TMath::Max(ThresCrossBinList[ncross]-NSB,IntMinBin);
      IntMaxBin= TMath::Min(ThresCrossBinList[ncross]+NSA-1,IntMaxBin);
      
      //no use for FixedThresCrossBin due to this for loop requiring ThresCrossBinList having atleast one element, so IntMin/MaxBin are always defined by comparison of their initial value to the threshold crossing bin -NSB or +NSA for min or max respectively
      
      // convert to pC, assume tcal is in ns, and 50ohm resistance
      Double_t pC_Conv = fSamples.tcal/50.;
      
      for(size_t i =IntMinBin ; i <IntMaxBin ; i++ ) {
	if(PeakFound && fSamples.samples_raw[i] < fSamples.ped+ThresVal) break; //truncate integration if pulse drops below thres before going through window of integration to avoid contribution from later pulse
	sped+=fSamples.ped*pC_Conv;
	sum+=fSamples.samples_raw[i]*pC_Conv;
	if ( i >= ThresCrossBinList[ncross] && !PeakFound) {
	  if (fSamples.samples_raw[i] > max) {
	    max = fSamples.samples_raw[i];
	  } else { //look for first decreasing sample after threshold crossing
	    PeakFound= kTRUE;
	    PeakBin = i-1;
	  }
	}
      }
      //
      //    std::cout << " Int = " << IntMinBin << " " << IntMaxBin<< " ThresCrossBin =   " << ThresCrossBin << " peak-found " << PeakFound << std::endl ;
      //
      Double_t VMid = (max+fSamples.ped)/2.;
      if (PeakFound) {
	for(size_t i =IntMinBin ; i <PeakBin+1 ; i++ ) {
	  if (VMid >= fSamples.samples_raw[i]  && VMid < fSamples.samples_raw[i+1]) {
	    FineTime = i+(VMid-fSamples.samples_raw[i])/(fSamples.samples_raw[i+1]-fSamples.samples_raw[i]);
	    /*
	    std::cout << "good finetime" << std::endl;
	  }else{
	    std::cout << "Vmid: " << VMid << " ped: " << fSamples.ped << " , fsamp i: " << fSamples.samples_raw[i] << " , fsamp i+1: " << fSamples.samples_raw[i+1] << " ,i is: " << i << std::endl;
	    */
	  }
	}
      }
      
      /*
	if (ThresCrossBin==vals.size()) {
	std::cout << " Threshold = " << fThresVal << " ped = " << fSamples.ped << " element = " << ThresCrossBin << " " << " Vmid = " << VMid << "  sum = " << sum << " max = " << max << " time = " << FineTime << " tcal = " << fSamples.tcal <<  std::endl;
	}
      
      if(FineTime == 0.0){
	std::cout << "time 0 at finetime and peakfound is: " << PeakFound << std::endl;
      }
      */
      fSamples.multipulse[ncross].integral.raw = sum;
      fSamples.multipulse[ncross].integral.val = (sum-sped)*fSamples.cal;
      fSamples.multipulse[ncross].time.raw = FineTime;
      fSamples.multipulse[ncross].time.val = (FineTime)*fSamples.tcal + fSamples.timeoffset;
      fSamples.multipulse[ncross].amplitude.raw = max;
      fSamples.multipulse[ncross].amplitude.val = (max-fSamples.ped)*fSamples.acal;
      if (max==0) fSamples.multipulse[ncross].amplitude.val=max;
      //
      //std::cout << "starting new threshold crossing" << std::endl;
    }
    fHasData = (vals.size() > 0);
  }
  
  void Waveform::Clear()
  {
    for(size_t i = 0; i < fSamples.samples.size(); i++) {
      fSamples.samples_raw[i] = fSamples.samples[i] = 0;
    }
    fHasData = false;
  }

}; // end SBSData
