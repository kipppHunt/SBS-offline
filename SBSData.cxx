#include "SBSData.h"
#include "TMath.h"
#include "TF1.h"
#include <algorithm>
#include "TGraph.h"
#include "TPad.h"
#include <iostream>
#include <vector>
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

  void Waveform::ProcessDecon(std::vector<Double_t> &vals)
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
    
    /*
    for(Int_t i = 0; i < NPedsum ; i++ ) {
      // Computing pedestal using first N samples
      pedsum_fsamps += fSamples.samples_raw[i];
      // Computing pedestal using last N samples
      pedsum_lsamps += fSamples.samples_raw[(totTimeSamps-1)-i];
    }
    */

    Double_t ped_sum_min = 100000;
    
    for(Int_t ii = 0; ii < vals.size() - NPedsum; ii++){
      Double_t ped_sum_curr = 0;
      for(Int_t jj = 0; jj < NPedsum; jj++){
	ped_sum_curr += fSamples.samples_raw[ ii + jj ];
      }
      if(ped_sum_curr < ped_sum_min) ped_sum_min = ped_sum_curr;
    }
    
    //SetPed( TMath::Min(pedsum_fsamps,pedsum_lsamps)/NPedsum );
    SetPed(ped_sum_min/NPedsum);
    for(Int_t i = 0; i < vals.size(); i++ ) {
      fSamples.samples[i] = (fSamples.samples_raw[i]-fSamples.ped)*fSamples.cal;
    }

    Double_t smooth_samp[vals.size()] = {0};
    
    for(Int_t i = 1; i < vals.size()-1; i++){
      Double_t currSamp = fSamples.samples[i];
      Double_t prevSamp = fSamples.samples[i-1];
      Double_t nextSamp = fSamples.samples[i+1];
      
      Double_t smooth_avg = (currSamp+prevSamp+nextSamp)/3.0;
      smooth_samp[i] = smooth_avg;
    }
    
    // Look for local maxima instead of just threshold crossings
    Double_t ThresVal=GetThres(); // mV
    //max betwen NPedsum-1 and 1 since we want to compare to a previous sample as well
    //should it be just NPedsum instead of subtracting 1?
    //UInt_t samp=TMath::Max(NPedsum-1,1);
    UInt_t samp = 2;
    std::vector<Int_t> LocalMaxBinList;
    std::vector<Int_t> ThresCrossBinList;
    //std::vector<Int_t> FirstMaxPerCross;
    Bool_t thres_cross_check = kFALSE;
    Bool_t local_max_check = kFALSE;

    std::vector<Double_t> DeconSamps;
    DeconSamps.push_back(0.0);
    DeconSamps.push_back(0.0);
    
    //while ( LocalMaxBin < vals.size() - 1) {
    
    Double_t ns_per_samp = 4.0;
    Double_t tau = 4.0/GetExamplePulseParams();
    Double_t period_tau_ratio = ns_per_samp/tau;

    Double_t decon_weight[3];
    
    decon_weight[0] = TMath::Exp(-1.0 + period_tau_ratio)/period_tau_ratio;
    decon_weight[1] = -2.0/(period_tau_ratio*TMath::Exp(1.0));
    decon_weight[2] = TMath::Exp(-1.0 - period_tau_ratio)/period_tau_ratio;

    Double_t sumsig2_decon = 0.0;
    Double_t Vpp = 2.0; //2V peak-to-peak
    Int_t nbits_ADC = 12;

    for(Int_t w = 0; w < 3; w++){
      sumsig2_decon += pow( decon_weight[w] * Vpp/pow(2,nbits_ADC), 2 );
    }

    Double_t sigma_samp_decon = sqrt(sumsig2_decon);
    Double_t thres_nsigma_decon = 3.0;
    Double_t thres_decon = thres_nsigma_decon*sigma_samp_decon * 1000; //1000 for V to mV conv.

    //std::cout << " thres_decon is " << thres_decon << std::endl;
    
    while( samp < vals.size() - 1 ){
      Double_t decon_samp = 0.0;
      for(Int_t decon_range = 0; decon_range < 3; decon_range++){
	decon_samp += decon_weight[decon_range]*(fSamples.samples_raw[samp - decon_range] - fSamples.ped);
      } 
      DeconSamps.push_back(decon_samp);
      samp++;
    }

    
    for(Int_t decon_ind = 1; decon_ind < vals.size()-1; decon_ind++){

      Double_t currSamp = DeconSamps[decon_ind], prevSamp = DeconSamps[decon_ind-1], nextSamp = DeconSamps[decon_ind+1];

      if( currSamp > prevSamp && currSamp > nextSamp && (currSamp > thres_decon || std::max(currSamp+prevSamp,currSamp+nextSamp) >= sqrt(2.0)*thres_decon) ){
	ThresCrossBinList.push_back(decon_ind);
	LocalMaxBinList.push_back(decon_ind+5);
      }

    }
    
    //loop over maxima and do peak erasing, taken from SBSGEMModule::find_clusters_1D_experimental

    /*
    
    std::vector<Int_t> peaks_to_erase;

    UInt_t n_peaks_erased = 0;
    
    for( auto maxima_bin_iter : LocalMaxBinList ){

      Double_t ADCmax = fSamples.samples_raw[maxima_bin_iter];
      Double_t prominence = ADCmax;
      Int_t samplo = maxima_bin_iter, samphi = maxima_bin_iter;
      Double_t ADCminright = ADCmax, ADCminleft = ADCmax;

      Bool_t higherpeakleft = false, higherpeakright = false;
      Int_t peakright = -1, peakleft = -1;

      //sigma_max should be based on the individual ch threshold which is defined as something like 3x the pedrms for that ch (atleast 3x for ecal, need to check hcal)
      Double_t sigma_max = ThresVal/3.0;

      //search for potentially overlapping pulses to the right of current local max
      //while(std::find(fSamples.samples_raw.begin(), fSamples.samples_raw.end(), samphi + 1) != fSamples.samples_raw.end()){
      while(samphi < fSamples.samples_raw.size() - 1){
	//while( fSamples.samples_raw.find( samphi + 1 ) != fSamples.samples_raw.end() ){
	
	samphi++;
	Double_t ADCtest = fSamples.samples_raw[samphi];

	//want to find lowest sample before next peak
	if( ADCtest < ADCminright && !higherpeakright ){
	  ADCminright = ADCtest;
	}

	//if this eval true then we have found a higher peak to the right of this local maxima so we will want to erase it if it doesnt pass prominence and time sep threshold
	//if( LocalMaxBinList.find(samphi) != LocalMaxBinList.end() && ADCtest > ADCmax ){
	if( std::find(LocalMaxBinList.begin(), LocalMaxBinList.end(), samphi) != LocalMaxBinList.end() && ADCtest > ADCmax ){
	  higherpeakright = true;
	  peakright = samphi;
	  break;
	}

      }

      //search for potentially overlapping pulses to the left of current local max
      //while( fSamples.samples_raw.find( samplo - 1 ) != fSamples.samples_raw.end() ){
      //while( std::find(fSamples.samples_raw.begin(), fSamples.samples_raw.end(), samplo - 1) != fSamples.samples_raw.end() ){
      while(samplo > 0){
	
	samplo--;
	Double_t ADCtest = fSamples.samples_raw[samplo];

	//want to find lowest sample before next peak
	if( ADCtest < ADCminleft && !higherpeakleft ){
	  ADCminleft = ADCtest;
	}

 	//if this eval true then we have found a higher peak to the left of this local maxima so we will want to erase it if it doesnt pass prominence threshold
	if( std::find(LocalMaxBinList.begin(), LocalMaxBinList.end(), samplo) != LocalMaxBinList.end() && ADCtest > ADCmax ){
	  //if( LocalMaxBinList.find(samplo) != LocalMaxBinList.end() && ADCtest > ADCmax ){
	  higherpeakleft = true;
	  peakleft = samplo;
	  break;
	}

    }

      Bool_t peak_close = false;
      if( !higherpeakleft ) ADCminleft = 0.0;
      if( !higherpeakright ) ADCminright = 0.0;

      Double_t sigma = ThresVal/3.0;
      Double_t prom_thres_sigma_factor = 1.0;

      UInt_t max_peak_time_sep = 6;
 
      //prominence is vertical diff btwn peaks, adc val diff
      //max separation here will be time difference btwn peaks, horizontal
      if( higherpeakright || higherpeakleft ){
	
	prominence = ADCmax - std::max ( ADCminleft, ADCminright );// define prom as diff btwn current max and highest valley to left or right


	//make right peak thres smaller since rise time for calo pulse is so fast compared to decay
	if( higherpeakleft && std::abs( peakleft - maxima_bin_iter ) <= 2*max_peak_time_sep ) peak_close = true;
	if( higherpeakright && std::abs( peakright - maxima_bin_iter ) <= 1*max_peak_time_sep ) peak_close = true;

	//for now just defining prom threshold, include fraction of current peak later
	if( peak_close && prominence < prom_thres_sigma_factor*sigma ){
	  n_peaks_erased++;
	  peaks_to_erase.push_back( maxima_bin_iter );
	}
      }

    }//end of loop over local max

    Int_t first_peak_bin;
    Bool_t thres_check = kFALSE;

    std::vector<Int_t> thres_to_erase;
    std::vector<Int_t> thres_to_insert;
    
    for( Int_t thres = 0; thres < ThresCrossBinList.size(); thres++){

      Int_t next_thres = 0;
      if(thres == ThresCrossBinList.size() - 1){
	next_thres = vals.size() - 1;
      }else{
	next_thres = ThresCrossBinList[thres+1];
      }

      Int_t thres_to_erase_temp = 0;
      Int_t thres_to_insert_temp = 0;
      Bool_t dont_erase_found = kFALSE;
	
      for( Int_t max = 0; max < LocalMaxBinList.size(); max++){

	Int_t next_max = 0;
	if(max == LocalMaxBinList.size() - 1){
	  next_max = vals.size() - 1;
	}else{
	  next_max = LocalMaxBinList[thres+1];
	}
	
	if(dont_erase_found) break;
	if(LocalMaxBinList[max] > ThresCrossBinList[thres] && LocalMaxBinList[max] < next_thres ){
	  
	  auto loc_max_ind_to_erase = std::find(peaks_to_erase.begin(), peaks_to_erase.end(), LocalMaxBinList[max]);
	  if(loc_max_ind_to_erase == peaks_to_erase.end()){
	    dont_erase_found = kTRUE;
	  }else{
	    thres_to_erase_temp = ThresCrossBinList[thres];
	    if(next_max < next_thres) thres_to_insert_temp = next_max - 5;
	  }
	  
	}	    
      }

      //std::cout << "erase this one " << thres_to_erase_temp << " and insert " << thres_to_insert_temp << std::endl;

      
      if(thres_to_erase_temp != 0) thres_to_erase.push_back(thres_to_erase_temp);
      if(thres_to_insert_temp != 0) thres_to_insert.push_back(thres_to_insert_temp);
      
    }
    
    //erase non-prominent or insignificant peaks
    for( Int_t peak_ind = 0; peak_ind < peaks_to_erase.size(); peak_ind++ ){
      Int_t peak_val = peaks_to_erase[peak_ind];

      std::cout << peak_val << " was erased " << std::endl;
      
      auto loc_max_ind_to_erase = std::find(LocalMaxBinList.begin(), LocalMaxBinList.end(), peak_val);
      //if(loc_max_ind_to_erase != LocalMaxBinList.end()) LocalMaxBinList.erase(loc_max_ind_to_erase);
      
    }
          
    //delete threshold crossings corresponding to now erased peaks
    for( Int_t thres_ind = 0; thres_ind < thres_to_erase.size(); thres_ind++ ){
      Int_t thres_val = thres_to_erase[thres_ind];

      //std::cout << thres_val << " thres was erased" << std::endl;
 
      auto thres_ind_to_erase = std::find(ThresCrossBinList.begin(), ThresCrossBinList.end(), thres_val);
      //if(thres_ind_to_erase != ThresCrossBinList.end()) ThresCrossBinList.erase(thres_ind_to_erase);
    }
    
    //insert new threshold crossings to account for new first peaks after a previously existing threshold crossing
    for( Int_t thres_ind = 0; thres_ind < thres_to_insert.size(); thres_ind++ ){
      Int_t new_thres_val = thres_to_insert[thres_ind];
      
      for(Int_t cross_ind = 0; cross_ind < ThresCrossBinList.size(); cross_ind++){

	Int_t next_thres = 0;
	if(thres_ind == ThresCrossBinList.size() - 1){
	  next_thres = vals.size() - 1;
	}else{
	  next_thres = ThresCrossBinList[thres_ind+1];
	}	
	
	if(thres_to_insert[thres_ind] > ThresCrossBinList[cross_ind] && thres_to_insert[thres_ind] < next_thres){
	  //ThresCrossBinList.insert(ThresCrossBinList.begin() + cross_ind, thres_to_insert[thres_ind]);
	  //std::cout << new_thres_val << " thres was inserted aft and bef: " << ThresCrossBinList[cross_ind] << " " << next_thres << std::endl;
	  break;
	}
      }
    }
    */
    fSamples.multipulse.resize(LocalMaxBinList.size());
    
    //now, similar to the SBSGEMModule method, we want to split remaining significant pulses and split any overlapping pulses contributions
    
    if(LocalMaxBinList.size() == 0) return;

    Double_t scaled_pulse_x[vals.size()] = {0};
    
    Double_t waveform_pulse_x[vals.size()] = {0};
    Double_t waveform_pulse_y[vals.size()] = {0};
    Double_t waveform_pulse_y_smooth[vals.size()] = {0};

    Double_t DeconSampsArr[vals.size()] = {0};
    
    for(UInt_t jj = 0; jj < vals.size(); jj++){
      scaled_pulse_x[jj] = jj;
      waveform_pulse_x[jj] = jj;
      waveform_pulse_y[jj] = fSamples.samples_raw[jj]-fSamples.ped;
      waveform_pulse_y_smooth[jj] = smooth_samp[jj];
      DeconSampsArr[jj] = DeconSamps[jj];
    }
    
    for(Int_t ncross = 0; ncross < ThresCrossBinList.size(); ncross++){
    //for(Int_t ncross = 0; ncross < 1; ncross++){

      UInt_t NSB = GetNSB();
      UInt_t NSA = GetNSA();
      UInt_t FixedThresCrossBin=GetFixThresBin();
      Double_t FineTime = 0;
      Double_t max  = 0;
      Double_t sum = 0;
      Double_t sped = 0;

      // convert to pC, assume tcal is in ns, and 50ohm resistance
      Double_t pC_Conv = fSamples.tcal/50.;
      
      Bool_t PeakFound= kFALSE;
      UInt_t PeakBin= 0;
      UInt_t IntMinBin= 0;
      UInt_t IntMaxBin= vals.size();
	  
      IntMinBin= TMath::Max(ThresCrossBinList[ncross]-NSB,IntMinBin);
      IntMaxBin= TMath::Min(ThresCrossBinList[ncross]+NSA-1,IntMaxBin);

      Double_t scaled_pulse_y[vals.size()] = {0};
      
      Double_t template_pulse_param = GetExamplePulseParams();
      template_pulse_scaled = new TF1("template_pulse_scaled","[0] + [1] * (x-[2]) * exp((x-[2])*(-[3]))");

      /*
      Double_t FineCrossTime = 0.0;
      Double_t VMid = DeconSamps[ThresCrossBinList[ncross]]/2.0;
      for(size_t i =IntMinBin ; i <LocalMaxBinList[ncross]+1 ; i++ ) {
	if (VMid >= DeconSamps[i]  && VMid < DeconSamps[i+1]) {
	  FineCrossTime = i+(VMid-DeconSamps[i])/(DeconSamps[i+1]-DeconSamps[i]);
	}
      }

      std::cout << ThresCrossBinList[ncross] << " and decon is " << DeconSamps[ThresCrossBinList[ncross]] << std::endl;
      */

      Double_t amp_from_decon = 0.0;
      Double_t time_from_decon = 0.0;
      for(Int_t decon_range = -2; decon_range <=2; decon_range++){
	Double_t decon_bin = ThresCrossBinList[ncross]+decon_range;
	Double_t decon_volt = DeconSamps[decon_bin];
	amp_from_decon += decon_volt;
	time_from_decon += decon_volt*decon_bin;
      }
      time_from_decon /= amp_from_decon;

      //std::cout << " time from decon calc is " << time_from_decon << std::endl;

      if(DeconSamps[ThresCrossBinList[ncross]] < thres_decon){
	delete template_pulse_scaled;
	continue;
      }
      
      template_pulse_scaled->FixParameter(0, 0.0);
      template_pulse_scaled->FixParameter(1, ( amp_from_decon ) * TMath::Exp(1) * template_pulse_param);
      template_pulse_scaled->FixParameter(2, Double_t(ThresCrossBinList[ncross]) - 1);
	  //template_pulse_scaled->FixParameter(2, FineCrossTime);
	  //std::cout << "old thres " << new_thres_cross << " shift is " << shift_par << std::endl;
      template_pulse_scaled->FixParameter(3, template_pulse_param);
	  
      
      for(Int_t fill_iter = IntMinBin; fill_iter < IntMaxBin; fill_iter++){
	Double_t template_pulse_val = template_pulse_scaled->Eval(fill_iter);
	scaled_pulse_y[fill_iter] = template_pulse_val;
	
	sped+=fSamples.ped*pC_Conv;
	
	if(template_pulse_val <= fSamples.samples_raw[fill_iter] - fSamples.ped){
	  sum += template_pulse_val;
	}else{
	  sum += fSamples.samples_raw[fill_iter] - fSamples.ped;
	}
	
      }

      fSamples.multipulse[ncross].integral.raw = sum+sped;
      fSamples.multipulse[ncross].integral.val = (sum)*fSamples.cal;
      fSamples.multipulse[ncross].time.raw = time_from_decon;
      fSamples.multipulse[ncross].time.val = (time_from_decon)*fSamples.tcal + fSamples.timeoffset;
      fSamples.multipulse[ncross].amplitude.raw = amp_from_decon+fSamples.ped;
      fSamples.multipulse[ncross].amplitude.val = amp_from_decon*fSamples.acal;
      if (amp_from_decon==0) fSamples.multipulse[ncross].amplitude.val=amp_from_decon;

      //std::cout << time_from_decon << " and calib " << (time_from_decon)*fSamples.tcal + fSamples.timeoffset << std::endl;
      /*
      auto waveform_plot = new TGraph(vals.size(), waveform_pulse_x, waveform_pulse_y);
      auto waveform_plot_smooth = new TGraph(vals.size(), waveform_pulse_x, waveform_pulse_y_smooth);
      auto scaled_waveform_plot = new TGraph(vals.size(), scaled_pulse_x, scaled_pulse_y);
      auto decon_plot = new TGraph(vals.size(), waveform_pulse_x, DeconSampsArr);
      
      if(ncross == 0){
      waveform_plot->Draw("ALP");
      waveform_plot->SetMarkerColor(kBlue+1);
      waveform_plot->SetLineWidth(3);
      
      decon_plot->SetLineColor(kRed+1);
      decon_plot->SetLineWidth(1);
      decon_plot->SetMarkerStyle(20); // Full circle
      decon_plot->SetMarkerSize(1); // Adjust size
      decon_plot->SetMarkerColor(kRed+1);
      decon_plot->Draw("SAME P");
      }

      scaled_waveform_plot->SetLineColor(kMagenta+2);
      scaled_waveform_plot->SetLineWidth(4);
      scaled_waveform_plot->Draw("SAME");
  
      gPad->Modified();
      gPad->Update();
      gPad->WaitPrimitive();
      */
      delete template_pulse_scaled;
      
    }
    /*
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
  
      Double_t VMid = (max+fSamples.ped)/2.;
      if (PeakFound) {
	for(size_t i =IntMinBin ; i <PeakBin+1 ; i++ ) {
	  if (VMid >= fSamples.samples_raw[i]  && VMid < fSamples.samples_raw[i+1]) {
	    FineTime = i+(VMid-fSamples.samples_raw[i])/(fSamples.samples_raw[i+1]-fSamples.samples_raw[i]);
	    
	  }
	}
      }
      
      
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
    */
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
    
    // Look for multiple threshold crossings
    Double_t ThresVal=GetThres(); // mV
    UInt_t ThresCrossBin=TMath::Max(NPedsum-1,0);
    std::vector<UInt_t> ThresCrossBinList;
    //    std::cout << " ped = " << fSamples.ped << " thres = " << ThresVal << std::endl ;
    while ( ThresCrossBin < vals.size() ) {
      //look for threshold crossings
      if( fSamples.samples_raw[ThresCrossBin] >= fSamples.ped+ThresVal ){
	//incase first sample is above threshold (is this necessary?)
	if( ThresCrossBin == 0 ){
	  ThresCrossBinList.push_back(ThresCrossBin);
	//otherwise just check that the previous sample was actually below threshold so we dont just add all samples above threshold
	}else if( fSamples.samples_raw[ThresCrossBin-1] < fSamples.ped+ThresVal ){
	  ThresCrossBinList.push_back(ThresCrossBin);
	}
      }
      ThresCrossBin++;
    }
    //
 
    fSamples.multipulse.resize(ThresCrossBinList.size());
    
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
