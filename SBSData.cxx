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

  //alternate idea from Andrew: find all maxima in waveform and sort them by amplitude (I think maybe good time could be way to sort also), then start with "best" peak and look forward and back in time to see if we go back below threshold before encountering another maxima
  //instead of threshold crossing list, would have maxima list (ROOT has some sort of built in peak finder method), then loop over samples before and after "best" peak until you hit drop below threshold and search for other samples in the peak list

  
  /* //comment out starting here
  void Waveform::ProcessMultiMaxima(std::vector<Double_t> &vals)
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
    
    // Get all maxima
    Double_t ThresVal=GetThres(); // mV
    std::vector<UInt_t> MaximaBinList;

    //TSpectrum is outdated but no real replacement in ROOT yet, and not in current halla analyzer build so maybe just try smoothing the waveform (by averaging adjacent points and replacing) and finding peaks by finding where first derivative goes sharply negative enough in conjunction with original sample value in that vecinity being above a certain "peak amplitude threshold"
    //std::vector<Double_t> smooth_waveform;
    //for(int i = 0; 
    
    TSpectrum *wavespectrum = new TSpectrum();
    Double_t dest[vals.size()];
    Int_t MaximaFound = wavespectrum->SearchHighRes(fSamples.samples_raw, dest, vals.size(), 2, 10, kFALSE, 1, kFALSE, 1);
    Double_t *MaximaPos = wavespectrum->GetPositionX();

    //the list of maxima returned from the Search function of the TSpectrum should rank them starting from highest amplitude and decreasing
    for(Int_t maxima_int = 0; maxima_int < MaximaFound; maxima_int++){
      Int_t maxima_bin = Int_t(MaximaPos[maxima_int]+0.5);
      if(fsamples.samples_raw[maxima_bin] >= fSamples.ped + ThresVal){
	MaximaBinList.push_back(maxima_bin);
      }
    }
 
    fSamples.multipulse.resize(MaximaBinList.size());
    
    //Loop over all maxima, perhaps now we should look to the left and right of the maximum in question until we find a sample going below threshold and if we find another maximum before going below threshold then we attempt to resolve the two pulses.
    for(int nmaxima = 0; nmaxima < MaximaBinList.size(); nmaxima++){

      UInt_t NSB = GetNSB();
      UInt_t NSA = GetNSA();
      UInt_t FixedThresCrossBin=GetFixThresBin();
      Double_t FineTime = 0;
      Double_t max  = 0;
      Double_t sumleft = 0, sumright = 0;
      Double_t spedleft = 0, spedright = 0;
      
      Bool_t PeakFound= kFALSE;
      UInt_t PeakBin= 0;
      UInt_t IntMinBin= 0;
      UInt_t IntMaxBin= vals.size();
      
      if (ThresCrossBinList[ncross] < vals.size()) {
	IntMinBin= TMath::Max(ThresCrossBinList[ncross]-NSB,IntMinBin);
	IntMaxBin= TMath::Min(ThresCrossBinList[ncross]+NSA-1,IntMaxBin);
      } else {
	IntMinBin = TMath::Max(FixedThresCrossBin-NSB,IntMinBin);
	IntMaxBin = TMath::Min(FixedThresCrossBin+NSA-1,IntMaxBin);
      }
      // convert to pC, assume tcal is in ns, and 50ohm resistance
      Double_t pC_Conv = fSamples.tcal/50.;

      Bool_t LeftThresFound = kFALSE, RightThresFound = kFALSE;
      Bool_t LeftOverLapFound = kFALSE, RightOverLapFound = kFALSE;
      UInt_t SampInc = 1;
      UInt_t SampsFromStartToPeak = 0;
      while(LeftThresFound == kFALSE && RightThresFound == kFALSE){
	UInt_t LeftInd = MaximaBinList[nmaxima]-SampInc;
	UInt_t RightInd = MaximaBinList[nmaxima]+SampInc;
	Double_t SampLeft = fSamples.samples_raw[LeftInd];
	Double_t SampRight = fSamples.samples_raw[RightInd];
	
	//if we encounter another maxima, we want to recalculate the sum and sped on whichever side of the maxima we are currently on, as well as go ahead and calculate the sum and sped for the opposing side of the other maxima which we encountered since we will need to split up the sample values between the two overlapping pulses
	
	if(LeftThresFound == kFALSE && LeftOverLapFound == kFALSE){
	  //check that we havent run into another maxima
	  for(Int_t MaximaCheck = 0; MaximaCheck < MaximaBinList.size(); MaximaCheck++){
	    if(MaximaBinList[MaximaCheck] == LeftInd){
	      LeftOverLapFound = kTRUE;
	      break;
	    }
	  }
	  //if we didnt run into maxima, gather sample info and sum into pulse info
	  if(SampLeft >= fSamples.ped + ThresVal){
	    spedleft+=(fSamples.ped*pC_Conv);
	    sumleft+=(SampLeft*pC_Conv);
	    SampsFromStartToPeak++;
	  }else{
	    LeftThresFound = kTRUE;
	  }
	}

	if(RightThresFound == kFALSE && RightOverLapFound == kFALSE){
	  //check that we havent run into another maxima
	  for(Int_t MaximaCheck = 0; MaximaCheck < MaximaBinList.size(); MaximaCheck++){
	    if(MaximaBinList[MaximaCheck] == RightInd){
	      RightOverLapFound = kTRUE;
	      spedright = 0.0;
	      sumright = 0.0;
	      break;
	    }
	  }
	  if(SampRight >= fSamples.ped + ThresVal){
	    spedright+=(fSamples.ped*pC_Conv);
	    sumright+=(SampRight*pC_Conv);
	  }else{
	    RightThresFound = kTRUE;
	  }
	}
	
	SampInc++;
	
      }

      //Start out with just splitting up the pulse based on the clean elastic "reference" pulse for this channel to find what portion of the samples should contribute to the current pulse
      //I guess the contribution to the adjacent pulse encountered should just be whats left over not contributing to the current pulse
      
      //Want to find first pulse in time such that we have a clean start for the rise and hopefully a good peak amplitude assuming the adjacent pulse isnt too close in time in which case the rise of the adjacent pulse may contribute to the peak amplitude of the first pulse
      //Assuming we have this good first pulse, we should fit the first pulse with our clean "reference" pulse for the given channel and use the fit to determine how much of the overlapping samples should contribute to the tail of the first pulse and then cascade in time from there

      if(RightOverLapFound == kTRUE){
	//starting assumption: if we can find a clean start to the rise of the first pulse in the waveform and its peak sample, then we can scale the "reference" pulse linearly to fit the shape determined by the starting point and peak of the pulse

	//find ratio of number of samps between start and peak of "reference" pulse and number of samps between start and peak of the actual pulse here to find horizontal scale factor
	//find ratio of peak amplitudes between "reference" pulse and actual pulse here to find vertical scale factor
	//then use these factors to scale the "reference" pulse to match the actual pulse here and then find the amplitude at the time-values (x-axis) coming after the peak sample

	//this is currently giving the sample values of the template pulses
	std::vector<Double_t> template_pulse = GetExamplePulseParams();
	
	TSpline3 *cubic_spline = new TSpline3("cubic_spline",0,template_pulse.size(),template_pulse.data(), template_pulse.size());

	Int_t template_dist_to_peak = std::max_element(template_pulse.begin(), template_pulse.end());
	Double_t template_peak = template_pulse[template_dist_to_peak];
	Double_t real_pulse_peak = fSamples.samples_raw[MaximaBinList[nmaxima]];

	Double_t h_scale = SampsFromStartToPeak / template_dist_to_peak;
	Double_t v_scale = real_pulse_peak / template_peak;

	Int_t max_samps_after_peak = 5;

	for(int samp_after_peak = SampsFromStartToPeak; samp_after_peak < max_samps_after_peak; samp_after_peak++){
	  Double_t scaled_samp = v_scale*cubic_spline->Eval(samp_after_peak*h_scale);
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
    fHasData = (vals.size() > 0);
  }
  */

  
  //filter function outline:
  
  //loop over all samples in noisey-wf excluding the ouside regions, like first 10 and last 10 samples
  // this will depend on the number of samples surrounding the peak of the ref-wf

  //in this loop over wf samples, loop over ref-wf
  //dot product: at the current sample of the noisey-wf, go left by the number of samples away from peak in ref-wf, and multiply the amp of the noisey-wf by the ref-wf, and then move over to the right until youve gone to the right-most end of the ref-wf. Sum up these multiplications and assign the summation to the matched filter function output value at the current sample in the noisey-wf. This dot product is the correlation btwn the noisey-wf and the ref-wf. NPS implements this as a convolution by time reversing the ref-wf (starting at the right most end of the ref-wf rather than the left most end)
  
  void Waveform::ProcessMulti(std::vector<Double_t> &vals)
  {
    UInt_t OverlapCount_full_wave = 0;
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

    // Get all maxima
    Double_t ThresVal=GetThres(); // mV

    //maximabinlist stuff EXPERIMENTAL
    std::vector<UInt_t> MaximaBinList;
    /*
    TSpectrum wavespectrum;
    Double_t dest[vals.size()];
    Int_t MaximaFound = wavespectrum.SearchHighRes(fSamples.samples_raw.data(), dest, vals.size(), 2, 10, kFALSE, 1, kFALSE, 1);
    Double_t *MaximaPos = wavespectrum.GetPositionX();

    //the list of maxima returned from the Search function of the TSpectrum should rank them starting from highest amplitude and decreasing
    for(Int_t maxima_int = 0; maxima_int < MaximaFound; maxima_int++){
      Int_t maxima_bin = Int_t(MaximaPos[maxima_int]+0.5);
      if(fSamples.samples_raw[maxima_bin] >= fSamples.ped + ThresVal){
	MaximaBinList.push_back(maxima_bin);
      }
    }
    */

    //thres_cross_check stuff is EXPERIMENTAL
    
    // Look for multiple threshold crossings
    UInt_t ThresCrossBin=TMath::Max(NPedsum-1,0);
    std::vector<UInt_t> ThresCrossBinList;
    Bool_t thres_cross_check = kFALSE;
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
	  thres_cross_check = kTRUE;
	}
	//EXPERIMENTAL
	/*
	//check for more than one peak being found after the threshold crossing before we go back below threshold
	if(thres_cross_check == kTRUE){
	  UInt_t n_max_found = 0;
	  for(UInt_t max_check = 0; max_check < MaximaFound; max_check++){
	    if( fSamples.samples_raw[ThresCrossBin] == MaximaBinList[max_check] ) n_max_found++;
	  }
	  if(n_max_found > 1){
	    fSamples.n_overlaps++;
	  }
	}
	*/
      }else{
	thres_cross_check = kFALSE;
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
      Bool_t OverlapFound= kFALSE;
      
      if (ThresCrossBinList[ncross] < vals.size()) {
	IntMinBin= TMath::Max(ThresCrossBinList[ncross]-NSB,IntMinBin);
	IntMaxBin= TMath::Min(ThresCrossBinList[ncross]+NSA-1,IntMaxBin);
      } else {
	IntMinBin = TMath::Max(FixedThresCrossBin-NSB,IntMinBin);
	IntMaxBin = TMath::Min(FixedThresCrossBin+NSA-1,IntMaxBin);
      }
      // convert to pC, assume tcal is in ns, and 50ohm resistance
      Double_t pC_Conv = fSamples.tcal/50.;
      
      for(size_t i =IntMinBin ; i <IntMaxBin ; i++ ) {
	sped+=fSamples.ped*pC_Conv;
	//accumulate sample sums in range of few samples before and after the threshold crossing (range determined by IntMinBin and IntMaxBin)
	sum+=fSamples.samples_raw[i]*pC_Conv;
	//if ( i >= ThresCrossBinList[ncross] && !PeakFound) {
	if ( i >= ThresCrossBinList[ncross] ) {
	  
	  if(fSamples.samples_raw[i] > max && !PeakFound) { //check current sample higher than max if peak hasnt already been found
	    max = fSamples.samples_raw[i];

	    //this check is only for overlapping pulses EXPERIMENTAL
	  } else if (PeakFound && i > PeakBin + 5 &&  fSamples.samples_raw[i] > fSamples.samples_raw[i-1] && fSamples.samples_raw[i-1] > fSamples.samples_raw[i-2]){ //if peak exists already, check 5 samples later if amplitude is increasing consistently
	    //the use of 5 samples later here is somewhat arbitrary but this should be defined as the "prominence" or the distance after finding a peak which we want to search for another potential overlapping pulse
	    OverlapFound = kTRUE;
	  } else if (!PeakFound && fSamples.samples_raw[i+1] < fSamples.samples_raw[i]){ //check that two samples after current max(at i and i+1, where peak is at i-1) are continuously decreasing, if so the current max is the peak
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
      //EXPERIMENTAL
      if(OverlapFound){
	OverlapCount_full_wave++;
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
    fSamples.n_overlaps = OverlapCount_full_wave;
  }
  
  void Waveform::Clear()
  {
    for(size_t i = 0; i < fSamples.samples.size(); i++) {
      fSamples.samples_raw[i] = fSamples.samples[i] = 0;
    }
    fHasData = false;
  }

}; // end SBSData
