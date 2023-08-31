// Author: AuroraPepino

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// GMCTGenerateWaveforms                                                      //
//                                                                            //
// Begin_Html <!--
/*-->

Description:
<p>
Generate waveforms in DCH
</p>

Usage:
<p>
No special features
</p>

Status:
<p>
In progress
</p>

<!--*/
// --> End_Html
//                                                                            //
// The event methods have been written by AuroraPepino.                       //
//                                                                            //
// Please note: The following information is only correct after executing     //
// the ROMEBuilder.                                                           //
//                                                                            //
// This task accesses the following folders :                                 //
//     DCHHit                                                                 //
//     WaveformData                                                           //
//                                                                            //
//                                                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

/* Generated header file containing necessary includes                        */
#include "generated/GMCTGenerateWaveformsGeneratedIncludes.h"

////////////////////////////////////////////////////////////////////////////////
/*  This header was generated by ROMEBuilder. Manual changes above the        *
 * following line will be lost next time ROMEBuilder is executed.             */
/////////////////////////////////////----///////////////////////////////////////

#include "TMath.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TRandom.h"
#include "TF1.h"
#include "TH1F.h"


#include "generated/GMCAnalyzer.h"
#include "tasks/GMCTGenerateWaveforms.h"
#include "util/XTRELTIME.h"
#include "ROMEiostream.h"
#include "generated/GMCWaveformData.h"
#include "generated/GMCDCHHit.h"
#include "TClonesArray.h"  
#include "GMCAnalyzer.h"
#include "generated/GMCTGenerateWaveforms_Base.h"
#include "util/Geometry.h"


// uncomment if you want to include headers of all folders
//#include "GMCAllFolders.h"

Double_t signalShape2(Double_t *x,Double_t *p);

const Double_t fDelayFromTrigger = 0.3;

ClassImp(GMCTGenerateWaveforms)

using namespace std;
//__________________________________________________________________________
void GMCTGenerateWaveforms::Init() {

  std::cout<<"inizializzazione waveforms"<<std::endl;
  // relazione x-t
  fxt=new XTRELTIME("XTData/X_T_Relation_8515.root");

  for (int i=0;i<6;i++) fSignalShapePar[i] = 0.;  

  fGeometry = Geometry::GetInstance();
  
  SetSignalParameters(); 
  
  InitPolyaFunction();

  //inizializzo shape2
  //  TF1 *shape = new TF1("shape",shapeSignal2,0.01,10.,7);     
}

//_________________________________________________________________________
void GMCTGenerateWaveforms::BeginOfRun() {


}

//_________________________________________________________________________
void GMCTGenerateWaveforms::Event() {

  TClonesArray *hits = gAnalyzer->GetDCHHits();
  //int n_hit=hits->GetEntries();
  cout<<"Numero totale di hit "<<hits->GetEntries()<<endl<<
    "___________________________________________________"<<endl;
  
//  Double_t *signR;   //segnale letto dalla parte destra
//  Double_t *signL;   //segnale letto dalla parte sinistra
  //  Double_t *time;   //vettore dei tempi               
                            
  //__________________________________________________________________________
     
  gAnalyzer->SetWaveformDataSize(hits->GetEntries());
  
  for (int j=0; j<hits->GetEntries(); j++){       
    
    GMCDCHHit *ah = (GMCDCHHit*)hits->At(j);
    Int_t num_cluster = ah->GetfNrClusters();//ritorna il num di cluster del singolo hit
    Int_t trk_id = ah->GetfTrkIDofHit();
    
    Int_t total_charge = (int)ah->GetfTotalCharge();
    Double_t impact_param = ah->GetfImpact();
    Double_t track_length = ah->GetfLength();
    cout <<"numero di cluster generati da questo hit "<<num_cluster<<endl<<
      " identificativo traccia di questo hit: " <<trk_id<<endl<<
      " con carica tot: "                       <<total_charge  <<" e-"<<endl<<
      " parametro di impatto: "                 <<impact_param  <<" mm"<<endl<<
      " lunghezza della traccia: "              <<track_length  <<" mm"<<endl;
    
    signR = new Double_t [fNrOfBins];
    signL = new Double_t [fNrOfBins];
    
    if (total_charge < 0) total_charge = 0;
    Double_t *peaks = new Double_t [total_charge];
    Double_t *timeseq = new Double_t [total_charge];

    Int_t idxOfPeaks[total_charge];
    Int_t ncntpk = 0;

    for (int l=0;l<fNrOfBins;l++) signR[l] = signL[l] = 0.; 
    for (int l=0;l<total_charge;l++) peaks[l] = timeseq[l]/* = peaksAmplitude[l]*/ = 0.;
    
    if (num_cluster > 0) {
      
      peaks = DefinePeakValues(ah);
      timeseq = DefineTimeValues(ah);
     
      for (int i=0; i<ah->GetfTotalCharge(); i++){

    	  //        CreateSignal(timeseq[i],10.*fVo*peaks[i]);
    	  Int_t idxpk = CreateSignal3(timeseq[i],fVo_1*peaks[i]);

    	  Bool_t found = false;
    	  int k = 0;
    	  while (k < ncntpk) {
    		  if (idxOfPeaks[k] == idxpk) found = true;
    		  k++;
    	  }

    	  if (!found) {
    		  idxOfPeaks[ncntpk] = idxpk;
    		  ncntpk++;
    	  }

      }

    } 
    cout<<"________________________________________________________________________"<<endl; 
    

    FEGain();
    
    if (gAnalyzer->GetGenerateWaveformsTaskBase()->GetSP()->GetNoiseOn()) AddNoise();
    
    
    if (gAnalyzer->GetGenerateWaveformsTaskBase()->GetSP()->GetDigitizationOn()) WaveDigitization();

    Double_t *peaksAmplitude = new Double_t [ncntpk];
    Int_t *peaksPos = new Int_t [ncntpk];
    for (int i=0; i<ncntpk; i++) {
      peaksPos[i] = idxOfPeaks[i];
      peaksAmplitude[i] = -signL[idxOfPeaks[i]];    
    }

    GMCWaveformData *wave_data = gAnalyzer->GetWaveformDataAt(j);
    wave_data->SetfCellId(ah->GetfCellId());
    wave_data->SetfNpoints(fNrOfBins);
    wave_data->SetfBinWidth(fTimeres);
    wave_data->SetfSignalL(signL);
    wave_data->SetfSignalR(signR);
    wave_data->SetfNofGenEle(total_charge);
    wave_data->SetfDTimeEle(timeseq);
    wave_data->SetfNrEleAvalanche(peaks);
    wave_data->SetfPeakSize(ncntpk);
    wave_data->SetfPeakBinPos(peaksPos);
    wave_data->SetfPeakAmplitude(peaksAmplitude);
  }
  
  //        TGraph *waveformR=new TGraph(Npoint,time,signR); 
  //        TCanvas *c1=new TCanvas();
  //        c1->cd();
  //        waveformR->Draw("ap");
  //        c1->SaveAs("eccomi2.gif");     
  //        tree->Write();
  //        my_waveform_file->Close();
  
}

//______________________________________________________________________________
void GMCTGenerateWaveforms::SetSignalParameters() {
  
  fTimeres = 0.001*gAnalyzer->GetGenerateWaveformsTaskBase()->GetSP()->GetBinTimeWidth(); //tempo in microsecondi
  
  fTimeWindow = 0.001*gAnalyzer->GetGenerateWaveformsTaskBase()->GetSP()->GetTimeWindow(); //tempo in microsecondi
  
  Double_t tauRUp = 0.001*gAnalyzer->GetGenerateWaveformsTaskBase()->GetSP()->GettauRumpUp(); //tempo in microsecondi
  
  Double_t tau = 0.001*gAnalyzer->GetGenerateWaveformsTaskBase()->GetSP()->GettauRumpDown();

  fSignalShapePar[2] = tau;
  fSignalShapePar[3] = tauRUp;
  fSignalShapePar[4] = gAnalyzer->GetGenerateWaveformsTaskBase()->GetSP()->GetFractionMixSign(); 
  fSignalShapePar[5] = 0.001*gAnalyzer->GetGenerateWaveformsTaskBase()->GetSP()->GettauRumpDown2(); 
  
  fConst1=fTimeres/tauRUp;
  fConst2=fTimeres/tau;
  
  Double_t Volt = gAnalyzer->GetGenerateWaveformsTaskBase()->GetSP()->GetHVTube();
  cout<<"V_filo: "<<Volt<<endl;
  
  Double_t mob = gAnalyzer->GetGenerateWaveformsTaskBase()->GetSP()->Getmobility();//(cm^2/us)*V^-1
  cout<<"mobilita: "<<mob<<endl;
  
  fGasGain = gAnalyzer->GetGenerateWaveformsTaskBase()->GetSP()->GetGasGain();
  cout<<"Gas gain: "<<fGasGain<<endl;
  
  Double_t ResistRR=gAnalyzer->GetGenerateWaveformsTaskBase()->GetSP()->GetResistRR();   
  cout<<"reisistenza di adattamento: "<<ResistRR<<endl;
  
  Double_t Ztube=fGeometry->GetTubeImpedance();
  cout<<"Ztube: "<<Ztube<<endl;
    
  Double_t Rtube=fGeometry->GetTubeSize();  //mm          
  Double_t Rwire=fGeometry->GetAnodeSize();  //mm
  cout<<"Tubo size (mm): "<<Rtube<<" Anodo size (mm): "<<Rwire<<endl;
    
  Double_t reflect = (ResistRR - Ztube)/(ResistRR + Ztube);
  Double_t logRat=TMath::Log(Rtube/Rwire);  
  Double_t ro=3.*Rwire*0.1*0.5; // distanza alla quale si forma la valanga (in cm)         
  Double_t to=pow(ro,2)*logRat/(2.*mob*Volt); //in sec
//  Double_t tMax = (TMath::Power(Rtube,2) - TMath::Power(Rwire,2))*logRat/
//    (2.*mob*Volt);
  cout<<"to: "<<to<<endl;        
  fVo = -0.5*TMath::Qe()/(2.*logRat*to)*ResistRR*(1.-reflect); // valore massimo del segnale
  fVo *= 1000.; // VO ESPRESSO in mV      
//  cout<<"Debug 1 fV0 = "<<fVo*fGasGain*1000.<<"  and to = "<<to*1.e6<<endl;  
  
  fVo_1 = (-0.5*TMath::Qe()*ResistRR*(1.-reflect))*1.e9; //10^9 per ottenere l'integrale in mV.us
  fNrOfBins =(int)(fTimeWindow/fTimeres);   //numero dei punti
  
  rumpUp = TMath::Nint(4.0*tauRUp/fTimeres);
  
  //define noise level
  fNoiseLevel = DefineNoiseLevel();
  
  fADCbin = TMath::Power(2.,gAnalyzer->GetGenerateWaveformsTaskBase()->GetSP()->GetADCResolution());
  fDigitStep  = gAnalyzer->GetGenerateWaveformsTaskBase()->GetSP()->GetADCFullRange();
  fDigitStep /= (double)fADCbin;
  
}

//______________________________________________________________________________
Double_t GMCTGenerateWaveforms::DefineNoiseLevel() {

  //define a noise level as fuction of Signal/noise ratio (mV)
  fSignalShapePar[0] = fVo_1*fGasGain;  
  fSignalShapePar[1] = 0.1;  //0.48
  Double_t tmpNoise = 999.;
  Int_t ii = 0;
  Double_t tmpx = 0.1;
  while (tmpNoise > signalShape2(&tmpx,fSignalShapePar)) {
    ii++;
    tmpx += ii*fTimeres;
    tmpNoise = signalShape2(&tmpx,fSignalShapePar);
  }

  //  tmpNoise *= gAnalyzer->GetGenerateWaveformsTaskBase()->GetSP()->GetFrontEndGain(); 
  tmpNoise /= gAnalyzer->GetGenerateWaveformsTaskBase()->GetSP()->GetSignalNoiseRatio();

  return tmpNoise;  

}

//______________________________________________________________________________
void GMCTGenerateWaveforms::InitPolyaFunction() {

    fPolya = new TF1("Polya","((1.+[2])/([1]*TMath::Gamma([2]+1)))*TMath::Power(x*(1.+ [2])/[1],[2])*TMath::Exp(-x*(1.+ [2])/[1])",1.e3,1.e7);


  fPolya->SetParameter(0,1.);          //Normalization
  fPolya->SetParameter(1,fGasGain);    //GasGain
  fPolya->SetParameter(2,0.5);         //Geometrical parameter for Drift Tube 
                                       //(see G.D. Alkhazov NIM A89 pag 155)
}

//______________________________________________________________________________
Double_t* GMCTGenerateWaveforms::DefinePeakValues(GMCDCHHit *ah) {

  //define the avalanche charge taking in account the statistic fluctuation by a Polya distribution

  Double_t  *pksVal = new Double_t [(int)ah->GetfTotalCharge()];
  Int_t nCls = 0; 
  for (int icls=0;icls<ah->GetfNrClusters();icls++) {

     Double_t carica_cluster = ah->GetfClsChargeAt(icls);
     for (int ichg=0;ichg<(int)carica_cluster;ichg++) {
     
        Double_t charge = fPolya->GetRandom();
        pksVal[nCls+ichg] = charge;    //10 -- is the FE gain
     }

     nCls += (int)carica_cluster;
  }
  
  return pksVal;
  
}

//______________________________________________________________________________
Double_t* GMCTGenerateWaveforms::DefineTimeValues(GMCDCHHit *ah) {

  //define the peak time: a gaussian smearing using longitudinal diffusion
  //is used

  Double_t  *timeVal = new Double_t [(int)ah->GetfTotalCharge()];
  Int_t nCls = 0; 
  for (int icls=0;icls<ah->GetfNrClusters();icls++) {

     Double_t carica_cluster = ah->GetfClsChargeAt(icls);
     Double_t distanza_drift = ah->GetfClsDistanceAt(icls);

     Float_t* x2time = fxt->spaceToTime((float)0.1*distanza_drift);
     
     //correct the right time for long distance
     if (distanza_drift >= 4.) x2time[0] = (3.094-TMath::Sqrt(3.094*3.094 + 
                                            4.*3.438*(0.0345-0.1*distanza_drift)))/(2.*3.438);
                                            
     for (int ichg=0;ichg<(int)carica_cluster;ichg++) {
     
        timeVal[nCls+ichg] = x2time[0] + gRandom->Gaus(0.,x2time[1]);
     }

     nCls += (int)carica_cluster;
  }
  
  return timeVal;
  
}

//______________________________________________________________________________
void GMCTGenerateWaveforms::CreateSignal(Double_t time, Double_t pValue) {

  //create an electrical signalfor each peak

  Int_t binCharge = (Int_t)((fDelayFromTrigger + time)/fTimeres);       

  if (binCharge >= fNrOfBins) return;  //condizione di controllo

  cout<<"bin charge: "<<binCharge<<"   "<<rumpUp<<endl;

  Float_t scoeff = 1;
  Int_t kBin = 0;   
  for (; kBin<=rumpUp; kBin++){  

    //PRIMA PARTE DEL SEGNALE	
    scoeff = (1.-TMath::Exp(-1.*((Float_t)kBin)*fConst1));
    
    signR[binCharge] += (pValue*scoeff);  // TMath::Log(1 + i*timeres/to);
    signL[binCharge] += (pValue*scoeff);               
    
    binCharge++;
  }

  kBin=0;
  pValue *= scoeff;
  while(binCharge < fNrOfBins) {
    
    scoeff = TMath::Exp(-1.*((Float_t)kBin+1.)*fConst2);
    signR[binCharge] += pValue*scoeff;     //tau/(to + i*timeres);
    signL[binCharge] += pValue*scoeff;     //tau/(to + i*timeres);
    binCharge++;
 
    kBin++;
  }
  
}

//______________________________________________________________________________
void GMCTGenerateWaveforms::CreateSignal2(Double_t time, Double_t pValue) {

  //scaling the pvalue from Voltage to Charge

  fSignalShapePar[0] = pValue;  
  fSignalShapePar[1] = fDelayFromTrigger+time;  //0.482

  Int_t binCharge = (Int_t)((fDelayFromTrigger + time)/fTimeres);       

  if (binCharge >= fNrOfBins) return;  //condizione di controllo

  Int_t kBin = binCharge;
  Double_t xtmp;   
  Double_t assetime[fNrOfBins];
  for (int iii=0;iii<fNrOfBins;iii++) assetime[iii] = iii*fTimeres;

  for (; kBin<fNrOfBins; kBin++){  

    xtmp = kBin*fTimeres;

    signR[kBin] += signalShape2(&xtmp,fSignalShapePar);
    signL[kBin] += signalShape2(&xtmp,fSignalShapePar);
  }               

}

//______________________________________________________________________________
Int_t GMCTGenerateWaveforms::CreateSignal3(Double_t time, Double_t pValue) {

  //scaling the pvalue from Voltage to Charge

  fSignalShapePar[0] = pValue;  
  fSignalShapePar[1] = fDelayFromTrigger+time;  //0.482

  Int_t binCharge = (Int_t)((fDelayFromTrigger + time)/fTimeres);       

  if (binCharge >= fNrOfBins) return fNrOfBins;  //condizione di controllo

  Int_t kBin = binCharge;
  Int_t nPeakMax = 0;
  Double_t xtmp, ampVal = 999.;   
  Double_t assetime[fNrOfBins];
  for (int iii=0;iii<fNrOfBins;iii++) assetime[iii] = iii*fTimeres;

  for (; kBin<fNrOfBins; kBin++){  

    xtmp = kBin*fTimeres;
    if (ampVal > signalShape2(&xtmp,fSignalShapePar)) {
      nPeakMax = kBin; 
      ampVal = signalShape2(&xtmp,fSignalShapePar);
    }

    signR[kBin] += signalShape2(&xtmp,fSignalShapePar);
    signL[kBin] += signalShape2(&xtmp,fSignalShapePar);
  }               

  return nPeakMax;

}

//______________________________________________________________________________
void GMCTGenerateWaveforms::AddNoise() {

  //Add gaussian noise: sigma = fNoiseLevel     
  for (Int_t i=0; i<fNrOfBins; i++){
      
    signL[i] += gRandom->Gaus(0.,fNoiseLevel);
    signR[i] += gRandom->Gaus(0.,fNoiseLevel);  

  }

}

//______________________________________________________________________________
void GMCTGenerateWaveforms::FEGain() {

  for (Int_t i=0; i<fNrOfBins; i++){
      
    signL[i] *= gAnalyzer->GetGenerateWaveformsTaskBase()->GetSP()->GetFrontEndGain();
    signR[i] *= gAnalyzer->GetGenerateWaveformsTaskBase()->GetSP()->GetFrontEndGain();

  }

}

//______________________________________________________________________________
void GMCTGenerateWaveforms::WaveDigitization() {

  Double_t  tmpVal;
  Double_t  Baseline = 50.;  //baseline 50 mv;
  for (Int_t i=0; i<fNrOfBins; i++){
    tmpVal = Baseline - signL[i];    //invert signal (to have a positive signal amplitude)
    signL[i] = (int)(0.001*tmpVal/fDigitStep);
    if (signR[i] > fADCbin) signL[i] = fADCbin;  

    tmpVal = Baseline - signR[i];
    signR[i] = (int)(0.001*tmpVal/fDigitStep);  
    if (signR[i] > fADCbin) signR[i] = fADCbin;  
    
  }

  //Riconverto per ottenere mV
  for (Int_t i=0; i<fNrOfBins; i++){
    signL[i] *= -(fDigitStep*1000);
    signL[i] += Baseline;
    signR[i] *= -(fDigitStep*1000);
    signR[i] += Baseline;
  }

}


//______________________________________________________________________________
void GMCTGenerateWaveforms::EndOfRun()
{
}

//______________________________________________________________________________
void GMCTGenerateWaveforms::Terminate()
{
}

//______________________________________________________________________________
Double_t signalShape2(Double_t *x, Double_t *p) {
  Double_t tmpx = x[0]-p[1];
  Double_t sign = 0.0;
  if (tmpx>0.0) {  
    sign = p[4]*TMath::Exp(-tmpx/p[2])/p[2];
    sign+= (1.0-p[4])*TMath::Exp(-tmpx/p[5])/p[5];
    sign*= (1.0-TMath::Exp(-tmpx/p[3]));
    sign*= (p[2]+p[3])*(p[3]+p[5]);
    sign/= (p[2]*p[5]+p[3]*p[5]+p[2]*p[3]*p[4]-p[3]*p[4]*p[5]);
    sign*=p[0];
  }

  //  sign+=p[6];
  return sign;

}



