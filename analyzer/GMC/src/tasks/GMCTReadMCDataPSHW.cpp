////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// GMCTReadMCDataPSHW                                                         //
//                                                                            //
// Begin_Html <!--
/*-->

<!--*/
// --> End_Html
//                                                                            //
//                                                                            //
// Please note: The following information is only correct after executing     //
// the ROMEBuilder.                                                           //
//                                                                            //
// This task accesses the following folders :                                 //
//     GeantTrack                                                             //
//     PSHWGeantStep                                                          //
//     PSHWHit                                                                //
//                                                                            //
//                                                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

/* Generated header file containing necessary includes                        */
#include "generated/GMCTReadMCDataPSHWGeneratedIncludes.h"

////////////////////////////////////////////////////////////////////////////////
/*  This header was generated by ROMEBuilder. Manual changes above the        *
 * following line will be lost next time ROMEBuilder is executed.             */
/////////////////////////////////////----///////////////////////////////////////

#include "generated/GMCAnalyzer.h"
#include "tasks/GMCTReadMCDataPSHW.h"
#include "ROMEiostream.h"

// uncomment if you want to include headers of all folders
//#include "GMCAllFolders.h"

#include "generated/GMCConfigToForm.h"
#include "generated/GMCMCRunHeader.h"
#include "generated/GMCGeantTrack.h"
#include "generated/GMCPSHWGeantStep.h"
#include "generated/GMCPSHWHit.h"

#include "generated/GMCEventLoop.h"

#include "ROGeometryHandle.hh"

ClassImp(GMCTReadMCDataPSHW)

//______________________________________________________________________________
void GMCTReadMCDataPSHW::Init()
{

  //instance the geometry
  fGeometry = GeometrySrvHndl::Instance();
  fGeometry->makeDetectors();

  _minimumEnergy = GetSP()->GetHitEnergyCut();
}

//______________________________________________________________________________
void GMCTReadMCDataPSHW::BeginOfRun()
{
}

//______________________________________________________________________________
void GMCTReadMCDataPSHW::Event()
{

  if (gAnalyzer->GetCurrentEventNumber()%100==0 ) { printf("\n ***** ReadMCDataPSHW Load event %lld ************************ \n",gAnalyzer->GetCurrentEventNumber()); }
  LoadEvent(gAnalyzer->GetCurrentEventNumber());

  gAnalyzer->SetPSHWHitSize(0);
  DefineHitsInPSHW();

}

//______________________________________________________________________________
void GMCTReadMCDataPSHW::LoadEvent(Int_t nev) {

  fBrDataPshw=gAnalyzer->GetPSHWGeantSteps();
  fBrDataTrk=gAnalyzer->GetGeantTracks();
}

//______________________________________________________________________________
void GMCTReadMCDataPSHW::DefineHitsInPSHW() {

  Int_t NrHits = 0;

  PshwChHitMap hitmap;
  hitmap.clear();
  fillHitMap( hitmap );

  for(PshwChHitMap::const_iterator idchn = hitmap.begin(); idchn != hitmap.end(); ++idchn) {
    unsigned long dchn_id = idchn->first;
//    fGeometry->GetCellHandle()->SelectCellDet(dchn_id);
//    CLHEP::Hep3Vector const& mid   = fGeometry->GetCellHandle()->GetCellCenter();
//    CLHEP::Hep3Vector const& w     = fGeometry->GetCellHandle()->GetCellDirection();

    std::vector<GMCPSHWGeantStep *> const& ihits = idchn->second;

    std::map<int,GMCPSHWHit *> hitPerTrk;
    std::map<int,GMCPSHWHit *>::iterator hitPerTrk_it;

//    std::map<int,TVector3 *> hitEndPerTrk;
//    std::map<int,TVector3 *>::iterator hitEndPerTrk_it;

    for( size_t i=0; i<ihits.size(); i++ ) {

      GMCPSHWGeantStep& hit = *(ihits.at(i));

      hitPerTrk_it = hitPerTrk.find(hit.GetfTrackID());
//      hitEndPerTrk_it = hitEndPerTrk.find(hit.GetfTrackID());

      if ( hitPerTrk_it == hitPerTrk.end() ) {
        gAnalyzer->SetPSHWHitSize(NrHits+1);
        GMCPSHWHit *ahit = gAnalyzer->GetPSHWHitAt(NrHits);
        ++NrHits;
        hitPerTrk.insert( std::pair<int,GMCPSHWHit *>( hit.GetfTrackID(), ahit ) );
//        hitEndPerTrk.insert( std::pair<int,TVector3 *>( hit.GetfTrackID(), hit.GetfPosEnding() ) );

        ahit->SetfChanId( dchn_id );
        ahit->SetfTrkIDofHit( hit.GetfTrackID() );
        ahit->SetfEntranceX( hit.GetfPos()->X() );
        ahit->SetfEntranceY( hit.GetfPos()->Y() );
        ahit->SetfEntranceZ( hit.GetfPos()->Z() );
        ahit->SetfEntranceMomX( hit.GetfMomentum()->X() );
        ahit->SetfEntranceMomY( hit.GetfMomentum()->Y() );
        ahit->SetfEntranceMomZ( hit.GetfMomentum()->Z() );
        ahit->SetfTotalEnergyLoss( hit.GetfEdep() );
        ahit->SetfGlobalTime( hit.GetfGlobalTime() );
        ahit->SetfToF( hit.GetfProperTime() );

      } else {
        GMCPSHWHit *ahit = hitPerTrk_it->second;
        ahit->SetfTotalEnergyLoss( ahit->GetfTotalEnergyLoss() + hit.GetfEdep() );
//        ahit->SetfLength( ahit->GetfLength() + hit.GetfStepLen() );
//        hitEndPerTrk_it->second=hit.GetfPos();
      }

    }

    //FIXME!!! Trackes Pile-up not taken into account!!!

//    for ( hitPerTrk_it=hitPerTrk.begin(); hitPerTrk_it!=hitPerTrk.end(); ++hitPerTrk_it ) {
//      GMCPSHWHit *ahit = hitPerTrk_it->second;
//    }

  }

  if (gAnalyzer->GetCurrentEventNumber()%100==0 ) { std::cout<<"PSHW nhits "<<gAnalyzer->GetPSHWHitSize()<<std::endl; }

}

//______________________________________________________________________________
void GMCTReadMCDataPSHW::fillHitMap ( PshwChHitMap& hitmap ) {
  int nentries=fBrDataPshw->GetEntries();
  for (int istep=0;istep<nentries;istep++) {
    GMCPSHWGeantStep *aStep = (GMCPSHWGeantStep *)fBrDataPshw->At(istep);
    if( (aStep->GetfEdep()-aStep->GetfNoIEdep())<_minimumEnergy ) continue; // Skip steps with very low energy deposition
    unsigned long unqChId = fGeometry->GetPSHWROChanHandle()->computeDet(aStep->GetfChamberNr(),aStep->GetfChannelNr());
    hitmap[unqChId].push_back(aStep);
  }
}

//______________________________________________________________________________
void GMCTReadMCDataPSHW::EndOfRun()
{
}

//______________________________________________________________________________
void GMCTReadMCDataPSHW::Terminate()
{
}
