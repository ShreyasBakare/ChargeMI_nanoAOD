#define nano9Ana_cxx
// The class definition in nano9Ana.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("nano9Ana.C")
// root> T->Process("nano9Ana.C","some options")
// root> T->Process("nano9Ana.C+")
//195479682  //78177 //270425
int evtno=270425;//195479682;  //enter the event number you want event display of
bool Mu =true;     //make false if you dont want goodMu event display
bool Ele =true;     //make false if you dont want goodElec event display
bool gen = true;   //make false if you dont want genPart event display

#include "nano9Ana.h"
#include <TH2.h>
#include <TStyle.h>
using namespace std;

int m =0;
int e =0;
int ssMu =0;
int gssMu =0;
int sse =0;
int gsse =0;
int bothe =0;
int bothMu =0;
int errorMu =0;
int errorElec =0;
void nano9Ana::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();
}

void nano9Ana::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();
  
  //Initialization of the counters:
  nEvtRan        = 0;
  nEvtTotal      = 0;
  //Other custom counters can be initialized here.

   _HstFile = new TFile(_HstFileName,"recreate");
  BookHistograms();
}

void nano9Ana::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
  
   _HstFile->Write();
  _HstFile->Close();

  //The following lines are displayed on the root prompt.
  cout<<"Total events ran = "<<nEvtRan<<endl;
  cout<<"Total good events = "<<nEvtTotal<<endl;

  //The following lines are written on the sum_<process name>.txt file
  ofstream fout(_SumFileName);
  fout<<"Total events ran = "<<nEvtRan<<endl;
  fout<<"Total good events  = "<<nEvtTotal<<endl;
}

void nano9Ana::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file. 
}

Bool_t nano9Ana::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.
  
  fReader.SetLocalEntry(entry);
  if(_data == 0)
    fReader_MC.SetLocalEntry(entry);
  if(_data == 1)
    fReader_Data.SetLocalEntry(entry);

  //Verbosity determines the number of processed events after which the root prompt is supposed to display a status update.
  if(_verbosity==0 && nEvtTotal%1000000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;      
  else if(_verbosity>0 && nEvtTotal%1000000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;
  
  //The following flags throws away some events based on unwanted properties (such as detector problems)
  GoodEvt2018 = (_year==2018 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && *Flag_BadPFMuonFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  GoodEvt2017 = (_year==2017 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && *Flag_BadPFMuonFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  GoodEvt2016 = (_year==2016 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && *Flag_BadPFMuonFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  
  GoodEvt = GoodEvt2018 && GoodEvt2017 && GoodEvt2016;
  
  nEvtRan++;                             //Total number of events containing everything (including the trash events).
  
  if(GoodEvt){
    nEvtTotal++;                         //Total number of events containing goodEvents
                                         //The analysis is done for these good events.


    //Construction of the arrays:
    
    //goodMu array :
    int nmu = 0;                         // This counts the number of muons in each event.
    goodMu.clear();                      // Make sure to empty the array from previous event.
    for(unsigned int i=0; i<(*nMuon); i++){
                                         // This loop runs over all the muon candidates. Some of them will pass our selection criteria.
                                         // These will be stored in the goodMu array.
      Lepton temp;                       // 'temp' is the i-th candidate.
      temp.v.SetPtEtaPhiM(Muon_pt[i],Muon_eta[i],Muon_phi[i],0.105); //the muon mass in GeV is 0.105
      temp.id = -13*Muon_charge[i];      //pdgID for mu- = 13, pdgID for mu+ = -13  
      temp.ind = i;

      //These are the flags the 'temp' object i.e. the muon candidate has to pass.
      bool passCuts = temp.v.Pt()>15 && fabs(temp.v.Eta())<2.4 && Muon_mediumId[i];
      passCuts = passCuts && Muon_pfRelIso04_all[i]<0.15;
      passCuts = passCuts && fabs(Muon_dxy[i])<0.05 && fabs(Muon_dz[i])<0.1;
      
      if(passCuts){
	goodMu.push_back(temp);          // If 'temp' satisfies all the conditions, it is pushed back into goodMu
	nmu++;                           // Everytime a 'temp' passes the flags, this counter increases by one.
      }
    }                                    // This 'for' loop has created a goodMu array.
    //Now we sort the goodMu in decreasing pT order
    Sort(0);                     
    
      
    //goodElec array :
    int nElec = 0;                         // This counts the number of Electrons in each event.
    goodElec.clear();                      // Make sure to empty the array from previous event.
    for(unsigned int i=0; i<(*nElectron); i++){
      // This loop runs over all the Electron candidates. Some of them will pass our selection criteria.
      // These will be stored in the goodElec array.
      Lepton temp;                       // 'temp' is the i-th candidate.
      temp.v.SetPtEtaPhiM(Electron_pt[i],Electron_eta[i],Electron_phi[i],0.000511); //the Electron mass in GeV is 0.000511
      temp.id = -11*Electron_charge[i];      //pdgID for e- = 11, pdgID for e+ = -11  
      temp.ind = i;
      
      //These are the flags the 'temp' object i.e. the electron candidate has to pass.
      bool passCuts = temp.v.Pt()>15 && fabs(temp.v.Eta())<2.4 && Electron_mvaFall17V2Iso_WP80[i];
      bool noMu = true;
      //remove temp closer to Mu
      for(int j=0; j<(int)goodMu.size(); j++){   
	if (goodMu.at(j).v.DeltaR(temp.v)<0.4){
	  noMu=false;
	}
      }
      
      passCuts = passCuts && noMu;
      //// passCuts = passCuts && fabs(Muon_dxy[i])<0.05 && fabs(Muon_dz[i])<0.1;
      
      if(passCuts){
	goodElec.push_back(temp);          // If 'temp' satisfies all the conditions, it is pushed back into goodElec
	nElec++;                           // Everytime a 'temp' passes the flags, this counter increases by one.
      }
    }                                    // This 'for' loop has created a goodElec array.
    
    //Now we sort the goodElec in decreasing pT order
    Sort(1);
    
    // goodPh array :
      int nPh = 0;                         // This counts the number of Photons in each event.
      goodPh.clear();                      // Make sure to empty the array from previous event.
      for(unsigned int i=0; i<(*nPhoton); i++){
                                         // This loop runs over all the Photon candidates. Some of them will pass our selection criteria.
                                         // These will be stored in the goodPh array.
	Lepton temp;                       // 'temp' is the i-th candidate.
	temp.v.SetPtEtaPhiM(Photon_pt[i],Photon_eta[i],Photon_phi[i],0); //the Photon mass in GeV is 0
	temp.id = 22;     //pdgID for ph = 22.
	temp.ind = i;

      //These are the flags the 'temp' object i.e. the Photon candidate has to pass.
	bool passCuts =  temp.v.Pt()>15 && fabs(temp.v.Eta())<2.4 && Photon_mvaID_WP80[i];
	////passCuts = passCuts && Muon_pfRelIso04_all[i]<0.15;
	//// passCuts = passCuts && fabs(Muon_dxy[i])<0.05 && fabs(Muon_dz[i])<0.1;

	bool noMu = true;
	bool noElec = true;
	//remove temp closer to Mu
	for(int j=0; j<(int)goodMu.size(); j++){   
	  if (goodMu.at(j).v.DeltaR(temp.v)<0.4){
	    noMu= false;
	  } 
	}
	//remove temp closer to Elec
	for(int j=0; j<(int)goodElec.size(); j++){   
	  if (goodElec.at(j).v.DeltaR(temp.v)<0.4){
	    noElec= false;
	  }
	}
	
	passCuts = passCuts && noElec && noMu;
	if(passCuts){
	  goodPh.push_back(temp);          // If 'temp' satisfies all the conditions, it is pushed back into goodPh
	  nPh++;                           // Everytime a 'temp' passes the flags, this counter increases by one.
	}
      }                                    // This 'for' loop has created a goodPh array.
    
    //Now we sort the goodPh in decreasing pT order
    Sort(2);


    
    // goodJet array :
      int nJets = 0;                         // This counts the number of Jets in each event.
      goodJet.clear();                      // Make sure to empty the array from previous event.
      for(unsigned int i=0; i<(*nJet); i++){
                                         // This loop runs over all the Jet candidates. Some of them will pass our selection criteria.
                                         // These will be stored in the goodJet array.
	Lepton temp;                       // 'temp' is the i-th candidate.
	temp.v.SetPtEtaPhiM(Jet_pt[i],Jet_eta[i],Jet_phi[i],Jet_mass[i]);
	//temp.id = N/A;     //pdgID for Jets = N/A
	temp.ind = i;

      //These are the flags the 'temp' object i.e. the Photon candidate has to pass.
	bool passCuts =  temp.v.Pt()>30 && fabs(temp.v.Eta())<5;
	////passCuts = passCuts && Muon_pfRelIso04_all[i]<0.15;
	//// passCuts = passCuts && fabs(Muon_dxy[i])<0.05 && fabs(Muon_dz[i])<0.1;
	
	bool noMu = true;
	bool noElec = true;
	bool noPh=true;
	//remove temp closer to Mu
	for(int j=0; j<(int)goodMu.size(); j++){   
	  if (goodMu.at(j).v.DeltaR(temp.v)<0.4){
	    noMu=false;
	  }
	}
	//remove temp closer to Elec
	for(int j=0; j<(int)goodElec.size(); j++){   
	  if (goodElec.at(j).v.DeltaR(temp.v)<0.4){
	    noElec=false;
	  }
	}
	//remove temp closer to Ph
	for(int j=0; j<(int)goodPh.size(); j++){   
	  if (goodPh.at(j).v.DeltaR(temp.v)<0.4){
	    noPh=false;
	  }
	}
	
	passCuts = passCuts && noElec && noMu && noPh;
	if(passCuts){
	  goodJet.push_back(temp);          // If 'temp' satisfies all the conditions, it is pushed back into goodJet
	  nJets++;                           // Everytime a 'temp' passes the flags, this counter increases by one.
	}
      }                                    // This 'for' loop has created a goodJet array.
    
      //Now we sort the goodJet in decreasing pT order
      Sort(3);                           
      
      
      // genPart array :
      int ngPart=0;
      int ngMu = 0;                         // This counts the number of genParticles in each event.
      int nge = 0;                         // This counts the number of genParticles in each event.
      int ngPh = 0;                         // This counts the number of genParticles in each event.
      genMu.clear();                      // Make sure to empty the array from previous event.
      genPart.clear();                      // Make sure to empty the array from previous event.
      genElec.clear();                      // Make sure to empty the array from previous event.
      genPh.clear();                      // Make sure to empty the array from previous event.

      if(nEvtTotal==evtno && gen) cout<<" ## goodMu.size()= "<<goodMu.size()<<" goodElec.size()= "<<goodElec.size()<<" goodPh.size()= "<<goodMu.size()<<" goodJet.size()= "<<goodJet.size()<<" genPart.size()= "<<genPart.size()<<endl<<"\nindex"<<"\tpdgID"<<"\tmomIndex "<<"momID"<<"\tstatus"<<"\teta"<<endl;
      
      for(unsigned int i=0; i<(*nGenPart); i++){
	// This loop runs over all the genParticles. Some of them will pass our selection criteria.
	// These will be stored in the genMu, genElec, genPh arrays.
	//########### Event display ###############
	if(nEvtTotal==evtno  && gen){
	  cout<<i<<"\t";
	  cout<<GenPart_pdgId[i]<<"\t";
	  cout<<GenPart_genPartIdxMother[i]<<"\t";
	  cout<<GenPart_pdgId[GenPart_genPartIdxMother[i]]<<"\t";
	  cout<<GenPart_status[i]<< "\t";
	  cout<<fabs(GenPart_eta[i])<<endl;
	}
	//#########################################
	Lepton temp;                       // 'temp' is the i-th Particle.
	temp.v.SetPtEtaPhiM(GenPart_pt[i],GenPart_eta[i],GenPart_phi[i],GenPart_mass[i]);
	temp.ind = i;
	temp.id = GenPart_pdgId[i];
	temp.momid  =GenPart_pdgId[GenPart_genPartIdxMother[i]];
	temp.momind=GenPart_genPartIdxMother[i];

	if (GenPart_status[i]==1){
	  genPart.push_back(temp);        // If 'temp' satisfies all the gPart conditions, it is pushed back into genPart
	  ngPart++;
	}
	
	if (fabs(temp.id)==13 && GenPart_status[i]==1 && fabs(temp.v.Eta())<2.4){
	  genMu.push_back(temp);          // If 'temp' satisfies all the gMu conditions, it is pushed back into genMu
	  ngMu++;
	}
	
	if (fabs(temp.id)==11 && GenPart_status[i]==1 && fabs(temp.v.Eta())<2.4){
	  genElec.push_back(temp);          // If 'temp' satisfies all the gElec conditions, it is pushed back into genElec
	  nge++;
	}
	
	if (fabs(temp.id)==22 && GenPart_status[i]==1 && fabs(temp.v.Eta())<2.4){
	  genPh.push_back(temp);          // If 'temp' satisfies all the gPh conditions, it is pushed back into genPh
	  ngPh++;
	}
      }
      //Now we sort the genPart in decreasing pT order (all in the same sort)
      Sort(4);      

      
      //##############
      // Analysis
      //##############
      // ############################################################################################################################################################     
      //Filling DeltaR for goodMu and closest genPart for DiMu
      if(goodMu.size()>1){
	h.Muprop[0]-> Fill((goodMu.at(0).v+goodMu.at(1).v).M());
	if((goodMu.at(0).v+goodMu.at(1).v).M()<=106 &&(goodMu.at(0).v+goodMu.at(1).v).M()>76){
	  for(int j=0; j<(int)goodMu.size(); j++){                                        //Find DRlgoodMug for all goodMu in same plot
	    float DRlgoodMug = 1000;
	    for(int i=0; i<(int)genPart.size(); i++){
	      if(goodMu.at(j).v.DeltaR(genPart.at(i).v) < DRlgoodMug){
		DRlgoodMug= goodMu.at(j).v.DeltaR(genPart.at(i).v);
		goodMu.at(j).genmatch = genPart.at(i).ind;
		goodMu.at(j).genId = genPart.at(i).id;
	      }
	    }
	    h.goodMugdeltaR[2]->Fill(DRlgoodMug);
	    if(j==0)  h.goodMugdeltaR[3]->Fill(DRlgoodMug);
	    if(j==1)  h.goodMugdeltaR[4]->Fill(DRlgoodMug);
	    h.Muprop[1]->Fill(goodMu.at(j).v.Pt());
	    h.Muprop[3]->Fill(goodMu.at(j).v.Eta());
	    if(goodMu.at(j).id ==13 &&  goodMu.at(j).genId == -13){
	      if(j!=0) errorMu++;
	      h.Muprop[2]->Fill(goodMu.at(j).v.Pt());
	      h.Muprop[4]->Fill(goodMu.at(j).v.Eta());
	    }
	    if(goodMu.at(j).id ==-13 &&  goodMu.at(j).genId == 13){
	      if(j!=0) errorMu++;
	      h.Muprop[2]->Fill(goodMu.at(j).v.Pt());
	      h.Muprop[4]->Fill(goodMu.at(j).v.Eta());
	    }
	  }
	  if((goodMu.at(0).id * goodMu.at(0).genId) < 0 && (goodMu.at(1).id * goodMu.at(1).genId) < 0){
	    bothMu++;
	  }
	  if(goodMu.at(0).id==goodMu.at(1).id){
	    ssMu++;
	    h.Muprop[12]->Fill(goodMu.at(0).v.Pt());
	    h.Muprop[14]->Fill(goodMu.at(0).v.Eta());
	  }
	  if(goodMu.at(0).genId==goodMu.at(1).genId){
	    // cout <<"$$$$$$$"<<nEvtTotal<<endl;
	    gssMu++;
	  }
	  
	  h.Zmu[0]->Fill(goodMu.at(0).genId);
	  h.Zmu[1]->Fill(goodMu.at(1).genId);
	}
	if((goodMu.at(0).v+goodMu.at(1).v).M()<=76 &&(goodMu.at(0).v+goodMu.at(1).v).M()>61){
	  if(goodMu.at(0).id==goodMu.at(1).id){
	    h.Muprop[13]->Fill(goodMu.at(0).v.Pt());
	    h.Muprop[15]->Fill(goodMu.at(0).v.Eta());
	  }
	}
	if((goodMu.at(0).v+goodMu.at(1).v).M()<=121 &&(goodMu.at(0).v+goodMu.at(1).v).M()>106){
	  if(goodMu.at(0).id==goodMu.at(1).id){
	    h.Muprop[13]->Fill(goodMu.at(0).v.Pt());
	    h.Muprop[15]->Fill(goodMu.at(0).v.Eta());
	  }
	}
      }
      
      if(nEvtTotal==evtno && Mu) cout<<"\ni-gMu"<<"\tsize"<<"\tindex"<<"\tpdgID"<<"\tgenID"<<endl;
      
      for(unsigned int i=0; i<(goodMu.size()); i++){
	//########### Event display ###############
	if(nEvtTotal==evtno && Mu){
	  cout<<i<<"\t";
	  cout<<goodMu.size()<<"\t";
	  cout<<goodMu.at(i).ind<<"\t";
	  cout<<goodMu.at(i).id<<"\t";
	  cout<<goodMu.at(i).genId<<endl;
	}
      }
    
      //Filling DeltaR for goodElec and closest genPart DiElec
      if(goodElec.size()>1){
	h.eprop[0]-> Fill((goodElec.at(0).v+goodElec.at(1).v).M());
	if((goodElec.at(0).v+goodElec.at(1).v).M()<106 &&(goodElec.at(0).v+goodElec.at(1).v).M()>76){
	  for(int j=0; j<(int)goodElec.size(); j++){                                        //Find DRlgoodeg for all goodElec in same plot
	    float DRlgoodeg = 1000;
	    for(int i=0; i<(int)genPart.size(); i++){
	      if(goodElec.at(j).v.DeltaR(genPart.at(i).v) < DRlgoodeg){
		DRlgoodeg= goodElec.at(j).v.DeltaR(genPart.at(i).v);
		goodElec.at(j).genmatch = genPart.at(i).ind;
		goodElec.at(j).genId = genPart.at(i).id;
	      }
	    }
	    h.goodegdeltaR[2]->Fill(DRlgoodeg);
	    if(j==0)  h.goodegdeltaR[3]->Fill(DRlgoodeg);
	    if(j==1)  h.goodegdeltaR[4]->Fill(DRlgoodeg);
	  }
	  for(int j=0; j<2; j++){  
	    h.eprop[1]->Fill(goodElec.at(j).v.Pt());
	    h.eprop[3]->Fill(goodElec.at(j).v.Eta());
	    if(goodElec.at(j).id == 11 &&  goodElec.at(j).genId == -11){
	      if(j!=0 && fabs(goodElec.at(0).v.Pt()-goodElec.at(1).v.Pt())>15){
		errorElec++;
	      }
	      //cout <<"#######-"<<nEvtTotal<<endl;
	      h.eprop[2]->Fill(goodElec.at(j).v.Pt());
	      h.eprop[4]->Fill(goodElec.at(j).v.Eta());
	    }
	    if(goodElec.at(j).id == -11 &&  goodElec.at(j).genId == 11){
	      if(j!=0 && fabs(goodElec.at(0).v.Pt()-goodElec.at(1).v.Pt())>15){
		errorElec++;
	      }
	      //cout <<"#######+"<<nEvtTotal<<endl;
	      h.eprop[2]->Fill(goodElec.at(j).v.Pt());
	      h.eprop[4]->Fill(goodElec.at(j).v.Eta());
	    }
	  }
	  if((goodElec.at(0).id * goodElec.at(0).genId) < 0 && (goodElec.at(1).id * goodElec.at(1).genId) < 0){
	    bothe++;
	  }
	  if(goodElec.at(0).id==goodElec.at(1).id){
	    // cout <<"$$$$$$$"<<nEvtTotal<<endl;
	    sse++;
	    h.eprop[12]->Fill(goodElec.at(0).v.Pt());
	    h.eprop[14]->Fill(goodElec.at(0).v.Eta());
	  }
	  if(goodElec.at(0).genId==goodElec.at(1).genId){
	    gsse++;
	    // cout <<"######"<<nEvtTotal<<endl;
	  }
	  
	  h.Ze[0]->Fill(goodElec.at(0).genId);
	  h.Ze[1]->Fill(goodElec.at(1).genId);
	}
	if((goodElec.at(0).v+goodElec.at(1).v).M()<=76 &&(goodElec.at(0).v+goodElec.at(1).v).M()>61){
	  if(goodElec.at(0).id==goodElec.at(1).id){
	    h.eprop[13]->Fill(goodElec.at(0).v.Pt());
	    h.eprop[15]->Fill(goodElec.at(0).v.Eta());
	  }
	}
	if((goodElec.at(0).v+goodElec.at(1).v).M()<=121 &&(goodElec.at(0).v+goodElec.at(1).v).M()>106){
	  if(goodElec.at(0).id==goodElec.at(1).id){
	    h.eprop[13]->Fill(goodElec.at(0).v.Pt());
	    h.eprop[15]->Fill(goodElec.at(0).v.Eta());
	  }
	}
      }
      if(nEvtTotal==evtno && Ele) cout<<"\ni-ge"<<"\tsize"<<"\tindex"<<"\tpdgID"<<"\tgenID"<<"\teeta"<<"\tept"<<endl;
      
      for(unsigned int i=0; i<(goodElec.size()); i++){
	//########### Event display ###############
	if(nEvtTotal==evtno && Ele){
	  cout<<i<<"\t";
	  cout<<goodElec.size()<<"\t";
	  cout<<goodElec.at(i).ind<<"\t";
	  cout<<goodElec.at(i).id<<"\t";
	  cout<<goodElec.at(i).genId<<"\t";
	  cout<<goodElec.at(i).v.Eta()<<"\t";
	  cout<<goodElec.at(i).v.Pt()<<endl;
	}
      }
      // ############################################################################################################################################################     
      if(goodMu.size()>2) m++;
      if(goodElec.size()>2) e++;
      if(nEvtTotal==evtno) cout<<"\ngM3"<<"\t"<<m<<"\tge3"<<"\t"<<e<<"\terrorMu"<<"\t"<<errorMu<<"\terrorElec"<<"\t"<<errorElec<<"\tssMu"<<"\t"<<ssMu<<"\tgssMu"<<"\t"<<gssMu<<"\tsse"<<"\t"<<sse<<"\tgsse"<<"\t"<<gsse<<"\tbothMu"<<"\t"<<bothMu<<"\tbothe"<<"\t"<<bothe<<endl;
      //########### ANALYSIS ENDS HERE ##############
  }//GoodEvt
  
  return kTRUE;
}




//######################################
//        USER DEFINED FUNCTIONS
//######################################
void nano9Ana::Sort(int opt)
{
  //This functions sorts an array in the decreasing order of pT.
  //For goodMu:
  if(opt==0){
    for(int i=0; i<(int)goodMu.size()-1; i++){
      for(int j=i+1; j<(int)goodMu.size(); j++){
	if( goodMu[i].v.Pt() < goodMu[j].v.Pt() ) swap(goodMu.at(i),goodMu.at(j));
      }
    }
  }
  
  //For goodElec:
  if(opt==1){
    for(int i=0; i<(int)goodElec.size()-1; i++){
      for(int j=i+1; j<(int)goodElec.size(); j++){
	if( goodElec[i].v.Pt() < goodElec[j].v.Pt() ) swap(goodElec.at(i),goodElec.at(j));
      }
    }
  }
  
  //For goodPh:
  if(opt==2){
    for(int i=0; i<(int)goodPh.size()-1; i++){
      for(int j=i+1; j<(int)goodPh.size(); j++){
	if( goodPh[i].v.Pt() < goodPh[j].v.Pt() ) swap(goodPh.at(i),goodPh.at(j));
      }
    }
  }
  
  //For goodJet:
  if(opt==3){
    for(int i=0; i<(int)goodJet.size()-1; i++){
      for(int j=i+1; j<(int)goodJet.size(); j++){
	if( goodJet[i].v.Pt() < goodJet[j].v.Pt() ) swap(goodJet.at(i),goodJet.at(j));
      }
    }
  }
    
  //For genPart:
  if(opt==4){
    for(int i=0; i<(int)genPart.size()-1; i++){
      for(int j=i+1; j<(int)genPart.size(); j++){
	if( genPart[i].v.Pt() < genPart[j].v.Pt() ) swap(genPart.at(i),genPart.at(j));
      }
    }
    for(int i=0; i<(int)genMu.size()-1; i++){
      for(int j=i+1; j<(int)genMu.size(); j++){
	if( genMu[i].v.Pt() < genMu[j].v.Pt() ) swap(genMu.at(i),genMu.at(j));
      }
    }
    for(int i=0; i<(int)genElec.size()-1; i++){
      for(int j=i+1; j<(int)genElec.size(); j++){
	if( genElec[i].v.Pt() < genElec[j].v.Pt() ) swap(genElec.at(i),genElec.at(j));
      }
    }
    for(int i=0; i<(int)genPh.size()-1; i++){
      for(int j=i+1; j<(int)genPh.size(); j++){
	if( genPh[i].v.Pt() < genPh[j].v.Pt() ) swap(genPh.at(i),genPh.at(j));
      }
    }
  }
  
  
  //Repeat this for the other arrays here.
}
int nano9Ana::GenMother(int ind, int mom_ind)
{
  int p_id = GenPart_pdgId[ind];
  int m_id = GenPart_pdgId[mom_ind];
  while(p_id==m_id){
    ind = mom_ind;
    mom_ind = GenPart_genPartIdxMother[ind];
    p_id = GenPart_pdgId[ind];
    m_id = GenPart_pdgId[mom_ind];
  }
  return m_id;
}

float nano9Ana::delta_phi(float phi1, float phi2)
{
  //The correct deltaPhi falls in the interval [0 , pi]
  phi1 = TVector2::Phi_0_2pi(phi1);
  phi2 = TVector2::Phi_0_2pi(phi2);
  float dphi = fabs(phi1 - phi2);
  if(dphi>TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
  return dphi;
}

float nano9Ana::transv_mass(float pt_lep, float MET, float dphi)
{
  //The inputs are the pt of the lepton, MET and dPhi between the lepton and MET
  float mT = sqrt(2* pt_lep * MET *(1-cos(dphi)));
  return mT;
}

void nano9Ana::BookHistograms()
{
  //The histograms are booked here.
  //Binning etc are done here.
  //These histograms are stored in the hst_<process name>.root file in the same order.

  //Example : new TH1F ("hst_name", "hst title", NBins, startVal, EndVal);
  
  h.Muprop[0] = new TH1F("muM","inv mass of Mu0 & Mu1",20000,0,2000);
  h.Muprop[1] = new TH1F("den_mupt","pt of all muons",200,0,200);
  h.Muprop[2] = new TH1F("num_mupt","pt of muons with charge MI",200,0,200);
  h.Muprop[3] = new TH1F("den_mueta","eta of all muons",400,-4,4);
  h.Muprop[4] = new TH1F("num_mueta","eta of muons with charge MI",400,-4,4);
  h.Muprop[12] = new TH1F("num_muptDATA","pt of muons with charge MI w/o genPart",200,0,200);
  h.Muprop[14] = new TH1F("num_muetaDATA","eta of muons with charge MI w/o genPart",400,-4,4);
  h.Muprop[13] = new TH1F("num_muptMinus","pt of muons with charge MI w/o genPart in 61-76 && 106-121",200,0,200);
  h.Muprop[15] = new TH1F("num_muetaMinus","eta of muons with charge MI w/o genPart in 61-76 && 106-121",400,-4,4);
  
  h.eprop[0] = new TH1F("eM","inv mass of e0 & e1",20000,0,2000);
  h.eprop[1] = new TH1F("den_ept","pt of all ELecs",200,0,200);
  h.eprop[2] = new TH1F("num_ept","pt of Elecs with charge MI",200,0,200);
  h.eprop[3] = new TH1F("den_eeta","eta of all Elecs",400,-4,4);
  h.eprop[4] = new TH1F("num_eeta","eta of Elecs with charge MI",400,-4,4);
  h.eprop[12] = new TH1F("num_eptDATA","pt of Elecs with charge MI w/o genPart",200,0,200);
  h.eprop[14] = new TH1F("num_eetaDATA","eta of Elecs with charge MI w/o genPart",400,-4,4);
  h.eprop[13] = new TH1F("num_eptMinus","pt of Elecs with charge MI w/o genPart in 61-76 && 106-12",200,0,200);
  h.eprop[15] = new TH1F("num_eetaMinus","eta of Elecs with charge MI w/o genPart in 61-76 && 106-12",400,-4,4);
  
  // h.ptreso[0] = new TH1F("ptMuReso","(goodMu Pt - genMu Pt)/goodMu Pt",600,-3,3);
  // h.ptreso[1] = new TH1F("ptElecReso","(goodElec Pt - genElec Pt)/goodElec Pt",600,-3,3);

  h.goodMugdeltaR[2] = new TH1F("2goodMu_gPart_deltaR","Solid angle differences of all goodMu and matching gPart",500,0,10);
  h.goodMugdeltaR[3] = new TH1F("goodMu0_gPart_deltaR","Delta R of Mu0 and matching gPart",500,0,10);
  h.goodMugdeltaR[4] = new TH1F("goodMu1_gPart_deltaR","Delta R of Mu1  and matching gPart",500,0,10);
  
  h.Zmu[0] =  new TH1F("Zmu0ID", "ID of closest genPart to Mu0",200,-100,100);
  h.Zmu[1] =  new TH1F("Zmu1ID", "ID of closest genPart to Mu1",200,-100,100);

  h.goodegdeltaR[2] = new TH1F("2goode_gPart_deltaR","Solid angle differences of all goodElec and matching gPart",500,0,10);
  h.goodegdeltaR[3] = new TH1F("goode0_gPart_deltaR","Delta R of e0 and matching gPart",500,0,10);
  h.goodegdeltaR[4] = new TH1F("goode1_gPart_deltaR","Delta R of e1  and matching gPart",500,0,10);
  
  h.Ze[0] =  new TH1F("Ze0ID", "ID of closest genPart to e0",200,-100,100);
  h.Ze[1] =  new TH1F("Ze1ID", "ID of closest genPart to e1",200,-100,100);
}
 
