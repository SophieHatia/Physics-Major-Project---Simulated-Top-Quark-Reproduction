#include "TQPOffline/TQPAnalysis.h"
#include "TTHbbConfiguration/GlobalConfiguration.h"
#include <math.h>
#include <iomanip>
// #include "RooUnfold/RooUnfold.h"

namespace TTHbb{


  // Here is where you put your global variables:
  int N_totalJets = 0;
  float SmallestJetPt = 100.0;

  //___________________________________________
  //
  TQPAnalysis::TQPAnalysis(std::string name) : ToolBase(true), m_outFile(0), m_hist_WlepPt(0), m_hist_WlepPx(0)
                                               {
    m_name = name;
    auto* config = TTHbb::GlobalConfiguration::get();


    
    //m_DeltaRCut_STDT=1.0;
    
    m_btaggingWP = (*config)("BTaggingWP");
    
    std::string muIDstring = (*config)("TQPAnalysis.SMTmuID"); // tight or medium+MI
    m_SMTmuID = -1;
    if(muIDstring=="medium+MI") m_SMTmuID = 0;
    if(muIDstring=="tight")     m_SMTmuID = 1;

    std::string btagSELstring = (*config)("TQPAnalysis.BTagSel"); // btagging selection (none, double-tag, one b outside, one b anywhere)
    m_BTagSel = 0; // none
    if(btagSELstring=="doubleTag")      m_BTagSel = 1; // double-tag
    if(btagSELstring=="oneTagOutside")  m_BTagSel = 2; // one b outside SMT jet
    if(btagSELstring=="oneTagAnywhere") m_BTagSel = 3; // one b anywhere
    
    m_tqph = new TQPHelper();
    m_tqph->setTaggingType(m_SMTmuID);//m_SMTmuID is the type of selection, 0 is momImbalance, 1 is Tight (default)
    
    // corrections
    m_corrections = new TQPCorrections();
    std::vector<int> dsIdList; dsIdList.clear();
    dsIdList.push_back(410525); // PH7
    dsIdList.push_back(410225); // aMCP8
    dsIdList.push_back(410501); // PP8
    dsIdList.push_back(410511); // PP8 rad Hi
    dsIdList.push_back(410512); // PP8 rad Lo
    
    // fragmentation function
//     m_corrections->InitFragFuncRew("../TQPOffline/data/MCweights_v6.root",dsIdList);
//     m_corrections->InitFragFuncRew("../TQPOffline/data/MCweights_v7.root",dsIdList);
//     m_corrections->InitFragFuncRew("../TQPOffline/data/MCweights_v7bis.root",dsIdList);
//     m_corrections->InitFragFuncRew("../TQPOffline/data/MCweights_v8.root",dsIdList);
//     m_corrections->InitFragFuncRew("../TQPOffline/data/MCweights_v8_cheat.root",dsIdList);
    //m_corrections->InitFragFuncRew("../TQPOffline/data/MCweights_v9b.root",dsIdList);
    m_corrections->InitFragFuncRew("../TQPOffline/data/MCweights_v12.root",dsIdList); //changed on 19/2/2018

    // ttbar mass
    m_corrections->InitTtMassRew("../TQPOffline/data/ttKinRew.root",dsIdList);    

    // B prod fractions and decay
//     m_corrections->InitProdFracRew("../TQPOffline/data/MCweights_v7.root",dsIdList);
//     m_corrections->InitBRsRew("../TQPOffline/data/MCweights_v7.root",dsIdList);
//     m_corrections->InitProdFracRew("../TQPOffline/data/MCweights_v7bis.root",dsIdList);
//     m_corrections->InitBRsRew("../TQPOffline/data/MCweights_v7bis.root",dsIdList);
    m_corrections->InitProdFracRew("../TQPOffline/data/MCweights_v12.root",dsIdList);  //changed on 19/2/2018
    m_corrections->InitBRsRew("../TQPOffline/data/MCweights_v12.root",dsIdList);  //changed on 19/2/2018
//     m_corrections->InitProdFracRew("../TQPOffline/data/MCweights_v4.root",dsIdList);
//     m_corrections->InitBRsRew("../TQPOffline/data/MCweights_v4.root",dsIdList);
    
    // Output file here and book histograms here
    m_outFile = new TFile("OUTPUTFILE.root","RECREATE"); 
    m_hist_WlepPt = new TH1F( "WlepPt_hist", "W lepton pT (GeV)", 100, 20., 200.);
    m_hist_WlepPx = new TH1F( "WlepPx_hist", "W lepton px (GeV)", 100, 20., 200.);


    m_outFile->cd();
   // ...
    
    m_selections=std::vector<TString>();
    m_selections.push_back("");
    m_selections.push_back("doubleTag");
    m_selections.push_back("oneTagOutside");
    m_selections.push_back("oneTagAnywhere");
    m_selections.push_back("doubleTag_pTmu5");
    m_selections.push_back("oneTagOutside_pTmu5");
    m_selections.push_back("oneTagAnywhere_pTmu5");
    m_selections.push_back("doubleTag_pTmu6");
    m_selections.push_back("oneTagOutside_pTmu6");
    m_selections.push_back("oneTagAnywhere_pTmu6");
    m_selections.push_back("doubleTag_pTjet40");
    m_selections.push_back("oneTagOutside_pTjet40");
    m_selections.push_back("oneTagAnywhere_pTjet40");
    m_selections.push_back("doubleTag_pTmujetRatioGr0p1");
    m_selections.push_back("oneTagOutside_pTmujetRatioGr0p1");
    m_selections.push_back("oneTagAnywhere_pTmujetRatioGr0p1");
    
    
    for(int i=0; i<(int)m_selections.size(); i++)
    {
		std::map<TString,long int> eventCounters;
        eventCounters["AllEvents"]=0;
        eventCounters["SMTtagged"]=0;
        eventCounters["SMTtagged_topbmu"]=0;
        eventCounters["SMTtagged_topbcmu"]=0;
        eventCounters["SMTtagged_topbtaumu"]=0;
        eventCounters["SMTtagged_topbctaumu"]=0;
        eventCounters["SMTtagged_notopbmu"]=0;
        eventCounters["SMTtagged_notopcmu"]=0;
        eventCounters["SMTtagged_notoptaumu"]=0;
        eventCounters["SMTtagged_notopWmu"]=0;
        eventCounters["SMTtagged_notop_unknown"]=0;
        eventCounters["SMTtagged_fake"]=0;   
        
        m_eventCounters[m_selections[i]]=eventCounters; 
    }
    
    
  }

  //___________________________________________
  //
  void TQPAnalysis::apply(TTHbb::Event* event){
	  
	//FIX ME:
    // 1) Remember to check carefully Mlmu_truth, for now this is filled starting from the truth muon linked from the reco, but we may want to add a truth-level selection
    
    m_tqph->setEvent(event);
    
    // get truth info for soft muons
    unsigned int Ntruthsmtmu = 0;
    m_truthsmtmu= m_tqph->getTruthSoftMuons(Ntruthsmtmu); //Ntruthsmtmu passed by reference and filled
    
    // get all the reco muons
    unsigned int Nsmtmu = 0;
    m_smtmu = m_tqph->getRecoSoftMuons(Nsmtmu);//event->m_customObj["smtmu"]; 
    //std::cout<< "ah "<<event->m_customObj["smtmu"].size()<<" "<<event->m_jets.size()<<std::endl;
    
   for(unsigned int i_jet=0;i_jet<event->m_jets.size();i_jet++)
   {
	    event->m_jets[i_jet]->charVariable("jet_isSMTagged") = false;
	    event->m_jets[i_jet]->charVariable("jet_isSMAntiTagged") = false;
   }

    int njetSMTtagged = 0;
    int highestPtSMTindex=-1;
    double highestPtSMT=-999.;
    std::shared_ptr<Particle> SMT_max_pT = std::make_shared<Particle>();

    // loop on Reco soft muons
    for(unsigned int i_smtmu=0;i_smtmu<Nsmtmu;i_smtmu++){
      //

      std::shared_ptr<Particle> SMT = m_smtmu[i_smtmu];
       
      if(!m_tqph->passSMTPreselection(SMT)) continue;   
      //float smtmu_momImbalance = m_tqph->getSMTMomentumImbalance(SMT);
      
      // find the associated jet
      float dRmin=999.;
      int jetIdx=m_tqph->getJetAssociatedWithSMT(SMT,dRmin); //dRmin is filled by the method

      // make sure to put here all the new decroations on smtmu, so that they are always filled before the continue is called
      SMT->intVariable("smtmu_jetIdx")    = jetIdx;
      SMT->floatVariable("smtmu_dRmin")   = dRmin;
      SMT->charVariable("smtmu_isTagged") = false;
      
      if(jetIdx<0) continue; // we skip the muon if this is not tagging anything
      if(!m_tqph->passSMTSelection(SMT))
      {
		   if(!event->m_jets[jetIdx]->charVariable("jet_isSMTagged")) //so we avoid to anti-tag something we already tagged 
		        event->m_jets[jetIdx]->charVariable("jet_isSMAntiTagged") = true;
		   continue;
	  }
      njetSMTtagged++;
      event->m_jets[jetIdx]->charVariable("jet_isSMTagged") = true;
      event->m_jets[jetIdx]->charVariable("jet_isSMAntiTagged") = false; //this is needed in case one jet is anti-tagged by one soft muon but tagged by another      
      
      if(SMT->p4().Pt()>highestPtSMT)
      {
		highestPtSMTindex=i_smtmu;
		highestPtSMT=SMT->p4().Pt();	
		SMT_max_pT=SMT;
	  }
    } //end of loop on soft muons

    // fill the histogram of the lepton pT:
    TLorentzVector Lep = event->m_leptons[0]->p4();
    float lep_pt     = Lep.Pt();
    //  std::cout<<"Lep pt: "<<lep_pt/1000<<std::endl; 
    m_hist_WlepPt->Fill(lep_pt/1000); //filling the W lepton

    float lep_px     = Lep.Px();
    m_hist_WlepPx->Fill(lep_px/1000);

    // looping over all jets
    std::cout<<"Njets: "<< event->m_jets.size() << std::endl; 
    N_totalJets += event->m_jets.size();
    for(unsigned int i_jet=0;i_jet<event->m_jets.size();i_jet++){
      TLorentzVector Jet = event->m_jets[i_jet]->p4();
      float Jet_pt = Jet.Pt();
      std::cout<<"looping over jets: "<< i_jet << " pt:" << Jet_pt/1000 << std::endl; 
      // check for the smallest jet:
      if (Jet_pt/1000 < SmallestJetPt){
	SmallestJetPt = Jet_pt/1000;
      }
    }
    

    
    this->saveEvent(event,SMT_max_pT);
    
    return;
  }//end of apply
  
  void TQPAnalysis::saveEvent(TTHbb::Event* event, std::shared_ptr<Particle> SMT_max_pT)
  {
	  
    bool isData = event->m_info->isData;
      bool isNominal = false;
      if(event->m_info->currentSyst=="nominal") isNominal=true;
      int selectedSMTjetIndex=-1;
      if(SMT_max_pT->E()>1.) selectedSMTjetIndex=SMT_max_pT->intVariable("smtmu_jetIdx");
      
      TLorentzVector Lep = event->m_leptons[0]->p4();
      TLorentzVector Met(0,0,0,0);
      Met.SetPtEtaPhiE(event->met_met,0,event->met_phi,event->met_met);
      float mtw = sqrt(2. * Lep.Pt() * Met.Pt() * (1. - cos( Lep.DeltaPhi(Met)) ) );
      float lep_pt     = Lep.Pt();
      float lep_eta    = Lep.Eta();
      float lep_phi    = Lep.Phi();
      float dPhi_lnu   = Lep.DeltaPhi(Met);
      int lep_charge   = event->m_leptons[0]->charge();
      int lep_truth_charge = 0;
      int ttbar_channel = 0; // 1: l+jets, 2: dilep   // FIXME: not clear how there can be all-had events in non-all-had sample... But this seems to happen...
     
      TLorentzVector truthLepPlus=m_tqph->getTruthLeptonTLV(1);
      TLorentzVector truthLepMinus=m_tqph->getTruthLeptonTLV(-1);
      TLorentzVector truthLepPlusBeforeRad=m_tqph->getTruthLeptonTLVBeforeRad(1);
      TLorentzVector truthLepMinusBeforeRad=m_tqph->getTruthLeptonTLVBeforeRad(-1);
      TLorentzVector truthLepPlusDressed=m_tqph->getTruthLeptonTLVDressed(1);
      TLorentzVector truthLepMinusDressed=m_tqph->getTruthLeptonTLVDressed(-1);
    
      if(truthLepPlus.E()>0){
        ttbar_channel++;
        if(Lep.DeltaR(truthLepPlus)<0.1)
        {
			  lep_truth_charge =  1;
		}
      }
      if(truthLepMinus.E()>0){
        ttbar_channel++;
        if(Lep.DeltaR(truthLepMinus)<0.1)
        {
			 lep_truth_charge = -1;
		}
      }
      int dsId = event->m_info->mcChannelNumber;
	  bool isNominalSample = false;
	  bool isSherpa=false;
	  if(dsId==410250 ||dsId==410251 || dsId==410252 ) isSherpa=true;
	  if(dsId==410501) isNominalSample = true;  //we store all systematics only for this
	  // set dsId to nominal in some cases
	  if(dsId==410493) dsId = 410501; // pp8 170
	  if(dsId==410494) dsId = 410501; // pp8 175
	  if(dsId==410248) dsId = 410501; // pp8 CF
	  if(dsId==410511) dsId = 410501; // pp8 radLow
	  if(dsId==410512) dsId = 410501; // pp8 radHi
	  if(dsId==410225) dsId = 410501; // aMC@NLO+pp8 
	  // and set it to PP6 in some other cases
	  if(dsId==410037) dsId = 410000; // pp6 170
	  if(dsId==410038) dsId = 410000; // pp6 171.5
	  if(dsId==410039) dsId = 410000; // pp6 173.5
	  if(dsId==410040) dsId = 410000; // pp6 175
	  if(dsId==410041) dsId = 410000; // pp6 177.5
	  
	  ///////////////////INITIALIZING EVERYITHING///////////////////////////////////////
	  float weight_smtmuSF=1.;
	  float weight_smtmuSF_SMT_SF_ID_STAT_UP=1.;
	  float weight_smtmuSF_SMT_SF_ID_STAT_DOWN=1.;
	  float weight_smtmuSF_SMT_SF_ID_SYST_UP=1.;
	  float weight_smtmuSF_SMT_SF_ID_SYST_DOWN=1.;
	  float weight_smtmuSF_SMT_SF_ID_STAT_LOWPT_UP=1.;
	  float weight_smtmuSF_SMT_SF_ID_STAT_LOWPT_DOWN=1.;
	  float weight_smtmuSF_SMT_SF_ID_SYST_LOWPT_UP=1.;
	  float weight_smtmuSF_SMT_SF_ID_SYST_LOWPT_DOWN=1.;
	  int nSMT=0;
	  int nSMA=0;
	  int nSMTandBTag=0;
	  int nSMTnoBTag=0;
	  int nSMAandBTag=0;
	  int nSMAnoBTag=0;
	  int nBTagOutsideSMT=0;
	  
	  float SMTmu_pt=-999.;
	  float SMTmu_eta=-999.;
	  float SMTmu_phi=-999.;
	  float SMTmu_dR=-999.;
	  int SMTmu_charge=-999;
	  int SMTmu_origin=-999;
	  int SMTmu_originFlag=0;
	  float SMTmu_d0sig=-999.;
	  float SMTmu_pTrel=-999.;
	  int SMTmu_smtIdx=-999;
	  
	  int isSTDTTruth=0;
	  float deltaR_lep_SMTmu=-999.;
	  
	  float Mlmu=-999.;
	  float Mlmu_truth=-999.;
	  float Mlmu_truth_beforeRad=-999.;
	  float Mlmu_truth_dressed=-999.;
	  
	  float top_pT=-999.;
	  float antitop_pT=-999.;
	  float ttbar_pT=-999.;
	  float ttbar_mass=-999.;
	  float SMTmu_xB=-999.;
	  
	  //int SMTmu_hadronOrig=-999;
	  int lastBmother_pdgid=0;
	  int lastCmother_pdgid=0;
	  int mother_pdgid=0;
	  
	  bool hasProblemInTTbarInfo=false;
	  float weight_ttbarMass = 1.;
	  float weight_fragFunc= 1.;
	  float weight_fragFunc_up= 1.;
	  float weight_fragFunc_down= 1.;
	  float weight_prodFrac          = 1.;
      float weight_prodFrac_Bup      = 1.;
      float weight_prodFrac_Bsup     = 1.;
      float weight_prodFrac_Baryup   = 1.;
      float weight_prodFrac_Bdown    = 1.;
      float weight_prodFrac_Bsdown   = 1.;
      float weight_prodFrac_Barydown = 1.;
      float weight_BRs               = 1.;
      float weight_BRs_bup           = 1.;
      float weight_BRs_bdown         = 1.;
      float weight_BRs_btocup        = 1.;
      float weight_BRs_btocdown      = 1.;
      float weight_BRs_cup           = 1.;
      float weight_BRs_cdown         = 1.;
      
      float SMTmu_B_pT= -999.;
      float SMTmu_b_pT= -999.;
      float SMTmu_b_pT_afterRad= -999.;
      float SMTmu_B_E=-999.;
      float SMTmu_b_E=-999.;
      float SMTmu_b_E_afterRad=-999.;
    
	  ////////////////////////////FILLING EVERYTHING/////////////////////////////////////
	  
	  for(unsigned int i_jet=0;i_jet<event->m_jets.size();i_jet++)
	  {
		 if(selectedSMTjetIndex>=0 &&i_jet!=selectedSMTjetIndex && event->m_jets[i_jet]->charVariable("isbtagged_"+m_btaggingWP) )
		 {
			 nBTagOutsideSMT++;
		 }
		 if(event->m_jets[i_jet]->charVariable("jet_isSMTagged"))
		 {
			  nSMT++; 
			  if(event->m_jets[i_jet]->charVariable("isbtagged_"+m_btaggingWP)) nSMTandBTag++;
			  else nSMTnoBTag++;
		 }
		 if(event->m_jets[i_jet]->charVariable("jet_isSMAntiTagged"))
		 {
			 nSMA++;
			 if(event->m_jets[i_jet]->charVariable("isbtagged_"+m_btaggingWP)) nSMAandBTag++;
			 else nSMAnoBTag++;
		 }
	  }
	  
	  if(selectedSMTjetIndex>=0) //it means we have a SMT
	  {
		 TLorentzVector SMTjet_tlv=event->m_jets[selectedSMTjetIndex]->p4();
		 TLorentzVector SMTplusLep=Lep+SMT_max_pT->p4();
		 Mlmu=SMTplusLep.M();
		 deltaR_lep_SMTmu=SMTjet_tlv.DeltaR(SMT_max_pT->p4());
		 SMTmu_pt=SMT_max_pT->p4().Pt();
		 SMTmu_eta=SMT_max_pT->p4().Eta();
		 SMTmu_phi=SMT_max_pT->p4().Phi();
		 SMTmu_dR=SMT_max_pT->p4().DeltaR(SMTjet_tlv);
		 SMTmu_charge=SMT_max_pT->intVariable("smtmu_q");
		 if(!isData) 
		 SMTmu_d0sig=SMT_max_pT->floatVariable("smtmu_d0_sig");
		 SMTmu_pTrel=m_tqph->getPtRel(SMT_max_pT,SMTjet_tlv);
		 SMTmu_smtIdx=SMT_max_pT->intVariable("smtmu_idx"); //this is now filled in TQPHelper::getRecoSoftMuons
		 if(!isData)
		  {
		   //get the corresponding truth muon (if available) for later use
		   int i_corr_truth=-1;
		   if(SMT_max_pT->checkIntVariable("smtmu_truthmuIndex")) i_corr_truth=SMT_max_pT->intVariable("smtmu_truthmuIndex");

		   std::shared_ptr<Particle> truthSMTmu = std::make_shared<Particle>();
		   if(i_corr_truth>=0)
		   {
			   truthSMTmu=m_truthsmtmu[i_corr_truth];
		       if(!isSherpa) isSTDTTruth=m_tqph->getSTDTTruth(truthSMTmu,lep_truth_charge);
		       SMTmu_origin=m_tqph->getSMTOriginCode(truthSMTmu);
		       //if(isSherpa) SMTmu_origin=truthSMTmu->intVariable("truthmu_origin");
		       
		       SMTmu_originFlag=truthSMTmu->intVariable("truthmu_originFlag");
		       
		       if(lep_truth_charge>0)
		       {
				   	Mlmu_truth=(truthSMTmu->p4()+truthLepPlus).M(); 
				   	Mlmu_truth_beforeRad=(m_tqph->getTruthMuonTLVBeforeRad(truthSMTmu)+truthLepPlusBeforeRad).M(); 
				   	Mlmu_truth_dressed=(m_tqph->getTruthMuonTLVDressed(truthSMTmu)+truthLepPlusDressed).M(); 
			   }
		       if(lep_truth_charge<0)
		       {
				   	Mlmu_truth=(truthSMTmu->p4()+truthLepMinus).M(); 
				   	Mlmu_truth_beforeRad=(m_tqph->getTruthMuonTLVBeforeRad(truthSMTmu)+truthLepMinusBeforeRad).M(); 
				   	Mlmu_truth_dressed=(m_tqph->getTruthMuonTLVDressed(truthSMTmu)+truthLepMinusDressed).M(); 
			   }
		       if(m_tqph->isTruthTTbarInfoAvailable() && !m_tqph->hasProblemInTTbarTruthInfo())
		       {
				   SMTmu_xB=m_tqph->getXB(truthSMTmu);
				   if(m_tqph->isMuonFromTopOrAntiTopDecayChain(truthSMTmu)) //excluded cases are for a muon not from top decay chain
				   {
				    TLorentzVector b_tlv, bltv_afterRad;
				    if(m_tqph->isMuonFromTopDecayChain(truthSMTmu))
				    {
						 b_tlv=m_tqph->getTruthBottomTLV();
						 bltv_afterRad=m_tqph->getTruthBottomAfterRadTLV();
					}
				    if(m_tqph->isMuonFromAntiTopDecayChain(truthSMTmu))
				    {
						  b_tlv=m_tqph->getTruthAntiBottomTLV();
						  bltv_afterRad=m_tqph->getTruthAntiBottomAfterRadTLV();
					}
				    SMTmu_b_pT=b_tlv.Pt();
                    SMTmu_b_pT_afterRad=bltv_afterRad.Pt();
                    SMTmu_b_E=b_tlv.E();
                    SMTmu_b_E_afterRad=bltv_afterRad.E();
				   }
			   }
			   
			   if(truthSMTmu->floatVariable("truthmu_lastBmother_pt")>0.) //the muon is coming in some way from a B
			   {
				   lastBmother_pdgid=truthSMTmu->intVariable("truthmu_lastBmother_pdgid");
			   }
			   
			   if(truthSMTmu->floatVariable("truthmu_lastCmother_pt")>0.) //the muon is coming in some way from a B
			   {
				   lastCmother_pdgid=truthSMTmu->intVariable("truthmu_firstCmother_pdgid"); //FIX ME: this must be lastCmother, missing for now in the ntuples but we will soon have it! (this should be a very small effect though)
			   }
			   
			   mother_pdgid=truthSMTmu->intVariable("truthmu_mother_pdgid");
			   
			   if(truthSMTmu->floatVariable("truthmu_lastBmother_pt")>0.) 
			   {
				   SMTmu_B_pT=truthSMTmu->floatVariable("truthmu_lastBmother_pt");
				   SMTmu_B_E=truthSMTmu->floatVariable("truthmu_lastBmother_E");
			   }

			   if(!isSherpa)
			   {
				   if(lastBmother_pdgid!=0){  //we have frag frac rwt only for muons coming from a B-hadron
					  weight_prodFrac            = m_corrections->GetProdFracRew(lastBmother_pdgid,dsId,"nominal");
					  if(isNominalSample){
						weight_prodFrac_Bup      = m_corrections->GetProdFracRew(lastBmother_pdgid,dsId,"Bup");
						weight_prodFrac_Bsup     = m_corrections->GetProdFracRew(lastBmother_pdgid,dsId,"Bsup");
						weight_prodFrac_Baryup   = m_corrections->GetProdFracRew(lastBmother_pdgid,dsId,"Baryup");
						weight_prodFrac_Bdown    = m_corrections->GetProdFracRew(lastBmother_pdgid,dsId,"Bdown");
						weight_prodFrac_Bsdown   = m_corrections->GetProdFracRew(lastBmother_pdgid,dsId,"Bsdown");
						weight_prodFrac_Barydown = m_corrections->GetProdFracRew(lastBmother_pdgid,dsId,"Barydown");
					  }
					}
					weight_BRs           = m_corrections->GetBRsRew(mother_pdgid,lastBmother_pdgid,lastCmother_pdgid,dsId,"nominal");
					if(isNominalSample){
						weight_BRs_bup       = m_corrections->GetBRsRew(mother_pdgid,lastBmother_pdgid,lastCmother_pdgid,dsId,"bup");
						weight_BRs_bdown     = m_corrections->GetBRsRew(mother_pdgid,lastBmother_pdgid,lastCmother_pdgid,dsId,"bdown");
						weight_BRs_btocup    = m_corrections->GetBRsRew(mother_pdgid,lastBmother_pdgid,lastCmother_pdgid,dsId,"btocup");
						weight_BRs_btocdown  = m_corrections->GetBRsRew(mother_pdgid,lastBmother_pdgid,lastCmother_pdgid,dsId,"btocdown");
						weight_BRs_cup       = m_corrections->GetBRsRew(mother_pdgid,lastBmother_pdgid,lastCmother_pdgid,dsId,"cup");
						weight_BRs_cdown     = m_corrections->GetBRsRew(mother_pdgid,lastBmother_pdgid,lastCmother_pdgid,dsId,"cdown");
					}
			  }
 
		       
		   }
		   if(SMT_max_pT->checkFloatVariable("smtmu_recoeff_tight_SF")) weight_smtmuSF=SMT_max_pT->floatVariable("smtmu_recoeff_tight_SF");
		   if(SMT_max_pT->checkFloatVariable("smtmu_recoeff_tight_syst_SF_ID_STAT_UP")) weight_smtmuSF_SMT_SF_ID_STAT_UP=SMT_max_pT->floatVariable("smtmu_recoeff_tight_syst_SF_ID_STAT_UP");
		   if(SMT_max_pT->checkFloatVariable("smtmu_recoeff_tight_syst_SF_ID_SYST_UP")) weight_smtmuSF_SMT_SF_ID_SYST_UP=SMT_max_pT->floatVariable("smtmu_recoeff_tight_syst_SF_ID_SYST_UP");
		   if(SMT_max_pT->checkFloatVariable("smtmu_recoeff_tight_syst_SF_ID_STAT_LOWPT_UP")) weight_smtmuSF_SMT_SF_ID_STAT_LOWPT_UP=SMT_max_pT->floatVariable("smtmu_recoeff_tight_syst_SF_ID_STAT_LOWPT_UP");
		   if(SMT_max_pT->checkFloatVariable("smtmu_recoeff_tight_syst_SF_ID_SYST_LOWPT_UP")) weight_smtmuSF_SMT_SF_ID_SYST_LOWPT_UP=SMT_max_pT->floatVariable("smtmu_recoeff_tight_syst_SF_ID_SYST_LOWPT_UP");
			   
		   if(SMT_max_pT->checkFloatVariable("smtmu_recoeff_tight_syst_SF_ID_STAT_DOWN")) weight_smtmuSF_SMT_SF_ID_STAT_DOWN=SMT_max_pT->floatVariable("smtmu_recoeff_tight_syst_SF_ID_STAT_DOWN");
		   if(SMT_max_pT->checkFloatVariable("smtmu_recoeff_tight_syst_SF_ID_SYST_DOWN")) weight_smtmuSF_SMT_SF_ID_SYST_DOWN=SMT_max_pT->floatVariable("smtmu_recoeff_tight_syst_SF_ID_SYST_DOWN");
		   if(SMT_max_pT->checkFloatVariable("smtmu_recoeff_tight_syst_SF_ID_STAT_LOWPT_DOWN")) weight_smtmuSF_SMT_SF_ID_STAT_LOWPT_DOWN=SMT_max_pT->floatVariable("smtmu_recoeff_tight_syst_SF_ID_STAT_LOWPT_DOWN");
		   if(SMT_max_pT->checkFloatVariable("smtmu_recoeff_tight_syst_SF_ID_SYST_LOWPT_DOWN")) weight_smtmuSF_SMT_SF_ID_SYST_LOWPT_DOWN=SMT_max_pT->floatVariable("smtmu_recoeff_tight_syst_SF_ID_SYST_LOWPT_DOWN");
		   		   
		  }//end of if(!isData)
      }//end of if(selectedSMTjetIndex>=0)
      
      //now we fill the ttbar truth information if available
      if(!isData)
      {
		   if(m_tqph->isTruthTTbarInfoAvailable() && !isSherpa)
		   {
			  hasProblemInTTbarInfo=m_tqph->hasProblemInTTbarTruthInfo(); //this may happen, for instance, in 4t events; more or less 0.05% of the events
			  TLorentzVector top_tlv=m_tqph->getTruthTopTLV();
			  TLorentzVector antitop_tlv=m_tqph->getTruthAntiTopTLV();
			  top_pT=top_tlv.Pt();
			  antitop_pT=antitop_tlv.Pt();
			  TLorentzVector ttbar_tlv=top_tlv+antitop_tlv;
			  ttbar_pT=ttbar_tlv.Pt();
			  ttbar_mass=ttbar_tlv.M();
			  
			  //now the reweightings
			  weight_ttbarMass = m_corrections->GetTtMassWeight(ttbar_mass/1e3,dsId);
			  
			  if(SMTmu_xB>0 && SMTmu_xB<100.)
			  {
				   weight_fragFunc = m_corrections->GetFragFuncWeight(SMTmu_xB,dsId);
			       if(isNominalSample){
						weight_fragFunc_up   = m_corrections->GetFragFuncWeight(SMTmu_xB,dsId,"up");
						weight_fragFunc_down = m_corrections->GetFragFuncWeight(SMTmu_xB,dsId,"down");
					}
			   }

		   }//end of if(m_tqph->isTruthTTbarInfoAvailable())
	  }//end of if(!isData)
	  ////////////////////////SAVING EVERYTHING//////////////////////////////////////////////////
	  event->intVariable("nSMT")                            = nSMT;  //these quantities are now jet based, instead of muon based
      event->intVariable("nSMA")                            = nSMA; //these quantities are now jet based, instead of muon based
      event->intVariable("nSMTandBTag"+m_btaggingWP)        = nSMTandBTag; //these quantities are now jet based, instead of muon based
      event->intVariable("nSMTnoBTag"+m_btaggingWP)         = nSMTnoBTag; //these quantities are now jet based, instead of muon based
      event->intVariable("nSMAandBTag"+m_btaggingWP)        = nSMAandBTag; //these quantities are now jet based, instead of muon based
      event->intVariable("nSMAnoBTag"+m_btaggingWP)         = nSMAnoBTag; //these quantities are now jet based, instead of muon based
      event->intVariable("nBTag"+m_btaggingWP+"OutsideSMT") = nBTagOutsideSMT; //these quantities are now jet based, instead of muon based
      event->intVariable("njetSMTagged")       = nSMT; //kept for backward compatibility for the time being
      
      if(!isData)
      {
       event->charVariable("hasProblemInTTbarInfo") = hasProblemInTTbarInfo;      
       event->intVariable("ttbar_channel")      = ttbar_channel;
       event->floatVariable("top_pT")           = top_pT;
       event->floatVariable("antitop_pT")           = antitop_pT;
       event->floatVariable("ttbar_pT")         = ttbar_pT;
       event->floatVariable("ttbar_mass")       = ttbar_mass;
       event->floatVariable("SMTmu_xB")         = SMTmu_xB;
		   //std::cout<<"debug: isTTbarInfoAvailable="<<m_tqph->isTruthTTbarInfoAvailable()<<" hasProblemInTTbarInfo="<<hasProblemInTTbarInfo<<" ttbar_channel="<<ttbar_channel<<" top_pT="<<top_pT<<" antitop_pT="<<antitop_pT<<" ttbar_pT="<<ttbar_pT<<" ttbar_mass="<<ttbar_mass<<" xB="<<SMTmu_xB<<std::endl;
		   
		event->floatVariable("SMTmu_B_pT")         = SMTmu_B_pT;
		event->floatVariable("SMTmu_b_pT")         = SMTmu_b_pT;
		event->floatVariable("SMTmu_b_pT_afterRad")= SMTmu_b_pT_afterRad;
		event->floatVariable("SMTmu_B_E")          = SMTmu_B_E;
		event->floatVariable("SMTmu_b_E")          = SMTmu_b_E;
		event->floatVariable("SMTmu_b_E_afterRad") = SMTmu_b_E_afterRad;
		
		//std::cout<<"truth tlv "<<event->floatVariable("SMTmu_B_pT")<<" "<<event->floatVariable("SMTmu_b_pT")<<" "<<event->floatVariable("SMTmu_B_E")<<" "<<event->floatVariable("SMTmu_b_E")<<" "<<event->floatVariable("SMTmu_b_pT_afterRad")<<" "<<event->floatVariable("SMTmu_b_E_afterRad")<<std::endl;

       event->floatVariable("weight_ttbarMass")     = weight_ttbarMass;
       event->floatVariable("weight_fragFunc")      = weight_fragFunc;
       if(isNominalSample){
		   event->floatVariable("weight_fragFunc_up")   = weight_fragFunc_up;
		   event->floatVariable("weight_fragFunc_down") = weight_fragFunc_down;
	   }
       event->floatVariable("weight_prodFrac")          = weight_prodFrac;
       if(isNominalSample){
		   event->floatVariable("weight_prodFrac_Bup")      = weight_prodFrac_Bup;
		   event->floatVariable("weight_prodFrac_Bsup")     = weight_prodFrac_Bsup;
		   event->floatVariable("weight_prodFrac_Baryup")   = weight_prodFrac_Baryup;
		   event->floatVariable("weight_prodFrac_Bdown")    = weight_prodFrac_Bdown;
		   event->floatVariable("weight_prodFrac_Bsdown")   = weight_prodFrac_Bsdown;
		   event->floatVariable("weight_prodFrac_Barydown") = weight_prodFrac_Barydown;
	   }
       event->floatVariable("weight_BRs")          = weight_BRs;
       if(isNominalSample){
		   event->floatVariable("weight_BRs_bup")      = weight_BRs_bup;
		   event->floatVariable("weight_BRs_bdown")    = weight_BRs_bdown;
		   event->floatVariable("weight_BRs_btocup")   = weight_BRs_btocup;
		   event->floatVariable("weight_BRs_btocdown") = weight_BRs_btocdown;
		   event->floatVariable("weight_BRs_cup")      = weight_BRs_cup;
		   event->floatVariable("weight_BRs_cdown")    = weight_BRs_cdown;
       }
      }//end of if(!isData)
      
      if(!isData) event->charVariable("isST_truth")         = (isSTDTTruth==1);
      if(!isData) event->charVariable("isDT_truth")         = (isSTDTTruth==(-1));
      if(!isData) event->charVariable("no_ST_DT_truth") = (isSTDTTruth==0);
      event->floatVariable("deltaR_lep_SMTmu") = deltaR_lep_SMTmu;

      event->floatVariable("mtw")              = mtw;
      event->floatVariable("lep_pt")           = lep_pt;
      event->floatVariable("lep_eta")          = lep_eta;
      event->floatVariable("lep_phi")          = lep_phi;
      event->intVariable("lep_charge")         = lep_charge;
      if(!isData) event->intVariable("lep_truth_charge")   = lep_truth_charge;
      event->floatVariable("dPhi_lnu")         = fabs(dPhi_lnu);
      event->floatVariable("Mlmu")             = Mlmu;
      if(!isData) event->floatVariable("Mlmu_truth")             = Mlmu_truth;
      if(!isData) event->floatVariable("Mlmu_truth_beforeRad")             = Mlmu_truth_beforeRad;
      if(!isData) event->floatVariable("Mlmu_truth_dressed")             = Mlmu_truth_dressed;
      
      //std::cout<<"reco "<<Mlmu<<" truth "<<Mlmu_truth<<" truth_beforeRad "<<Mlmu_truth_beforeRad<<" truth_dressed "<<Mlmu_truth_dressed<<std::endl;
     
      event->floatVariable("SMTmu_pt")         = SMTmu_pt;
      event->floatVariable("SMTmu_eta")        = SMTmu_eta;
      event->floatVariable("SMTmu_phi")        = SMTmu_phi;
      event->floatVariable("SMTmu_dR")         = SMT_max_pT->floatVariable("smtmu_dRmin");
      event->intVariable("SMTmu_charge")       = SMTmu_charge;
      if(!isData) event->intVariable("SMTmu_origin")       = SMTmu_origin;
      if(!isData) event->intVariable("SMTmu_originFlag")       = SMTmu_originFlag;
      if(!isData) event->intVariable("SMTmu_mother_pdgid")       = mother_pdgid;
      if(!isData) event->intVariable("SMTmu_Bmother_pdgid")       =lastBmother_pdgid;
      if(!isData) event->intVariable("SMTmu_Cmother_pdgid")       =lastCmother_pdgid;
      event->floatVariable("SMTmu_d0sig")      = SMTmu_d0sig;
      event->floatVariable("SMTmu_pTrel")      = SMTmu_pTrel;
      event->intVariable("SMTmu_smtIdx")      = SMTmu_smtIdx;
      event->intVariable("SMTmu_jetIdx")       = selectedSMTjetIndex;
      
      //std::cout<<" -----"<<std::endl;
      //std::cout<<" smtmu pt "<<event->floatVariable("SMTmu_pt")<<" origin "<<event->intVariable("SMTmu_origin")<<" originFlag "<<event->intVariable("SMTmu_originFlag")<<" mother "<<event->intVariable("SMTmu_mother_pdgid")<<" Bmother "<<event->intVariable("SMTmu_Bmother_pdgid")<<" Cmother "<<event->intVariable("SMTmu_Cmother_pdgid")<<std::endl;
      //std::cout<<" weight ttbar mass "<<event->floatVariable("weight_ttbarMass")<<std::endl;
      //std::cout<<" weight frag func "<<event->floatVariable("weight_fragFunc")<<" up "<<event->floatVariable("weight_fragFunc_up")<<" down "<<event->floatVariable("weight_fragFunc_down")<<std::endl;
      //std::cout<<" weight prod frac "<<event->floatVariable("weight_prodFrac")<<" syst "<<event->floatVariable("weight_prodFrac_Bup")<<" "<<event->floatVariable("weight_prodFrac_Bsup")<<" "<<event->floatVariable("weight_prodFrac_Baryup")<<" "<<event->floatVariable("weight_prodFrac_Bdown")<<" "<<event->floatVariable("weight_prodFrac_Bsdown")<<" "<<event->floatVariable("weight_prodFrac_Barydown")<<std::endl;
      //std::cout<<" BR "<<event->floatVariable("weight_BRs")<<" syst "<<event->floatVariable("weight_BRs_bup")<<" "<<event->floatVariable("weight_BRs_bdown")<<" "<<event->floatVariable("weight_BRs_btocup")<<" "<<event->floatVariable("weight_BRs_btocdown")<<" "<<event->floatVariable("weight_BRs_cup")<<" "<<event->floatVariable("weight_BRs_cdown")<<std::endl;

	  if(!isData){
		  event->floatVariable("weight_smtmuSF")                           = weight_smtmuSF;
		  if(isNominal) //to save space
		  {
			  event->floatVariable("weight_smtmuSF_SMT_SF_ID_STAT_UP")         = weight_smtmuSF_SMT_SF_ID_STAT_UP;
			  event->floatVariable("weight_smtmuSF_SMT_SF_ID_STAT_DOWN")       = weight_smtmuSF_SMT_SF_ID_STAT_DOWN;
			  event->floatVariable("weight_smtmuSF_SMT_SF_ID_SYST_UP")         = weight_smtmuSF_SMT_SF_ID_SYST_UP;
			  event->floatVariable("weight_smtmuSF_SMT_SF_ID_SYST_DOWN")       = weight_smtmuSF_SMT_SF_ID_SYST_DOWN;
			  event->floatVariable("weight_smtmuSF_SMT_SF_ID_STAT_LOWPT_UP")   = weight_smtmuSF_SMT_SF_ID_STAT_LOWPT_UP;
			  event->floatVariable("weight_smtmuSF_SMT_SF_ID_STAT_LOWPT_DOWN") = weight_smtmuSF_SMT_SF_ID_STAT_LOWPT_DOWN;
			  event->floatVariable("weight_smtmuSF_SMT_SF_ID_SYST_LOWPT_UP")   = weight_smtmuSF_SMT_SF_ID_SYST_LOWPT_UP;
			  event->floatVariable("weight_smtmuSF_SMT_SF_ID_SYST_LOWPT_DOWN") = weight_smtmuSF_SMT_SF_ID_SYST_LOWPT_DOWN;
			  
			  //std::cout<<"weights SF "<<event->floatVariable("weight_smtmuSF")<<" "<<event->floatVariable("weight_smtmuSF_SMT_SF_ID_STAT_UP")<<" "<<event->floatVariable("weight_smtmuSF_SMT_SF_ID_STAT_DOWN")<<" "<<event->floatVariable("weight_smtmuSF_SMT_SF_ID_SYST_UP")<<" "<<event->floatVariable("weight_smtmuSF_SMT_SF_ID_SYST_DOWN")<<" "<<event->floatVariable("weight_smtmuSF_SMT_SF_ID_STAT_LOWPT_UP")<<" "<<event->floatVariable("weight_smtmuSF_SMT_SF_ID_STAT_LOWPT_DOWN")<<" "<<event->floatVariable("weight_smtmuSF_SMT_SF_ID_SYST_LOWPT_UP")<<" "<<event->floatVariable("weight_smtmuSF_SMT_SF_ID_SYST_LOWPT_DOWN")<<std::endl;
		  }
	  }//end of if(!isData)
	  
	  
	//we save some event counters for quick cross-checks
	for(int i=0; i<(int)m_selections.size(); i++)
	{

		 m_eventCounters[m_selections[i]]["AllEvents"]++;
		 if(SMTmu_pt<100.) continue; //it means the event has no SMT tagged jets
		 if(m_selections[i].Contains("doubleTag") && nSMTandBTag<1) continue;
		 if(m_selections[i].Contains("oneTagOutside") && nBTagOutsideSMT<1) continue;
		 if(m_selections[i].Contains("oneTagAnywhere") && (nBTagOutsideSMT<1 && nSMTandBTag<1)) continue;
		 if(m_selections[i].Contains("pTmu5") && SMTmu_pt<5000.) continue;
		 if(m_selections[i].Contains("pTmu6") && SMTmu_pt<6000.) continue;
		 if(m_selections[i].Contains("pTjet40") && event->m_jets[selectedSMTjetIndex]->p4().Pt()<40000.) continue;
		 if(m_selections[i].Contains("pTmujetRatioGr0p1") && (SMTmu_pt/event->m_jets[selectedSMTjetIndex]->p4().Pt())<0.1) continue;
		 	 
		 m_eventCounters[m_selections[i]]["SMTtagged"]++;
		 if(abs(SMTmu_originFlag)==0) m_eventCounters[m_selections[i]]["SMTtagged_fake"]++;
		 else if(abs(SMTmu_originFlag)==1) m_eventCounters[m_selections[i]]["SMTtagged_topbmu"]++;
		 else if(abs(SMTmu_originFlag)==10) m_eventCounters[m_selections[i]]["SMTtagged_topbcmu"]++;
		 else if(abs(SMTmu_originFlag)==100) m_eventCounters[m_selections[i]]["SMTtagged_topbtaumu"]++;
		 else if(abs(SMTmu_originFlag)==110) m_eventCounters[m_selections[i]]["SMTtagged_topbctaumu"]++;
		 else if(abs(SMTmu_originFlag)==1001) m_eventCounters[m_selections[i]]["SMTtagged_notopbmu"]++;
		 else if(abs(SMTmu_originFlag)==1002||abs(SMTmu_originFlag)==1003) m_eventCounters[m_selections[i]]["SMTtagged_notopcmu"]++;
		 else if(abs(SMTmu_originFlag)==1004||abs(SMTmu_originFlag)==10041||abs(SMTmu_originFlag)==10042||abs(SMTmu_originFlag)==10043) m_eventCounters[m_selections[i]]["SMTtagged_notoptaumu"]++;
		 else if(abs(SMTmu_originFlag)==2000 ||abs(SMTmu_originFlag)==3000 ) m_eventCounters[m_selections[i]]["SMTtagged_notopWmu"]++;
		 else m_eventCounters[m_selections[i]]["SMTtagged_notop_unknown"]++;
		 
	}

	  
	  return;	  
  }//end of saveEvent

  //___________________________________________
  // 
  void TQPAnalysis::finalise(){
	  
	  
    //	return; //not to print these for the time being
	
	std::vector<TString> outord;
	outord.push_back("topbmu");
	outord.push_back("topbcmu");
	outord.push_back("topbtaumu");
	outord.push_back("topbctaumu");
	outord.push_back("notopbmu");
	outord.push_back("notopcmu");
	outord.push_back("notoptaumu");
	outord.push_back("notopWmu");
	outord.push_back("notop_unknown");
	outord.push_back("fake");
	
	std::cout<<std::setprecision(2);
	std::cout<<std::fixed;
	
	std::cout<<"--------Printing event counts--------"<<std::endl;
	double allEv=(double) m_eventCounters[""]["AllEvents"];
	std::cout<<"AllEvents = "<<allEv<<std::endl<<std::endl;
	std:: cout<<"Total Number of jets: " << N_totalJets <<std::endl<<std::endl;
	std:: cout<<"Smallest pT out of all the jets: " << SmallestJetPt <<std::endl<<std::endl;

	m_outFile->cd();
	m_outFile->Write();
	
	return;
  }//end of finalise
    
}
