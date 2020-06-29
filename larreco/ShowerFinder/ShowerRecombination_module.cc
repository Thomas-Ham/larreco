#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//Larsoft includes 
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Shower.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
//#include "sbndcode/RecoUtils/RecoUtils.h"

//Root Includes
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
//C++ Includes 
#include <vector>
#include <iostream>

namespace ana {
  class ShowerRecombination;
}

class ana::ShowerRecombination : public art::EDAnalyzer {
public:
  
  ShowerRecombination(const fhicl::ParameterSet& pset);
  
  void analyze(const art::Event& evt);
  void endJob();
  void beginJob();

private:
  //fcl Parameters
 
  std::string fGenieGenModuleLabel;
  std::string fLArGeantModuleLabel;
  std::string fHitsModuleLabel;
  std::string fTrackModuleLabel;
  
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;

  TH1D* energyDepHist;
  TH1D* energyColHist;
  TH1D* energyRecoHist;
  TH1D* energyRecoColHist;
  TH1D* energyRecoDepHist;

  TH1D* recoRatioCol;
  TH1D* recoRatioDep;
  TH1D* recombinationAll;
  TH1D* recombinationReco;

  TH1D* recoEffDep;
  TH1D* recoEffCol;

  TH1D* recombinationHist;

  TH2D* chargeEnergyHist;
  TH2D* energyChargeHist;
  TH2D* recombination2Hist;
  TH2D* recombination2HistLow;
  TH2D* recombination2HistLowest;
  TH2D* reconstruction2Hist;
  TH2D* primary2Hist;
  TH2D* secondary2Hist;
  TH2D* other2Hist;
  TH2D* channel2Hist;
  TH2D* channel2HistAverage;

  TH2D* recombinationIDEHist;
  TH2D* recombinationIDEHistPrimary;
  TH2D* recombinationIDEHistSecondary;
  TH2D* recombinationIDEHistOther;

  TH2D* end2Hist;
  TH2D* endIDEHist;
  TH2D* endIDEHistPrimary;
  TH2D* endIDEHistSecondary;
  TH2D* endIDEHistOther;

  TH1D* recombinationHistStart;
  TH2D* recombination2HistStart;
  TH2D* recombinationIDEHistStart;

  TH1D* recombinationHistStartPrimary;
  TH2D* recombination2HistStartPrimary;
  TH2D* recombinationIDEHistStartPrimary;
 
  TH1D* recombinationHistStartSecondary;
  TH2D* recombination2HistStartSecondary;
  TH2D* recombinationIDEHistStartSecondary;

  std::ofstream mytestfile;
  std::ofstream myfile1;
  std::ofstream myfile2;
  //std::ofstream primaries;
  //std::ofstream secondaries;
  //std::ofstream others; 

  int evtcounter;
};

ana::ShowerRecombination::ShowerRecombination(const fhicl::ParameterSet& pset) : EDAnalyzer(pset){
 
  fGenieGenModuleLabel   = pset.get<std::string>("GenieGenModuleLabel");
  fLArGeantModuleLabel   = pset.get<std::string>("LArGeantModuleLabel");
  fHitsModuleLabel       = pset.get<std::string>("HitsModuleLabel");
  fTrackModuleLabel      = pset.get<std::string>("TrackModuleLabel");
 }


void ana::ShowerRecombination::beginJob(){
  
  evtcounter=0;
  
  mytestfile.open("MyTestFile.txt");
  myfile1.open ("MyLow.txt");
  myfile2.open ("MyLower.txt");
  //primaries.open ("MyPrimaries.txt");
  //secondaries.open ("MySecondaries.txt");
  //others.open ("MyOthers.txt");

  art::ServiceHandle<art::TFileService> tfs;
 
  std::string fGenieGenModuleLabel;
  std::string fLArGeantModuleLabel;
  std::string fHitsModuleLabel;
  std::string fTrackModuleLabel;

  energyDepHist = tfs->make<TH1D>("energyDepHist","Histogram of Total Energy Deposited",800,0,1.6e3);
  energyColHist = tfs->make<TH1D>("energyColHist","Histogram of Total Energy Collected",800,0,1.6e3);
  energyRecoHist = tfs->make<TH1D>("energyRecoHist","Histogram of Total Energy Reconstructed",800,0,1.6e3);
  energyRecoColHist  =  tfs->make<TH1D>("energyRecoColHist","Histogram of Collected Energy in Reconstructed Hits",800,0,1.6e3);
  energyRecoDepHist = tfs->make<TH1D>("energyRecoDepHist","Histogram of Deposited Energy in Reconstrucetd Hits",800,0,1.6e3);

  recombinationHist =  tfs->make<TH1D>("recombinationHist","Histogram of Recombination",2000,0,1);

  recoRatioCol = tfs->make<TH1D>("recoRatioCol","Ratio of Collected Energy in Reconstructed Hits / All IDEs",100,0,1);  
  recoRatioDep = tfs->make<TH1D>("recoRatioDep","Ratio of Deposited Energy in Reconstructed Hits / All IDEs",100,0,1);
  recombinationAll = tfs->make<TH1D>("recombinationAll","Recombination Factor: All IDEs",75,0,0.75);
  recombinationReco = tfs->make<TH1D>("recombinationReco","Recombination Factor: Reconstructed IDEs",75,0,0.75);

  recoEffDep = tfs->make<TH1D>("recoEffDep","Reconstruction Efficiency: Energy Deposited",60,0.45,0.75);
  recoEffCol = tfs->make<TH1D>("recoEffCol","Reconstruction Efficiency: Energy Collected",60,0.8,1.1);

  chargeEnergyHist = tfs->make<TH2D>("chargeEnergyHist","Reconstructed Energy vs Deposited Energy",500,0,5000,500,0,5000);
  energyChargeHist = tfs->make<TH2D>("energyChargeHist","Deposited Energy vs Reconstructed Charge",1000,0,2e8,500,0,5000);
  recombination2Hist = tfs->make<TH2D>("recombination2Hist","Energy at Wire vs Deposited Energy",501,-0.1,50,501,-0.1,50);
  recombination2HistLow = tfs->make<TH2D>("recombination2HistLow","Energy at Wire vs Deposited Energy",501,-0.1,50,501,-0.1,50);
  recombination2HistLowest = tfs->make<TH2D>("recombination2HistLowest","Energy at Wire vs Deposited Energy",501,-0.1,50,501,-0.1,50);
  reconstruction2Hist = tfs->make<TH2D>("reconstruction2Hist","Reconstructed Energy vs Deposited Energy",501,-0.1,50,501,-0.1,50);
  primary2Hist= tfs->make<TH2D>("primary2Hist","Reconstructed Energy vs Deposited Energy: TrackID=1",501,-0.1,50,500,-0.1,50);
  secondary2Hist= tfs->make<TH2D>("secondary2Hist","Reconstructed Energy vs Deposited Energy: TrackID=-1",501,-0.1,50,500,-0.1,50);
  other2Hist= tfs->make<TH2D>("other2Hist","Reconstructed Energy vs Deposited Energy: TrackID!=+-1",501,-0.1,50,500,-0.1,50);
  channel2Hist = tfs->make<TH2D>("channel2Hist","Channel vs Recombination Factor",11214,0,11213,100,0,1);
  channel2HistAverage = tfs->make<TH2D>("channel2HistAverage","Channel vs Recombination Factor",11214,0,11213,100,0,1);
  recombinationIDEHist = tfs->make<TH2D>("recombinationIDEHist","Recombination vs Deposited Energy",1501,-0.01,15,100,0,1);
  recombinationIDEHistPrimary = tfs->make<TH2D>("recombinationIDEHistPrimary","Recombination vs Deposited Energy: Primary",1501,-0.01,15,100,0,1);
  recombinationIDEHistSecondary = tfs->make<TH2D>("recombinationIDEHistSecondary","Recombination vs Deposited Energy: Secondary",1501,-0.01,15,100,0,1);
  recombinationIDEHistOther = tfs->make<TH2D>("recombinationIDEHistOther","Recombination vs Deposited Energy: Other",1501,-0.01,15,100,0,1);

  end2Hist = tfs->make<TH2D>("end2Hist","Energy at Wire vs Deposited Energy: End of Track",501,-0.1,50,501,-0.1,50);
  endIDEHist =  tfs->make<TH2D>("endIDEHist","Recombination vs Deposited Energy: End",1501,-0.01,15,100,0,1);
  endIDEHistPrimary =  tfs->make<TH2D>("endIDEHistPrimary","Recombination vs Deposited Energy: End: Primary",1501,-0.01,15,100,0,1);
  endIDEHistSecondary =  tfs->make<TH2D>("endIDEHistSecondary","Recombination vs Deposited Energy: End: Secondary",1501,-0.01,15,100,0,1);
  endIDEHistOther =  tfs->make<TH2D>("endIDEHistOther","Recombination vs Deposited Energy: End: Other",1501,-0.01,15,100,0,1);

  recombinationHistStart =  tfs->make<TH1D>("recombinationHistStart","Histogram of Recombination: Start",2000,0,1);
  recombination2HistStart = tfs->make<TH2D>("recombination2HistStart","Energy at Wire vs Deposited Energy: Start",501,-0.1,50,501,-0.1,50);
  recombinationIDEHistStart = tfs->make<TH2D>("recombinationIDEHistStart","Recombination vs Deposited Energy: Start",1501,-0.01,15,100,0,1);

  recombinationHistStartPrimary =  tfs->make<TH1D>("recombinationHistStartPrimary","Histogram of Recombination: Start: Primary",2000,0,1);
  recombination2HistStartPrimary = tfs->make<TH2D>("recombination2HistStartPrimary","Energy at Wire vs Deposited Energy: Start: Primary",
						   501,-0.1,50,501,-0.1,50);
  recombinationIDEHistStartPrimary = tfs->make<TH2D>("recombinationIDEHistStartPrimary","Recombination vs Deposited Energy: Start: Primary",
						     1501,-0.01,15,100,0,1);

  recombinationHistStartSecondary =  tfs->make<TH1D>("recombinationHistStartSecondary","Histogram of Recombination: Start: Secondary",2000,0,1);
  recombination2HistStartSecondary = tfs->make<TH2D>("recombination2HistStartSecondary","Energy at Wire vs Deposited Energy: Start: Secondary",
						   501,-0.1,50,501,-0.1,50);
  recombinationIDEHistStartSecondary = tfs->make<TH2D>("recombinationIDEHistStartSecondary","Recombination vs Deposited Energy: Start: Secondary",
						     1501,-0.01,15,100,0,1);

}

void ana::ShowerRecombination::analyze(const art::Event& evt){
  
  evtcounter++;
  
  art::ServiceHandle<cheat::BackTrackerService> backtracker;

  //Getting  MC truth information                                                                                                                                                                                                            
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  if(evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
    {art::fill_ptr_vector(mclist, mctruthListHandle);}
  

  //Getting the SimWire Information                                                                                                                                                                                                          
  //Get the SimChannels so that we can find the IDEs deposited on them.                                                                                                                                                                      
  art::Handle<std::vector<sim::SimChannel> > simChannelHandle;
  std::vector<art::Ptr<sim::SimChannel> > simchannels;
  if(evt.getByLabel(fLArGeantModuleLabel,simChannelHandle))
    {art::fill_ptr_vector(simchannels, simChannelHandle);}

  //Getting the Hit Information                                                                                                                                                               #include "DetectorPropertiesStandard.h"    
  art::Handle<std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hits;
  if(evt.getByLabel(fHitsModuleLabel,hitListHandle))
    {art::fill_ptr_vector(hits, hitListHandle);}
  //std::cout << "hits.size(): " << hits.size() << std::endl;

  //List the particles in the event                                                                                                                                                                                                          
  const sim::ParticleList& particles = particleInventory->ParticleList();

  //Loop over the particles                                                                                                                                                                                                                  
  std::map<int,const simb::MCParticle*> trueParticles;

  double averageRecombination=0;
  double averageDep=0;
  double averageCol=0;
  double averageRecoCharge=0;
  int counter=0;
  //Make a map of Track id and pdgcode                                                                                                                                                                                                        
  for (sim::ParticleList::const_iterator particleIt = particles.begin(); particleIt != particles.end(); ++particleIt){
    const simb::MCParticle *particle = particleIt->second;
    trueParticles[particle->TrackId()] = particle;

    if (particle->Mother()==0){
      mytestfile<<"Event: "<<evtcounter<<", Track End Position: "<<particle->EndZ();		
    }
    /*
    if (particle->Mother()==1) {
      mytestfile<<"Event: "<<evtcounter<<" TrackID: "<<particle->TrackId()<<" particle: "<<particle->PdgCode()<<" process "<<
	particle->Process()<<" end z: "<<particle->EndZ()<<"\n";
      }
    */
  }
  
  double ADCtoElectrons[3] = {0.02354, 0.02130, 0.02354};
  //double ElectronstoADC = 1./0.02354;
  double finalz = 0;
    
  for (std::vector<art::Ptr<sim::SimChannel>>::iterator channelIt = simchannels.begin(); channelIt!=simchannels.end(); ++channelIt){
      
      
      auto tdc_ide_map = (*channelIt)->TDCIDEMap(); 
      
      for(auto const& tdc_ide_pair : tdc_ide_map) {
	double channelAverage=0;
	int channelCounter=0;
	auto const& ide_v = tdc_ide_pair.second;
	for(auto const& ide : ide_v) {
	  
	  //std::cout<<(*channelIt)->Channel()<<std::endl;
	  double energy = ide.energy;
	  double numberOfElectrons = ide.numElectrons;
	  double recombination = (numberOfElectrons*23.6)/(energy*1e6);
	  channelAverage += recombination;
	  channelCounter++;
	  
	  //if (numberOfElectrons>2*ElectronstoADC)
	  //{
	  //Fill Histograms;
          recombinationHist->Fill(recombination);

	  recombination2Hist->Fill(energy,numberOfElectrons*23.6/1e6);
	  averageRecombination+=recombination;

	  averageDep+=energy;
	  averageCol+=numberOfElectrons;

	  counter++;
	  
          recombinationIDEHist->Fill(energy,recombination);
	  
	  const simb::MCParticle *particle =  trueParticles[abs(ide.trackID)];//need abs to account for -1 track id particles
	  
	  if (ide.trackID==1) {
	    primary2Hist->Fill(energy,numberOfElectrons*23.6/1e6);
	    recombinationIDEHistPrimary->Fill(energy,recombination);
            //primaries<<"Recombination: "<<recombination<<" Energy Deposited: "<<energy<<", Electrons Energy: "<<numberOfElectrons<<" Track ID: "<<(*IDEIt).trackID<< " SimChannel: "<<(*channelIt)->Channel()<<" PDG Code: "<<particle->PdgCode()<<"\n";
	    if (finalz< ide.z)
	      {
		finalz=ide.z;
	      }
	  
	  }
	  else if (ide.trackID==-1) {
	    secondary2Hist->Fill(energy,numberOfElectrons*23.6/1e6);
	    recombinationIDEHistSecondary->Fill(energy,recombination);
            //secondaries<<"Recombination: "<<recombination<<" Energy Deposited: "<<energy<<", Electrons Energy: "<<numberOfElectrons<<" Track ID: "<<(*IDEIt).trackID<< " SimChannel: "<<(*channelIt)->Channel()<<"\n";

	  }
	  else  {
	    other2Hist->Fill(energy,numberOfElectrons*23.6/1e6);
            recombinationIDEHistOther->Fill(energy,recombination);
            //others<<"Recombination: "<<recombination<<" Energy Deposited: "<<energy<<", Electrons Energy: "<<numberOfElectrons<<" Track ID: "<<(*IDEIt).trackID<< " SimChannel: "<<(*channelIt)->Channel()<<" PDG Code: "<<particle->PdgCode()<<"\n";

	  }	    
	  
	  channel2Hist->Fill((*channelIt)->Channel(),recombination);

	  //const simb::MCParticle *particle = trueParticles.at((*IDEIt).trackID);
	  //if (recombination<0.1) myfile<<"Recombination: "<<recombination<<" Energy Deposited: "<<energy<<", Electrons Energy: "<<numberOfElectrons<<" Track ID: "<<(*IDEIt).trackID<<" Particle PDG: "<<particle->PdgCode()<<"\n";

          if (recombination>0.22 && recombination<0.38 && energy>2 && energy<12){
	    myfile1<<"Recombination: "<<recombination<<" Energy Deposited: "<<energy<<", Electrons Energy: "<<numberOfElectrons<<" Track ID: "<<ide.trackID<< " SimChannel: "<<(*channelIt)->Channel()<<" PDG Code: "<<particle->PdgCode()<<"\n";
	    recombination2HistLow->Fill(energy,numberOfElectrons*23.6/1e6);
	    myfile1<<"Event: "<<evtcounter<<", Track End Position: "<<particle->EndZ()<<", IDE position: "<<ide.z<<", PDG Code: "<<particle->PdgCode()<<"\n";

	  }
	  if (recombination<0.22 && energy>4){
	    myfile2<<"Recombination: "<<recombination<<" Energy Deposited: "<<energy<<", Electrons Energy: "<<numberOfElectrons<<" Track ID: "<<ide.trackID<< " SimChannel: "<<(*channelIt)->Channel()<<" PDG Code: "<<particle->PdgCode()<<"\n";
	    myfile2<<"Event: "<<evtcounter<<", Track End Position: "<<particle->EndZ()<<", IDE position: "<<ide.z<<"\n" ;

	    recombination2HistLowest->Fill(energy,numberOfElectrons*23.6/1e6);
	  }
	  
	  if ((particle->EndZ()-ide.z)<0.3) {
	    end2Hist->Fill(energy,numberOfElectrons*23.6/1e6);
	    endIDEHist->Fill(energy,recombination);
	    
	    //mytestfile<<"Recombination: "<<recombination<<" Energy Deposited: "<<energy<<", Electrons Energy: "<<numberOfElectrons<<
	    //    " Track ID: "<<ide.trackID<< " SimChannel: "<<(*channelIt)->Channel()<<" PDG Code: "<<particle->PdgCode()<<"\n";
	    //mytestfile<<"Event: "<<evtcounter<<", Track End Position: "<<particle->EndZ()<<", IDE position: "<<ide.z<<"\n" ;
	  
	  }
          if ((particle->EndZ()-ide.z)<0.3 && ide.trackID==1) {
            endIDEHistPrimary->Fill(energy,recombination);
          }

          if ((particle->EndZ()-ide.z)<0.3 && ide.trackID==-1) {
            endIDEHistSecondary->Fill(energy,recombination);
          }

          if ((particle->EndZ()-ide.z)<0.3 && abs(ide.trackID)!=1 && particle->PdgCode()==11) {
            endIDEHistOther->Fill(energy,recombination);
          }

          if (ide.z < (0.8*particle->EndZ())) {
	      recombinationHistStart->Fill(recombination);
	      recombination2HistStart->Fill(energy,numberOfElectrons*23.6/1e6);
	      recombinationIDEHistStart->Fill(energy,recombination);
	  }

          if (ide.z < (0.8*particle->EndZ()) && ide.trackID==1) {
	    recombinationHistStartPrimary->Fill(recombination);
	    recombination2HistStartPrimary->Fill(energy,numberOfElectrons*23.6/1e6);
	    recombinationIDEHistStartPrimary->Fill(energy,recombination);
          }

          if (ide.z < (0.8*particle->EndZ()) && ide.trackID==-1) {
            recombinationHistStartSecondary->Fill(recombination);
            recombination2HistStartSecondary->Fill(energy,numberOfElectrons*23.6/1e6);
            recombinationIDEHistStartSecondary->Fill(energy,recombination);
          }
	  

	  //}
	}
	//mytestfile<<", and energy: "<<energysum;
	channel2HistAverage->Fill((*channelIt)->Channel(),channelAverage/channelCounter);
      }      
    }

  //mytestfile<<" Penultimate momentum: "<<finalMomentum<<"\n";

  double averageRecoCol = 0;
  double averageRecoDep=0;
 
  for (std::vector<art::Ptr<recob::Hit>>::iterator hitIt = hits.begin(); hitIt != hits.end();++hitIt)
    { 
      double recoCharge = (*hitIt)->Integral();
      int view = (int)(*hitIt)->View() ;
      double recoElectrons = recoCharge/ADCtoElectrons[view];
      averageRecoCharge += recoElectrons;

      std::vector<sim::TrackIDE> trackIDEs = backtracker->HitToTrackIDEs(*hitIt);
      for (std::vector<sim::TrackIDE>::iterator IDEIt = trackIDEs.begin(); IDEIt!=trackIDEs.end();++IDEIt)
	{          
	  double IDEenergy = (*IDEIt).energy;
	  double IDEnumberOfElectrons = (*IDEIt).numElectrons;
	      
	  averageRecoCol += IDEnumberOfElectrons;
	  averageRecoDep += IDEenergy;
	}
    }
  
  //mytestfile<<"Dep: "<<averageDep*1e6/3.<<" Col: "<<averageCol/3.<<"\n";

  energyDepHist->Fill(averageDep/3.);
  energyColHist->Fill(averageCol*23.6/3./1e6);
  energyRecoHist->Fill(averageRecoCharge*23.6/3./1e6);

  energyRecoColHist->Fill(averageRecoCol*23.6/3./1e6);
  energyRecoDepHist->Fill(averageRecoDep/3.);
  
  recoRatioCol->Fill(averageRecoCol/averageCol);
  recoRatioDep->Fill(averageRecoDep/averageDep);
  recombinationAll->Fill((averageCol*23.6/1e6)/averageDep);
  recombinationReco->Fill((averageRecoCol*23.6/1e6)/averageRecoDep);

  recoEffDep->Fill((averageRecoCharge*23.6/1e6)/averageDep);
  recoEffCol->Fill(averageRecoCharge/averageCol);

  energyDepHist->SetTitle("IDE Energy Deposited");
  energyColHist->SetTitle("IDE Energy Collected at Wire");
  energyRecoHist->SetTitle("Reconstructed Energy");

  energyDepHist->GetXaxis()->SetTitle("IDE Energy Deposited (MeV)");
  energyColHist->GetXaxis()->SetTitle("IDE Energy Collected at Wire (MeV)");
  energyRecoHist->GetXaxis()->SetTitle("Reconstructed Energy (MeV)");


  //double electrondEdx = ;
  //double muondEdx = ;
  
  //double electronR = 0.7/(1+0.0486*electrondEdx/(0.5*1.3954));
  //double muonR = 0.7/(1+0.0486*muondEdx/(0.5*1.3954));

  //std::cout<<electronR<<" and "<<muonR<<std::endl;

  chargeEnergyHist->Fill(averageDep/3.,averageRecoCharge*23.6/1e6/3.);
  energyChargeHist->Fill(averageRecoCharge/3.,averageDep/3.);

  chargeEnergyHist->GetYaxis()->SetTitle("Total Energy Collected on Wires (MeV)");
  chargeEnergyHist->GetXaxis()->SetTitle("Total Energy Deposited by Track (MeV)");
  energyChargeHist->GetYaxis()->SetTitle("Total Energy Deposited by Track (MeV)");
  energyChargeHist->GetXaxis()->SetTitle("Total Charge Collected on Wires (ADC*Tick)");
  chargeEnergyHist->GetYaxis()->SetTitleOffset(1);  
  energyChargeHist->GetYaxis()->SetTitleOffset(1);

  recombination2Hist->GetYaxis()->SetTitle("Energy Collected on Wires (MeV)");
  recombination2Hist->GetXaxis()->SetTitle("Energy Deposited by IDE (MeV)");
  recombination2Hist->GetYaxis()->SetTitleOffset(1);

  primary2Hist->GetYaxis()->SetTitle("Energy Collected on Wires (MeV)");
  primary2Hist->GetXaxis()->SetTitle("Energy Deposited by IDE (MeV)");
  primary2Hist->GetYaxis()->SetTitleOffset(1);

  secondary2Hist->GetYaxis()->SetTitle("Energy Collected on Wires (MeV)");
  secondary2Hist->GetXaxis()->SetTitle("Energy Deposited by IDE (MeV)");
  secondary2Hist->GetYaxis()->SetTitleOffset(1);

  reconstruction2Hist->GetYaxis()->SetTitle("Reconstructed Energy Collected on Wires (MeV)");
  reconstruction2Hist->GetXaxis()->SetTitle("Energy Deposited by Hit (MeV)");
  reconstruction2Hist->GetYaxis()->SetTitleOffset(1); 

  channel2Hist->GetYaxis()->SetTitle("IDE Energy at Wire/IDE Energy Deposited");
  channel2Hist->GetXaxis()->SetTitle("Channel Number");
  channel2Hist->GetYaxis()->SetTitleOffset(1);

  recombinationIDEHist->GetYaxis()->SetTitle("Recombination Factor R");
  recombinationIDEHist->GetXaxis()->SetTitle("Deposited Energy (MeV)");
  recombinationIDEHist->GetYaxis()->SetTitleOffset(1);

  recombinationIDEHistPrimary->GetYaxis()->SetTitle("Recombination Factor R");
  recombinationIDEHistPrimary->GetXaxis()->SetTitle("Deposited Energy (MeV)");
  recombinationIDEHistPrimary->GetYaxis()->SetTitleOffset(1);

  recombinationIDEHistSecondary->GetYaxis()->SetTitle("Recombination Factor R");
  recombinationIDEHistSecondary->GetXaxis()->SetTitle("Deposited Energy (MeV)");
  recombinationIDEHistSecondary->GetYaxis()->SetTitleOffset(1);

  recombinationIDEHistOther->GetYaxis()->SetTitle("Recombination Factor R");
  recombinationIDEHistOther->GetXaxis()->SetTitle("Deposited Energy (MeV)");
  recombinationIDEHistOther->GetYaxis()->SetTitleOffset(1);

  endIDEHist->GetYaxis()->SetTitle("Recombination Factor R");
  endIDEHist->GetXaxis()->SetTitle("Deposited Energy (MeV)");
  endIDEHist->GetYaxis()->SetTitleOffset(1);

  endIDEHistPrimary->GetYaxis()->SetTitle("Recombination Factor R");
  endIDEHistPrimary->GetXaxis()->SetTitle("Deposited Energy (MeV)");
  endIDEHistPrimary->GetYaxis()->SetTitleOffset(1);

  endIDEHistOther->GetYaxis()->SetTitle("Recombination Factor R");
  endIDEHistOther->GetXaxis()->SetTitle("Deposited Energy (MeV)");
  endIDEHistOther->GetYaxis()->SetTitleOffset(1);

  recombinationIDEHistStart->GetYaxis()->SetTitle("Recombination Factor R");
  recombinationIDEHistStart->GetXaxis()->SetTitle("Deposited Energy (MeV)");
  recombinationIDEHistStart->GetYaxis()->SetTitleOffset(1);

  recombinationIDEHistStartSecondary->GetYaxis()->SetTitle("Recombination Factor R");
  recombinationIDEHistStartSecondary->GetXaxis()->SetTitle("Deposited Energy (MeV)");
  recombinationIDEHistStartSecondary->GetYaxis()->SetTitleOffset(1);

  recombinationIDEHistStartPrimary->GetYaxis()->SetTitle("Recombination Factor R");
  recombinationIDEHistStartPrimary->GetXaxis()->SetTitle("Deposited Energy (MeV)");
  recombinationIDEHistStartPrimary->GetYaxis()->SetTitleOffset(1);


}

void ana::ShowerRecombination::endJob(){
  mytestfile.close();
  myfile1.close();
  myfile2.close();
}


DEFINE_ART_MODULE(ana::ShowerRecombination)





