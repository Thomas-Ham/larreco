//############################################################################
//### Name:        ShowerEnergyCheater                                     ###
//### Author:      Ed Tyley                                                ###
//### Date:        16.07.19                                                ###
//### Description: Cheating tool using truth for shower energy             ###
//############################################################################

#include "larreco/ShowerFinder/ShowerTools/IShowerTool.h"

//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft Includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larreco/RecoAlg/TRACSAlg.h"
#include "larreco/RecoAlg/TRACSCheatingAlg.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "nusimdata/SimulationBase/MCTruth.h"

//ROOT Inclides
#include "TFile.h"
#include "TTree.h"

//C++ Includes
#include <iostream>
#include <cmath>
#include <fstream>

//Root Includes

namespace ShowerRecoTools {

  class ShowerEnergyCheater:IShowerTool {

    public:

      ShowerEnergyCheater(const fhicl::ParameterSet& pset);

      ~ShowerEnergyCheater();

      //Calculate Cheating StartPosition 
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
			   art::Event& Event, 
			   reco::shower::ShowerElementHolder& ShowerEleHolder) override;

    private:

        void FindLargestShower(std::vector<std::tuple<unsigned int, unsigned int, unsigned int, int, double>> energy_largest_shower);     

        //Algorithm functions
        shower::TRACSCheatingAlg fTRACSCheatingAlg;

        //FCL
        art::InputTag fPFParticleModuleLabel;
        art::InputTag fHitModuleLabel;
 
        unsigned int showernum = 0;
        art::SubRunNumber_t subRunN;                                                                                             
        art::EventNumber_t EventN; 
        int hitsize;
        double trueEnergy;

	std::vector<std::tuple<unsigned int, unsigned int, unsigned int, int, double>> n_hit_energy;

        TFile *m_file;    ///< TFile for saving event information                                                               
        TTree *fOutputTree;  
 };


    ShowerEnergyCheater::ShowerEnergyCheater(const fhicl::ParameterSet& pset) :
        IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
        fTRACSCheatingAlg(pset.get<fhicl::ParameterSet>("TRACSCheatingAlg")),
        fPFParticleModuleLabel(pset.get<art::InputTag>("PFParticleModuleLabel","")),
        fHitModuleLabel(pset.get<art::InputTag>("HitModuleLabel"))
    {
        // The TFileService lets us define a tree and takes care of writing it to file          
        art::ServiceHandle<art::TFileService> tfs;                                                                                                                       
        m_file      = new TFile("EnergyfromNumElectrons.root", "UPDATE");                          
        fOutputTree = tfs->make<TTree>("energy_cheater","Energy Cheater");                                                          

        //add branches                                                                      
        fOutputTree->Branch("Subrun", &subRunN, "Subrun/i");                                          
        fOutputTree->Branch("Event", &EventN, "Event/i");                                      
        fOutputTree->Branch("ShowerN", &showernum, "ShowerN/i");                                        
        fOutputTree->Branch("NHits", &hitsize, "NHits/I");                      
        fOutputTree->Branch("Energy", &trueEnergy, "Energy/d");   
    }

    ShowerEnergyCheater::~ShowerEnergyCheater()
    {
        m_file->cd();                                                                                   
        fOutputTree->CloneTree()->Write("energy_cheater", TObject::kOverwrite); 

        bool find_largest_shower = true;                                                                                         
        if(find_largest_shower){                                                                            
            FindLargestShower(n_hit_energy);                                                                                                    
        }     

        m_file->Close();
    }

    int ShowerEnergyCheater::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle, 
						   art::Event& Event, 
						   reco::shower::ShowerElementHolder& ShowerEleHolder){

    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Shower Energy Cheater ~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

    showernum = ShowerEleHolder.GetShowerNumber();
		
    // get subrun number
    subRunN = Event.subRun();
    // get event number 
    EventN = Event.id().event();

    //Could store these in the shower element holder and just calculate once?
    std::map<int,const simb::MCParticle*> trueParticles = fTRACSCheatingAlg.GetTrueParticleMap();
    std::map<int,std::vector<int> > showersMothers = fTRACSCheatingAlg.GetTrueChain(trueParticles);

    //Get the hits from the shower:
    art::Handle<std::vector<recob::PFParticle> > pfpHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
      throw cet::exception("ShowerEnergyCheater") << "Could not get the pandora pf particles. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
      return 1;
    }

    //Get the clusters
    art::Handle<std::vector<recob::Cluster> > clusHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, clusHandle)){
      throw cet::exception("ShowerEnergyCheater") << "Could not get the pandora clusters. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
      return 1;
    }

    art::FindManyP<recob::Cluster> fmc(pfpHandle, Event, fPFParticleModuleLabel);
    std::vector<art::Ptr<recob::Cluster> > clusters = fmc.at(pfparticle.key());

    //Get the hit association
    art::FindManyP<recob::Hit> fmhc(clusHandle, Event, fPFParticleModuleLabel);

    std::vector<art::Ptr<recob::Hit> > showerHits;
    //std::map<geo::View_t, std::vector<art::Ptr<recob::Hit> >> showerHits;	

    for(auto const& cluster: clusters){
        // Get the view.
        geo::View_t view = cluster->View();

        //Get the hits
        showerHits = fmhc.at(cluster.key());

        //showerHits.insert(showerHits.end(),hits.begin(),hits.end());
        if(view == 2){
	    hitsize = showerHits.size();
     	}
    }

    //Get the true particle from the shower - number specifies plane
    // Gets (shower id, deposited shower energy) in the plane
    //std::cout << "hitsize = " << hitsize << std::endl;
    std::pair<int,double> ShowerTrackInfo = fTRACSCheatingAlg.TrueParticleIDFromTrueChain(showersMothers,showerHits,2);
       
    // Get the energy 
    if(ShowerTrackInfo.first==-99999) {
        mf::LogError("ShowerEnergyCheater") << "True Shower Not Found. Setting hits and energy to -999.";
        hitsize    = -999;
        trueEnergy = -999;
        //return 1;
    }
    else{
        trueEnergy = 0;
        // Get true energy of showering particle 
        const simb::MCParticle* trueParticle = trueParticles[ShowerTrackInfo.first];
        trueEnergy = trueParticle->E() * 1000; // change to MeV
    }
    std::cout << "Subrun = " << subRunN << "    Event = " << EventN << "    shower num = " << showernum << "    hits = " << hitsize <<  "    true energy = " << trueEnergy << std::endl; 
		
    //Holder for the final product
    std::vector<double> ShowerEnergyCheater;
    ShowerEnergyCheater.push_back(trueEnergy);

    std::vector<double> EnergyError = {-999,-999,-999};
	 
    ShowerEleHolder.SetElement(ShowerEnergyCheater,EnergyError,"ShowerEnergyCheater");

    // check showers are coming from a photon/electron        
    // if(abs(trueParticle->PdgCode()) == 11 || trueParticle->PdgCode() == 22){ 
        n_hit_energy.push_back(std::make_tuple(subRunN ,EventN, showernum, hitsize, trueEnergy));
    // }

    fOutputTree->Fill();
	
    return 0;
    }


// Function to find only the largest shower
void ShowerEnergyCheater::FindLargestShower(std::vector<std::tuple<unsigned int, unsigned int, unsigned int, int, double>> energy_largest_shower){
    // Cut showers so we are only left with the biggest (most hits) from each event.
    // Feel like there should be a more ROOT-y way to do this...
    unsigned int i = 0;                                                                 
    while(i < energy_largest_shower.size()){  

        // If the shower at (i-1)'th position from the same event has fewer hits, delete it
        if(std::get<2>(energy_largest_shower[i]) != 0 && std::get<3>(energy_largest_shower[i]) > std::get<3>(energy_largest_shower[(i-1)])){   

            energy_largest_shower.erase(energy_largest_shower.begin() + (i-1));   
        }
        // Delete any remaining i'th non primary shower (i.e. the non primary shower which has fewer hits than the one at (i-1))
        else if(std::get<2>(energy_largest_shower[i]) != 0){ 
    
            energy_largest_shower.erase(energy_largest_shower.begin() + (i));
        }

        else{
            i++;
        }
    }


    // Add new tree to root file which has only the largest shower from each event
    TTree *energy_cheater_LS = new TTree("energy_cheater_LS", "Energy Cheater Largest Shower");

    energy_cheater_LS->Branch("Subrun", &subRunN, "Subrun/i");                                                                              
    energy_cheater_LS->Branch("Event", &EventN, "Event/i");    
    energy_cheater_LS->Branch("ShowerN", &showernum, "ShowerN/i");
    energy_cheater_LS->Branch("NHits", &hitsize, "NHits/I");                                                     
    energy_cheater_LS->Branch("Energy", &trueEnergy, "Energy/d"); 

    for(unsigned int i = 0; i < energy_largest_shower.size(); i++){
        subRunN    = std::get<0>(energy_largest_shower[i]);
        EventN     = std::get<1>(energy_largest_shower[i]);
        showernum  = std::get<2>(energy_largest_shower[i]);
        hitsize    = std::get<3>(energy_largest_shower[i]);
        trueEnergy = std::get<4>(energy_largest_shower[i]);

        energy_cheater_LS->Fill();
    }

    m_file->Write("", TObject::kOverwrite); // save only the new version of the tree - without the arguments was duplicating original tree

}

} // namespace 

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerEnergyCheater)
