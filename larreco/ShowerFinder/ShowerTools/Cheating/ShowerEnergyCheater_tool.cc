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

      //Algorithm functions
      shower::TRACSCheatingAlg fTRACSCheatingAlg;

      //FCL
      art::InputTag fPFParticleModuleLabel;
      art::InputTag fHitModuleLabel;
 
 
			std::vector<double> plane2_energy;

			int showernum = 0;
			std::vector<std::pair<int,std::pair< int, unsigned int>>> hit_energy_vec;
			//std::vector<std::pair<int, std::pair<int, double>>> n_hit_energy;
			std::vector<std::tuple<int, int, int, int, double>> n_hit_energy;

 };


  ShowerEnergyCheater::ShowerEnergyCheater(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fTRACSCheatingAlg(pset.get<fhicl::ParameterSet>("TRACSCheatingAlg")),
    fPFParticleModuleLabel(pset.get<art::InputTag>("PFParticleModuleLabel","")),
    fHitModuleLabel(pset.get<art::InputTag>("HitModuleLabel"))
  {
  }

  ShowerEnergyCheater::~ShowerEnergyCheater()
  {
  }

  int ShowerEnergyCheater::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle, 
						   art::Event& Event, 
						   reco::shower::ShowerElementHolder& ShowerEleHolder){


    //const simb::MCParticle* trueParticle;

    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Shower Energy Cheater ~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

		showernum = ShowerEleHolder.GetShowerNumber();
		
		// get subrun number
		art::SubRunNumber_t subRunN = Event.subRun();
		// get event number 
		art::EventNumber_t EventN = Event.id().event();


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

      //if(view==2){
      

      //Get the hits
        std::vector<art::Ptr<recob::Hit> > hits = fmhc.at(cluster.key());

			

       showerHits.insert(showerHits.end(),hits.begin(),hits.end());
			 if(view == 2){
				std::cout << "shower num: " << showernum << "		view: " << view << "   shower hits: " << showerHits.size() << std::endl; // note, prints the view in seemingly random order 
     	 }
		 					
		
	}

		//Get the true particle from the shower - number specifies plane
    std::pair<int,double> ShowerTrackInfo = fTRACSCheatingAlg.TrueParticleIDFromTrueChain(showersMothers,showerHits,2);

		std::cout << "First:" << ShowerTrackInfo.first << " Second: " << ShowerTrackInfo.second << std::endl;
		
	
		
		const simb::MCParticle* trueParticle = trueParticles[ShowerTrackInfo.first];
		//if(/*trueParticle->PdgCode() == 11 ||*/ trueParticle->PdgCode() == 22){ 
		double trueEnergy = trueParticle->E() * 1000; // change to MeV
		std::cout << "shower num = " << showernum << "	true energy = " << trueEnergy << std::endl; 
		//std::cout << "trueParticles = " << *trueParticle << std::endl;
//	}
	
		
		if(ShowerTrackInfo.first==-99999) {
    	mf::LogError("ShowerEnergyCheater") << "True Shower Not Found";
      return 1;
    }
	

		//Holder for the final product
		std::vector<double> ShowerEnergyCheater;
		ShowerEnergyCheater.push_back(trueEnergy);

		std::vector<double> EnergyError = {-999,-999,-999};
		 
		ShowerEleHolder.SetElement(ShowerEnergyCheater,EnergyError,"ShowerEnergyCheater");



		double hitSize = showerHits.size();
		//plane2_energy.push_back(ShowerTrackInfo.second);
		//plane2_energy.push_back(std::make_pair(showernum, trueEnergy));
		n_hit_energy.push_back(std::make_tuple(subRunN ,EventN, showernum, hitSize, trueEnergy));
		// Make a .txt file with the truth energies from the plane
		std::ofstream outfile2;
		outfile2.open("truth_hit_energy_vec_spacecharge.txt");
		for (auto i = n_hit_energy.begin(); i != n_hit_energy.end(); ++i){
			//std::cout << i->first << "	" i->second->first << std::endl;

			outfile2 << std::get<0>(*i) << "		" << std::get<1>(*i) << "		" << std::get<2>(*i) << "		" << std::get<3>(*i) << "		" << std::get<4>(*i) << std::endl;
    }
	
	return 0;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerEnergyCheater)
