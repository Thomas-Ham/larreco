//############################################################################
//### Name:        ShowerRecoEnergyfromNumElectrons                        ###
//### Author:      Tom Ham			                           ###
//### Date:        01/04/2020                                              ###
//### Description: Tool for finding the Energy of the shower by going      ###
//###              from number of hits -> number of electrons -> energy.   ###
//###              Derived from the linear energy algorithm, written for   ###
//###              the EMShower_module.cc                                  ###
//############################################################################

#include "larreco/ShowerFinder/ShowerTools/IShowerTool.h"

//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art_root_io/TFileService.h"

//LArSoft Includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h" 

//C++ Includes
#include <iostream>
#include <vector>
#include <tuple>

namespace ShowerRecoTools {

  class ShowerRecoEnergyfromNumElectrons:IShowerTool {

    public:

    ShowerRecoEnergyfromNumElectrons(const fhicl::ParameterSet& pset);
    
    ~ShowerRecoEnergyfromNumElectrons(); 
        
    //Physics Function. Calculate the shower Energy.
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
    art::Event& Event,
    reco::shower::ShowerElementHolder& ShowerElementHolder
    ) override;

    private:
    
    double CalculateEnergy(std::vector<art::Ptr<recob::Hit> >& hits, geo::View_t& view);
    double ModBoxRecombination(double EField);
    void FindLargestShower(std::vector<std::tuple<unsigned int, unsigned int, unsigned int, int, double>> recoenergy_largest_shower);    
 
    art::InputTag fPFParticleModuleLabel;

    //Services
    detinfo::DetectorProperties const* detprop = nullptr;
    art::ServiceHandle<geo::Geometry> fGeom;
    calo::CalorimetryAlg              fCalorimetryAlg;    
   
    //fcl params 
    std::string fShowerEnergyfromNumElectrons; 
    bool fSCECorrectEField;
    double fNominalModBoxdEdx;
    double fNominalRecombinationFactor;

    // Declare variables etc.
    double Energy                = 0;
    double nominal_Efield        = 0;
    double localEfield           = 0;
    double localEfield_cweighted = 0;
    
    unsigned int showernum = 0;
    art::SubRunNumber_t subRunN;
    art::EventNumber_t EventN;
    int hitsize;

    std::vector<std::tuple<int, double>> hit_key_EField_sp;

    // vec to store subrun #, event #, shower #, # of hits and energy
    std::vector<std::tuple<unsigned int, unsigned int, unsigned int, int, double>> n_hit_energy; // more useful when making plots

    bool write_to_file = false;

    TFile *m_file;    ///< TFile for saving event information  
    TTree *fOutputTree;

}; // class

    ShowerRecoEnergyfromNumElectrons::ShowerRecoEnergyfromNumElectrons(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fPFParticleModuleLabel(pset.get<art::InputTag>("PFParticleModuleLabel","")),
    detprop(lar::providerFrom<detinfo::DetectorPropertiesService>()),
    fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
    fShowerEnergyfromNumElectrons(pset.get<std::string>("ShowerEnergyfromNumElectrons")),
    fSCECorrectEField(pset.get<bool>("SCECorrectEField")),
    fNominalModBoxdEdx(pset.get<double>("NominalModBoxdEdx")),
    fNominalRecombinationFactor(pset.get<double>("NominalRecombinationFactor"))
    {

        if(write_to_file){
            // The TFileService lets us define a tree and takes care of writing it to file
            art::ServiceHandle<art::TFileService> tfs;
            m_file      = new TFile("Energy_files/Anode_sample_SCE.root", "UPDATE"); 
            fOutputTree = tfs->make<TTree>("recoenergy_SCE_correction_ModBox","Reco Energy SCE correction ModBox");

            //add branches                                                                          
            fOutputTree->Branch("Subrun", &subRunN, "Subrun/i");
            fOutputTree->Branch("Event", &EventN, "Event/i");
            fOutputTree->Branch("ShowerN", &showernum, "ShowerN/i");
            fOutputTree->Branch("NHits", &hitsize, "NHits/I");
            fOutputTree->Branch("Energy", &Energy, "Eventd");
        }
    }

    ShowerRecoEnergyfromNumElectrons::~ShowerRecoEnergyfromNumElectrons()
    {
        if(write_to_file){
            //store output tree                                                                                            
            m_file->cd();                                                  
            fOutputTree->CloneTree()->Write("recoenergy_SCE_correction_ModBox", TObject::kOverwrite);

            // Find only the largest shower from each event and write to the output file
            // To find the largest shower need to compare with the previous one and since the tool iterates over them one at a time 
            // we need to run at the end once we've looked at all the showers (hence we're in the destructor since it runs at the end).
            bool find_largest_shower = true;
            if(find_largest_shower){
                FindLargestShower(n_hit_energy);
            }

            m_file->Close();
        }
    }

    int ShowerRecoEnergyfromNumElectrons::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
    art::Event& Event,
    reco::shower::ShowerElementHolder& ShowerEleHolder
    ){

    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Shower Reco Energy Tool ~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

    // get shower number per event
    showernum = ShowerEleHolder.GetShowerNumber();

    // get subrun number
    subRunN = Event.subRun();

    // get event number 
    EventN = Event.id().event();
    
    //ShowerEleHolder.PrintElements();

    // Get the number of planes
    unsigned int numPlanes = fGeom->Nplanes();

    //Holder for the final product
    std::vector<double> ShowerRecoEnergyfromNumElectrons(numPlanes, -999);
    std::vector<double> ShowerRecoEnergyfromNumElectronsError(numPlanes, -999);

    // Get the assocated pfParicle vertex PFParticles
    art::Handle<std::vector<recob::PFParticle> > pfpHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
        throw cet::exception("ShowerRecoEnergyfromNumElectrons") << "Could not get the pandora pf particles. Something is not configured correctly Please give the correct pandora module label. Stopping";
    return 1;
    }

    std::map<geo::View_t, std::vector<art::Ptr<recob::Hit>>> view_hits;

    //Get the clusters
    art::Handle<std::vector<recob::Cluster> > clusHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, clusHandle)){
      throw cet::exception("ShowerRecoEnergyfromNumElectrons") << "Could not get the pandora clusters. Something is not configured correctly Please give the correct pandora module label. Stopping";
      return 1;
    }
    art::FindManyP<recob::Cluster> fmc(pfpHandle, Event, fPFParticleModuleLabel);
    std::vector<art::Ptr<recob::Cluster> > clusters = fmc.at(pfparticle.key());


    //Get the sapcepoints
    art::Handle<std::vector<recob::SpacePoint> > spHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, spHandle)){
        throw cet::exception("thShowerSpacePoints") << "Could not get the pandora space points. Something is not cofingured coreectly Please give the correct pandoa module     label. Stopping";
    return 1;
    }
    art::FindManyP<recob::SpacePoint> fmsp(pfpHandle, Event, fPFParticleModuleLabel);
    std::vector<art::Ptr<recob::SpacePoint> > spacepoints = fmsp.at(pfparticle.key());

    // Get E-field -  Electric Field in the drift region in KV/cm
    nominal_Efield  = detprop->Efield(); // Nominal E-field 

    // Get hit association for spacepoints
    int hits_total = 0;
    hit_key_EField_sp.clear();
    art::FindManyP<recob::Hit> fmhsp(spHandle, Event, fPFParticleModuleLabel);
    for(auto const& spacepoint : spacepoints){
        //Get the hits
        std::vector<art::Ptr<recob::Hit> > hits = fmhsp.at(spacepoint.key());
        hits_total += hits.size();

        // get local Efield at SP_pos
        TVector3 SP_pos = IShowerTool::GetTRACSAlg().SpacePointPosition(spacepoint); 
        double localEfield_SP = IShowerTool::GetTRACSAlg().SCECorrectEField(nominal_Efield, SP_pos);

        for(auto const& h : hits){
            hit_key_EField_sp.push_back(std::make_tuple(h.key(), localEfield_SP));
        }

                
    }
    std::cout << "hits associated with SP: " << hits_total << std::endl;
 
    for(auto i = hit_key_EField_sp.begin(); i != hit_key_EField_sp.end(); i++){
        std::cout << std::get<0>(*i) << "   " << std::get<1>(*i) << std::endl;
    }

    // Not all hits from the clusters will match to Space points. If that's the case let's use some sort of shower centre
    // Get shower centre 
    TVector3 showercentre = IShowerTool::GetTRACSAlg().ShowerCentre(spacepoints);
  
    // Get the charge weighted shower centre - sum(charge*position)/sum(charge)
    // Charge is obtained from the lifetime corrected hits and hits are obtained from the spacepoints 
    TVector3 chargecentre = IShowerTool::GetTRACSAlg().ShowerCentre(spacepoints, fmhsp);

    // Get space charge corrected E-field
    if(fSCECorrectEField){
        // Check the shower centre exists and is sensible
        // Think the charge weighted centre can take nan values which causes the E-field calculation to fail
        if(showercentre.Mag() >= 0 && chargecentre.Mag() >=0){ 
            localEfield = IShowerTool::GetTRACSAlg().SCECorrectEField(nominal_Efield, showercentre);
            localEfield_cweighted = IShowerTool::GetTRACSAlg().SCECorrectEField(nominal_Efield, chargecentre);
        }
        else{
            mf::LogWarning("ShowerRecoEnergyfromNumElectrons") << "Shower centre calculation doesn't look to be sensible. Reconstruction is probably dodgy." << std::endl;
        }

    }

    //Get the hit association for clusters
    art::FindManyP<recob::Hit> fmhc(clusHandle, Event, fPFParticleModuleLabel);

    std::vector<std::vector<art::Ptr<recob::Hit> > > trackHits;
    trackHits.resize(numPlanes);
    
    //Loop over the clusters in the plane and get the hits
    for(auto const& cluster: clusters){
        std::cout << "cluster key: " << cluster.key() << std::endl;
        //Get the hits
        std::vector<art::Ptr<recob::Hit> > hits = fmhc.at(cluster.key());

        if(hits.size() == 0){
            mf::LogWarning("ShowerRecoEnergyfromNumElectrons") << "No hit for the cluster. This suggest the find many is wrong."<< std::endl;
            continue;
        }

        //Get the view.
        geo::View_t view = cluster->View();

        view_hits[view].insert(view_hits[view].end(),hits.begin(),hits.end());
    }




    std::map<unsigned int, double > view_energies;
    std::vector<art::Ptr<recob::Hit>> hits;
    //Accounting for events crossing the cathode.
    for(auto const& view_hit_iter: view_hits){
 	
        hits = view_hit_iter.second;
        geo::View_t view = view_hit_iter.first;
 	
        //Calculate the Energy
        Energy = CalculateEnergy(hits,view);
        
        hitsize = hits.size();
		
        // Print out the energy for each plane
        std::cout << "Subrun: " << subRunN << " Event: " << EventN << " ShowerNum: " << showernum << " View: "<< view << "  hitsize: " << hitsize  <<  " energy: " << Energy << std::endl;

       
        unsigned int viewNum = view;
        view_energies[viewNum] = Energy;
        
        ShowerRecoEnergyfromNumElectrons.at(view) = Energy;   

    }

    //TODO think of a better way of doing this
    for (unsigned int plane=0; plane<numPlanes; ++plane) {
       
        try{
            Energy = view_energies.at(plane);
            if (Energy<0){
                mf::LogWarning("ShowerLinearEnergy") << "Negative shower energy: "<<Energy << ". Setting the energy to -999." << std::endl;
                Energy=-999;
             }
            if(plane == 2){
                n_hit_energy.push_back(std::make_tuple(subRunN ,EventN, showernum, hitsize, Energy)); //save info for collection plane
            }   
        } 

        catch(...){
            mf::LogWarning("ShowerLinearEnergy") <<"No energy calculation for plane "<< plane << ". Setting the energy to -999." << std::endl;
            // if there's no calculation, set the energy to -999.
            Energy = -999;
            if(plane == 2){
                n_hit_energy.push_back(std::make_tuple(subRunN, EventN, showernum, hitsize, Energy)); //save info for collection plane
            }
        }              	
    }
 
    if(ShowerRecoEnergyfromNumElectrons.size() == 0){
        throw cet::exception("ShowerLinearEnergy") << "Energy Vector is empty";
        return 1;
        }

	 
    ShowerEleHolder.SetElement(ShowerRecoEnergyfromNumElectrons,ShowerRecoEnergyfromNumElectronsError,fShowerEnergyfromNumElectrons);
    if(write_to_file){
        fOutputTree->Fill();
    }

    return 0;

}


// function to calculate the reco energy
double ShowerRecoEnergyfromNumElectrons::CalculateEnergy(std::vector<art::Ptr<recob::Hit> >& hits, geo::View_t& view){
 
    double totalCharge = 0;
    double totalEnergy = 0;
    double nElectrons = 0;
    double recombination = 0;
    
    // double A = 0.8;
    // double nominal_k_dEdx = ((A / fNominalRecombinationFactor) - 1) * nominal_Efield; // Birk's method, see https://arxiv.org/pdf/1306.1712.pdf

    // Loop over the hits
    for(auto const& h : hits){
        const art::Ptr<recob::Hit> hit = hits[0]; 
        if(fSCECorrectEField){
            auto it = std::find_if(hit_key_EField_sp.begin(), hit_key_EField_sp.end(), [&](const std::tuple<unsigned int, double>& e){
                return std::get<0>(e) == h.key();
            });
    
            // do a per hit SCE correction since we have associated space points
            if(it != hit_key_EField_sp.end()){
                std::cout << "hit has SP " << std::endl;
                recombination = ModBoxRecombination(std::get<1>(*it));
            }

            // no associated space points :( so correct these hits using the shower centre position 
            else{
                std::cout << "hit has no SP " << std::endl;
                
                recombination = ModBoxRecombination(localEfield_cweighted);
            }
        }

        else{
            recombination = fNominalRecombinationFactor;
            std::cout << "Using the nominal Efield of " << nominal_Efield << "kV/cm and the nominal recombination factor of " << recombination << "." << std::endl;
        }   
 

    totalCharge += (hit)->Integral() * fCalorimetryAlg.LifetimeCorrection((hit)->PeakTime()) / recombination; // obtain charge and correct for lifetime and recombination
    }

    // calculate # of electrons and the corresponding energy
    nElectrons = fCalorimetryAlg.ElectronsFromADCArea(totalCharge, view);
    totalEnergy = (nElectrons / util::kGeVToElectrons) * 1000; // energy in MeV
    return totalEnergy;
 
}

// ModBox function to calculate recomobination
double ShowerRecoEnergyfromNumElectrons::ModBoxRecombination(double EField){
    double rho   = detprop->Density();
    double Alpha = util::kModBoxA;
    double Beta  = util::kModBoxB/(rho * EField);

    double recombination = 0;
    recombination = std::log(Alpha + Beta * fNominalModBoxdEdx)/(Beta * fNominalModBoxdEdx);
    std::cout << "Using the charge weighted local Efield of " << localEfield_cweighted << " kV/cm and a recombination factor of " << recombination << "." << std::endl;

    return recombination; 
}

// Function to find only the largest shower
void ShowerRecoEnergyfromNumElectrons::FindLargestShower(std::vector<std::tuple<unsigned int, unsigned int, unsigned int, int, double>> recoenergy_largest_shower){
    // Cut showers so we are only left with the biggest (most hits) from each event.
    // Feel like there should be a more ROOT-y way to do this...
    unsigned int i = 0;                                                                 
    while(i < recoenergy_largest_shower.size()){  

        // If the shower at (i-1)'th position from the same event has fewer hits, delete it
        if(std::get<2>(recoenergy_largest_shower[i]) != 0 && std::get<3>(recoenergy_largest_shower[i]) > std::get<3>(recoenergy_largest_shower[(i-1)])){   

            recoenergy_largest_shower.erase(recoenergy_largest_shower.begin() + (i-1));   
        }
        // Delete any remaining i'th non primary shower (i.e. the non primary shower which has fewer hits than the one at (i-1))
        else if(std::get<2>(recoenergy_largest_shower[i]) != 0){ 
    
            recoenergy_largest_shower.erase(recoenergy_largest_shower.begin() + (i));
        }

        else{
            i++;
        }
    }


    // Add new tree to root file which has only the largest shower from each event
    TTree *recoenergy_LS = new TTree("recoenergy_SCE_correction_ModBoxLS", "Reco Energy SCE Correction ModBox Largest Shower");
    //recoenergy_LS    = tfs->make<TTree>("recoenergy_LS","Reco Energy Largest Shower");

    recoenergy_LS->Branch("Subrun", &subRunN, "Subrun/i");                                                                              
    recoenergy_LS->Branch("Event", &EventN, "Event/i");    
    recoenergy_LS->Branch("ShowerN", &showernum, "ShowerN/i");
    recoenergy_LS->Branch("NHits", &hitsize, "NHits/I");                                                     
    recoenergy_LS->Branch("Energy", &Energy, "Energy/d"); 

    for(unsigned int i = 0; i < recoenergy_largest_shower.size(); i++){
        subRunN   = std::get<0>(recoenergy_largest_shower[i]);
        EventN    = std::get<1>(recoenergy_largest_shower[i]);
        showernum = std::get<2>(recoenergy_largest_shower[i]);
        hitsize   = std::get<3>(recoenergy_largest_shower[i]);
        Energy    = std::get<4>(recoenergy_largest_shower[i]);

        recoenergy_LS->Fill();
    }

    m_file->Write("", TObject::kOverwrite); // save only the new version of the tree - without the arguments was duplicating original tree
}

} // namespace

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerRecoEnergyfromNumElectrons)










