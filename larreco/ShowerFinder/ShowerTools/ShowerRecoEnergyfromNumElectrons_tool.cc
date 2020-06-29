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

//#include "larevt/SpaceChargeServices/SpaceChargeService.h"

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
    
    art::InputTag fPFParticleModuleLabel;

    //Services
    detinfo::DetectorProperties const* detprop = nullptr;
    art::ServiceHandle<geo::Geometry> fGeom;
    calo::CalorimetryAlg              fCalorimetryAlg;    

    bool fSCECorrectEField;

    // Declare variables etc.
    double Energy                = 0;
    double nominal_Efield        = 0;
    double localEfield           = 0;
    double localEfield_cweighted = 0;

    // vec to store subrun #, event #, shower #, # of hits and energy
    std::vector<std::tuple<int, int, int, int, double>> n_hit_energy; // more useful when making plots
		
    int showernum = 0;
	
};

    ShowerRecoEnergyfromNumElectrons::ShowerRecoEnergyfromNumElectrons(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fPFParticleModuleLabel(pset.get<art::InputTag>("PFParticleModuleLabel","")),
    detprop(lar::providerFrom<detinfo::DetectorPropertiesService>()),
    fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
    fSCECorrectEField(pset.get<bool>("SCECorrectEField"))
    //fSCE(lar::providerFrom<spacecharge::SpaceChargeService>())
    { 
    }

    ShowerRecoEnergyfromNumElectrons::~ShowerRecoEnergyfromNumElectrons()
    { 
    }

    int ShowerRecoEnergyfromNumElectrons::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
    art::Event& Event,
    reco::shower::ShowerElementHolder& ShowerEleHolder
    ){


    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Shower Reco Energy Tool ~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
		
    // get shower number per event
    showernum = ShowerEleHolder.GetShowerNumber();

    // get subrun number
    art::SubRunNumber_t subRunN = Event.subRun();

    // get event number 
    art::EventNumber_t EventN = Event.id().event();
		
    //ShowerEleHolder.PrintElements();

    //Holder for the final product
    std::vector<double> ShowerRecoEnergyfromNumElectrons;
    
    // Get the number of planes
    unsigned int numPlanes = fGeom->Nplanes();

    // Get the assocated pfParicle vertex PFParticles
    art::Handle<std::vector<recob::PFParticle> > pfpHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
        throw cet::exception("ShowerRecoEnergyfromNumElectrons") << "Could not get the pandora pf particles. Something is not configured correctly Please give the correct pandora module label. Stopping";
    return 1;
    }

    std::map<geo::View_t, std::vector<art::Ptr<recob::Hit> > > view_hits;

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

    /*
    // variables to hold position of spacepoints
    double X = 0, Y = 0, Z = 0;
    // loop over the space points 
    for (auto const &sps : spacepoints){
        // sum up all the XYZ positions of the spacepoints 
        X += sps->XYZ()[0];
        Y += sps->XYZ()[1];
        Z += sps->XYZ()[2];
    }

    // Get the average spacepoint postions - (consider charge weighting this in the future)
    X = X/spacepoints.size();
    Y = Y/spacepoints.size();
    Z = Z/spacepoints.size();
    std::cout << "number of spacepoints = " << spacepoints.size() << std::endl;
    std::cout << "X = " << X << "       Y = " << Y << "     Z = " << Z << std::endl;
    */
    
    // Get shower centre (remove above at some point)
    TVector3 showercentre = IShowerTool::GetTRACSAlg().ShowerCentre(spacepoints);
    std::cout << "shower centre = (" << showercentre[0] << ", " << showercentre[1] << ", " << showercentre[2] << ")" << std::endl;    

     // Get E-field -  Electric Field in the drift region in KV/cm
    nominal_Efield  = detprop->Efield(); // Nominal E-field
    std::cout << "EField = " << nominal_Efield << std::endl;
    
    // Get hit association for spacepoints
    art::FindManyP<recob::Hit> fmhsp(spHandle, Event, fPFParticleModuleLabel);

    // Get the charge weighted shower centre - sum(charge*position)/sum(charge)
    // Charge is obtained from the lifetime corrected hits and hits are obtained from the spacepoints 
    TVector3 chargecentre = IShowerTool::GetTRACSAlg().ShowerCentre(spacepoints, fmhsp);
    std::cout << "charge centre = " << chargecentre[0] << ", " << chargecentre[1] << ", " << chargecentre[2] << std::endl;

    // Enable SCE 
    fSCECorrectEField = true;
    // Get space charge corrected E-field
    if(fSCECorrectEField){
        // Check the shower centre exists and is sensible
        // Think the charge weighted centre can take nan values which causes the E-field calculation to fail
        if(showercentre.Mag() >= 0 && chargecentre.Mag() >=0){ 
            localEfield = IShowerTool::GetTRACSAlg().SCECorrectEField(nominal_Efield, showercentre);
            std::cout << "local efield = " << localEfield << std::endl;
            localEfield_cweighted = IShowerTool::GetTRACSAlg().SCECorrectEField(nominal_Efield, chargecentre);
            std::cout << "local efield weighted = " << localEfield_cweighted << std::endl;
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

    //Get the hits
    std::vector<art::Ptr<recob::Hit> > hits = fmhc.at(cluster.key());
    if(hits.size() == 0){
        mf::LogWarning("ShowerRecoEnergyfromNumElectrons") << "No hit for the cluster. This suggest the find many is wrong."<< std::endl;
        continue;
    }

    //Get the view.
    geo::View_t view = cluster->View();
    //std::cout << "view = " << view << "	hits = " << hits.size() << std::endl;

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
        //std::cout << "hits = " << hits.size() << std::endl;
		
        // Print out the energy for each plane
        std::cout<<"View: "<< view <<  " and energy: "<<Energy<<std::endl;;

       
        unsigned int viewNum = view;
        view_energies[viewNum] = Energy;

    }

    //TODO think of a better way of doing this
    for (unsigned int plane=0; plane<numPlanes; ++plane) {
       
    try{
        Energy = view_energies.at(plane);
        if (Energy<0){
            mf::LogWarning("ShowerLinearEnergy") << "Negative shower energy: "<<Energy;
            Energy=-999;
        }
        if(plane == 2){
            n_hit_energy.push_back(std::make_tuple(subRunN ,EventN, showernum, hits.size(), Energy)); //save info for collection plane
        }   
        } 

    catch(...){
    mf::LogWarning("ShowerLinearEnergy") <<"No energy calculation for plane "<<plane<<std::endl;
    // if there's no calculation, set the energy to -999.
    Energy = -999;
    if(plane == 2){
        n_hit_energy.push_back(std::make_tuple(subRunN, EventN, showernum, hits.size(), Energy)); //save info for collection plane
    }
    }              	
    ShowerRecoEnergyfromNumElectrons.push_back(Energy);   
    }
 
    if(ShowerRecoEnergyfromNumElectrons.size() == 0){
        throw cet::exception("ShowerLinearEnergy") << "Energy Vector is empty";
        return 1;
        }

    std::vector<double> EnergyError = {-999,-999,-999};
	 
    ShowerEleHolder.SetElement(ShowerRecoEnergyfromNumElectrons,EnergyError,"ShowerEnergy");
	

    bool write_to_file = true; 
    // Make a .txt file with the subrun, event number, showernum, hits and energies from the plane
    // Useful info when making plots
    if(write_to_file){
        std::ofstream outfile;
        outfile.open("reco_hit_energy_vec_corrected.txt");
        for (auto i = n_hit_energy.begin(); i != n_hit_energy.end(); ++i){	
            outfile << std::get<0>(*i) << "		" << std::get<1>(*i) << "		" << std::get<2>(*i) << "		" << std::get<3>(*i) << "		" << std::get<4>(*i) << std::endl;
        }
    }
    return 0;

}


// function to calculate the reco energy
double ShowerRecoEnergyfromNumElectrons::CalculateEnergy(std::vector<art::Ptr<recob::Hit> >& hits, geo::View_t& view){
 
    double totalCharge = 0;
    double totalEnergy = 0;
    double correctedtotalCharge = 0;
    double nElectrons = 0;
    double nominal_recombination = 0.64; //constant factor for every shower (study done by Ed)	
    double recombination = 0;
    
    double A = 0.8;
    double nominal_k_dEdx = ((A / nominal_recombination) - 1) * nominal_Efield; // see https://arxiv.org/pdf/1306.1712.pdf


    if(fSCECorrectEField){
        recombination = A / (1 + (nominal_k_dEdx / localEfield_cweighted));
    }
    else{
        recombination = nominal_recombination;
    }   
    std::cout << "recombination = " << recombination << std::endl;
    
    for (art::PtrVector<recob::Hit>::const_iterator hit = hits.begin(); hit != hits.end(); ++hit){
        totalCharge += (*hit)->Integral() * fCalorimetryAlg.LifetimeCorrection((*hit)->PeakTime()); // obtain charge and correct for lifetime
    }
 		
    // correct charge due to recombination
    correctedtotalCharge = totalCharge / recombination;
    // calculate # of electrons and the corresponding energy
    nElectrons = fCalorimetryAlg.ElectronsFromADCArea(correctedtotalCharge, view);
    totalEnergy = (nElectrons / util::kGeVToElectrons) * 1000; // energy in MeV
    return totalEnergy;
 
    }


}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerRecoEnergyfromNumElectrons)


