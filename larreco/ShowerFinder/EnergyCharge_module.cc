// Framework includes                                                                                                                                                             
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
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/BeamInfo.h"
// #include "larcore/Utilities/AssociationUtil.h"                                                                                                                             
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
//#include "lardata/RecoAlg/TrackMomentumCalculator.h"                                                                                                                            
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"

//Root Includes
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"

//C++ Includes 
#include <vector>
#include <iostream>
#include <string>

#include "larcorealg/Geometry/ChannelMapAlg.h"

namespace ana {
  class EvsQ;
}

class ana::EvsQ : public art::EDAnalyzer {
public:

  EvsQ(const fhicl::ParameterSet& pset);

  void analyze(const art::Event& evt);
  void endJob();
  void beginJob();
private:


  std::string fHitsModuleLabel;
  std::string fTrackModuleLabel;
  std::string fLArGeantModuleLabel;
//std::string fSimulationProducerLabel;
  //Define Variables Heere 

  //  gInterpreter->GenerateDictionary("std::map<geo::PlaneID,TGraph*>","TGraph.h;map");
  std::map<geo::PlaneID,TGraph*> ChangeVsEnergyGraph_map;
	art::ServiceHandle<geo::Geometry> geom;

};

ana::EvsQ::EvsQ(const fhicl::ParameterSet& pset) : EDAnalyzer(pset){
  
  fHitsModuleLabel = pset.get<std::string>("HitsModuleLabel");
  fTrackModuleLabel = pset.get<std::string>("TrackModuleLabel");
  fLArGeantModuleLabel =  pset.get<std::string>("LArGeantModuleLabel");
  //fSimulationProducerLabel = pset.get<std::string>("SimulationProducerLabel");

  TFile output_file("dEdX.root","RECREATE");

}



void ana::EvsQ::beginJob() {

  //Get the geometry                                                                                                                                                            
  for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){ChangeVsEnergyGraph_map[plane_id] = new TGraph(0); std::cout << plane_id << std::endl;}

}


void ana::EvsQ::analyze(const art::Event& evt) {

  //Service handles 
  art::ServiceHandle<cheat::BackTrackerService> backtracker;
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

  // Getting the Track Information                                                                                                                                          
  art::Handle< std::vector<recob::Track> > trackListHandle; 
  std::vector<art::Ptr<recob::Track> > tracklist;                                                      

  // Filling the tracklist from the tracklistHandle                                                                                                                  
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    {art::fill_ptr_vector(tracklist, trackListHandle);}


  // Getting the Hit Information                                                                                                                                           
  art::Handle< std::vector<recob::Hit> > hitListHandle;                                                           
  std::vector<art::Ptr<recob::Hit> > hits;                                                                                



	// Filling the hitlist from the hitlistHandle                                                                                                                           
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    {art::fill_ptr_vector(hits, hitListHandle);}

  //Association between Tracks and 2d Hits                                                                                                                             
  art::FindManyP<recob::Track>       fmtk(hitListHandle,   evt, fTrackModuleLabel);

  std::map<geo::PlaneID,double> charge_sum_map;
  double charge;
  double charge_deposited;
	geo::PlaneID PlaneID;
	geo::WireID wireid;
  //Loop over the hits and the charge (in ADC) found at the wires  
  for(std::vector<art::Ptr<recob::Hit> >::iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt){

    art::Ptr<recob::Hit> hit = *hitIt;

    //Calculate the charge deposited on the wire then extrapolate back to the the charge deposited in the tpc
    //charge = hit->Integral();
    charge = hit->Integral();
    charge_deposited = charge*TMath::Exp((detprop->SamplingRate()*(hit->PeakTime())) / (detprop->ElectronLifetime()*1e3));
	
    //Find the plane id
    wireid = hit->WireID();
    //std::cout << "wire id = " << wireid.Plane << std::endl;
	  PlaneID = wireid.planeID();



   charge_sum_map[PlaneID] += charge_deposited; 

    //Split the Hit into its IDE for each geant track it associates with.                                                  
    std::vector<sim::TrackIDE> trackIDs = backtracker->HitToTrackIDEs(hit);

    ////Loop over the geant track id's to the hit and find the energy that trakc deposited.   
    //    for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
    //      G4Tracks_IDE[trackIDs.at(idIt).trackID] = trackIDs.at(idIt).energy;
    // std::cout << "Track ID: " << trackIDs.at(idIt).trackID << " Energy: " << trackIDs.at(idIt).energy << std::endl;
    //}
 }


	art::Handle< std::vector<sim::SimChannel> > simchArray;
  evt.getByLabel(fLArGeantModuleLabel,simchArray);
  if(!simchArray.isValid())
    throw cet::exception(__PRETTY_FUNCTION__)
      << "Did not find sim::SimChannel with a label: " << fLArGeantModuleLabel.c_str() << std::endl;

	// total eneegy for each plane
  double total_energy0 = 0, total_energy1 = 0, total_energy2 = 0;

  // Loop over SimChannel 
	for(size_t simch_index=0; simch_index<simchArray->size(); ++simch_index) {
    const art::Ptr<sim::SimChannel> simch_ptr(simchArray,simch_index);
    
		// Each channel has a map inside it that connects a time slice to energy deposits in the detector.
    auto tdc_ide_map = simch_ptr->TDCIDEMap();

		// get id associated with channel
		auto const channel = simch_ptr->Channel();

		// map channel id to wire in tpc
		auto const* geom2 = lar::providerFrom<geo::Geometry>();
		auto wireIDs = geom2->ChannelToWire(channel);
		auto const& wireID = wireIDs.at(0);

		// get the deposited energy for each of the planes
		if(wireID.Plane == 2){
			for(auto const& tdc_ide_pair : tdc_ide_map) {

				//auto const& tdc   = tdc_ide_pair.first;
      	auto const& ide_v = tdc_ide_pair.second;

      	for(auto const& ide : ide_v) {
				total_energy2 +=  ide.energy;
      	}
			}
		}

		if(wireID.Plane == 1){
			for(auto const& tdc_ide_pair : tdc_ide_map) {

				//auto const& tdc   = tdc_ide_pair.first;
      	auto const& ide_v = tdc_ide_pair.second;

      	for(auto const& ide : ide_v) {
				total_energy1 +=  ide.energy;
      	}
			}
		}
	 
		if(wireID.Plane == 0){
			for(auto const& tdc_ide_pair : tdc_ide_map) {

				//auto const& tdc   = tdc_ide_pair.first;
      	auto const& ide_v = tdc_ide_pair.second;

      	for(auto const& ide : ide_v) {
				total_energy0 +=  ide.energy;
      	}
			}
		}
	 
}	 

	std::cout << "total energy 2 = " << total_energy2 << std::endl;

	//Fill the plane graphs with the Charge vs Energy info.
  for(std::map<geo::PlaneID,double>::iterator plane_iter=charge_sum_map.begin(); plane_iter!=charge_sum_map.end(); ++plane_iter){
		//std::cout << "plane_iter =" << plane_iter->first.Plane << std::endl;
		if(plane_iter->first.Plane == 0){
			ChangeVsEnergyGraph_map[plane_iter->first]->SetPoint(ChangeVsEnergyGraph_map[plane_iter->first]->GetN(),plane_iter->second,total_energy0);
		}

		if(plane_iter->first.Plane == 1){
			ChangeVsEnergyGraph_map[plane_iter->first]->SetPoint(ChangeVsEnergyGraph_map[plane_iter->first]->GetN(),plane_iter->second,total_energy1);
		}

		if(plane_iter->first.Plane == 2){
			ChangeVsEnergyGraph_map[plane_iter->first]->SetPoint(ChangeVsEnergyGraph_map[plane_iter->first]->GetN(),plane_iter->second,total_energy2);
		}

	}

	return; 
}



void ana::EvsQ::endJob() {

  TFile output_file("energyVschargeGraphs.root","RECREATE");
  

  //Fit a line to the energy depoisted on each plane                                                                                                                              
  std::cout << "****************************************************************" << std::endl;
  for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){
    
    if(ChangeVsEnergyGraph_map[plane_id]->GetN() < 3){continue;}
  
    TF1 *QvsE_fit = new TF1("QvsE_fit","pol1");
    ChangeVsEnergyGraph_map[plane_id]->Fit(QvsE_fit, "rob=0.75Q");

    //Get the paramters we want for the shower algorithm.                                                                                                                         
    double coefficient = QvsE_fit->GetParameter(0);
    double gradient    = QvsE_fit->GetParameter(1);

    std::cout << "** Plane:" << plane_id << " ** coefficient: " << coefficient << " ** gradient: " << gradient << " ** " << std::endl;


    //Create and Draw the canvas graphs                                                                                                                                          
    std::string title_string;
    std::stringstream sstm_QvsE_Canvas;
    sstm_QvsE_Canvas << "Charge Deposited Vs Energy Despoisted Plane: " << plane_id;
    title_string = sstm_QvsE_Canvas.str();
    const char* graph_title = title_string.c_str();
    TCanvas *c = new TCanvas(graph_title,graph_title,900,600);
    c->cd();
    ChangeVsEnergyGraph_map[plane_id]->Draw("AP");
    ChangeVsEnergyGraph_map[plane_id]->SetMarkerStyle(8);
    ChangeVsEnergyGraph_map[plane_id]->SetTitle(graph_title);
    ChangeVsEnergyGraph_map[plane_id]->GetXaxis()->SetTitle("Deposited Charge");
    ChangeVsEnergyGraph_map[plane_id]->GetYaxis()->SetTitle("Deposited Energy");

    c->Write();

    //Delete canvas' and fit                                                                                                                                                     
    delete c;
    QvsE_fit->Delete();
  }
  std::cout << "****************************************************************" << std::endl;
  std::cout << "********** These result in an Shower Energy of MeV *************" << std::endl;
  std::cout << "****************************************************************" << std::endl;

}




DEFINE_ART_MODULE(ana::EvsQ)
