//############################################################################
//### Name:        ShowerTrackDirection                                    ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the shower direction using the         ###
//###              initial track method.                                   ###
//############################################################################

#include "larreco/ShowerFinder/ShowerTools/IShowerTool.h"

//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft Includes 
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larreco/RecoAlg/SBNShowerAlg.h"

//C++ Includes 
#include <iostream>
#include <cmath>

//Root Includes 
#include "TVector3.h"
#include "TMath.h"
#include "TH1.h" 

#include "TFile.h"

namespace ShowerRecoTools {


	class ShowerTrackDirection:IShowerTool {

		public:

			ShowerTrackDirection(const fhicl::ParameterSet& pset);

			~ShowerTrackDirection(); 

			//Generic Direction Finder
			int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
					art::Event& Event,
					reco::shower::ShowerElementHolder& ShowerEleHolder
					) override;

		private:

			// Define standard art tool interface
			void configure(const fhicl::ParameterSet& pset) override;

			//Algoritma
			shower::SBNShowerAlg       fSBNShowerAlg;

			//fcl
			bool fUsePandoraVertex;
			bool fUseStartPosition;
			bool fDebugEVD; 

			//services
			detinfo::DetectorProperties const* fDetProp;

	};


	ShowerTrackDirection::ShowerTrackDirection(const fhicl::ParameterSet& pset)
		:fSBNShowerAlg(pset.get<fhicl::ParameterSet>("SBNShowerAlg")),
		fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>())

	{
		configure(pset);
	}

	ShowerTrackDirection::~ShowerTrackDirection()
	{
	}

	void ShowerTrackDirection::configure(const fhicl::ParameterSet& pset)
	{
		fUsePandoraVertex         = pset.get<bool>("UsePandoraVertex");
		fUseStartPosition         = pset.get<bool>("UseStartPosition");
		fDebugEVD                 = pset.get<bool>("DebugEVD");  

	}

	int ShowerTrackDirection::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
			art::Event& Event,
			reco::shower::ShowerElementHolder& ShowerEleHolder){

		std::cout << "ShowerTrackDirection_tool" << std::endl;

		//Check the Track has been defined
		if(!ShowerEleHolder.CheckElement("InitialTrack")){
			mf::LogError("ShowerTrackDirection") << "Initial track not set, returning"<< std::endl;
			return 1;
		}

		//Check the start position is set.
		if(!ShowerEleHolder.CheckElement("ShowerStartPosition") && fUsePandoraVertex){
			mf::LogError("ShowerTrackDirection") << "Start position not set, returning "<< std::endl;
			return 1;
		}

		recob::Track InitialTrack;
		ShowerEleHolder.GetElement("InitialTrack",InitialTrack);

		TH1D* XHist = new TH1D("XHist","XHist",30,-1,1);
		TH1D* YHist = new TH1D("YHist","YHist",30,-1,1);
		TH1D* ZHist = new TH1D("ZHist","ZHist",30,-1,1);

		if(fUseStartPosition){ 
			geo::Point_t StartPosition;
			if(fUsePandoraVertex){
				TVector3 StartPosition_vec = {-999, -999, -999};
				ShowerEleHolder.GetElement("ShowerStartPosition",StartPosition_vec);
				StartPosition.SetCoordinates(StartPosition_vec.X(),StartPosition_vec.Y(),StartPosition_vec.Z()); 
			}
			else{ 
				StartPosition = InitialTrack.Start();
			}

			//Loop over the trajectory points and average the direction
			for(unsigned int traj=0; traj<InitialTrack.NumberTrajectoryPoints(); ++traj){
				geo::Point_t  TrajPosition = InitialTrack.LocationAtPoint(traj);
				//std::cout << "Traj Position X: " << TrajPosition.X() << " Y: " << TrajPosition.Y() << " Z: " <<  TrajPosition.Z() << std::endl; 
				//std::cout << "Start Position X: " <<  StartPosition.X() <<" Y: " << StartPosition.Y() << " Z: " << StartPosition.Z() << std::endl;
				geo::Vector_t Direction    = (TrajPosition - StartPosition).Unit();

				XHist->Fill(Direction.X());
				YHist->Fill(Direction.Y());
				ZHist->Fill(Direction.Z());
			}

			TFile *MyFile = new TFile("test.root","RECREATE");
			XHist->Write();
			YHist->Write();
			ZHist->Write();
			MyFile->Close();

			//geo::Vector_t Mean = {XHist->GetMean(), YHist->GetMean(), ZHist->GetMean()};
			//geo::Vector_t RMS  = {XHist->GetRMS(), YHist->GetRMS(), ZHist->GetRMS()};

			float N = 0.;
			TVector3 Direction_Mean = {0,0,0};
			//Remove trajectory points from the mean that are not with one sigma.
			for(unsigned int traj=0; traj< InitialTrack.NumberTrajectoryPoints(); ++traj){
				geo::Point_t TrajPosition = InitialTrack.LocationAtPoint(traj);
				geo::Vector_t Direction   = (TrajPosition - StartPosition).Unit();
				// if((TMath::Abs((Direction-Mean).X()) < 2*RMS.X()) && 
				//    (TMath::Abs((Direction-Mean).Y()) < 2*RMS.Y()) && 
				//    (TMath::Abs((Direction-Mean).Z()) < 2*RMS.Z())){
				TVector3 Direction_vec = {Direction.X(),Direction.Y(),Direction.Z()};
				Direction_Mean += Direction_vec;
				++N;
				//}
			}


			//Take the mean value
			if(N > 0){
				Direction_Mean *= 1/N;
				TVector3 Direction = Direction_Mean.Unit();
				TVector3 DirectionErr = {-999,-999,-999};
				ShowerEleHolder.SetElement(Direction,DirectionErr,"ShowerDirection");
			}
			else{
				mf::LogError("ShowerDirection") << "None of the points are within 2 sigma"<< std::endl;
				return 1;
			}

			if (fDebugEVD){
				std::cout << "running mad" << std::endl;
				fSBNShowerAlg.DebugEVD(pfparticle,Event,ShowerEleHolder);
			}
			return 0;

		}	else { //If not using start position
			//Loop over the trajectory points and average the direction
			for(unsigned int traj=0; traj<InitialTrack.NumberTrajectoryPoints(); ++traj){
				geo::Vector_t  Direction = InitialTrack.DirectionAtPoint(traj);
				XHist->Fill(Direction.X());
				YHist->Fill(Direction.Y());
				ZHist->Fill(Direction.Z());
			}

			geo::Vector_t Mean = {XHist->GetMean(), YHist->GetMean(), ZHist->GetMean()};
			geo::Vector_t RMS  = {XHist->GetRMS(), YHist->GetRMS(), ZHist->GetRMS()};

			//Remove trajectory points from the mean that are not with one sigma.
			float N = 0.;
			TVector3 Direction_Mean = {0,0,0};
			for(unsigned int traj=0; traj<InitialTrack.NumberTrajectoryPoints(); ++traj){
				geo::Vector_t Direction = InitialTrack.DirectionAtPoint(traj);
				if((TMath::Abs((Direction-Mean).X()) < 2*RMS.X()) && 
						(TMath::Abs((Direction-Mean).Y())< 2*RMS.Y()) && 
						(TMath::Abs((Direction-Mean).Z()) < 2*RMS.Z())){
					TVector3 Direction_vec = {Direction.X(),Direction.Y(),Direction.Z()};
					Direction_Mean += Direction_vec;
					++N;
				}
			}

			//Take the mean value
			if(N>0){
				Direction_Mean *= 1/N;
				TVector3 Direction = Direction_Mean.Unit();
				TVector3 DirectionErr = {-999,-999,-999};
				ShowerEleHolder.SetElement(Direction,DirectionErr,"ShowerDirection");
			}
			else{
				mf::LogError("ShowerDirection") << "None of the points are within 2 sigma"<< std::endl;
				return 1;
			}

			if (fDebugEVD){
				std::cout << "running mad" << std::endl;
				fSBNShowerAlg.DebugEVD(pfparticle,Event,ShowerEleHolder);
			}

			return 0;
		}
	}
}


DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerTrackDirection)

