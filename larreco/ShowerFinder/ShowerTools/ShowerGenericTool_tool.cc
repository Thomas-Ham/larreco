//############################################################################
//### Name:        ShowerGenericTool                                       ###
//### Author:      You                                                     ###
//### Date:        13.05.19                                                ###
//### Description: Generic form of the shower tools                        ###
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
#include "lardataobj/RecoBase/PFParticle.h"

namespace ShowerRecoTools {


  class ShowerGenericTool: public IShowerTool {

    public:

      ShowerGenericTool(const fhicl::ParameterSet& pset);

      //Generic Direction Finder
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
          art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder
          ) override;

    private:

      //Function to add the assoctions
      int AddAssociations(const art::Ptr<recob::PFParticle>& pfpPtr, art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder) override;
 
      std::string fShowerRecoEnergy;
      std::string fBestPlaneOutputLabel;   
  };


  ShowerGenericTool::ShowerGenericTool(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fShowerRecoEnergy(pset.get<std::string>("ShowerEnergyfromNumElectrons")),
    fBestPlaneOutputLabel(pset.get<std::string>("BestPlane")) 
  {
  }

  int ShowerGenericTool::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event, reco::shower::ShowerElementHolder& ShowerEleHolder){
    std::cout << "generic tool" << std::endl;
    std::vector<double> RecoEnergy;
      ShowerEleHolder.PrintElements(); 
      ShowerEleHolder.CheckElement(fShowerRecoEnergy);         
      ShowerEleHolder.GetElement(fShowerRecoEnergy, RecoEnergy);

    for(unsigned int i = 0; i < RecoEnergy.size(); i++){
        std::cout << RecoEnergy[i] << std::endl;
     }

      return 0;
  }

  int ShowerGenericTool::AddAssociations(const art::Ptr<recob::PFParticle>& pfpPtr, art::Event& Event,
      reco::shower::ShowerElementHolder& ShowerEleHolder
      ){
    return 0;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerGenericTool)
