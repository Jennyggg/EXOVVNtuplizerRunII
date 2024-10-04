#ifndef HELPER_H
#include "../interface/helper.h"
#include "fastjet/contrib/Njettiness.hh"
#endif
#include "../interface/PPSTrackingEfficiency.h"
#include <TRandom3.h>
#ifndef InstantonNtuplizer_H
#define InstantonNtuplizer_H
using namespace std;
class InstantonNtuplizer : public CandidateNtuplizer {

public:
  InstantonNtuplizer(  edm::EDGetTokenT<pat::MuonCollection>    muonToken   ,
                       edm::EDGetTokenT<reco::VertexCollection> verticeToken,
                       edm::EDGetTokenT<reco::BeamSpot>             beamToken ,
                       edm::EDGetTokenT<pat::PackedCandidateCollection> packedpfcandidatesToken,
                       edm::EDGetTokenT<edm::TriggerResults> triggertoken,
                       edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerobject,
                       edm::EDGetTokenT<reco::GenParticleCollection> genptoken,
                       edm::EDGetTokenT<pat::PackedGenParticleCollection> packedgenptoken,
                       edm::EDGetTokenT<pat::JetCollection>     jettoken,
                       edm::EDGetTokenT<std::vector<reco::VertexCompositePtrCandidate> > svToken,
                       edm::EDGetTokenT<edm::SortedCollection<CaloTower>> CaloTowerCollection ,
                       edm::EDGetTokenT<reco::ForwardProtonCollection> ctppsProton_single_rp,
                       edm::EDGetTokenT<reco::ForwardProtonCollection> ctppsProton_multi_rp,
                       std::map< std::string, bool >& runFlags,
                       std::map< std::string, double >& runValues,
                       std::map< std::string, std::string >& runStrings,
                       Json::Value& instan_lumi,
                       NtupleBranches* nBranches);
  ~InstantonNtuplizer(void);
  bool fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );


private:
   edm::EDGetTokenT<pat::MuonCollection>    muonToken_   ;
   edm::EDGetTokenT<reco::VertexCollection> verticeToken_   ;
   edm::EDGetTokenT<reco::BeamSpot> bsToken_   ;
   edm::EDGetTokenT<pat::PackedCandidateCollection>             packedpfcandidatesToken_;
   edm::EDGetTokenT<edm::TriggerResults>                     HLTtriggersToken_;
   edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection>  triggerObjects_;
   edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
   edm::EDGetTokenT<pat::PackedGenParticleCollection> packedgenParticlesToken_;
   edm::EDGetTokenT<pat::JetCollection>     jetInputToken_    ;
   edm::EDGetTokenT<std::vector<reco::VertexCompositePtrCandidate> > svToken_ ;
   edm::Handle<pat::MuonCollection>                            muons_;
   edm::Handle< reco::VertexCollection >  vertices_;
   edm::Handle< reco::BeamSpot >  beamspot_;
   edm::Handle< std::vector<pat::PackedCandidate> > packedpfcandidates_   ;
   edm::Handle< edm::TriggerResults>                         HLTtriggers_;
   edm::Handle<pat::TriggerObjectStandAloneCollection>       triggerObjects;
   edm::Handle< reco::GenParticleCollection >  genParticles_;
   edm::Handle< std::vector<pat::PackedGenParticle> >  packedgenParticles_;
   edm::Handle<pat::JetCollection>              jets_           ;
   edm::Handle<std::vector<reco::VertexCompositePtrCandidate>> svs_;
   edm::Handle<std::vector<reco::ForwardProton>> ctppsProton_single_rp_;
   edm::Handle<std::vector<reco::ForwardProton>> ctppsProton_multi_rp_;
   edm::EDGetTokenT<edm::SortedCollection<CaloTower>> CaloTowerCollection_;
   edm::EDGetTokenT<reco::ForwardProtonCollection> ctppsProton_single_rpToken_;
   edm::EDGetTokenT<reco::ForwardProtonCollection> ctppsProton_multi_rpToken_; 
   edm::ESHandle<TransientTrackBuilder> builder;
   const MagneticField                 *fMagneticField;


   bool runOnMC_;
   bool runOnMCPU_;
   bool runOnZeroBias_;
   bool runOnHerwigInstanton_;
   bool runOnSherpaInstanton_;
   bool runOnMCTest_;
   bool verbose_;
   bool doTrack_;
   bool doTrackJet_;
   bool doForwardProton_;
   bool trackrandomdrop_;


   bool flag_fill = false;
   helper aux;
   float chi = 0.;
   float ndf = 0.;
   Json::Value instan_lumi_;
   std::map< std::string, std::string > runStrings_; 
   std::vector<std::string> detNames = { "EB", "EE", "HB", "HE", "HFm", "HFp", "unknown" };
   std::vector<std::pair<double,double>> detLimits = {
     {0    , 1.479 }, // EB (Exclude transition region between calo barrel and endcap)
     {1.479, 3.0   }, // EE (Exclude transition region between calo barrel and endcap)
     {0    , 1.3   }, // HB
     {1.3  , 3.0   }, // HE
     {3.0  , 5.2   }, // HFm    // sleontsi
     {3.0  , 5.2   }, // HFp*/
   };
   enum CaloType { kEB, kEE, kHB, kHE, kHFp, kHFm, nCaloTypes };
   CaloType GetTowerSubdetHad(double&);
   CaloType GetTowerSubdetEm(double&);
   const std::map<CaloType, std::string> caloName = {
     { kEB  , "EB"  },
     { kEE  , "EE"  },
     { kHB  , "HB"  },
     { kHE  , "HE"  },
     { kHFp , "HFp" },
     { kHFm , "HFm" },
   };
   const std::map<CaloType, double> noiseThreshold = {
     { kEB  , 0.7 },
     { kEE  , 7.5 },
     { kHB  , 2.8  },
     { kHE  , 2.4  },
     { kHFp , 7.2  },
     { kHFm , 7.5  },
   };
   const double maxEtaEB = detLimits.at(0).second;
   const double minEtaEE = detLimits.at(1).first;
   const double maxEtaEE = detLimits.at(1).second;
 
   const double maxEtaHB = detLimits.at(2).second;
   const double minEtaHE = detLimits.at(3).first;
   const double maxEtaHE = detLimits.at(3).second;
 
   const double minEtaHF = abs(detLimits.at(4).first);
   const double maxEtaHF = abs(detLimits.at(4).second);

   PPSTrackingEfficiency* forwardProtonEff_;
};
#endif
