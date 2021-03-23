#ifndef HELPER_H
#include "../interface/helper.h"
#endif

#ifndef InstantonNtuplizer_H
#define InstantonNtuplizer_H

#include <TRandom3.h>

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
                       std::map< std::string, bool >& runFlags,
                       std::map< std::string, double >& runValues,
                       std::map< std::string, std::string >& runStrings,
                       Json::Value& instan_lumi,
                       NtupleBranches* nBranches );
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
   edm::Handle<pat::MuonCollection>                            muons_                  ;
   edm::Handle< reco::VertexCollection >  vertices_;
   edm::Handle< reco::BeamSpot >  beamspot_;
   edm::Handle< std::vector<pat::PackedCandidate> > packedpfcandidates_   ;
   edm::Handle< edm::TriggerResults>                         HLTtriggers_;
   edm::Handle<pat::TriggerObjectStandAloneCollection>       triggerObjects;
   edm::Handle< reco::GenParticleCollection >  genParticles_;
   edm::Handle< std::vector<pat::PackedGenParticle> >  packedgenParticles_;
   edm::Handle<pat::JetCollection>              jets_           ;
   edm::ESHandle<TransientTrackBuilder> builder;
   const MagneticField                 *fMagneticField;


   bool runOnMC_;
   bool verbose_;


   bool flag_fill = false;
   helper aux;
   float chi = 0.;
   float ndf = 0.;
   Json::Value instan_lumi_;
};
#endif
