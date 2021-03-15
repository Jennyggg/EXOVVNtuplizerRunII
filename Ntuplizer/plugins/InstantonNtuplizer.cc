#include "../interface/InstantonNtuplizer.h"
#include "TMatrixTBase.h"
#include "TMatrixDSymEigen.h"
#include "TVector3.h"
#include "../interface/helper.h"

void calcS_T_B(std::vector<pat::PackedCandidate>& momenta, Float_t& Spherocity,Float_t& Thrust,Float_t& Broaden, TVector3& taxis, int nSeed){
  if (momenta.size()==0) {Spherocity=-1;Thrust=-1;Broaden=-1;return;}
  vector<TLorentzVector> p;
  TLorentzVector p_total;
  for(size_t trkindex=0; trkindex<momenta.size(); trkindex++){
    TLorentzVector momenta_p4(momenta[trkindex].px(),momenta[trkindex].py(),momenta[trkindex].pz(),momenta[trkindex].energy());
    p.push_back(momenta_p4);
    p_total+=momenta_p4;
  }
  TVector3 v_boost = -p_total.BoostVector();
  float p2_sum=0;
  float p_sum=0;
  TMatrixDSym Sph_matrix(3);
//Calculate the Spherocity in the COM frame
  for (size_t j=0; j< momenta.size(); j++) {
    p[j].Boost(v_boost);
    p2_sum += p[j].P()*p[j].P();
    p_sum += p[j].P();
    Sph_matrix(0,0)+=(p[j][0])*(p[j][0]);
    Sph_matrix(1,1)+=(p[j][1])*(p[j][1]);
    Sph_matrix(2,2)+=(p[j][2])*(p[j][2]);
    Sph_matrix(0,1)+=(p[j][0])*(p[j][1]);
    Sph_matrix(0,2)+=(p[j][0])*(p[j][2]);
    Sph_matrix(1,2)+=(p[j][1])*(p[j][2]);
        
  }
  Sph_matrix(1,0) = Sph_matrix(0,1);
  Sph_matrix(2,0) = Sph_matrix(0,2);
  Sph_matrix(2,1) = Sph_matrix(1,2);
  TMatrixDSymEigen eigen(Sph_matrix);
  TVectorD eigenVal=eigen.GetEigenValues();
  Spherocity=1.5*(eigenVal.Sum()-eigenVal.Max())/p2_sum;

//Calculate the Thrust

  std::sort(p.begin(),p.end(),[](TLorentzVector a, TLorentzVector b){
    return a.P()>b.P();
  });
  Thrust = 0;
  Float_t Thr;
  vector<Float_t> Thr_val;
  vector<TVector3> Thr_vec;
  for (unsigned int j=0; j < pow(2,nSeed); j++){
    int sign=j;
    TVector3 foo(0,0,0);
    for (int k = 0 ; k < nSeed ; k++) {
      (sign % 2) == 1 ? foo += p[k].Vect() : foo -= p[k].Vect();
      sign /= 2;
    }
    foo=foo.Unit();
    double diff=999.;
    while (diff > 1e-5) {
      TVector3 foobar(0,0,0);
      for (size_t k=0 ; k<momenta.size() ; k++)
        foo.Dot(p[k].Vect())>0 ? foobar+=p[k].Vect() : foobar-=p[k].Vect();
      diff = (foo-foobar.Unit()).Mag();
      foo=foobar.Unit();
    }
    Thr=0.;
    for (size_t k=0 ; k<momenta.size() ; k++)
      Thr+=abs(foo.Dot(p[k].Vect()));
    Thr_val.push_back(Thr);
    Thr_vec.push_back(foo);
  }
  Thrust=0.;
  for (unsigned int j=0; j < Thr_val.size(); j++){
    if (Thr_val[j]>Thrust){
      Thrust=Thr_val[j];
      taxis=Thr_vec[j];
    }
  }
  Thrust=1-Thrust/p_sum;
  Broaden=0.;
//Calculate the Broaden
  for(size_t j=0; j< momenta.size(); j++)
    Broaden+=abs(taxis.Cross(p[j].Vect()).Mag());
  Broaden/=p_sum*2;
}

//===================================================================================================================
InstantonNtuplizer::InstantonNtuplizer(  edm::EDGetTokenT<pat::MuonCollection>    muonToken   ,
                                         edm::EDGetTokenT<reco::VertexCollection> verticeToken,
                                         edm::EDGetTokenT<reco::BeamSpot>             beamToken,
                                         edm::EDGetTokenT<pat::PackedCandidateCollection> packedpfcandidatesToken,
                                         edm::EDGetTokenT<edm::TriggerResults> triggertoken,
                                         edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerobject,
                                         edm::EDGetTokenT<reco::GenParticleCollection> genptoken,
                                         edm::EDGetTokenT<pat::PackedGenParticleCollection> packedgenptoken,
                                         edm::EDGetTokenT<pat::JetCollection>     jettoken,
                                         std::map< std::string, bool >& runFlags,
                                         std::map< std::string, double >& runValues,
                                         std::map< std::string, std::string >& runStrings,
                                         NtupleBranches* nBranches )
: CandidateNtuplizer ( nBranches )
  , muonToken_          ( muonToken )
  , verticeToken_          ( verticeToken )
  , bsToken_ (beamToken)
  , packedpfcandidatesToken_(packedpfcandidatesToken)
  , HLTtriggersToken_   ( triggertoken )
  , triggerObjects_     ( triggerobject )
  , genParticlesToken_( genptoken )
  , packedgenParticlesToken_( packedgenptoken )
  , jetInputToken_     ( jettoken     )
  , runOnMC_   (runFlags["runOnMC"])
  , verbose_   (runFlags["verbose"])
{
  if(verbose_){
    std::cout << "[InstantonNtuplizer] runOnMC    = " << runOnMC_ << std::endl;
  }
}
InstantonNtuplizer::~InstantonNtuplizer(void){

}
using namespace fastjet;
bool InstantonNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){
  if(verbose_) std::cout << "[InstantonNtuplizer] ---------------- event, run, lumi = " << event.id().event() << " " << event.id().run() << " " << event.id().luminosityBlock() << "----------------" << std::endl;
//TrackJet calculating parameters
  double R=0.4;
  JetDefinition jet_def(antikt_algorithm, R);
  float TrackJet_PtCut=3.0; //Pt Cut > 3.0 GeV
  int nSeed_ThrustCalculation=3;

//Get the Tokens from MINIAOD
  if(runOnMC_){
    event.getByToken(genParticlesToken_ , genParticles_);
    event.getByToken(packedgenParticlesToken_ , packedgenParticles_);//packed particles all have status 1 (final state)
  }

  event.getByToken(verticeToken_   , vertices_     );
  event.getByToken(bsToken_, beamspot_);
  event.getByToken(muonToken_   , muons_    );
  event.getByToken(triggerObjects_  , triggerObjects);
  event.getByToken(HLTtriggersToken_,HLTtriggers_);
  std::vector<pat::Muon> muoncollection;
  muoncollection.clear();

  event.getByToken( packedpfcandidatesToken_               , packedpfcandidates_      );
  event.getByToken(jetInputToken_      , jets_    );
  std::vector<pfcand_struct> pfcands;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
//  std::vector<int> alltracks_idx;
// Count number of tracks and displaced tracks (transverse impact paramter > 0.02 cm)
  int N_Trk_total=0;
  int N_Trk_Displaced_total=0;
  int N_Trk_goodDisplaced_total=0;
  int PV_N_good_hasgoodMuon = 0;
  int PV_N_good_hasgoodPFMuon = 0;
  int N_Jet_total=jets_->size();
  std::vector<int> N_Trk(vertices_->size(),0);
  std::vector<int> N_Trk_highPurity(vertices_->size(),0);
  std::vector<int> N_Trk_highPurity_pt05(vertices_->size(),0);
  std::vector<int> N_Trk_highPurity_pt08(vertices_->size(),0);
  std::vector<int> N_Trk_highPurity_pt1(vertices_->size(),0);
  std::vector<int> N_Trk_PVAssociationQuality0(vertices_->size(),0);
  std::vector<int> N_Trk_PVAssociationQuality1(vertices_->size(),0);
  std::vector<int> N_Trk_PVAssociationQuality4(vertices_->size(),0);
  std::vector<int> N_Trk_PVAssociationQuality5(vertices_->size(),0);
  std::vector<int> N_Trk_PVAssociationQuality6(vertices_->size(),0);
  std::vector<int> N_Trk_PVAssociationQuality7(vertices_->size(),0);
  std::vector<int> N_Trk_PVAssociationQualityLeq4_highPurity(vertices_->size(),0);
  std::vector<int> N_Trk_PVAssociationQualityLeq4_highPurity_pt05(vertices_->size(),0);
  std::vector<int> N_Trk_PVAssociationQualityLeq4_highPurity_pt08(vertices_->size(),0);
  std::vector<int> N_Trk_PVAssociationQualityLeq4_highPurity_pt1(vertices_->size(),0);
  std::vector<int> N_Trk_Displaced(vertices_->size(),0);
  std::vector<int> N_Trk_Displaced_highPurity(vertices_->size(),0);
  std::vector<int> N_Trk_Displaced_highPurity_pt05(vertices_->size(),0);
  std::vector<int> N_Trk_Displaced_highPurity_pt08(vertices_->size(),0);
  std::vector<int> N_Trk_Displaced_highPurity_pt1(vertices_->size(),0);
  std::vector<int> N_Trk_Displaced_PVAssociationQuality0(vertices_->size(),0);
  std::vector<int> N_Trk_Displaced_PVAssociationQuality1(vertices_->size(),0);
  std::vector<int> N_Trk_Displaced_PVAssociationQuality4(vertices_->size(),0);
  std::vector<int> N_Trk_Displaced_PVAssociationQuality5(vertices_->size(),0);
  std::vector<int> N_Trk_Displaced_PVAssociationQuality6(vertices_->size(),0);
  std::vector<int> N_Trk_Displaced_PVAssociationQuality7(vertices_->size(),0);
  std::vector<int> N_Trk_Displaced_PVAssociationQualityLeq4_highPurity(vertices_->size(),0);
  std::vector<int> N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt05(vertices_->size(),0);
  std::vector<int> N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt08(vertices_->size(),0);
  std::vector<int> N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt1(vertices_->size(),0);
  std::vector<int> N_Trk_goodDisplaced(vertices_->size(),0);
  std::vector<int> N_Trk_goodDisplaced_highPurity(vertices_->size(),0);
  std::vector<int> N_Trk_goodDisplaced_highPurity_pt05(vertices_->size(),0);
  std::vector<int> N_Trk_goodDisplaced_highPurity_pt08(vertices_->size(),0);
  std::vector<int> N_Trk_goodDisplaced_highPurity_pt1(vertices_->size(),0);
  std::vector<int> N_Trk_goodDisplaced_PVAssociationQuality0(vertices_->size(),0);
  std::vector<int> N_Trk_goodDisplaced_PVAssociationQuality1(vertices_->size(),0);
  std::vector<int> N_Trk_goodDisplaced_PVAssociationQuality4(vertices_->size(),0);
  std::vector<int> N_Trk_goodDisplaced_PVAssociationQuality5(vertices_->size(),0);
  std::vector<int> N_Trk_goodDisplaced_PVAssociationQuality6(vertices_->size(),0);
  std::vector<int> N_Trk_goodDisplaced_PVAssociationQuality7(vertices_->size(),0);
  std::vector<int> N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity(vertices_->size(),0);
  std::vector<int> N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt05(vertices_->size(),0);
  std::vector<int> N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt08(vertices_->size(),0);
  std::vector<int> N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt1(vertices_->size(),0);
  std::vector<float> Trk_mass(vertices_->size(),0);
  std::vector<int> N_Jet(vertices_->size(),0);
  std::vector<int> vtx_N_goodMuon(vertices_->size(),0);
  std::vector<int> vtx_N_goodPFMuon(vertices_->size(),0);
  std::vector<bool> vtx_isBPHtrigger_fromMuon(vertices_->size(),false); 
  std::vector<bool> vtx_isBPHtrigger_fromPFMuon(vertices_->size(),false);
  std::vector<float> vtx_goodMuonIP(vertices_->size(),999.);
  std::vector<float> vtx_goodPFMuon_PVDistance(vertices_->size(),999.);
  std::vector<std::vector<reco::TransientTrack>> alltracks(vertices_->size());
  std::vector<std::vector<pat::PackedCandidate>> alltracks_pf(vertices_->size());
//  for(reco::VertexCollection::const_iterator vtx = vertices_->begin(); vtx != vertices_->end(); ++vtx){
//    cout<<"PV "<<vtx->position().X()<<" "<<vtx->position().Y()<<" "<<vtx->position().Z()<<endl;
//  }
  std::vector<pat::PackedCandidate> alltracks_pfmuon;
  for( size_t ii = 0; ii < packedpfcandidates_->size(); ++ii ){
    pat::PackedCandidate pf = (*packedpfcandidates_)[ii];
      if(pf.pdgId()==13||pf.pdgId()==-13){ 
        alltracks_pfmuon.push_back(pf);
      }
    if(pf.hasTrackDetails()){
      auto it1 = pf.bestTrack();
      if(it1==nullptr){
        std::cout << "This is null!!" << std::endl;
      }
      bool isDisplaced=false;
      bool isgoodDisplaced=false;
      bool ishighPurity=false;
      bool ishighPurity_pt05=false;
      bool ishighPurity_pt08=false;
      bool ishighPurity_pt1=false;
      N_Trk_total++;
//associate the track to a primary vertex 
      int vtxindex=pf.vertexRef().key();
//      for(reco::VertexCollection::const_iterator vtx = vertices_->begin(); vtx != vertices_->end(); ++vtx, ++vtxindex){
//        if(pf.vertexRef()->z()==vtx->position().Z()){
      N_Trk[vtxindex]++;
//Refer to: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookChapter7#Selecting_Good_Quality_Tracks
      if(it1->quality(reco::Track::highPurity)){
        ishighPurity=true;
        N_Trk_highPurity[vtxindex]++;
        if(pf.pt()>0.5) {ishighPurity_pt05=true; N_Trk_highPurity_pt05[vtxindex]++;}
        if(pf.pt()>0.8) {ishighPurity_pt08=true; N_Trk_highPurity_pt08[vtxindex]++;}
        if(pf.pt()>1) {ishighPurity_pt1=true;  N_Trk_highPurity_pt1[vtxindex]++;}
      }
      reco::TransientTrack  _track_ = (*builder).build(pf.pseudoTrack());
      _track_.setBeamSpot(*beamspot_);
      GlobalPoint vert(pf.vertexRef()->x(), pf.vertexRef()->y(), pf.vertexRef()->z());
      TrajectoryStateClosestToPoint  traj = _track_.trajectoryStateClosestToPoint(vert);
      if(traj.perigeeParameters().transverseImpactParameter() > 0.02){
//          math::XYZPoint vert(pf.vertexRef()->x(), pf.vertexRef()->y(), pf.vertexRef()->z());
//          if(pf.dxy(vert)> 0.002){
        N_Trk_Displaced[vtxindex]++;
        N_Trk_Displaced_highPurity[vtxindex]+=ishighPurity;
        N_Trk_Displaced_highPurity_pt05[vtxindex]+=ishighPurity_pt05;
        N_Trk_Displaced_highPurity_pt08[vtxindex]+=ishighPurity_pt08;
        N_Trk_Displaced_highPurity_pt1[vtxindex]+=ishighPurity_pt1;
        N_Trk_Displaced_total++;
        isDisplaced=true;
        if(traj.perigeeParameters().transverseImpactParameter()/traj.perigeeError().transverseImpactParameterError()>5.){
          N_Trk_goodDisplaced[vtxindex]++;
          N_Trk_goodDisplaced_highPurity[vtxindex]+=ishighPurity;
          N_Trk_goodDisplaced_highPurity_pt05[vtxindex]+=ishighPurity_pt05;
          N_Trk_goodDisplaced_highPurity_pt08[vtxindex]+=ishighPurity_pt08;
          N_Trk_goodDisplaced_highPurity_pt1[vtxindex]+=ishighPurity_pt1;
          N_Trk_goodDisplaced_total++;
          isgoodDisplaced=true;
        }     
      }
      switch(pf.pvAssociationQuality()){
        case pat::PackedCandidate::PVAssociationQuality::NotReconstructedPrimary: N_Trk_PVAssociationQuality0[vtxindex]++;
                N_Trk_Displaced_PVAssociationQuality0[vtxindex]+=(int)isDisplaced;
                N_Trk_goodDisplaced_PVAssociationQuality0[vtxindex]+=(int)isgoodDisplaced;
                break;
        case pat::PackedCandidate::PVAssociationQuality::OtherDeltaZ: N_Trk_PVAssociationQuality1[vtxindex]++;
                N_Trk_Displaced_PVAssociationQuality1[vtxindex]+=(int)isDisplaced;
                N_Trk_goodDisplaced_PVAssociationQuality1[vtxindex]+=(int)isgoodDisplaced;
                break;
        case pat::PackedCandidate::PVAssociationQuality::CompatibilityBTag: N_Trk_PVAssociationQuality4[vtxindex]++;
                N_Trk_Displaced_PVAssociationQuality4[vtxindex]+=(int)isDisplaced;
                N_Trk_goodDisplaced_PVAssociationQuality4[vtxindex]+=(int)isgoodDisplaced;
                N_Trk_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)ishighPurity;
                N_Trk_PVAssociationQualityLeq4_highPurity_pt05[vtxindex]+=(int)ishighPurity_pt05;
                N_Trk_PVAssociationQualityLeq4_highPurity_pt08[vtxindex]+=(int)ishighPurity_pt08;
                N_Trk_PVAssociationQualityLeq4_highPurity_pt1[vtxindex]+=(int)ishighPurity_pt1;
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(ishighPurity&&isDisplaced);
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt05[vtxindex]+=(int)(ishighPurity_pt05&&isDisplaced);
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt08[vtxindex]+=(int)(ishighPurity_pt08&&isDisplaced);
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt1[vtxindex]+=(int)(ishighPurity_pt1&&isDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(ishighPurity&&isgoodDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt05[vtxindex]+=(int)(ishighPurity_pt05&&isgoodDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt08[vtxindex]+=(int)(ishighPurity_pt08&&isgoodDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt1[vtxindex]+=(int)(ishighPurity_pt1&&isgoodDisplaced);
                break;
        case pat::PackedCandidate::PVAssociationQuality::CompatibilityDz: N_Trk_PVAssociationQuality5[vtxindex]++;
                N_Trk_Displaced_PVAssociationQuality5[vtxindex]+=(int)isDisplaced;
                N_Trk_goodDisplaced_PVAssociationQuality5[vtxindex]+=(int)isgoodDisplaced;
                N_Trk_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)ishighPurity;
                N_Trk_PVAssociationQualityLeq4_highPurity_pt05[vtxindex]+=(int)ishighPurity_pt05;
                N_Trk_PVAssociationQualityLeq4_highPurity_pt08[vtxindex]+=(int)ishighPurity_pt08;
                N_Trk_PVAssociationQualityLeq4_highPurity_pt1[vtxindex]+=(int)ishighPurity_pt1;
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(ishighPurity&&isDisplaced);
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt05[vtxindex]+=(int)(ishighPurity_pt05&&isDisplaced);
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt08[vtxindex]+=(int)(ishighPurity_pt08&&isDisplaced);
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt1[vtxindex]+=(int)(ishighPurity_pt1&&isDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(ishighPurity&&isgoodDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt05[vtxindex]+=(int)(ishighPurity_pt05&&isgoodDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt08[vtxindex]+=(int)(ishighPurity_pt08&&isgoodDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt1[vtxindex]+=(int)(ishighPurity_pt1&&isgoodDisplaced);
                break;
        case pat::PackedCandidate::PVAssociationQuality::UsedInFitLoose: N_Trk_PVAssociationQuality6[vtxindex]++;
                N_Trk_Displaced_PVAssociationQuality6[vtxindex]+=(int)isDisplaced;
                N_Trk_goodDisplaced_PVAssociationQuality6[vtxindex]+=(int)isgoodDisplaced;
                N_Trk_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)ishighPurity;
                N_Trk_PVAssociationQualityLeq4_highPurity_pt05[vtxindex]+=(int)ishighPurity_pt05;
                N_Trk_PVAssociationQualityLeq4_highPurity_pt08[vtxindex]+=(int)ishighPurity_pt08;
                N_Trk_PVAssociationQualityLeq4_highPurity_pt1[vtxindex]+=(int)ishighPurity_pt1;
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(ishighPurity&&isDisplaced);
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt05[vtxindex]+=(int)(ishighPurity_pt05&&isDisplaced);
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt08[vtxindex]+=(int)(ishighPurity_pt08&&isDisplaced);
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt1[vtxindex]+=(int)(ishighPurity_pt1&&isDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(ishighPurity&&isgoodDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt05[vtxindex]+=(int)(ishighPurity_pt05&&isgoodDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt08[vtxindex]+=(int)(ishighPurity_pt08&&isgoodDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt1[vtxindex]+=(int)(ishighPurity_pt1&&isgoodDisplaced);
                break;
        case pat::PackedCandidate::PVAssociationQuality::UsedInFitTight: N_Trk_PVAssociationQuality7[vtxindex]++;
                N_Trk_Displaced_PVAssociationQuality7[vtxindex]+=(int)isDisplaced;
                N_Trk_goodDisplaced_PVAssociationQuality7[vtxindex]+=(int)isgoodDisplaced;
                N_Trk_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)ishighPurity;
                N_Trk_PVAssociationQualityLeq4_highPurity_pt05[vtxindex]+=(int)ishighPurity_pt05;
                N_Trk_PVAssociationQualityLeq4_highPurity_pt08[vtxindex]+=(int)ishighPurity_pt08;
                N_Trk_PVAssociationQualityLeq4_highPurity_pt1[vtxindex]+=(int)ishighPurity_pt1;
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(ishighPurity&&isDisplaced);
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt05[vtxindex]+=(int)(ishighPurity_pt05&&isDisplaced);
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt08[vtxindex]+=(int)(ishighPurity_pt08&&isDisplaced);
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt1[vtxindex]+=(int)(ishighPurity_pt1&&isDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(ishighPurity&&isgoodDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt05[vtxindex]+=(int)(ishighPurity_pt05&&isgoodDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt08[vtxindex]+=(int)(ishighPurity_pt08&&isgoodDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt1[vtxindex]+=(int)(ishighPurity_pt1&&isgoodDisplaced);
                break;
      }
      alltracks[vtxindex].push_back(_track_);
      alltracks_pf[vtxindex].push_back(pf);

//          break;
//        }
//      }
    }
  }

//Calculate total mass of the tracks for each PV
  for( size_t j = 0; j < vertices_->size(); j++){
    math::XYZTLorentzVector p4(0,0,0,0);
    for(int k = 0; k < N_Trk[j]; k++)
      p4 = p4 + alltracks_pf[j][k].p4();
    Trk_mass[j] = p4.M();
  }

//Calculate TrackJets
  for( size_t j = 0; j < vertices_->size(); j++){
    std::vector<PseudoJet> particles;
    std::vector<TLorentzVector> TrackJet_P4_thisvertex;
    for( size_t trkindex = 0; trkindex < alltracks_pf[j].size(); trkindex++)
      particles.push_back(PseudoJet(alltracks_pf[j][trkindex].px(),alltracks_pf[j][trkindex].py(),alltracks_pf[j][trkindex].pz(),alltracks_pf[j][trkindex].energy()));

    ClusterSequence cs(particles, jet_def);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
    Selector select_pt       = SelectorPtMin(TrackJet_PtCut);
    vector<PseudoJet> selected_jets = select_pt(jets);
//    cout << "Clustered with " << jet_def.description() << endl;
    N_Jet[j] = selected_jets.size();
    N_Jet_total += selected_jets.size();
    for(size_t jetindex=0; jetindex < selected_jets.size(); jetindex++){
      TLorentzVector Jet_P4(selected_jets[jetindex].px(),selected_jets[jetindex].py(),selected_jets[jetindex].pz(),selected_jets[jetindex].E());
      TrackJet_P4_thisvertex.push_back(Jet_P4);
    }
    nBranches_->TrackJet_P4.push_back(TrackJet_P4_thisvertex);
  }

//Calculate event shape observables: Spherocity, Thrust, Broadening

  Float_t Spherocity, Thrust, Broaden;
  TVector3 taxis;

  for(size_t j = 0; j < vertices_->size(); j++){
    calcS_T_B(alltracks_pf[j], Spherocity,Thrust,Broaden, taxis, min(nSeed_ThrustCalculation,(int)alltracks_pf[j].size()));
    nBranches_->Instanton_Trk_spherocity.push_back(Spherocity);
    nBranches_->Instanton_Trk_thrust.push_back(Thrust);
    nBranches_->Instanton_Trk_broaden.push_back(Broaden);
  }

/*  if(runOnMC_){
    cout<<"number of packedgenparticles "<<packedgenParticles_->size()<<", number of genparticles "<<genParticles_->size()<<", number of tracks "<<alltracks_pf.size()<<endl;
    int N_finalgenparticles=0;
    for(size_t kk=0; kk< packedgenParticles_->size();kk++){
    if (!(*packedgenParticles_)[kk].isPromptFinalState()||!((*packedgenParticles_)[kk].pdgId()==211||(*packedgenParticles_)[kk].pdgId()==-211)) continue;
    N_finalgenparticles++;
    bool matched=false;
    for(size_t ivertex=0; ivertex< alltracks_pf.size(); ivertex++){
      for (size_t ipf=0; ipf< alltracks_pf[ivertex].size(); ipf++){
        if( (*packedgenParticles_)[kk].pdgId()==alltracks_pf[ivertex][ipf].pdgId() ){
          Float_t dr = reco::deltaR((*packedgenParticles_)[kk].eta(),(*packedgenParticles_)[kk].phi(),
                                     alltracks_pf[ivertex][ipf].eta(),alltracks_pf[ivertex][ipf].phi());
          Float_t ptratio = (*packedgenParticles_)[kk].pt()/alltracks_pf[ivertex][ipf].bestTrack()->pt();
          if(ptratio > 0.85 && ptratio < 1.15 && dr<1.){
          //  break;
          }
        }
      }
    }
    if(matched) cout<<"matched to a track"<<endl;
    else cout<<"this particle unmatched, pdgId "<<(*packedgenParticles_)[kk].pdgId()<<endl;
    }  
    cout<<"number of final genparticle " <<N_finalgenparticles<<endl;

  }
*/
  std::vector<string> FiredTriggers;
  std::vector<string> FiredBPHTriggers;
  bool isBPHdata=false;
  const edm::TriggerNames& trigNames = event.triggerNames(*HLTtriggers_);
  for (unsigned int i = 0, n = HLTtriggers_->size(); i < n; ++i) {
    if(HLTtriggers_->accept(i)){
      FiredTriggers.push_back(trigNames.triggerName(i));
//    std::cout << "Trigger " << trigNames.triggerName(i) << ": " << (HLTtriggers_->accept(i) ? "PASS" : "fail (or not run)") << std::endl;
      if(isBPHtrigger(trigNames.triggerName(i))){
        isBPHdata=true;
        FiredBPHTriggers.push_back(trigNames.triggerName(i));
      }
    }
  }
  for(size_t imuon = 0; imuon < muons_->size(); ++ imuon){
    const pat::Muon & muon = (*muons_)[imuon];
    if(!(muon.track().isNonnull())) continue;
//    cout<<"muon pt="<<muon.pt()<<", eta="<<muon.eta()<<endl; 
    if(muon.pt() < 4) continue;
    if(TMath::Abs(muon.eta()) > 2.4) continue;
/*    for (pat::TriggerObjectStandAlone obj : *triggerObjects) {    
      obj.unpackPathNames(trigNames);
      obj.unpackFilterLabels(event, *HLTtriggers_);
      std::vector<std::string> pathNamesAll  = obj.pathNames(false);
      if (!pathNamesAll.size()) continue;
      if(!imuon){
        bool triggermuon=false;
        for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
          if((obj.pdgId()==13||obj.pdgId()==-13)&&isBPHtrigger(pathNamesAll[h])){
            cout<<"pathname for this obj: "<<pathNamesAll[h]<<endl;
            triggermuon=true;
          }
        }
        if(triggermuon)
          cout<<"trigger muon pt = "<<obj.pt()<<",phi = "<<obj.phi()<<", eta = "<<obj.eta()<<endl;
      }
    }
*/
//    cout<<"muon pt = "<<muon.pt()<<", phi = "<<muon.phi()<<", eta = "<<muon.eta()<<endl;

//assign good muons to vtx

    const reco::TrackRef track_muon = muon.muonBestTrack();
    reco::TransientTrack transienttrack_muon = (*builder).build(track_muon);
    double closestIP=999.;
    int closestvtxindex=-1;
    for( size_t nn = 0; nn < vertices_->size(); nn++){
      std::pair<bool,Measurement1D> IP3D = IPTools::absoluteImpactParameter3D(transienttrack_muon,(*vertices_)[nn]);
      if(IP3D.second.value() < closestIP){
        closestIP = IP3D.second.value();
        closestvtxindex = nn;
      }
      vtx_goodMuonIP[nn]=std::min(vtx_goodMuonIP[nn],(float)IP3D.second.value());
    }

    if(closestvtxindex!=-1)  vtx_N_goodMuon[closestvtxindex]++;
    else {
      std::cout<<"Warning: muon vertex unfound"<<std::endl;
      continue;
    }
    int pfvtxindex= -1;
    for(size_t ipfmuon = 0; ipfmuon < alltracks_pfmuon.size(); ++ ipfmuon){
      if(muon.pdgId()==alltracks_pfmuon[ipfmuon].pdgId()){
        double ptratio = muon.pt()/alltracks_pfmuon[ipfmuon].pt();
        double deltaR = reco::deltaR(muon.eta(), muon.phi(), alltracks_pfmuon[ipfmuon].eta(), alltracks_pfmuon[ipfmuon].phi());
        if(ptratio>0.95&&ptratio<1.05&&deltaR<0.015){
          pfvtxindex = alltracks_pfmuon[ipfmuon].vertexRef().key();
          vtx_N_goodPFMuon[pfvtxindex]++;
          for( size_t nn = 0; nn < vertices_->size(); nn++){
            vtx_goodPFMuon_PVDistance[nn] = std::min(vtx_goodPFMuon_PVDistance[nn],(float)std::sqrt(std::pow((*vertices_)[nn].position().x()-(*vertices_)[pfvtxindex].position().x(),2)+std::pow((*vertices_)[nn].position().y()-(*vertices_)[pfvtxindex].position().y(),2)+std::pow((*vertices_)[nn].position().z()-(*vertices_)[pfvtxindex].position().z(),2)));

          }
          break;
        }
      }
    }
//label vertice if the muons coming out of it fired the trigger
//    if(muon.triggered("HLT_Mu7_IP4_part0_v2"))
//    cout<<"triggered !"<<endl;
    if(!isBPHdata) continue;
    for(size_t mm=0; mm < FiredBPHTriggers.size(); mm++){
      if(muon.triggered(FiredBPHTriggers[mm].c_str())){// cout<<" this muon fired the BPH trigger"<<endl;
        vtx_isBPHtrigger_fromMuon[closestvtxindex] = true;
        if(pfvtxindex!= -1) vtx_isBPHtrigger_fromPFMuon[pfvtxindex] = true;
        break;
      }
    }

  }
  for( size_t nn = 0; nn < vertices_->size(); nn++){
    if((*vertices_)[nn].ndof()>=4. && (*vertices_)[nn].position().Rho()<=2.0 && fabs((*vertices_)[nn].position().Z())<=24.0){
      if(vtx_isBPHtrigger_fromMuon[nn]) PV_N_good_hasgoodMuon++;
      if(vtx_isBPHtrigger_fromPFMuon[nn]) PV_N_good_hasgoodPFMuon++;
    }
  }
  nBranches_->Instanton_N_Trk_total = N_Trk_total;
  nBranches_->Instanton_N_Trk_Displaced_total = N_Trk_Displaced_total;
  nBranches_->Instanton_N_Trk_goodDisplaced_total = N_Trk_goodDisplaced_total;
  nBranches_->Instanton_N_TrackJet_total = N_Jet_total;
  nBranches_->PV_N_good_hasgoodMuon = PV_N_good_hasgoodMuon;
  nBranches_->PV_N_good_hasgoodPFMuon = PV_N_good_hasgoodPFMuon;
  for( size_t jj = 0; jj < vertices_->size(); ++jj ){
    nBranches_->Instanton_N_Trk.push_back(N_Trk[jj]);
    nBranches_->Instanton_N_Trk_highPurity.push_back(N_Trk_highPurity[jj]);
    nBranches_->Instanton_N_Trk_highPurity_pt05.push_back(N_Trk_highPurity_pt05[jj]);
    nBranches_->Instanton_N_Trk_highPurity_pt08.push_back(N_Trk_highPurity_pt08[jj]);
    nBranches_->Instanton_N_Trk_highPurity_pt1.push_back(N_Trk_highPurity_pt1[jj]);
    nBranches_->Instanton_N_Trk_PVAssociationQuality0.push_back(N_Trk_PVAssociationQuality0[jj]);
    nBranches_->Instanton_N_Trk_PVAssociationQuality1.push_back(N_Trk_PVAssociationQuality1[jj]);
    nBranches_->Instanton_N_Trk_PVAssociationQuality4.push_back(N_Trk_PVAssociationQuality4[jj]);
    nBranches_->Instanton_N_Trk_PVAssociationQuality5.push_back(N_Trk_PVAssociationQuality5[jj]);
    nBranches_->Instanton_N_Trk_PVAssociationQuality6.push_back(N_Trk_PVAssociationQuality6[jj]);
    nBranches_->Instanton_N_Trk_PVAssociationQuality7.push_back(N_Trk_PVAssociationQuality7[jj]);
    nBranches_->Instanton_N_Trk_PVAssociationQualityLeq4_highPurity.push_back(N_Trk_PVAssociationQualityLeq4_highPurity[jj]);
    nBranches_->Instanton_N_Trk_PVAssociationQualityLeq4_highPurity_pt05.push_back(N_Trk_PVAssociationQualityLeq4_highPurity_pt05[jj]);
    nBranches_->Instanton_N_Trk_PVAssociationQualityLeq4_highPurity_pt08.push_back(N_Trk_PVAssociationQualityLeq4_highPurity_pt08[jj]);
    nBranches_->Instanton_N_Trk_PVAssociationQualityLeq4_highPurity_pt1.push_back(N_Trk_PVAssociationQualityLeq4_highPurity_pt1[jj]);
    nBranches_->Instanton_N_Trk_Displaced.push_back(N_Trk_Displaced[jj]);
    nBranches_->Instanton_N_Trk_Displaced_highPurity.push_back(N_Trk_Displaced_highPurity[jj]);
    nBranches_->Instanton_N_Trk_Displaced_highPurity_pt05.push_back(N_Trk_Displaced_highPurity_pt05[jj]);
    nBranches_->Instanton_N_Trk_Displaced_highPurity_pt08.push_back(N_Trk_Displaced_highPurity_pt08[jj]);
    nBranches_->Instanton_N_Trk_Displaced_highPurity_pt1.push_back(N_Trk_Displaced_highPurity_pt1[jj]);
    nBranches_->Instanton_N_Trk_Displaced_PVAssociationQuality0.push_back(N_Trk_Displaced_PVAssociationQuality0[jj]);
    nBranches_->Instanton_N_Trk_Displaced_PVAssociationQuality1.push_back(N_Trk_Displaced_PVAssociationQuality1[jj]);
    nBranches_->Instanton_N_Trk_Displaced_PVAssociationQuality4.push_back(N_Trk_Displaced_PVAssociationQuality4[jj]);
    nBranches_->Instanton_N_Trk_Displaced_PVAssociationQuality5.push_back(N_Trk_Displaced_PVAssociationQuality5[jj]);
    nBranches_->Instanton_N_Trk_Displaced_PVAssociationQuality6.push_back(N_Trk_Displaced_PVAssociationQuality6[jj]);
    nBranches_->Instanton_N_Trk_Displaced_PVAssociationQuality7.push_back(N_Trk_Displaced_PVAssociationQuality7[jj]);
    nBranches_->Instanton_N_Trk_Displaced_PVAssociationQualityLeq4_highPurity.push_back(N_Trk_Displaced_PVAssociationQualityLeq4_highPurity[jj]);
    nBranches_->Instanton_N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt05.push_back(N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt05[jj]);
    nBranches_->Instanton_N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt08.push_back(N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt08[jj]);
    nBranches_->Instanton_N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt1.push_back(N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt1[jj]);
    nBranches_->Instanton_N_Trk_goodDisplaced.push_back(N_Trk_goodDisplaced[jj]);
    nBranches_->Instanton_N_Trk_goodDisplaced_highPurity.push_back(N_Trk_goodDisplaced_highPurity[jj]);
    nBranches_->Instanton_N_Trk_goodDisplaced_highPurity_pt05.push_back(N_Trk_goodDisplaced_highPurity_pt05[jj]);
    nBranches_->Instanton_N_Trk_goodDisplaced_highPurity_pt08.push_back(N_Trk_goodDisplaced_highPurity_pt08[jj]);
    nBranches_->Instanton_N_Trk_goodDisplaced_highPurity_pt1.push_back(N_Trk_goodDisplaced_highPurity_pt1[jj]);
    nBranches_->Instanton_N_Trk_goodDisplaced_PVAssociationQuality0.push_back(N_Trk_goodDisplaced_PVAssociationQuality0[jj]);
    nBranches_->Instanton_N_Trk_goodDisplaced_PVAssociationQuality1.push_back(N_Trk_goodDisplaced_PVAssociationQuality1[jj]);
    nBranches_->Instanton_N_Trk_goodDisplaced_PVAssociationQuality4.push_back(N_Trk_goodDisplaced_PVAssociationQuality4[jj]);
    nBranches_->Instanton_N_Trk_goodDisplaced_PVAssociationQuality5.push_back(N_Trk_goodDisplaced_PVAssociationQuality5[jj]);
    nBranches_->Instanton_N_Trk_goodDisplaced_PVAssociationQuality6.push_back(N_Trk_goodDisplaced_PVAssociationQuality6[jj]);
    nBranches_->Instanton_N_Trk_goodDisplaced_PVAssociationQuality7.push_back(N_Trk_goodDisplaced_PVAssociationQuality7[jj]);
    nBranches_->Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity.push_back(N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity[jj]);
    nBranches_->Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt05.push_back(N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt05[jj]);
    nBranches_->Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt08.push_back(N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt08[jj]);
    nBranches_->Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt1.push_back(N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt1[jj]);
    nBranches_->Instanton_Trk_mass.push_back(Trk_mass[jj]);
    nBranches_->Instanton_N_TrackJet.push_back(N_Jet[jj]);
    nBranches_->Instanton_vtx_N_goodMuon.push_back(vtx_N_goodMuon[jj]);
    nBranches_->Instanton_vtx_N_goodPFMuon.push_back(vtx_N_goodPFMuon[jj]);
    nBranches_->Instanton_vtx_isBPHtrigger_fromMuon.push_back(vtx_isBPHtrigger_fromMuon[jj]);
    nBranches_->Instanton_vtx_isBPHtrigger_fromPFMuon.push_back(vtx_isBPHtrigger_fromPFMuon[jj]);
    nBranches_->Instanton_vtx_goodMuonIP.push_back(vtx_goodMuonIP[jj]);
    nBranches_->Instanton_vtx_goodPFMuon_PVDistance.push_back(vtx_goodPFMuon_PVDistance[jj]);
  }
  return true;
}
