#include "../interface/InstantonNtuplizer.h"
#include "TMatrixTBase.h"
#include "TMatrixDSymEigen.h"
#include "TVector3.h"
#include "../interface/helper.h"
#include "../interface/PPSCut.h"
#include "../interface/PPSEff.h"
#include "math.h"
#include "../interface/PPSTrackingEfficiency.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "CondFormats/RunInfo/interface/LHCInfo.h"
#include "CondFormats/DataRecord/interface/LHCInfoRcd.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "DataFormats/CTPPSReco/interface/TotemRPRecHit.h"
#include "DataFormats/CTPPSReco/interface/TotemRPUVPattern.h"
#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"
#include "DataFormats/ProtonReco/interface/ForwardProton.h"
#include "DataFormats/ProtonReco/interface/ForwardProtonFwd.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLiteFwd.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "RecoPPS/ProtonReconstruction/interface/ProtonReconstructionAlgorithm.h"
#include "CondFormats/DataRecord/interface/CTPPSInterpolatedOpticsRcd.h"
#include "CondFormats/PPSObjects/interface/LHCInterpolatedOpticalFunctionsSetCollection.h"

#include "Geometry/Records/interface/VeryForwardRealGeometryRecord.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/CTPPSGeometry.h"

#include <cstring>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include <stdlib.h>
using namespace edm;
using namespace std;
using namespace reco;
using namespace pat;
template <typename T>
//void calcS_T_B(std::vector<pat::PackedCandidate>& momenta, Float_t& Spherocity,Float_t& Thrust,Float_t& Broaden, TVector3& taxis, int nSeed){
void calcS_T_B(std::vector<T>& momenta, Float_t& Spherocity,Float_t& Thrust,Float_t& Broaden,Float_t& TransverseSpherocity,Float_t& TransverseThrust, Float_t& TransverseSpherocity_pTunweighted, TVector3& taxis, TVector3& transtaxis, int nSeed){
  float PI = 3.141592653;
  if (momenta.size()<2) {Spherocity=-1;Thrust=-1;Broaden=-1;TransverseSpherocity=-1;return;}
  vector<TLorentzVector> p;
  vector<TVector3> units;
  vector<TVector3> pt;
  float pt_scalarsum = 0.;
  TLorentzVector p_total;
  for(size_t trkindex=0; trkindex<momenta.size(); trkindex++){
    TLorentzVector momenta_p4(momenta[trkindex].px(),momenta[trkindex].py(),momenta[trkindex].pz(),momenta[trkindex].energy());
    units.push_back(TVector3(momenta[trkindex].px(),momenta[trkindex].py(),0.0).Unit());
    pt.push_back(TVector3(momenta[trkindex].px(),momenta[trkindex].py(),0.0));
    p.push_back(momenta_p4);
    p_total+=momenta_p4;
    pt_scalarsum += momenta[trkindex].pt();
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

//Calculate the Transverse Spherocity
  TVector3 tsaxis;
  float  tsph = 99999.;
  for (const TVector3& u : units){
    float s = 0.0;
    for (size_t k=0; k<pt.size(); k++) s += fabs(pt[k].Cross(u).Mag());
    if (s < tsph){
      tsph = s;
      tsaxis = u;
    }
  }
  TransverseSpherocity = 0.25*PI*PI*tsph*tsph/pt_scalarsum/pt_scalarsum;

//Calculate the pT unweighted Transverse Spherocity 
  TVector3 tsaxis_pTunweight;
  float  tsph_pTunweight = 99999.;
  for (const TVector3& u : units){
    float s = 0.0;
    for (size_t k=0; k<units.size(); k++) s += fabs(units[k].Cross(u).Mag());
    if (s < tsph_pTunweight){
      tsph_pTunweight = s;
      tsaxis_pTunweight = u;
    }
  }
  TransverseSpherocity_pTunweighted = 0.25*PI*PI*tsph_pTunweight*tsph_pTunweight/units.size()/units.size();
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
//Calculate the transverse thrust
  std::sort(pt.begin(),pt.end(),[](TVector3 a, TVector3 b){
    return a.Pt()>b.Pt();
  });
  TransverseThrust = 0;
  Float_t TransThr;
  vector<Float_t> TransThr_val;
  vector<TVector3> TransThr_vec;
  for (unsigned int j=0; j < pow(2,nSeed); j++){
    int sign=j;
    TVector3 foo(0,0,0);
    for (int k = 0 ; k < nSeed ; k++) {
      (sign % 2) == 1 ? foo += pt[k] : foo -= pt[k];
      sign /= 2;
    }
    foo=foo.Unit();
    double diff=999.;
    while (diff > 1e-5) {
      TVector3 foobar(0,0,0);
      for (size_t k=0 ; k<momenta.size() ; k++)
        foo.Dot(pt[k])>0 ? foobar+=pt[k] : foobar-=pt[k];
      diff = (foo-foobar.Unit()).Mag();
      foo=foobar.Unit();
    }
    TransThr=0.;
    for (size_t k=0 ; k<momenta.size() ; k++)
      TransThr+=abs(foo.Dot(pt[k]));
    TransThr_val.push_back(TransThr);
    TransThr_vec.push_back(foo);
  }
  TransverseThrust=0.;
  for (unsigned int j=0; j < TransThr_val.size(); j++){
    if (TransThr_val[j]>TransverseThrust){
      TransverseThrust=TransThr_val[j];
      transtaxis=TransThr_vec[j];
    }
  }
  TransverseThrust=1-TransverseThrust/pt_scalarsum;


  Broaden=0.;
//Calculate the Broaden
  for(size_t j=0; j< momenta.size(); j++)
    Broaden+=abs(taxis.Cross(p[j].Vect()).Mag());
  Broaden/=p_sum*2;
}

unsigned int Search_LumiBlock(Json::Value& Run, int Block, string RunNumber){
  if(Run[RunNumber][0u][0u]>Block||Run[RunNumber][Run[RunNumber].size()-1u][0u]<Block) return 0;
  unsigned int index_low = 0;
  unsigned int index_high = Run[RunNumber].size()-1;
  unsigned int index_mid = index_high/2;
  while(Run[RunNumber][index_mid][0u] != Block && index_mid != index_low && index_mid != index_high){
    if (Run[RunNumber][index_mid][0u] < Block) index_low=index_mid;
    else index_high=index_mid;
    index_mid = (index_low+index_high)/2;
  }
  if(Run[RunNumber][index_mid][0u] == Block) return index_mid;
  else {cout<<"lumi block unfound"<<endl;return 0;}
}
template <typename T>
void calc_rapidity(std::vector<T>& momenta, Float_t& rapiditygap,Float_t& xi_x,Float_t& xi_y){
  if (momenta.size()<2) {rapiditygap=-1;xi_x=-1;xi_y=-1;return;}
  vector<TLorentzVector> p;
  TLorentzVector p_total;
  for(size_t trkindex=0; trkindex<momenta.size(); trkindex++){
    TLorentzVector momenta_p4(momenta[trkindex].px(),momenta[trkindex].py(),momenta[trkindex].pz(),momenta[trkindex].energy());
    p.push_back(momenta_p4);
    p_total+=momenta_p4;
  }
  std::sort(p.begin(),p.end(),[](TLorentzVector a, TLorentzVector b){
    return a.Eta()>b.Eta();
  });
  Int_t gap_index=0;
  Float_t gap_max=0;
  for (size_t i=0; i< p.size()-1; i++) {
    if(p[i].Eta()-p[i+1].Eta()>gap_max){
      gap_max=p[i].Eta()-p[i+1].Eta();
      gap_index=i;
    }
  }
  rapiditygap=p[gap_index].Eta()-p[gap_index+1].Eta();
  TLorentzVector p_total_X;
  TLorentzVector p_total_Y;
  for(size_t ix=0; ix<=(size_t)gap_index; ix++)
    p_total_X+=p[ix];
  for(size_t iy=gap_index+1; iy<p.size(); iy++)
    p_total_Y+=p[iy];
  xi_x=p_total_X.M2()/p_total.M2();
  xi_y=p_total_Y.M2()/p_total.M2();
}



void FindDaughter(const reco::Candidate* particle, vector<const reco::Candidate*> &daughter){
  if(particle->status()==1&&particle->numberOfDaughters()==0){
//    cout<<"final state found"<<endl;
    daughter.push_back(particle);
    return;
  }
  else if(particle->numberOfDaughters()==0){
//    cout<<"0 daughters, status "<<particle->status()<<endl;
    return;
  }
  else{
//    cout<<particle->numberOfDaughters()<<" daughters"<<endl;
//    cout<<"first daughter "<<particle->daughter(0)->pdgId()<<endl;
    for(unsigned int d=0; d<particle->numberOfDaughters(); d++){
//      cout<<"daughter "<<d<<endl;
      FindDaughter(particle->daughter(d),daughter);
    }
  }
}

int cal_IPsign(const reco::Candidate* particle,double pv_x, double pv_y, double pv_z){
  if(particle->px()*(particle->vx()-pv_x)+particle->py()*(particle->vy()-pv_y)+particle->pz()*(particle->vz()-pv_z))
    return 1;
  else
    return -1;
}

float cal_dxy(const reco::Candidate* particle, double pv_x, double pv_y, double pv_z){
  return cal_IPsign(particle,pv_x,pv_y,pv_z)*fabs((particle->vx()-pv_x)*particle->py()-(particle->vy()-pv_y)*particle->px())/particle->pt();
}

float cal_dz(const reco::Candidate* particle, double pv_x, double pv_y, double pv_z){
  return cal_IPsign(particle,pv_x,pv_y,pv_z)*fabs(particle->vz()-pv_z-particle->pz()*(particle->px()*(particle->vx()-pv_x)+particle->py()*(particle->vy()-pv_y))/(particle->pt()*particle->pt()));
}

float cal_IP3D(const reco::Candidate* particle, double pv_x, double pv_y, double pv_z){
  return cal_IPsign(particle,pv_x,pv_y,pv_z)*sqrt((pow((particle->vx()-pv_x)*particle->py()-(particle->vy()-pv_y)*particle->px(),2)+pow((particle->vx()-pv_x)*particle->pz()-(particle->vz()-pv_z)*particle->px(),2)+pow((particle->vy()-pv_y)*particle->pz()-(particle->vz()-pv_z)*particle->py(),2))/(pow(particle->pt(),2)+pow(particle->pz(),2)));

}

void cal_p4_pt(std::vector<std::vector<pat::PackedCandidate>> pf_byvtx, vector<float> &mass, vector<float> &pt_scalarsum,vector<float> &pt,vector<float> &pz,vector<float> &eta){
  mass.clear();
  pt_scalarsum.clear();
  pt.clear();
  pz.clear();
  eta.clear();
  for(size_t ivtx=0; ivtx<pf_byvtx.size(); ivtx++){
    math::XYZTLorentzVector p4_vertex(0,0,0,0);
    float pt_vertex = 0.;
    for(size_t itrk=0; itrk<pf_byvtx[ivtx].size(); itrk++){
      p4_vertex += pf_byvtx[ivtx][itrk].p4();
      pt_vertex += pf_byvtx[ivtx][itrk].pt();
    }
    mass.push_back(p4_vertex.M());
    pt.push_back(p4_vertex.pt());
    pz.push_back(p4_vertex.pz());
    eta.push_back(p4_vertex.Eta());
    pt_scalarsum.push_back(pt_vertex);
  }
}
template <typename T>
void cal_TrackJet(std::vector<std::vector<T>> pf_byvtx,vector<int> &N_TrackJet, vector<vector<float>> &TrackJet_px,vector<vector<float>> &TrackJet_py,vector<vector<float>> &TrackJet_pz,vector<vector<float>> &TrackJet_pt,vector<vector<float>> &TrackJet_E,vector<vector<float>> &TrackJet_mass,vector<vector<float>> &TrackJet_phi,vector<vector<float>> &TrackJet_eta, vector<float> &TrackJet_deltaR, vector<float> &TrackJet_deltaphi, float TrackJet_PtCut,vector<float>& Jettiness1,vector<float>& Jettiness2,vector<float>& Jettiness3,vector<float>& Jettiness4,vector<float>& Jettiness5,vector<float>& Jettiness6,vector<float>& Jettiness7,vector<float>& Jettiness8,float fastjet_beta=1){
//  N_TrackJet.clear();
  using namespace fastjet;
  using namespace fastjet::contrib;
  double R=0.4;
  JetDefinition jet_def(antikt_algorithm, R);
  for(size_t ivtx=0; ivtx<pf_byvtx.size(); ivtx++){
    std::vector<PseudoJet> particles;
    for( size_t trkindex = 0; trkindex < pf_byvtx[ivtx].size(); trkindex++)
      particles.push_back(PseudoJet(pf_byvtx[ivtx][trkindex].px(),pf_byvtx[ivtx][trkindex].py(),pf_byvtx[ivtx][trkindex].pz(),pf_byvtx[ivtx][trkindex].energy()));

    ClusterSequence cs(particles, jet_def);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
    Selector select_pt       = SelectorPtMin(TrackJet_PtCut);
    vector<PseudoJet> selected_jets = select_pt(jets);
//    float jet_E_vtx=0,jet_pt_vtx=0,jet_px_vtx=0,jet_py_vtx=0,jet_pz_vtx=0;
//    for(size_t jetindex=0; jetindex < selected_jets.size(); jetindex++){
//      jet_E_vtx+=selected_jets[jetindex].E();
//      jet_pt_vtx+=selected_jets[jetindex].pt();
//      jet_px_vtx+=selected_jets[jetindex].px();
//      jet_py_vtx+=selected_jets[jetindex].py();
//      jet_pz_vtx+=selected_jets[jetindex].pz();
//    }
//    jet_E.push_back(jet_E_vtx);
//    jet_pt.push_back(jet_pt_vtx);
//    jet_px.push_back(jet_px_vtx);
//    jet_py.push_back(jet_py_vtx);
//    jet_pz.push_back(jet_pz_vtx);
    N_TrackJet.push_back((int)selected_jets.size());
    for(size_t jetindex=0; jetindex < selected_jets.size(); jetindex++){
      TrackJet_px[ivtx].push_back(selected_jets[jetindex].px());
      TrackJet_py[ivtx].push_back(selected_jets[jetindex].py());
      TrackJet_pz[ivtx].push_back(selected_jets[jetindex].pz());
      TrackJet_pt[ivtx].push_back(selected_jets[jetindex].pt());
      TrackJet_E[ivtx].push_back(selected_jets[jetindex].E());
      TrackJet_mass[ivtx].push_back(selected_jets[jetindex].m());
      TrackJet_phi[ivtx].push_back(selected_jets[jetindex].phi());
      TrackJet_eta[ivtx].push_back(selected_jets[jetindex].eta());
    }
    if(selected_jets.size()>1){
      TVector3 v1(selected_jets[0].px(),selected_jets[0].py(),selected_jets[0].pz());
      TVector3 v2(selected_jets[1].px(),selected_jets[1].py(),selected_jets[1].pz());
      TrackJet_deltaR.push_back(v1.Angle(v2));
      TVector3 v1t(selected_jets[0].px(),selected_jets[0].py(),0);
      TVector3 v2t(selected_jets[1].px(),selected_jets[1].py(),0);
      TrackJet_deltaphi.push_back(v1t.Angle(v2t));
    }
    else{
      TrackJet_deltaR.push_back(-1.0);
      TrackJet_deltaphi.push_back(-1.0);
    }
    Njettiness njettiness = Njettiness(OnePass_KT_Axes(),UnnormalizedMeasure(fastjet_beta));
    Jettiness1.push_back(njettiness.getTau(1,particles));
    Jettiness2.push_back(njettiness.getTau(2,particles));
    Jettiness3.push_back(njettiness.getTau(3,particles));
    Jettiness4.push_back(njettiness.getTau(4,particles));
    Jettiness5.push_back(njettiness.getTau(5,particles));
    Jettiness6.push_back(njettiness.getTau(6,particles));
    Jettiness7.push_back(njettiness.getTau(7,particles));
    Jettiness8.push_back(njettiness.getTau(8,particles));
  }
  return;
}

//===================================================================================================================
InstantonNtuplizer::CaloType InstantonNtuplizer::GetTowerSubdetHad(double&eta)
{
  if(eta > -maxEtaHF && eta < -minEtaHF) return kHFm;
  if(eta >  minEtaHF && eta <  maxEtaHF) return kHFp;

  if(fabs(eta) > 0 && fabs(eta) < maxEtaHB) return kHB;

  if(fabs(eta) > minEtaHE && fabs(eta) < maxEtaHE) return kHE;

  return nCaloTypes;
}

InstantonNtuplizer::CaloType InstantonNtuplizer::GetTowerSubdetEm(double&eta)
{

  if(fabs(eta) > 0 && fabs(eta) < maxEtaEB) return kEB;

  if(fabs(eta) > minEtaEE && fabs(eta) < maxEtaEE) return kEE;

  return nCaloTypes;
}

InstantonNtuplizer::InstantonNtuplizer(  edm::EDGetTokenT<pat::MuonCollection>    muonToken   ,
                       edm::EDGetTokenT<reco::VertexCollection> verticeToken,
                       edm::EDGetTokenT<reco::BeamSpot>             beamToken ,
                       edm::EDGetTokenT<pat::PackedCandidateCollection> packedpfcandidatesToken,
                       edm::EDGetTokenT<edm::TriggerResults> triggertoken,
                       edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerobject,
                       edm::EDGetTokenT<reco::GenParticleCollection> genptoken,
                       edm::EDGetTokenT<pat::PackedGenParticleCollection> packedgenptoken,
                       edm::EDGetTokenT<pat::JetCollection>     jettoken,
                       edm::EDGetTokenT<std::vector<reco::VertexCompositePtrCandidate> > svToken,
                       edm::EDGetTokenT<edm::SortedCollection<CaloTower>> CaloTowerCollection,
                       edm::EDGetTokenT<reco::ForwardProtonCollection> ctppsProton_single_rpToken,
                       edm::EDGetTokenT<reco::ForwardProtonCollection> ctppsProton_multi_rpToken,
                       std::map< std::string, bool >& runFlags,
                       std::map< std::string, double >& runValues,
                       std::map< std::string, std::string >& runStrings,
                       Json::Value& instan_lumi,
                       NtupleBranches* nBranches)
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
  , svToken_   (svToken)
  , CaloTowerCollection_ ( CaloTowerCollection )
  , ctppsProton_single_rpToken_ (ctppsProton_single_rpToken) 
  , ctppsProton_multi_rpToken_ (ctppsProton_multi_rpToken)
  , runOnMC_   (runFlags["runOnMC"])
  , runOnMCPU_ (runFlags["runOnMCPU"])
  , runOnZeroBias_ (runFlags["runOnZeroBias"])
  , runOnHerwigInstanton_ (runFlags["runOnHerwigInstanton"])
  , runOnSherpaInstanton_ (runFlags["runOnSherpaInstanton"])
  , runOnMCTest_ (runFlags["runOnMCTest"])
  , verbose_   (runFlags["verbose"])
  , doTrack_   (runFlags["doTrack"])
  , doTrackJet_   (runFlags["doTrackJet"])
  , doForwardProton_   (runFlags["doForwardProton"])
  , trackrandomdrop_ (runFlags["TrackRandomDrop"])
  , instan_lumi_ (instan_lumi)
  , runStrings_ (runStrings)
{
  int year=2018;
  if(runStrings_["run"].find("18") != std::string::npos) year = 2018;
  else if(runStrings_["run"].find("17") != std::string::npos) year = 2017;
  else if(runStrings_["run"].find("16") != std::string::npos) year = 2016;
  forwardProtonEff_ = new PPSTrackingEfficiency(year); 
  if(verbose_){
    std::cout << "[InstantonNtuplizer] runOnMC    = " << runOnMC_ << std::endl;
  }
}
InstantonNtuplizer::~InstantonNtuplizer(void){

}
bool InstantonNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){
  if(verbose_) std::cout << "[InstantonNtuplizer] ---------------- event, run, lumi = " << event.id().event() << " " << event.id().run() << " " << event.id().luminosityBlock() << "----------------" << std::endl;

  if(!runOnMC_){
    unsigned int lumi_index = Search_LumiBlock(instan_lumi_, (unsigned int)event.id().luminosityBlock(),std::to_string(event.id().run()));
    std::stringstream mean_ss;
    mean_ss << instan_lumi_[std::to_string(event.id().run())][lumi_index][3u];
//    nBranches_->Instan_Lumi_per_bunch_mean = std::stof(mean_ss.str());
    nBranches_->Instan_Lumi_per_bunch_mean = std::atof(mean_ss.str().c_str());


    std::stringstream RMS_ss;
    RMS_ss << instan_lumi_[std::to_string(event.id().run())][lumi_index][2u];
//    nBranches_->Instan_Lumi_per_bunch_RMS = std::stof(RMS_ss.str());
    nBranches_->Instan_Lumi_per_bunch_RMS = std::atof(RMS_ss.str().c_str());
  }

//TrackJet calculating parameters
//  double R=0.4;
//  JetDefinition jet_def(antikt_algorithm, R);
  float TrackJet_PtCut=0.0; //Pt Cut > 3.0 GeV
  int nSeed_ThrustCalculation=3;

//Get the Tokens from MINIAOD
  double PV_X = 0;
  double PV_Y = 0;
  double PV_Z = 0;
//  std::cout<<"runOnMC_"<<(int)runOnMC_<<", runOnMCPU_"<<(int)runOnMCPU_<<std::endl;
  vector<pat::PackedGenParticle> ChargedFSParticle;
  vector<pat::PackedGenParticle> ChargedFSParticle_eta2p4;
  vector<pat::PackedGenParticle> ChargedFSParticle_eta2p4pt05;
  vector<pat::PackedGenParticle> ChargedFSParticle_eta2p4pt035;
  vector<pat::PackedGenParticle> ChargedFSParticle_eta05pt05;
  vector<pat::PackedGenParticle> ChargedFSParticle_eta1pt05;
  vector<pat::PackedGenParticle> ChargedFSParticle_eta1p5pt05;
  vector<pat::PackedGenParticle> ChargedFSParticle_eta2pt05;
  vector<pat::PackedGenParticle> ChargedFSParticle_eta05to1pt05;
  vector<pat::PackedGenParticle> ChargedFSParticle_eta1to1p5pt05;
  vector<pat::PackedGenParticle> ChargedFSParticle_eta1p5to2pt05;
  vector<pat::PackedGenParticle> ChargedFSParticle_eta2to2p4pt05;
  vector<vector<pat::PackedGenParticle>> ChargedFSParticle_eta2p4pt05_jet;
  if(runOnMC_){
    event.getByToken(genParticlesToken_ , genParticles_);
    event.getByToken(packedgenParticlesToken_ , packedgenParticles_);//packed particles all have status 1 (final state)

//To find the instanton positions
//For Herwig events, it searches for the gluons with number of daughters > 7, so the selected gluons go through instanton processes
    if(runOnHerwigInstanton_){
      float ins_px=0, ins_py=0, ins_pz=0, ins_E=0;
      for( unsigned p=0; p<genParticles_->size(); ++p ){
        if((*genParticles_)[p].numberOfDaughters()>7&&(*genParticles_)[p].pdgId()==21){
          ins_px+=(*genParticles_)[p].px();
          ins_py+=(*genParticles_)[p].py();
          ins_pz+=(*genParticles_)[p].pz();
          ins_E+=(*genParticles_)[p].energy();
          PV_X=(*genParticles_)[p].vx();
          PV_Y=(*genParticles_)[p].vy();
          PV_Z=(*genParticles_)[p].vz();
        }
        nBranches_->Instanton_VX=PV_X;
        nBranches_->Instanton_VY=PV_Y;
        nBranches_->Instanton_VZ=PV_Z;
        nBranches_->Instanton_PX=ins_px;
        nBranches_->Instanton_PY=ins_py;
        nBranches_->Instanton_PZ=ins_pz;
        nBranches_->Instanton_E=ins_E;
        nBranches_->Instanton_MASS=sqrt(ins_E*ins_E-ins_px*ins_px-ins_py*ins_py-ins_pz*ins_pz);
        nBranches_->Instanton_RAPIDITY=0.5*log((ins_E+sqrt(ins_px*ins_px+ins_py*ins_py+ins_pz*ins_pz))/(ins_E-sqrt(ins_px*ins_px+ins_py*ins_py+ins_pz*ins_pz)));
      }

    }
//Sherpa simulates the instantons in a different way: the instanton is considered as a particle (pid=999) instead of a process. So we search for the generated instantons by its pid and store its position.
    else if(runOnSherpaInstanton_){
//      bool findins = false;
      for( unsigned p=0; p<genParticles_->size(); ++p ){
        if((*genParticles_)[p].pdgId()==999){
//          findins = true;
          PV_X=(*genParticles_)[p].vx();
          PV_Y=(*genParticles_)[p].vy();
          PV_Z=(*genParticles_)[p].vz();
          nBranches_->Instanton_VX=(*genParticles_)[p].vx();
          nBranches_->Instanton_VY=(*genParticles_)[p].vy();
          nBranches_->Instanton_VZ=(*genParticles_)[p].vz();
          nBranches_->Instanton_PX=(*genParticles_)[p].px();
          nBranches_->Instanton_PY=(*genParticles_)[p].py();
          nBranches_->Instanton_PZ=(*genParticles_)[p].pz();
          nBranches_->Instanton_E=(*genParticles_)[p].energy();
          nBranches_->Instanton_MASS=(*genParticles_)[p].mass();
          nBranches_->Instanton_RAPIDITY=(*genParticles_)[p].rapidity();
          break;
        }
      }
//      if(findins == false) cout<<"Cannot find the generated instanton !"<<endl;
    }
    else if(runOnMCTest_){
      for( unsigned p=0; p<genParticles_->size(); ++p ){
//        if((*genParticles_)[p].numberOfMothers()==2 && (*genParticles_)[p].mother(0)->pdgId()==2212 && (*genParticles_)[p].mother(1)->pdgId()==2212){
        if((*genParticles_)[p].numberOfDaughters()>7&&(*genParticles_)[p].pdgId()==21){
//          cout<<"pdgId="<<(*genParticles_)[p].pdgId()<<", vx="<<(*genParticles_)[p].vx()<<", vy="<<(*genParticles_)[p].vy()<<", vz="<<(*genParticles_)[p].vz()<<endl;
//          nBranches_->Instanton_VX=(*genParticles_)[p].vx();
//          nBranches_->Instanton_VY=(*genParticles_)[p].vy();
//          nBranches_->Instanton_VZ=(*genParticles_)[p].vz();
//          break;
        }
        if((*genParticles_)[p].pdgId()==999) cout<<"Instanton position "<<(*genParticles_)[p].vx()<<" "<<(*genParticles_)[p].vy()<<" "<<(*genParticles_)[p].vz()<<endl;
      }
    }
/*
    if(runOnHerwigInstanton_||runOnSherpaInstanton_){
      vector<reco::GenParticle> genParticles_fscharge;
      for( unsigned p=0; p<genParticles_->size(); ++p ){
        if((*genParticles_)[p].status()==1&&(*genParticles_)[p].charge()!=0&&(*genParticles_)[p].eta()>-2.4&&(*genParticles_)[p].eta()<2.4&&(*genParticles_)[p].pt()>=0.5){
          genParticles_fscharge.push_back((*genParticles_)[p]);
//          cout<<"pdgId="<<(*genParticles_)[p].pdgId()<<" charge="<<(*genParticles_)[p].charge()<<" px="<<(*genParticles_)[p].px()<<" py="<<(*genParticles_)[p].py()<<" pz="<<(*genParticles_)[p].pz()<<" pt="<<(*genParticles_)[p].pt()<<endl;
        }
      }
      Float_t gen_sph,gen_thrust,gen_broaden;
      TVector3 gen_taxis;
      calcS_T_B(genParticles_fscharge, gen_sph,gen_thrust,gen_broaden, gen_taxis, min(nSeed_ThrustCalculation,(int)genParticles_fscharge.size()));
//      cout<<gen_sph<<" "<<gen_thrust<<" "<<gen_broaden<<" "<<(int)genParticles_fscharge.size()<<endl;
    }
*/
    double gen_sum_eT_HFp=0;
    double gen_sum_eT_HFm=0;
    double gen_sum_energy_HFp=0;
    double gen_sum_energy_HFm=0;
    double gen_leading_energy_HFp=0;
    double gen_leading_energy_HFm=0;
    for( unsigned p=0; p<packedgenParticles_->size(); ++p ){
      if((*packedgenParticles_)[p].status()==1){
        nBranches_->Instanton_N_gen_FSParticle++;
        double eta = (*packedgenParticles_)[p].eta();
        CaloType subdetHad = GetTowerSubdetHad(eta);
        if(subdetHad == kHFp){
          gen_sum_eT_HFp += (*packedgenParticles_)[p].et();
          gen_sum_energy_HFp += (*packedgenParticles_)[p].energy();
          if((*packedgenParticles_)[p].energy()>gen_leading_energy_HFp) gen_leading_energy_HFp=(*packedgenParticles_)[p].energy();
        }
        if(subdetHad == kHFm){
          gen_sum_eT_HFm += (*packedgenParticles_)[p].et();
          gen_sum_energy_HFm += (*packedgenParticles_)[p].energy();
          if((*packedgenParticles_)[p].energy()>gen_leading_energy_HFm) gen_leading_energy_HFm=(*packedgenParticles_)[p].energy();
        }
        if((*packedgenParticles_)[p].charge()!=0){
          nBranches_->Instanton_N_gen_ChargedFSParticle++;
          ChargedFSParticle.push_back((*packedgenParticles_)[p]);
          if((*packedgenParticles_)[p].eta()<2.4&&(*packedgenParticles_)[p].eta()>-2.4){ 
//          if((*packedgenParticles_)[p].eta()<4&&(*packedgenParticles_)[p].eta()>-4){ //TODO FOR DEBUGGING. DO NOT USE IN SAMPLE PRODUCTION!!! TODO
            nBranches_->Instanton_N_gen_ChargedFSParticle_eta2p4++;
            nBranches_->ChargedFSParticle_eta2p4_px.push_back((*packedgenParticles_)[p].px());
            nBranches_->ChargedFSParticle_eta2p4_py.push_back((*packedgenParticles_)[p].py());
            nBranches_->ChargedFSParticle_eta2p4_pz.push_back((*packedgenParticles_)[p].pz());
            nBranches_->ChargedFSParticle_eta2p4_p0.push_back((*packedgenParticles_)[p].energy());
            nBranches_->ChargedFSParticle_eta2p4_pt.push_back((*packedgenParticles_)[p].pt());
            nBranches_->ChargedFSParticle_eta2p4_mass.push_back((*packedgenParticles_)[p].mass());
            nBranches_->ChargedFSParticle_eta2p4_eta.push_back((*packedgenParticles_)[p].eta());
            nBranches_->ChargedFSParticle_eta2p4_phi.push_back((*packedgenParticles_)[p].phi());
            nBranches_->ChargedFSParticle_eta2p4_pdgId.push_back((*packedgenParticles_)[p].pdgId());
            nBranches_->ChargedFSParticle_eta2p4_charge.push_back((*packedgenParticles_)[p].charge());
            nBranches_->ChargedFSParticle_eta2p4_matched.push_back(0);
            ChargedFSParticle_eta2p4.push_back((*packedgenParticles_)[p]);
            if((*packedgenParticles_)[p].pt()>0.5) {
              nBranches_->Instanton_N_gen_ChargedFSParticle_eta2p4pt05++;
              ChargedFSParticle_eta2p4pt05.push_back((*packedgenParticles_)[p]);  
              if((*packedgenParticles_)[p].eta()<2&&(*packedgenParticles_)[p].eta()>-2){
                  nBranches_->Instanton_N_gen_ChargedFSParticle_eta2pt05++;
                  ChargedFSParticle_eta2pt05.push_back((*packedgenParticles_)[p]);
                 
                if((*packedgenParticles_)[p].eta()<1.5&&(*packedgenParticles_)[p].eta()>-1.5){
                    nBranches_->Instanton_N_gen_ChargedFSParticle_eta1p5pt05++;
                    ChargedFSParticle_eta1p5pt05.push_back((*packedgenParticles_)[p]);
                  if((*packedgenParticles_)[p].eta()<1&&(*packedgenParticles_)[p].eta()>-1){
                    nBranches_->Instanton_N_gen_ChargedFSParticle_eta1pt05++;
                    ChargedFSParticle_eta1pt05.push_back((*packedgenParticles_)[p]);
                    if((*packedgenParticles_)[p].eta()<0.5&&(*packedgenParticles_)[p].eta()>-0.5){
                      nBranches_->Instanton_N_gen_ChargedFSParticle_eta05pt05++;
                      ChargedFSParticle_eta05pt05.push_back((*packedgenParticles_)[p]);
                    }
                    else{
                      nBranches_->Instanton_N_gen_ChargedFSParticle_eta05to1pt05++;
                      ChargedFSParticle_eta05to1pt05.push_back((*packedgenParticles_)[p]);
                    }
                  }
                  else{
                    nBranches_->Instanton_N_gen_ChargedFSParticle_eta1to1p5pt05++;
                    ChargedFSParticle_eta1to1p5pt05.push_back((*packedgenParticles_)[p]);
                  }
                }
                else{
                  nBranches_->Instanton_N_gen_ChargedFSParticle_eta1p5to2pt05++;
                  ChargedFSParticle_eta1p5to2pt05.push_back((*packedgenParticles_)[p]);
                }
              }
              else{
                nBranches_->Instanton_N_gen_ChargedFSParticle_eta2to2p4pt05++;
                ChargedFSParticle_eta2to2p4pt05.push_back((*packedgenParticles_)[p]);
              }
            }
            if((*packedgenParticles_)[p].pt()>0.35) {
              nBranches_->Instanton_N_gen_ChargedFSParticle_eta2p4pt035++;
              ChargedFSParticle_eta2p4pt035.push_back((*packedgenParticles_)[p]);
            }
          }
        }
      }
    }
    ChargedFSParticle_eta2p4pt05_jet.push_back(ChargedFSParticle_eta2p4pt05);
    nBranches_->gen_sum_energy_HFp = gen_sum_energy_HFp;
    nBranches_->gen_sum_energy_HFm = gen_sum_energy_HFm;
    nBranches_->gen_sum_eT_HFp = gen_sum_eT_HFp;
    nBranches_->gen_sum_eT_HFm = gen_sum_eT_HFm;
    nBranches_->gen_leading_energy_HFp = gen_leading_energy_HFp;
    nBranches_->gen_leading_energy_HFm = gen_leading_energy_HFm;
    std::vector<std::vector<float>> chjet_eta2p4pt05_px(1);
    std::vector<std::vector<float>> chjet_eta2p4pt05_py(1);
    std::vector<std::vector<float>> chjet_eta2p4pt05_pz(1);
    std::vector<std::vector<float>> chjet_eta2p4pt05_pt(1);
    std::vector<std::vector<float>> chjet_eta2p4pt05_E(1);
    std::vector<std::vector<float>> chjet_eta2p4pt05_mass(1);
    std::vector<std::vector<float>> chjet_eta2p4pt05_phi(1);
    std::vector<std::vector<float>> chjet_eta2p4pt05_eta(1);
    std::vector<int> N_chjet_eta2p4pt05;
    std::vector<float> ch_eta2p4pt05_deltaR;
    std::vector<float> ch_eta2p4pt05_deltaphi;
    std::vector<float> ch_1jettiness;
    std::vector<float> ch_2jettiness;
    std::vector<float> ch_3jettiness;
    std::vector<float> ch_4jettiness;
    std::vector<float> ch_5jettiness;
    std::vector<float> ch_6jettiness;
    std::vector<float> ch_7jettiness;
    std::vector<float> ch_8jettiness;
    cal_TrackJet(ChargedFSParticle_eta2p4pt05_jet,N_chjet_eta2p4pt05, chjet_eta2p4pt05_px,chjet_eta2p4pt05_py,chjet_eta2p4pt05_pz,chjet_eta2p4pt05_pt,chjet_eta2p4pt05_E,chjet_eta2p4pt05_mass,chjet_eta2p4pt05_phi,chjet_eta2p4pt05_eta,ch_eta2p4pt05_deltaR,ch_eta2p4pt05_deltaphi, TrackJet_PtCut,ch_1jettiness,ch_2jettiness,ch_3jettiness,ch_4jettiness,ch_5jettiness,ch_6jettiness,ch_7jettiness,ch_8jettiness);
    nBranches_->Instanton_N_gen_ChargedFSParticle_eta2p4pt05_jet = N_chjet_eta2p4pt05[0];
    nBranches_->ChargedFSParticle_eta2p4pt05_jet_deltaR = ch_eta2p4pt05_deltaR[0];
    nBranches_->ChargedFSParticle_eta2p4pt05_jet_deltaphi = ch_eta2p4pt05_deltaphi[0];
    nBranches_->Instanton_gen_1Jettiness_eta2p4pt05 = ch_1jettiness[0];
    nBranches_->Instanton_gen_2Jettiness_eta2p4pt05 = ch_2jettiness[0];
    nBranches_->Instanton_gen_3Jettiness_eta2p4pt05 = ch_3jettiness[0];
    nBranches_->Instanton_gen_4Jettiness_eta2p4pt05 = ch_4jettiness[0];
    nBranches_->Instanton_gen_5Jettiness_eta2p4pt05 = ch_5jettiness[0];
    nBranches_->Instanton_gen_6Jettiness_eta2p4pt05 = ch_6jettiness[0];
    nBranches_->Instanton_gen_7Jettiness_eta2p4pt05 = ch_7jettiness[0];
    nBranches_->Instanton_gen_8Jettiness_eta2p4pt05 = ch_8jettiness[0];
    if(doTrackJet_){
      for(int ijet=0; ijet< N_chjet_eta2p4pt05[0];ijet++){
        nBranches_->ChargedFSParticle_eta2p4pt05_jet_px.push_back(chjet_eta2p4pt05_px[0][ijet]);
        nBranches_->ChargedFSParticle_eta2p4pt05_jet_py.push_back(chjet_eta2p4pt05_py[0][ijet]);
        nBranches_->ChargedFSParticle_eta2p4pt05_jet_pz.push_back(chjet_eta2p4pt05_pz[0][ijet]);
        nBranches_->ChargedFSParticle_eta2p4pt05_jet_pt.push_back(chjet_eta2p4pt05_pt[0][ijet]);
        nBranches_->ChargedFSParticle_eta2p4pt05_jet_E.push_back(chjet_eta2p4pt05_E[0][ijet]);
        nBranches_->ChargedFSParticle_eta2p4pt05_jet_mass.push_back(chjet_eta2p4pt05_mass[0][ijet]);
        nBranches_->ChargedFSParticle_eta2p4pt05_jet_phi.push_back(chjet_eta2p4pt05_phi[0][ijet]);
        nBranches_->ChargedFSParticle_eta2p4pt05_jet_eta.push_back(chjet_eta2p4pt05_eta[0][ijet]);
      }
    }
    Float_t gen_spherocity,gen_thrust,gen_broaden,gen_transversespherocity,gen_transversethrust,gen_transversespherocity_pTunweighted;
    Float_t gen_spherocity_eta2p4,gen_thrust_eta2p4,gen_broaden_eta2p4,gen_transversespherocity_eta2p4,gen_transversethrust_eta2p4,gen_transversespherocity_pTunweighted_eta2p4;
    Float_t gen_spherocity_eta2p4pt05,gen_thrust_eta2p4pt05,gen_broaden_eta2p4pt05,gen_transversespherocity_eta2p4pt05,gen_transversethrust_eta2p4pt05,gen_transversespherocity_pTunweighted_eta2p4pt05;
    Float_t gen_spherocity_eta2p4pt035,gen_thrust_eta2p4pt035,gen_broaden_eta2p4pt035,gen_transversespherocity_eta2p4pt035,gen_transversethrust_eta2p4pt035,gen_transversespherocity_pTunweighted_eta2p4pt035;
    Float_t gen_spherocity_eta05pt05,gen_thrust_eta05pt05,gen_broaden_eta05pt05,gen_transversespherocity_eta05pt05,gen_transversethrust_eta05pt05,gen_transversespherocity_pTunweighted_eta05pt05;
    Float_t gen_spherocity_eta1pt05,gen_thrust_eta1pt05,gen_broaden_eta1pt05,gen_transversespherocity_eta1pt05,gen_transversethrust_eta1pt05,gen_transversespherocity_pTunweighted_eta1pt05;
    Float_t gen_spherocity_eta2pt05,gen_thrust_eta2pt05,gen_broaden_eta2pt05,gen_transversespherocity_eta2pt05,gen_transversethrust_eta2pt05,gen_transversespherocity_pTunweighted_eta2pt05;
    Float_t gen_spherocity_eta05to1pt05,gen_thrust_eta05to1pt05,gen_broaden_eta05to1pt05,gen_transversespherocity_eta05to1pt05,gen_transversethrust_eta05to1pt05,gen_transversespherocity_pTunweighted_eta05to1pt05;
    Float_t gen_spherocity_eta1to1p5pt05,gen_thrust_eta1to1p5pt05,gen_broaden_eta1to1p5pt05,gen_transversespherocity_eta1to1p5pt05,gen_transversethrust_eta1to1p5pt05,gen_transversespherocity_pTunweighted_eta1to1p5pt05;
    Float_t gen_spherocity_eta1p5to2pt05,gen_thrust_eta1p5to2pt05,gen_broaden_eta1p5to2pt05,gen_transversespherocity_eta1p5to2pt05,gen_transversethrust_eta1p5to2pt05,gen_transversespherocity_pTunweighted_eta1p5to2pt05;
    Float_t gen_spherocity_eta2to2p4pt05,gen_thrust_eta2to2p4pt05,gen_broaden_eta2to2p4pt05,gen_transversespherocity_eta2to2p4pt05,gen_transversethrust_eta2to2p4pt05,gen_transversespherocity_pTunweighted_eta2to2p4pt05;
    TVector3 gen_taxis;
    TVector3 gen_taxis_eta2p4;
    TVector3 gen_taxis_eta2p4pt05;
    TVector3 gen_taxis_eta2p4pt035;
    TVector3 gen_taxis_eta05pt05;
    TVector3 gen_taxis_eta1pt05;
    TVector3 gen_taxis_eta2pt05;
    TVector3 gen_taxis_eta05to1pt05;
    TVector3 gen_taxis_eta1to1p5pt05;
    TVector3 gen_taxis_eta1p5to2pt05;
    TVector3 gen_taxis_eta2to2p4pt05;
    TVector3 gen_transtaxis;
    TVector3 gen_transtaxis_eta2p4;
    TVector3 gen_transtaxis_eta2p4pt05;
    TVector3 gen_transtaxis_eta2p4pt035;
    TVector3 gen_transtaxis_eta05pt05;
    TVector3 gen_transtaxis_eta1pt05;
    TVector3 gen_transtaxis_eta2pt05;
    TVector3 gen_transtaxis_eta05to1pt05;
    TVector3 gen_transtaxis_eta1to1p5pt05;
    TVector3 gen_transtaxis_eta1p5to2pt05;
    TVector3 gen_transtaxis_eta2to2p4pt05;
    calcS_T_B(ChargedFSParticle, gen_spherocity,gen_thrust,gen_broaden,gen_transversespherocity,gen_transversethrust,gen_transversespherocity_pTunweighted, gen_taxis,gen_transtaxis, min(nSeed_ThrustCalculation,(int)ChargedFSParticle.size())); 
    calcS_T_B(ChargedFSParticle_eta2p4, gen_spherocity_eta2p4,gen_thrust_eta2p4,gen_broaden_eta2p4,gen_transversespherocity_eta2p4,gen_transversethrust_eta2p4,gen_transversespherocity_pTunweighted_eta2p4, gen_taxis_eta2p4,gen_transtaxis_eta2p4, min(nSeed_ThrustCalculation,(int)ChargedFSParticle_eta2p4.size()));
    calcS_T_B(ChargedFSParticle_eta2p4pt05, gen_spherocity_eta2p4pt05,gen_thrust_eta2p4pt05,gen_broaden_eta2p4pt05,gen_transversespherocity_eta2p4pt05,gen_transversethrust_eta2p4pt05,gen_transversespherocity_pTunweighted_eta2p4pt05, gen_taxis_eta2p4pt05,gen_transtaxis_eta2p4pt05, min(nSeed_ThrustCalculation,(int)ChargedFSParticle_eta2p4pt05.size()));
    calcS_T_B(ChargedFSParticle_eta2p4pt035, gen_spherocity_eta2p4pt035,gen_thrust_eta2p4pt035,gen_broaden_eta2p4pt035,gen_transversespherocity_eta2p4pt035,gen_transversethrust_eta2p4pt035,gen_transversespherocity_pTunweighted_eta2p4pt035, gen_taxis_eta2p4pt035,gen_transtaxis_eta2p4pt035, min(nSeed_ThrustCalculation,(int)ChargedFSParticle_eta2p4pt035.size()));
    calcS_T_B(ChargedFSParticle_eta05pt05, gen_spherocity_eta05pt05,gen_thrust_eta05pt05,gen_broaden_eta05pt05,gen_transversespherocity_eta05pt05,gen_transversethrust_eta05pt05,gen_transversespherocity_pTunweighted_eta05pt05, gen_taxis_eta05pt05,gen_transtaxis_eta05pt05, min(nSeed_ThrustCalculation,(int)ChargedFSParticle_eta05pt05.size()));
    calcS_T_B(ChargedFSParticle_eta1pt05, gen_spherocity_eta1pt05,gen_thrust_eta1pt05,gen_broaden_eta1pt05,gen_transversespherocity_eta1pt05,gen_transversethrust_eta1pt05,gen_transversespherocity_pTunweighted_eta1pt05, gen_taxis_eta1pt05,gen_transtaxis_eta1pt05, min(nSeed_ThrustCalculation,(int)ChargedFSParticle_eta1pt05.size()));
    calcS_T_B(ChargedFSParticle_eta2pt05, gen_spherocity_eta2pt05,gen_thrust_eta2pt05,gen_broaden_eta2pt05,gen_transversespherocity_eta2pt05,gen_transversethrust_eta2pt05,gen_transversespherocity_pTunweighted_eta2pt05, gen_taxis_eta2pt05,gen_transtaxis_eta2pt05, min(nSeed_ThrustCalculation,(int)ChargedFSParticle_eta2pt05.size()));
    calcS_T_B(ChargedFSParticle_eta05to1pt05, gen_spherocity_eta05to1pt05,gen_thrust_eta05to1pt05,gen_broaden_eta05to1pt05,gen_transversespherocity_eta05to1pt05,gen_transversethrust_eta05to1pt05,gen_transversespherocity_pTunweighted_eta05to1pt05, gen_taxis_eta05to1pt05,gen_transtaxis_eta05to1pt05, min(nSeed_ThrustCalculation,(int)ChargedFSParticle_eta05to1pt05.size()));
    calcS_T_B(ChargedFSParticle_eta1to1p5pt05, gen_spherocity_eta1to1p5pt05,gen_thrust_eta1to1p5pt05,gen_broaden_eta1to1p5pt05,gen_transversespherocity_eta1to1p5pt05,gen_transversethrust_eta1to1p5pt05,gen_transversespherocity_pTunweighted_eta1to1p5pt05,gen_taxis_eta1to1p5pt05,gen_transtaxis_eta1to1p5pt05, min(nSeed_ThrustCalculation,(int)ChargedFSParticle_eta1to1p5pt05.size()));
    calcS_T_B(ChargedFSParticle_eta1p5to2pt05, gen_spherocity_eta1p5to2pt05,gen_thrust_eta1p5to2pt05,gen_broaden_eta1p5to2pt05,gen_transversespherocity_eta1p5to2pt05,gen_transversethrust_eta1p5to2pt05,gen_transversespherocity_pTunweighted_eta1p5to2pt05, gen_taxis_eta1p5to2pt05,gen_transtaxis_eta1p5to2pt05, min(nSeed_ThrustCalculation,(int)ChargedFSParticle_eta1p5to2pt05.size()));
    calcS_T_B(ChargedFSParticle_eta2to2p4pt05, gen_spherocity_eta2to2p4pt05,gen_thrust_eta2to2p4pt05,gen_broaden_eta2to2p4pt05,gen_transversespherocity_eta2to2p4pt05,gen_transversethrust_eta2to2p4pt05,gen_transversespherocity_pTunweighted_eta2to2p4pt05, gen_taxis_eta2to2p4pt05,gen_transtaxis_eta2to2p4pt05, min(nSeed_ThrustCalculation,(int)ChargedFSParticle_eta2to2p4pt05.size()));
    nBranches_->Instanton_gen_spherocity = gen_spherocity;
    nBranches_->Instanton_gen_thrust = gen_thrust;
    nBranches_->Instanton_gen_broaden = gen_broaden;
    nBranches_->Instanton_gen_transversespherocity = gen_transversespherocity;
    nBranches_->Instanton_gen_transversethrust = gen_transversethrust;
    nBranches_->Instanton_gen_transversespherocity_pTunweighted = gen_transversespherocity_pTunweighted;

    nBranches_->Instanton_gen_spherocity_eta2p4 = gen_spherocity_eta2p4;
    nBranches_->Instanton_gen_thrust_eta2p4 = gen_thrust_eta2p4;
    nBranches_->Instanton_gen_broaden_eta2p4 = gen_broaden_eta2p4;
    nBranches_->Instanton_gen_transversespherocity_eta2p4 = gen_transversespherocity_eta2p4;
    nBranches_->Instanton_gen_transversethrust_eta2p4 = gen_transversethrust_eta2p4;
    nBranches_->Instanton_gen_transversespherocity_pTunweighted_eta2p4 = gen_transversespherocity_pTunweighted_eta2p4;

    nBranches_->Instanton_gen_spherocity_eta2p4pt05 = gen_spherocity_eta2p4pt05;
    nBranches_->Instanton_gen_spherocity_eta05pt05 = gen_spherocity_eta05pt05;
    nBranches_->Instanton_gen_spherocity_eta1pt05 = gen_spherocity_eta1pt05;
    nBranches_->Instanton_gen_spherocity_eta2pt05 = gen_spherocity_eta2pt05;
    nBranches_->Instanton_gen_thrust_eta2p4pt05 = gen_thrust_eta2p4pt05;
    nBranches_->Instanton_gen_broaden_eta2p4pt05 = gen_broaden_eta2p4pt05;
    nBranches_->Instanton_gen_transversespherocity_eta2p4pt05 = gen_transversespherocity_eta2p4pt05;
    nBranches_->Instanton_gen_transversespherocity_eta05pt05 = gen_transversespherocity_eta05pt05;
    nBranches_->Instanton_gen_transversespherocity_eta1pt05 = gen_transversespherocity_eta1pt05;
    nBranches_->Instanton_gen_transversespherocity_eta2pt05 = gen_transversespherocity_eta2pt05;
    nBranches_->Instanton_gen_transversespherocity_eta05to1pt05 = gen_transversespherocity_eta05to1pt05;
    nBranches_->Instanton_gen_transversespherocity_eta1to1p5pt05 = gen_transversespherocity_eta1to1p5pt05;
    nBranches_->Instanton_gen_transversespherocity_eta1p5to2pt05 = gen_transversespherocity_eta1p5to2pt05;
    nBranches_->Instanton_gen_transversespherocity_eta2to2p4pt05 = gen_transversespherocity_eta2to2p4pt05;
    nBranches_->Instanton_gen_transversethrust_eta2p4pt05 = gen_transversethrust_eta2p4pt05;
    nBranches_->Instanton_gen_transversethrust_eta05pt05 = gen_transversethrust_eta05pt05;
    nBranches_->Instanton_gen_transversethrust_eta1pt05 = gen_transversethrust_eta1pt05;
    nBranches_->Instanton_gen_transversethrust_eta2pt05 = gen_transversethrust_eta2pt05;
    nBranches_->Instanton_gen_transversethrust_eta05to1pt05 = gen_transversethrust_eta05to1pt05;
    nBranches_->Instanton_gen_transversethrust_eta1to1p5pt05 = gen_transversethrust_eta1to1p5pt05;
    nBranches_->Instanton_gen_transversethrust_eta1p5to2pt05 = gen_transversethrust_eta1p5to2pt05;
    nBranches_->Instanton_gen_transversethrust_eta2to2p4pt05 = gen_transversethrust_eta2to2p4pt05;

    nBranches_->Instanton_gen_transversespherocity_pTunweighted_eta2p4pt05 = gen_transversespherocity_pTunweighted_eta2p4pt05;

    nBranches_->Instanton_gen_spherocity_eta2p4pt035 = gen_spherocity_eta2p4pt035;
    nBranches_->Instanton_gen_thrust_eta2p4pt035 = gen_thrust_eta2p4pt035;
    nBranches_->Instanton_gen_broaden_eta2p4pt035 = gen_broaden_eta2p4pt035;
    nBranches_->Instanton_gen_transversespherocity_eta2p4pt035 = gen_transversespherocity_eta2p4pt035;
    nBranches_->Instanton_gen_transversethrust_eta2p4pt035 = gen_transversethrust_eta2p4pt035;
    nBranches_->Instanton_gen_transversespherocity_pTunweighted_eta2p4pt035 = gen_transversespherocity_pTunweighted_eta2p4pt035;

    Float_t gen_rapiditygap_eta2p4pt05,gen_xix_eta2p4pt05,gen_xiy_eta2p4pt05;
    calc_rapidity(ChargedFSParticle_eta2p4pt05, gen_rapiditygap_eta2p4pt05,gen_xix_eta2p4pt05,gen_xiy_eta2p4pt05);
    nBranches_-> Instanton_gen_rapiditygap_eta2p4pt05 =  gen_rapiditygap_eta2p4pt05;
    nBranches_-> Instanton_gen_xix_eta2p4pt05 =  gen_xix_eta2p4pt05;
    nBranches_-> Instanton_gen_xiy_eta2p4pt05 =  gen_xiy_eta2p4pt05;

    Float_t gen_rapiditygap_eta2p4pt035,gen_xix_eta2p4pt035,gen_xiy_eta2p4pt035;
    calc_rapidity(ChargedFSParticle_eta2p4pt035, gen_rapiditygap_eta2p4pt035,gen_xix_eta2p4pt035,gen_xiy_eta2p4pt035);
    nBranches_-> Instanton_gen_rapiditygap_eta2p4pt035 =  gen_rapiditygap_eta2p4pt035;
    nBranches_-> Instanton_gen_xix_eta2p4pt035 =  gen_xix_eta2p4pt035;
    nBranches_-> Instanton_gen_xiy_eta2p4pt035 =  gen_xiy_eta2p4pt035;

    TLorentzVector gen_eta2p4pt05;
    for(size_t i=0; i<ChargedFSParticle_eta2p4pt05.size(); ++i){
      TLorentzVector momenta_p4(ChargedFSParticle_eta2p4pt05[i].px(),ChargedFSParticle_eta2p4pt05[i].py(),ChargedFSParticle_eta2p4pt05[i].pz(),ChargedFSParticle_eta2p4pt05[i].energy());
      gen_eta2p4pt05+=momenta_p4;
    }
    TLorentzVector gen_eta2p4pt035;
    for(size_t i=0; i<ChargedFSParticle_eta2p4pt035.size(); ++i){
      TLorentzVector momenta_p4(ChargedFSParticle_eta2p4pt035[i].px(),ChargedFSParticle_eta2p4pt035[i].py(),ChargedFSParticle_eta2p4pt035[i].pz(),ChargedFSParticle_eta2p4pt035[i].energy());
      gen_eta2p4pt035+=momenta_p4;
    }
    nBranches_-> Instanton_gen_mass_eta2p4pt05=gen_eta2p4pt05.M();
    nBranches_-> Instanton_gen_pt_eta2p4pt05=gen_eta2p4pt05.T();
    nBranches_-> Instanton_gen_pz_eta2p4pt05=gen_eta2p4pt05.Z();
    nBranches_-> Instanton_gen_eta_eta2p4pt05=gen_eta2p4pt05.Eta();
    nBranches_-> Instanton_gen_mass_eta2p4pt035=gen_eta2p4pt035.M();
    nBranches_-> Instanton_gen_pt_eta2p4pt035=gen_eta2p4pt035.T();
    nBranches_-> Instanton_gen_pz_eta2p4pt035=gen_eta2p4pt035.Z();
    nBranches_-> Instanton_gen_eta_eta2p4pt035=gen_eta2p4pt035.Eta();
    if(runOnMCPU_){
      bool PVfound=false;
      for( unsigned p=0; p<genParticles_->size(); ++p ){
        if((*genParticles_)[p].numberOfMothers()==1&&(*genParticles_)[p].mother(0)->pdgId()==2212&&((*genParticles_)[p].pdgId()==4||(*genParticles_)[p].pdgId()==5||(*genParticles_)[p].pdgId()==-4||(*genParticles_)[p].pdgId()==-5)) nBranches_->Instanton_N_genPromptHeavyQ++;
        if(PVfound==false&&(*genParticles_)[p].numberOfMothers()==1&&(*genParticles_)[p].mother(0)->pdgId()==2212){
          PV_X = (*genParticles_)[p].vx();
          PV_Y = (*genParticles_)[p].vy();
          PV_Z = (*genParticles_)[p].vz();
          PVfound=true;
        } 
      } 
      long unsigned int len_FS=0;
      vector<const reco::Candidate*> finalstates;
//      cout<<"begin this event"<<endl;
      for( unsigned p=0; p<genParticles_->size(); ++p ){
//        if(((*genParticles_)[p].numberOfMothers()>=1&&(*genParticles_)[p].mother(0)->pdgId()==2212)||((*genParticles_)[p].numberOfMothers()>=2&&(*genParticles_)[p].mother(1)->pdgId()==2212)){cout<<"pdgId="<<(*genParticles_)[p].pdgId()<<", vx="<<(*genParticles_)[p].vx()<<", vy="<<(*genParticles_)[p].vy()<<", vz="<<(*genParticles_)[p].vz()<<endl;}
//        if((*genParticles_)[p].status()==4){cout<<"pdgId="<<(*genParticles_)[p].pdgId()<<", vx="<<(*genParticles_)[p].vx()<<", vy="<<(*genParticles_)[p].vy()<<", vz="<<(*genParticles_)[p].vz()<<endl;}
//        if((*genParticles_)[p].pdgId()==999){cout<<"Instanton vx="<<(*genParticles_)[p].vx()<<", vy="<<(*genParticles_)[p].vy()<<", vz="<<(*genParticles_)[p].vz()<<endl;}
//        if((((*genParticles_)[p].pdgId()<=5&&(*genParticles_)[p].pdgId()>=-5)||(*genParticles_)[p].pdgId()==21)&&(*genParticles_)[p].numberOfMothers()==2&&(*genParticles_)[p].mother(0)->pdgId()==21&&(*genParticles_)[p].mother(0)->numberOfDaughters()>7){cout<<"pdgId="<<(*genParticles_)[p].pdgId()<<", vx="<<(*genParticles_)[p].vx()<<", vy="<<(*genParticles_)[p].vy()<<", vz="<<(*genParticles_)[p].vz()<<", mother1 pdgId="<<(*genParticles_)[p].mother(0)->pdgId()<<" #child="<<(*genParticles_)[p].mother(0)->numberOfDaughters()<<", mother2 pdgId="<<(*genParticles_)[p].mother(1)->pdgId()<<" #child="<<(*genParticles_)[p].mother(1)->numberOfDaughters()<<endl;}
//        if((*genParticles_)[p].pdgId()==21&&(*genParticles_)[p].numberOfDaughters()>7){cout<<"vx="<<(*genParticles_)[p].vx()<<", vy="<<(*genParticles_)[p].vy()<<", vz="<<(*genParticles_)[p].vz()<<", #child="<<(*genParticles_)[p].numberOfDaughters()<<", #mothers="<<(*genParticles_)[p].numberOfMothers()<<endl;}
        if((*genParticles_)[p].numberOfMothers()==1&&(*genParticles_)[p].mother(0)->pdgId()==2212){
          FindDaughter(&(*genParticles_)[p],finalstates);
          if((*genParticles_)[p].pdgId()==4||(*genParticles_)[p].pdgId()==5||(*genParticles_)[p].pdgId()==-4||(*genParticles_)[p].pdgId()==-5){
            for (long unsigned int iFS=len_FS; iFS<finalstates.size(); iFS++){
              nBranches_->genParticle_isfromheavyq.push_back(1);
              nBranches_->genParticle_dxy.push_back(cal_dxy(finalstates[iFS],PV_X,PV_Y,PV_Z));
              nBranches_->genParticle_dz.push_back(cal_dz(finalstates[iFS],PV_X,PV_Y,PV_Z));
              nBranches_->genParticle_IP3D.push_back(cal_IP3D(finalstates[iFS],PV_X,PV_Y,PV_Z));
//              nBranches_->Instanton_N_gen_FSParticle++;
              nBranches_->genParticle_dvxy.push_back((float)sqrt(pow((*genParticles_)[p].vx()-PV_X,2)+pow((*genParticles_)[p].vy()-PV_Y,2)));
              nBranches_->genParticle_dvz.push_back((float)fabs((*genParticles_)[p].vz()-PV_Z));
              if((float)sqrt(pow((*genParticles_)[p].vx()-PV_X,2)+pow((*genParticles_)[p].vy()-PV_Y,2))>0.002)
                nBranches_->Instanton_N_gen_ChargedFS_dvxyp002++;
//              if(finalstates[iFS]->charge()!=0) nBranches_->Instanton_N_gen_ChargedFSParticle++;
              if(cal_dxy(finalstates[iFS],PV_X,PV_Y,PV_Z)>0.02) {
                nBranches_->Instanton_N_gen_Displaced++;
                if(finalstates[iFS]->charge()!=0) nBranches_->Instanton_N_gen_ChargedDisplaced++;
              }
            }
          }
          else{
            for (long unsigned int iFS=len_FS; iFS<finalstates.size(); iFS++){
              nBranches_->genParticle_isfromheavyq.push_back(0);
              nBranches_->genParticle_dxy.push_back(cal_dxy(finalstates[iFS],PV_X,PV_Y,PV_Z));
              nBranches_->genParticle_dz.push_back(cal_dz(finalstates[iFS],PV_X,PV_Y,PV_Z));
              nBranches_->genParticle_IP3D.push_back(cal_IP3D(finalstates[iFS],PV_X,PV_Y,PV_Z));
              nBranches_->Instanton_N_gen_FSParticle++;
              nBranches_->genParticle_dvxy.push_back((float)sqrt(pow((*genParticles_)[p].vx()-PV_X,2)+pow((*genParticles_)[p].vy()-PV_Y,2)));
              nBranches_->genParticle_dvz.push_back((float)fabs((*genParticles_)[p].vz()-PV_Z));
              if((float)sqrt(pow((*genParticles_)[p].vx()-PV_X,2)+pow((*genParticles_)[p].vy()-PV_Y,2))>0.002)
                nBranches_->Instanton_N_gen_ChargedFS_dvxyp002++;
              if(finalstates[iFS]->charge()!=0) nBranches_->Instanton_N_gen_ChargedFSParticle++;
              if(cal_dxy(finalstates[iFS],PV_X,PV_Y,PV_Z)>0.02) {
                nBranches_->Instanton_N_gen_Displaced++;
                if(finalstates[iFS]->charge()!=0) nBranches_->Instanton_N_gen_ChargedDisplaced++;
              }
            }
          }
          len_FS=finalstates.size();
        }
      }
    }   

  }

  event.getByToken(verticeToken_   , vertices_     );
  event.getByToken(bsToken_, beamspot_);
  event.getByToken(muonToken_   , muons_    );
  event.getByToken(triggerObjects_  , triggerObjects);
  event.getByToken(HLTtriggersToken_,HLTtriggers_);
  std::vector<pat::Muon> muoncollection;
  muoncollection.clear();


//cout<<"vertices_ size "<<vertices_->size()<<endl;


  event.getByToken( packedpfcandidatesToken_               , packedpfcandidates_      );
  event.getByToken( ctppsProton_single_rpToken_, ctppsProton_single_rp_);
  event.getByToken( ctppsProton_multi_rpToken_, ctppsProton_multi_rp_);
  event.getByToken(jetInputToken_      , jets_    );
  event.getByToken(svToken_, svs_);
  std::vector<pfcand_struct> pfcands;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);

//cout<<"N sv "<<svs_->size()<<endl;
//for (size_t i=0; i<svs_->size(); i++) cout<<i<<", N mother "<<(*svs_)[i].numberOfMothers()<<endl;
  nBranches_->Instanton_N_SV_default=(int)svs_->size();
//  std::vector<int> alltracks_idx;
// Count number of tracks and displaced tracks (transverse impact paramter > 0.02 cm)
  int N_Trk_total=0;
  int N_Trk_Displaced_total=0;
  int N_Trk_goodDisplaced_total=0;
  int PV_N_good_hasgoodMuon = 0;
  int PV_N_good_hasgoodPFMuon = 0;
//  int N_Jet_total=0;
//  int N_Jet_TrkCut_total=0;
//  int N_Jet_TrkCut_pt05_total=0;
//  int N_Jet_TrkCut_pt1_total=0;
  std::vector<int> N_Trk(vertices_->size(),0);
  std::vector<int> N_Trk_highPurity(vertices_->size(),0);
  std::vector<int> N_Trk_highPurity_pt05(vertices_->size(),0);
  std::vector<int> N_Trk_highPurity_pt05_lowIP(vertices_->size(),0);
  std::vector<int> N_Trk_highPurity_pt05_matched(vertices_->size(),0);
  std::vector<int> N_Trk_highPurity_pt08(vertices_->size(),0);
  std::vector<int> N_Trk_highPurity_pt035(vertices_->size(),0);
  std::vector<int> N_Trk_highPurity_pt035_lowIP(vertices_->size(),0);
  std::vector<int> N_Trk_highPurity_eta05pt05(vertices_->size(),0);
  std::vector<int> N_Trk_highPurity_eta1pt05(vertices_->size(),0);
  std::vector<int> N_Trk_highPurity_eta1p5pt05(vertices_->size(),0);
  std::vector<int> N_Trk_highPurity_eta2pt05(vertices_->size(),0);
  std::vector<int> N_Trk_highPurity_eta05to1pt05(vertices_->size(),0);
  std::vector<int> N_Trk_highPurity_eta1to1p5pt05(vertices_->size(),0);
  std::vector<int> N_Trk_highPurity_eta1p5to2pt05(vertices_->size(),0);
  std::vector<int> N_Trk_highPurity_eta2to2p4pt05(vertices_->size(),0);
//  std::vector<int> N_Trk_PVAssociationQuality0(vertices_->size(),0);
//  std::vector<int> N_Trk_PVAssociationQuality1(vertices_->size(),0);
//  std::vector<int> N_Trk_PVAssociationQuality4(vertices_->size(),0);
//  std::vector<int> N_Trk_PVAssociationQuality5(vertices_->size(),0);
//  std::vector<int> N_Trk_PVAssociationQuality6(vertices_->size(),0);
//  std::vector<int> N_Trk_PVAssociationQuality7(vertices_->size(),0);
  std::vector<int> N_Trk_PVAssociationQualityLeq4_highPurity(vertices_->size(),0);
  std::vector<int> N_Trk_PVAssociationQualityLeq4_highPurity_pt05(vertices_->size(),0);
  std::vector<int> N_Trk_PVAssociationQualityLeq4_highPurity_pt08(vertices_->size(),0);
  std::vector<int> N_Trk_Displaced(vertices_->size(),0);
  std::vector<int> N_Trk_Displaced_highPurity(vertices_->size(),0);
  std::vector<int> N_Trk_Displaced_highPurity_pt05(vertices_->size(),0);
  std::vector<int> N_Trk_Displaced_highPurity_pt08(vertices_->size(),0);
//  std::vector<int> N_Trk_Displaced_PVAssociationQuality0(vertices_->size(),0);
//  std::vector<int> N_Trk_Displaced_PVAssociationQuality1(vertices_->size(),0);
//  std::vector<int> N_Trk_Displaced_PVAssociationQuality4(vertices_->size(),0);
//  std::vector<int> N_Trk_Displaced_PVAssociationQuality5(vertices_->size(),0);
//  std::vector<int> N_Trk_Displaced_PVAssociationQuality6(vertices_->size(),0);
//  std::vector<int> N_Trk_Displaced_PVAssociationQuality7(vertices_->size(),0);
  std::vector<int> N_Trk_Displaced_PVAssociationQualityLeq4_highPurity(vertices_->size(),0);
  std::vector<int> N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt05(vertices_->size(),0);
  std::vector<int> N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt08(vertices_->size(),0);
//  std::vector<int> N_Trk_sig2Displaced_highPurity(vertices_->size(),0);
//  std::vector<int> N_Trk_sig3Displaced_highPurity(vertices_->size(),0);
//  std::vector<int> N_Trk_sig5Displaced_highPurity(vertices_->size(),0);
//  std::vector<int> N_Trk_sig2Displaced_PVAssociationQualityLeq4_highPurity(vertices_->size(),0);
//  std::vector<int> N_Trk_sig3Displaced_PVAssociationQualityLeq4_highPurity(vertices_->size(),0);
//  std::vector<int> N_Trk_sig5Displaced_PVAssociationQualityLeq4_highPurity(vertices_->size(),0);

//  std::vector<int> N_Trk_fromPV2_highPurity(vertices_->size(),0);
//  std::vector<int> N_Trk_Displaced_fromPV2_highPurity(vertices_->size(),0);
//  std::vector<int> N_Trk_Displaced_fromPV2_highPurity_pt05(vertices_->size(),0);
//  std::vector<int> N_Trk_Displaced_fromPV2_highPurity_pt08(vertices_->size(),0);
//  std::vector<int> N_Trk_Displaced_fromPV2_highPurity_pt1(vertices_->size(),0);
//  std::vector<int> N_Trk_goodDisplaced_fromPV2_highPurity_sig3(vertices_->size(),0);
//  std::vector<int> N_Trk_goodDisplaced_fromPV2_highPurity_sig3_pt05(vertices_->size(),0);
//  std::vector<int> N_Trk_goodDisplaced_fromPV2_highPurity_sig3_pt08(vertices_->size(),0);
//  std::vector<int> N_Trk_goodDisplaced_fromPV2_highPurity_sig3_pt1(vertices_->size(),0);
  std::vector<int> N_Trk_formSV_sig2_pt05(vertices_->size(),0);
  std::vector<int> N_Trk_formSV_sig3_pt05(vertices_->size(),0);
  std::vector<int> N_Trk_formSV_sig5_pt05(vertices_->size(),0);
  std::vector<int> N_Trk_formSV_dxyp02_pt05(vertices_->size(),0);
  std::vector<int> N_Trk_PVAQLeq4_formSV_dxyp02_pt05(vertices_->size(),0);
  std::vector<int> N_SV_sig2_pt05(vertices_->size(),0);
  std::vector<int> N_SV_sig3_pt05(vertices_->size(),0);
  std::vector<int> N_SV_sig5_pt05(vertices_->size(),0);
  std::vector<int> N_Trk_DisplaceVp02_pt05(vertices_->size(),0);
  std::vector<int> N_Trk_formSV_dxyp02_sig2_pt05(vertices_->size(),0);
  std::vector<int> N_Trk_goodDisplaced(vertices_->size(),0);
  std::vector<int> N_Trk_goodDisplaced_highPurity(vertices_->size(),0);
  std::vector<int> N_Trk_goodDisplaced_highPurity_pt05(vertices_->size(),0);
  std::vector<int> N_Trk_goodDisplaced_highPurity_pt08(vertices_->size(),0);
//  std::vector<int> N_Trk_goodDisplaced_PVAssociationQuality0(vertices_->size(),0);
//  std::vector<int> N_Trk_goodDisplaced_PVAssociationQuality1(vertices_->size(),0);
//  std::vector<int> N_Trk_goodDisplaced_PVAssociationQuality4(vertices_->size(),0);
//  std::vector<int> N_Trk_goodDisplaced_PVAssociationQuality5(vertices_->size(),0);
//  std::vector<int> N_Trk_goodDisplaced_PVAssociationQuality6(vertices_->size(),0);
//  std::vector<int> N_Trk_goodDisplaced_PVAssociationQuality7(vertices_->size(),0);
  std::vector<int> N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity(vertices_->size(),0);
  std::vector<int> N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt05(vertices_->size(),0);
  std::vector<int> N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt08(vertices_->size(),0);
  std::vector<float> Trk_masssum(vertices_->size(),0);
  std::vector<float> Trk_ptsum(vertices_->size(),0);
  std::vector<float> Trk_ptmiss(vertices_->size(),0);
  std::vector<float> Trk_etasum(vertices_->size(),0);
  std::vector<float> Trk_pzsum(vertices_->size(),0);
  std::vector<float> Trk_TrkCut_mass(vertices_->size(),0);
  std::vector<float> Trk_TrkCut_pt05_mass(vertices_->size(),0);
  std::vector<float> Trk_TrkCut_pt(vertices_->size(),0);
  std::vector<float> Trk_TrkCut_pt05_pt(vertices_->size(),0);
  std::vector<float> Trk_TrkCut_ptmiss(vertices_->size(),0);
  std::vector<float> Trk_TrkCut_pt05_ptmiss(vertices_->size(),0);
  std::vector<float> Trk_TrkCut_etasum(vertices_->size(),0);
  std::vector<float> Trk_TrkCut_pt05_etasum(vertices_->size(),0);
  std::vector<float> Trk_TrkCut_pzsum(vertices_->size(),0);
  std::vector<float> Trk_TrkCut_pt05_pzsum(vertices_->size(),0);
  std::vector<float> Trk_Purity_pt05_mass(vertices_->size(),0);
  std::vector<float> Trk_Purity_pt05_pt(vertices_->size(),0);
  std::vector<float> Trk_Purity_pt05_ptmiss(vertices_->size(),0);
  std::vector<float> Trk_Purity_pt05_etasum(vertices_->size(),0);
  std::vector<float> Trk_Purity_pt05_pzsum(vertices_->size(),0);
  std::vector<float> Trk_Purity_pt05_lowIP_mass(vertices_->size(),0);
  std::vector<float> Trk_Purity_pt05_lowIP_pt(vertices_->size(),0);
  std::vector<float> Trk_Purity_pt05_lowIP_ptmiss(vertices_->size(),0);
  std::vector<float> Trk_Purity_pt05_lowIP_etasum(vertices_->size(),0);
  std::vector<float> Trk_Purity_pt05_lowIP_pzsum(vertices_->size(),0);
  std::vector<float> Trk_Purity_pt035_mass(vertices_->size(),0);
  std::vector<float> Trk_Purity_pt035_pt(vertices_->size(),0);
  std::vector<float> Trk_Purity_pt035_ptmiss(vertices_->size(),0);
  std::vector<float> Trk_Purity_pt035_etasum(vertices_->size(),0);
  std::vector<float> Trk_Purity_pt035_pzsum(vertices_->size(),0);
  std::vector<float> Trk_Purity_pt035_lowIP_mass(vertices_->size(),0);
  std::vector<float> Trk_Purity_pt035_lowIP_pt(vertices_->size(),0);
  std::vector<float> Trk_Purity_pt035_lowIP_ptmiss(vertices_->size(),0);
  std::vector<float> Trk_Purity_pt035_lowIP_etasum(vertices_->size(),0);
  std::vector<float> Trk_Purity_pt035_lowIP_pzsum(vertices_->size(),0);
  std::vector<int> N_Jet(vertices_->size(),0);
//  std::vector<int> N_Jet_TrkCut(vertices_->size(),0);
//  std::vector<int> N_Jet_TrkCut_pt05(vertices_->size(),0);
//  std::vector<int> N_Jet_TrkCut_pt05_tjpt10(vertices_->size(),0);
//  std::vector<int> N_Jet_TrkCut_pt1(vertices_->size(),0);
  std::vector<int> vtx_N_goodMuon(vertices_->size(),0);
  std::vector<int> vtx_N_goodPFMuon(vertices_->size(),0);
  std::vector<bool> vtx_isBPHtrigger_fromMuon(vertices_->size(),false); 
  std::vector<bool> vtx_isBPHtrigger_fromPFMuon(vertices_->size(),false);
  std::vector<float> vtx_goodMuonIP(vertices_->size(),999.);
  std::vector<float> vtx_goodPFMuon_PVDistance(vertices_->size(),999.);
  std::vector<float> vtx_genvertex_PVDistance(vertices_->size(),999.);
  std::vector<std::vector<reco::TransientTrack>> alltracks(vertices_->size());
  std::vector<std::vector<pat::PackedCandidate>> alltracks_pf(vertices_->size());
  std::vector<std::vector<reco::TransientTrack>> alltracks_purity_pt05(vertices_->size());
  std::vector<std::vector<pat::PackedCandidate>> alltracks_pf_purity_pt05(vertices_->size());
  std::vector<std::vector<reco::TransientTrack>> alltracks_purity_pt05_lowIP(vertices_->size());
  std::vector<std::vector<pat::PackedCandidate>> alltracks_pf_purity_pt05_lowIP(vertices_->size());
  std::vector<std::vector<reco::TransientTrack>> alltracks_purity_eta05pt05(vertices_->size());
  std::vector<std::vector<pat::PackedCandidate>> alltracks_pf_purity_eta05pt05(vertices_->size());
  std::vector<std::vector<reco::TransientTrack>> alltracks_purity_eta1pt05(vertices_->size());
  std::vector<std::vector<pat::PackedCandidate>> alltracks_pf_purity_eta1pt05(vertices_->size());
  std::vector<std::vector<reco::TransientTrack>> alltracks_purity_eta1p5pt05(vertices_->size());
  std::vector<std::vector<pat::PackedCandidate>> alltracks_pf_purity_eta1p5pt05(vertices_->size());
  std::vector<std::vector<reco::TransientTrack>> alltracks_purity_eta2pt05(vertices_->size());
  std::vector<std::vector<pat::PackedCandidate>> alltracks_pf_purity_eta2pt05(vertices_->size());
  std::vector<std::vector<reco::TransientTrack>> alltracks_purity_eta05to1pt05(vertices_->size());
  std::vector<std::vector<pat::PackedCandidate>> alltracks_pf_purity_eta05to1pt05(vertices_->size());
  std::vector<std::vector<reco::TransientTrack>> alltracks_purity_eta1to1p5pt05(vertices_->size());
  std::vector<std::vector<pat::PackedCandidate>> alltracks_pf_purity_eta1to1p5pt05(vertices_->size());
  std::vector<std::vector<reco::TransientTrack>> alltracks_purity_eta1p5to2pt05(vertices_->size());
  std::vector<std::vector<pat::PackedCandidate>> alltracks_pf_purity_eta1p5to2pt05(vertices_->size());
  std::vector<std::vector<reco::TransientTrack>> alltracks_purity_eta2to2p4pt05(vertices_->size());
  std::vector<std::vector<pat::PackedCandidate>> alltracks_pf_purity_eta2to2p4pt05(vertices_->size());
  std::vector<std::vector<reco::TransientTrack>> alltracks_TrkCut(vertices_->size());
  std::vector<std::vector<pat::PackedCandidate>> alltracks_pf_TrkCut(vertices_->size());
  std::vector<std::vector<reco::TransientTrack>> alltracks_TrkCut_pt05(vertices_->size());
  std::vector<std::vector<pat::PackedCandidate>> alltracks_pf_TrkCut_pt05(vertices_->size());
  std::vector<std::vector<reco::TransientTrack>> alltracks_purity_pt035(vertices_->size());
  std::vector<std::vector<pat::PackedCandidate>> alltracks_pf_purity_pt035(vertices_->size());
  std::vector<std::vector<reco::TransientTrack>> alltracks_purity_pt035_lowIP(vertices_->size());
  std::vector<std::vector<pat::PackedCandidate>> alltracks_pf_purity_pt035_lowIP(vertices_->size());
  std::vector<std::vector<float>> Trk_px(vertices_->size());
  std::vector<std::vector<float>> Trk_py(vertices_->size());
  std::vector<std::vector<float>> Trk_pz(vertices_->size());
  std::vector<std::vector<float>> Trk_p0(vertices_->size());
  std::vector<std::vector<float>> Trk_pt(vertices_->size());
  std::vector<std::vector<float>> Trk_mass(vertices_->size());
  std::vector<std::vector<float>> Trk_eta(vertices_->size());
  std::vector<std::vector<float>> Trk_phi(vertices_->size());
  std::vector<std::vector<int>>   Trk_PVAssociationQuality(vertices_->size());
  std::vector<std::vector<int>>  Trk_isHighPurity(vertices_->size());
  std::vector<std::vector<float>> Trk_vx(vertices_->size());
  std::vector<std::vector<float>> Trk_vy(vertices_->size());
  std::vector<std::vector<float>> Trk_vz(vertices_->size());
  std::vector<std::vector<float>> Trk_IP(vertices_->size());
  std::vector<std::vector<float>> Trk_IPsignificance(vertices_->size());
  std::vector<std::vector<float>> Trk_z0(vertices_->size());
  std::vector<std::vector<float>> Trk_z0significance(vertices_->size());
  std::vector<std::vector<int>> Trk_matched(vertices_->size());
  std::vector<std::vector<int>>   Trk_pdgId(vertices_->size());
  std::vector<std::vector<int>>   Trk_charge(vertices_->size());
  std::vector<std::vector<float>> TrackJet_px(vertices_->size());
  std::vector<std::vector<float>> TrackJet_py(vertices_->size());
  std::vector<std::vector<float>> TrackJet_pt(vertices_->size());
  std::vector<std::vector<float>> TrackJet_pz(vertices_->size());
  std::vector<std::vector<float>> TrackJet_E(vertices_->size());
  std::vector<std::vector<float>> TrackJet_mass(vertices_->size());
  std::vector<std::vector<float>> TrackJet_phi(vertices_->size());
  std::vector<std::vector<float>> TrackJet_eta(vertices_->size());
//  for(reco::VertexCollection::const_iterator vtx = vertices_->begin(); vtx != vertices_->end(); ++vtx){
//    cout<<"PV "<<vtx->position().X()<<" "<<vtx->position().Y()<<" "<<vtx->position().Z()<<endl;
//  }

//  int ijet = 0;
  int njets_pt20=0;
  int njets_pt20_MVAPU05=0;
 // cout<<"test jet"<<endl;
  for (const pat::Jet &j : *jets_) {
    if (j.pt() < 20 || j.eta()>2.5 || j.eta()<-2.5) continue;
    njets_pt20++;
    if(j.userFloat("pileupJetId:fullDiscriminant")>0.5) njets_pt20_MVAPU05++;
//    printf("jet  with pt %5.1f (raw pt %5.1f), eta %+4.2f, btag CSV %.3f, CISV %.3f, pileup mva disc %+.2f\n",j.pt(), j.pt()*j.jecFactor("Uncorrected"), j.eta(), std::max(0.f,j.bDiscriminator("combinedSecondaryVertexBJetTags")), std::max(0.f,j.bDiscriminator("combinedInclusiveSecondaryVertexBJetTags")), j.userFloat("pileupJetId:fullDiscriminant"));
//    if ((++ijet) == 1) {
//      std::vector daus(j.daughterPtrVector());
//      std::sort(daus.begin(), daus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2->pt(); }); // the joys of C++11
//      for (unsigned int i2 = 0, n = daus.size(); i2 < n; ++i2) {
//        const pat::PackedCandidate &cand = dynamic_cast<const pat::PackedCandidate &>(*daus[i2]);
//        printf("         constituent %3d: vtxid %d,  pt %6.2f, dz(pv) %+.3f, pdgId %+3d\n", i2,cand.vertexRef().key(),cand.pt(),cand.dz(),cand.pdgId());
//      }
//    }

  }
  nBranches_->Instanton_N_PFAK4Jet_pt20_total = njets_pt20;
  nBranches_->Instanton_N_PFAK4Jet_pt20_PUMVA05_total = njets_pt20_MVAPU05;

  std::vector<pat::PackedCandidate> alltracks_pfmuon;
  std::string period = "";
  if(runStrings_["run"].find("18") != std::string::npos) period = "2018";
  else if(runStrings_["run"].find("17") != std::string::npos) period = "2017_postTS2";
  else if(runStrings_["run"].find("16") != std::string::npos) period = "2016_postTS2";
  edm::ESHandle<LHCInfo> hLHCInfo;
  std::string lhcInfoLabel("");
  iSetup.get<LHCInfoRcd>().get(lhcInfoLabel, hLHCInfo);
  double beamXangle=hLHCInfo->crossingAngle();
  PPSEff *PPS_eff_ = new PPSEff();
  //for( size_t ii = 0; ii < ctppsProton_single_rp_->size(); ++ii ){
  //  reco::ForwardProton fp = (*ctppsProton_single_rp_)[ii];
  //  cout<<"Forwand proton single rp "<<ii<<", pt="<<fp.pt()<<", vx="<<fp.vx()<<", vy="<<fp.vy()<<", vz="<<fp.vz()<<", chi2="<<fp.chi2()<<", ndof="<<fp.ndof()<<", thetaX="<<fp.thetaX()<<", thetaY="<<fp.thetaY()<<endl;
  //}
  //cout<<"doForwardProton_"<<doForwardProton_<<endl;
  if(doForwardProton_){
    int N_ForwardProton = 0;
    for( size_t ii = 0; ii < ctppsProton_multi_rp_->size(); ++ii ){
      reco::ForwardProton fp = (*ctppsProton_multi_rp_)[ii];
      //cout<<"Forwand proton multi rp "<<ii<<", pt="<<fp.pt()<<", vx="<<fp.vx()<<", vy="<<fp.vy()<<", vz="<<fp.vz()<<", chi2="<<fp.chi2()<<", ndof="<<fp.ndof()<<", thetaX="<<fp.thetaX()<<", thetaY="<<fp.thetaY()<<", xi="<<fp.xi()<<endl;
      //cout<<"validFit="<<fp.validFit()<<endl;
      //cout<<"len "<<fp.contributingLocalTracks().size()<<endl;
      //cout<<(fp.contributingLocalTracks()[0])->x()<<", "<<(fp.contributingLocalTracks()[0])->y()<<endl;
      if(!fp.validFit()) continue;
      CTPPSDetId detid( (*(fp.contributingLocalTracks().begin()))->rpId() );
      //cout<<"period="<<period<<", arm="<<std::to_string(detid.arm())<<endl;
      if(!IsWithinAperture(period, "arm"+std::to_string(detid.arm()), beamXangle, fp.xi(), fp.thetaX())) continue;
      bool isPixel = (detid.station()==2);
      bool VetoTrack = false;
      if (isPixel){
        for(int itrk=0; itrk<(int)fp.contributingLocalTracks().size(); itrk++){
          CTPPSDetId rpId(fp.contributingLocalTracks()[itrk]->rpId());
          if(rpId.station()!=2) continue;
          if(!PPS_eff_->InFiducial(detid.arm(), 2, runOnMC_?319313:event.id().run(), fp.contributingLocalTracks()[itrk]->x(), fp.contributingLocalTracks()[itrk]->y())) VetoTrack = true;
        }  
      }
      if(!VetoTrack) {
        N_ForwardProton ++;
        nBranches_->ForwardProton_vx.push_back((float)fp.vx());
        nBranches_->ForwardProton_vxError.push_back((float)fp.vxError());
        nBranches_->ForwardProton_vy.push_back((float)fp.vy());
        nBranches_->ForwardProton_vyError.push_back((float)fp.vyError());
        nBranches_->ForwardProton_vz.push_back((float)fp.vz());
        nBranches_->ForwardProton_px.push_back((float)fp.px());
        nBranches_->ForwardProton_py.push_back((float)fp.py());
        nBranches_->ForwardProton_pz.push_back((float)fp.pz());
        nBranches_->ForwardProton_pt.push_back((float)fp.pt());
        nBranches_->ForwardProton_chi2.push_back((float)fp.chi2());
        nBranches_->ForwardProton_ndof.push_back((float)fp.ndof());
        nBranches_->ForwardProton_xi.push_back((float)fp.xi());
        nBranches_->ForwardProton_xiError.push_back((float)fp.xiError());
        nBranches_->ForwardProton_thetaXError.push_back((float)fp.thetaXError());
        nBranches_->ForwardProton_thetaX.push_back((float)fp.thetaX());
        nBranches_->ForwardProton_thetaY.push_back((float)fp.thetaY());
        nBranches_->ForwardProton_thetaYError.push_back((float)fp.thetaYError());
        double pps_eff_value = forwardProtonEff_->getEfficiency(fp, runOnMC_?319313:event.id().run());
        double pps_eff_statunc = forwardProtonEff_->getEfficiencyunc_stat(fp,runOnMC_?319313:event.id().run());
        double pps_eff_showerineff = forwardProtonEff_->getEfficiency_showerineff(fp);
        nBranches_->ForwardProton_eff.push_back((float)pps_eff_value);
        nBranches_->ForwardProton_eff_statunc.push_back((float)pps_eff_statunc);
        nBranches_->ForwardProton_showerineff.push_back((float)pps_eff_showerineff);
      }
    }
    nBranches_->Instanton_N_ForwardProton = N_ForwardProton;
  }
  double calo_sum_eT_HFp=0;
  double calo_sum_eT_HFm=0;
  double calo_sum_energy_HFp=0;
  double calo_sum_energy_HFm=0;
  double calo_leading_energy_HFp=0;
  double calo_leading_energy_HFm=0;
  for( size_t ii = 0; ii < packedpfcandidates_->size(); ++ii ){
    pat::PackedCandidate pf = (*packedpfcandidates_)[ii];
    double eta = pf.eta();
    CaloType subdetHad = GetTowerSubdetHad(eta);
    //if((subdetHad == kHFp)||(subdetHad == kHFm)){
    //  cout<<"Find PF Candidate in HF region"<<endl;
    //  cout<<"eta="<<pf.eta()<<", pdgId="<<pf.pdgId()<<", energy="<<pf.energy()<<", et="<<pf.et()<<", pt="<<pf.pt()<<", hcalFraction="<<pf.hcalFraction()<<", caloFraction="<<pf.caloFraction()<<endl;
    //  cout<<"rawHcalFraction="<<pf.rawHcalFraction()<<", rawCaloFraction="<<pf.rawCaloFraction()<<endl;
    //}
    if(subdetHad == kHFp){
      calo_sum_eT_HFp += pf.et();
      calo_sum_energy_HFp += pf.energy();
      if(pf.energy()>calo_leading_energy_HFp) calo_leading_energy_HFp=pf.energy();
    }
    if(subdetHad == kHFm){
      calo_sum_eT_HFm += pf.et();
      calo_sum_energy_HFm += pf.energy();
      if(pf.energy()>calo_leading_energy_HFm) calo_leading_energy_HFm=pf.energy();
    }
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
//      bool isgoodDisplaced_sig3=false;
      bool ishighPurity=false;
      bool isLowIP = false;
      bool ishighPurity_eta05pt05=false;
      bool ishighPurity_eta1pt05=false;
      bool ishighPurity_eta1p5pt05=false;
      bool ishighPurity_eta2pt05=false;
      bool ishighPurity_pt05=false;
      bool ishighPurity_pt08=false;
      bool ishighPurity_pt035=false;
//      bool issig2Displaced=false;
//      bool issig3Displaced=false;
//      bool issig5Displaced=false;
      int PVAssociationQuality = -1;
      if (trackrandomdrop_){ 
         double random_drop=(double)rand()/(double)RAND_MAX;
         if ((pf.pt()<=20&&random_drop>=0.979)||(pf.pt()>20&&random_drop>=0.99)) continue;
      }
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
        if(pf.pt()>0.5&&pf.eta()<2.4&&pf.eta()>-2.4) {
            ishighPurity_pt05=true;
            N_Trk_highPurity_pt05[vtxindex]++;
          if(pf.eta()<2&&pf.eta()>-2){
            ishighPurity_eta2pt05=true; 
            N_Trk_highPurity_eta2pt05[vtxindex]++;
            if(pf.eta()<1.5&&pf.eta()>-1.5){
              ishighPurity_eta1p5pt05=true;
              N_Trk_highPurity_eta1p5pt05[vtxindex]++;
              if(pf.eta()<1&&pf.eta()>-1){
                ishighPurity_eta1pt05=true;
                N_Trk_highPurity_eta1pt05[vtxindex]++;
                if(pf.eta()<0.5&&pf.eta()>-0.5){
                  ishighPurity_eta05pt05=true;
                  N_Trk_highPurity_eta05pt05[vtxindex]++;
                }
                else N_Trk_highPurity_eta05to1pt05[vtxindex]++;
              }
              else N_Trk_highPurity_eta1to1p5pt05[vtxindex]++;
            }
            else N_Trk_highPurity_eta1p5to2pt05[vtxindex]++;
          }
          else N_Trk_highPurity_eta2to2p4pt05[vtxindex]++;
        }
        if(pf.pt()>0.8) {ishighPurity_pt08=true; N_Trk_highPurity_pt08[vtxindex]++;}
        if(pf.pt()>0.35&&pf.eta()<2.4&&pf.eta()>-2.4) {ishighPurity_pt035=true;  N_Trk_highPurity_pt035[vtxindex]++;}

      }
      reco::TransientTrack  _track_ = (*builder).build(pf.pseudoTrack());
      _track_.setBeamSpot(*beamspot_);
      GlobalPoint vert(pf.vertexRef()->x(), pf.vertexRef()->y(), pf.vertexRef()->z());
      TrajectoryStateClosestToPoint  traj = _track_.trajectoryStateClosestToPoint(vert);
      if(fabs((float)traj.perigeeParameters().transverseImpactParameter()/traj.perigeeError().transverseImpactParameterError())<8.0&&fabs((float)traj.perigeeParameters().longitudinalImpactParameter()/traj.perigeeError().longitudinalImpactParameterError())<13.0) isLowIP=true;
      if(ishighPurity_pt05&&isLowIP) N_Trk_highPurity_pt05_lowIP[vtxindex]++;
      if(ishighPurity_pt035&&isLowIP) N_Trk_highPurity_pt035_lowIP[vtxindex]++;
//      if(traj.perigeeParameters().transverseImpactParameter()/traj.perigeeError().transverseImpactParameterError()>2.){
//        issig2Displaced=true;
//        N_Trk_sig2Displaced_highPurity[vtxindex]+=(int)ishighPurity;
//      }
//      if(traj.perigeeParameters().transverseImpactParameter()/traj.perigeeError().transverseImpactParameterError()>3.){
//        issig3Displaced=true;
//        N_Trk_sig3Displaced_highPurity[vtxindex]+=(int)ishighPurity;
//      }
//      if(traj.perigeeParameters().transverseImpactParameter()/traj.perigeeError().transverseImpactParameterError()>5.){
//        issig5Displaced=true;
//        N_Trk_sig5Displaced_highPurity[vtxindex]+=(int)ishighPurity;
//      }

      if(traj.perigeeParameters().transverseImpactParameter() > 0.02){
//          math::XYZPoint vert(pf.vertexRef()->x(), pf.vertexRef()->y(), pf.vertexRef()->z());
//          if(pf.dxy(vert)> 0.002){
        N_Trk_Displaced[vtxindex]++;
        N_Trk_Displaced_highPurity[vtxindex]+=ishighPurity;
        N_Trk_Displaced_highPurity_pt05[vtxindex]+=ishighPurity_pt05;
        N_Trk_Displaced_highPurity_pt08[vtxindex]+=ishighPurity_pt08;
        N_Trk_Displaced_total++;
        isDisplaced=true;
        if(traj.perigeeParameters().transverseImpactParameter()/traj.perigeeError().transverseImpactParameterError()>5.){
          N_Trk_goodDisplaced[vtxindex]++;
          N_Trk_goodDisplaced_highPurity[vtxindex]+=ishighPurity;
          N_Trk_goodDisplaced_highPurity_pt05[vtxindex]+=ishighPurity_pt05;
          N_Trk_goodDisplaced_highPurity_pt08[vtxindex]+=ishighPurity_pt08;
          N_Trk_goodDisplaced_total++;
          isgoodDisplaced=true;
        }
//        if(traj.perigeeParameters().transverseImpactParameter()/traj.perigeeError().transverseImpactParameterError()>3.) {
//          isgoodDisplaced_sig3=true;
//        }    
      }
/*      if(pf.fromPV()==2){
        N_Trk_Displaced_fromPV2_highPurity[vtxindex]+=(int)(isDisplaced&&ishighPurity);
        N_Trk_Displaced_fromPV2_highPurity_pt05[vtxindex]+=(int)(isDisplaced&&ishighPurity_pt05);
        N_Trk_Displaced_fromPV2_highPurity_pt08[vtxindex]+=(int)(isDisplaced&&ishighPurity_pt08);
        N_Trk_Displaced_fromPV2_highPurity_pt1[vtxindex]+=(int)(isDisplaced&&ishighPurity_pt1);
        N_Trk_goodDisplaced_fromPV2_highPurity_sig3[vtxindex]+=(int)(isgoodDisplaced_sig3&&ishighPurity);
        N_Trk_goodDisplaced_fromPV2_highPurity_sig3_pt05[vtxindex]+=(int)(isgoodDisplaced_sig3&&ishighPurity_pt05);
        N_Trk_goodDisplaced_fromPV2_highPurity_sig3_pt08[vtxindex]+=(int)(isgoodDisplaced_sig3&&ishighPurity_pt08);
        N_Trk_goodDisplaced_fromPV2_highPurity_sig3_pt1[vtxindex]+=(int)(isgoodDisplaced_sig3&&ishighPurity_pt1);
        N_Trk_fromPV2_highPurity[vtxindex]+=(int)ishighPurity;
      }
*/
      switch(pf.pvAssociationQuality()){
        case pat::PackedCandidate::PVAssociationQuality::NotReconstructedPrimary:
//                N_Trk_PVAssociationQuality0[vtxindex]++;
//                N_Trk_Displaced_PVAssociationQuality0[vtxindex]+=(int)isDisplaced;
//                N_Trk_goodDisplaced_PVAssociationQuality0[vtxindex]+=(int)isgoodDisplaced;
                PVAssociationQuality=0;
                break;
        case pat::PackedCandidate::PVAssociationQuality::OtherDeltaZ: 
//                N_Trk_PVAssociationQuality1[vtxindex]++;
//                N_Trk_Displaced_PVAssociationQuality1[vtxindex]+=(int)isDisplaced;
//                N_Trk_goodDisplaced_PVAssociationQuality1[vtxindex]+=(int)isgoodDisplaced;
                PVAssociationQuality=1;
                break;
        case pat::PackedCandidate::PVAssociationQuality::CompatibilityBTag: 
//                N_Trk_PVAssociationQuality4[vtxindex]++;
//                N_Trk_Displaced_PVAssociationQuality4[vtxindex]+=(int)isDisplaced;
//                N_Trk_goodDisplaced_PVAssociationQuality4[vtxindex]+=(int)isgoodDisplaced;
                N_Trk_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)ishighPurity;
                N_Trk_PVAssociationQualityLeq4_highPurity_pt05[vtxindex]+=(int)ishighPurity_pt05;
                N_Trk_PVAssociationQualityLeq4_highPurity_pt08[vtxindex]+=(int)ishighPurity_pt08;
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(ishighPurity&&isDisplaced);
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt05[vtxindex]+=(int)(ishighPurity_pt05&&isDisplaced);
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt08[vtxindex]+=(int)(ishighPurity_pt08&&isDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(ishighPurity&&isgoodDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt05[vtxindex]+=(int)(ishighPurity_pt05&&isgoodDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt08[vtxindex]+=(int)(ishighPurity_pt08&&isgoodDisplaced);
//                N_Trk_sig2Displaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(issig2Displaced&&ishighPurity);
//                N_Trk_sig3Displaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(issig3Displaced&&ishighPurity);
//                N_Trk_sig5Displaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(issig5Displaced&&ishighPurity);
                PVAssociationQuality=4;
                break;
        case pat::PackedCandidate::PVAssociationQuality::CompatibilityDz: 
//                N_Trk_PVAssociationQuality5[vtxindex]++;
//                N_Trk_Displaced_PVAssociationQuality5[vtxindex]+=(int)isDisplaced;
//                N_Trk_goodDisplaced_PVAssociationQuality5[vtxindex]+=(int)isgoodDisplaced;
                N_Trk_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)ishighPurity;
                N_Trk_PVAssociationQualityLeq4_highPurity_pt05[vtxindex]+=(int)ishighPurity_pt05;
                N_Trk_PVAssociationQualityLeq4_highPurity_pt08[vtxindex]+=(int)ishighPurity_pt08;
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(ishighPurity&&isDisplaced);
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt05[vtxindex]+=(int)(ishighPurity_pt05&&isDisplaced);
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt08[vtxindex]+=(int)(ishighPurity_pt08&&isDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(ishighPurity&&isgoodDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt05[vtxindex]+=(int)(ishighPurity_pt05&&isgoodDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt08[vtxindex]+=(int)(ishighPurity_pt08&&isgoodDisplaced);
//                N_Trk_sig2Displaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(issig2Displaced&&ishighPurity);
//                N_Trk_sig3Displaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(issig3Displaced&&ishighPurity);
//                N_Trk_sig5Displaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(issig5Displaced&&ishighPurity);
                PVAssociationQuality=5;
                break;
        case pat::PackedCandidate::PVAssociationQuality::UsedInFitLoose:
//                N_Trk_PVAssociationQuality6[vtxindex]++;
//                N_Trk_Displaced_PVAssociationQuality6[vtxindex]+=(int)isDisplaced;
//                N_Trk_goodDisplaced_PVAssociationQuality6[vtxindex]+=(int)isgoodDisplaced;
                N_Trk_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)ishighPurity;
                N_Trk_PVAssociationQualityLeq4_highPurity_pt05[vtxindex]+=(int)ishighPurity_pt05;
                N_Trk_PVAssociationQualityLeq4_highPurity_pt08[vtxindex]+=(int)ishighPurity_pt08;
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(ishighPurity&&isDisplaced);
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt05[vtxindex]+=(int)(ishighPurity_pt05&&isDisplaced);
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt08[vtxindex]+=(int)(ishighPurity_pt08&&isDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(ishighPurity&&isgoodDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt05[vtxindex]+=(int)(ishighPurity_pt05&&isgoodDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt08[vtxindex]+=(int)(ishighPurity_pt08&&isgoodDisplaced);
//                N_Trk_sig2Displaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(issig2Displaced&&ishighPurity);
//                N_Trk_sig3Displaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(issig3Displaced&&ishighPurity);
//                N_Trk_sig5Displaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(issig5Displaced&&ishighPurity);
                PVAssociationQuality=6;
                break;
        case pat::PackedCandidate::PVAssociationQuality::UsedInFitTight:
//                N_Trk_PVAssociationQuality7[vtxindex]++;
//                N_Trk_Displaced_PVAssociationQuality7[vtxindex]+=(int)isDisplaced;
//                N_Trk_goodDisplaced_PVAssociationQuality7[vtxindex]+=(int)isgoodDisplaced;
                N_Trk_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)ishighPurity;
                N_Trk_PVAssociationQualityLeq4_highPurity_pt05[vtxindex]+=(int)ishighPurity_pt05;
                N_Trk_PVAssociationQualityLeq4_highPurity_pt08[vtxindex]+=(int)ishighPurity_pt08;
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(ishighPurity&&isDisplaced);
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt05[vtxindex]+=(int)(ishighPurity_pt05&&isDisplaced);
                N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt08[vtxindex]+=(int)(ishighPurity_pt08&&isDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(ishighPurity&&isgoodDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt05[vtxindex]+=(int)(ishighPurity_pt05&&isgoodDisplaced);
                N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt08[vtxindex]+=(int)(ishighPurity_pt08&&isgoodDisplaced);
//                N_Trk_sig2Displaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(issig2Displaced&&ishighPurity);
//                N_Trk_sig3Displaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(issig3Displaced&&ishighPurity);
//                N_Trk_sig5Displaced_PVAssociationQualityLeq4_highPurity[vtxindex]+=(int)(issig5Displaced&&ishighPurity);
                PVAssociationQuality=7;
                break;
      }
      alltracks[vtxindex].push_back(_track_);
      alltracks_pf[vtxindex].push_back(pf);
      if(ishighPurity_pt05){
        alltracks_purity_pt05[vtxindex].push_back(_track_);
        alltracks_pf_purity_pt05[vtxindex].push_back(pf);
        if(isLowIP){
          alltracks_purity_pt05_lowIP[vtxindex].push_back(_track_);
          alltracks_pf_purity_pt05_lowIP[vtxindex].push_back(pf);
        }
      }
      if(ishighPurity_eta05pt05){
        alltracks_purity_eta05pt05[vtxindex].push_back(_track_);
        alltracks_pf_purity_eta05pt05[vtxindex].push_back(pf);
      }
      if(ishighPurity_eta1pt05){
        alltracks_purity_eta1pt05[vtxindex].push_back(_track_);
        alltracks_pf_purity_eta1pt05[vtxindex].push_back(pf);
      }
      if(ishighPurity_eta1p5pt05){
        alltracks_purity_eta1p5pt05[vtxindex].push_back(_track_);
        alltracks_pf_purity_eta1p5pt05[vtxindex].push_back(pf);
      }
      if(ishighPurity_eta2pt05){
        alltracks_purity_eta2pt05[vtxindex].push_back(_track_);
        alltracks_pf_purity_eta2pt05[vtxindex].push_back(pf);
      }
      if(!ishighPurity_eta05pt05&&ishighPurity_eta1pt05){
        alltracks_purity_eta05to1pt05[vtxindex].push_back(_track_);
        alltracks_pf_purity_eta05to1pt05[vtxindex].push_back(pf);
      }
      if(!ishighPurity_eta1pt05&&ishighPurity_eta1p5pt05){
        alltracks_purity_eta1to1p5pt05[vtxindex].push_back(_track_);
        alltracks_pf_purity_eta1to1p5pt05[vtxindex].push_back(pf);
      }
      if(!ishighPurity_eta1p5pt05&&ishighPurity_eta2pt05){
        alltracks_purity_eta1p5to2pt05[vtxindex].push_back(_track_);
        alltracks_pf_purity_eta1p5to2pt05[vtxindex].push_back(pf);
      }
      if(!ishighPurity_eta2pt05&&ishighPurity_pt05){
        alltracks_purity_eta2to2p4pt05[vtxindex].push_back(_track_);
        alltracks_pf_purity_eta2to2p4pt05[vtxindex].push_back(pf);
      }
      if(PVAssociationQuality>=4&&ishighPurity){ 
        alltracks_TrkCut[vtxindex].push_back(_track_);
        alltracks_pf_TrkCut[vtxindex].push_back(pf);
        if(ishighPurity_pt05){
          alltracks_TrkCut_pt05[vtxindex].push_back(_track_);
          alltracks_pf_TrkCut_pt05[vtxindex].push_back(pf);
        }
      }
      if(ishighPurity_pt035){
        alltracks_purity_pt035[vtxindex].push_back(_track_);
        alltracks_pf_purity_pt035[vtxindex].push_back(pf);
        if(isLowIP){
          alltracks_purity_pt035_lowIP[vtxindex].push_back(_track_);
          alltracks_pf_purity_pt035_lowIP[vtxindex].push_back(pf);
        }
      }
      if(doTrack_){
        Trk_px[vtxindex].push_back((float)pf.px());
        Trk_py[vtxindex].push_back((float)pf.py());
        Trk_pz[vtxindex].push_back((float)pf.pz());
        Trk_pt[vtxindex].push_back((float)pf.pt());
        Trk_p0[vtxindex].push_back((float)pf.energy());
        Trk_mass[vtxindex].push_back((float)pf.mass());
        Trk_eta[vtxindex].push_back((float)pf.eta());
        Trk_phi[vtxindex].push_back((float)pf.phi());
        Trk_PVAssociationQuality[vtxindex].push_back(PVAssociationQuality);
        Trk_isHighPurity[vtxindex].push_back((int)ishighPurity);
        Trk_vx[vtxindex].push_back((float)pf.vx());
        Trk_vy[vtxindex].push_back((float)pf.vy());
        Trk_vz[vtxindex].push_back((float)pf.vz());
        Trk_IP[vtxindex].push_back((float)traj.perigeeParameters().transverseImpactParameter());
        Trk_IPsignificance[vtxindex].push_back((float)traj.perigeeParameters().transverseImpactParameter()/traj.perigeeError().transverseImpactParameterError());
        Trk_z0[vtxindex].push_back((float)traj.perigeeParameters().longitudinalImpactParameter());
        Trk_z0significance[vtxindex].push_back((float)traj.perigeeParameters().longitudinalImpactParameter()/traj.perigeeError().longitudinalImpactParameterError());
        Trk_charge[vtxindex].push_back(pf.charge());
        Trk_pdgId[vtxindex].push_back(pf.pdgId());
        if(runOnMC_) Trk_matched[vtxindex].push_back(0);
      }
    }
  }


    nBranches_->calo_sum_energy_HFp = calo_sum_energy_HFp;
    nBranches_->calo_sum_energy_HFm = calo_sum_energy_HFm;
    nBranches_->calo_sum_eT_HFp = calo_sum_eT_HFp;
    nBranches_->calo_sum_eT_HFm = calo_sum_eT_HFm;
    nBranches_->calo_leading_energy_HFp = calo_leading_energy_HFp;
    nBranches_->calo_leading_energy_HFm = calo_leading_energy_HFm;

//Calculate total mass of the tracks for each PV
/*  for( size_t j = 0; j < vertices_->size(); j++){
    math::XYZTLorentzVector p4(0,0,0,0);
    math::XYZTLorentzVector p4_TrkCut(0,0,0,0);
    math::XYZTLorentzVector p4_TrkCut05(0,0,0,0);
    math::XYZTLorentzVector p4_TrkCut1(0,0,0,0);
    float pt_TrkCut_sum=0;
    float pt_TrkCut_pt05_sum=0;
    float pt_TrkCut_pt1_sum=0;
    for(int k = 0; k < N_Trk[j]; k++)
      p4 = p4 + alltracks_pf[j][k].p4();
    for(int k = 0; k < N_Trk_PVAssociationQualityLeq4_highPurity[j]; k++){
      p4_TrkCut = p4_TrkCut + alltracks_pf_TrkCut[j][k].p4();
      pt_TrkCut_sum += alltracks_pf_TrkCut[j][k].pt();
    }
    for(int k = 0; k < N_Trk_PVAssociationQualityLeq4_highPurity_pt05[j]; k++){
      p4_TrkCut05 = p4_TrkCut05 + alltracks_pf_TrkCut_pt05[j][k].p4();
      pt_TrkCut_pt05_sum += alltracks_pf_TrkCut_pt05[j][k].pt();
    }
    for(int k = 0; k < N_Trk_PVAssociationQualityLeq4_highPurity_pt1[j]; k++){
      p4_TrkCut1 = p4_TrkCut1 + alltracks_pf_TrkCut_pt1[j][k].p4();
      pt_TrkCut_pt1_sum += alltracks_pf_TrkCut_pt1[j][k].pt();
    }
    Trk_masssum[j] = p4.M();
    Trk_TrkCut_mass[j] = p4_TrkCut.M();
    Trk_TrkCut_pt05_mass[j] = p4_TrkCut05.M();
    Trk_TrkCut_pt1_mass[j] = p4_TrkCut1.M();
    Trk_TrkCut_pt[j] = pt_TrkCut_sum;
    Trk_TrkCut_pt05_pt[j] = pt_TrkCut_pt05_sum;
    Trk_TrkCut_pt1_pt[j] = pt_TrkCut_pt1_sum;
    Trk_TrkCut_ptmiss[j] = p4_TrkCut.pt();
    Trk_TrkCut_pt05_ptmiss[j] = p4_TrkCut05.pt();
    Trk_TrkCut_pt1_ptmiss[j] = p4_TrkCut1.pt();
  }
*/
  cal_p4_pt(alltracks_pf, Trk_masssum, Trk_ptsum, Trk_ptmiss, Trk_pzsum,Trk_etasum);
  cal_p4_pt(alltracks_pf_purity_pt05, Trk_Purity_pt05_mass,Trk_Purity_pt05_pt,Trk_Purity_pt05_ptmiss,Trk_Purity_pt05_pzsum,Trk_Purity_pt05_etasum);
  cal_p4_pt(alltracks_pf_purity_pt05_lowIP, Trk_Purity_pt05_lowIP_mass,Trk_Purity_pt05_lowIP_pt,Trk_Purity_pt05_lowIP_ptmiss,Trk_Purity_pt05_lowIP_pzsum,Trk_Purity_pt05_lowIP_etasum);
  cal_p4_pt(alltracks_pf_TrkCut, Trk_TrkCut_mass,Trk_TrkCut_pt,Trk_TrkCut_ptmiss,Trk_TrkCut_pzsum,Trk_TrkCut_etasum);
  cal_p4_pt(alltracks_pf_TrkCut_pt05, Trk_TrkCut_pt05_mass,Trk_TrkCut_pt05_pt,Trk_TrkCut_pt05_ptmiss,Trk_TrkCut_pt05_pzsum,Trk_TrkCut_pt05_etasum);
  cal_p4_pt(alltracks_pf_purity_pt035, Trk_Purity_pt035_mass,Trk_Purity_pt035_pt,Trk_Purity_pt035_ptmiss,Trk_Purity_pt035_pzsum,Trk_Purity_pt035_etasum);
  cal_p4_pt(alltracks_pf_purity_pt035_lowIP, Trk_Purity_pt035_lowIP_mass,Trk_Purity_pt035_lowIP_pt,Trk_Purity_pt035_lowIP_ptmiss,Trk_Purity_pt035_lowIP_pzsum,Trk_Purity_pt035_lowIP_etasum);





    // Calo and HF 

/*
    edm::Handle<edm::SortedCollection<CaloTower>> CaloTowerHandle;
    event.getByToken(CaloTowerCollection_, CaloTowerHandle);


    double sum_energy_HFp = 0;
    double sum_energy_HFm = 0;
    double sum_eT_HFp = 0;
    double sum_eT_HFm = 0;
    double leading_energy_HFp = 0;
    double leading_energy_HFm = 0;
    int nTowerHFp = 0;
    int nTowerHFm = 0;
    for (edm::SortedCollection<CaloTower>::const_iterator calo = CaloTowerHandle->begin(); calo != CaloTowerHandle->end(); ++calo) {
      double eta;
      double phi;
      double energy_HFp = 0;
      double energy_HFm = 0;
      double eT_HFp = 0;
      double eT_HFm = 0;


      phi=calo->phi();
      eta=calo->eta();

      CaloType subdetHad = GetTowerSubdetHad(eta);
      
      if(subdetHad==kHFp && (calo->hadEnergy() || calo->hadEt()) ){
        cout<<"HTp"<<endl;
        eT_HFp = calo->hadEt();
        energy_HFp = calo->hadEnergy();
        cout<<"eT_HFp "<<eT_HFp<<", energy_HFp"<<endl;
        sum_eT_HFp += eT_HFp;
        sum_energy_HFp += energy_HFp;
        nTowerHFp += 1;
        if (energy_HFp>leading_energy_HFp) leading_energy_HFp = energy_HFp;
        nBranches_->HF_hit_map->Fill(eta,phi);
        nBranches_->HFp_eta.push_back(eta);
        nBranches_->HFp_phi.push_back(phi);
        nBranches_->HFp_energy.push_back(calo->hadEnergy());
      }

      if(subdetHad==kHFm && (calo->hadEnergy() || calo->hadEt()) ){
        cout<<"HTm"<<endl;
        eT_HFm = calo->hadEt();
        energy_HFm = calo->hadEnergy();
        cout<<"eT_HFm "<<eT_HFp<<", energy_HFm"<<endl;
        sum_eT_HFm += eT_HFm;
        sum_energy_HFm += energy_HFm;
        nTowerHFm += 1;
        if (energy_HFm>leading_energy_HFm) leading_energy_HFm = energy_HFm;
        nBranches_->HF_hit_map->Fill(eta,phi);
        nBranches_->HFm_eta.push_back(eta);
        nBranches_->HFm_phi.push_back(phi);
        nBranches_->HFm_energy.push_back(calo->hadEnergy());
      }
    }
    cout<<"sum_energy_HFp"<<sum_energy_HFp<<endl;
    cout<<"sum_energy_HFm"<<sum_energy_HFm<<endl;
    cout<<"leading_energy_HFp"<<leading_energy_HFp<<endl;
    cout<<"leading_energy_HFm"<<leading_energy_HFm<<endl;
    nBranches_->calo_sum_energy_HFp = sum_energy_HFp;
    nBranches_->calo_sum_energy_HFm = sum_energy_HFm;
    nBranches_->calo_sum_eT_HFp = sum_eT_HFp;
    nBranches_->calo_sum_eT_HFm = sum_eT_HFm;
    nBranches_->calo_nTowerHFp = nTowerHFp;
    nBranches_->calo_nTowerHFm = nTowerHFm;
    nBranches_->calo_leading_energy_HFp = leading_energy_HFp;
    nBranches_->calo_leading_energy_HFm = leading_energy_HFm;
*/

//find secondary vertices by fitting the trajectory of track pairs.

  std::vector<reco::TransientTrack> trkpair(2);
  for( size_t j = 0; j < vertices_->size(); j++){
    std::vector<bool> isfromSV_sig2_pt05(alltracks_purity_pt05[j].size(),false);
    std::vector<bool> isfromSV_sig3_pt05(alltracks_purity_pt05[j].size(),false);
    std::vector<bool> isfromSV_sig5_pt05(alltracks_purity_pt05[j].size(),false);
    std::vector<bool> isfromSV_dxyp02_pt05(alltracks_purity_pt05[j].size(),false);
    std::vector<bool> isfromSV_dxyp02_pt05_PVAQLeq4(alltracks_purity_pt05[j].size(),false);
    TVector3 pvtx;
    pvtx.SetXYZ((*vertices_)[j].x(),(*vertices_)[j].y(),(*vertices_)[j].z());
   for( int trkindex1 = 0; trkindex1 < (int)alltracks_purity_pt05[j].size()-1; trkindex1++){
      trkpair[0]=alltracks_purity_pt05[j][trkindex1];
      for( size_t trkindex2 = trkindex1+1; trkindex2 < alltracks_purity_pt05[j].size(); trkindex2++){
        trkpair[1]=alltracks_purity_pt05[j][trkindex2];
        KalmanVertexFitter kalman_fitter;
        TransientVertex vtx12 = kalman_fitter.vertex(trkpair);
//        cout<<"track "<<trkindex1<<", "<<trkindex2<<endl;
        if (vtx12.isValid()&&TMath::Prob(vtx12.totalChiSquared(), int(vtx12.degreesOfFreedom()))>0.5){
          VertexDistanceXY vdistXY1;
          Measurement1D distXY1;
          distXY1 = vdistXY1.distance(vtx12.vertexState(), (*vertices_)[j]); 
//          TVector3 svtx;
//          svtx.SetXYZ(vtx12.position().x(),vtx12.position().y(),vtx12.position().z());
//          TVector3 vdiff = svtx-pvtx;
//          TVector3 p1,p2;
//          p1.SetXYZ(alltracks_pf_TrkCut[j][trkindex1].px(),alltracks_pf_TrkCut[j][trkindex1].py(),alltracks_pf_TrkCut[j][trkindex1].pz());
//          p2.SetXYZ(alltracks_pf_TrkCut[j][trkindex2].px(),alltracks_pf_TrkCut[j][trkindex2].py(),alltracks_pf_TrkCut[j][trkindex2].pz());
//          if(vdiff.Dot(p1)<0||vdiff.Dot(p2)<0) continue;
//          cout<<"distance "<<distXY1.value()<<", error "<<distXY1.error()<<", prob "<<TMath::Prob(vtx12.totalChiSquared(), int(vtx12.degreesOfFreedom()))<<endl;

          if(distXY1.value()/distXY1.error()>2.) {
            N_SV_sig2_pt05[j]++;
            isfromSV_sig2_pt05[trkindex1]=isfromSV_sig2_pt05[trkindex2]=true;
          }
          if(distXY1.value()/distXY1.error()>3.) {
            N_SV_sig3_pt05[j]++;
            isfromSV_sig3_pt05[trkindex1]=isfromSV_sig3_pt05[trkindex2]=true;
          }
          if(distXY1.value()/distXY1.error()>5.) {
            N_SV_sig5_pt05[j]++;
            isfromSV_sig5_pt05[trkindex1]=isfromSV_sig5_pt05[trkindex2]=true;
         }
          if(distXY1.value()>0.02){
            isfromSV_dxyp02_pt05[trkindex1]=true;
            isfromSV_dxyp02_pt05[trkindex2]=true;
            isfromSV_dxyp02_pt05_PVAQLeq4[trkindex1]=(alltracks_pf_purity_pt05[j][trkindex1].pvAssociationQuality()>=4);
            isfromSV_dxyp02_pt05_PVAQLeq4[trkindex2]=(alltracks_pf_purity_pt05[j][trkindex2].pvAssociationQuality()>=4);
          }
        }
      }
    }
//find tracks with vertices displaced by 0.02cm to the PV in the transverse plane
    for( size_t trkindex = 0; trkindex < alltracks_pf_purity_pt05[j].size(); trkindex++){
      if(sqrt(pow(alltracks_pf_purity_pt05[j][trkindex].vx()-(*vertices_)[j].x(),2)+pow(alltracks_pf_purity_pt05[j][trkindex].vy()-(*vertices_)[j].y(),2))>0.02){
        TVector3 ptrk;
        ptrk.SetXYZ(alltracks_pf_purity_pt05[j][trkindex].px(),alltracks_pf_purity_pt05[j][trkindex].py(),alltracks_pf_purity_pt05[j][trkindex].pz());
        TVector3 vtrk;
        vtrk.SetXYZ(alltracks_pf_purity_pt05[j][trkindex].vx(),alltracks_pf_purity_pt05[j][trkindex].vy(),alltracks_pf_purity_pt05[j][trkindex].vz());
        TVector3 vdiff = vtrk-pvtx;
        if(vdiff.Dot(ptrk)<0) continue;
        N_Trk_DisplaceVp02_pt05[j]++;
        if(isfromSV_sig2_pt05[(int)trkindex]) N_Trk_formSV_dxyp02_sig2_pt05[j]++;
      }
    }
    N_Trk_formSV_sig2_pt05[j]=count(isfromSV_sig2_pt05.begin(),isfromSV_sig2_pt05.end(),true);
    N_Trk_formSV_sig3_pt05[j]=count(isfromSV_sig3_pt05.begin(),isfromSV_sig3_pt05.end(),true);
    N_Trk_formSV_sig5_pt05[j]=count(isfromSV_sig5_pt05.begin(),isfromSV_sig5_pt05.end(),true);
    N_Trk_formSV_dxyp02_pt05[j]=count(isfromSV_dxyp02_pt05.begin(),isfromSV_dxyp02_pt05.end(),true);
    N_Trk_PVAQLeq4_formSV_dxyp02_pt05[j]=count(isfromSV_dxyp02_pt05_PVAQLeq4.begin(),isfromSV_dxyp02_pt05_PVAQLeq4.end(),true);
  }

//  for( size_t j = 0; j < vertices_->size(); j++){
//    for( size_t trkindex = 0; trkindex < alltracks_pf_TrkCut[j].size(); trkindex++){
//      if(sqrt(pow(alltracks_pf_TrkCut[j][trkindex].vx()-(*vertices_)[j].position().X(),2)+pow(alltracks_pf_TrkCut[j][trkindex].vy()-(*vertices_)[j].position().Y(),2))>0.002){
//        N_Trk_DisplaceVp002[j]++;
//        N_Trk_DisplaceVp002_SVsig2[j]++;
//      }
//    }
//  }
//Calculate TrackJets
/*
  for( size_t j = 0; j < vertices_->size(); j++){
    std::vector<PseudoJet> particles;
//    std::vector<TLorentzVector> TrackJet_P4_thisvertex;
    for( size_t trkindex = 0; trkindex < alltracks_pf[j].size(); trkindex++)
      particles.push_back(PseudoJet(alltracks_pf[j][trkindex].px(),alltracks_pf[j][trkindex].py(),alltracks_pf[j][trkindex].pz(),alltracks_pf[j][trkindex].energy()));

    ClusterSequence cs(particles, jet_def);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
    Selector select_pt       = SelectorPtMin(TrackJet_PtCut);
    vector<PseudoJet> selected_jets = select_pt(jets);
//    cout << "Clustered with " << jet_def.description() << endl;
    N_Jet[j] = selected_jets.size();
    N_Jet_total += selected_jets.size();
//    for(size_t jetindex=0; jetindex < selected_jets.size(); jetindex++){
//      TLorentzVector Jet_P4(selected_jets[jetindex].px(),selected_jets[jetindex].py(),selected_jets[jetindex].pz(),selected_jets[jetindex].E());
//      TrackJet_P4_thisvertex.push_back(Jet_P4);
//    }
//    nBranches_->TrackJet_P4.push_back(TrackJet_P4_thisvertex);
  }

  for( size_t j = 0; j < vertices_->size(); j++){
    std::vector<PseudoJet> particles;
    float TrackJet_TrkCut_E_thisvertex=0;
    float TrackJet_TrkCut_pt_thisvertex=0;
    float TrackJet_TrkCut_px_thisvertex=0;
    float TrackJet_TrkCut_py_thisvertex=0;
    float TrackJet_TrkCut_pz_thisvertex=0;
    for( size_t trkindex = 0; trkindex < alltracks_pf_TrkCut[j].size(); trkindex++)
      particles.push_back(PseudoJet(alltracks_pf_TrkCut[j][trkindex].px(),alltracks_pf_TrkCut[j][trkindex].py(),alltracks_pf_TrkCut[j][trkindex].pz(),alltracks_pf_TrkCut[j][trkindex].energy()));

    ClusterSequence cs(particles, jet_def);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
    Selector select_pt       = SelectorPtMin(TrackJet_PtCut);
    vector<PseudoJet> selected_jets = select_pt(jets);
    N_Jet_TrkCut[j] = selected_jets.size();
    N_Jet_TrkCut_total += selected_jets.size();
    for(size_t jetindex=0; jetindex < selected_jets.size(); jetindex++){
      TrackJet_TrkCut_E_thisvertex+=selected_jets[jetindex].E();
      TrackJet_TrkCut_pt_thisvertex+=selected_jets[jetindex].pt();
      TrackJet_TrkCut_px_thisvertex+=selected_jets[jetindex].px();
      TrackJet_TrkCut_py_thisvertex+=selected_jets[jetindex].py();
      TrackJet_TrkCut_pz_thisvertex+=selected_jets[jetindex].pz();
    }
    nBranches_->Instanton_TrackJet_TrkCut_E.push_back(TrackJet_TrkCut_E_thisvertex);
    nBranches_->Instanton_TrackJet_TrkCut_pt.push_back(TrackJet_TrkCut_pt_thisvertex);
    nBranches_->Instanton_TrackJet_TrkCut_px.push_back(TrackJet_TrkCut_px_thisvertex);
    nBranches_->Instanton_TrackJet_TrkCut_py.push_back(TrackJet_TrkCut_py_thisvertex);
    nBranches_->Instanton_TrackJet_TrkCut_pz.push_back(TrackJet_TrkCut_pz_thisvertex);
  }

  for( size_t j = 0; j < vertices_->size(); j++){
    std::vector<PseudoJet> particles;
    float TrackJet_TrkCut_pt05_E_thisvertex=0;
    float TrackJet_TrkCut_pt05_pt_thisvertex=0;
    float TrackJet_TrkCut_pt05_px_thisvertex=0;
    float TrackJet_TrkCut_pt05_py_thisvertex=0;
    float TrackJet_TrkCut_pt05_pz_thisvertex=0;
    for( size_t trkindex = 0; trkindex < alltracks_pf_TrkCut_pt05[j].size(); trkindex++)
      particles.push_back(PseudoJet(alltracks_pf_TrkCut_pt05[j][trkindex].px(),alltracks_pf_TrkCut_pt05[j][trkindex].py(),alltracks_pf_TrkCut_pt05[j][trkindex].pz(),alltracks_pf_TrkCut_pt05[j][trkindex].energy()));

    ClusterSequence cs(particles, jet_def);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
    Selector select_pt       = SelectorPtMin(TrackJet_PtCut);
    vector<PseudoJet> selected_jets = select_pt(jets);
    Selector select_pt10       = SelectorPtMin(10);
    vector<PseudoJet> selected_jets_tjpt10 = select_pt10(jets);
    N_Jet_TrkCut_pt05[j] = selected_jets.size();
    N_Jet_TrkCut_pt05_tjpt10[j] = selected_jets_tjpt10.size();
    N_Jet_TrkCut_pt05_total += selected_jets.size();
    for(size_t jetindex=0; jetindex < selected_jets.size(); jetindex++){
      TrackJet_TrkCut_pt05_E_thisvertex+=selected_jets[jetindex].E();
      TrackJet_TrkCut_pt05_pt_thisvertex+=selected_jets[jetindex].pt();
      TrackJet_TrkCut_pt05_px_thisvertex+=selected_jets[jetindex].px();
      TrackJet_TrkCut_pt05_py_thisvertex+=selected_jets[jetindex].py();
      TrackJet_TrkCut_pt05_pz_thisvertex+=selected_jets[jetindex].pz();
    }
    nBranches_->Instanton_TrackJet_TrkCut_pt05_E.push_back(TrackJet_TrkCut_pt05_E_thisvertex);
    nBranches_->Instanton_TrackJet_TrkCut_pt05_pt.push_back(TrackJet_TrkCut_pt05_pt_thisvertex);
    nBranches_->Instanton_TrackJet_TrkCut_pt05_px.push_back(TrackJet_TrkCut_pt05_px_thisvertex);
    nBranches_->Instanton_TrackJet_TrkCut_pt05_py.push_back(TrackJet_TrkCut_pt05_py_thisvertex);
    nBranches_->Instanton_TrackJet_TrkCut_pt05_pz.push_back(TrackJet_TrkCut_pt05_pz_thisvertex);
  }

  for( size_t j = 0; j < vertices_->size(); j++){
    std::vector<PseudoJet> particles;
    float TrackJet_TrkCut_pt1_E_thisvertex;
    float TrackJet_TrkCut_pt1_pt_thisvertex;
    float TrackJet_TrkCut_pt1_px_thisvertex;
    float TrackJet_TrkCut_pt1_py_thisvertex;
    float TrackJet_TrkCut_pt1_pz_thisvertex;
    for( size_t trkindex = 0; trkindex < alltracks_pf_TrkCut_pt1[j].size(); trkindex++)
      particles.push_back(PseudoJet(alltracks_pf_TrkCut_pt1[j][trkindex].px(),alltracks_pf_TrkCut_pt1[j][trkindex].py(),alltracks_pf_TrkCut_pt1[j][trkindex].pz(),alltracks_pf_TrkCut_pt1[j][trkindex].energy()));

    ClusterSequence cs(particles, jet_def);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
    Selector select_pt       = SelectorPtMin(TrackJet_PtCut);
    vector<PseudoJet> selected_jets = select_pt(jets);
    N_Jet_TrkCut_pt1[j] = selected_jets.size();
    N_Jet_TrkCut_pt1_total += selected_jets.size();
    for(size_t jetindex=0; jetindex < selected_jets.size(); jetindex++){
      TrackJet_TrkCut_pt1_E_thisvertex+=selected_jets[jetindex].E();
      TrackJet_TrkCut_pt1_pt_thisvertex+=selected_jets[jetindex].pt();
      TrackJet_TrkCut_pt1_px_thisvertex+=selected_jets[jetindex].px();
      TrackJet_TrkCut_pt1_py_thisvertex+=selected_jets[jetindex].py();
      TrackJet_TrkCut_pt1_pz_thisvertex+=selected_jets[jetindex].pz();
    }
    nBranches_->Instanton_TrackJet_TrkCut_pt1_E.push_back(TrackJet_TrkCut_pt1_E_thisvertex);
    nBranches_->Instanton_TrackJet_TrkCut_pt1_pt.push_back(TrackJet_TrkCut_pt1_pt_thisvertex);
    nBranches_->Instanton_TrackJet_TrkCut_pt1_px.push_back(TrackJet_TrkCut_pt1_px_thisvertex);
    nBranches_->Instanton_TrackJet_TrkCut_pt1_py.push_back(TrackJet_TrkCut_pt1_py_thisvertex);
    nBranches_->Instanton_TrackJet_TrkCut_pt1_pz.push_back(TrackJet_TrkCut_pt1_pz_thisvertex);
  }
*/
  //cal_TrackJet(alltracks_pf,nBranches_->Instanton_N_TrackJet, TrackJet_PtCut);
  //nBranches_->Instanton_N_TrackJet_total = std::accumulate(nBranches_->Instanton_N_TrackJet.begin(), nBranches_->Instanton_N_TrackJet.end(), 0);
  //cal_TrackJet(alltracks_pf_TrkCut,nBranches_->Instanton_N_TrackJet_TrkCut, TrackJet_PtCut);
  //nBranches_->Instanton_N_TrackJet_TrkCut_total = std::accumulate(nBranches_->Instanton_N_TrackJet_TrkCut.begin(), nBranches_->Instanton_N_TrackJet_TrkCut.end(), 0);
  //cal_TrackJet(alltracks_pf_TrkCut_pt05,nBranches_->Instanton_N_TrackJet_TrkCut_pt05, TrackJet_PtCut);
  //nBranches_->Instanton_N_TrackJet_TrkCut_pt05_total = std::accumulate(nBranches_->Instanton_N_TrackJet_TrkCut_pt05.begin(), nBranches_->Instanton_N_TrackJet_TrkCut_pt05.end(), 0);
  //cal_TrackJet(alltracks_pf_TrkCut_pt05,nBranches_->Instanton_N_TrackJet_TrkCut_pt05_tjpt10, 10.);
  cal_TrackJet(alltracks_pf_purity_pt05_lowIP,nBranches_->Instanton_N_TrackJet_Purity_pt05,TrackJet_px,TrackJet_py,TrackJet_pz,TrackJet_pt,TrackJet_E,TrackJet_mass,TrackJet_phi,TrackJet_eta,nBranches_->Instanton_TrackJet_deltaR,nBranches_->Instanton_TrackJet_deltaphi, TrackJet_PtCut,nBranches_->Instanton_1Jettiness_Purity_pt05,nBranches_->Instanton_2Jettiness_Purity_pt05,nBranches_->Instanton_3Jettiness_Purity_pt05,nBranches_->Instanton_4Jettiness_Purity_pt05,nBranches_->Instanton_5Jettiness_Purity_pt05,nBranches_->Instanton_6Jettiness_Purity_pt05,nBranches_->Instanton_7Jettiness_Purity_pt05,nBranches_->Instanton_8Jettiness_Purity_pt05);
  for(size_t j = 0; j < vertices_->size(); j++){
    nBranches_->Instanton_TrackJet_px.push_back(TrackJet_px[j]);
    nBranches_->Instanton_TrackJet_py.push_back(TrackJet_py[j]);
    nBranches_->Instanton_TrackJet_pz.push_back(TrackJet_pz[j]);
    nBranches_->Instanton_TrackJet_pt.push_back(TrackJet_pt[j]);
    nBranches_->Instanton_TrackJet_E.push_back(TrackJet_E[j]);
    nBranches_->Instanton_TrackJet_mass.push_back(TrackJet_mass[j]);
    nBranches_->Instanton_TrackJet_phi.push_back(TrackJet_phi[j]);
    nBranches_->Instanton_TrackJet_eta.push_back(TrackJet_eta[j]);
  }
  //cal_TrackJet(alltracks_pf_purity_pt05,nBranches_->Instanton_N_TrackJet_Purity_pt05_tjpt10, 10.);
  //cal_TrackJet(alltracks_pf_TrkCut_pt1,nBranches_->Instanton_N_TrackJet_TrkCut_pt1, TrackJet_PtCut);
  //nBranches_->Instanton_N_TrackJet_TrkCut_pt1_total = std::accumulate(nBranches_->Instanton_N_TrackJet_TrkCut_pt1.begin(), nBranches_->Instanton_N_TrackJet_TrkCut_pt1.end(), 0);
//Calculate event shape (bservables: Spherocity, Thrust, Broadening



  Float_t Spherocity, Thrust, Broaden, TransverseSpherocity,TransverseThrust,TransverseSpherocity_pTunweighted;
  TVector3 taxis;
  TVector3 transtaxis;
  Float_t Rapiditygap,Xix,Xiy;
  for(size_t j = 0; j < vertices_->size(); j++){
    calcS_T_B(alltracks_pf_purity_pt05_lowIP[j], Spherocity,Thrust,Broaden,TransverseSpherocity,TransverseThrust,TransverseSpherocity_pTunweighted, taxis,transtaxis, min(nSeed_ThrustCalculation,(int)alltracks_pf_purity_pt05_lowIP[j].size()));
    nBranches_->Instanton_Trk_spherocity.push_back(Spherocity);
    nBranches_->Instanton_Trk_thrust.push_back(Thrust);
    nBranches_->Instanton_Trk_broaden.push_back(Broaden);
    nBranches_->Instanton_Trk_transversespherocity.push_back(TransverseSpherocity);
    nBranches_->Instanton_Trk_transversethrust.push_back(TransverseThrust);
    nBranches_->Instanton_Trk_transversespherocity_pTunweighted.push_back(TransverseSpherocity_pTunweighted);
    calc_rapidity(alltracks_pf_purity_pt05_lowIP[j], Rapiditygap,Xix,Xiy);
    nBranches_->Instanton_Trk_rapiditygap.push_back(Rapiditygap);
    nBranches_->Instanton_Trk_xix.push_back(Xix);
    nBranches_->Instanton_Trk_xiy.push_back(Xiy);
  }

  Float_t Spherocity_pt035, Thrust_pt035, Broaden_pt035, TransverseSpherocity_pt035,TransverseThrust_pt035,TransverseSpherocity_pTunweighted_pt035;
  TVector3 taxis_pt035;
  TVector3 transtaxis_pt035;
  Float_t Rapiditygap_pt035,Xix_pt035,Xiy_pt035;
  for(size_t j = 0; j < vertices_->size(); j++){ 
    calcS_T_B(alltracks_pf_purity_pt035_lowIP[j], Spherocity_pt035,Thrust_pt035,Broaden_pt035,TransverseSpherocity_pt035,TransverseThrust_pt035,TransverseSpherocity_pTunweighted_pt035, taxis_pt035,transtaxis_pt035, min(nSeed_ThrustCalculation,(int)alltracks_pf_purity_pt035_lowIP[j].size()));
    nBranches_->Instanton_Trk_spherocity_pt035.push_back(Spherocity_pt035);
    nBranches_->Instanton_Trk_thrust_pt035.push_back(Thrust_pt035);
    nBranches_->Instanton_Trk_broaden_pt035.push_back(Broaden_pt035);
    nBranches_->Instanton_Trk_transversespherocity_pt035.push_back(TransverseSpherocity_pt035);
    nBranches_->Instanton_Trk_transversethrust_pt035.push_back(TransverseThrust_pt035);
    nBranches_->Instanton_Trk_transversespherocity_pTunweighted_pt035.push_back(TransverseSpherocity_pTunweighted_pt035);
    calc_rapidity(alltracks_pf_purity_pt035_lowIP[j], Rapiditygap_pt035,Xix_pt035,Xiy_pt035);
    nBranches_->Instanton_Trk_rapiditygap_pt035.push_back(Rapiditygap_pt035);
    nBranches_->Instanton_Trk_xix_pt035.push_back(Xix_pt035);
    nBranches_->Instanton_Trk_xiy_pt035.push_back(Xiy_pt035);
  }


  Float_t Spherocity_eta1, Thrust_eta1, Broaden_eta1, TransverseSpherocity_eta1,TransverseThrust_eta1, TransverseSpherocity_pTunweighted_eta1;
  TVector3 taxis_eta1;
  TVector3 transtaxis_eta1;
  for(size_t j = 0; j < vertices_->size(); j++){
    calcS_T_B(alltracks_pf_purity_eta1pt05[j], Spherocity_eta1,Thrust_eta1,Broaden_eta1,TransverseSpherocity_eta1,TransverseThrust_eta1,TransverseSpherocity_pTunweighted_eta1, taxis_eta1,transtaxis_eta1, min(nSeed_ThrustCalculation,(int)alltracks_pf_purity_eta1pt05[j].size()));
    nBranches_->Instanton_Trk_spherocity_eta1.push_back(Spherocity_eta1);
    //nBranches_->Instanton_Trk_thrust_eta1.push_back(Thrust_eta1);
    //nBranches_->Instanton_Trk_broaden_eta1.push_back(Broaden_eta1);
    nBranches_->Instanton_Trk_transversespherocity_eta1.push_back(TransverseSpherocity_eta1);
    nBranches_->Instanton_Trk_transversethrust_eta1.push_back(TransverseThrust_eta1);
  }

  Float_t Spherocity_eta2, Thrust_eta2, Broaden_eta2, TransverseSpherocity_eta2,TransverseThrust_eta2,TransverseSpherocity_pTunweighted_eta2;
  TVector3 taxis_eta2;
  TVector3 transtaxis_eta2;
  for(size_t j = 0; j < vertices_->size(); j++){
    calcS_T_B(alltracks_pf_purity_eta2pt05[j], Spherocity_eta2,Thrust_eta2,Broaden_eta2,TransverseSpherocity_eta2,TransverseThrust_eta2, TransverseSpherocity_pTunweighted_eta2, taxis_eta2,transtaxis_eta2, min(nSeed_ThrustCalculation,(int)alltracks_pf_purity_eta2pt05[j].size()));
    nBranches_->Instanton_Trk_spherocity_eta2.push_back(Spherocity_eta2);
    //nBranches_->Instanton_Trk_thrust_eta2.push_back(Thrust_eta2);
    //nBranches_->Instanton_Trk_broaden_eta2.push_back(Broaden_eta2);
    nBranches_->Instanton_Trk_transversespherocity_eta2.push_back(TransverseSpherocity_eta2);
    nBranches_->Instanton_Trk_transversethrust_eta2.push_back(TransverseThrust_eta2);
  }

  Float_t Spherocity_eta05, Thrust_eta05, Broaden_eta05, TransverseSpherocity_eta05,TransverseThrust_eta05,TransverseSpherocity_pTunweighted_eta05;
  TVector3 taxis_eta05;
  TVector3 transtaxis_eta05;
  for(size_t j = 0; j < vertices_->size(); j++){
    calcS_T_B(alltracks_pf_purity_eta05pt05[j], Spherocity_eta05,Thrust_eta05,Broaden_eta05,TransverseSpherocity_eta05,TransverseThrust_eta05,TransverseSpherocity_pTunweighted_eta05, taxis_eta05,transtaxis_eta05, min(nSeed_ThrustCalculation,(int)alltracks_pf_purity_eta05pt05[j].size()));
    nBranches_->Instanton_Trk_spherocity_eta05.push_back(Spherocity_eta05);
    nBranches_->Instanton_Trk_transversespherocity_eta05.push_back(TransverseSpherocity_eta05);
    nBranches_->Instanton_Trk_transversethrust_eta05.push_back(TransverseThrust_eta05);
  }

  Float_t Spherocity_eta05to1, Thrust_eta05to1, Broaden_eta05to1, TransverseSpherocity_eta05to1,TransverseThrust_eta05to1,TransverseSpherocity_pTunweighted_eta05to1;
  TVector3 taxis_eta05to1;
  TVector3 transtaxis_eta05to1;
  for(size_t j = 0; j < vertices_->size(); j++){
    calcS_T_B(alltracks_pf_purity_eta05to1pt05[j], Spherocity_eta05to1,Thrust_eta05to1,Broaden_eta05to1,TransverseSpherocity_eta05to1,TransverseThrust_eta05to1,TransverseSpherocity_pTunweighted_eta05to1, taxis_eta05to1,transtaxis_eta05to1, min(nSeed_ThrustCalculation,(int)alltracks_pf_purity_eta05to1pt05[j].size()));
    nBranches_->Instanton_Trk_transversespherocity_eta05to1.push_back(TransverseSpherocity_eta05to1);
    nBranches_->Instanton_Trk_transversethrust_eta05to1.push_back(TransverseThrust_eta05to1);
  }

  Float_t Spherocity_eta1to1p5, Thrust_eta1to1p5, Broaden_eta1to1p5, TransverseSpherocity_eta1to1p5,TransverseThrust_eta1to1p5,TransverseSpherocity_pTunweighted_eta1to1p5;
  TVector3 taxis_eta1to1p5;
  TVector3 transtaxis_eta1to1p5;
  for(size_t j = 0; j < vertices_->size(); j++){
    calcS_T_B(alltracks_pf_purity_eta1to1p5pt05[j], Spherocity_eta1to1p5,Thrust_eta1to1p5,Broaden_eta1to1p5,TransverseSpherocity_eta1to1p5,TransverseThrust_eta1to1p5,TransverseSpherocity_pTunweighted_eta1to1p5, taxis_eta1to1p5,transtaxis_eta1to1p5, min(nSeed_ThrustCalculation,(int)alltracks_pf_purity_eta1to1p5pt05[j].size()));
    nBranches_->Instanton_Trk_transversespherocity_eta1to1p5.push_back(TransverseSpherocity_eta1to1p5);
    nBranches_->Instanton_Trk_transversethrust_eta1to1p5.push_back(TransverseThrust_eta1to1p5);
  }


  Float_t Spherocity_eta1p5to2, Thrust_eta1p5to2, Broaden_eta1p5to2, TransverseSpherocity_eta1p5to2,TransverseThrust_eta1p5to2,TransverseSpherocity_pTunweighted_eta1p5to2;
  TVector3 taxis_eta1p5to2;
  TVector3 transtaxis_eta1p5to2;
  for(size_t j = 0; j < vertices_->size(); j++){
    calcS_T_B(alltracks_pf_purity_eta1p5to2pt05[j], Spherocity_eta1p5to2,Thrust_eta1p5to2,Broaden_eta1p5to2,TransverseSpherocity_eta1p5to2,TransverseThrust_eta1p5to2,TransverseSpherocity_pTunweighted_eta1p5to2, taxis_eta1p5to2,transtaxis_eta1p5to2, min(nSeed_ThrustCalculation,(int)alltracks_pf_purity_eta1p5to2pt05[j].size()));
    nBranches_->Instanton_Trk_transversespherocity_eta1p5to2.push_back(TransverseSpherocity_eta1p5to2);
    nBranches_->Instanton_Trk_transversethrust_eta1p5to2.push_back(TransverseThrust_eta1p5to2);
  } 

  Float_t Spherocity_eta2to2p4, Thrust_eta2to2p4, Broaden_eta2to2p4, TransverseSpherocity_eta2to2p4,TransverseThrust_eta2to2p4,TransverseSpherocity_pTunweighted_eta2to2p4;
  TVector3 taxis_eta2to2p4;
  TVector3 transtaxis_eta2to2p4;
  for(size_t j = 0; j < vertices_->size(); j++){
    calcS_T_B(alltracks_pf_purity_eta2to2p4pt05[j], Spherocity_eta2to2p4,Thrust_eta2to2p4,Broaden_eta2to2p4,TransverseSpherocity_eta2to2p4,TransverseThrust_eta2to2p4,TransverseSpherocity_pTunweighted_eta2to2p4, taxis_eta2to2p4,transtaxis_eta2to2p4, min(nSeed_ThrustCalculation,(int)alltracks_pf_purity_eta2to2p4pt05[j].size()));
    nBranches_->Instanton_Trk_transversespherocity_eta2to2p4.push_back(TransverseSpherocity_eta2to2p4);
    nBranches_->Instanton_Trk_transversethrust_eta2to2p4.push_back(TransverseThrust_eta2to2p4);
  }


  Float_t Spherocity_TrkCut, Thrust_TrkCut, Broaden_TrkCut, TransverseSpherocity_TrkCut,TransverseThrust_TrkCut,TransverseSpherocity_pTunweighted_TrkCut;
  TVector3 taxis_TrkCut;
  TVector3 transtaxis_TrkCut;

  for(size_t j = 0; j < vertices_->size(); j++){
    calcS_T_B(alltracks_pf_TrkCut[j], Spherocity_TrkCut,Thrust_TrkCut,Broaden_TrkCut,TransverseSpherocity_TrkCut,TransverseThrust_TrkCut, TransverseSpherocity_pTunweighted_TrkCut, taxis_TrkCut,transtaxis_TrkCut, min(nSeed_ThrustCalculation,(int)alltracks_pf_TrkCut[j].size()));
    nBranches_->Instanton_Trk_TrkCut_spherocity.push_back(Spherocity_TrkCut);
    nBranches_->Instanton_Trk_TrkCut_thrust.push_back(Thrust_TrkCut);
    nBranches_->Instanton_Trk_TrkCut_broaden.push_back(Broaden_TrkCut);
    nBranches_->Instanton_Trk_TrkCut_transversespherocity.push_back(TransverseSpherocity_TrkCut);
    nBranches_->Instanton_Trk_TrkCut_transversethrust.push_back(TransverseThrust_TrkCut);
  }

  Float_t Spherocity_TrkCut_pt05, Thrust_TrkCut_pt05, Broaden_TrkCut_pt05, TransverseSpherocity_TrkCut_pt05,TransverseThrust_TrkCut_pt05,TransverseSpherocity_pTunweighted_TrkCut_pt05;
  TVector3 taxis_TrkCut_pt05;
  TVector3 transtaxis_TrkCut_pt05;
  for(size_t j = 0; j < vertices_->size(); j++){
    calcS_T_B(alltracks_pf_TrkCut_pt05[j], Spherocity_TrkCut_pt05,Thrust_TrkCut_pt05,Broaden_TrkCut_pt05,TransverseSpherocity_TrkCut_pt05,TransverseThrust_TrkCut_pt05,TransverseSpherocity_pTunweighted_TrkCut_pt05, taxis_TrkCut_pt05,transtaxis_TrkCut_pt05, min(nSeed_ThrustCalculation,(int)alltracks_pf_TrkCut_pt05[j].size()));
    nBranches_->Instanton_Trk_TrkCut_pt05_spherocity.push_back(Spherocity_TrkCut_pt05);
    nBranches_->Instanton_Trk_TrkCut_pt05_thrust.push_back(Thrust_TrkCut_pt05);
    nBranches_->Instanton_Trk_TrkCut_pt05_broaden.push_back(Broaden_TrkCut_pt05);
    nBranches_->Instanton_Trk_TrkCut_pt05_transversespherocity.push_back(TransverseSpherocity_TrkCut_pt05);
    nBranches_->Instanton_Trk_TrkCut_pt05_transversethrust.push_back(TransverseThrust_TrkCut_pt05);



/*    cout<<"reco tracks="<<(int)alltracks_pf_TrkCut_pt05[j].size()<<" sph="<<Spherocity_TrkCut_pt05<<" thrust="<<Thrust_TrkCut_pt05<<" broaden="<<Broaden_TrkCut_pt05<<endl;
    if(j==0){
      for(size_t trkindex=0; trkindex<alltracks_pf_TrkCut_pt05[j].size(); trkindex++)
        cout<<"trk px="<<alltracks_pf_TrkCut_pt05[j][trkindex].px()<<" py="<<alltracks_pf_TrkCut_pt05[j][trkindex].py()<<" pz="<<alltracks_pf_TrkCut_pt05[j][trkindex].pz()<<" pt="<<alltracks_pf_TrkCut_pt05[j][trkindex].pt()<<endl;
    }
*/
  }


  if(runOnMC_){
    for(size_t kk=0; kk< ChargedFSParticle_eta2p4.size();kk++){
      bool matched=false;
      int ivertex_match=0;
      int ipf_match=0;
      float dr_min = 100.;
      for(size_t ivertex=0; ivertex< alltracks_pf.size(); ivertex++){
        for (size_t ipf=0; ipf < alltracks_pf[ivertex].size(); ipf++){
          if( ChargedFSParticle_eta2p4[kk].charge()==alltracks_pf[ivertex][ipf].charge() ){
            Float_t dr = reco::deltaR(ChargedFSParticle_eta2p4[kk].eta(),ChargedFSParticle_eta2p4[kk].phi(),
                                      alltracks_pf[ivertex][ipf].eta(),alltracks_pf[ivertex][ipf].phi());
            Float_t ptratio = alltracks_pf[ivertex][ipf].bestTrack()->pt()/ChargedFSParticle_eta2p4[kk].pt();
            if(ptratio > 0. && ptratio < 2. && dr<0.3){
            //if(dr<0.3){
              if(!matched){
                ivertex_match = ivertex;
                ipf_match = ipf;
                dr_min=dr;
                matched = true;
              }
              else{
                if(dr<dr_min){
                  ivertex_match = ivertex;
                  ipf_match = ipf;
                  dr_min=dr;
                }
              }
            }
          }
        }
      }
      if(matched){
        Trk_matched[ivertex_match][ipf_match] = 1;
        nBranches_->ChargedFSParticle_eta2p4_matched[kk] = 1;
        if(nBranches_->ChargedFSParticle_eta2p4_pt[kk]>0.5) nBranches_->Instanton_N_gen_ChargedFSParticle_eta2p4pt05_matched++;
        if(Trk_pt[ivertex_match][ipf_match]>0.5&&Trk_eta[ivertex_match][ipf_match]<2.4&&Trk_eta[ivertex_match][ipf_match]>-2.4&&Trk_isHighPurity[ivertex_match][ipf_match])   N_Trk_highPurity_pt05_matched[ivertex_match]++;
      }
    }  
    
  }

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
// Filter out the events not firing the ZeroBias trigger if running on the ZeroBias data
//  std::cout<<"on MC?"<<runOnMC_<<std::endl;
//  std::cout<<"on ZeroBias?"<<runOnZeroBias_<<std::endl;
//  std::cout<<"is HLT_ZeroBias_v6 triggered?"<<(std::find(FiredTriggers.begin(),FiredTriggers.end(),"HLT_ZeroBias_v6")!=FiredTriggers.end())<<std::endl;
  if(!runOnMC_&&runOnZeroBias_&&(std::find(FiredTriggers.begin(),FiredTriggers.end(),"HLT_ZeroBias_v6")==FiredTriggers.end())) return false;
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
  if(runOnMC_&&(runOnHerwigInstanton_||runOnSherpaInstanton_)){
  //  std::cout<<"calculate distance to PV"<<std::endl;
    for( size_t nn = 0; nn < vertices_->size(); nn++){
      vtx_genvertex_PVDistance[nn] = (float)std::sqrt(std::pow((*vertices_)[nn].position().x()-PV_X,2)+std::pow((*vertices_)[nn].position().y()-PV_Y,2)+std::pow((*vertices_)[nn].position().z()-PV_Z,2));
//      std::cout<<"distance "<<vtx_genvertex_PVDistance[nn]<<std::endl;
      nBranches_->Instanton_vtx_genvertex_PVDistance.push_back(vtx_genvertex_PVDistance[nn]);
    }
  }
  nBranches_->Instanton_N_Trk_total = N_Trk_total;
  nBranches_->Instanton_N_Trk_Displaced_total = N_Trk_Displaced_total;
  nBranches_->Instanton_N_Trk_goodDisplaced_total = N_Trk_goodDisplaced_total;
//  nBranches_->Instanton_N_TrackJet_total = N_Jet_total;
//  nBranches_->Instanton_N_TrackJet_TrkCut_total = N_Jet_TrkCut_total;
 // nBranches_->Instanton_N_TrackJet_TrkCut_pt05_total = N_Jet_TrkCut_pt05_total;
//  nBranches_->Instanton_N_TrackJet_TrkCut_pt1_total = N_Jet_TrkCut_pt1_total; 
  nBranches_->PV_N_good_hasgoodMuon = PV_N_good_hasgoodMuon;
  nBranches_->PV_N_good_hasgoodPFMuon = PV_N_good_hasgoodPFMuon;

//cout<<"Event processed. Begin to write the file "<<endl;
  for( size_t jj = 0; jj < vertices_->size(); ++jj ){
    nBranches_->Instanton_N_Trk.push_back(N_Trk[jj]);
    nBranches_->Instanton_N_Trk_highPurity.push_back(N_Trk_highPurity[jj]);
    nBranches_->Instanton_N_Trk_highPurity_pt05.push_back(N_Trk_highPurity_pt05[jj]);
    nBranches_->Instanton_N_Trk_highPurity_pt05_lowIP.push_back(N_Trk_highPurity_pt05_lowIP[jj]);
    nBranches_->Instanton_N_Trk_highPurity_pt05_matched.push_back(N_Trk_highPurity_pt05_matched[jj]);
    nBranches_->Instanton_N_Trk_highPurity_pt035.push_back(N_Trk_highPurity_pt035[jj]);
    nBranches_->Instanton_N_Trk_highPurity_pt035_lowIP.push_back(N_Trk_highPurity_pt035_lowIP[jj]);
    nBranches_->Instanton_N_Trk_highPurity_pt08.push_back(N_Trk_highPurity_pt08[jj]);
    nBranches_->Instanton_N_Trk_highPurity_eta1pt05.push_back(N_Trk_highPurity_eta1pt05[jj]);
    nBranches_->Instanton_N_Trk_highPurity_eta2pt05.push_back(N_Trk_highPurity_eta2pt05[jj]);
    nBranches_->Instanton_N_Trk_highPurity_eta05pt05.push_back(N_Trk_highPurity_eta05pt05[jj]);
    nBranches_->Instanton_N_Trk_highPurity_eta1p5pt05.push_back(N_Trk_highPurity_eta1p5pt05[jj]);
    nBranches_->Instanton_N_Trk_highPurity_eta05to1pt05.push_back(N_Trk_highPurity_eta05to1pt05[jj]);
    nBranches_->Instanton_N_Trk_highPurity_eta1to1p5pt05.push_back(N_Trk_highPurity_eta1to1p5pt05[jj]);
    nBranches_->Instanton_N_Trk_highPurity_eta1p5to2pt05.push_back(N_Trk_highPurity_eta1p5to2pt05[jj]);
    nBranches_->Instanton_N_Trk_highPurity_eta2to2p4pt05.push_back(N_Trk_highPurity_eta2to2p4pt05[jj]);
//    nBranches_->Instanton_N_Trk_PVAssociationQuality0.push_back(N_Trk_PVAssociationQuality0[jj]);
//    nBranches_->Instanton_N_Trk_PVAssociationQuality1.push_back(N_Trk_PVAssociationQuality1[jj]);
//    nBranches_->Instanton_N_Trk_PVAssociationQuality4.push_back(N_Trk_PVAssociationQuality4[jj]);
//    nBranches_->Instanton_N_Trk_PVAssociationQuality5.push_back(N_Trk_PVAssociationQuality5[jj]);
//    nBranches_->Instanton_N_Trk_PVAssociationQuality6.push_back(N_Trk_PVAssociationQuality6[jj]);
//    nBranches_->Instanton_N_Trk_PVAssociationQuality7.push_back(N_Trk_PVAssociationQuality7[jj]);
    nBranches_->Instanton_N_Trk_PVAssociationQualityLeq4_highPurity.push_back(N_Trk_PVAssociationQualityLeq4_highPurity[jj]);
    nBranches_->Instanton_N_Trk_PVAssociationQualityLeq4_highPurity_pt05.push_back(N_Trk_PVAssociationQualityLeq4_highPurity_pt05[jj]);
    nBranches_->Instanton_N_Trk_PVAssociationQualityLeq4_highPurity_pt08.push_back(N_Trk_PVAssociationQualityLeq4_highPurity_pt08[jj]);
    nBranches_->Instanton_N_Trk_Displaced.push_back(N_Trk_Displaced[jj]);
    nBranches_->Instanton_N_Trk_Displaced_highPurity.push_back(N_Trk_Displaced_highPurity[jj]);
    nBranches_->Instanton_N_Trk_Displaced_highPurity_pt05.push_back(N_Trk_Displaced_highPurity_pt05[jj]);
    nBranches_->Instanton_N_Trk_Displaced_highPurity_pt08.push_back(N_Trk_Displaced_highPurity_pt08[jj]);
//    nBranches_->Instanton_N_Trk_Displaced_PVAssociationQuality0.push_back(N_Trk_Displaced_PVAssociationQuality0[jj]);
//    nBranches_->Instanton_N_Trk_Displaced_PVAssociationQuality1.push_back(N_Trk_Displaced_PVAssociationQuality1[jj]);
//    nBranches_->Instanton_N_Trk_Displaced_PVAssociationQuality4.push_back(N_Trk_Displaced_PVAssociationQuality4[jj]);
//    nBranches_->Instanton_N_Trk_Displaced_PVAssociationQuality5.push_back(N_Trk_Displaced_PVAssociationQuality5[jj]);
//    nBranches_->Instanton_N_Trk_Displaced_PVAssociationQuality6.push_back(N_Trk_Displaced_PVAssociationQuality6[jj]);
//    nBranches_->Instanton_N_Trk_Displaced_PVAssociationQuality7.push_back(N_Trk_Displaced_PVAssociationQuality7[jj]);
    nBranches_->Instanton_N_Trk_Displaced_PVAssociationQualityLeq4_highPurity.push_back(N_Trk_Displaced_PVAssociationQualityLeq4_highPurity[jj]);
    nBranches_->Instanton_N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt05.push_back(N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt05[jj]);
    nBranches_->Instanton_N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt08.push_back(N_Trk_Displaced_PVAssociationQualityLeq4_highPurity_pt08[jj]);
//    nBranches_->Instanton_N_Trk_sig2Displaced_highPurity.push_back(N_Trk_sig2Displaced_highPurity[jj]);
//    nBranches_->Instanton_N_Trk_sig3Displaced_highPurity.push_back(N_Trk_sig3Displaced_highPurity[jj]);
//    nBranches_->Instanton_N_Trk_sig5Displaced_highPurity.push_back(N_Trk_sig5Displaced_highPurity[jj]);
//    nBranches_->Instanton_N_Trk_sig2Displaced_PVAssociationQualityLeq4_highPurity.push_back(N_Trk_sig2Displaced_PVAssociationQualityLeq4_highPurity[jj]);
//    nBranches_->Instanton_N_Trk_sig3Displaced_PVAssociationQualityLeq4_highPurity.push_back(N_Trk_sig3Displaced_PVAssociationQualityLeq4_highPurity[jj]);
//    nBranches_->Instanton_N_Trk_sig5Displaced_PVAssociationQualityLeq4_highPurity.push_back(N_Trk_sig5Displaced_PVAssociationQualityLeq4_highPurity[jj]);
//    nBranches_->Instanton_N_Trk_fromPV2_highPurity.push_back(N_Trk_fromPV2_highPurity[jj]);
//    nBranches_->Instanton_N_Trk_Displaced_fromPV2_highPurity.push_back(N_Trk_Displaced_fromPV2_highPurity[jj]);
//    nBranches_->Instanton_N_Trk_Displaced_fromPV2_highPurity_pt05.push_back(N_Trk_Displaced_fromPV2_highPurity_pt05[jj]);
//    nBranches_->Instanton_N_Trk_Displaced_fromPV2_highPurity_pt08.push_back(N_Trk_Displaced_fromPV2_highPurity_pt08[jj]);
//    nBranches_->Instanton_N_Trk_Displaced_fromPV2_highPurity_pt1.push_back(N_Trk_Displaced_fromPV2_highPurity_pt1[jj]);
//    nBranches_->Instanton_N_Trk_goodDisplaced_fromPV2_highPurity_sig3.push_back(N_Trk_goodDisplaced_fromPV2_highPurity_sig3[jj]);
//    nBranches_->Instanton_N_Trk_goodDisplaced_fromPV2_highPurity_sig3_pt05.push_back(N_Trk_goodDisplaced_fromPV2_highPurity_sig3_pt05[jj]);
//    nBranches_->Instanton_N_Trk_goodDisplaced_fromPV2_highPurity_sig3_pt08.push_back(N_Trk_goodDisplaced_fromPV2_highPurity_sig3_pt08[jj]);
//    nBranches_->Instanton_N_Trk_goodDisplaced_fromPV2_highPurity_sig3_pt1.push_back(N_Trk_goodDisplaced_fromPV2_highPurity_sig3_pt1[jj]);
    nBranches_->Instanton_N_SV_sig2_pt05.push_back(N_SV_sig2_pt05[jj]);
    nBranches_->Instanton_N_SV_sig3_pt05.push_back(N_SV_sig3_pt05[jj]);
    nBranches_->Instanton_N_SV_sig5_pt05.push_back(N_SV_sig5_pt05[jj]);
    nBranches_->Instanton_N_Trk_formSV_sig2_pt05.push_back(N_Trk_formSV_sig2_pt05[jj]);
    nBranches_->Instanton_N_Trk_formSV_sig3_pt05.push_back(N_Trk_formSV_sig3_pt05[jj]);
    nBranches_->Instanton_N_Trk_formSV_sig5_pt05.push_back(N_Trk_formSV_sig5_pt05[jj]);
    nBranches_->Instanton_N_Trk_formSV_dxyp02_pt05.push_back(N_Trk_formSV_dxyp02_pt05[jj]);
    nBranches_->Instanton_N_Trk_PVAQLeq4_formSV_dxyp02_pt05.push_back(N_Trk_PVAQLeq4_formSV_dxyp02_pt05[jj]);
    nBranches_->Instanton_N_Trk_DisplaceVp02_pt05.push_back(N_Trk_DisplaceVp02_pt05[jj]);
    nBranches_->Instanton_N_Trk_formSV_dxyp02_sig2_pt05.push_back(N_Trk_formSV_dxyp02_sig2_pt05[jj]);
    nBranches_->Instanton_N_Trk_goodDisplaced.push_back(N_Trk_goodDisplaced[jj]);
    nBranches_->Instanton_N_Trk_goodDisplaced_highPurity.push_back(N_Trk_goodDisplaced_highPurity[jj]);
    nBranches_->Instanton_N_Trk_goodDisplaced_highPurity_pt05.push_back(N_Trk_goodDisplaced_highPurity_pt05[jj]);
    nBranches_->Instanton_N_Trk_goodDisplaced_highPurity_pt08.push_back(N_Trk_goodDisplaced_highPurity_pt08[jj]);
//    nBranches_->Instanton_N_Trk_goodDisplaced_PVAssociationQuality0.push_back(N_Trk_goodDisplaced_PVAssociationQuality0[jj]);
//    nBranches_->Instanton_N_Trk_goodDisplaced_PVAssociationQuality1.push_back(N_Trk_goodDisplaced_PVAssociationQuality1[jj]);
//    nBranches_->Instanton_N_Trk_goodDisplaced_PVAssociationQuality4.push_back(N_Trk_goodDisplaced_PVAssociationQuality4[jj]);
//    nBranches_->Instanton_N_Trk_goodDisplaced_PVAssociationQuality5.push_back(N_Trk_goodDisplaced_PVAssociationQuality5[jj]);
//    nBranches_->Instanton_N_Trk_goodDisplaced_PVAssociationQuality6.push_back(N_Trk_goodDisplaced_PVAssociationQuality6[jj]);
//    nBranches_->Instanton_N_Trk_goodDisplaced_PVAssociationQuality7.push_back(N_Trk_goodDisplaced_PVAssociationQuality7[jj]);
    nBranches_->Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity.push_back(N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity[jj]);
    nBranches_->Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt05.push_back(N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt05[jj]);
    nBranches_->Instanton_N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt08.push_back(N_Trk_goodDisplaced_PVAssociationQualityLeq4_highPurity_pt08[jj]);
    nBranches_->Instanton_Trk_mass.push_back(Trk_masssum[jj]);
    nBranches_->Instanton_Trk_TrkCut_mass.push_back(Trk_TrkCut_mass[jj]);
    nBranches_->Instanton_Trk_TrkCut_pt05_mass.push_back(Trk_TrkCut_pt05_mass[jj]);
    nBranches_->Instanton_Trk_Purity_pt05_lowIP_mass.push_back(Trk_Purity_pt05_lowIP_mass[jj]);
    nBranches_->Instanton_Trk_Purity_pt035_mass.push_back(Trk_Purity_pt035_mass[jj]);
    nBranches_->Instanton_Trk_Purity_pt035_lowIP_mass.push_back(Trk_Purity_pt035_lowIP_mass[jj]);
    nBranches_->Instanton_Trk_TrkCut_pt.push_back(Trk_TrkCut_pt[jj]);
    nBranches_->Instanton_Trk_TrkCut_pt05_pt.push_back(Trk_TrkCut_pt05_pt[jj]);
    nBranches_->Instanton_Trk_Purity_pt05_pt.push_back(Trk_Purity_pt05_pt[jj]);
    nBranches_->Instanton_Trk_Purity_pt05_lowIP_pt.push_back(Trk_Purity_pt05_lowIP_pt[jj]);
    nBranches_->Instanton_Trk_Purity_pt035_pt.push_back(Trk_Purity_pt035_pt[jj]);
    nBranches_->Instanton_Trk_Purity_pt035_lowIP_pt.push_back(Trk_Purity_pt035_lowIP_pt[jj]);
    nBranches_->Instanton_Trk_TrkCut_ptmiss.push_back(Trk_TrkCut_ptmiss[jj]);
    nBranches_->Instanton_Trk_TrkCut_pt05_ptmiss.push_back(Trk_TrkCut_pt05_ptmiss[jj]);
    nBranches_->Instanton_Trk_Purity_pt05_lowIP_ptmiss.push_back(Trk_Purity_pt05_lowIP_ptmiss[jj]);
    nBranches_->Instanton_Trk_Purity_pt035_ptmiss.push_back(Trk_Purity_pt035_ptmiss[jj]);
    nBranches_->Instanton_Trk_Purity_pt035_lowIP_ptmiss.push_back(Trk_Purity_pt035_lowIP_ptmiss[jj]);
    nBranches_->Instanton_Trk_Purity_pt05_eta.push_back(Trk_Purity_pt05_etasum[jj]);
    nBranches_->Instanton_Trk_Purity_pt05_lowIP_eta.push_back(Trk_Purity_pt05_lowIP_etasum[jj]);
    nBranches_->Instanton_Trk_Purity_pt035_eta.push_back(Trk_Purity_pt035_etasum[jj]);
    nBranches_->Instanton_Trk_Purity_pt035_lowIP_eta.push_back(Trk_Purity_pt035_lowIP_etasum[jj]);
    nBranches_->Instanton_Trk_Purity_pt05_pz.push_back(Trk_Purity_pt05_pzsum[jj]);
    nBranches_->Instanton_Trk_Purity_pt05_lowIP_pz.push_back(Trk_Purity_pt05_lowIP_pzsum[jj]);
    nBranches_->Instanton_Trk_Purity_pt035_pz.push_back(Trk_Purity_pt035_pzsum[jj]);
    nBranches_->Instanton_Trk_Purity_pt035_lowIP_pz.push_back(Trk_Purity_pt035_lowIP_pzsum[jj]);
//    nBranches_->Instanton_N_TrackJet.push_back(N_Jet[jj]);
//    nBranches_->Instanton_N_TrackJet_TrkCut.push_back(N_Jet_TrkCut[jj]);
//    nBranches_->Instanton_N_TrackJet_TrkCut_pt05.push_back(N_Jet_TrkCut_pt05[jj]);
//    nBranches_->Instanton_N_TrackJet_TrkCut_pt05_tjpt10.push_back(N_Jet_TrkCut_pt05_tjpt10[jj]);
//    nBranches_->Instanton_N_TrackJet_TrkCut_pt1.push_back(N_Jet_TrkCut_pt1[jj]);
    nBranches_->Instanton_vtx_N_goodMuon.push_back(vtx_N_goodMuon[jj]);
    nBranches_->Instanton_vtx_N_goodPFMuon.push_back(vtx_N_goodPFMuon[jj]);
    nBranches_->Instanton_vtx_isBPHtrigger_fromMuon.push_back(vtx_isBPHtrigger_fromMuon[jj]);
    nBranches_->Instanton_vtx_isBPHtrigger_fromPFMuon.push_back(vtx_isBPHtrigger_fromPFMuon[jj]);
    nBranches_->Instanton_vtx_goodMuonIP.push_back(vtx_goodMuonIP[jj]);
    nBranches_->Instanton_vtx_goodPFMuon_PVDistance.push_back(vtx_goodPFMuon_PVDistance[jj]);
  }
  if(doTrack_){
    for( size_t jj = 0; jj < vertices_->size(); ++jj ){
      nBranches_->Trk_px.push_back(Trk_px[jj]);
      nBranches_->Trk_py.push_back(Trk_py[jj]);
      nBranches_->Trk_pz.push_back(Trk_pz[jj]);
      nBranches_->Trk_p0.push_back(Trk_p0[jj]);
      nBranches_->Trk_pt.push_back(Trk_pt[jj]);
      nBranches_->Trk_mass.push_back(Trk_mass[jj]);
      nBranches_->Trk_eta.push_back(Trk_eta[jj]);
      nBranches_->Trk_phi.push_back(Trk_phi[jj]);
      nBranches_->Trk_PVAssociationQuality.push_back(Trk_PVAssociationQuality[jj]);
      nBranches_->Trk_isHighPurity.push_back(Trk_isHighPurity[jj]);
      nBranches_->Trk_vx.push_back(Trk_vx[jj]);
      nBranches_->Trk_vy.push_back(Trk_vy[jj]);
      nBranches_->Trk_vz.push_back(Trk_vz[jj]);
      nBranches_->Trk_IP.push_back(Trk_IP[jj]);
      nBranches_->Trk_IPsignificance.push_back(Trk_IPsignificance[jj]);
      nBranches_->Trk_z0.push_back(Trk_z0[jj]);
      nBranches_->Trk_z0significance.push_back(Trk_z0significance[jj]);
      nBranches_->Trk_matched.push_back(Trk_matched[jj]);
      nBranches_->Trk_pdgId.push_back(Trk_pdgId[jj]);
      nBranches_->Trk_charge.push_back(Trk_charge[jj]);
    }
  }
//  cout<<"finish this event"<<endl;
  return true;
}
