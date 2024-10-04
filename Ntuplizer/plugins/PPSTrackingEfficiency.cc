#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "DataFormats/ProtonReco/interface/ForwardProton.h"
#include "TFile.h"
#include "TH2D.h"
#include <string>
#include <cmath>
#include "../interface/PPSTrackingEfficiency.h"

PPSTrackingEfficiency::PPSTrackingEfficiency(int year)
 : theYear(year) {
  init();
}

std::string PPSTrackingEfficiency::theEra(int runNumber){
                if(runNumber>=297050&& runNumber<=299329) return "2017B";
                if(runNumber>=299368&& runNumber<=300780) return "2017C1";
                if(runNumber>=300806&& runNumber<=302029) return "2017C2";
                if(runNumber>=302031&& runNumber<=302663) return "2017D";
                if(runNumber>=303825&& runNumber<=304797) return "2017E";
                if(runNumber>=305044&& runNumber<=305114) return "2017F1";
                if(runNumber>=305178&& runNumber<=305902) return "2017F2";
                if(runNumber>=305967&& runNumber<=306460) return "2017F3";
                if(runNumber>=306936&& runNumber<=307082) return "2017H";
                if(runNumber>=315252&& runNumber<=316995) return "2018A";
                if(runNumber>=316998&& runNumber<=317696) return "2018B1";
                if(runNumber>=318622&& runNumber<=319312) return "2018B2";
                if(runNumber>=319313&& runNumber<=320393) return "2018C";
                if(runNumber>=320394&& runNumber<=322633) return "2018D1";
                if(runNumber>=323363&& runNumber<=325273) return "2018D2";
                return "";
}

void PPSTrackingEfficiency::init() {
  efficiencyFile_multiRP = new TFile(pixelEfficiency_multiRP.c_str(),"read");
  efficiencyFile_singleRP_far = new TFile(pixelEfficiency_singleRP.c_str(),"read");
  std::vector<std::string> eras;
  if (theYear == 2017){
    efficiencyFile_singleRP_near = new TFile(stripsEfficiency_radiation.c_str(),"read");
    eras = {"2017B","2017C1","2017C2","2017D","2017E","2017F1","2017F2","2017F3","2017H"};
  }
  else if (theYear == 2018){
    efficiencyFile_singleRP_near = new TFile(pixelEfficiency_radiation.c_str(),"read");
    eras = {"2018A","2018B1","2018B2","2018C","2018D1","2018D2"};
  }
  else {
    std::cout
        << "WARNING: no efficiency file for near RP, this might lead to errors";
    return;
  }
  for (auto sector : sectors) {
    std::string rpid = (sector.first == 0) ? "2" : "102";
    std::map<std::string,std::string> pixelString_near;
    std::map<std::string,std::string> pixelString_far;
    std::map<std::string,std::string> stripsString_near;
    std::map<std::string,std::string> stripsString_multiTrack;
    for(std::string era : eras){
      pixelString_near[era] = "Pixel/" + std::to_string(theYear) + "/" +
                                era + "/" + "h" +
                               sector.second + "_210_" +
                               era + "_all_2D";
      pixelString_far[era] = "Pixel/" + std::to_string(theYear) + "/" +
                              era + "/" + "h" +
                              sector.second + "_220_" +
                              era + "_all_2D";
      stripsString_near[era] = "Strips/" + std::to_string(theYear) + "/" +
                                era + "/h" +
                                sector.second + "_" + std::to_string(theYear) +
                                era + "_RP" + rpid + "_all_2D";
      stripsString_multiTrack[era] = "Strips/" + std::to_string(theYear) + "/" +
                                era + "/h" + sector.second + "multitrackeff_" +
                                era + "_avg_RP" + rpid;
    }

    detId_far_ptr = new CTPPSDetId(4, sector.first, 2, 3);
    detId_near_ptr = 0;
    if (theYear == 2017) {
      detId_near_ptr = new CTPPSDetId(3, sector.first, 0, 3);
      bool begin=true;
      for(std::string era : eras){
        if ((TH2D *)efficiencyFile_singleRP_near->Get(stripsString_near[era].c_str())) {
          if(begin)
            efficiencyMaps_singleRP[*detId_near_ptr] = {{era, (TH2D *)efficiencyFile_singleRP_near->Get(stripsString_near[era].c_str())}};
          else
            efficiencyMaps_singleRP[*detId_near_ptr].insert(std::pair<std::string,TH2D *>(era,(TH2D *)efficiencyFile_singleRP_near->Get(stripsString_near[era].c_str())));
        }
        else{
          std::cout << "ERROR: Radiation efficiency map " << stripsString_near[era]
                    << " not found. I will probably crash!" << std::endl;
        }
        if ((TH2D *)efficiencyFile_singleRP_near->Get(stripsString_multiTrack[era].c_str())) {
          double multitrackEff =
              ((TH2D *)efficiencyFile_singleRP_near->Get(stripsString_multiTrack[era].c_str()))
                  ->GetBinContent(1);
          efficiencyMaps_singleRP[*detId_near_ptr][era]->Scale(multitrackEff);
        }
        else{
          std::cout << "ERROR: Multitracking efficiency histogram "
                    << stripsString_multiTrack[era]
                    << " not found. I will probably crash!" << std::endl;
        }
        begin=false;
      }
      inefficiency_shower[*detId_near_ptr] = 0.0;
    }
    else if (theYear == 2018) {
      detId_near_ptr = new CTPPSDetId(4, sector.first, 0, 3);
      bool begin=true;
      for(std::string era : eras){
        if (efficiencyFile_singleRP_near->Get(pixelString_near[era].c_str())) {
          if(begin)
            efficiencyMaps_singleRP[*detId_near_ptr] = {{era, (TH2D *)efficiencyFile_singleRP_near->Get(pixelString_near[era].c_str())}};
          else
            efficiencyMaps_singleRP[*detId_near_ptr].insert(std::pair<std::string,TH2D *>(era,(TH2D *)efficiencyFile_singleRP_near->Get(pixelString_near[era].c_str())));
        }
        else{
          std::cout << "ERROR: Radiation efficiency map " << pixelString_near[era]
                    << " not found. I will probably crash!" << std::endl;
        }
        begin=false;
      }
      inefficiency_shower[*detId_near_ptr] = sector.first?0.017:0.015;
    }
    bool begin=true;
    for(std::string era : eras){
      if (efficiencyFile_multiRP->Get(pixelString_far[era].c_str())) {
        if(begin)
          efficiencyMaps_multiRP[*detId_near_ptr] = {{era, (TH2D *)efficiencyFile_multiRP->Get(pixelString_far[era].c_str())}};
        else
          efficiencyMaps_multiRP[*detId_near_ptr].insert(std::pair<std::string,TH2D *>(era,(TH2D *)efficiencyFile_multiRP->Get(pixelString_far[era].c_str())));
      }
      else{
        std::cout << "ERROR: MultiRP efficiency map " << pixelString_far[era]
                  << " not found. I will probably crash!" << std::endl;
      }
      if (efficiencyFile_singleRP_far->Get(pixelString_far[era].c_str())) {
        if(begin)
          efficiencyMaps_singleRP[*detId_far_ptr] = {{era, (TH2D *)efficiencyFile_singleRP_far->Get(pixelString_far[era].c_str())}};
        else
          efficiencyMaps_singleRP[*detId_far_ptr].insert(std::pair<std::string,TH2D *>(era,(TH2D *)efficiencyFile_singleRP_far->Get(pixelString_far[era].c_str())));
      }
      else{
        std::cout << "ERROR: SingleRP efficiency map " << pixelString_far[era]
                  << " not found. I will probably crash!" << std::endl;
      }
      begin=false;
    }
  }
}

double PPSTrackingEfficiency::getEfficiency(reco::ForwardProton proton, int runNumber) {
    CTPPSLocalTrackLite track = *(proton.contributingLocalTracks().at(0));
    CTPPSDetId detId = CTPPSDetId(track.rpId());
  if (proton.method() == reco::ForwardProton::ReconstructionMethod::singleRP) {
    int bin =
        efficiencyMaps_singleRP[detId][theEra(runNumber)]->FindBin(track.x(), track.y());
    return efficiencyMaps_singleRP[detId][theEra(runNumber)]->GetBinContent(bin);
  } else if (proton.method() ==
             reco::ForwardProton::ReconstructionMethod::multiRP) {
    int bin_near =
        efficiencyMaps_singleRP[detId][theEra(runNumber)]->FindBin(track.x(), track.y());
    double efficiency_near =
        efficiencyMaps_singleRP[detId][theEra(runNumber)]->GetBinContent(bin_near);
    int bin_multiRP =
        efficiencyMaps_multiRP[detId][theEra(runNumber)]->FindBin(track.x(), track.y());
    double efficiency_multiRP =
        efficiencyMaps_multiRP[detId][theEra(runNumber)]->GetBinContent(bin_multiRP);
    return efficiency_near * efficiency_multiRP;
  } else {
    std::cout << "WARNING: proton reconstruction method invalid. Returning "
                 "efficiency = 0"
              << std::endl;
    return 0;
  }
}

double PPSTrackingEfficiency::getEfficiencyunc_stat(reco::ForwardProton proton, int runNumber) {
    CTPPSLocalTrackLite track = *(proton.contributingLocalTracks().at(0));
    CTPPSDetId detId = CTPPSDetId(track.rpId());
  if (proton.method() == reco::ForwardProton::ReconstructionMethod::singleRP) {
    int bin =
        efficiencyMaps_singleRP[detId][theEra(runNumber)]->FindBin(track.x(), track.y());
    return efficiencyMaps_singleRP[detId][theEra(runNumber)]->GetBinError(bin);
  } else if (proton.method() ==
             reco::ForwardProton::ReconstructionMethod::multiRP) {
    int bin_near =
        efficiencyMaps_singleRP[detId][theEra(runNumber)]->FindBin(track.x(), track.y());
    double efficiencyunc_near =
        efficiencyMaps_singleRP[detId][theEra(runNumber)]->GetBinError(bin_near);
    double efficiency_near =
        efficiencyMaps_singleRP[detId][theEra(runNumber)]->GetBinContent(bin_near);
    int bin_multiRP =
        efficiencyMaps_multiRP[detId][theEra(runNumber)]->FindBin(track.x(), track.y());
    double efficiencyunc_multiRP =
        efficiencyMaps_multiRP[detId][theEra(runNumber)]->GetBinError(bin_multiRP);
    double efficiency_multiRP =
        efficiencyMaps_multiRP[detId][theEra(runNumber)]->GetBinContent(bin_multiRP);
    return sqrt(pow(efficiencyunc_near * efficiency_multiRP,2)+pow(efficiency_near * efficiencyunc_multiRP,2));
  } else {
    std::cout << "WARNING: proton reconstruction method invalid. Returning "
                 "efficiency = 0"
              << std::endl;
    return 0;
  }
}

double PPSTrackingEfficiency::getEfficiency_showerineff(reco::ForwardProton proton) {
    CTPPSLocalTrackLite track = *(proton.contributingLocalTracks().at(0));
    CTPPSDetId detId = CTPPSDetId(track.rpId());
    return inefficiency_shower[detId];
}

PPSTrackingEfficiency::~PPSTrackingEfficiency(void){
    inefficiency_shower.clear();
    for(auto &itr: efficiencyMaps_multiRP){
      for(auto &itritr : itr.second)
        delete (itritr.second);
      itr.second.clear();
    }
    efficiencyMaps_multiRP.clear(); 
   for(auto &itr : efficiencyMaps_singleRP){
      for(auto &itritr : itr.second)
        delete (itritr.second);
      itr.second.clear();
    }
    efficiencyMaps_singleRP.clear();
    delete detId_far_ptr;
    delete detId_near_ptr;
    efficiencyFile_multiRP->Close();
    efficiencyFile_singleRP_far->Close();
    efficiencyFile_singleRP_near->Close();
    delete efficiencyFile_multiRP;
    delete efficiencyFile_singleRP_far;
    delete efficiencyFile_singleRP_near;
    sectors.clear();
}
