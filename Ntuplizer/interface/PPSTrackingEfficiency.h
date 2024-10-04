#ifndef PPSTrackingEfficiency_H
#define PPSTrackingEfficiency_H

#include <memory>

#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "DataFormats/ProtonReco/interface/ForwardProton.h"
#include "TFile.h"
#include "TH2D.h"
#include <string>
#include <cmath>
class PPSTrackingEfficiency {
public:
  PPSTrackingEfficiency(int year);
  ~PPSTrackingEfficiency();
  double getEfficiency(reco::ForwardProton proton, int runNumber);
  double getEfficiencyunc_stat(reco::ForwardProton proton, int runNumber);
  double getEfficiency_showerineff(reco::ForwardProton proton);
private:
  void init();
  int theYear;
  TFile* efficiencyFile_multiRP;
  TFile* efficiencyFile_singleRP_far;
  TFile* efficiencyFile_singleRP_near;
  std::map<int, std::string> sectors = {{0, "45"}, {1, "56"}};
  CTPPSDetId* detId_far_ptr;
  CTPPSDetId* detId_near_ptr;
  std::map<CTPPSDetId,double> inefficiency_shower;
  std::map<CTPPSDetId, std::map<std::string,TH2D*>> efficiencyMaps_multiRP;
  std::map<CTPPSDetId, std::map<std::string,TH2D*>> efficiencyMaps_singleRP;

  std::string pixelEfficiency_singleRP =
      "./data/ctpps/"
      "pixelEfficiencies_singleRP_reMiniAOD.root";
  std::string pixelEfficiency_multiRP =
      "./data/ctpps/"
      "pixelEfficiencies_multiRP_reMiniAOD.root";
  std::string pixelEfficiency_radiation =
      "./data/ctpps/"
      "pixelEfficiencies_radiation.root";
  std::string stripsEfficiency_radiation =
      "./data/ctpps/"
      "PreliminaryEfficiencies_July132020_1D2DMultiTrack.root";
  std::string theEra(int runNumber);
};
#endif
