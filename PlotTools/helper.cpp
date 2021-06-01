#include <iostream>
#include <sstream>
#include <string>
using namespace std;
void splitstring(string origin, vector<string>& split ){
  istringstream iss(origin);
  string word;
  while (iss >> word){
    cout<<word<<"\n";
    split.push_back(word);
  }
  cout<<endl;
}
void rootfilelist(string origin, vector<string>& list){
  vector<string> split;
  splitstring(origin,split);
  for(int i=0; i<split.size(); i++){
    if(split[i].find(".root")==std::string::npos){
      FILE *proc = popen(("/bin/ls -1 "+split[i]+"/*.root").c_str(),"r");
      char buf[1024];
      while ( !feof(proc) && fgets(buf,sizeof(buf),proc) )
      {
//        printf("Line read: %s",buf);
          std::string s(buf);
          list.push_back(s.substr(0,s.size()-1));
          cout<<s<<endl;
      }
    }
    else{
      list.push_back(split[i]);
      cout<<split[i]<<endl;
    }
  }
}

void HistError(TChain * tree,string PUWeight_tree, string Branch_v, string Branch_error,string selection, TH1F * Hist_origin, int nbins, float min, float max){
  tree->AddFriend("tree",PUWeight_tree.c_str());
  TH1F *hist = new TH1F("hist","hist",nbins,min,max);
  hist->Sumw2(kTRUE);
  tree->Draw((Branch_v+">>hist").c_str(),("("+selection+")*"+Branch_error).c_str(),"goff");
  for (Int_t i=0;i<=nbins;i++) {
    Hist_origin->SetBinError(i,sqrt(Hist_origin->GetBinError(i)*Hist_origin->GetBinError(i)+hist->GetBinError(i)*hist->GetBinError(i)));
  }
}
