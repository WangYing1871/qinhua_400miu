//--------------------------------Stamp-------------------------------
//^-^ Author: Zhi Heng            Email: wy187110@mail.ustc.edu.cn     
//^-^ Time: 2024-12-30 10:29:37   Posi: Hefei
//^-^ File: noise.cpp
//--------------------------------------------------------------------
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <filesystem>
#include <stdexcept>
#include <map>
#include <memory>
#include <set>

#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1I.h"
#include "TGraph.h"
#include "TF1.h"

#include "util.hpp"
#include "data_strcut_cint.h"

constexpr static const std::size_t c_l= 4;
constexpr static const std::size_t c_c= 64;
struct noise_store_t{
  long counts=0;
  double mean[c_l][c_c];
  double rms[c_l][c_c];

  void clear(){
    counts=0;
    std::memset(mean,0.,sizeof(mean));
    std::memset(rms,0.,sizeof(rms));
  }
};
typedef std::map<uint16_t,uint16_t> dcm2fec_map_t;
void read_map(std::string const fname, dcm2fec_map_t& map){
  namespace fs = std::filesystem;
  map.clear();
  if (!fs::exists(fname)){
    std::stringstream sstr("");
    sstr<<"invalid file name: ["<<fname<<"], please confirm your path";
    throw std::invalid_argument(sstr.str().c_str()); }
  std::ifstream fin(fname.c_str());
  while(!fin.eof()){
    std::string sbuf; std::getline(fin,sbuf);
    if (!sbuf.empty() && sbuf[0] != '#'){
      std::stringstream sstr(sbuf.c_str());
      uint16_t a,b;
      sstr>>a>>b;
      map[a] = b;
    }
  }
 fin.close();
}

void draw_wave(entry_new const* data, std::string const& name){
  //save_wave<TGraph>("event1/wave10",)
}
using std::cout; using std::endl; using std::string; using std::vector; 
int main(int argc, char* argv[]){
  if (argc<2){ std::cerr<<"input file name needed\n"; return 0; }
  std::string v1_map_name="calc/v1_map.txt";
  dcm2fec_map_t v1_map;
  std::string v2_map_name="calc/v2_map.txt";
  dcm2fec_map_t v2_map;
  read_map(v1_map_name,v1_map);
  read_map(v2_map_name,v2_map);

  info_out(v1_map.size());
  info_out(v2_map.size());

  auto* fin = new TFile(argv[1]);
  auto* tree = (TTree*)fin->Get("CollectionTree");
  entry_new* data = new entry_new;
  tree->SetBranchAddress("data",std::addressof(data));
  auto entries = tree->GetEntries();
  info_out(entries);

  

  auto* rfout0 = new TFile("noiseL0.root","recreate");
  auto* tree0_out = new TTree("noise_tree","noise_tree");
  noise_store_t noise0_store;
  tree0_out->Branch("total_cnt",&noise0_store.counts,"&noise0_store.count/L");
  tree0_out->Branch("mean",noise0_store.mean,"noise0_store.mean[4][64]/D");
  tree0_out->Branch("rms",noise0_store.rms,"noise0_store.rms[4][64]/D");

  //TODO factory need!
  // functional
  // register
  // manager
  auto* rfout1 = new TFile("noiseL1.root","recreate");
  auto* tree1_out = new TTree("noise_tree","noise_tree");
  noise_store_t noise1_store;
  tree1_out->Branch("total_cnt",&noise1_store.counts,"&noise1_store.count/L");
  tree1_out->Branch("mean",noise1_store.mean,"noise1_store.mean[4][64]/D");
  tree1_out->Branch("rms",noise1_store.rms,"noise1_store.rms[4][64]/D");

  //std::unordered_map<uint16_t,TH1I*> baseline_his;
  std::unordered_map<uint16_t,std::shared_ptr<TH1I>> baseline_his;
  for(int i=0; i<entries; ++i){
    tree->GetEntry(i);
    if (i>200){
      info_out("[WARN] too much entries for baseline meaningless. 200 used");
      break; }
    for(std::size_t index=0; auto&& x : data->global_ids){
      auto const& adcs = data->adcs[index++];
      if (baseline_his.find(x)==baseline_his.end()){
        std::stringstream sstr("");
        sstr<<"baseline-"<<x;
        //baseline_his.emplace(x,new TH1I(sstr.str().c_str(),sstr.str().c_str(),300,500,800));
        baseline_his.emplace(x,std::make_shared<TH1I>(sstr.str().c_str(),sstr.str().c_str(),500,500,1000));
        baseline_his[x].get()->SetDirectory(nullptr);
      }
      for (auto&& y : adcs){
        baseline_his[x].get()->Fill(y);
      }
      //int fec_id = x>>8;
      //int channel_id = x&0xFF;
      //switch(fec_id){
      //  case 0:
      //    if (v1_map.find(channel_id)!=v1_map.end()){
      //      
      //    }
      //    break;
      //}
    }
  }

  auto const& remap = [&](uint16_t gid)->std::pair<std::size_t,std::size_t>{
    constexpr static const std::size_t max_sz = (std::numeric_limits<std::size_t>::max)();
    uint8_t lid= gid>>8;
    uint8_t cid = gid&0xFF;
    switch(lid){
      case 0:
        if (v1_map.find(cid)!=v1_map.end()) return std::make_pair(0,v1_map[cid]);
        if (v2_map.find(cid)!=v2_map.end()) return std::make_pair(2,v2_map[cid]);
        break;
      case 1:
        if (v1_map.find(cid)!=v1_map.end()) return std::make_pair(1,v1_map[cid]);
        if (v2_map.find(cid)!=v2_map.end()) return std::make_pair(3,v2_map[cid]);
        break;
      case 2:
        if (v1_map.find(cid)!=v1_map.end()) return std::make_pair(0,v1_map[cid]);
        if (v2_map.find(cid)!=v2_map.end()) return std::make_pair(2,v2_map[cid]);
        break;
      case 3:
        if (v1_map.find(cid)!=v1_map.end()) return std::make_pair(1,v1_map[cid]);
        if (v2_map.find(cid)!=v2_map.end()) return std::make_pair(3,v2_map[cid]);
        break;
      default:
        break;
    }
    return std::make_pair(max_sz,max_sz);
  };

  noise0_store.clear();

  std::multiset<double> view_set;
  for (auto&& [x,y] : baseline_his){
    TF1 f_gaus("f_gaus","gaus"
        ,y.get()->GetMean()-3*y.get()->GetRMS()
        ,y.get()->GetMean()+3*y.get()->GetRMS()
        );
    y.get()->Fit(&f_gaus,"RQ");
    if ((x>>8)<2){
      auto idx = remap(x);
      noise0_store.mean[idx.first][idx.second]=f_gaus.GetParameter(1);
      noise0_store.rms[idx.first][idx.second]=f_gaus.GetParameter(2);
      noise0_store.counts++;

    }
    if ((x>>8)>=2){
      auto idx = remap(x);
      noise1_store.mean[idx.first][idx.second]=f_gaus.GetParameter(1);
      noise1_store.rms[idx.first][idx.second]=f_gaus.GetParameter(2);
      noise1_store.counts++;
    }
    

    
  }
  for(auto&& x : view_set)
    info_out(x);

  fin->Close();
  tree0_out->Fill();
  rfout0->Write();
  rfout0->Close();
  tree1_out->Fill();
  rfout1->Write();
  rfout1->Close();

  //auto* his_file = new TFile("his.root","recreate");
  //his_file->cd();
  //his_file->Write();
  //his_file->Close();

  return 0;
}
