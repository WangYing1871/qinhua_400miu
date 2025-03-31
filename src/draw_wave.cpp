//--------------------------------Stamp-------------------------------
//^-^ Author: Zhi Heng            Email: wy187110@mail.ustc.edu.cn     
//^-^ Time: 2024-10-31 15:55:10   Posi: Hefei
//^-^ File: draw_wave.cpp
//--------------------------------------------------------------------
#ifndef info_out
#define info_out(X) std::cout<<"==> "<<__LINE__<<" "<<#X<<" |"<<(X)<<"|\n"
#endif
#include <algorithm>
#include <numeric>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include <filesystem>

#include "TTree.h"
#include "TFile.h"
#include "TFolder.h"
#include "TGraph.h"
#include "TH1F.h"
#include "data_strcut_cint.h"

//std::map<int,std::vector<uint16_t>> pix_map;

struct point_t{
  float x;
  float y;
  float z;
  //point_t(float a, float b, float c):x(a),y(b),z(c) {}
};

TFolder* waves_view(std::string const& name, entry_new* data){
  auto* folder = new TFolder(name.c_str(),name.c_str());
  auto& ids = data->global_ids;
  auto& adcs = data->adcs;

  std::map<std::string,TGraph*> grps;
  for (int i=0; i<ids.size(); ++i){
    std::stringstream sstr("");
    auto layer_id = uint8_t(ids[i]>>8);
    auto channel_id = uint8_t(ids[i]&0xFF);
    sstr<<"wave"<<(int)layer_id<<"-"<<(int)channel_id;
    auto* graph = new TGraph();
    graph->SetName(sstr.str().c_str());
    graph->SetTitle(sstr.str().c_str());
    auto max_adc = *std::max_element(std::begin(adcs[i]),std::end(adcs[i]));
   // pix_map[ids[i]].emplace_back(max_adc);
    for (std::size_t index=0; auto&& x : adcs[i])
      graph->SetPoint(index,index,x), index++;
    grps[sstr.str()] = graph;
    //folder->Add(graph);
  }
  for (auto&& [x,y] : grps) folder->Add(y);
  return folder;
}

template <class _iter>
void draw_graph(TGraph& g, std::string const& name, _iter begin, _iter end){
  g.SetName(name.c_str());
  g.SetTitle(name.c_str());
  int distance = std::distance(begin,end);
  for (int i=0; i<distance; ++i)
    g.SetPoint(i,i,*std::next(begin,i));
}

void demo00(uint64_t evt_id, entry_tm* data){
}


using std::cout; using std::endl; using std::string; using std::vector;
int main(int argc, char* argv[]){
  {
    std::string fname = argv[1];
    std::string fout_name = "wave.root";
    if (argc>=3) fout_name = argv[2];
    auto* fin = new TFile(fname.c_str());
    auto* tree = (TTree*)fin->Get("CollectionTree");
    entry_new* data = new entry_new;
    tree->SetBranchAddress("data",std::addressof(data));

    auto* fout = new TFile(fout_name.c_str(),"recreate");
    fout->cd();

    auto entries = tree->GetEntries();

    //for (int i=0; i<tree->GetEntries(); ++i){
    for (int i=0; i<(entries>200 ? 200 : entries); ++i){
    //for (int i=200; i<300; ++i){
      tree->GetEntry(i);
      std::stringstream sstr("");
      sstr<<"Event-"<<i;
      waves_view(sstr.str(),data)->Write();
    }
    fout->Write(); fout->Close();
    fin->Close();
    return 0;
  }
  //{
  //  TFile* file = new TFile("temp.root","recreate");
  //  file->cd();
  //  for (int i=0; i<100000; ++i){
  //    std::stringstream sstr;
  //    sstr<<"his"<<i;
  //    TH1F f1(sstr.str().c_str(),sstr.str().c_str(),1000,-5,5);
  //    f1.FillRandom("gaus",10000);
  //    f1.Write();
  //    if (i%1000==0){
  //      info_out(std::filesystem::file_size("temp.root"));
  //    }
  //  }
  //  file->Write();
  //  file->Close();
  //}
  {
    std::string fname = argv[1];
    auto* fin = new TFile(fname.c_str());
    auto* tree = (TTree*)fin->Get("CollectionTree");
    auto* data = new entry_tm;
    tree->SetBranchAddress("data",std::addressof(data));
    info_out(tree->GetEntries());
    auto* fout = new TFile("wave.root","recreate");
    fout->cd();
    std::map<std::string,TGraph> grp_map;
    for (long long i=0; i<tree->GetEntries(); ++i){
      grp_map.clear();
      tree->GetEntry(i);
      std::stringstream sstr("");
      sstr<<"Event"<<i;
      TFolder folder(sstr.str().c_str(),sstr.str().c_str());
      for (std::size_t index=0; auto&& x : data->det_ids){
        sstr = std::stringstream("");
        sstr<<"wave-"<<(x>>16)<<"-"<<(x&0xFFFF);
        auto& adc = data->adcs[index++];
        TGraph g;
        draw_graph(g,sstr.str(),std::begin(adc),std::end(adc));
        //g.Write();
        grp_map.emplace(sstr.str(),g);
      }
      for (auto&& [x,y] : grp_map) folder.Add(&y);
      folder.Write();
      if (i%1000==0)
        info_out(std::filesystem::file_size("wave.root"));
    }
    fout->Write();
    fout->Close();



  }
  return 0;
}
