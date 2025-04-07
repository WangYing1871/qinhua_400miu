//--------------------------------Stamp-------------------------------
//^-^ Author: Zhi Heng            Email: wy187110@mail.ustc.edu.cn     
//^-^ Time: 2024-12-21 16:05:28   Posi: Hefei
//^-^ File: time_sort.cpp
//--------------------------------------------------------------------
#ifndef info_out
# define info_out(X) std::cout<<"==> "<<__LINE__<<" "<<#X<<" |"<<(X)<<"|\n"
#endif
#define OFFSET(X,Y) std::next(std::begin(X),Y)
#define ROFFSET(X,Y) std::prev(std::end(X),Y)
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <unordered_map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <filesystem>
#include <ctime>

#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TFolder.h"

#include "clusters.hpp"
#include "data_strcut_cint.h"
namespace util{
std::string trim_first_digits(std::string const& str){
  auto iter = std::find_if(std::begin(str),std::end(str),[](auto a){return std::isdigit(a);});
  if (iter != std::end(str)){
    std::string ret("");
    for(auto it = iter; it != std::end(str); it++){
      if (!std::isdigit(*it)) return ret; 
      ret.push_back(*it); }
  }
  std::cerr<<"[Warning] No time_stamp info get. use now time point\n";
  return std::to_string(std::time(0));
}

}
namespace helper{
int re_channel_id(int in){
  return in<11 ? in :
    in>=11 && in<22 ? in+1 : 
     in>=22 && in<45 ? in+2 :
        in>=45 && in<56 ? in+3 :
          in>=56 ? in+4 : in;
}
}

namespace wave_handler{
struct wave_t{
  
};

struct wave_filter{

};
}

struct aget_raw_data_t{
  static constexpr const std::size_t s_max_channels = 256;
  long time_stamp;
  int trigger_id;
  int nhits;
  int chip_id[s_max_channels];
  int chn_id[s_max_channels];
  //int adc_data[s_max_channels][1024];
  int adc_data[s_max_channels][512];

  uint8_t is_valid[s_max_channels];

  void clear(){

  }
};

// Temporary solution //FIXME 2024-12-26 16:28:16 by Zhi Heng 
struct layer_store_t{
  typedef layer_store_t self_t;
  std::string m_name;
  TFile* m_file;
  TTree* m_data_tree;
  aget_raw_data_t m_data_buf;

public:
  layer_store_t()=default;
  ~layer_store_t() noexcept = default;
  layer_store_t(std::string const& v):m_name(v) {}



  template <class _tp>
  self_t* set(_tp* v){
    if constexpr(std::is_same<TTree,_tp>::value){ m_data_tree = v; return this; }
    if constexpr(std::is_same<TFile,_tp>::value){ m_file = v; return this; }
    return this;
  }
  void write_tree(){
    if (!m_file || !m_data_tree) return;
    m_file->cd(); m_data_tree->Write(); }

  aget_raw_data_t& get_store() &{return m_data_buf; }
  void fill() {if (m_data_tree) m_data_tree->Fill();}
  void cd() {if (m_file) m_file->cd();}
  void write_close() { write_tree(); if (m_file) m_file->Write(), m_file->Close();}

  bool draw_wave(std::string const& folder_name){
    //auto* folder = new TFolder(folder_name.c_str(),folder_name.c_str());
    if (!m_file || m_file->IsZombie()) return false;
    if (m_data_buf.nhits==0) return false;
    auto folder = std::make_shared<TFolder>(folder_name.c_str(),folder_name.c_str());
    for (int i=0; i<m_data_buf.nhits; ++i){
      std::stringstream sstr("");
      sstr<<m_data_buf.chip_id[i]<<"-"<<m_data_buf.chn_id[i];
      auto* graph = new TGraph();
      graph->SetName(sstr.str().c_str());
      graph->SetTitle(m_data_buf.is_valid[i] ? "O" : "X");
      //graph->SetDirectory(nullptr);
      for (std::size_t index=0; auto&& x : m_data_buf.adc_data[i]) {
        graph->SetPoint(index,index,x), index++;
      }
      folder->Add(graph);
    }
    m_file->cd();
    folder->Write();
    return true;
  }

  std::ostream& display(std::ostream& os = std::cout) const{
    info_out(m_data_buf.trigger_id);

    for (int i=0; i<m_data_buf.nhits; ++i){
      for (int j=0; j<512; ++j)
        os << m_data_buf.adc_data[i][j]<<std::endl;
      info_out("");
    }
    return os;
  }
};

layer_store_t make_store(std::string const& name
    ,std::string const& fname, std::string const& tname){
  layer_store_t s(name);
  s.m_file = new TFile(fname.c_str(),"recreate");
  s.m_data_tree = new TTree(tname.c_str(),tname.c_str());
  s.m_data_tree->Branch("time_stamp",std::addressof(s.m_data_buf.time_stamp)
      ,"s.m_data_buf.time_stamp/L");
  s.m_data_tree->Branch("trigger_id",std::addressof(s.m_data_buf.trigger_id)
      ,"s.m_data_buf.trigger_id/I");
  s.m_data_tree->Branch("nhits",std::addressof(s.m_data_buf.nhits),"s.m_data_buf.nhits/I");
  s.m_data_tree->Branch("chip_id",s.m_data_buf.chip_id,"chip_id[s.m_data_buf.nhits]/I");
  s.m_data_tree->Branch("chn_id",s.m_data_buf.chn_id,"chn_id[s.m_data_buf.nhits]/I");
  s.m_data_tree->Branch("adc_data",s.m_data_buf.adc_data,"adc_data[s.m_data_buf.nhits][512]/I");
  //s.m_data_tree->Branch("adc_data",s.m_data_buf.adc_data,"adc_data[s.m_data_buf.nhits][1024]/I");
  return s;
}

//void save_wave
typedef std::map<uint16_t,uint16_t> dcm2fec_map_t;
void read_map(std::string const fname, dcm2fec_map_t& map){
  namespace fs = std::filesystem;
  map.clear();
  if (!fs::exists(fname)){
    throw std::invalid_argument("Invalid file name"); }
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



bool wave_filter(std::vector<uint16_t> const& adcs  
    ,std::pair<float,float> ms, int scale, int comp0, int comp1,int width=20){
  auto min_iter = std::min_element(adcs.begin(),adcs.end());
  auto max_iter = std::max_element(adcs.begin(),adcs.end());
  if (*min_iter<ms.first-comp0*ms.second) return false;
  if (*max_iter<ms.first+comp1*ms.second) return false;
  int distance = (int)std::distance(adcs.begin(),max_iter);
  if (distance<width || distance>adcs.size()-width) return false;
  uint32_t sum_max=0;
  //info_out(distance); info_out(width);
  size_t beg = distance-width<=0 ? 0 : distance-width;
  size_t end = distance+width>adcs.size() ? adcs.size(): distance+width;
  size_t range = end-beg;
  //info_out(beg); info_out(end);
  for (size_t i=beg; i<end; ++i){
    sum_max += adcs[i];
  }
  float peak_sz = (float)sum_max/range;
  //info_out(max_iter);
  if (peak_sz<ms.first+scale*ms.second) return false;
  if (peak_sz<0.8**max_iter) return false;
  return true;
}
int main(int argc, char* argv[]){

  std::string v1_map_name="calc/v1_map.txt";
  dcm2fec_map_t v1_map;
  std::string v2_map_name="calc/v2_map.txt";
  dcm2fec_map_t v2_map;
  read_map(v1_map_name,v1_map);
  read_map(v2_map_name,v2_map);



  std::string fin_name = argv[1];
  std::string baseline_name = argv[2];
  auto* fin = new TFile(fin_name.c_str());
  auto ts_str = util::trim_first_digits(fin_name.substr(fin_name.find_last_of('/')+1));
  auto nm0  = "run_"+ts_str+"out_mt-0.root";
  auto nm1  = "run_"+ts_str+"out_mt-1.root";
  //info_out(nm0);
  //exit(0);
  namespace fs = std::filesystem;
  if (!fs::exists(fin_name) || !fs::is_regular_file(fin_name))
    throw std::invalid_argument("invalid input root file");
  auto* tree = (TTree*)fin->Get("CollectionTree");
  entry_new* data = new entry_new;
  tree->SetBranchAddress("data",std::addressof(data));
  auto entries = tree->GetEntries();

  //auto* rfout = new TFile("time_sort.root","recreate");
  //auto* tree_out = new TTree("fec_origin_data","fec_origin_data");
  //aget_raw_data_t store_buf;
  //tree_out->Branch("time_stamp",std::addressof(store_buf.time_stamp),"store_buf.time_stamp/L");
  //tree_out->Branch("trigger_id",std::addressof(store_buf.trigger_id),"store_buf.trigger_id/I");
  //tree_out->Branch("nhits",std::addressof(store_buf.nhits),"store_buf.nhits/I");
  //tree_out->Branch("chip_id",store_buf.chip_id,"chip_id[store_buf.nhits]/I");
  //tree_out->Branch("chn_id",store_buf.chn_id,"chn_id[store_buf.nhits]/I");
  //tree_out->Branch("adc_data",store_buf.adc_data,"adc_data[store_buf.nhits][1024]/I");
 // std::unordered_map<std::string,TTree*>

  std::map<std::pair<int,int>,std::pair<float,float>> baseline_map;
  std::ifstream prestal_fin(baseline_name.c_str());
  while(!prestal_fin.eof()){
    std::string sbuf;
    std::getline(prestal_fin,sbuf);
    if(!sbuf.empty() && sbuf[0] != '#'){
      std::stringstream sstr(sbuf.c_str());
      int a, d; float b, c,e;
      sstr>>a>>d>>b>>c>>e;
      baseline_map[std::make_pair(a,d)]=std::make_pair(b,c);
    }
  }

  typedef struct{
    //uint16_t max;
    uint16_t adcs[1024];
    uint64_t evt_id;
    int glb_id;
    bool is_valid = true;
  }max_id_t;
  std::multimap<uint64_t,max_id_t> tm_vs_adc;
  //double ww=5000;
  double ww=500;
  std::multiset<std::size_t> hits_no_map;
  int trigger_id_idx=-1;
  int ww0=0, ww1=0;
#ifdef DEBUG
  std::vector<std::vector<int>> strip_ids;

#endif
  




  auto storeL0 = make_store("layer0",nm0,"fec_origin_data");
  auto storeL1 = make_store("layer1",nm1,"fec_origin_data");
  std::size_t total_wv = 0;
  std::size_t valid_wv = 0;

  std::size_t min_strip_id=3;
  auto const& store = [&](){
    if (tm_vs_adc.size()<=10) return;
    auto& store_buf = storeL0.get_store();
    auto& store1_buf = storeL1.get_store();
    auto iter_front = std::begin(tm_vs_adc);
    for (auto iter=tm_vs_adc.begin(); iter != ROFFSET(tm_vs_adc,1);){
      auto iter_next = std::next(iter);
      if (iter_next->first-iter->first>ww
          ){
        if (true //TODO
            && 
            (std::distance(iter_front,iter_next)>aget_raw_data_t::s_max_channels
            || std::distance(iter_front,iter_next)<min_strip_id)
            )
        {
          info_out(std::distance(iter_front,iter_next));
          iter = iter_front = iter_next;
          continue;
        }
        store_buf.nhits=0;
        store1_buf.nhits=0;
        //info_out(std::distance(iter_front,iter_next));

        for(auto iter_tmp = iter_front; iter_tmp != iter_next; ++iter_tmp){
          auto id = iter_tmp->second.glb_id;
          uint8_t adm_layer_id = id>>8;
          uint8_t adm_channel_id = id&0xFF;

          if (adm_layer_id==0){
            if(v1_map.find(adm_channel_id)!=v1_map.end()){
              store_buf.chip_id[store_buf.nhits]=0;
              store_buf.chn_id[store_buf.nhits]=helper::re_channel_id(v1_map[adm_channel_id]);
              for (int i=0; i<1024; i+=2)
                store_buf.adc_data[store_buf.nhits][i/2] = 
                  (iter_tmp->second.adcs[i]+iter_tmp->second.adcs[i+1])/2;
              store_buf.is_valid[store_buf.nhits] = iter_tmp->second.is_valid;
              store_buf.nhits++;
            }
            else if(v2_map.find(adm_channel_id)!=v2_map.end()){
              store_buf.chip_id[store_buf.nhits]=2;
              store_buf.chn_id[store_buf.nhits]=helper::re_channel_id(v2_map[adm_channel_id]);
              for (int i=0; i<1024; i+=2)
                store_buf.adc_data[store_buf.nhits][i/2] = 
                  (iter_tmp->second.adcs[i]+iter_tmp->second.adcs[i+1])/2;
              store_buf.is_valid[store_buf.nhits] = iter_tmp->second.is_valid;
              store_buf.nhits++;
            }
          }else if(adm_layer_id==1){
            if(v1_map.find(adm_channel_id)!=v1_map.end()){
              store_buf.chip_id[store_buf.nhits]=1;
              store_buf.chn_id[store_buf.nhits]=helper::re_channel_id(v1_map[adm_channel_id]);
              for (int i=0; i<1024; i+=2)
                store_buf.adc_data[store_buf.nhits][i/2] = 
                  (iter_tmp->second.adcs[i]+iter_tmp->second.adcs[i+1])/2;
              store_buf.is_valid[store_buf.nhits] = iter_tmp->second.is_valid;
              store_buf.nhits++;
            }
            else if(v2_map.find(adm_channel_id)!=v2_map.end()){
              store_buf.chip_id[store_buf.nhits]=3;
              store_buf.chn_id[store_buf.nhits]=helper::re_channel_id(v2_map[adm_channel_id]);
              for (int i=0; i<1024; i+=2)
                store_buf.adc_data[store_buf.nhits][i/2] = 
                  (iter_tmp->second.adcs[i]+iter_tmp->second.adcs[i+1])/2;
              store_buf.is_valid[store_buf.nhits] = iter_tmp->second.is_valid;
              store_buf.nhits++;
            }
          }
          
          if (adm_layer_id==2){
            if(v1_map.find(adm_channel_id)!=v1_map.end()){
              store1_buf.chip_id[store1_buf.nhits]= 0;
              store1_buf.chn_id[store1_buf.nhits]= helper::re_channel_id(v1_map[adm_channel_id]);
              for (int i=0; i<1024; i+=2)
                store1_buf.adc_data[store1_buf.nhits][i/2] = 
                  (iter_tmp->second.adcs[i]+iter_tmp->second.adcs[i+1])/2;
              store1_buf.nhits++;
            }
            else if(v2_map.find(adm_channel_id)!=v2_map.end()){
              store1_buf.chip_id[store1_buf.nhits]= 2;
              store1_buf.chn_id[store1_buf.nhits]= helper::re_channel_id(v2_map[adm_channel_id]);
              for (int i=0; i<1024; i+=2)
                store1_buf.adc_data[store1_buf.nhits][i/2] = 
                  (iter_tmp->second.adcs[i]+iter_tmp->second.adcs[i+1])/2;
              store1_buf.nhits++;
            }
          }
          else if(adm_layer_id==3){
            if(v1_map.find(adm_channel_id)!=v1_map.end()){
              store1_buf.chip_id[store1_buf.nhits]= 1;
              store1_buf.chn_id[store1_buf.nhits]= helper::re_channel_id(v1_map[adm_channel_id]);
              for (int i=0; i<1024; i+=2)
                store1_buf.adc_data[store1_buf.nhits][i/2] = 
                  (iter_tmp->second.adcs[i]+iter_tmp->second.adcs[i+1])/2;
              store1_buf.nhits++;
            }
            else if(v2_map.find(adm_channel_id)!=v2_map.end()){
              store1_buf.chip_id[store1_buf.nhits]= 3;
              store1_buf.chn_id[store1_buf.nhits]= helper::re_channel_id(v2_map[adm_channel_id]);
              for (int i=0; i<1024; i+=2)
                store1_buf.adc_data[store1_buf.nhits][i/2] = 
                  (iter_tmp->second.adcs[i]+iter_tmp->second.adcs[i+1])/2;
              store1_buf.nhits++;
            }
          }
        }

        if (store_buf.nhits>0 || store1_buf.nhits>0){
        //if (store_buf.nhits>0 && store1_buf.nhits>0){ //Just for test
          trigger_id_idx++;
          store_buf.time_stamp = std::find_if(
              iter_front
              ,iter_next
              ,[](auto v){
                auto id = v.second.glb_id;
                uint8_t adm_layer_id = id>>8;
                if (adm_layer_id==0 || adm_layer_id==1) return true;
                return false;
              }
              )->first;
          store_buf.trigger_id = trigger_id_idx;
          storeL0.fill();


          store1_buf.time_stamp = std::find_if(
              iter_front
              ,iter_next
              ,[](auto v){
                auto id = v.second.glb_id;
                uint8_t adm_layer_id = id>>8;
                if (adm_layer_id==2 || adm_layer_id==3) return true;
                return false;
              }
              )->first;
          store1_buf.trigger_id = trigger_id_idx;
          storeL1.fill();

          //if (ww0<200 || ww1<200){
          //  std::stringstream sstr(""); sstr<<"Event-"<<trigger_id_idx;
          //  if (storeL0.draw_wave(sstr.str())) ww0++;
          //  if (storeL1.draw_wave(sstr.str())) ww1++;
          //}
        }
        iter = iter_front = iter_next;
      }else iter++;



    }
    tm_vs_adc.clear();
  };



  info_out(entries);
  for(int i=0; i<entries; ++i){
    if (i%1000==0) info_out("progress");
    tree->GetEntry(i);
    std::map<int,uint64_t> ts_map;
    for(std::size_t index=0; auto&& x : data->fec_ids) ts_map[x] = data->time_stamps[index++];
    for(std::size_t index=0; auto&& x : data->global_ids){
      auto const& adcs = data->adcs[index++];
      int fec_id = x>>8;
      int channel_id = x&0xFF;
      auto ms = baseline_map.at(std::make_pair(fec_id,channel_id));
      if (ts_map.find(fec_id)==ts_map.end()) continue;
      uint64_t ts = ts_map.at(fec_id);
      auto max = std::max_element(std::begin(adcs),std::end(adcs));
      uint16_t max_position = std::distance(std::begin(adcs),max);
      uint64_t ts_peak = std::round(ts*8.33) + (max_position-624)*25;
      max_id_t max_id;
      //max_id.max = *max;
      std::copy(OFFSET(adcs,0),ROFFSET(adcs,0),OFFSET(max_id.adcs,0));
      max_id.glb_id = x;
      max_id.evt_id = i;
      max_id.is_valid = wave_filter(adcs,ms,4,5,8,20);
      total_wv++;
      if (max_id.is_valid) valid_wv++;
      //////if (is_valid) save_wave(adcs,"1/channel0");
      tm_vs_adc.emplace(ts_peak,max_id);
    }
    if (tm_vs_adc.size()>5000){
      store();
    }
  }


  //for (auto iter = OFFSET(tm_vs_adc,0); iter != OFFSET(tm_vs_adc,10); ++iter){

  //  info_out("CCCC");
  //  for (auto&& x : iter->second.adcs) std::cout<<x<<" ";
  //  info_out("");

  //}
  store();
  info_out(total_wv);
  info_out(valid_wv);

  //info_out(tm_vs_adc.size());
  //info_out(std::count_if(
  //      std::begin(tm_vs_adc),std::end(tm_vs_adc),[]<class _tp>(_tp v){return v.second.is_valid;}));


  fin->Close();

  storeL0.write_close();
  storeL1.write_close();

  //rfout->cd();
  //tree_out->Write();
  //rfout->Write();
  //rfout->Close();

  
  //{
  //  uint64_t prev_value = 0;
  //  std::multiset<uint64_t> tm_dis;
  //for (auto&& [x,y] : tm_vs_adc){
  // // std::cout<<x<<" "<<x-prw
  //  if (prev_value==0){ prev_value = x; continue;}
  //  tm_dis.insert(x-prev_value); prev_value = x;
  //}
  //for (auto x : tm_dis)
  //  std::cout<<x<<"\n";
  //}
  {
    //for (auto&&x : hits_no_map) info_out(x);
  }
  return 0;
}
