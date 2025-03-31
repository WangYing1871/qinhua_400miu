#include <map>
#include <set>
#include <future>
#include <stack>
#include "TTree.h"
#include "TFile.h"
#include "TH1I.h"
#include "TF1.h"

#include "util.hpp"
#include "data_strcut_cint.h"
#include "unpack.h"

namespace{

//template <class _tp>


struct event_record_t{
  uint32_t m_evt_id;
  std::set<uint8_t> m_boards_id;

public:
  event_record_t(uint32_t a, uint8_t b):
    m_evt_id(a) {m_boards_id.insert(b);}
  ~event_record_t() noexcept = default;

public:
  void clear() {m_evt_id = 0; m_boards_id.clear();}
  void insert_boards(uint8_t v) {m_boards_id.insert(v);}
  unsigned long erase_board(uint8_t v) {return m_boards_id.erase(v);}
  bool has_board(uint8_t v) const {return m_boards_id.find(v)!=m_boards_id.end();}
  bool empty() const {return m_boards_id.empty();}
  uint32_t event_id() const {return m_evt_id;}
#ifdef DEBUG
public:
  std::ostream& display(std::ostream& os=std::cout){
    os<<"event_id: "<<m_evt_id<<" [";
    for (auto&& x : m_boards_id) os<<(int)x<<" ";
    return os<<"]\n";
  }
#endif
};


}
waveform_by_entry::waveform_by_entry(std::string const& n):parse_base_t(n){
}

/* TODO  how to write valid ??*/ 
bool waveform_by_entry::valid(head_t const& head){
  //if (head.start_tag!=cs_start_tag) return false;
  //auto ps = head.package_size&0x1FFF;
  //if ((head.package_size&0x1FFF) != 20) return false;
  //if (head.package.crc32... /* TODO*/)
  return true;
}

bool waveform_by_entry::valid(body_t const& body){
  //if (body.start_tag!=cs_start_tag) return false;
  //if(body.adc_package_size&0x1FFF!=2060 
  //    && body.adc_package_size&0x1FFF!=14) return false;
  bool ret = 
    body.reserved2==0;
  return ret;
}

bool waveform_by_entry::valid(tail_t const& tail){
  //if (tail.start_tag!=cs_start_tag) return false;
  //if((tail.package_size&0x1FFF)!=12) return false;
  return true;
}

void waveform_by_entry::store(unit_t const& unit){
  //std::cout<<std::dec<<get_event_id(unit.heads[0])<<": "<<unit.heads.size()<<" "
  //  <<unit.bodys.size()<<" "<<unit.tails.size()<<std::endl;




  if (!m_store_ref || !m_tree_ref) return;
  m_store_ref->event_id = get_event_id(unit.heads[0]);
  m_store_ref->fec_ids.resize(unit.heads.size());
  m_store_ref->time_stamps.resize(unit.heads.size());
  for (int i=0; i<unit.heads.size(); ++i){
    uint64_t ts = get_timestamp(unit.heads[i]);
    ts &= 0xFFFFFFFFFFFF;
    m_store_ref->time_stamps[i] = ts;
  }
  std::transform(std::begin(unit.heads),std::end(unit.heads),std::begin(m_store_ref->fec_ids)
      ,[&](head_t const& v){return this->get_fec_id(v);});
  m_store_ref->global_ids.resize(unit.bodys.size());

  m_store_ref->adcs.resize(unit.bodys.size());
  std::size_t index=0;
  for (auto&& x : unit.bodys){
    m_store_ref->global_ids[index] = ((uint16_t)get_fec_id(x)<<8) + (uint16_t)get_channel_id(x);
    m_store_ref->adcs[index].resize(1024);
    std::size_t j=0;
    for (auto&& y : x.adc) m_store_ref->adcs[index][j++] = (y&0x0FFF);
    index++;
  }
  m_tree_ref->Fill();
}
void waveform_by_entry::store(){
  if (!m_store_ref || !m_tree_ref) return;
  m_store_ref->fec_ids.resize(m_unit.heads.size());
  std::transform(std::begin(m_unit.heads),std::end(m_unit.heads),std::begin(m_store_ref->fec_ids)
      ,[&](head_t const& v){return this->get_fec_id(v);});

  m_store_ref->global_ids.resize(m_unit.bodys.size());
  m_store_ref->adcs.resize(m_unit.bodys.size());
  std::size_t index=0;
  for (auto&& x : m_unit.bodys){
    m_store_ref->global_ids[index] = ((uint16_t)get_fec_id(x)<<8) + (uint16_t)get_channel_id(x);
    m_store_ref->adcs[index].resize(1024);
    std::size_t j=0;
    for (auto&& y : x.adc) m_store_ref->adcs[index][j++] = (y&0x0FFF);
    index++;
  }
  m_tree_ref->Fill();
}


bool waveform_by_entry::parse1(char*& iter, char* const& end){
#ifdef DEBUG
  struct seq_t{
    char m_tag;
    int m_info;
    int64_t m_pos;
    seq_t() = default;
    seq_t(char a, int b, int64_t c):m_tag(a),m_info(b),m_pos(c) {}
    ~seq_t() noexcept = default;
  };
  std::vector<seq_t> seq_infos;
#endif
  //typedef struct{
  //  uint32_t m_evt_id;
  //  std::set<uint8_t> m_boards_id;
  //} event_record_t;

  std::stack<::event_record_t> event_table;


  using namespace util;
  auto const& parse_head = [&](head_t& head)->bool{
    read_int(head.start_tag,iter);
    read_int(head.package_size,iter);
    if ((head.package_size&0x1FFF)!=20) return false;
    read_int(head.fec_id,iter);
    read_int<uint64_t,6>(head.time_stamp,iter);
    read_int(head.event_id,iter);
    read_int(head.hit_channel_no,iter);
    read_int(head.reserved,iter);
    read_int(head.crc32,iter);
    return true;
  };
  auto const& parse_body = [&](body_t& bd)->bool{
    read_int(bd.start_tag,iter);
    read_int(bd.adc_package_size,iter);
    if((bd.adc_package_size&0x1FFF)!=2060 
        && (bd.adc_package_size&0x1FFF)!=14)
      return false;
    read_int(bd.reserved0,iter);
    read_int(bd.channel_id,iter);
    read_int(bd.reserved1,iter);
    for (auto&& x : bd.adc) read_int(x,iter);
    read_int(bd.reserved2,iter);
    read_int(bd.crc32,iter);
    return true;
  };
  auto const& parse_tail = [&](tail_t& tail)->bool{
    read_int(tail.start_tag,iter);
    read_int(tail.package_size,iter);
    if ((tail.package_size&0x1FFF)!=12) return false;
    read_int(tail.reserved,iter);
    read_int(tail.event_size,iter);
    read_int(tail.crc32,iter);
    return true;
  };
  char* begin = iter;
  char* back = iter;

  size_t index = 0;
  int tail_count=0;

  //uint32_t current_event_idx=0;

  size_t bytes = std::distance(iter,end);
  ui::Xslider slider(0,std::distance(iter,end));
  slider.set_label("unpack");
  auto ft0 = std::async(&ui::Xslider::start,std::ref(slider));

  std::list<body_t> invalid_bodys;
  std::list<tail_t> invalid_tails;

  while(std::distance(iter,end)>0){
    slider.progress(bytes-std::distance(iter,end));
    if(iter[0]==cs_start_tag){
      head_t head; back = iter;
      if (parse_head(head) && valid(head)){
#ifdef DEBUG
        seq_infos.emplace_back(seq_t{'H',(int)get_event_id(head),get_fec_id(head)});
#endif
        auto fec_id = get_fec_id(head);
        auto evt_id = get_event_id(head);
        if (in_memory.find(evt_id)==in_memory.end()){
          in_memory[evt_id];
#ifdef DEBUG
          info_out("head get");
          //info_out(std::distance(begin,iter));
#endif

          //event_table.push(event_record_t{evt_id,fec_id});
        }
        if (event_table.empty() || event_table.top().m_evt_id != evt_id){
          event_table.push(event_record_t{evt_id,fec_id});
        }else{
          event_table.top().insert_boards(fec_id);
        }
        //std::cout<<"head get: "<<evt_id<<std::endl;
        //for (int i=-20; i<0; ++i)
        //  std::cout<<std::hex<<(int)*(iter+i)<<std::endl;
        //std::cout<<std::hex<<(int)head.start_tag<<std::dec<<std::endl;
        //std::cout<<std::hex<<(int)head.package_size<<std::dec<<std::endl;
        //std::cout<<std::hex<<(int)head.fec_id<<std::dec<<std::endl;
        //std::cout<<std::hex<<(uint64_t)head.time_stamp<<std::dec<<std::endl;
        //std::cout<<std::hex<<(int)head.event_id<<std::dec<<std::endl;
        //info_out(std::distance(iter,end));
        //info_out(std::distance(begin,iter));
        in_memory.at(evt_id).heads.push_back(head);
        continue;
      }else{
        iter=back;
      }

      body_t body; back = iter;
      if (parse_body(body) && valid(body)){
        auto fec_id = get_fec_id(body);
#ifdef DEBUG
        seq_infos.emplace_back(seq_t{'B',fec_id,std::distance(begin,iter)});
        info_out("body get");
        info_out(event_table.size());
#endif
        //info_out(event_table.size());
        if (!event_table.empty()){
#ifdef DEBUG
          //event_table.top().display();
#endif
          if (event_table.top().has_board(fec_id)){
            in_memory.at(event_table.top().event_id()).bodys.push_back(body);
          }else{
            invalid_bodys.push_back(body);
          }
          //auto& ref_curr_heads = in_memory.at(current_event_idx).heads;
          

          //in_memory.at(current_event_idx).bodys.push_back(body);
        }
        //std::cout<<"b"<<body.reserved2<<std::endl;
        continue;
      }else{
        iter = back;
      }
      
      tail_t tail; back = iter;
      if (parse_tail(tail) && valid(tail)){
        auto fec_id = get_fec_id(tail);
#ifdef DEBUG
        info_out("tail get");
        seq_infos.emplace_back(seq_t{'T',fec_id,std::distance(begin,iter)});
#endif
        if(!event_table.empty()){
          if (event_table.top().has_board(fec_id)){
            event_table.top().erase_board(fec_id);
            in_memory.at(event_table.top().event_id()).tails.emplace_back(tail);
            if (event_table.top().empty()){
              event_table.pop();
              if (!event_table.empty()){
#ifdef DEBUG
                info_out("evt-get");
                //if (invalid_bodys.size() !=0){
                //  info_out(event_table.size());
                //  info_out(invalid_bodys.size());
                //  for (auto&& x: invalid_bodys) info_out((int)get_fec_id(x));
                //  event_table.top().display();
                //  exit(0);

                //}
#endif
                for (auto iter = invalid_bodys.begin(); iter != invalid_bodys.end();){
                  auto iter_buf = std::next(iter);
                  if (event_table.top().has_board(get_fec_id(*iter))){
                    invalid_bodys.erase(iter);
                    in_memory.at(event_table.top().event_id()).bodys.emplace_back(*iter);
                  }
                  iter = iter_buf;
                }
                int loop_time=0;
                while(invalid_tails.size()) {
#ifdef DEBUG
                      for(auto&& x : invalid_tails)
                        info_out((int)get_fec_id(x));
                      info_out(event_table.top().event_id());
#endif
                      
                  loop_time++;
                  size_t table_sz = event_table.size();
                  size_t tails_sz = invalid_tails.size();
                  for (auto iter = invalid_tails.begin(); iter != invalid_tails.end(); ){
                    auto iter_buf = std::next(iter);
                    if (event_table.top().has_board(get_fec_id(*iter))){
                      invalid_tails.erase(iter);
                      event_table.top().erase_board(get_fec_id(*iter));
                      in_memory.at(event_table.top().event_id()).tails.emplace_back(*iter);
                      if (event_table.top().empty()){
                        event_table.pop();
                      }
                    }
                    iter = iter_buf;
                  }
                  if (table_sz==event_table.size() || tails_sz==invalid_tails.size())
                    break;
                    
#ifdef DEBUG
                //  if (invalid_tails.size()){
                //info_out(invalid_bodys.size());
                //info_out(invalid_tails.size());
                //info_out(event_table.top().event_id());
                ////event_table.pop();
                ////info_out(event_table.top().event_id());
                //info_out(event_table.size());
                //info_out(in_memory.size());
                //info_out(std::prev(std::end(in_memory))->first);
                //info_out(in_memory.at(499388).bodys.size());
                //info_out(in_memory.at(499278).bodys.size());
                //info_out(in_memory.at(499441).bodys.size());
                //exit(0);
                //}
#endif

                    if (loop_time>=10){

                      for(auto&& x : invalid_tails)
                        info_out((int)get_fec_id(x));
                info_out(invalid_bodys.size());
                info_out(invalid_tails.size());
                info_out(event_table.top().event_id());
                //event_table.pop();
                //info_out(event_table.top().event_id());
                info_out(event_table.size());
                info_out(in_memory.size());
                info_out(std::prev(std::end(in_memory))->first);
                info_out(in_memory.at(605644).bodys.size());
                info_out(in_memory.at(605650).bodys.size());
                info_out(std::distance(begin,iter));


                      info_out("ERROR!");
                      exit(0);
                    }
                }
#ifdef DEBUG
                //if (_debug_tag){
                //info_out(invalid_bodys.size());
                //info_out(event_table.size());
                //info_out(in_memory.size());
                //info_out(std::prev(std::end(in_memory))->first);
                //info_out(in_memory.at(26211).bodys.size());
                //info_out(in_memory.at(26212).bodys.size());
                //exit(0);
                //}
#endif

                
              }
            }

          }else{
#ifdef DEBUG
            //info_out("invalid tail get!");
            //    info_out(in_memory.at(57332).bodys.size());
            //    info_out(invalid_bodys.size());
            //    info_out(event_table.size());
            //    info_out(in_memory.size());
            //    info_out(std::prev(std::end(in_memory))->first);
            //exit(0);
#endif
            invalid_tails.emplace_back(tail);
          }



          //in_memory.at(current_event_idx).tails.emplace_back(tail);
        }
        continue;
      }else{
        iter = back;
      }

      iter++;
    }else
      iter++;
  }
  slider.stop();


  dump();

  info_out(invalid_bodys.size());
  info_out(invalid_tails.size());

#ifdef DEBUG
  info_out(invalid_bodys.size());
  info_out(invalid_tails.size());
  for(auto&& x : seq_infos)
    std::cout<<x.m_tag<<" "<<x.m_info<<" "<<x.m_pos<<std::endl;
#endif
  return true;
}
void waveform_by_entry::dump(){
  if (in_memory.size()<5000 && m_mode==1 && m_first_invoke==true) return;
  m_first_invoke = true;
  size_t evt_get = 0;
  //info_out("==============");
  //for (auto&& [x,y] : in_memory)
  //  info_out(x);
  //info_out("==============");
  for (auto iter = in_memory.begin(); iter != in_memory.end(); ++iter){
    //std::cout<<get_event_id(iter->second.heads[0])
    //  <<" "<<iter->second.heads.size()
    //  <<" "<<iter->second.bodys.size()
    //  <<" "<<iter->second.tails.size()
    //  <<std::endl;
    std::set<uint8_t> head_fecs;
    for (auto x : iter->second.heads) head_fecs.insert(get_fec_id(x));
    std::set<uint8_t> tail_fecs;
    for (auto x : iter->second.tails) tail_fecs.insert(get_fec_id(x));
    if (
        true
        //&& iter->second.heads.size()==m_fec_count
        //&& iter->second.tails.size()==m_fec_count
        //&& head_fecs==tail_fecs
        && iter->second.heads.size()==iter->second.tails.size()
        ){
      store(iter->second);
      evt_get++;
    }
  }
  std::cout<<"event parsed: "<<in_memory.size()<<" event store: "<<evt_get<<std::endl;
  in_memory.clear();

}

bool waveform_by_entry::parse(char*& iter, char* const& end){
  using namespace util;
  auto const& parse_head = [&](head_t& head)->bool{
    read_int(head.start_tag,iter);
    read_int(head.package_size,iter);
    if ((head.package_size&0x1FFF)!=20) return false;
    read_int(head.fec_id,iter);
    read_int<uint64_t,6>(head.time_stamp,iter);
    read_int(head.event_id,iter);
    read_int(head.hit_channel_no,iter);
    read_int(head.reserved,iter);
    read_int(head.crc32,iter);
    return true;
  };
  auto const& parse_body = [&](body_t& bd)->bool{
    read_int(bd.start_tag,iter);
    read_int(bd.adc_package_size,iter);
    if((bd.adc_package_size&0x1FFF)!=2060 
        && (bd.adc_package_size&0x1FFF)!=14)
      return false;
    read_int(bd.reserved0,iter);
    read_int(bd.channel_id,iter);
    read_int(bd.reserved1,iter);
    for (auto&& x : bd.adc) read_int(x,iter);
    read_int(bd.reserved2,iter);
    read_int(bd.crc32,iter);
    return true;
  };
  auto const& parse_tail = [&](tail_t& tail)->bool{
    read_int(tail.start_tag,iter);
    read_int(tail.package_size,iter);
    if ((tail.package_size&0x1FFF)!=12) return false;
    read_int(tail.reserved,iter);
    read_int(tail.event_size,iter);
    read_int(tail.crc32,iter);
    return true;
  };

  enum class e_state : uint8_t{
    k_unknow
    ,k_have_head
    ,k_in_body
    ,k_have_tail
  };

  int head_count = 0;
  char* begin = iter;

  e_state stream_state = e_state::k_unknow;
  char* back = iter;
  while(std::distance(iter,end)>0){
    if(iter[0]==cs_start_tag){
      back = iter;
      if (stream_state==e_state::k_unknow 
           || stream_state==e_state::k_have_head
           || stream_state==e_state::k_in_body
           || stream_state==e_state::k_have_tail){
        m_unit.heads.resize(m_unit.heads.size()+1);
        auto& h_ref = m_unit.heads.back();
        if (parse_head(h_ref) && valid(h_ref)){
          //info_out("head get");
          //info_out(std::distance(begin,back));
          //info_out((int)get_fec_id(h_ref));
          head_count++;
          //info_out(head_count);
          stream_state=e_state::k_have_head;
          continue;
        }else{
          iter = back;
          m_unit.heads.erase(std::prev(m_unit.heads.end()));
        }
      }
        
      back = iter;
      if (stream_state==e_state::k_have_head
          || stream_state==e_state::k_in_body
          || stream_state==e_state::k_have_tail){
        char* back = iter;
        m_unit.bodys.resize(m_unit.bodys.size()+1);
        auto& b_ref = m_unit.bodys.back();
        if (parse_body(b_ref) && valid(b_ref)){
          //info_out("body get");
          //info_out(std::distance(begin,back));
          stream_state=e_state::k_in_body;
          continue;
        }
        else{
          iter = back;
          m_unit.bodys.erase(std::prev(m_unit.bodys.end()));
        }
      }

      back = iter;
      if ((stream_state==e_state::k_have_head
          || stream_state==e_state::k_in_body
          || stream_state==e_state::k_have_tail)){
        m_unit.tails.resize(m_unit.tails.size()+1);
        auto& t_ref = m_unit.tails.back();
        if (parse_tail(t_ref) && valid(t_ref)){
          //info_out("tail get");
          //info_out(std::distance(begin,back));
          //info_out((int)get_fec_id(t_ref));
          stream_state=e_state::k_have_tail;
          head_count--;
          //info_out(head_count);
          if (head_count==0) { 
            //display(std::cout);
            store(); clear();
            stream_state=e_state::k_unknow;
          }
          continue;
        }else{
          iter = back;
          m_unit.tails.erase(std::prev(m_unit.tails.end()));
        }
      }

      iter++;
    }else {
      iter++;
    }
  }
  return true;
  
}

std::ostream& waveform_by_entry::display(std::ostream& os) const{
  os
    <<m_unit.heads.size()<<"\n";
    os<<"\t";
    for (auto&& x : m_unit.heads) os<<(int)get_fec_id(x)<<" ";
  os
    <<"\n"
    <<m_unit.bodys.size()<<"\n";
    for (auto&& x : m_unit.bodys) os<<(int)get_fec_id(x)<<" ";
  os
    <<"\n"
    <<m_unit.tails.size()<<"\n";
  return os; }

//---------------------------------------------------------------------
namespace{
std::string get_title(std::string const& name){
  auto p0 = name.find_last_of('.');
  if (p0 != std::string::npos){
    auto p1 = name.find_last_of('/');
    if (p1!=std::string::npos){
      return name.substr(p1+1,p0-p1-1);
    }
    return name.substr(0,p0);
  }
  return name;
}

inline uint16_t make_id(uint8_t fec_id, uint8_t channel_id){ return (fec_id<<8) + channel_id; }
inline std::pair<uint8_t,uint8_t> get_fec_id(uint16_t id){
  return std::make_pair(uint16_t(id>>8), uint16_t(id&0xFF)); }
inline uint32_t make_det_id(uint16_t isx, uint16_t channel_id){ return (isx<<16) + channel_id; }
inline std::pair<uint16_t,uint16_t> get_dec_id(uint32_t id){
  return std::make_pair(uint16_t(id>>16), uint16_t(id&0xFFFF)); }

}
int main(int argc, char* argv []){
  if (argc<3){
    std::cerr<<"Usage: unpack your/path/of/dat/file mode(0:baselie 1:trigger)"<<std::endl;
    return 0;
  }
  std::string dat_name = argv[1];
  int mode = std::stoi(argv[2]);
  TTree* data_tree = new TTree("CollectionTree","CollectionTree");
  entry_new entry_buffer;
  data_tree->Branch("data",std::addressof(entry_buffer));
  std::string entry_out_file = dat_name.substr(
      0,dat_name.find_last_of("."))+"_entry.root";
  TFile* fout = new TFile(entry_out_file.c_str(),"recreate");
  
  //TODO !!! read-write-queue needed! //FIXME
  std::ifstream fin(dat_name.c_str(),std::ios::binary);
  fin.seekg(0,std::ios_base::end);
  size_t fsz = fin.tellg();
  info_out(fsz);
  waveform_by_entry wf;
  //wf.fec_count(fec_count);
  wf.set_store(entry_buffer);
  wf.set_tree(data_tree);
  wf.set_mode(mode);

  std::size_t oneG = (std::size_t)1024*1024*1024;
  for(int i=0; i<std::ceil(fsz/(float)oneG); ++i){
  //for(int i=0; i<4; ++i){
    fin.seekg(i*oneG,std::ios_base::beg);
    char* data = new char[oneG];
    fin.read(data,oneG);
    auto readed = fin.gcount();
    info_out(readed);
    char* iter_beg = data;
    wf.do_parse(iter_beg,iter_beg+readed);
    delete[] data;
  }

  fin.close();
  fout->cd();
  data_tree->Write(); 
  fout->Write(); fout->Close(); 
  std::cout<<"Raw Root Store: "<<entry_out_file<<"\n";
  typedef typename util::terminal_color tc;
  using namespace util;
  std::cout<<terminal_color(
      tc::display_mode::k_underline
      ,tc::f_color::k_white
      ,tc::b_color::k_blue)
    <<"---->UNPACK DONE<----" <<util::terminal_reset() <<std::endl;

  if (mode==0){
    std::unordered_map<std::string,std::string> argv_map;
    std::unordered_map<int,float> compres_table;
    util::tt_map_t tt_map;
    util::mis_param_t mp;
    //std::string entry_out_file="";
    std::string prestal_name ="prestal.txt";
    std::string baseline_name = "baseline.root";
    prestal_name = "calc/"+get_title(dat_name)+"_prestal.txt";
    baseline_name = "calc/"+get_title(dat_name)+"_baseline.root";
    if (!std::filesystem::exists("calc"))
      std::filesystem::create_directory("calc");
    TFile* rfin = new TFile(entry_out_file.c_str());
    auto* fout = new TFile(baseline_name.c_str(),"recreate");
    auto* data_tree  = static_cast<TTree*>(rfin->Get("CollectionTree"));
    entry_new* entry_buffer_ptr = new entry_new;
    data_tree->SetBranchAddress("data",std::addressof(entry_buffer_ptr));
    std::map<int,TH1I*> baseline_map;
    for (long long i=0; i<data_tree->GetEntries(); ++i){
      data_tree->GetEntry(i);
      for (int j=0; j<entry_buffer_ptr->global_ids.size(); ++j){
        int x = entry_buffer_ptr->global_ids.at(j);
        if (baseline_map.find(x)==baseline_map.end()){
          std::stringstream sstr("");
          sstr<<"baseline-"<<(int)(x>>8)<<"_"<<(int)(x&0xFF);
          auto const& ref = entry_buffer_ptr->adcs[j];
          uint16_t min = *util::min_element(std::begin(ref),std::end(ref));
          uint16_t max = *util::max_element(std::begin(ref),std::end(ref));
          min -= 10; max += 10;
          baseline_map[x] = new TH1I(sstr.str().c_str(),sstr.str().c_str()
              ,max-min,min,max);
        }
        for (auto&& adc : entry_buffer_ptr->adcs[j]){ baseline_map[x]->Fill(adc); }
      }
    }
    fout->cd();
    std::unordered_map<int,std::pair<float,float>> mean_and_rms;


    info_out(prestal_name);
    std::ofstream prestal_txt(prestal_name.c_str());
    prestal_txt<<"#"<<"Prestal File V0.0.1\t"<<util::time_to_str()<<"\n";
    prestal_txt<<"#"<<"Generate by "<<entry_out_file<<"\n";
    prestal_txt<<"#"<<"fecid channelid mean sigma chi2/ndf\n";
    for (auto iter = baseline_map.begin(); iter != baseline_map.end(); ++iter) {
      TF1 f_gaus("f_gaus","gaus"
          ,iter->second->GetMean()-3*iter->second->GetRMS()
          ,iter->second->GetMean()+3*iter->second->GetRMS());
      iter->second->Fit(&f_gaus,"RQ");
      mean_and_rms[iter->first].first=iter->second->GetMean();
      mean_and_rms[iter->first].second=iter->second->GetRMS();

      prestal_txt<<(int)::get_fec_id(iter->first).first<<
        " "<<(int)::get_fec_id(iter->first).second
        <<" "<<f_gaus.GetParameter(1)<<" "
        <<f_gaus.GetParameter(2)
        <<" "<<f_gaus.GetChisquare()/f_gaus.GetNDF()
        <<"\n";
    }
    prestal_txt.close();

    TH1I mean_dis("mean","mean",3000,0,3000);
    TH1F rms_dis("rms","rms",3000,0,3000);
    for (auto [x,y] : mean_and_rms){
      mean_dis.SetBinContent(x,y.first);
      rms_dis.SetBinContent(x,y.second); }
    mean_dis.Write(); rms_dis.Write();
    rfin->Close(); fout->Write(); fout->Close();
    delete entry_buffer_ptr;

    typedef typename util::terminal_color tc;
    std::cout<<util::terminal_color(
        tc::display_mode::k_underline
        ,tc::f_color::k_white
        ,tc::b_color::k_blue)
      <<"---->PRESTAL<----" <<util::terminal_reset() <<std::endl;
    std::cout<<"Baseline Root Store: "<<baseline_name<<"\n";
  }

  return 0;
}
