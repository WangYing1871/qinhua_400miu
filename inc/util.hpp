#ifndef util_HPP
#define util_HPP 1

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <map>
#include <vector>
#include <queue>

#define info_out(X) std::cout<<"==> "<<__LINE__<<" "<<#X<<" |"<<(X)<<"|\n"
//#ifndef __CINT__
//#include <boost/type_index.hpp>
//#include <boost/timer/timer.hpp>
//#define type_name(x) std::cout<<std::dec<<boost::typeindex::type_id<decltype(x)>().pretty_name() \
//  <<" sizeof:"<<sizeof(decltype(x)) \
//  <<" alignof:"<<alignof(decltype(x))<<std::endl;
//#endif
namespace util{
std::string time_to_str();
 
namespace math{
double gaus(double* x,double* p);
}

//template <class _iter_t>
//bool is_peak(_iter_t first, _iter_t last){
//  int distance = std::distance(first,last);
//  if (std::abs(distance)<3) return false;
//  size_t side_sz = std::floor(distance/3.);
//  size_t middle_sz = distance-2*side_sz;
//  long long side0 = std::accumulate(first,std::next(first,side_sz),0);
//  long long middle = std::accumulate(std::next(first,side_sz),std::prev(last,side_sz),0);
//  long long side1 = std::accumulate(std::prev(last,side_sz),last,0);
//  return (middle>=side1 && middle>=side0 && side0>=0 && side1>=0)
//  || (middle<=side1 && middle<=side0 && side0<=0 && side1<=0);
//}

//template <std::size_t _vs,class _iter>
//std::vector<_iter> first_max_element(_iter begin, _iter end){
//  std::vector<_iter> ret;
//  if (std::distance(first,last)<_vs) return ret;
//  ret.resize(_vs);
//
//
//
//
//}

template <class _iter>
_iter max_element(_iter first, _iter last){
  if (std::distance(first,last)==0) return first;
  _iter ret=first;
  auto me = *first;
  for (;first != last; ++first) if(*first>me) {me = *first; ret=first;}
  return ret;
}
template <class _iter>
_iter min_element(_iter first, _iter last){
  if (std::distance(first,last)==0) return first;
  _iter ret=first;
  auto me = *first;
  for (;first != last; ++first) if(*first<me) {me = *first; ret=first;}
  return ret;
}

inline void trim_space(std::string& s){
  for (auto iter = s.begin(); iter != s.end(); ++iter)
    if (!std::isspace(*iter))
      for(auto riter = s.rbegin(); riter != s.rend(); ++riter)
        if (!std::isspace(*riter)){
          s = s.substr(
              std::distance(s.begin(),iter)
              ,std::distance(riter,s.rend())-std::distance(s.begin(),iter)
              ); return; }
  s.clear();
}

namespace _meta{
template <int _idx_v, class... _args>
auto get_v(_args&&... params){
  auto tuple = std::make_tuple<_args...>(std::forward<_args>(params)...);
  if constexpr (_idx_v>=sizeof...(_args)) return std::nullopt;
  else return std::make_optional(std::get<_idx_v>(tuple));
}

template <class _tp ,int _idx_v,class... _args>
_tp get_param(_args&&... params){
  _tp ret{};
  auto tuple = std::make_tuple<_args...>(std::forward<_args>(params)...);
  if constexpr (sizeof...(_args)>=(_idx_v+1)){
    auto p = std::get<_idx_v>(tuple);
    if constexpr (std::is_constructible
        <_tp,typename std::remove_reference<decltype(p)>::type>::value){
      ret = p;
    }
  }
  return ret;
}
}
}

namespace util{
std::unordered_map<std::string,std::string> read_argv(std::string const& pg);

template <class _tp
  ,int _bytes_v = sizeof(_tp)
  ,bool _is_big_end_v = false
  ,class = typename std::enable_if<std::is_integral<_tp>::value>::type>
inline void read_int(_tp& data, char*& iter){
  static_assert(_bytes_v<=sizeof(_tp));
  char rt[_bytes_v];
  if constexpr (_is_big_end_v) for (int i=0; i<_bytes_v; ++i) rt[i] = *iter++;
  else for (int i=_bytes_v-1; i>=0; --i) rt[i] = *iter++;
  data = *reinterpret_cast<_tp*>(rt); }
inline void read_int(uint8_t& data, char*& iter){ data = *iter++; }

enum class e_command_address : uint16_t{
  k_FEID = 0x1806
  ,k_ChannelID = 0x1808
  ,k_Threshold_Tag= 0x180B
  ,k_Threshold_Value = 0x180A
  ,k_Config = 0x180C
  ,k_BaseLine = 0x1112 /*0x1112-0x1116*/
  ,k_TThreshold = 0x1814
  ,k_Confirm = 0x0001
  ,k_Reset = 0x0000
  ,k_IIRFilter = 0x1100
  ,k_TrigDelayCycle = 0x1804
  ,k_HitWidthCycle = 0x180e
  ,k_TrigRiseStep = 0x1810
  ,k_NHitChannel = 0x1802
};
inline void head_write(char* positon, e_command_address cmd, char base=0x00){
  uint16_t value = uint16_t(cmd);
  positon[0] = base+0x00 + (value&0x0F);
  positon[1] = base+0x10 + ((value>>4)&0x0F);
  positon[2] = base+0x20 + ((value>>8)&0x0F);
  positon[3] = base+0x30 + ((value>>12)&0x0F);
};
inline void head_write(char* positon, uint16_t cmd, char base = 0x00){
  uint16_t value = uint16_t(cmd);
  positon[0] = base+0x00 + (value&0x0F);
  positon[1] = base+0x10 + ((value>>4)&0x0F);
  positon[2] = base+0x20 + ((value>>8)&0x0F);
  positon[3] = base+0x30 + ((value>>12)&0x0F);
};
template <class... _args>
void generate_frame(e_command_address command, char* positon
    , _args&&... params){
  if (!positon) return;
  switch(command){
    case (e_command_address::k_FEID):{
      head_write(positon,command);
      uint8_t fec_id = _meta::get_param<uint8_t,0>(std::forward<_args>(params)...);
      positon[4] = 0x40+fec_id%16;
      positon[5] = 0x50+fec_id/16;
      positon[6] = 0x60;
      positon[7] = 0x70;
      positon[8] = 0x83;
      break;
    }
    case (e_command_address::k_ChannelID):{
      head_write(positon,command);
      uint8_t cid = _meta::get_param<uint8_t,0>(std::forward<_args>(params)...);
      positon[4] = 0x40+cid%16;
      positon[5] = 0x50+cid/16;
      positon[6] = 0x60;
      positon[7] = 0x70;
      positon[8] = 0x83;
      break;
    }
    case (e_command_address::k_Threshold_Tag):{
      head_write(positon,command);
      positon[4] = 0x40;
      positon[5] = 0x50;
      positon[6] = 0x60;
      positon[7] = 0x70;
      positon[8] = 0x83;
      break;
    }
    case (e_command_address::k_Threshold_Value):{
      head_write(positon,command);
      uint32_t mean=_meta::get_param<uint32_t,0>(std::forward<_args>(params)...);
      uint32_t sigma=_meta::get_param<uint32_t,1>(std::forward<_args>(params)...);
      float compres=_meta::get_param<float,2>(std::forward<_args>(params)...);
      uint64_t value = mean+compres*sigma;
      uint16_t value_u16 = value;
      if (value>0x0FFF){
        std::cerr<<"Warning!! threshold set overflow to 12 bits max! Set to '0x0FFF'"<<"\n";
        value_u16 = 0x0FFF; }

      //std::cout<<mean<<" "<<sigma<<" "<<compres<<" "<<std::hex<<value_u16<<std::dec<<std::endl;

      positon[4] = 0x40+(value_u16&0x0F);
      positon[5] = 0x50+((value_u16>>4)&0x0F);
      positon[6] = 0x60+((value_u16>>8)&0x0F);
      positon[7] = 0x70+((value_u16>>12)&0x0F);
      positon[8] = 0x83;
      break;
    }
    case (e_command_address::k_Config):{
      head_write(positon,command);
      positon[4] = 0x41;
      positon[5] = 0x50;
      positon[6] = 0x60;
      positon[7] = 0x70;
      positon[8] = 0x83;
      break;
    }
    case (e_command_address::k_BaseLine):{
      //std::vector<uint16_t> bl = _meta::get_param<std::vector<uint16_t>,0>(std::forward<_args>(params)...);
      //std::vector<uint16_t> bl{};
      //if constexpr(std::is_same<decltype(opt),std::nullopt_t>::value);
      //else bl = opt.value();
      

      break;
    }
    default:
      break;
  }
}


void generate_frame_trigger(char*& positon,uint16_t b, uint16_t f, uint16_t t);
std::size_t generate_frame_baseline(char* positon, std::vector<uint16_t> const& bl);
std::size_t generate_frame_irrfilter(char* positon, std::vector<uint16_t> const& bl);
std::size_t generate_frame_trig_delay_cycle(char* positon, uint16_t v);
std::size_t generate_frame_hitwidth_cycle(char* positon, uint8_t v);
std::size_t generate_frame_trig_rise_step(char* positon, uint8_t v);
std::size_t generate_frame_nhit_channel(char* positon, uint8_t v);

}

namespace util{
struct terminal_color{
  enum class display_mode : uint8_t{
     k_default_ = 0
    ,k_highlight
    ,k_underline
    ,k_flash
    ,k_anti
    ,k_not_visible
    ,k_anti_bold
    ,k_anti_underline
    ,k_anti_flash
    ,k_anti_anti
    ,k_visible
  };

  enum class f_color : uint8_t{
    k_black = 30
    ,k_red
    ,k_green
    ,k_yellow
    ,k_blue
    ,k_pink
    ,k_cyan
    ,k_white
  };

  enum class b_color : uint8_t{
    k_black = 40
    ,k_red
    ,k_green
    ,k_yellow
    ,k_blue
    ,k_pink
    ,k_cyan
    ,k_white
  };

  display_mode mode;
  f_color fc;
  b_color bc;

  inline terminal_color(uint8_t a, uint8_t b, uint8_t c){
    mode = display_mode(a); fc = f_color(b); bc = b_color(c); }
  inline terminal_color(display_mode a=display_mode::k_underline
      , f_color b=f_color::k_white, b_color c=b_color::k_green){
    mode = a, fc = b, bc = c; }
  std::string gen() const;
};
inline std::ostream& operator<<(std::ostream& os, terminal_color const& b){ return os<<b.gen().c_str(); }
inline std::ostream& operator<<(std::ostream& os, terminal_color&& b){ return os<<b.gen().c_str(); }
struct terminal_reset{ };
inline std::ostream& operator<<(std::ostream& os, terminal_reset&& b){ return os<<"\033[0m"; }
}

#include "data_strcut_cint.h"
namespace util{
void cluster_comp(cluster& v);

struct lcs_t{
  int layer_id=0;
  int channel_id=0;
  float compres=3;
};
struct lc_adc_t{
  uint16_t layer_id=0;
  uint16_t channel_id=0;
  uint16_t adc_t=0;
};

struct mis_param_t{
 uint8_t s_rise_step = 0x0;
 uint16_t s_delay_time = 0x258;
 uint8_t s_wait_cycle = 10;
 uint8_t s_nhit_strips = 0x02;
};

lcs_t parse_lcs(std::string str);
lc_adc_t parse_lc_adc(std::string str);

typedef std::unordered_map<int,std::pair<float,float>> mean_and_rms_t;
typedef std::unordered_map<uint16_t,uint16_t> tt_map_t;
void generate_configs(
  std::string const& path
  ,std::unordered_map<int,float> const& compres_table
  ,mean_and_rms_t const& maps
  ,float compres
  );
void generate_mis(std::string const& path,mis_param_t para);
void generate_tt(std::string const& path,tt_map_t const& map);

}
#include <atomic>
#include <mutex>
namespace util{
namespace _ui{
struct Xslider{
  typedef Xslider self_t;

public:
  Xslider() = default;
  Xslider(std::size_t begin, std::size_t end):begin(begin),end(end) {}
  ~Xslider() noexcept { stop(); }
  inline void progress(size_t p){ current.store(p); }
  inline self_t& set_label(std::string const& v) {label=v; return *this;}
  inline self_t& set_char(char a, char b) {
    std::lock_guard<std::mutex> lock(m_mutex);
    chars = std::make_pair(a,b);
    return *this; }
  void start();
  void stop();

private:
  std::mutex m_mutex;
  std::atomic<bool> stop_tag=false;
  std::pair<char,char> chars{'=','-'};
  size_t updatedel=100;
  std::size_t begin=0;
  std::size_t end=100;
  //std::size_t length=100; //FIXME
  std::atomic<size_t> current=0;
  std::string label="progress";
};

}
}

namespace ui = util::_ui;

#include <boost/program_options.hpp>
namespace util{
struct getopt{
  std::string project_name = "JHS@JV\nTPC\nOpenSource Analysis Codes";
  boost::program_options::options_description desc;
  boost::program_options::variables_map vm;
  int run(int argc, char* argv[]);

  template <class _tp>
  bool count(std::string const& key,_tp& v){
    if (vm.count(key.c_str())){
      v = vm[key].as<_tp>();
      return true;
    }else return false;
  }
};
}

#endif
