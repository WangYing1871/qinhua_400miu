#include "util.hpp"
#include <chrono>
namespace{
auto const& format00 = [](std::size_t v)->std::string{
  return v<10 ? std::string("0")+std::to_string(v) : std::to_string(v); };
}

namespace util{

namespace math{
double gaus(double* x,double* p){
  return p[0]*std::exp(-(x[0]-p[1])*(x[0]-p[1])/2/p[2]/p[2]); }
}

std::unordered_map<std::string,std::string> read_argv(std::string const& pg){
  std::ifstream fin("argv.ini");
  if (!fin){
    throw std::invalid_argument("No InI file");
    exit(0);
  }
  std::string sbuf;
  std::unordered_map<std::string,std::string> ret;
  std::stringstream pgname(""); pgname<<"["<<pg<<"]";
  while(!fin.eof()){
    std::getline(fin,sbuf);
    if (sbuf==pgname.str()){
      //std::cout<<"WY\n";
      while(!fin.eof()){
        auto back = fin.tellg();
        std::getline(fin,sbuf);
        if (auto iter = sbuf.find_first_of('#'); iter != std::string::npos){
          bool tag=true;
          for (std::size_t i=0; i<iter; ++i) if (!std::isspace(sbuf[i])){ tag=false; break;}
          if (tag) continue; }
        if (sbuf.find_first_of('[') != std::string::npos) { fin.seekg(back); break;}
        if (!sbuf.empty() && sbuf.find_first_of('=') != std::string::npos){
          std::string key = sbuf.substr(0,sbuf.find_first_of('='));
          std::string value =  sbuf.substr(sbuf.find_first_of('=')+1);
          trim_space(key); trim_space(value);
          if (!key.empty() && !value.empty()) ret[key]=value;
        }
      }
    }
  }
  fin.close();
  return ret;
}
void generate_frame_trigger(char*& positon,uint16_t b, uint16_t f, uint16_t t){
  //head_write(positon,e_command_address::k_TThreshold);
  //head_write(positon+4,0x0001);
  t &= 0x0FFF;
  head_write(positon,e_command_address::k_FEID);
  positon[4] = 0x40 + (b&0x000F);
  positon[5] = 0x50 + ((b>>4)&0x000F);
  positon[6] = 0x60 + ((b>>8)&0x000F);
  positon[7] = 0x70 + ((b>>12)&0x000F);
  positon[8] = 0x83;
  head_write(positon+9,e_command_address::k_ChannelID);
  positon[13] = 0x40 + (f&0x000F);
  positon[14] = 0x50 + ((f>>4)&0x000F);
  positon[15] = 0x60 + ((f>>8)&0x000F);
  positon[16] = 0x70 + ((f>>12)&0x000F);
  positon[17] = 0x83;
  head_write(positon+18,e_command_address::k_Threshold_Value);
  positon[22] = 0x40 + (t&0x000F);
  positon[23] = 0x50 + ((t>>4)&0x000F);
  positon[24] = 0x60 + ((t>>8)&0x000F);
  positon[25] = 0x70 + ((t>>12)&0x000F);
  positon[26] = 0x83;
  head_write(positon+27,e_command_address::k_Config);
  head_write(positon+31,e_command_address::k_Confirm,0x40);
  positon[35] = 0x83;
  positon+=36;
}

std::size_t generate_frame_baseline(char* positon, std::vector<uint16_t> const& bl){
  char* begin = positon;
  assert(bl.size()>=3);
  for (size_t i=0; i<bl.size(); ++i){
    head_write(positon,(uint16_t)e_command_address::k_BaseLine+2*i);
    head_write(positon+4,bl[i],0x40);
    positon[8] = 0x83; positon += 9;
  }
  return std::distance(begin,positon);
}
std::size_t generate_frame_irrfilter(char* positon, std::vector<uint16_t> const& bl){
  char* begin = positon;
  assert(bl.size()>=8);
  for (size_t i=0; i<bl.size(); ++i){
    head_write(positon,(uint16_t)e_command_address::k_IIRFilter+2*i);
    head_write(positon+4,bl[i],0x40);
    positon[8] = 0x83; positon += 9;
  }
  return std::distance(begin,positon);
}
std::size_t generate_frame_trig_delay_cycle(char* positon, uint16_t v){
  head_write(positon,e_command_address::k_TrigDelayCycle);
  head_write(positon+4,v,0x40);
  positon[8] = 0x83;
  return 9;
}
std::size_t generate_frame_hitwidth_cycle(char* positon, uint8_t v){
  head_write(positon,e_command_address::k_HitWidthCycle);
  head_write(positon+4,v,0x40);
  positon[8] = 0x83;
  return 9;
}
std::size_t generate_frame_trig_rise_step(char* positon, uint8_t v){
  head_write(positon,e_command_address::k_TrigRiseStep);
  head_write(positon+4,v,0x40);
  positon[8] = 0x83;
  return 9;
}
std::size_t generate_frame_nhit_channel(char* positon, uint8_t v){
  head_write(positon,e_command_address::k_NHitChannel);
  head_write(positon+4,v,0x40);
  positon[8] = 0x83;
  return 9; }

std::string terminal_color::gen() const{
  return  
    std::string{"\033"}
    + std::string{"["}
    + std::to_string(std::underlying_type<display_mode>::type(mode))
    + std::string{";"}
    + std::to_string(std::underlying_type<f_color>::type(fc))
    + std::string{";"}
    + std::to_string(std::underlying_type<b_color>::type(bc))
    + std::string{"m"}
    ;
}
std::string time_to_str(){
  using namespace std::chrono;
  system_clock::time_point tp_now = system_clock::now();
  system_clock::duration dsp = tp_now.time_since_epoch();
  time_t msp = duration_cast<microseconds>(dsp).count();
  time_t sse = msp/1000000;
  std::tm ct = *std::localtime(&sse);
  std::stringstream sstr;
  sstr<<1900 + ct.tm_year << "-"<<1+ ct.tm_mon << "-" << ct.tm_mday
    <<"_"<<format00(ct.tm_hour)<<"-"<<format00(ct.tm_min)<<"-"<<format00(ct.tm_sec)<<"_"<< msp%1000000;
  return sstr.str(); }


}

#include <numeric>
namespace util{
void cluster_comp(cluster& v){
  v.strips = v.range.second-v.range.first;
  v.sum_amp = std::accumulate(v.amp.begin(),v.amp.end(),0); }
}

namespace util{
lcs_t parse_lcs(std::string str){
  lcs_t rt;
  int index=0;
  while(str.find_first_of('>')!=std::string::npos){
    std::string v_str = str.substr(str.find_first_of('<')+1,str.find_first_of('>')-1);
    if (index==0) rt.layer_id=std::stoi(v_str.c_str());
    else if (index==1) rt.channel_id= std::stoi(v_str.c_str());
    else if (index==2) rt.compres=std::stof(v_str.c_str());
    str = str.substr(str.find_first_of('>')+1);
    index++; }
  return rt; }

lc_adc_t parse_lc_adc(std::string str){
  lc_adc_t rt;
  int index=0;
  while(str.find_first_of('>')!=std::string::npos){
    std::string v_str = str.substr(str.find_first_of('<')+1,str.find_first_of('>')-1);
    if (index==0) rt.layer_id=std::stoul(v_str.c_str());
    else if (index==1) rt.channel_id= std::stoul(v_str.c_str());
    else if (index==2){ std::stringstream sstr(""); sstr<<v_str; sstr>>std::hex>>rt.adc_t; }
    str = str.substr(str.find_first_of('>')+1);
    index++; }
  return rt; }

void generate_configs(
  std::string const& path
  ,std::unordered_map<int,float> const& compres_table
  ,mean_and_rms_t const& maps
  ,float compres
  ){
  float compres_back = compres;
  using namespace util;
  char* frame_buf = new char[9];
  std::vector<char> frame(0);
  for (auto [x,y] : maps){
    uint16_t gid = x;
    if (compres_table.find(gid) != compres_table.end()) compres=compres_table.at(gid);
    else compres=compres_back;
    uint8_t fec_id= uint8_t(x>>8);
    uint8_t channel_id = uint8_t(x&0xFF);
    std::stringstream sstr("");
    sstr<<"L"<<(int)fec_id<<"CompressThreshold_"<<(int)channel_id<<".dat";
    if (!std::filesystem::exists(path.c_str()))
      std::filesystem::create_directory(path.c_str());
    std::string fname=path+"/"+sstr.str();
    std::ofstream fout(fname.c_str(),std::ios::out|std::ios::binary);

    head_write(frame_buf,e_command_address::k_TThreshold);
    head_write(frame_buf+4,e_command_address::k_Reset,0x40);
    frame_buf[8] = 0x83;
    for (int i=0; i<9; ++i) frame.push_back(frame_buf[i]);

    generate_frame(e_command_address::k_FEID,frame_buf,fec_id);
    for (int i=0; i<9; ++i) frame.push_back(frame_buf[i]);
    generate_frame(e_command_address::k_ChannelID,frame_buf,channel_id);
    for (int i=0; i<9; ++i) frame.push_back(frame_buf[i]);
    generate_frame(e_command_address::k_Threshold_Tag,frame_buf);
    for (int i=0; i<9; ++i) frame.push_back(frame_buf[i]);
    for (int i=0; i<9; ++i) frame_buf[i]=0;
    generate_frame(e_command_address::k_Threshold_Value,frame_buf
        ,(uint16_t)std::floor(y.first)
        ,(uint16_t)std::floor(y.second)
        ,compres);
    for (int i=0; i<9; ++i) frame.push_back(frame_buf[i]);
    for (int i=0; i<9; ++i) frame_buf[i]=0;
    generate_frame(e_command_address::k_Config,frame_buf);
    for (int i=0; i<9; ++i) frame.push_back(frame_buf[i]);
    fout.write(frame.data(),frame.size());
    frame.clear();
    fout.close();
  }
  delete[] frame_buf;
}

void generate_tt(std::string const& path
    ,tt_map_t const& map){
   using namespace util;
  //std::string fname = "HitThreshold.dat";
  //fname = path+"/"+fname;
  auto* data = new char[map.size()*36+9];
  char* begin = data;
  auto const& helper = [&](uint16_t fec_id, uint16_t channel_id
      ,uint16_t y){
    char data[36+9];
    head_write(data,e_command_address::k_TThreshold);
    head_write(data+4,e_command_address::k_Confirm,0x40);
    data[8] = 0x83;
    char* begin = &data[9];
    generate_frame_trigger(begin,fec_id,channel_id,y);
    std::stringstream sstr("");
    sstr<<path+"/"+"HitThreshold"<<fec_id<<"_"<<channel_id<<".dat";
    std::ofstream fout(sstr.str().c_str(),std::ios::out|std::ios::binary);
    fout.write(data,45);
    fout.close();
  };
  data[8] = 0x83;
  data+=9;
  for (auto [x,y] : map){
    uint8_t fec_id= uint8_t(x>>8);
    uint8_t channel_id = uint8_t(x&0xFF);
    helper(fec_id,channel_id,y);
  }
  delete[] begin;
}

void generate_mis(std::string const& path
    ,mis_param_t mp){
  using namespace util;
  char* data = new char[100];
  std::string fname = "Baseline.dat";
  {
  fname = path+"/"+fname;
  std::ofstream fout(fname.c_str(),std::ios::binary|std::ios::out);
  fout.write(data,generate_frame_baseline(data,{0,0,0}));
  fout.close();
  }

  {
  fname = "HitWidthCycle.dat";
  fname = path+"/"+fname;
  std::ofstream fout(fname.c_str(),std::ios::binary|std::ios::out);
  fout.write(data,generate_frame_hitwidth_cycle(data,mp.s_wait_cycle));
  fout.close();
  }

  {
  fname = "IIRFilter.dat";
  fname = path+"/"+fname;
  std::ofstream fout(fname.c_str(),std::ios::binary|std::ios::out);
  fout.write(data,generate_frame_irrfilter(data,{0,0,0,0,0,0,0,0}));
  fout.close();
  }

  {
  fname = "TrigRiseStep.dat";
  fname = path+"/"+fname;
  std::ofstream fout(fname.c_str(),std::ios::binary|std::ios::out);
  fout.write(data,generate_frame_trig_rise_step(data,mp.s_rise_step));
  fout.close();
  }

  {
  fname = "TrigDelayCycle.dat";
  fname = path+"/"+fname;
  std::ofstream fout(fname.c_str(),std::ios::binary|std::ios::out);
  //std::cout<<__LINE__<<" "<<delay_time<<std::endl;
  fout.write(data,generate_frame_trig_delay_cycle(data,mp.s_delay_time));
  fout.close();
  }

  {
  fname = "NHitChannel.dat";
  fname = path+"/"+fname;
  std::ofstream fout(fname.c_str(),std::ios::binary|std::ios::out);
  fout.write(data,generate_frame_nhit_channel(data,mp.s_nhit_strips));
  fout.close();
  }
  delete[] data; }

}

#include <thread>
namespace util::_ui{
void Xslider::start(){
  while(!stop_tag.load()){
    float schedule = (current.load()-begin)/(float)(end-begin);
    schedule*=100.; schedule=std::round(schedule);
    schedule = schedule>=100. ? 100. : schedule;
    std::cout<<"\r"<<label<<"\t|";
    for (size_t i=0; i<schedule; ++i) std::cout.put(chars.first);
    for (size_t i=schedule; i<100; ++i) std::cout.put(chars.second);
    std::cout<<"|\t"<<std::fixed<<std::setprecision(1) <<schedule<<" %"<<std::flush;
    std::this_thread::sleep_for(std::chrono::milliseconds(updatedel));
  }
  std::cout<<"\n";
}
void Xslider::stop(){ stop_tag.store(true); }

}

namespace util{
int getopt::run(int argc, char* argv[]){
  namespace po = boost::program_options;
  //desc = boost::program_options::options_description{project_name.c_str()};
  desc.add_options()
    ("help", "JHS@JV analysis")
    //("config,c",po::value<std::string>()
    //  ->default_value("argv.ini"),"Input param by Ini file")
    ("config,c","Input param by Ini file")
    ("file,f",po::value<std::string>(),"raw dat file")
    ("rfile,r",po::value<std::string>()->default_value("<unpack>"),"unpacked root file")
    ("board,b",po::value<int>()->default_value(6),"FEE board number")
    ("map,m",po::value<std::string>(),"fec to det map")
    ("prestal,p",po::value<std::string>(),"prestal file")
    ("mode,M",po::value<int>()->default_value(1),"Handle mode\n\t[0]: BaseLine\n\t[1]: Signal")
    ("rawfile,R",po::value<std::string>()->default_value("<unpack>")
             ,"Raw root file to generate prestal")
    ("prestal_file,P",po::value<std::string>()->default_value("prestal.txt")
             ,"prestal file to sotre path")
    ("view,VV",po::value<std::string>()->default_value("filter.root")
             ,"filter root file store")
    ;
  try{
    po::store(po::parse_command_line(argc,argv,desc),vm);
    po::notify(vm);
  }catch(std::exception& e){
    std::cerr<<"[ERROR]: "<<e.what()<<"\n"
      <<"HELP: \n"<<desc
      <<"\n";
    return -1;

  }catch(...){
    std::cerr<< "Exception of unknown type!\n";
    return -1;
  }
  if (vm.count("help")){
    std::cout<<desc<<"\n";
    return 0;
  }
  return -1;
}
}

