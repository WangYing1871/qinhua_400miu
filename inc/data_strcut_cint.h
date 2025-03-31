#ifndef data_strcut_cint_H
#define data_strcut_cint_H 1 
#include "TObject.h"
#include <array>
#include <list>
#include <vector>
#include <cstdint>
struct entry_new : public TObject{
  /* raw_id */ uint32_t event_id;
  std::vector<uint8_t> fec_ids;
  std::vector<uint8_t> hit_channel_nos;
  std::vector<uint64_t> time_stamps;
  std::vector<uint16_t> global_ids;
  std::vector<std::vector<uint16_t>> adcs;

public:
  ClassDef(entry_new,1)
};

struct cluster : public TObject{
  std::pair<int,int> range;
  std::vector<uint16_t> amp;
  std::vector<int> amp_time;
  std::size_t strips;
  uint16_t sum_amp;

public:
  //cluster() = default;
  //~cluster() noexcept = default;
public:

  ClassDef(cluster,1)
};

struct entry_tm : public TObject{
  std::vector<uint32_t> det_ids;
  std::vector<uint64_t> start_times;
  std::vector<std::vector<uint16_t>> adcs;

public:
  ClassDef(entry_tm,1)
};



#endif
