#ifndef unpack_H
#define unpack_H 1 
#include <list>
#include <cstdint>
struct parse_base_t{
  std::string m_name;

public:
  parse_base_t() = default;
  parse_base_t(std::string const& v):m_name(v) {}
  ~parse_base_t() noexcept = default;

  virtual bool parse(char*&, char* const&) = 0;
  virtual void clear() {}
  virtual std::ostream& display(std::ostream& os = std::cout) const {return os;}
};
#include "data_strcut_cint.h"
#include "TTree.h"
struct waveform_pack_head{
  uint8_t start_tag;
  //bit-field defined method [deprecated]
  //uint64_t package_size:16;
  //uint64_t time_stamp:48;
  uint16_t package_size;
  uint8_t fec_id;
  uint64_t time_stamp;
  uint32_t event_id;
  uint8_t hit_channel_no;
  uint8_t reserved;
  uint32_t crc32;
  //waveform_pack_head& operator=(waveform_pack_head&) = default;

};
struct waveform_pack_body{
  uint8_t start_tag;
  uint16_t adc_package_size;
  uint8_t reserved0;
  uint8_t channel_id;
  uint8_t reserved1;
  uint16_t adc[1024];
  uint16_t reserved2;
  uint32_t crc32;

};

struct waveform_pack_tail{
  uint8_t start_tag;
  uint16_t package_size;
  uint8_t reserved;
  uint32_t event_size;
  uint32_t crc32;
};


struct waveform_by_entry : public parse_base_t{
  typedef waveform_by_entry self_t;
  typedef parse_base_t base_t;
  typedef waveform_pack_head head_t;
  typedef waveform_pack_body body_t;
  typedef waveform_pack_tail tail_t;

  waveform_by_entry() = default;
  waveform_by_entry(std::string const&);
  ~waveform_by_entry() noexcept = default;


  virtual bool parse(char*& iter, char* const& end) override;
  bool parse1(char*& iter, char* const& end);
  virtual void clear() override {m_unit.bodys.clear(); m_unit.heads.clear(); m_unit.tails.clear();}
  virtual std::ostream& display(std::ostream& os = std::cout) const override;

public:
  inline void set_store(entry_new& ref) {m_store_ref=std::addressof(ref);}
  inline void set_tree(TTree* ref) {m_tree_ref=ref;}


private:
  bool valid(head_t const&);
  bool valid(body_t const&);
  bool valid(tail_t const&);

private:
  entry_new* m_store_ref=nullptr;
  TTree* m_tree_ref=nullptr;
  struct unit_t{
    std::vector<head_t> heads;
    std::list<body_t> bodys;
    std::vector<tail_t> tails;
  };
  unit_t m_unit;
  std::map<uint32_t,unit_t> in_memory;
  void store();
  void store(unit_t const&);
  constexpr static uint8_t const cs_start_tag = 0x5a;

  std::size_t board_no = 1;
public:
  inline uint32_t get_event_id(head_t const& head) const{ return head.event_id;}
  inline uint8_t get_fec_id(head_t const& v) const{ return v.fec_id&0x3F;}
  inline uint8_t get_fec_id(body_t const& v) const{ return v.reserved0&0x3F;}
  inline uint8_t get_fec_id(tail_t const& v) const{ return v.reserved&0x3F;}

  inline uint64_t get_timestamp(head_t const& v) const { return v.time_stamp;}
  inline uint8_t get_channel_id(body_t const& v) const{ return v.channel_id & 0x7F; }

  bool do_parse(char*& iter, char* const& end){
    //parse(iter,end);
    return parse1(iter,end);
    //return parse(iter,end);
  }

  void dump();
  int m_mode = 1;
  bool m_first_invoke = false;

public:
  void set_mode(int v) {m_mode=v;}


  struct event_id_recordor{
    //TODO

  };

  uint8_t m_fec_count=1;

  void fec_count(uint8_t v) {m_fec_count=v;}
public:
  static constexpr std::index_sequence<8,8> const glb_id_idx={};
  static constexpr float const s_ts_unit = 8.33; /* ns */
};
#endif
