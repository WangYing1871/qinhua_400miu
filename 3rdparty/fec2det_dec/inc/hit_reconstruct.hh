#pragma once
#include <iostream>
#include <cstring>
#include <vector>
#include <ctime>   
#include <iomanip>
// root lib
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TGraph.h>
// Self-define lib
#include "hit_dec.hh"
#include "constant.hh"
#include "PositionData.hh" // 

class hit_reconstruct
{
public:
	hit_reconstruct(string save_filename
      , string l_filename
      , bool is_write_txt = false
      , string enc_filename = "./encoding_list_64to512.csv"
      , string a_idx_filename = "./group_lists.csv"
      , string enc_filename2 = ""
      , bool is_allow_dec_discon = true);
	~hit_reconstruct();
	bool fill_data(string data_filename);
	bool fill_noise(string noise_filename, double max_rms = 200.);
	void complete_adding();
	void set_noise_sca(int st, int sp);
	void set_data_sca(int st, int sp);
	void set_sigma(int s, int);
	bool is_display_warning = false;
	/// @brief Error msg for main function to show
	string error_msg = "";
	bool set_unmultiplexing_mapping(string filename, int offset = 0);
	bool is_init_done = false;
	

private:
	// For data out file
	bool init_dataout_file(string save_filename, bool is_write_txt);
	bool is_noise_fill = false;
	bool is_write_txt_file = false;
	hit_dec* dec_task;
	hit_dec* dec_task_cg; // Encoding with Complete graph
	hit_dec* dec_task_bg; // Encoding with Bipartite graph

	int noise_sca_start = 256;
	int noise_sca_stop = 511;
	double sigma = 5.;
	int data_sca_start = 0;
	int data_sca_stop = 255;

	bool is_save_sca = false;

	string encoding_filename;
	bool load_encoding_idx(string idx_filename);
	vector<vector<int>> encoding_idx;
	vector<int> encoding_scheme_select;

	TFile* fec_hit_file;
	TTree* fec_hit_tree;
	int nhits;
	int trigger_id;
	int hit_strips[cMAX_DEC_HITS];
	int hit_asic_chn[cMAX_DEC_HITS];
	int hit_amp[cMAX_DEC_HITS];
	double hit_time[cMAX_DEC_HITS];
	double hit_max_position[cMAX_DEC_HITS];
	int dimension_idx[cMAX_DEC_HITS];
	int cluster_number;
	int cluster_size[cMAX_DEC_HITS];
	int cluster_holed_num[cMAX_DEC_HITS];
	//int chip_id[cMAX_DEC_HITS];
	//int chn_id[cMAX_DEC_HITS];
	//int adc_data[cMAX_HITS][cSCA_NUM];

	// For data in file
	int init_datain_tree(string);
	TFile* fec_data_file;
	TTree* fec_data_tree;
	//long time_stamp_data[cELINK_NUM];
	int nhits_data;
	int trigger_id_data;
	int chip_id_data[cMAX_HITS];
	int chn_id_data[cMAX_HITS];
	int adc_data_all[cMAX_HITS][cSCA_NUM];
	//TBranch* b_time_stamp_data[cELINK_NUM];
	TBranch* b_nhits_data;
	TBranch* b_trigger_id_data;
	TBranch* b_chip_id_data;
	TBranch* b_chn_id_data;
	TBranch* b_adc_data_all;

	ofstream txt_file;
	string txt_filename;

	void fill_txt_file();

	string log_filename;
	ofstream log_file;

	// Get chn data over threshold
	vector<int> single_event_hit_chn[cCHIP_NUM];
	vector<int> single_event_hit_amp[cCHIP_NUM];
	vector<double> single_event_hit_time[cCHIP_NUM];
	vector<double> single_event_max_position[cCHIP_NUM];
	
	bool board_data_cmp();
	vector<double> single_chn_cmp(int* chn_data);
	vector<double> single_chn_cmp(int* chn_data, int vth, double mean);
	double fec_noise_mean[cCHIP_NUM][cCHN_NUM];
	double fec_noise_rms[cCHIP_NUM][cCHN_NUM];
	int fec_vth[cCHIP_NUM][cCHN_NUM];
	int max_rms = 100;

	void hit_data_dec();

	bool is_unmultiplexing = false;
	string unmultiplexing_mapping_filename = "./default_unmultiplexing_mapping.csv";
	int unmultiplexing_offset = 0;
	int unmultiplexing_mapping_idx[cCHN_NUM];
	bool load_unmultiplexing_mappting_file(string filename);
	void unmultiplexing_mapping();
	
	// Leading edge fit for micro-TPC
	// 'y = f0+f_max/(1 + exp(-(t-t_s)/tau))','independent','t','coefficients',{'f0','f_max','t_s','tau'});
	static Double_t sigmoid_fun(Double_t* x, Double_t* par);
	vector<double> leading_fit(double* x, double* y, int len, double* init_par);
};


