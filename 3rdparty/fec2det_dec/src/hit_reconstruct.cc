#define info_out(X) std::cout<<"==> "<<__LINE__<<" "<<#X<<" |"<<(X)<<"|\n"
#include "hit_reconstruct.hh"

hit_reconstruct::hit_reconstruct(
    string save_filename
    ,string l_filename
    ,bool is_write_txt
    ,string enc_filename
    ,string a_idx_filename, string enc_filename2, bool is_allow_dec_discon)
{
	log_filename = l_filename + "_log.txt";
	log_file.open(log_filename, ios_base::app);
	init_dataout_file(save_filename, is_write_txt);
	// 初始化解码程序，需要先载入解码表，后续直接运行即可，可以加快程序运行速度
	dec_task_cg = new hit_dec(enc_filename, 40, 512, 64);
	dec_task_cg->is_allow_discontinue = is_allow_dec_discon;
	if (dec_task_cg->is_list_get)
	{
		log_file << "Encoding lists " << enc_filename << " get successfully" << endl;
	}
	else
	{
		log_file << "[Error]: Encoding lists " << enc_filename << " get failed!" << endl;
		is_init_done = false;
		error_msg = "Encoding lists " + enc_filename + " get failed!";
		return;
	}
	// 如果有第二个解码表，载入第二个解码表，这段代码是当存在两个编码方案混合测试的时候才需要的
	if (!enc_filename2.empty())
	{
		dec_task_bg = new hit_dec(enc_filename2, 40, 512, 64);
		dec_task_bg->is_allow_discontinue = is_allow_dec_discon;
		if (dec_task_bg->is_list_get)
		{
			log_file << "Encoding lists " << enc_filename2 << " get successfully" << endl;
		}
		else
		{
			log_file << "[Error]: Encoding lists " << enc_filename2 << " get failed!" << endl;
			is_init_done = false;
			error_msg = "Encoding lists " + enc_filename2 + " get failed!";
			return;
		}
	}

	if (!load_encoding_idx(a_idx_filename))
	{
		log_file << "[Warning]: Not specify the ASIC index file. Exit!!!" << endl;
		error_msg = "Not specify the ASIC index file. Exit!!!";
		return;
	}
	time_t now = time(0);
	for (int i = 0; i < 20; i++)
	{
		log_file << "--";
	}
	log_file << endl
			 << "Excute time: " << ctime(&now) << endl;
	is_init_done = true;
}

hit_reconstruct::~hit_reconstruct()
{
	if (fec_data_file->IsOpen())
	{
		fec_data_file->Close();
		delete fec_data_file;
	}
	if (fec_hit_file->IsOpen())
	{
		fec_hit_file->Close();
		delete fec_hit_tree;
		delete fec_hit_file;
	}
	if (txt_file.is_open())
	{
		txt_file.close();
	}
	if (log_file.is_open())
	{
		log_file.close();
	}
}

bool hit_reconstruct::fill_data(string data_filename)
{
	int event_num = init_datain_tree(data_filename);
	if (event_num <= 0)
	{
		error_msg = "Event number is " + to_string(event_num) + ". Return!";
		return false;
	}
	int one_percent = event_num / 100;
	one_percent = (one_percent == 0 ? 1 : one_percent);
	//cout << "Processing \t 0%" << setprecision(4);
	setbuf(stdout, NULL);
	for (int i = 0; i < event_num; i++)
	{
		// The following code is used to show the progress of the program
		if (i % one_percent == 0 && i != 0)
		{

			cout << cCLR_LINE;
			cout << "Processing "
				 << "\t" << double(i) / event_num * 100 << "% ";
		}

		fec_data_tree->GetEntry(i);
		trigger_id = trigger_id_data;

		// nhits is the number of hits in one event
		nhits = 0;
		cluster_number = 0;
		if (!board_data_cmp())
		{
			// If the board data is not over the threshold, still Fill it inorder to keep the trigger ID synchronized
			fec_hit_tree->Fill();
			continue;
		}

		hit_data_dec();
		fec_hit_tree->Fill();
		if (is_write_txt_file)
		{
			fill_txt_file();
		}
	}
	cout << "\b\b\b\b\b\b\b\b 100%    " << endl;
	return true;
}

bool hit_reconstruct::fill_noise(string noise_filename, double max_rms)
{
	TFile *noise_file;
	TTree *noise_tree;
	TBranch *b_mean;
	TBranch *b_rms;
	double mean[4][68];
	double rms[4][68];
	noise_file = new TFile(noise_filename.c_str());
	if (noise_file->IsZombie() || !noise_file->IsOpen())
	{
		cout << cRED << noise_filename << " does not exist!!!" << cRESET << endl;
		log_file << "[Error]: " << noise_filename << " does not exist!!!" << endl;
		delete noise_file;
		is_noise_fill = false;
		return false;
	}
	noise_tree = (TTree *)noise_file->Get("noise_tree");
	noise_tree->SetBranchAddress("mean", mean, &b_mean);
	noise_tree->SetBranchAddress("rms", rms, &b_rms);
	if (noise_tree->GetEntries() == 0)
	{
		cout << cRED << "Noise tree in " << noise_filename << " does not have data" << cRESET << endl;
		log_file << "[Error]: "
				 << "Noise tree in " << noise_filename << " does not have data" << endl;
		is_noise_fill = false;
		return false;
	}
	noise_tree->GetEntry(0);
	for (int j = 0; j < cCHIP_NUM; j++)
	{
		for (int k = 0; k < cCHN_NUM; k++)
		{
			fec_noise_mean[j][k] = mean[j][k];
			if (rms[j][k] < max_rms)
			{
				fec_noise_rms[j][k] = rms[j][k];
				fec_vth[j][k] = mean[j][k] + sigma * rms[j][k];
			}
			else
			{
				fec_noise_rms[j][k] = max_rms;
				fec_vth[j][k] = mean[j][k] + sigma * max_rms;
			}
		}
	}
	delete noise_file;
	cout << cGREEN << "Noise file in " << noise_filename << " read done!" << cRESET << endl;
	is_noise_fill = true;
	return true;
}

void hit_reconstruct::complete_adding()
{
	fec_hit_file->cd();
	fec_hit_tree->Write("", TObject::kOverwrite);
	fec_hit_file->Flush();
	fec_hit_file->Close();
	log_file.flush();
	log_file.close();
}

void hit_reconstruct::set_noise_sca(int st, int sp)
{
	if (st > sp || st < 0 || sp > 511)
	{
		cout << cYELLOW << "[Warning]: Start point is larger than the stop point!" << cRESET << endl;
		log_file << "[Warning]: Start point is larger than the stop point!" << endl;
		noise_sca_start = 256;
		noise_sca_stop = 511;
	}
	else
	{
		noise_sca_start = st;
		noise_sca_stop = sp;
	}
}

void hit_reconstruct::set_data_sca(int st, int sp)
{
	if (st > sp || st < 0 || sp > 511)
	{
		cout << cYELLOW << "[Warning]: Start point is larger than the stop point!" << cRESET << endl;
		log_file << "[Warning]: Start point is larger than the stop point!" << endl;
		data_sca_start = 0;
		data_sca_stop = 256;
	}
	else
	{
		data_sca_start = st;
		data_sca_stop = sp;
		cout << cCYAN << "Data sca start: " << data_sca_start << " stop: " << data_sca_stop << cRESET << endl;
	}
}

void hit_reconstruct::set_sigma(int s, int m)
{
	if (s < 0)
	{
		cout << cYELLOW << "[Warning]:  Sigma is lower than zero!" << cRESET << endl;
		log_file << "[Warning]:  Sigma is lower than zero!" << endl;
		sigma = 5;
	}
	else
	{
		sigma = s;
	}
	if (m < 0)
	{
		cout << cYELLOW << "[Warning]: Max vth is lower than 0. Set to 100!" << cRESET << endl;
		log_file << "[Warning]: Max vth is lower than 0. Set to 100!" << endl;
		m = 100;
	}
	else
	{
		max_rms = m;
		if (m < 20)
		{
			cout << cYELLOW << "Max vth is quite small, and this may introduce noise" << cRESET << endl;
			log_file << "Max vth is quite small, and this may introduce noise" << endl;
		}
	}
}

bool hit_reconstruct::set_unmultiplexing_mapping(string filename, int offset)
{
	is_unmultiplexing = true;
	if (!load_unmultiplexing_mappting_file(filename))
	{
		return false;
	}
	unmultiplexing_offset = offset;
	cout << cGREEN << "Set multiplexing mapping file to " << filename << " with offset " << offset << cRESET << endl;
	return true;
}

bool hit_reconstruct::init_dataout_file(string save_filename, bool is_write_txt)
{
	string root_filename = save_filename + ".root";
	fec_hit_file = new TFile(root_filename.c_str(), "recreate");
	fec_hit_tree = new TTree("asic_hit_data", "FEC HIT");

	fec_hit_tree->Branch("nhits", &nhits, "nhits/I");
	fec_hit_tree->Branch("trigger_id", &trigger_id, "trigger_id/I");
	fec_hit_tree->Branch("hit_strips", hit_strips, "hit_strips[nhits]/I");
	fec_hit_tree->Branch("hit_asic_chn", hit_asic_chn, "hit_asic_chn[nhits]/I");
	fec_hit_tree->Branch("hit_amp", hit_amp, "hit_amp[nhits]/I");
	fec_hit_tree->Branch("hit_time", hit_time, "hit_time[nhits]/D");
	fec_hit_tree->Branch("hit_max_position", hit_max_position, "hit_max_position[nhits]/D");
	fec_hit_tree->Branch("dimension_idx", dimension_idx, "dimension_idx[nhits]/I");
	fec_hit_tree->Branch("cluster_number", &cluster_number, "cluster_number/I");
	fec_hit_tree->Branch("cluster_size", cluster_size, "cluster_size[cluster_number]/I");
	fec_hit_tree->Branch("cluster_holed_num", cluster_holed_num, "cluster_holed_num[cluster_number]/I");
	// fec_hit_tree->Branch("chip_id", chip_id, "chip_id[nhits]/I");
	// fec_hit_tree->Branch("chn_id", chn_id, "chn_id[nhits]/I");
	// fec_hit_tree->Branch("adc_data", adc_data, "adc_data[nhits][512]"); // The second idx of adc_data (i.e. 512) must be equal to cSCA_NUM
	cout << cCYAN << "TTree in TFile " << cMAGENTA << root_filename << cCYAN << " init done" << cRESET << endl;
	log_file << "TTree in TFile " << root_filename << " init done" << endl;

	if (is_write_txt)
	{
		txt_filename = save_filename + ".txt";
		txt_file.open(txt_filename);
		if (txt_file.is_open())
		{
			cout << cCYAN << "Save hit position txt file to " << cMAGENTA << txt_filename << cRESET << endl;
			log_file << "Save hit position txt file to " << txt_filename << endl;
			is_write_txt_file = true;
		}
		else
		{
			cout << cRED << "Cannot create file: " << txt_filename << cRESET << endl;
			log_file << "[Error]: "
					 << "Cannot create file: " << txt_filename << endl;
			is_write_txt_file = false;
		}
	}
	return true;
}

bool hit_reconstruct::load_encoding_idx(string idx_filename)
{
	ifstream idx_file;
	idx_file.open(idx_filename);
	if (idx_file.is_open())
	{
		string list_s;
		vector<int> idx_s;
		int line = 0;
		while (getline(idx_file, list_s))
		{
			string s;
			istringstream list_ss(list_s);
			while (list_ss)
			{
				if (!getline(list_ss, s, ',') || s == "" || s == "\r")
				{
					break;
				}
				// Check to insure that only 4 ASIC are used

				if (stoi(s) < 0 || stoi(s) > 4)
				{
					cout << cRED << "Index file " << idx_filename << " content error at line " << line << cRESET << endl;
					log_file << "[Error]: "
							 << "Index file " << idx_filename << " content error at line " << line << endl;
					return false;
				}
				idx_s.push_back(stoi(s));
			}
			if (idx_s.size() > 5 || line > 3)
			{
				cout << cRED << "Index file " << idx_filename << " content error at line " << line << cRESET << endl;
				log_file << "[Error]: "
						 << "Index file " << idx_filename << " content error at line " << line << endl;
				return false;
			}
			vector<int> idx_tmp;
			idx_tmp.insert(idx_tmp.end(), idx_s.begin(), idx_s.end() - 1);
			encoding_idx.push_back(idx_tmp);
			encoding_scheme_select.push_back(idx_s.back());
			idx_s.clear();
			line++;
		}
	}
	else
	{
		cout << cRED << "Cannot read ASIC index file-- " << idx_filename << cRESET << endl;
		log_file << "[Error]: "
				 << "Cannot read ASIC index file-- " << idx_filename << endl;
		return false;
	}
	log_file << "Load encoding index file " << idx_filename << endl;

  //info_out("encoding_idx:");
  //for(auto&& x : encoding_idx){
  //  for (auto&& y : x) std::cout<<y<<" ";
  //  std::cout<<std::endl;
  //}
	return true;
}

int hit_reconstruct::init_datain_tree(string data_filename)
{
	fec_data_file = new TFile(data_filename.c_str());
	if (fec_data_file->IsZombie() || !fec_data_file->IsOpen())
	{
		cout << cRED << "File" << data_filename << " open failed. Please check the file" << cRESET << endl;
		log_file << "[Error]: "
				 << "File" << data_filename << " open failed. Please check the file" << endl;
		return -1;
	}
	fec_data_tree = (TTree *)fec_data_file->Get("fec_origin_data");
	if (fec_data_tree == nullptr)
	{
		cout << cRED << "Cannot get TTree* fec_origin_data in file \"" << data_filename << "\"" << cRESET << endl;
		log_file << "[Error]: "
				 << "Cannot get TTree * fec_origin_data in file \"" << data_filename << "\"" << endl;
		return -1;
	}
	fec_data_tree->SetBranchAddress("trigger_id", &trigger_id_data, &b_trigger_id_data);
	fec_data_tree->SetBranchAddress("nhits", &nhits_data, &b_nhits_data);
	fec_data_tree->SetBranchAddress("chip_id", chip_id_data, &b_chip_id_data);
	fec_data_tree->SetBranchAddress("chn_id", chn_id_data, &b_chn_id_data);
	fec_data_tree->SetBranchAddress("adc_data", &adc_data_all[0][0], &b_adc_data_all);
	cout << cGREEN << "Good! TTree in TFile: " << data_filename << " get done!" << cRESET << endl;
	log_file << "TTree in TFile: " << data_filename << " get done!" << endl;
	return fec_data_tree->GetEntries();
}

void hit_reconstruct::fill_txt_file()
{
	if (!txt_file.is_open())
	{
		cout << cRED << "txt file is not in open status. Please check the code" << cRESET << endl;
		return;
	}
	txt_file << "{" << endl;
	txt_file << "trigger_id : " << trigger_id << endl
			 << "nhits : " << nhits << endl
			 << "cluster_number: " << cluster_number << endl;
	int i = 0;
	for (int k = 0; k < cluster_number; k++)
	{
		txt_file << "-------- Cluster " << k << "--------" << endl;
		int base = i;
		txt_file << "    cluster size:\t" << cluster_size[k] << endl;
		txt_file << "    dimension_idx:\t" << dimension_idx[base] << endl;
		txt_file << "    hit_strips:\t\t";
		for (int j = 0; j < cluster_size[k]; j++)
		{
			txt_file << hit_strips[j + base];
			if (j != cluster_size[k] - 1)
			{
				txt_file << ","
						 << "\t";
			}
			i++;
		}
		txt_file << endl;
		txt_file << "    hit_asic_chn:\t";
		for (int j = 0; j < cluster_size[k]; j++)
		{
			txt_file << hit_asic_chn[j + base];
			if (j != cluster_size[k] - 1)
			{
				txt_file << ","
						 << "\t";
			}
		}
		txt_file << endl;
		txt_file << "    hit_amp:\t";
		for (int j = 0; j < cluster_size[k]; j++)
		{
			txt_file << hit_amp[j + base];
			if (j != cluster_size[k] - 1)
			{
				txt_file << ","
						 << "\t";
			}
		}
		txt_file << endl;
		txt_file << "    hit_time:\t";
		for (int j = 0; j < cluster_size[k]; j++)
		{
			txt_file << hit_time[j + base];
			if (j != cluster_size[k] - 1)
			{
				txt_file << ","
						 << "\t";
			}
		}
		txt_file << endl;
		txt_file << "    hit_max_position:\t";
		for (int j = 0; j < cluster_size[k]; j++)
		{
			txt_file << hit_max_position[j + base];
			if (j != cluster_size[k] - 1)
			{
				txt_file << ","
						 << "\t";
			}
		}
		txt_file << endl;
	}
	txt_file << endl
			 << "}" << endl;
}

bool hit_reconstruct::board_data_cmp()
{
	for (int i = 0; i < cCHIP_NUM; i++)
	{
		single_event_hit_chn[i].clear();
		single_event_hit_amp[i].clear();
		single_event_hit_time[i].clear();
	}
	bool is_hit = false;
	for (int i = 0; i < nhits_data; i++)
	{
		int chip = chip_id_data[i];
		int chn = chn_id_data[i];
		if (chip >= cCHIP_NUM || chn >= cCHN_NUM || chip < 0 || chn < 0)
		{
			cout << cRED << "Chip or chn ID error. Chip: " << chip << " chn: " << chn
				 << " trigger ID: " << trigger_id << " i: " << i << cRESET << endl;
			log_file << "[Error]: "
					 << "Chip or chn ID error. Chip: " << chip << " chn: " << chn
					 << " trigger ID: " << trigger_id << " i: " << i << endl;
			return false;
		}
		if (chn == 11 || chn == 22 || chn == 45 || chn == 56)
		{
			continue;
		}

		// hit_info 用于存储单次比较后的通道信息，按照顺序是：幅度，时间，最大值位置
		vector<double> hit_info;
		if (is_noise_fill)
		{
			hit_info = single_chn_cmp(&adc_data_all[i][0], fec_vth[chip][chn], fec_noise_mean[chip][chn]);
		}
		else
		{
			hit_info = single_chn_cmp(&adc_data_all[i][0]);
		}
		if (chn > 11 && chn < 22)
		{
			chn -= 1;
		}
		else if (chn > 22 && chn < 45)
		{
			chn -= 2;
		}
		else if (chn > 45 && chn < 56)
		{
			chn -= 3;
		}
		else if (chn > 56)
		{
			chn -= 4;
		}
		if (hit_info.size() != 0)
		{
			single_event_hit_chn[chip].push_back(chn);
			single_event_hit_amp[chip].push_back(hit_info[0]);
			single_event_hit_time[chip].push_back(hit_info[1]);
			single_event_max_position[chip].push_back(hit_info[2]);
			is_hit = true;
		}
	}
	return is_hit;
}

vector<double> hit_reconstruct::single_chn_cmp(int *chn_data)
{
	TH1F *chn_noise = new TH1F("Noise", "single_chn_noise", 100, 0, 600);
	for (int i = noise_sca_start; i < noise_sca_stop; i++)
	{
		chn_noise->Fill(chn_data[i]);
	}
	double pedestal = chn_noise->GetMean();
	double rms = chn_noise->GetStdDev();
	if (rms > max_rms)
	{
		rms = max_rms;
	}
	double vth = pedestal + sigma * rms;
	delete chn_noise;
	return single_chn_cmp(chn_data, vth, pedestal);
}

vector<double> hit_reconstruct::single_chn_cmp(int *chn_data, int vth, double mean)
{
	int max_value = 0;
	int max_position = 0;
	vector<double> hit_info;
	for (int i = data_sca_start; i < data_sca_stop; i++)
	{
		if (chn_data[i] > vth)
		{
			if (chn_data[i] > max_value)
			{
				max_value = chn_data[i];
				max_position = i;
			}
		}
	}
	if (max_value == 0)
	{
		return hit_info;
	}
	double half_maximun_position = 0;
	double half_maximun_value = max_value / 2;
	for (int i = data_sca_start; i < data_sca_stop - 1; i++)
	{
		if (chn_data[i] < half_maximun_value && chn_data[i + 1] > half_maximun_value)
		{
			half_maximun_position = i;
			break;
		}
	}
	// 这段代码暂时注释掉，是为了用来提取波形的前沿信息，用于Micro-TPC方法使用的，
	// 目前先采用最大值的半高位置作为时间信息，即前沿定时 half_maximun_position
	/* int start_position = 0;
	for (int i = max_position; i >= data_sca_start; i--)
	{
		if (chn_data[i] > mean && chn_data[i - 1] < mean)
		{
			start_position = i;
			break;
		}
	}
	int data_len = max_position - start_position + 10;
	double* xx = new double[data_len];
	double* yy = new double[data_len];
	for (int i = 0; i < data_len; i++)
	{
		xx[i] = i + start_position - 5;
		yy[i] = chn_data[i + start_position - 5];
	}
	double* init_par = new double[4];
	init_par[0] = mean;
	init_par[1] = max_value;
	init_par[2] = (max_position + start_position) / 2.;
	init_par[3] = 2;
	vector<double> fit_result = leading_fit(xx, yy, data_len, init_par);
	if (fit_result[6] < 0.7)
	{
		return hit_info;
	}*/
	hit_info.push_back(max_value - mean);
	hit_info.push_back(half_maximun_position);
	hit_info.push_back(max_position);
	return hit_info;
}

void hit_reconstruct::hit_data_dec()
{
	// 每次运行前清空nhits 和 cluser_number，以免后续程序判断错误
	nhits = 0;
	cluster_number = 0;
	// Clear the data stored in root for not confusing while in TBrowser()
	for (int i = 0; i < cMAX_DEC_HITS; i++)
	{
		hit_strips[i] = 0;
		hit_asic_chn[i] = 0;
		hit_amp[i] = 0;
		hit_time[i] = 0;
		hit_max_position[i] = 0;
		dimension_idx[i] = 0;
		cluster_size[i] = 0;
		cluster_holed_num[i] = 0;
	}

	// encoding_idx存放的是同一个维度下的ASIC的索引，在测试中可能是两个芯片对应同一个
	// 维度，甚至是三个芯片对应同一个维度,encoding_idx指定了ASIC的顺序，在下面的处理
	// 程序中将同一个维度对应的芯片的数据合并，并且指定是第几个芯片的数据
  //info_out(encoding_idx.size()); exit(0);
	for (int i = 0; i < encoding_idx.size(); i++)
	{
		vector<int> sub_idx = encoding_idx[i];
		vector<int> aget_chn;
		vector<int> aget_amp;
		vector<double> signal_time;
		vector<double> signal_max_position;
		vector<int> asic_idx;
		for (int j = 0; j < sub_idx.size(); j++) {
			aget_chn.insert(aget_chn.end(),
          single_event_hit_chn[sub_idx[j]].begin()
          ,single_event_hit_chn[sub_idx[j]].end());
			aget_amp.insert(aget_amp.end()
          ,single_event_hit_amp[sub_idx[j]].begin()
          ,single_event_hit_amp[sub_idx[j]].end());
			signal_time.insert(signal_time.end()
          ,single_event_hit_time[sub_idx[j]].begin()
          ,single_event_hit_time[sub_idx[j]].end());
			signal_max_position.insert(signal_max_position.end()
          ,single_event_max_position[sub_idx[j]].begin()
          , single_event_max_position[sub_idx[j]].end());
			for (int k = 0; k < single_event_hit_chn[sub_idx[j]].size(); k++) asic_idx.push_back(j);
    }

    //for (auto&& x : asic_idx)  info_out(x);

    // ??
		if (encoding_scheme_select[i] == 0) dec_task = dec_task_cg;
		else dec_task = dec_task_bg; 
		if (dec_task->run_dec(aget_chn, asic_idx))
		{
			int base = nhits;
			int dec_strip_cnt = 0;
			for (int j = 0; j < dec_task->cluster_num.size(); j++)
			{
				// Remove the oxo or oxoxo or oxoxoxo... cluster, where o is the hitted strips and x is the hole in the cluster
				if ((dec_task->cluster_num[j] - dec_task->cluster_holed[j] == 1) 
            && dec_task->cluster_num[j] != 1)
				{
					dec_strip_cnt += dec_task->cluster_num[j];
					continue;
				}
				for (int k = 0; k < dec_task->cluster_num[j]; k++)
				{
					hit_strips[nhits] = dec_task->hit_strip[k + dec_strip_cnt];
					hit_asic_chn[nhits] = aget_chn[dec_task->hit_seq[dec_strip_cnt + k]];
					hit_amp[nhits] = aget_amp[dec_task->hit_seq[k + dec_strip_cnt]];
					hit_time[nhits] = signal_time[dec_task->hit_seq[k + dec_strip_cnt]];
					hit_max_position[nhits] = signal_max_position[dec_task->hit_seq[k + dec_strip_cnt]];
					dimension_idx[nhits] = i;
					nhits++;
				}
				
				cluster_size[cluster_number] = dec_task->cluster_num[j];
				cluster_holed_num[cluster_number] = dec_task->cluster_holed[j];
				cluster_number++;
				dec_strip_cnt += dec_task->cluster_num[j];
			}
		}
	}
	
}

bool hit_reconstruct::load_unmultiplexing_mappting_file(string filename)
{
	fstream mapping_file;
	mapping_file.open(filename);
	if (!mapping_file.is_open())
	{
		cout << cRED << "Cannot open mapping file: " << filename << cRESET << endl;
		log_file << "[Error]: "
				 << "Cannot open mapping file: " << filename << endl;
		return false;
	}
	string line;
	int line_num = 0;
	while (getline(mapping_file, line) && line_num < cCHN_NUM)
	{
		// check if line is int
		if (line.find_first_not_of("0123456789") != string::npos)
		{
			cout << cRED << "Mapping file: " << filename << " content error at line " << line_num << cRESET << endl;
			log_file << "[Error]: "
					 << "Mapping file: " << filename << " content error at line " << line_num << endl;
			return false;
		}
		unmultiplexing_mapping_idx[line_num] = stoi(line);
	}
	return true;
}

void hit_reconstruct::unmultiplexing_mapping()
{
	for (int i = 0; i < cCHIP_NUM; i++)
	{
	}
}

Double_t hit_reconstruct::sigmoid_fun(Double_t *x, Double_t *par)
{
	return par[0] + par[1] / (1 + TMath::Exp(-(x[0] - par[2]) / par[3]));
}

vector<double> hit_reconstruct::leading_fit(double *x, double *y, int len, double *init_par)
{
	TF1 *fit_fun = new TF1("sigmoid_fun", sigmoid_fun, x[0] - 1, x[len - 1] + 1, 4);
	fit_fun->SetParameters(init_par[0], init_par[1], init_par[2], init_par[3]);
	fit_fun->SetParNames("f_0", "f_max", "t0", "tau");
	TGraph *gr = new TGraph(len, x, y);
	TFitResultPtr r = gr->Fit(fit_fun, "SQ");
	double *fit_result = fit_fun->GetParameters();
	double chi = fit_fun->GetChisquare();
	double ndf = fit_fun->GetNDF();
	double corr = gr->GetCorrelationFactor();
	vector<double> result;
	for (int i = 0; i < 4; i++)
	{
		result.push_back(fit_result[i]);
	}
	result.push_back(chi);
	result.push_back(ndf);
	result.push_back(corr);
	return result;
}
