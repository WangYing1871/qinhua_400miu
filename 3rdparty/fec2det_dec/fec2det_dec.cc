#include "fec2det_dec.hh"

using namespace std;

int main(int argc, char **argv)
{
	if (argc == 1)
	{
		cout << cRED << "At least -f parameter is needed!" << cRESET << endl;
		show_info(argv[0]);
		return -1;
	}
	bool arg_hit = false;
	// -f
	string data_filename;
	// -e
	string encoding_list_filename = "./encoding_list_64to512.csv";
	// -E
	string encoding_list_filename2 = "";
	// -i
	string idx_filename = "./group_lists.csv";
	// -n -N
	int noise_sca_start = 256;
	int noise_sca_stop = 511;
	// -s
	double sigma = 5.;
	// -d -D
	int data_sca_start = 0;
	int data_sca_stop = 255;
	// -T
	bool is_write_txt_file = false;
	// -o
	string save_file_suffix = "dec";
	// -O save_path_set
	string save_path_set = "data_dec";
	// -m
	int max_rms = 100;
	// -v
	bool is_display_warning = false;
	// -F Noise file
	string noise_filename;
	// -C
	bool is_allow_dec_discontinue = true;
	// -U
	bool is_unencoding = false;
	// -Z
	int uncoding_offset = 0;
	while ((c = getopt(argc, argv, "f:e:E:i:n:N:s:d:D:t:To:O:m:vF:CUZ:")) != EOF)
	{
		arg_hit = true;
		switch (c)
		{
		case 'f':
			data_filename = optarg;
			break;
		case 'e':
			encoding_list_filename = optarg;
			break;
		case 'E':
			encoding_list_filename2 = optarg;
			break;
		case 'i':
			idx_filename = optarg;
			break;
		case 'n':
			noise_sca_start = atoi(optarg);
			break;
		case 'N':
			noise_sca_stop = atoi(optarg);
			break;
		case 's':
			sigma = atof(optarg);
			break;
		case 'd':
			data_sca_start = atoi(optarg);
			break;
		case 'D':
			data_sca_stop = atoi(optarg);
			break;
		case 'T':
			is_write_txt_file = true;
			break;
		case 'o':
			save_file_suffix = optarg;
			break;
		case 'O':
			save_path_set = optarg;
			break;
		case 'm':
			max_rms = atoi(optarg);
			break;
		case 'v':
			is_display_warning = true;
			break;
		case 'F':
			noise_filename = optarg;
			break;
		case 'C':
			is_allow_dec_discontinue = false;
			break;
		case 'U':
			is_unencoding = true;
			break;
		case 'Z':
			uncoding_offset = atoi(optarg);
			break;
		case '?':
			cout << cRED << "<<----------\?\?---------->>" << endl;
			cout << "unknow parameters " << argv[optind] << endl;
			cout << "<<----------\?\?---------->>" << endl;
			cout << "Or forget break in the last case???" << cRESET << endl;
			show_info(argv[0]);
			return -1;
			break;
		default:
			cout << cRED << "<<----------\?\?---------->>" << endl;
			cout << "Parameter error at No. " << optind << ", value: " << argv[optind - 1] << " " << argv[optind] << endl;
			cout << "<<----------\?\?---------->>" << cRESET << endl;
			show_info(argv[0]);
			return -1;
			break;
		}
	}
	if (!arg_hit)
	{
		cout << cRED << "No parameter" << cRESET << endl;
		show_info(argv[0]);
		return -1;
	}
	if (data_filename.empty())
	{
		cout << cRED << "No data file" << cRESET << endl;
		show_info(argv[0]);
		return -1;
	}
	vector<string> data_filename_path = split_string(data_filename, "/");
	string save_path;
	if (data_filename_path[1].empty())
	{
		save_path = "./" + save_path_set;
	}
	else
	{
		save_path = data_filename_path[1] + "/" + save_path_set;
	}
	if (stat(save_path.c_str(), &info) == -1)
	{
		cout << cYELLOW << "Save path : " << save_path << " does not exist, create it." << cRESET << endl;
		fs::create_directories(save_path);
	}
	string relative_filename = data_filename_path[0];
	vector<string> data_file_prefix = split_string(relative_filename, ".root");

	if (data_file_prefix.size() != 2)
	{
		cout << cRED << "Filename does not contain .root" << cRESET << endl;
		return -1;
	}
	if (is_allow_dec_discontinue)
	{
		cout << cYELLOW << "Allow one channel discontinue in position dec" << cRESET << endl;
		cout << cYELLOW << "允许一条解码不连续" << cRESET << endl;
		save_file_suffix = save_file_suffix + "_discon_";
	}
	else
	{
		cout << cBLUE << "The position dec result must be continuous" << cRESET << endl;
		cout << cBLUE << "解码必须连续" << cRESET << endl;
		save_file_suffix = save_file_suffix + "_con_";
	}
	vector<string> name_tmp = split_string(data_file_prefix[1], "_mt");
	string save_filename = save_path + name_tmp[1] + "_" + save_file_suffix + to_string(int(sigma)) + "sigma" + name_tmp[0];

	// Get the name part and remove the numnber after it "-"
	vector<string> name_tmp2 = split_string(name_tmp[0], "-");
	string log_filename = save_path + name_tmp2[1] + save_file_suffix + to_string(int(sigma)) + "sigma" + name_tmp[0];
	hit_reconstruct *fec2det = new hit_reconstruct(
      save_filename
      ,log_filename
      ,is_write_txt_file
      ,encoding_list_filename
      ,idx_filename
      ,encoding_list_filename2
      ,is_allow_dec_discontinue);
	if (!fec2det->is_init_done)
	{
		cout << cRED << "Error: " << fec2det->error_msg << cRESET << endl;
		return -1;
	}
	if (!noise_filename.empty())
	{
		if (!fec2det->fill_noise(noise_filename))
			return -1;
	}
	else
	{
		cout << cYELLOW << "Warning: ";
		cout << "No noise file, use default noise. In some peaking time config this may cause wrong result." << cRESET << endl;
	}
	fec2det->set_noise_sca(noise_sca_start, noise_sca_stop);
	fec2det->set_data_sca(data_sca_start, data_sca_stop);
	fec2det->set_sigma(sigma, max_rms);
	fec2det->is_display_warning = is_display_warning;
	if (!fec2det->fill_data(data_filename))
	{
		cout << cRED << "Error: " << fec2det->error_msg << cRESET << endl;
		return -1;
	}
	fec2det->complete_adding();
	cout << "Done" << endl;
	for (int i = 0; i < 10; i++)
	{
		cout << "--";
	}
	cout << endl;
}

void show_info(char *name)
{
	for (int i = 0; i < 30; i++)
		cout << "-";
	cout << endl;
	cout << cGREEN << "Usage: " << name << " -[parameter]" << endl
		 << cBLUE << "Parameter lists:" << cGREEN << endl
		 << "    -f [First data filename] (required) <待处理文件名>" << endl
		 << cCYAN
		 << "    -F [Threshold filename (*.root). If not specificated, the vth will be calculated from the input data] <阈值文件名，如果不指定，则从数据中计算>" << endl
		 << "    -e [Encoding table] (Default: ./encoding_list_64to512.csv) <编码表>" << endl
		 << "    -E [Secondary encoding table] (Default: ./encoding_list_64to512.csv) <二级编码表>" << endl
		 << "    -i [ASIC index file] (Default: ./group_lists.csv) <ASIC索引表>" << endl
		 << "    -n [Noise SCA start. If -t is set, this parameter will be ignored] (Default: 256) <计算噪声的SCA起始位置，如果设置了-t参数，则忽略此参数>" << endl
		 << "    -N [Noise SCA end. If -t is set, this parameter will be ignored] (Default: 511) <计算噪声的SCA结束位置，如果设置了-t参数，则忽略此参数>" << endl
		 << "    -s [Sigma level. 5 sigma indicates a confident level of 99.9999%] (Default: 5) <计算阈值的sigma值，5 sigma表示99.9999%的置信度>" << endl
		 << "    -d [Data SCA start] (Default: 0) <信号的SCA起始位置>" << endl
		 << "    -D [Data SCA stop] (Default: 255) <信号的SCA结束位置>" << endl
		 << "    -T <Write data to txt file> <将数据写入txt文件>" << endl
		 << "    -o [Save file prefix] (Default: dec) <保存文件名补充后缀>" << endl
		 << "    -O [Save path set] (Default: data_dec) <保存路径>" << endl
		 << "    -m [Max vth value] (Default: 1000) <最大阈值>" << endl
		 << "    -v <Show decoding warning> <显示解码警告>" << endl
		 << "    -C <Do not allow discontinue position dec> <不允许通道间断解码>" << endl
		 << cRESET << endl;
}
vector<string> split_string(const string &str, string delims, int inc)
{
	std::vector<std::string> result;
	size_t i = str.rfind(delims, str.length());
	if (i == str.length() - 1 || i == string::npos)
	{
       
		result.push_back(str);
	}
	else
	{
		result.push_back(str.substr(i + inc, str.length() - i));
		result.push_back(str.substr(0, i));
	}
	return result;
}
