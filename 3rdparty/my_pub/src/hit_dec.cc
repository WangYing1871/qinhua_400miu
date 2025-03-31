#include <algorithm>
#include "hit_dec.hh"
#ifndef info_out
#define info_out(X) std::cout<<"==> "<<__LINE__<<" "<<#X<<" |"<<(X)<<"|\n"
#endif

hit_dec::hit_dec(
    string encoding_lists_filename
    , int max_chn // 40
    , int det_chn // 512
    , int electronic_chn_num //64
    , bool is_dis_warning)
{
    lists_filename = encoding_lists_filename;
    max_hit_chn = max_chn;
    det_chn_num = det_chn;
    is_list_get = get_encoded_lists(lists_filename, electronic_chn_num);
    is_display_warning = is_dis_warning;
    cout << cGREEN << "Using encoding list: " << lists_filename << endl << cRESET;
}

bool hit_dec::run_dec(vector<int> chn, vector<int> idx)
{
    if (!is_list_get)
    {
        cout << cRED << "Hit lists not get please check the code!!! " << cRESET << endl;
        return false;
    }
    if (chn.size() == 0)
    {
        return false;
    }
    if (chn.size() > max_hit_chn)
    {
        if (is_display_warning)
        {
            cout << cYELLOW << "Total hit number of channels is too large: " << chn.size()
                << " (set max to " << max_hit_chn << "). Considering spark event" << endl << cRESET;
        }
        return false;
    }
    if (chn.size() != idx.size())
    {
        cout << cRED << "Variable size error. chn size: " << chn.size() << " idx size: " << idx.size() << endl << cRESET << endl;
        return false;
    }
    clear();
    /*
    if (chn.size() == 1) 
    {
        hit_strip.push_back(-chn[0]);
        hit_seq.push_back(0);
        hit_idx.push_back(idx[0]);
        return true;
    }
    */
    return hit_dec_total(chn, idx) != 0;
}

hit_dec::~hit_dec()
{
}

void hit_dec::clear()
{
    hit_strip.clear();
    hit_seq.clear();
    hit_idx.clear();
    cluster_num.clear();
    cluster_holed.clear();
    strip_all.clear();
}

bool hit_dec::get_encoded_lists(string filename, int electronic_chn_num)
{
    ifstream detector_list_file;
    detector_list_file.open(filename);
    string list_s;
    vector<string> list_array;
    if (detector_list_file.is_open())
    {
        // Read the title line
        //getline(detector_list_file, list_s);
        int line = 0;
        while (getline(detector_list_file, list_s))
        {
            istringstream list_ss(list_s);
            info_out(list_s);
            while (list_ss)
            {
                string s;
                if (!getline(list_ss, s, ',') || s == "" || s == "\r")
                    break;
                list_array.push_back(s);
            }
            for (int i = 1; i < list_array.size(); i++)
            {
                encoded_lists[line].push_back(stoi(list_array[i]));
            }
            line++;
            list_array.clear();
        }
        if (line != electronic_chn_num)
        {
            cout << cYELLOW << filename << " content error!!!" << cRESET << endl;
            detector_list_file.close();
            return false;
        }
    }
    else
    {
        cout << "File " << filename << " open error!!!" << endl;
        //detector_list_file.close();
        return false;
    }
    detector_list_file.close();

    //for (auto&& x : encoded_lists){
    //  info_out(x.size());
    //  if (x.size());
    //  for (auto&& y : x) std::cout<<y<<" ";
    //  std::cout<<std::endl;
    //}
    return true;
}

int hit_dec::hit_dec_total(vector<int> chn, vector<int> idx)
{
  info_out(chn.size());
  for(auto&& x : chn) info_out(x);
  /*
  if (*std::max_element(chn.begin(),chn.end())>383
        || *std::min_element(chn.begin(),chn.end())<0){
    //...
    return 0;
  }
  std::vector<int> idx_tmp;
  std::vector<int> seq_tmp;
  for (int i=0; i<chn.size(); ++i){
    for (int j=0; j<encoded_lists[chn[i]].size(); ++j){
      strip_all.emplace_back(encoded_lists[chn[i]][j]+idx[i]*det_chn_num);
      seq_tmp.emplace_back(i);
      idx_tmp.emplace_back(idx[i]);
    }
  }
  */


  //info_out(*max_element(chn.begin(),chn.end()));
  //info_out(*min_element(chn.begin(),chn.end()));
  // ??? 
    if ((*max_element(chn.begin(), chn.end()) > 383) || (*min_element(chn.begin(), chn.end()) < 0))
    {
        cout << cRED << "Channel reconstrut error" << cRESET << endl;
        return 0;
    }
    
    vector<int> idx_tmp;
    vector<int> seq_tmp;
    for (int i = 0; i < chn.size(); i++)
    {
      //info_out(encoded_lists[chn[i]].size());
        for (int j = 0; j < encoded_lists[chn[i]].size(); j++)
        {
          //info_out(idx[i]);
            // Fill det_chn with all the possible encoding chn
            strip_all.push_back(encoded_lists[chn[i]][j] + idx[i] * det_chn_num);
            seq_tmp.push_back(i);
            idx_tmp.push_back(idx[i]);
        }
    }

    // 如果是单根条击中的事例，将所有可能的解码结果返回，因此cluster_num设置为1，自然cluster_holed也为0
    if (chn.size() == 1)
    {
        hit_strip.insert(hit_strip.end(), strip_all.begin(), strip_all.end());
        hit_idx.insert(hit_idx.end(), idx_tmp.begin(), idx_tmp.end());
        hit_seq.insert(hit_seq.end(), seq_tmp.begin(), seq_tmp.end());
        for (int k = 0; k < strip_all.size(); k++)
        {
            cluster_num.push_back(1);
            cluster_holed.push_back(0);
        }
        return cluster_num.size();
    }

    vector<int> strip_seq = bubble_sort_chn(strip_all);
    for(auto&& x : strip_all){
      info_out(x);
    }
    resort(idx_tmp, strip_seq);
    resort(seq_tmp, strip_seq);
    for(auto&& x : idx_tmp) info_out(x);
    for(auto&& x : seq_tmp) info_out(x);
    vector<int> continuous_num = calc_continuation(strip_all);
    vector<int> nno_posi;// Record the non-one element position of vector continuous_num
    int cnt = 0;
    // continuous_num.size() is 
    //for (auto&& x : continuous_num){
    //  info_out(x);
    //}
    while (cnt < continuous_num.size())
    {
        if (continuous_num[cnt] == 1)
        {
            cnt++;
            continue;
        }
        else
        {
            nno_posi.push_back(cnt);
            cnt += continuous_num[cnt];
        }
    }
    if (nno_posi.size() == 0)
    {
        return 0;
    }
    //for(auto&& x : nno_posi){
    //  info_out(x);
    //}
    //exit(0);

    // Store the continuous hit strips
    int holed_num = 0;
    for (int i = 0; i < nno_posi.size(); i++)
    {
        int posi_start = nno_posi[i];
        int last_strip = strip_all[posi_start];
        //info_out(last_strip); exit(0);
        for (int j = 0; j < continuous_num[posi_start]; j++) 
        {
          info_out(posi_start);
          info_out(strip_all[posi_start+j]-1);
          info_out(last_strip);
            if (j != 0)
            {
                if (last_strip != (strip_all[posi_start + j] - 1))
                {
                    holed_num++;
                }
                last_strip = strip_all[posi_start + j];
            }
            hit_strip.push_back(strip_all[posi_start + j]);
            hit_idx.push_back(idx_tmp[posi_start + j]);
            hit_seq.push_back(seq_tmp[posi_start + j]);
        }
        //info_out(continuous_num[posi_start]); exit(0);
        cluster_num.push_back(continuous_num[posi_start]);
        cluster_holed.push_back(holed_num);
        holed_num = 0;
    }
    return cluster_num.size();
   
}

/// <summary>
/// Bubble sort the input array and return the index of the sorted array
/// </summary>
/// <param name="d_chn">The array that need be sorted, with reference parammeter input so that this array will be sorted</param>
/// <returns>The index of the sorted result and can be used to change the order of other array</returns>
vector<int> hit_dec::bubble_sort_chn(vector<int>& d_chn)
{
    int len = d_chn.size();
    vector<int> _idx;
    for (int i = 0; i < len; i++)
    {
        _idx.push_back(i);
    }
    int tmp;
    for (int i = 0; i < len - 1; i++)
    {
        for (int j = 0; j < len - i - 1; j++)
        {
            if (d_chn[j] > d_chn[j + 1])
            {
                tmp = d_chn[j];
                d_chn[j] = d_chn[j + 1];
                d_chn[j + 1] = tmp;
                tmp = _idx[j];
                _idx[j] = _idx[j + 1];
                _idx[j + 1] = tmp;
            }
        }
    }
    return _idx;
}

vector<int> hit_dec::calc_continuation(vector<int> d_chn)
{
    if (d_chn.size() == 0)
    {
        return vector<int>();
    }
    int k = 0;
    int len = d_chn.size();
    vector<int> continuous_num;
    for (int i = 0; i < len; i++)
    {
        bool is_continue = true;
        int j = i;
        continuous_num.push_back(1);
        while (is_continue && (j < (len - 1)))
        {
            if (((d_chn[j] + 1) == d_chn[j + 1]) 
                || (is_allow_discontinue && ((d_chn[j] + 2) == d_chn[j + 1])))
            {
                continuous_num[i]++;
                j++;
            }
            else
            {
                is_continue = false;
            }
        }
    }
    continuous_num.push_back(1);
    return continuous_num;
}
