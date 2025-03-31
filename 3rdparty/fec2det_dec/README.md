# 电子学击中位置到探测器击中位置的解码

> 如果代码中中文乱码了，请使用gbk编码打开
> If the Chinese characters in the code are garbled, please open with gbk encoding

程序实现从读出的电子学数据中挑选出击中的通道，并且对其进行解码的功能，最终只存储解码成功的数据

## 输入参数
 + -f 待解码的数据文件名，需要经过[fec_data_dec程序处理](..\fec_data_dec)
 + -n 阈值文件名，如果不指定的话，会从数据中选一段来计算噪声水平。建议采用指定的方式，否则在成形时间较长的情况下，计算的阈值会有问题。指定的文件需要在采数前空采一段数据，然后使用[fec_noise_calc](..\fec_noise_calc)程序计算得到。
 + -n 如果不指定阈值文件，那么就从数据中选取一段来计算噪声水平，-n表示这段数据的起始SCA位置
 + -N 如果不指定阈值文件，那么就从数据中选取一段来计算噪声水平，-N表示这段数据的结束SCA位置
 + -e 指定采用的编码表，目前默认采用的是完全图的欧拉回路构造的编码表，文件名为encoding_list_64to512.csv
 + -E 指定第二个编码表，默认为完全二分图的欧拉回路构造的编码表，文件名为encoding_list_64to512_bg.csv
 + -i ASIC参与编码的顺序，在大面积探测器中，可能是两个ASIC读出探测器的一个维度，因此需要使用-i参数指定参与编码的ASIC的顺序
 + -s 噪声的 sigma 值，用于计算阈值，默认为 3
 + -d 有效数据的起始SCA位置，默认为 0
 + -D 有效数据的结束SCA位置，默认为 255
 + -T 输出文本文件，最初用于条使用，后续将去掉
 + -o 输出文件名后缀，用于区分是否经过解码，默认为 _dec
 + -m 最大的rms值，在有些情况下，可能计算得到的rms值特别大，用这个来限制一下。这种情况一般是系统没有接好，后续这个功能将取消。
 + -v 是否显示warning信息，默认为不显示
 + -C 解码的过程是否允许不连续点的存在，默认为允许，并输出warning

## 程序初始化
### 创建输出文件 
+ 最开始先创建输出文件，包括如下的内容
```c++
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
```
### 基本运行参数设置
#### 初始化解码函数
目前允许同时使用两个解码函数，以指针方式初始化
在当前的设计中

#### 读取ASIC芯片的顺序表
在前面提到的 -i 参数指定的芯片顺序中，通常情况下是使用两个芯片来读取一个探测器一个维度的信息。目前默认的芯片顺序表是 [group_lists_20220715_link0.csv](./group_lists_20220715_link0.csv)， 文件内容如下。包含两行，每一行表示一个编码的组合，每一行的最后一个数字表示使用的编码表的类型。0 表示 -e 参数指定的第一个编码表，1 表示 -E 参数指定的编码表。每行除最后一个数字外，其余的数字表示参与编码的芯片号，比如说第一行表示的是芯片号为 0 和 1 的芯片参与编码，使用的编码表为 -e 参数指定的第一个编码表；第二行表示的是芯片号为 2 和 3 的芯片参与编码，使用的编码表为 -e 参数指定的第一个编码表。

| 0    | 1    | 0    |
| ---- | ---- | ---- |
| 2    | 3    | 0    |

再举一个例子，如下的一个芯片分组顺序表表示的是第 0 和 1 号芯片使用 -e 参数指定的第一个编码表，第 3 号芯片使用 -E 参数指定的第二个编码表，而第 2 号芯片不参与编码。
| 0    | 1    | 0    |
| ---- | ---- | ---- |
| 3    | 1 |


## 程序运行
### 数据比较
+ 从每个板获取到的数据会被送到board_data_cmp()中进行比较，这个函数将会返回一个vector<double>类型的数据。如果没有数据在指定的范围内过阈，那么返回一个长度为0的数组；如果在数据指定的位置过阈了，返回击中信号的最大值，时间，位置信息
+ 上述的比较结果将会存到single_event_hit_chn, single_event_hit_amp, single_event_hit_time, single_event_hit_max_position中，用于解码程序处理
+ 若一次事例中有通道信号过阈，那么就会调用解码函数hit_data_dec()

### 解码
+ 在解码的过程中，先将要合到一起解码的通道存到一个vector<int> aget_chn中，注意，这些通道可能是分布在不同的ASIC中的。如现在采用的$400~\mathrm{mm}$探测器，一个维度的通道数为1000路，需要两个ASIC进行编码读出，因此在解码的过程中需要考虑两个ASIC的解码，尤其是击中的位置在两个ASIC的交界处的情况。
+ 接下来调用hit_dec->run_dec()，文件在[hit_dec.cc](..\..\my_pub\src\hit_dec.cc)中
+ 程序会进行一些击中通道的信息检查，其中一项就是输入的通道数数目，如果超过了设定的max_hit_chn， 那么就会输出warning信息，但是程序会继续运行，这个参数的设定在[hit_dec.hh](..\..\my_pub\lib\hit_dec.hh)中，其余的检查不再赘述
+ 之后将调用hit_dec_total函数，这个函数采用遍历的办法
  + 首先先将所有的击中通道可能对应的探测器通道存到一个数组中
  + 然后按照探测器通道号从小到大的顺序进行排列
  + 由于在探测器侧的击中通道号是连续的，因此直接筛选出连续的通道即可
    + 筛选过程中需要同时将探测器通道对应的电子学通道和相应的赋值都要存下来
    + 之后计算每个连续点中间有多少个 hole，这个hole的定义是：如果两个相邻的电子学通道号相差大于1，那么就认为中间有一个hole
+ 然后在主程序中，将解码得到的信息存到相应的变量中，最后将这些变量存到TTree中
