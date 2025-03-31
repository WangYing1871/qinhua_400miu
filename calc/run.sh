#!/bin/bash

usb_path="/home/wangying/work/400/400moun_system/data"
#echo ${usb_path}
rawdata_path=""
rawdata_filename="run_20240806"
dec_out="dec_out"
data_dec="data_dec"
append="out"
noise_file="noise_20240806-0_out"
encode_csv_path_8="./encoding_list_64to512_real.csv"
encode_csv_path_16="./encoding_list_64to1024_real.csv"
group_csv_path_3="./group_lists_jw_link3.csv"
group_csv_path="./group_lists_jw.csv"
position_csv_path="./link_position_jw_6layer.csv"
zposition_csv_path="./z_position_jw_6layer.csv"

#cd ./usb_data_separate/build/
#for i in {0..10}
#do 
#	./usb_data_separate -f /home/wangying/work/400/400moun_system/data/run_20240806-${i}.dat -O /home/wangying/work/400/400moun_system/data/dec_out/
#done
#
#cd ../../
#
#echo ${usb_path}

#exit 0

#cd ./fec_data_dec/build/
#for file in {0..10}
#do
#	for num in {0..5}
#	do 
##echo "./fec_data_dec -f ${usb_path}/${rawdata_path}/${dec_out}/${rawdata_filename}-${file}_mt-${num}.dat"
#    ./fec_data_dec -f ${usb_path}/${rawdata_path}/${dec_out}/${rawdata_filename}-${file}_mt-${num}.dat 
#	done
#done


#cd ../../

#cd ./fec_noise_calc/build/
#for num in {0..5}
#do
#	./fec_noise_calc -f ${usb_path}/${rawdata_path}/${dec_out}/${noise_file}_mt-${num}.root -s 10 -o
#done
#
#exit 0

#cd ../../

cd ./fec2det_dec/build/
for file in {0..0}
do
	for num in {0..1}
	do
    ./fec2det_dec -e ${encode_csv_path_8} -i ${group_csv_path} -F ${usb_path}/${rawdata_path}/${dec_out}/result_${noise_file}_mt-${num}.root -f ${usb_path}/${rawdata_path}/${dec_out}/${rawdata_filename}-${file}_out_mt-${num}.root -d 200 -D 300 -C -s 10
	done
done

cd ../../

cd ./hit_position_reconstrcut/build/
for file in {0..20}
do
	for num in {0..6}
	do
	    	./hit_position_reconstrcut -p ${position_csv_path} -z ${zposition_csv_path} -L 0x3F -l 6 -f ${usb_path}/${rawdata_path}/${dec_out}/${data_dec}/${rawdata_filename}-${file}_out_dec_con_10sigma_mt-${num}.root -V
	done
done

#cd ../../

cd ./tracker_seek/build/
for file in {0..20}
do
	./tracker_seek -f ${usb_path}/${rawdata_path}/${dec_out}/${data_dec}/${rawdata_filename}-${file}_out_dec_con_10sigma_mt_hit_data_link0x3F.root -L 5 -D 6 -l 0x3F -r ../offset_test.txt -j -T
done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#hadd 合并所有tracker_seek文件，作为画图的输入
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#画图,和计算效率与位置分辨率

#calc_and_draw()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#cd ../../

#cd ./adas_muxbrd_fit/build/
#for file in {0..119}
#do 
#	./adas_muxbrd_sync -f ${usb_path}/${rawdata_path}/${dec_out}/${data_dec}/${rawdata_filename}-${file}_out_5sigma_discon_5sigma_mt_hit_data_link0x3F.root -O result
#done
