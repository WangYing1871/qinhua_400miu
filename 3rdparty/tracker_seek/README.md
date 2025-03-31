# best_tracker
这个程序主要用来给 [hit_position_reconstruction](../hit_position_reconstruction) 处理后的所有可能的击中位置
程序会遍历每一层探测器的所有可能的击中位置，然后进行拟合，根据拟合的 RMSE 以及是否满足用户设定的条件选择最优的拟合结果存在 root 文件和 txt 文件中

offset_initial.txt  每行对应一个探测器，每行的前三列对应xyz方向的坐标，后三列对应xyz方向的旋转

## 程序输入
需要设置
