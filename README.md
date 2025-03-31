清华400 缪子初始数据转换

语言
  cmake
  c++

依赖
  CERN-ROOT
  Boost（选用）
  c++ 标准库

本工程分两步：
   解包: unpack
   时间排序(电子学逻辑导致): time_sort
                              
make:
  mkdir your/path/of/build && cd your/path/of/build && cmake path/of/this_project && cmake --build , -j4
  make出两个可执行文件
  [cmake /path/of/this/project -DBUILD_CALC_AND_DRAW=ON] enable calcu_and_draw function 

用法
  解包:两个参数： unpack your/path/data.dat 1
    p0: 初始数据文件
    p1: 模式(0: 基线 1：数据)
  时间排序: time_sort your/path/of/unpack.root your/path/of/baseline.txt

使用例：
  ./unpack ../rawdat/20241220162844_baseline.dat 0
  ./unpack ../rawdat/20241220162911_trigger.dat 1
  ./time_sort ../rawdat/20241220162911_trigger_entry.root calc/20241220162844_baseline_prestal.txt
  
测试数据:
  测试数据位于rawdat中

ana:
  ref aget codes
