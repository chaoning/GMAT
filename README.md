# 使用说明-中文版

## 1. 联系方式
ningchao(at)sdau(dot)edu(dot)cn

## 2. 安装
所有模块均已在64位linux系统下编译完成,下载并修改为可执行权限即可使用。  
修改文件权限命令:chmod 777 文件

## 3.功能模块

### 3.1 关系矩阵构建 
#### （1）参数详解: 

```
 -b, --bfile [FILE]                                      Plink Binary PED files前缀
 -d, --dfile [FILE]                                      Dosage基因型文件全名
 -a, --fam [FILE]                                        PLINK格式的样本信息文件全名
 -o, --out [FILE]                                        输出文件前缀
 -g, --grm [agrm, dgrm_as, dgrm_gs, aagrm, 
            adgrm_as, adgrm_gs, ddgrm_as, 
            ddgrm_gs, agrm_dosage]                       基因组关系矩阵类型
            agrm:                                        加性基因组关系矩阵
            dgrm_as:                                     显性基因组关系矩阵,生物学意义显性矩阵,可用于关联分析
            dgrm_gs:                                     显性基因组关系矩阵,育种学意义显性矩阵,用于显性效应剖分、基因组选择
            aagrm:                                       加加上位基因组关系矩阵
            adgrm_as:                                    加显上位基因组关系矩阵,其中显性矩阵是生物学意义显性矩阵
            adgrm_gs:                                    加显上位基因组关系矩阵,其中显性矩阵是育种学意义显性矩阵
            ddgrm_as:                                    显显上位基因组关系矩阵,其中显性矩阵是生物学意义显性矩阵
            ddgrm_gs:                                    显显上位基因组关系矩阵,其中显性矩阵是育种学意义显性矩阵   
            agrm_dosage:                                 加性基因组关系矩阵,读取Dosage基因型文件计算   
 -n, --npart [num, default is 20]                        将全基因组标记等分,分部分读取计算,节省内存
 -f, --fmt [num, default is 2]                           关系矩阵输出格式:0,“矩阵”格式;1,“行号、列号、值”格式;2,“个体号、个体号、值”格式
 -i, --inv [No argument]                                 求逆矩阵
 -m, --mkl_num_threads [num, default is max]             MKL库线程数，默认为当前服务器最大线程数。
 -v, --val [float, Default is 0.001]                     加到矩阵对角线的值,保证矩阵正定
 -s, --skipcols [int, Default is 1]                      读取文件时跳过的列数,仅对Dosage基因型文件(--dfile)有效
 -h, --help                                              打印此帮助信息

```
#### （2）输入文件：  

Plink Binary PED files (*.bed, *.bim, *.fam)：参考PLINK软件说明书  

Dosage基因型文件：第一行为文件头，后面每行表示一个位点的信息，前几列是位点位置信息，后面是个体的基因型编码值（0、1、2）或者基因型剂量（0-2之间的连续值，此外可以是基因表达量、微生物丰度等。读取时，需设置--skipcols参数，跳过位点信息。示例如下，前6列为位点位置信息。

| CHR  | SNP                | (C)M   | POS  | COUNTED | ALT  | 12659_14462 | 12659_14463 | 12659_14464 | 12659_14465 | 12659_14466 | 12659_14467 | 12659_14468 |
| ---- | :----------------- | ------ | ---- | ------- | ---- | ----------- | ----------- | ----------- | ----------- | ----------- | ----------- | ----------- |
| 1    | mCV25266528        | 10.977 | 0    | A       | T    | 1           | 1           | 1           | 0           | 1           | 1           | 1           |
| 1    | gnf01.037.906      | 20.135 | 0    | A       | T    | 1           | 2           | 1           | 0           | 1           | 1           | 1           |
| 1    | petAF067836-350A-1 | 29.195 | 0    | T       | A    | 1           | 1           | 1           | 1           | 1           | 1           | 1           |
| 1    | mCV23057534        | 38.008 | 0    | A       | T    | 2           | 2           | 1           | 1           | 1           | 1           | 2           |
| 1    | rs13475932         | 39.214 | 0    | A       | T    | 2           | 2           | 2           | 2           | 2           | 2           | 2           |

PLINK格式的样本信息文件：即Plink Binary PED files 文件中的*.fam文件，Dosage基因型文件出现时需要提供此文件，用以提供个体信息。

#### （3）输出文件：

关系矩阵文件：*\.\<grm\>\.\<fmt\>_fmt，依据给出--fmt参数的不同，输出不同格式的基因组关系矩阵，--fmt 0输出矩阵格式，--fmt 1输出行号、列号、值形式，--fmt 2输出个体号、个体号、值形式。

关系矩阵的逆矩阵文件：*.\<giv\>.\<fmt\>_fmt，依据给出--fmt参数的不同，输出不同格式的基因组关系矩阵的逆矩阵。

个体号文件：*.id，与.fam文件的第二列相同

#### （4）加性基因组关系矩阵

读取Plink Binary PED files (plink.bed, plink.bim, plink.fam),根据需求,可生成不同输出格式的加性基因组关系矩阵及逆矩阵。 
参考文献:VanRaden PM. Efficient methods to compute genomic predictions. J Dairy Sci. 2008 Nov;91(11):4414-23. 

```
#示例数据在mouse文件夹下
# 1. 默认参数输出
gmatrix --bfile plink --grm agrm --out test 
# 输出结果包含两个文件: test.id文件包含个体id号;test.agrm.id_fmt包含加性基因组关系矩阵信息, 为三列:个体号、个体号、值。

# 2. 输出格式为行号、列号、值
gmatrix --bfile plink --grm agrm --fmt 1 --out test 
# 输出结果包含两个文件: test.id文件包含个体id号;test.agrm.ind_fmt包含加性基因组关系矩阵信息, 为三列:行号、列号、值。

# 3. 输出格式为矩阵形式
gmatrix --bfile plink --grm agrm --fmt 0 --out test 
# 输出结果包含两个文件: test.id文件包含个体id号;test.agrm.mat_fmt包含加性基因组关系矩阵信息, 矩阵形式。

# 4. 输出逆矩阵
gmatrix --bfile plink --grm agrm --inv --out test
# 输出结果包含三个文件：test.id文件包含个体id号；test.agrm.id_fmt包含加性基因组关系矩阵信息, 为三列:个体号、个体号、值；test.agiv.id_fmt包含加性基因组关系矩阵的逆矩阵信息, 为三列:个体号、个体号、值。

# 5. 将标记分部分读取、依次计算以节省内存(增加--npart参数值,默认为20)
gmatrix --bfile plink --grm agrm --npart 30 --out test 

# 6. 设置mkl线程数(默认为当前服务器最大线程数)
gmatrix --bfile plink --grm agrm --mkl_num_threads 10 --out test 

# 7. 矩阵对角线上加值,保证正定(默认为0.001)
gmatrix --bfile plink --grm agrm --val 0.01 --out test 

# 8. 读取Dosage基因型文件,计算加性基因组关系矩阵
gmatrix --dfile plink.traw --fam plink.fam --skipcols 6 --grm agrm_dosage --out test 

```
#### （5）显性基因组关系矩阵
参考文献:Zulma G Vitezica, Luis Varona, Andres Legarra, On the Additive and Dominant Variance and Covariance of Individuals Within the Genomic Selection Scope, Genetics, Volume 195, Issue 4, 1 December 2013, Pages 1223–1230. 
```
#示例数据在mouse文件夹下，运行方式与（4）加性基因组关系矩阵构建类似
# 1. 构建生物学意义的显性基因组关系矩阵
gmatrix --bfile plink --grm dgrm_as --out test

# 2. 构建育种学意义的显性基因组关系矩阵
gmatrix --bfile plink --grm dgrm_gs --out test

```
#### （6）加加上位基因组关系矩阵
依据加性和加性基因组关系矩阵直积计算得到
```
#示例数据在mouse文件夹下，运行方式与（4）加性基因组关系矩阵构建类似
gmatrix --bfile plink --grm aagrm --out test
```
#### （7）加显上位基因组关系矩阵
依据加性和显性基因组关系矩阵直积计算得到
```
#示例数据在mouse文件夹下，运行方式与（4）加性基因组关系矩阵构建类似
# 1. 显性矩阵为生物学意义的关系矩阵
gmatrix --bfile plink --grm adgrm_as --out test
# 2. 显性矩阵为育种学意义的关系矩阵
gmatrix --bfile plink --grm adgrm_gs --out test

```
#### （8）显显上位基因组关系矩阵

依据显性和显性基因组关系矩阵直积计算得到
```
#示例数据在mouse文件夹下，运行方式与（4）加性基因组关系矩阵构建类似
# 1. 显性矩阵为生物学意义的关系矩阵
gmatrix --bfile plink --grm ddgrm_as --out test
# 2. 显性矩阵为育种学意义的关系矩阵
gmatrix --bfile plink --grm ddgrm_gs --out test
```
### 3.2 单性状全基因组关联分析

#### （1）参数详解: 

```
 -d, --data <FILE>                         数据文件,包含表头(列名),第一列必须是个体号,缺失值用NA表示
 -b, --bfile <PLINK>                       Plink Binary PED files 前缀
 -o, --out <FILE>                          输出文件前缀
 -g, --grm <FILE>                          加性基因组关系矩阵前缀
 -t, --trait <str>                         待分析性状在数据文件中的列名
 --repeat                                  使用重复力模型(有重复测量值的性状,例如,猪产仔数)
 --covar <str1,str2,...,str(n)>            放入模型中作为固定效应的连续变量列名
 --class <str1,str2,...,str(n)>            放入模型中作为固定效应的分类变量列名
 --snp_set <start_snp,end_snp>             起始和终止SNP名字,对此区段的SNP进行关联分析
 --continue                                估计方差组分时利用上一轮方差组分结果继续迭代
 --condition <SNP1,SNP2,...,SNP(n)>        加入到固定效应中的SNP名字
 --maxiter <200>                           不加SNP标记时估计方差组分的最大迭代次数
 --maxiter0 <30>                           对每个标记逐个检验时的方差组分迭代次数
 --cc_par <1.0e-8>                         方差组分收敛标准,方差组分变化差值的平方和的平方根
 --cc_gra <1.0e-6>                         方差组分收敛标准,一阶偏导的平方和的平方根
 --npart <20>                              将全基因组标记等分,分部分进行关联分析,以节省内存。
 --exact                                   利用精确的方法进行全基因组关联分析,每个标记均进行方差组分估计。
 --p_cut <1.0e-3>                          利用近似方式进行关联分析时,低于此P值的SNP利用精确方法重新计算P值。
 --parallel <ith,total>                    全基因组标记等分为total份,仅对第i份进行分析。
 --mkl_num_threads <max>                   MKL线程数，默认为当前服务器最大线程数。
 -h, --help                                打印此帮助信息并退出。
```
#### （2）输入文件

数据文件：包含文件头，第一列必须是个体号，后面是固定效应和表型，缺失值用***NA***表示（必须是大写的***NA***）。

| id    | mean | sex  | age  | treat | trait |
| ----- | ---- | ---- | ---- | ----- | ----- |
| 14462 | 1    | 0    | 126  | 0     | 0.58  |
| 14463 | 1    | 0    | 91   | 1     | 0.39  |
| 14464 | 1    | 1    | 126  | 0     | 0.37  |
| 14465 | 1    | 0    | 91   | 1     | NA    |
| 14466 | 1    | 0    | 91   | 1     | 0.84  |
| 14467 | 1    | 0    | 91   | 1     | 0.61  |
| 14468 | 1    | 1    | 91   | 1     | NA    |
| 14469 | 1    | 1    | 91   | 1     | 0.59  |
| 14470 | 1    | 1    | 124  | 0     | 0.13  |

加性基因组关系矩阵文件：gmatrix模块构建的加性基因组关系矩阵，包含关系矩阵和个体号文件，关系矩阵为个体号、个体号、值格式。可调用gmatrix --bfile \<PLINK\> --grm agrm --out \<FILE\> 实现。

Plink Binary PED files (*.bed, *.bim, *.fam)：参考PLINK软件说明书 

#### （3）输出文件：

主要结果文件：prefix.起始标记位置_终止标记位置.res；共13列，依次是标记顺序（order）（从0开始）、染色体号（chro）、标记号（snp_id）、遗传距离（cm）、物理距离（base）、次要等位基因（allele0）、主要等位基因（allele1）、主要等位基因频率（freq）、方差组分变化差值的平方和的平方根（cc_par)、一阶偏导的平方和的平方根（cc_gra)、标记效应值（eff）、标记效应值的标准误（se）、检验p值（p）。  

cc_par=1000.0、cc_gra=1000.0表示该标记没有重新估计方差组分，而是利用null model得到方差组分进行P值计算。

| order | chro | snp_id             | cm     | base | allele0 | allele1 | freq     | cc_par   | cc_gra   | eff      | se       | p        |
| ----- | ---- | ------------------ | ------ | ---- | ------- | ------- | -------- | -------- | -------- | -------- | -------- | -------- |
| 0     | 1    | mCV25266528        | 10.977 | 0    | A       | T       | 0.549463 | 2.02E-10 | 3.15E-07 | 0.038012 | 0.03837  | 0.321844 |
| 1     | 1    | gnf01.037.906      | 20.135 | 0    | A       | T       | 0.530675 | 1.20E-10 | 1.87E-07 | -0.01398 | 0.030733 | 0.649209 |
| 2     | 1    | petAF067836-350A-1 | 29.195 | 0    | T       | A       | 0.50115  | 4.62E-10 | 7.17E-07 | 0.015869 | 0.032372 | 0.623985 |
| 3     | 1    | mCV23057534        | 38.008 | 0    | A       | T       | 0.57362  | 5.13E-10 | 7.99E-07 | -0.04683 | 0.042799 | 0.27387  |
| 4     | 1    | rs13475932         | 39.214 | 0    | A       | T       | 0.555598 | 7.69E-11 | 1.22E-07 | 0.029993 | 0.046017 | 0.51455  |
| 5     | 1    | gnf01.075.385      | 39.786 | 0    | A       | T       | 0.554064 | 7.97E-11 | 1.27E-07 | 0.041027 | 0.045727 | 0.369612 |
| 6     | 1    | mCV23431007        | 41.28  | 0    | A       | T       | 0.59931  | 9.46E-11 | 1.47E-07 | -0.01466 | 0.044237 | 0.740378 |
| 7     | 1    | mCV23433457        | 41.735 | 0    | A       | T       | 0.617331 | 5.54E-10 | 8.55E-07 | -0.03208 | 0.043901 | 0.464992 |
| 8     | 1    | UT-1-92.862916     | 45.802 | 0    | A       | T       | 0.562883 | 4.08E-10 | 6.40E-07 | 0.032283 | 0.038055 | 0.396257 |
| 9     | 1    | CEL-1-111503693    | 50.839 | 0    | T       | A       | 0.512653 | 1.09E-10 | 1.71E-07 | 0.046876 | 0.054495 | 0.38968  |
| 10    | 1    | mCV22824651        | 50.839 | 0    | T       | A       | 0.509586 | 2.10E-10 | 3.27E-07 | 0.031347 | 0.054163 | 0.562755 |

方差组分文件：prefix.var；包含文件头，共四列；第一列是随机效应编码，第二列是随机效应方差协方差的行号、第三列是随机效应方差协方差的列号、第四列是方差组分值。单性状分析不涉及方差协方差结构，第二行为加性遗传方差、第三行为残差方差。

| random | cov_row | cov_col | value    |
| ------ | ------- | ------- | -------- |
| 1      | 1       | 1       | 0.087966 |
| 2      | 1       | 1       | 0.14003  |

方差组分迭代文件：prefix.ai_mat、prefix.em_mat、prefix.wemai_mat、prefix.fd_mat分别是最后一次迭代的AI矩阵、EM矩阵、AI和EM加权阵、一阶偏导向量。

#### （4）示例代码

小鼠数据：mouse文件夹下

```
# 计算加性基因组关系矩阵
gmatrix --bfile plink --grm agrm --out test

# 仅估计方差组分，sex和treat作为分类变量放入固定效应，age作为协变量放入固定效应
uvlmm --data pheno --grm test --trait trait --class sex,treat --covar age --out test

# 全基因组关联分析
uvlmm --bfile plink --data pheno --grm test --trait trait --class sex,treat --covar age --out test
# test.0_1407.res为主要输出结果文件

# 将标记等分为3份进行全基因组关联分析，并行提交，缩短运算时间
uvlmm --bfile plink --data pheno --grm test --trait trait --class sex,treat --covar age --out test --parallel 1,3
uvlmm --bfile plink --data pheno --grm test --trait trait --class sex,treat --covar age --out test --parallel 2,3
uvlmm --bfile plink --data pheno --grm test --trait trait --class sex,treat --covar age --out test --parallel 3,3
# test.0_469.res、test.469_938.res和test.938_1407.res为主要输出结果文件，合并后为全基因组分析结果

# 对一区段内标记进行关联分析，mCV25266528为起始标记，rs3663706为终止标记
uvlmm --bfile plink --data pheno --grm test --trait trait --class sex,treat --covar age --out test --snp_set mCV25266528,rs3663706

# 加入Top SNP进行条件GWAS分析，可放入多个Top SNP，逗号隔开
uvlmm --bfile plink --data pheno --grm test --trait trait --class sex,treat --covar age --out test_condition --condition rs13477831

# 利用准确但相对耗时的方法进行全基因组关联分析，即全基因组标记均进行方差组分估计
uvlmm --exact --bfile plink --data pheno --grm test --trait trait --class sex,treat --covar age --out test 

# 下列代码同样实现上述功能，即p值低于1的标记均重新利用精确算法重新计算p值，不推荐
uvlmm --bfile plink --data pheno --grm test --trait trait --class sex,treat --covar age --out test --p_cut 1
```
酵母数据：yeast文件夹下，有重复测量值，利用重复力模型进行全基因组关联分析

```
# 计算加性基因组关系矩阵
gmatrix --bfile yeast --grm agrm --out test

# 仅估计方差组分，batch为分类变量、无协变量
uvlmm --repeat --data pheno --grm test --trait trait --class batch --out test

# 全基因组关联分析
uvlmm --repeat --bfile yeast --data pheno --grm test --trait trait --class batch --out test
```

### 3.3多关系矩阵模块

此模块可将多个基因组关系矩阵同时放入模块，剖分不同随机效应的方差组分，用以解析不同随机效应占表型方差的比例，输入和输出文件格式与uvlmm估计方差组分结果相同。此模块与remma仅估计方差组分的功能相同

#### （1）参数详解

```
 -d, --data <FILE>                         Plink Binary PED files前缀
 -o, --out <FILE>                          输出文件前缀
 -g, --grm <FILE1,FILE2,...,FILE(n)>       多个基因组关系矩阵全名
 -t, --trait <str>                         待分析性状在数据文件中的列名
 --repeat                                  使用重复力模型(有重复测量值的性状,例如,猪产仔数)
 --covar <str1,str2,...,str(n)>            放入模型中作为固定效应的连续变量列名
 --class <str1,str2,...,str(n)>            放入模型中作为固定效应的分类变量列名
 --continue                                估计方差组分时利用上一轮方差组分结果继续迭代
 --maxiter <200>                           最大迭代次数
 --cc_par <1.0e-8>                         方差组分收敛标准,方差组分变化差值的平方和的平方根
 --cc_gra <1.0e-6>                         方差组分收敛标准,一阶偏导的平方和的平方根
 --mkl_num_threads <max>                   MKL线程数，默认为当前服务器最大线程数。
 -h, --help                                打印此帮助信息并退出。
```

#### （2）示例代码

```
# 小鼠数据：mouse文件夹下，解析不同类型遗传效应的方差组分
# 构建关系矩阵
gmatrix --bfile plink --grm agrm --out test  #加性
gmatrix --bfile plink --grm dgrm_gs --out test #显性
gmatrix --bfile plink --grm aagrm --out test #加加互作
gmatrix --bfile plink --grm adgrm_gs --out test #加显互作
gmatrix --bfile plink --grm ddgrm_gs --out test #显显互作

# 估计方差组分
multigrm --data pheno --grm test.agrm.id_fmt,test.dgrm_as.id_fmt,test.aagrm.id_fmt,test.adgrm_as.id_fmt,test.ddgrm_as.id_fmt --class sex,treat --covar age --trait trait --out test
```



### 3.4多性状全基因组关联分析

#### （1）参数详解

```
 --balance                               当成平衡数据分析，同一个体不同性状的测量值有缺失，将删除此个体.
 --unbalance                             当成非平衡数据分析，允许同一个体不同性状的测量值有缺失.
 -d, --data <FILE>                       数据文件,包含表头(列名),第一列必须是个体号,缺失值用NA表示
 -b, --bfile <FILE>                      Plink Binary PED files 前缀
 -o, --out <FILE>                        输出文件前缀
 -g, --grm <FILE>                        加性基因组关系矩阵前缀.
 -t, --trait <str1,str2,...,str(n)>      待分析性状的列名
 --covar <str1,str2,...,str(n)>          放入固定效应的连续变量列名
 --class <str1,str2,...,str(n)>          放入固定效应的分类变量列名
 --snp_set <start_snp,end_snp>           起始和终止SNP名字,对此区段的SNP进行关联分析
 --continue                              估计方差组分时利用上一轮方差组分结果继续迭代
 --condition <SNP1,SNP2,...,SNP(n)>      加入固定效应中的条件SNP名字
 --maxiter <200>                         不加SNP标记时估计方差组分的最大迭代次数
 --maxiter0 <30>                         对每个标记逐个检验时的方差组分迭代次数
 --cc_par <1.0e-8>                       方差组分收敛标准,方差组分变化差值的平方和的平方根
 --cc_gra <1.0e-6>                       方差组分收敛标准,一阶偏导的平方和的平方根
 -n --npart <20>                         将全基因组标记等分,分部分进行关联分析,以节省内存
 --exact                                 利用精确的方法进行全基因组关联分析,每个标记均进行方差组分估计
 --p_cut <1.0e-3>                        利用近似方式进行关联分析时,低于此P值的SNP利用精确方法重新计算P值
 --parallel <ith,total>                  全基因组标记等分为total份,仅对第i份进行分析
 --mkl_num_threads <max>                 MKL线程数，默认为当前服务器最大线程数
 -h, --help                              打印此帮助信息
```

#### （2）输入文件

数据文件：包含文件头，第一列必须是个体号，后面是固定效应和多个性状表型值，缺失值用***NA***表示（必须是大写的***NA***）。

| ID    | trait1   | trait2   |
| ----- | -------- | -------- |
| 01_01 | -2.25383 | NA       |
| 01_02 | -1.88775 | NA       |
| 01_03 | 2.415068 | -0.3207  |
| 01_04 | 2.417437 | NA       |
| 01_06 | -0.63451 | -1.44897 |
| 01_07 | 1.923285 | 1.545476 |
| 01_09 | 0.940961 | NA       |
| 01_10 | -0.05571 | 0.26933  |

加性基因组关系矩阵文件：gmatrix模块构建的加性基因组关系矩阵，包含关系矩阵和个体号文件，关系矩阵为个体号、个体号、值格式。可调用gmatrix --bfile \<PLINK\> --grm agrm --out \<FILE\> 实现。

Plink Binary PED files (\*.bed, \*.bim, \*.fam)：参考PLINK软件说明书 

#### （3）输出文件

主要结果文件：prefix.起始标记位置\_终止标记位置.res；依据性状数（n）不同而不同，依次是标记顺序（order）（从0开始）、染色体号（chro）、标记号（snp\_id）、遗传距离（cm）、物理距离（base）、次要等位基因（allele0）、主要等位基因（allele1）、主要等位基因频率（freq）、方差组分变化差值的平方和的平方根（cc\_par)、一阶偏导的平方和的平方根（cc\_gra)、标记效应值[eff\_(\*)，共n列]、标记效应估计值的方差或协方差[eff_var_(\*）\_(\*)，共(n+1)\*n/2列]、检验p值（p）。  

| order | chro | snp_id  | cm   | base  | allele0 | allele1 | freq     | cc_par | cc_gra | eff_0    | eff_1    | eff_var_0_0 | eff_var_1_0 | eff_var_1_1 | p        |
| ----- | ---- | ------- | ---- | ----- | ------- | ------- | -------- | ------ | ------ | -------- | -------- | ----------- | ----------- | ----------- | -------- |
| 67    | 10   | 5866308 | 0    | 33847 | T       | G       | 0.500662 | 1000   | 1000   | -0.0572  | -0.03159 | 0.001968    | 0.001767    | 0.001922    | 0.243204 |
| 68    | 10   | 5866473 | 0    | 34012 | A       | G       | 0.501103 | 1000   | 1000   | -0.06544 | -0.03818 | 0.001969    | 0.001768    | 0.001923    | 0.179184 |
| 69    | 10   | 5866608 | 0    | 34147 | A       | G       | 0.501544 | 1000   | 1000   | -0.06831 | -0.04005 | 0.001969    | 0.001768    | 0.001923    | 0.155641 |
| 70    | 10   | 5866899 | 0    | 34438 | G       | A       | 0.502426 | 1000   | 1000   | -0.07561 | -0.04937 | 0.001972    | 0.001771    | 0.001926    | 0.140538 |
| 71    | 10   | 5866926 | 0    | 34465 | A       | G       | 0.502867 | 1000   | 1000   | -0.07919 | -0.05279 | 0.001952    | 0.001751    | 0.001906    | 0.122015 |
| 72    | 10   | 5867142 | 0    | 34681 | G       | A       | 0.502867 | 1000   | 1000   | -0.07949 | -0.05302 | 0.001949    | 0.001749    | 0.001904    | 0.12013  |
| 73    | 10   | 5867298 | 0    | 34837 | T       | C       | 0.502867 | 1000   | 1000   | -0.06886 | -0.04319 | 0.001947    | 0.001747    | 0.001902    | 0.176801 |
| 74    | 10   | 5867379 | 0    | 34918 | G       | A       | 0.502867 | 1000   | 1000   | -0.06425 | -0.03881 | 0.00195     | 0.001749    | 0.001905    | 0.204452 |
| 75    | 10   | 5867496 | 0    | 35035 | A       | T       | 0.503308 | 1000   | 1000   | -0.06838 | -0.04386 | 0.001947    | 0.001746    | 0.001902    | 0.190853 |
| 76    | 10   | 5868024 | 0    | 35563 | T       | C       | 0.504191 | 1000   | 1000   | -0.08424 | -0.05801 | 0.001948    | 0.001747    | 0.001903    | 0.102217 |

方差组分文件：prefix.var；包含文件头，共四列；第一列是随机效应编码，第二列是随机效应方差协方差的行号、第三列是随机效应方差协方差的列号、第四列是方差组分值。random=1表示加性遗传方差协方差结构，random=2表示残差方差协方差结构。

| random | cov_row | cov_col | value    |
| ------ | ------- | ------- | -------- |
| 1      | 1       | 1       | 0.240966 |
| 1      | 2       | 1       | 0.237031 |
| 1      | 2       | 2       | 0.233331 |
| 2      | 1       | 1       | 0.558848 |
| 2      | 2       | 1       | 0.170067 |
| 2      | 2       | 2       | 0.578756 |

方差组分迭代文件：prefix.ai_mat、prefix.em_mat、prefix.wemai_mat、prefix.fd_mat分别是最后一次迭代的AI矩阵、EM矩阵、AI和EM加权阵、一阶偏导向量。

#### （4）平衡数据示例

同一个体的所有性状均有测量值，存在缺失的个体将被删掉

```
# 示例数据在mouse_long文件夹下
# 计算加性基因组关系矩阵
gmatrix --bfile plink --grm agrm --out test

# 仅估计方差组分，三性状模型，sex作为分类变量放入固定效应中
mvlmm --balance --data phe.balance.txt --grm test --trait trait1,trait2,trait3 --class sex --out test

# 全基因组关联分析
mvlmm --balance --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3 --class sex --out test
# test.0_11833.res为主要输出结果文件

# 将标记等分为3份进行全基因组关联分析，并行提交，缩短运算时间
mvlmm --balance --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3 --class sex --out test --parallel 1,3
mvlmm --balance --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3 --class sex --out test --parallel 2,3
mvlmm --balance --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3 --class sex --out test --parallel 3,3
# test.0_3944.res、test.3944_7888.res和test.7888_11833.res为主要输出结果文件，合并后为全基因组分析结果

# 对一区段内标记进行关联分析
mvlmm --balance --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3 --class sex --out test --snp_set UNC010515443,UNC172447

#加入Top SNP进行条件GWAS分析
mvlmm --balance --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3 --class sex --out test_condition --condition UNC18848064
```

#### （5）非平衡数据示例

同一个个体的多个性状可存在缺失测量值

```
# 示例数据在yeast文件夹下

# 计算加性基因组关系矩阵
gmatrix --bfile yeast --grm agrm --out test

# 仅估计方差组分，两性状模型
mvlmm --unbalance --data pheno2 --grm test --trait trait1,trait2 --out test

# 全基因组关联分析
mvlmm --unbalance --bfile yeast --data pheno2 --grm test --trait trait1,trait2 --out test
# test.0_28220.res为主要输出结果文件

# 将标记等分为3份进行全基因组关联分析，并行提交，缩短运算时间
mvlmm --unbalance --data pheno2 --grm test --trait trait1,trait2 --out test  # 估计方差组分
mvlmm --unbalance --bfile yeast --data pheno2 --grm test --trait trait1,trait2 --out test --continue --maxiter 0 --parallel 1,3
mvlmm --unbalance --bfile yeast --data pheno2 --grm test --trait trait1,trait2 --out test --continue --maxiter 0 --parallel 2,3
mvlmm --unbalance --bfile yeast --data pheno2 --grm test --trait trait1,trait2 --out test --continue --maxiter 0 --parallel 3,3
# 加入"--continue --maxiter 0"表示以null模型估计的方差组分为初值，不迭代，直接用于全基因组关联分析，避免多次估计方差组分耗时
# test.0_9406.res、test.9406_18812.res和test.18812_28220.res为主要结果文件，合并后为全基因组分析结果

# 对一区段内标记进行关联分析
mvlmm --unbalance --bfile yeast --data pheno2 --grm test --trait trait1,trait2 --out test --snp_set 7808236,8836653

# 加入Top SNP进行条件GWAS分析
mvlmm --unbalance --bfile yeast --data pheno2 --grm test --trait trait1,trait2 --out test_condition --condition 7895795
```

### 3.5纵向性状全基因组关联分析

Cite：

Chao Ning, Dan Wang, Lei Zhou, Julong Wei, Yuanxin Liu, Huimin Kang, Shengli Zhang, Xiang Zhou, Shizhong Xu and Jian-Feng Liu. Efficient Multivariate Analysis Algorithms for Longitudinal Genome-wide Association Studies. *Bioinformatics*, 2019, 35(23): 4879-4885.

Chao Ning†, Dan Wang†, Xianrui Zheng†, Qin Zhang, Shengli Zhang, Raphael Mrode, Jian-Feng Liu. Eigen decomposition expedites longitudinal genome-wide association studies for milk production traits in Chinese Holstein. *Genetics Selection Evolution*, 2018, 50(1): 12

Chao Ning†, Huimin Kang†, Lei Zhou, Dan Wang, Haifei Wang, Aiguo Wang, Jinluan Fu, Shengli Zhang, Jian-Feng Liu. Performance gains in genome-wide association studies for longitudinal traits via modeling time-varied effects. *Scientific Reports*, 2017, 7(1): 590.

#### （1）参数详解

```
 --balance                                 当成平衡数据分析，同一个体不同时间点的测量值有缺失，将删除此个体.
 --unbalance                               当成非平衡数据分析，允许同一个体不同时间点的测量值有缺失.
 -d, --data <FILE>                         数据文件,包含表头(列名),第一列必须是个体号,缺失值用NA表示
 -b, --bfile <FILE>                        Plink Binary PED files 前缀
 -o, --out <FILE>                          输出文件前缀
 -g, --grm <FILE>                          加性基因组关系矩阵前缀
 -t, --trait <str1,str2,...,str(n)>        待分析性状列名
 --tpoint <t1,t2,...,t(n) or str>          与分析性状对应的测量时间点（平衡数据）或者测量时间点的列名（非平衡数据）
 --tcovar <str1,str2,...,str(n)>           放入固定效应的随时间变化的协变量列名
 --tclass <str1,str2,...,str(n)>           放入固定效应的随时间变化的分类变量列名
 --covar <str1,str2,...,str(n)>            放入固定效应的不随时间变化的协变量列名
 --class <str1,str2,...,str(n)>            放入固定效应的不随时间变化的分类变量列名
 --snp_set <start_snp,end_snp>             起始和终止SNP名字,对此区段的SNP进行关联分析
 --order <3,3,3>                           固定效应、加性遗传效应和永久环境效应的勒让德多项式阶数，平衡数据时三者必须相等
 --continue                                估计方差组分时利用上一轮方差组分结果继续迭代
 --condition <SNP1,SNP2,...,SNP(n)>        加入固定效应中的条件SNP名字
 --maxiter <200>                           不加SNP标记时估计方差组分的最大迭代次数
 --maxiter0 <30>                           对每个标记逐个检验时的方差组分迭代次数
 --cc_par <1.0e-8>                         方差组分收敛标准,方差组分变化差值的平方和的平方根
 --cc_gra <1.0e-6>                         方差组分收敛标准,一阶偏导的平方和的平方根
 -n --npart <20>                           将全基因组标记等分,分部分进行关联分析,以节省内存
 --exact                                   利用精确的方法进行全基因组关联分析,每个标记均进行方差组分估计
 --p_cut <1.0e-3>                          利用近似方式进行关联分析时,低于此P值的SNP利用精确方法重新计算P值
 --parallel <ith,total>                    全基因组标记等分为total份,仅对第i份进行分析
 --mkl_num_threads <max>                   MKL线程数
 -h, --help                                打印此帮助信息并退出
```

#### （2）输入文件

数据文件：包含文件头，第一列必须是个体号，缺失值用***NA***表示（必须是大写的***NA***）。   

对于平衡数据，每个个体一行，顺次放入不同时间点的测量值，--tpoint提供的测量时间点应与性状出现的顺序相同，示例如下：

| ID   | sex  | trait1 | trait2 | trait3 | trait4 | trait5 | trait6 | trait7 | trait8 | trait9 | trait10 | trait11 | trait12 | trait13 | trait14 | trait15 | trait16 |
| ---- | ---- | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------- | ------- | ------- | ------- | ------- | ------- | ------- |
| 1419 | 1    | 4.9    | 8.4    | 10.5   | 11.9   | 13.8   | 15.7   | 16.7   | 18.2   | 17.5   | 18.1    | 18.4    | 18.7    | 19.1    | 19.5    | 19.5    | 21.4    |
| 1422 | 2    | 4.6    | 8.3    | 10.9   | 14     | 16.3   | 17.2   | 18     | 21.3   | 19.8   | 20.6    | 21.7    | 20.4    | 22.2    | 22.8    | 23.2    | 22.1    |
| 1433 | 2    | 8      | 9.7    | 11     | 13.7   | 17.3   | 17.5   | 17.7   | 18.4   | 19.1   | 19.4    | 19.3    | 19.6    | 20.2    | 21.6    | 21.3    | 21.6    |
| 1441 | 2    | 5.1    | 9.7    | 13.1   | 15.7   | 18.5   | 19.1   | 19.9   | 19.3   | 20     | 20.8    | 20.2    | 21.7    | 22      | 21.7    | 22.2    | 23.4    |
| 1457 | 1    | 4.9    | 8.4    | 10.6   | 11.5   | 13.7   | 14.9   | 15.6   | 15.4   | 15.5   | 16.4    | 15.9    | 16.4    | 16.5    | 17.4    | 17.9    | 18.5    |
| 1464 | 1    | 6.2    | 10.1   | 13.7   | 15.8   | 18.2   | 18.1   | 18.4   | 20     | 20.5   | 21.6    | 20.7    | 20.2    | 21.6    | 22.1    | 23.3    | 23.6    |
| 1471 | 2    | 5.2    | 9.2    | 12.2   | 16.1   | 21.3   | 20.8   | 22     | 22.4   | 23.6   | 23.9    | 23.6    | 24.1    | 23.9    | 24.2    | 25.2    | 26      |
| 1479 | 2    | 4.3    | 7.6    | 9.1    | 10.9   | 13.3   | 14.1   | 15.8   | 15.3   | 16.2   | 17.4    | 17.3    | 16.9    | 18.2    | 18.7    | 19.4    | 19.4    |
| 1487 | 2    | 5      | 8.7    | 11.6   | 13.2   | 17.8   | 17.2   | 17.8   | 18.2   | 19.5   | 20.6    | 20.9    | 21.1    | 21.3    | 22.2    | 22.2    | 22.5    |
| 1495 | 1    | 5.7    | 10.7   | 13.2   | 15.7   | 17.7   | 19     | 20.2   | 19.6   | 20.8   | 20.9    | 20.3    | 20.4    | 21.3    | 21.9    | 21.9    | 21.9    |

对于非平衡数据，每个个体多行，测量时间点与性状测量值对应放入文件中，示例如下：

| ID   | Sex  | weak | trait |
| ---- | ---- | ---- | ----- |
| 1419 | 1    | 1    | 4.9   |
| 1419 | 1    | 2    | 8.4   |
| 1419 | 1    | 3    | 10.5  |
| 1419 | 1    | 4    | 11.9  |
| 1419 | 1    | 5    | 13.8  |
| 1419 | 1    | 6    | 15.7  |
| 1419 | 1    | 7    | 16.7  |
| 1419 | 1    | 8    | 18.2  |
| 1419 | 1    | 9    | 17.5  |
| 1419 | 1    | 10   | 18.1  |

加性基因组关系矩阵文件：gmatrix模块构建的加性基因组关系矩阵，包含关系矩阵和个体号文件，关系矩阵为个体号、个体号、值格式。可调用gmatrix --bfile <PLINK> --grm agrm --out <FILE> 实现。

Plink Binary PED files (\*.bed, \*.bim, \*.fam)：参考PLINK软件说明书 

#### （3）输出文件

主要结果文件：prefix.起始标记位置\_终止标记位置.res；依据拟合标记的勒让德多项式阶数不同（df）不同而不同，依次是标记顺序（order）（从0开始）、染色体号（chro）、标记号（snp\_id）、遗传距离（cm）、物理距离（base）、次要等位基因（allele0）、主要等位基因（allele1）、主要等位基因频率（freq）、方差组分变化差值的平方和的平方根（cc\_par)、一阶偏导的平方和的平方根（cc\_gra)、拟合标记效应不同阶数的多项式系数值[eff\_(\*)，共df+1列]、估计值的方差或协方差[eff_var_(\*）\_(\*)，共(df+2)\*(df+1)/2列]、检验p值（p）。  



#### （4）平衡数据

不同个体在所有时间点均有测量值,无缺失

```
# mouse_long 数据，1-16周，每周测量一次
# 计算加性基因组关系矩阵
gmatrix --bfile plink --grm agrm --out test

# 仅估计方差组分，sex为随时间变化的分类效应
longwas --balance --data phe.balance.txt --grm test --trait trait1,trait2,trait3,trait4,trait5,trait6,trait7,trait8,trait9,trait10,trait11,trait12,trait13,trait14,trait15,trait16 --tpoint 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 --tclass sex --out test

# 全基因组关联分析
longwas --balance --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3,trait4,trait5,trait6,trait7,trait8,trait9,trait10,trait11,trait12,trait13,trait14,trait15,trait16 --tpoint 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 --tclass sex --out test
# test.0_11833.res为主要输出结果文件


# 将标记等分为3份进行全基因组关联分析，并行提交，缩短运算时间
longwas --balance --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3,trait4,trait5,trait6,trait7,trait8,trait9,trait10,trait11,trait12,trait13,trait14,trait15,trait16 --tpoint 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 --tclass sex --out test --parallel 1,3
longwas --balance --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3,trait4,trait5,trait6,trait7,trait8,trait9,trait10,trait11,trait12,trait13,trait14,trait15,trait16 --tpoint 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 --tclass sex --out test --parallel 2,3
longwas --balance --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3,trait4,trait5,trait6,trait7,trait8,trait9,trait10,trait11,trait12,trait13,trait14,trait15,trait16 --tpoint 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 --tclass sex --out test --parallel 3,3
# test.0_3944.res、test.3944_7888.res和test.7888_11833.res为主要输出结果文件，合并后全基因组关联分析结果

# 对一区段内标记进行关联分析
longwas --balance --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3,trait4,trait5,trait6,trait7,trait8,trait9,trait10,trait11,trait12,trait13,trait14,trait15,trait16 --tpoint 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 --tclass sex --out test --snp_set UNC010515443,UNC172447

#加入Top SNP进行条件GWAS分析
longwas --balance --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3,trait4,trait5,trait6,trait7,trait8,trait9,trait10,trait11,trait12,trait13,trait14,trait15,trait16 --tpoint 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 --tclass sex --out test_condition --condition UNC18848064
```

#### （5）非平衡数据

同一个体在不同时间点可能存在缺失测量值

```
# mouse_long 数据

# 计算加性基因组关系矩阵
gmatrix --bfile plink --grm agrm --out test

# 仅估计方差组分，sex为随时间变化的分类效应
longwas --unbalance --data phe.unbalance.txt --grm test --trait trait --tpoint weak --tclass sex --out test_un

# 全基因组关联分析
mvlmm --unbalance --bfile yeast --data pheno2 --grm test --trait trait1,trait2 --out test
# test.0_28220.res为主要输出结果文件

# 将标记等分为3份进行全基因组关联分析，并行提交，缩短运算时间
mvlmm --unbalance --data pheno2 --grm test --trait trait1,trait2 --out test  # 估计方差组分
mvlmm --unbalance --bfile yeast --data pheno2 --grm test --trait trait1,trait2 --out test --continue --maxiter 0 --parallel 1,3
mvlmm --unbalance --bfile yeast --data pheno2 --grm test --trait trait1,trait2 --out test --continue --maxiter 0 --parallel 2,3
mvlmm --unbalance --bfile yeast --data pheno2 --grm test --trait trait1,trait2 --out test --continue --maxiter 0 --parallel 3,3
# 加入"--continue --maxiter 0"表示以null模型估计的方差组分为初值，不迭代，直接用于全基因组关联分析，避免多次估计方差组分耗时
# test.0_9406.res、test.9406_18812.res和test.18812_28220.res为主要结果文件，合并后为全基因组分析结果

# 对一区段内标记进行关联分析
mvlmm --unbalance --bfile yeast --data pheno2 --grm test --trait trait1,trait2 --out test --snp_set 7808236,8836653

# 加入Top SNP进行条件GWAS分析
mvlmm --unbalance --bfile yeast --data pheno2 --grm test --trait trait1,trait2 --out test_condition --condition 7895795
```

### 3.6非加性效应检验

Cite:

Dan Wang, Hui Tang, Jian-Feng Liu, Shizhong Xu, Qin Zhang* and Chao Ning*. Rapid Epistatic Mixed Model Association Studies by Controlling Multiple Polygenic Effects. *Bioinformatics*, 2020.

Chao Ning, Dan Wang, Huimin Kang, Raphael Mrode, Lei Zhou, Shizhong Xu, Jian-Feng Liu. A rapid epistatic mixed-model association analysis by linear retransformations of genomic estimated values. *Bioinformatics*, 2018, 34(11): 1817-1825.

#### （1）参数详解

```
 -d, --data <FILE>                       数据文件,包含表头(列名),第一列必须是个体号,缺失值用NA表示
 -b, --bfile <PLINK>                     Plink Binary PED files 前缀
 -o, --out <FILE>                        输出文件前缀
 -g, --grm <FILE1,FILE2,...,FILE(n)>     多个基因组关系矩阵全名
 -t, --trait <str>                       待分析性状在数据文件中的列名
 --epis <abya,abyd,dbyd,add,dom>         检验遗传效应类型
        abya                             加加上位
        abyd                             加显上位
        dbyd                             显显上位
        add                              加性
        dom                              显性
 --repeat                                重复力模型
 --covar <str1,str2,...,str(n)>          放入模型中作为固定效应的连续变量列名
 --class <str1,str2,...,str(n)>          放入模型中作为固定效应的分类变量列名
 --continue                              估计方差组分时利用上一轮方差组分结果继续迭代
 --condition <SNP1,SNP2,...,SNP(n)>      加入到固定效应中的SNP名字
 --maxiter <1000>                        不加SNP标记时估计方差组分的最大迭代次数
 --cc_par <1.0e-8>                       方差组分收敛标准,方差组分变化差值的平方和的平方根
 --cc_gra <1.0e-6>                       方差组分收敛标准,一阶偏导的平方和的平方根
 --npart <20>                            将全基因组标记等分,分部分进行关联分析,以节省内存。
 --exact                                 利用精确的方法进行全基因组关联分析
 --num_pair <100000>                     用以计算效应估计值近似标准误的随机标记对数
 --parallel <ith,total>                  全基因组标记等分为total份,仅对第i份进行分析
 --snp_set <start_snp,end_snp>           起始和终止SNP名字,此区段内的标记与全基因组标记进行互作检验.
 --p_cut <1.0e-4>                        利用近似方式进行关联分析时,低于此P值的SNP利用精确方法重新计算P值
 --mkl_num_threads <max>                 MKL线程数，默认为当前服务器最大线程数
 -h, --help                              打印此帮助信息并退出
```

#### （2）输入文件

数据文件：包含文件头，第一列必须是个体号，后面是固定效应和表型，缺失值用***NA***表示（必须是大写的***NA***）。

| id    | mean | sex  | age  | treat | trait |
| ----- | ---- | ---- | ---- | ----- | ----- |
| 14462 | 1    | 0    | 126  | 0     | 0.58  |
| 14463 | 1    | 0    | 91   | 1     | 0.39  |
| 14464 | 1    | 1    | 126  | 0     | 0.37  |
| 14465 | 1    | 0    | 91   | 1     | NA    |
| 14466 | 1    | 0    | 91   | 1     | 0.84  |
| 14467 | 1    | 0    | 91   | 1     | 0.61  |
| 14468 | 1    | 1    | 91   | 1     | NA    |
| 14469 | 1    | 1    | 91   | 1     | 0.59  |
| 14470 | 1    | 1    | 124  | 0     | 0.13  |

加性基因组关系矩阵文件：gmatrix模块构建的加性基因组关系矩阵，包含关系矩阵和个体号文件，关系矩阵为个体号、个体号、值格式。可调用gmatrix --bfile \<PLINK\> --grm agrm --out \<FILE\> 实现。

Plink Binary PED files (*.bed, *.bim, *.fam)：参考PLINK软件说明书 

#### （3）输出文件

preix.遗传效应类型.起始标记位置\_终止标记位置.res：最终结果文件，0-7列为第一个标记的信息，8-15为第二个标记的信息，eff为效应估计值，se为标准误，p为显著性p值。

| order0 | chro0 | snp_id0         | cm0    | base0    | allele0_0 | allele1_0 | freq0    | order1 | chro1 | snp_id1          | cm1    | base1    | allele0_1 | allele1_1 | freq1    | eff      | se       | p           |
| ------ | ----- | --------------- | ------ | -------- | --------- | --------- | -------- | ------ | ----- | ---------------- | ------ | -------- | --------- | --------- | -------- | -------- | -------- | ----------- |
| 17     | 1     | mCV23509126     | 84.458 | 0        | T         | A         | 0.603144 | 1184   | 15    | rs13482580       | 20.003 | 53022566 | T         | A         | 0.508436 | -0.09924 | 0.025462 | 9.72E-05    |
| 18     | 1     | CEL-1-182091301 | 84.931 | 0        | T         | A         | 0.619248 | 1181   | 15    | rs3660290        | 19.176 | 50857204 | T         | A         | 0.533742 | -0.10317 | 0.025842 | 6.54E-05    |
| 18     | 1     | CEL-1-182091301 | 84.931 | 0        | T         | A         | 0.619248 | 1182   | 15    | rs3692040        | 19.262 | 51075548 | T         | A         | 0.533359 | -0.10315 | 0.02584  | 6.55E-05    |
| 18     | 1     | CEL-1-182091301 | 84.931 | 0        | T         | A         | 0.619248 | 1183   | 15    | rs6165881        | 19.453 | 51817647 | T         | A         | 0.521856 | -0.1063  | 0.025636 | 3.37E-05    |
| 18     | 1     | CEL-1-182091301 | 84.931 | 0        | T         | A         | 0.619248 | 1184   | 15    | rs13482580       | 20.003 | 53022566 | T         | A         | 0.508436 | -0.11057 | 0.025723 | 1.72E-05    |
| 44     | 1     | rs3663706       | 20.389 | 42013488 | T         | A         | 0.500767 | 1260   | 16    | rs4219239        | 53.25  | 91670164 | T         | A         | 0.519172 | 0.092843 | 0.024599 | 0.000160505 |
| 116    | 1     | rs6189020       | 54.7   | 1.25E+08 | T         | A         | 0.550613 | 811    | 10    | CEL-10-113177617 | 61.02  | 0        | T         | A         | 0.560583 | -0.11908 | 0.028256 | 2.51E-05    |
| 117    | 1     | rs3687720       | 54.7   | 1.25E+08 | T         | A         | 0.55023  | 811    | 10    | CEL-10-113177617 | 61.02  | 0        | T         | A         | 0.560583 | -0.11915 | 0.028258 | 2.48E-05    |
| 118    | 1     | rs13476094      | 54.7   | 1.25E+08 | T         | A         | 0.55023  | 811    | 10    | CEL-10-113177617 | 61.02  | 0        | T         | A         | 0.560583 | -0.11915 | 0.028258 | 2.48E-05    |
| 119    | 1     | rs13476095      | 54.7   | 1.26E+08 | T         | A         | 0.550613 | 811    | 10    | CEL-10-113177617 | 61.02  | 0        | T         | A         | 0.560583 | -0.11922 | 0.02826  | 2.46E-05    |

preix.遗传效应类型.起始标记位置\_终止标记位置.approx：为近似检验结果文件

preix.遗传效应类型.起始标记位置\_终止标记位置.random_pair：随机生成的SNP对

preix.遗传效应类型.起始标记位置\_终止标记位置.random_pair.res：随机生成的SNP对的检验结果

#### （4）示例代码

**加加互作检验**

注：仅检验加加互作时，模型中必含加性和加加互作关系矩阵，其他类型矩阵（显性、加显互作、显显互作）非必须

```
# 示例数据在mouse文件夹下
# 加性和加加上位关系矩阵构建
gmatrix --bfile plink --grm agrm --out test
gmatrix --bfile plink --grm aagrm --out test

# 仅估计方差组分
remma --data pheno --grm test.agrm.id_fmt,test.aagrm.id_fmt --class sex,treat --covar age --trait trait --out test

# 近似检验：随机筛选100000对标记进行检验，获取效应估计值近似标准误，应用于全基因组标记进行近似检验，P值低于1.0e-4的互作对将重新计算P值
remma --epis abya --bfile plink --data pheno --grm test.agrm.id_fmt,test.aagrm.id_fmt --class sex,treat --covar age --trait trait --out test
# test.abya.0_1407.res 为主要结果文件

# 精确检验，仅保留P值小于1.0e-4的位点，速度慢，不推荐
remma --epis abya --bfile plink --data pheno --grm test.agrm.id_fmt,test.aagrm.id_fmt --class sex,treat --covar age --trait trait --out test --exact
# test.abya.0_1407.res 为主要结果文件

# 等分为3份进行全基因组关联分析，并行提交，缩短运算时间
remma --data pheno --grm test.agrm.id_fmt,test.aagrm.id_fmt --class sex,treat --covar age --trait trait --out test # 仅估计方差组分
remma --epis abya --bfile plink --data pheno --grm test.agrm.id_fmt,test.aagrm.id_fmt --class sex,treat --covar age --trait trait --continue --maxiter 0 --out test --parallel 1,3
remma --epis abya --bfile plink --data pheno --grm test.agrm.id_fmt,test.aagrm.id_fmt --class sex,treat --covar age --trait trait --continue --maxiter 0 --out test --parallel 2,3
remma --epis abya --bfile plink --data pheno --grm test.agrm.id_fmt,test.aagrm.id_fmt --class sex,treat --covar age --trait trait --continue --maxiter 0 --out test --parallel 3,3
# 加入"--continue --maxiter 0"表示以null模型估计的方差组分为初值，不迭代，直接用于全基因组关联分析，避免多次估计方差组分耗时
# test.abya.0_259.res、test.abya.259_595.res和test.abya.595_1407.res为结果文件，合并后为全基因组分析结果

******

# 示例数据在yeast下，包含重复测量值
# 加性和加加上位关系矩阵构建
gmatrix --bfile yeast --grm agrm --out test
gmatrix --bfile yeast --grm aagrm --out test

# 仅估计方差组分
remma --repeat --data pheno --grm test.agrm.id_fmt,test.aagrm.id_fmt --class batch --trait trait --out test

# 近似检验：随机筛选100000对标记进行检验，获取效应估计值近似标准误，应用于全基因组标记进行近似检验，P值低于1.0e-4的互作对将重新计算P值
remma --epis abya --repeat --bfile yeast --data pheno --grm test.agrm.id_fmt,test.aagrm.id_fmt --class batch --trait trait --out test
# test.abya.0_28220.res 为主要结果文件

# 精确检验，仅保留P值小于1.0e-4的位点，速度慢，不推荐
remma --epis abya --repeat --bfile yeast --data pheno --grm test.agrm.id_fmt,test.aagrm.id_fmt --class batch --trait trait --out test --exact
# test.abya.0_28220.res 为主要结果文件

# 等分为5份进行全基因组关联分析，并行提交，缩短运算时间，推荐
remma --repeat --data pheno --grm test.agrm.id_fmt,test.aagrm.id_fmt --class batch --trait trait --out test # 仅估计方差组分
remma --epis abya --repeat --bfile yeast --data pheno --grm test.agrm.id_fmt,test.aagrm.id_fmt --class batch --trait trait --continue --maxiter 0 --out test --parallel 1,5
remma --epis abya --repeat --bfile yeast --data pheno --grm test.agrm.id_fmt,test.aagrm.id_fmt --class batch --trait trait --continue --maxiter 0 --out test --parallel 2,5
remma --epis abya --repeat --bfile yeast --data pheno --grm test.agrm.id_fmt,test.aagrm.id_fmt --class batch --trait trait --continue --maxiter 0 --out test --parallel 3,5
remma --epis abya --repeat --bfile yeast --data pheno --grm test.agrm.id_fmt,test.aagrm.id_fmt --class batch --trait trait --continue --maxiter 0 --out test --parallel 4,5
remma --epis abya --repeat --bfile yeast --data pheno --grm test.agrm.id_fmt,test.aagrm.id_fmt --class batch --trait trait --continue --maxiter 0 --out test --parallel 5,5
# 加入"--continue --maxiter 0"表示以null模型估计的方差组分为初值，不迭代，直接用于全基因组关联分析，避免多次估计方差组分耗时
# test.abya.0_5644.res、test.abya.5644_11288.res、test.abya.11288_16932.res、test.abya.16932_22576.res和test.abya.22576_28220.res为结果文件，合并后为全基因组分析结果
```

**加显互作检验**

注：仅检验加加互作时，模型需包含加性、显性、加加互作、加显互作、显显互作矩阵

```
# 示例数据在mouse文件夹下
# 构建关系矩阵
gmatrix --bfile plink --grm agrm --out test
gmatrix --bfile plink --grm dgrm_as --out test
gmatrix --bfile plink --grm aagrm --out test
gmatrix --bfile plink --grm adgrm_as --out test
gmatrix --bfile plink --grm ddgrm_as --out test

# 仅估计方差组分
remma --data pheno --grm test.agrm.id_fmt,test.dgrm_as.id_fmt,test.aagrm.id_fmt,test.adgrm_as.id_fmt,test.ddgrm_as.id_fmt --class sex,treat --covar age --trait trait --out test

# 近似检验：随机筛选100000对标记进行检验，获取效应估计值近似标准误，应用于全基因组标记进行近似检验，P值低于1.0e-4的互作对将重新计算P值
remma --epis abyd --bfile plink --data pheno --grm test.agrm.id_fmt,test.dgrm_as.id_fmt,test.aagrm.id_fmt,test.adgrm_as.id_fmt,test.ddgrm_as.id_fmt --class sex,treat --covar age --trait trait --out test
# test.abyd.0_1407.res 为主要结果文件

# 精确检验，仅保留P值小于1.0e-4的位点，速度慢，不推荐
remma --epis abyd --bfile plink --data pheno --grm test.agrm.id_fmt,test.dgrm_as.id_fmt,test.aagrm.id_fmt,test.adgrm_as.id_fmt,test.ddgrm_as.id_fmt --class sex,treat --covar age --trait trait --out test --exact
# test.abyd.0_1407.res 为主要结果文件

# 等分为3份进行全基因组关联分析，并行提交，缩短运算时间
remma --data pheno --grm test.agrm.id_fmt,test.dgrm_as.id_fmt,test.aagrm.id_fmt,test.adgrm_as.id_fmt,test.ddgrm_as.id_fmt --class sex,treat --covar age --trait trait --out test # 仅估计方差组分
remma --epis abyd --bfile plink --data pheno --grm test.agrm.id_fmt,test.dgrm_as.id_fmt,test.aagrm.id_fmt,test.adgrm_as.id_fmt,test.ddgrm_as.id_fmt --class sex,treat --covar age --trait trait --continue --maxiter 0 --out test --parallel 1,3
remma --epis abyd --bfile plink --data pheno --grm test.agrm.id_fmt,test.dgrm_as.id_fmt,test.aagrm.id_fmt,test.adgrm_as.id_fmt,test.ddgrm_as.id_fmt --class sex,treat --covar age --trait trait --continue --maxiter 0 --out test --parallel 2,3
remma --epis abyd --bfile plink --data pheno --grm test.agrm.id_fmt,test.dgrm_as.id_fmt,test.aagrm.id_fmt,test.adgrm_as.id_fmt,test.ddgrm_as.id_fmt --class sex,treat --covar age --trait trait --continue --maxiter 0 --out test --parallel 3,3
# 加入"--continue --maxiter 0"表示以null模型估计的方差组分为初值，不迭代，直接用于全基因组关联分析，避免多次估计方差组分耗时
# test.abyd.0_469.res、test.abyd.469_938.res和test.abyd.938_1407.res为结果文件，合并后为全基因组分析结果

```

**显显互作检验**

注：仅检验加加互作时，模型需包含加性、显性、加加互作、加显互作、显显互作矩阵

```
# 示例数据在mouse文件夹下
# 构建关系矩阵
gmatrix --bfile plink --grm agrm --out test
gmatrix --bfile plink --grm dgrm_as --out test
gmatrix --bfile plink --grm aagrm --out test
gmatrix --bfile plink --grm adgrm_as --out test
gmatrix --bfile plink --grm ddgrm_as --out test

# 仅估计方差组分
remma --data pheno --grm test.agrm.id_fmt,test.dgrm_as.id_fmt,test.aagrm.id_fmt,test.adgrm_as.id_fmt,test.ddgrm_as.id_fmt --class sex,treat --covar age --trait trait --out test

# 近似检验：随机筛选100000对标记进行检验，获取效应估计值近似标准误，应用于全基因组标记进行近似检验，P值低于1.0e-4的互作对将重新计算P值
remma --epis dbyd --bfile plink --data pheno --grm test.agrm.id_fmt,test.dgrm_as.id_fmt,test.aagrm.id_fmt,test.adgrm_as.id_fmt,test.ddgrm_as.id_fmt --class sex,treat --covar age --trait trait --out test
# test.dbyd.0_1407.res 为主要结果文件

# 精确检验，仅保留P值小于1.0e-4的位点，速度慢，不推荐
remma --epis dbyd --bfile plink --data pheno --grm test.agrm.id_fmt,test.dgrm_as.id_fmt,test.aagrm.id_fmt,test.adgrm_as.id_fmt,test.ddgrm_as.id_fmt --class sex,treat --covar age --trait trait --out test --exact
# test.dbyd.0_1407.res 为主要结果文件

# 等分为3份进行全基因组关联分析，并行提交，缩短运算时间
remma --data pheno --grm test.agrm.id_fmt,test.dgrm_as.id_fmt,test.aagrm.id_fmt,test.adgrm_as.id_fmt,test.ddgrm_as.id_fmt --class sex,treat --covar age --trait trait --out test # 仅估计方差组分
remma --epis dbyd --bfile plink --data pheno --grm test.agrm.id_fmt,test.dgrm_as.id_fmt,test.aagrm.id_fmt,test.adgrm_as.id_fmt,test.ddgrm_as.id_fmt --class sex,treat --covar age --trait trait --continue --maxiter 0 --out test --parallel 1,3
remma --epis dbyd --bfile plink --data pheno --grm test.agrm.id_fmt,test.dgrm_as.id_fmt,test.aagrm.id_fmt,test.adgrm_as.id_fmt,test.ddgrm_as.id_fmt --class sex,treat --covar age --trait trait --continue --maxiter 0 --out test --parallel 2,3
remma --epis dbyd --bfile plink --data pheno --grm test.agrm.id_fmt,test.dgrm_as.id_fmt,test.aagrm.id_fmt,test.adgrm_as.id_fmt,test.ddgrm_as.id_fmt --class sex,treat --covar age --trait trait --continue --maxiter 0 --out test --parallel 3,3
# 加入"--continue --maxiter 0"表示以null模型估计的方差组分为初值，不迭代，直接用于全基因组关联分析，避免多次估计方差组分耗时
# test.dbyd.0_259.res、test.dbyd.259_595.res和test.dbyd.595_1407.res为结果文件，合并后为全基因组分析结果
```

**加性检验**

注：至少包含加性基因组关系矩阵，可包含多个，建议使用uvlmm模块进行加性检验

```
# 示例数据在mouse文件夹下
# 构建关系矩阵
gmatrix --bfile plink --grm agrm --out test

# 仅估计方差组分
remma --data pheno --grm test.agrm.id_fmt --class sex,treat --covar age --trait trait --out test

# 检验
remma --epis add --bfile plink --data pheno --grm test.agrm.id_fmt --class sex,treat --covar age --trait trait --out test
# test.add.0_1407.res 为主要结果文件
```

**显性检验**

注：至少包含加性和显性基因组关系矩阵，可包含多个

```
# 示例数据在mouse文件夹下
# 构建关系矩阵
gmatrix --bfile plink --grm agrm --out test
gmatrix --bfile plink --grm dgrm_as --out test

# 仅估计方差组分
remma --data pheno --grm test.agrm.id_fmt,test.dgrm_as.id_fmt --class sex,treat --covar age --trait trait --out test

# 检验
remma --epis dom --bfile plink --data pheno --grm test.agrm.id_fmt,test.dgrm_as.id_fmt --class sex,treat --covar age --trait trait --out test
# test.dom.0_1407.res 为主要结果文件
```



