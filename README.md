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
 -b, --bfile [FILE]			Plink Binary PED files前缀
 -d, --dfile [FILE] 		Dosage基因型文件全名
 -a, --fam [FILE]            PLINK格式的样本信息文件全名
 -o, --out [FILE] 			输出文件前缀
 -g, --grm [string] 		基因组关系矩阵类型
      agrm: 				加性基因组关系矩阵
      dgrm_as: 				显性基因组关系矩阵,生物学意义显性矩阵,可用于关联分析
      dgrm_gs: 				显性基因组关系矩阵,育种学意义显性矩阵,用于显性效应剖分、基因组选择
      aagrm: 				加加上位基因组关系矩阵
      adgrm_as: 			加显上位基因组关系矩阵,其中显性矩阵是生物学意义显性矩阵
      adgrm_gs: 			加显上位基因组关系矩阵,其中显性矩阵是育种学意义显性矩阵
      ddgrm_as: 			显显上位基因组关系矩阵,其中显性矩阵是生物学意义显性矩阵        
      ddgrm_gs: 			显显上位基因组关系矩阵,其中显性矩阵是育种学意义显性矩阵   
      agrm_dosage: 			加性基因组关系矩阵,读取Dosage基因型文件计算   
 -n, --npart [num, default is 20] 			将全基因组标记等分,分部分读取计算,节省内存
 -f, --fmt [num, default is 2] 				关系矩阵输出格式:0,“矩阵”格式;1,“行号、列号、值”格式;2,“个体号、个体号、值”格式
 -i, --inv [No argument] 					求逆矩阵
 -m, --mkl [num, default is max] 			MKL库线程数，默认为当前服务器最大线程数。
 -v, --val [float, Default is 0.001] 		加到矩阵对角线的值,保证矩阵正定
 -s, --skipcols [int, Default is 1] 		读取文件时跳过的列数,仅对Dosage基因型文件(--dfile)有效
 -h, --help 								打印此帮助信息
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

关系矩阵文件：*.<grm>.<fmt>_fmt，依据给出--fmt参数的不同，输出不同格式的基因组关系矩阵，--fmt 0输出矩阵格式，--fmt 1输出行号、列号、值形式，--fmt 2输出个体号、个体号、值形式。

关系矩阵的逆矩阵文件：*.<giv>.<fmt>_fmt，依据给出--fmt参数的不同，输出不同格式的基因组关系矩阵的逆矩阵。

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
 -g, --grm <FILE>                          加性基因组关系矩阵前缀,
 -t, --trait <str>                         待分析性状在数据文件中的列名
 --repeat                              		使用重复力模型(有重复测量值的性状,例如,猪产仔数)
 --covar <str1,str2,...,str(n)>            	放入模型中作为固定效应的连续变量列名
 --class <str1,str2,...,str(n)>            	放入模型中作为固定效应的分类变量列名
 --snp_set <start_snp,end_snp>            	起始和终止SNP名字,对此区段的SNP进行关联分析
 --continue                       			估计方差组分时利用上一轮方差组分结果继续迭代
 --condition <SNP1,SNP2,...,SNP(n)>        	加入到固定效应中的SNP名字
 --maxiter <200>                     	    不加SNP标记时估计方差组分的最大迭代次数
 --maxiter0 <30>                    		对每个标记逐个检验时的方差组分迭代次数
 --cc_par <1.0e-8>                  		方差组分收敛标准,方差组分变化差值的平方和的平方根
 --cc_gra <1.0e-6>                  		方差组分收敛标准,一阶偏导的平方和的平方根
 --npart <20>                     			将全基因组标记等分,分部分进行关联分析,以节省内存。
 --exact                       				利用精确的方法进行全基因组关联分析,每个标记均进行方差组分估计。
 --p_cut <1.0e-3>                   		利用近似方式进行关联分析时,低于此P值的SNP利用精确方法重新计算P值。
 --parallel <ith,total>               		全基因组标记等分为total份,仅对第i份进行分析。
 --mkl_num_threads <max>                	MKL线程数，默认为当前服务器最大线程数。
 -h, --help                       			打印此帮助信息并退出。
```
#### （2）输入文件：

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

加性基因组关系矩阵文件：gmatrix模块构建的加性基因组关系矩阵，包含关系矩阵和个体号文件，关系矩阵为个体号、个体号、值格式。可调用gmatrix --bfile <PLINK> --grm agrm --out <FILE> 实现。

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
uvlmm --bfile plink --data pheno --grm test --trait trait --class sex,treat --covar age --out test --condition rs13477831

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



### 3.3 多性状全基因组关联分析
**平衡数据:不同性状在所有个体中均有测量值,无缺失数据**  
参数详解:

```
 -d, --data [FILE]                        数据文件,包含表头(列名),第一列必须是个体号,缺失值用NA表示
 -b, --bfile [FILE]                        Plink Binary PED files 前缀
 -o, --out [FILE]                          输出文件前缀
 -g, --grm [FILE]                          加性基因组关系矩阵前缀,
 -t, --trait [str1,str2,...,str(n)]              待分析性状列名
 --repeat                              使用重复力模型(有重复测量值的性状,例如,猪产仔数)
 --covar=str1,str2,...,str(n)            放入固定效应的连续变量列名
 --class=str1,str2,...,str(n)            放入固定效应的分类变量列名
 --snp_set=start_snp,end_snp            起始和终止SNP名字,对此区段的SNP进行关联分析
 --continue                       估计方差组分时利用上一轮方差组分结果继续迭代
 --condition=SNP1,SNP2,...,SNP(n)        加入固定效应中的条件SNP名字
 --maxiter=200                    不加SNP标记时估计方差组分的最大迭代次数
 --maxiter0=30                    对每个标记逐个检验时的方差组分迭代次数
 --cc_par=1.0e-8                  方差组分收敛标准,方差组分变化差值的平方和的平方根
 --cc_gra=1.0e-3                  方差组分收敛标准,一阶偏导的平方和的平方根
 --npart=20                     将全基因组标记等分,分部分进行关联分析,以节省内存。
 --exact                       利用精确的方法进行全基因组关联分析,每个标记均进行方差组分估计。
 --p_cut=1.0e-3                   利用近似方式进行关联分析时,低于此P值的SNP利用精确方法重新计算P值。
 --parallel ith,total               全基因组标记等分为total份,仅对第i份进行分析。
 --mkl_num_threads=10                MKL线程数
 -h, --help                       打印此帮助信息并退出。
```

```
# mouse_long 数据
# 计算加性基因组关系矩阵
gmatrix --bfile plink --grm agrm --out test

# 仅估计方差组分
bmvlmm --data phe.balance.txt --grm test --trait trait1,trait2,trait3 --class sex --out test

# 全基因组关联分析
bmvlmm --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3 --class sex --out test

# 将标记等分为3份进行全基因组关联分析
bmvlmm --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3 --class sex --out test --parallel 1,3
bmvlmm --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3 --class sex --out test --parallel 2,3
bmvlmm --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3 --class sex --out test --parallel 3,3

# 对一区段内标记进行关联分析
bmvlmm --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3 --class sex --out test --snp_set UNC010515443,UNC172447

#加入Top SNP进行条件GWAS分析
bmvlmm --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3 --class sex --out test --condition UNC18848064

```

**非平衡数据:同一个个体的多个性状可能存在缺失测量值**  
参数详解:

```
 -d, --data [FILE]                        数据文件,包含表头(列名),第一列必须是个体号,缺失值用NA表示
 -b, --bfile [FILE]                        Plink Binary PED files 前缀
 -o, --out [FILE]                          输出文件前缀
 -g, --grm [FILE]                          加性基因组关系矩阵前缀,
 -t, --trait [str1,str2,...,str(n)]              待分析性状列名
 --repeat                              使用重复力模型(有重复测量值的性状,例如,猪产仔数)
 --covar=str1,str2,...,str(n)            放入固定效应的连续变量列名
 --class=str1,str2,...,str(n)            放入固定效应的分类变量列名
 --snp_set=start_snp,end_snp            起始和终止SNP名字,对此区段的SNP进行关联分析
 --continue                       估计方差组分时利用上一轮方差组分结果继续迭代
 --condition=SNP1,SNP2,...,SNP(n)        加入固定效应中的条件SNP名字
 --maxiter=200                    不加SNP标记时估计方差组分的最大迭代次数
 --cc_par=1.0e-8                  方差组分收敛标准,方差组分变化差值的平方和的平方根
 --cc_gra=1.0e-3                  方差组分收敛标准,一阶偏导的平方和的平方根
 --npart=20                     将全基因组标记等分,分部分进行关联分析,以节省内存。
 --p_cut=1.0e-3                   利用近似方式进行关联分析时,低于此P值的SNP利用精确方法重新计算P值。
 --parallel ith,total               全基因组标记等分为total份,仅对第i份进行分析。
 --mkl_num_threads=10                MKL线程数
 -h, --help                       打印此帮助信息并退出。
```

```
# mouse_long 数据
# 计算加性基因组关系矩阵
gmatrix --bfile plink --grm agrm --out test

# 仅估计方差组分
unmvlmm --data phe.balance.txt --grm test --trait trait1,trait2,trait3 --class sex --out test

# 全基因组关联分析
unmvlmm --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3 --class sex --out test
# 将标记等分为3份进行全基因组关联分析
unmvlmm --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3 --class sex --out test --parallel 1,3
unmvlmm --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3 --class sex --out test --parallel 2,3
unmvlmm --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3 --class sex --out test --parallel 3,3

# 对一区段内标记进行关联分析
unmvlmm --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3 --class sex --out test --snp_set UNC010515443,UNC172447

#加入Top SNP进行条件GWAS分析
unmvlmm --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3 --class sex --out test --condition UNC18848064
```

### 3.4 纵向性状全基因组关联分析
**平衡数据:不同个体在所有时间点均有测量值,无缺失**  
参数详解:

```
 -d, --data [FILE]                        数据文件,包含表头(列名),第一列必须是个体号,缺失值用NA表示
 -b, --bfile [FILE]                        Plink Binary PED files 前缀
 -o, --out [FILE]                          输出文件前缀
 -g, --grm [FILE]                          加性基因组关系矩阵前缀,
 -t, --trait [str1,str2,...,str(n)]              待分析性状列名
 --tpoint=t1,t2,...,t(n)                测量时间点
 --tcovar=str1,str2,...,str(n)           随时间变化的协变量列名
 --tclass=str1,str2,...,str(n)          随时间变化的分类变量列名
 --covar=str1,str2,...,str(n)            放入固定效应的连续变量列名
 --class=str1,str2,...,str(n)            放入固定效应的分类变量列名
 --snp_set=start_snp,end_snp            起始和终止SNP名字,对此区段的SNP进行关联分析
 --order=num								勒让德多项式阶数
 --continue                       估计方差组分时利用上一轮方差组分结果继续迭代
 --condition=SNP1,SNP2,...,SNP(n)        加入固定效应中的条件SNP名字
 --maxiter=200                    不加SNP标记时估计方差组分的最大迭代次数
 --maxiter0=30                    对每个标记逐个检验时的方差组分迭代次数
 --cc_par=1.0e-8                  方差组分收敛标准,方差组分变化差值的平方和的平方根
 --cc_gra=1.0e-3                  方差组分收敛标准,一阶偏导的平方和的平方根
 --npart=20                     将全基因组标记等分,分部分进行关联分析,以节省内存。
 --exact                       利用精确的方法进行全基因组关联分析,每个标记均进行方差组分估计。
 --p_cut=1.0e-3                   利用近似方式进行关联分析时,低于此P值的SNP利用精确方法重新计算P值。
 --parallel ith,total               全基因组标记等分为total份,仅对第i份进行分析。
 --mkl_num_threads=10                MKL线程数
 -h, --help                       打印此帮助信息并退出。
```

```
# mouse_long 数据
# 计算加性基因组关系矩阵
gmatrix --bfile plink --grm agrm --out test

# 仅估计方差组分
blongwas --data phe.balance.txt --grm test --trait trait1,trait2,trait3,trait4,trait5,trait6,trait7,trait8,trait9,trait10,trait11,trait12,trait13,trait14,trait15,trait16 --tpoint 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 --tclass sex --out test

# 全基因组关联分析
blongwas --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3,trait4,trait5,trait6,trait7,trait8,trait9,trait10,trait11,trait12,trait13,trait14,trait15,trait16 --tpoint 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 --tclass sex --out test

# 将标记等分为3份进行全基因组关联分析
blongwas --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3,trait4,trait5,trait6,trait7,trait8,trait9,trait10,trait11,trait12,trait13,trait14,trait15,trait16 --tpoint 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 --tclass sex --out test --parallel 1,3
blongwas --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3,trait4,trait5,trait6,trait7,trait8,trait9,trait10,trait11,trait12,trait13,trait14,trait15,trait16 --tpoint 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 --tclass sex --out test --parallel 2,3
blongwas --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3,trait4,trait5,trait6,trait7,trait8,trait9,trait10,trait11,trait12,trait13,trait14,trait15,trait16 --tpoint 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 --tclass sex --out test --parallel 3,3

# 对一区段内标记进行关联分析
blongwas --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3,trait4,trait5,trait6,trait7,trait8,trait9,trait10,trait11,trait12,trait13,trait14,trait15,trait16 --tpoint 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 --tclass sex --out test --snp_set UNC010515443,UNC172447

#加入Top SNP进行条件GWAS分析
blongwas --bfile plink --data phe.balance.txt --grm test --trait trait1,trait2,trait3,trait4,trait5,trait6,trait7,trait8,trait9,trait10,trait11,trait12,trait13,trait14,trait15,trait16 --tpoint 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 --tclass sex --out test --condition UNC18848064
```

**非平衡数据:同一个体在不同时间点可能存在缺失测量值**  
参数详解:

```
-d, --data [FILE]                        数据文件,包含表头(列名),第一列必须是个体号,缺失值用NA表示
 -b, --bfile [FILE]                        Plink Binary PED files 前缀
 -o, --out [FILE]                          输出文件前缀
 -g, --grm [FILE]                          加性基因组关系矩阵前缀,
 -t, --trait str              					待分析性状列名
 --tpoint str                			测量时间点列名
 --tcovar=str1,str2,...,str(n)           随时间变化的协变量列名
 --tclass=str1,str2,...,str(n)          随时间变化的分类变量列名
 --covar=str1,str2,...,str(n)            放入固定效应的连续变量列名
 --class=str1,str2,...,str(n)            放入固定效应的分类变量列名
 --snp_set=start_snp,end_snp            起始和终止SNP名字,对此区段的SNP进行关联分析
 --order=num1,num2,num3						固定回归、加性遗传效应随机回归勒让德多项式阶数
 --continue                       估计方差组分时利用上一轮方差组分结果继续迭代
 --condition=SNP1,SNP2,...,SNP(n)        加入固定效应中的条件SNP名字
 --maxiter=200                    不加SNP标记时估计方差组分的最大迭代次数
 --cc_par=1.0e-8                  方差组分收敛标准,方差组分变化差值的平方和的平方根
 --cc_gra=1.0e-3                  方差组分收敛标准,一阶偏导的平方和的平方根
 --npart=20                     将全基因组标记等分,分部分进行关联分析,以节省内存。
 --p_cut=1.0e-3                   利用近似方式进行关联分析时,低于此P值的SNP利用精确方法重新计算P值。
 --parallel ith,total               全基因组标记等分为total份,仅对第i份进行分析。
 --mkl_num_threads=10                MKL线程数
 -h, --help                       打印此帮助信息并退出。
```

```
# mouse_long 数据
# 计算加性基因组关系矩阵
gmatrix --bfile plink --grm agrm --out test

# 仅估计方差组分
unlongwas --data phe.unbalance.txt --grm test --trait trait --tpoint weak --tclass Sex --out test

# 全基因组关联分析
unlongwas --bfile plink --data phe.unbalance.txt --grm test --trait trait --tpoint weak --tclass Sex --out test

# 将标记等分为3份进行全基因组关联分析
unlongwas --bfile plink --data phe.unbalance.txt --grm test --trait trait --tpoint weak --tclass Sex --out test --parallel 1,3
unlongwas --bfile plink --data phe.unbalance.txt --grm test --trait trait --tpoint weak --tclass Sex --out test --parallel 2,3
unlongwas --bfile plink --data phe.unbalance.txt --grm test --trait trait --tpoint weak --tclass Sex --out test --parallel 3,3

# 对一区段内标记进行关联分析
unlongwas --bfile plink --data phe.unbalance.txt --grm test --trait trait --tpoint weak --tclass Sex --out test --snp_set UNC010515443,UNC172447

#加入Top SNP进行条件GWAS分析
unlongwas --bfile plink --data phe.unbalance.txt --grm test --trait trait --tpoint weak --tclass Sex --out test --condition UNC18848064
```
