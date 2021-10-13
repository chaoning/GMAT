# 使用说明-中文版

## 1. 联系方式
ningchao(at)sdau(dot)edu(dot)cn

## 2. 安装
所有模块均已在64位linux系统下编译完成,下载并修改为可执行权限即可使用。 
修改文件权限命令:chmod 777 文件

## 3.功能模块



### 3.1 关系矩阵构建
参数详解: 
```
 -b, --bfile [FILE]			Plink Binary PED files 前缀
 -d, --dfile [FILE] 		Dosage基因型文件全名
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
 -n, --npart [num, default is 10] 			将全基因组标记等分,分部分读取计算,节省内存
 -f, --fmt [num, default is 2] 				关系矩阵输出格式:0,“矩阵”格式;1,“行号、列号、值”格式;2,“个体号、个体号、值”格式
 -i, --inv [No argument] 					求逆矩阵
 -m, --mkl [num, default is 10] 			MKL库线程数
 -v, --val [float, Default is 0.001] 		加到矩阵对角线的值,保证矩阵正定
 -s, --skipcols [int, Default is 1] 		读取文件时跳过的列数,仅对Dosage基因型文件(--dfile)有效
 -h, --help 								打印此帮助信息
```
#### 加性基因组关系矩阵
读取Plink Binary PED files (plink.bed, plink.bim, plink.fam),根据需求,可生成不同输出格式的加性基因组关系矩阵。 
参考文献:VanRaden PM. Efficient methods to compute genomic predictions. J Dairy Sci. 2008 Nov;91(11):4414-23. 
```
# 1. 默认参数输出
gmatrix --bfile plink --grm agrm --out test 
# 输出结果包含两个文件: test.id文件包含个体id号;test.agrm.id_fmt包含加性基因组关系矩阵信息, 为三列:个体号、个体号、值。

# 2. 输出格式为行号、列号、值
gmatrix --bfile plink --grm agrm --fmt 1 --out test 

# 3. 输出格式为矩阵形式
gmatrix --bfile plink --grm agrm --fmt 0 --out test 

# 4. 输出逆矩阵
gmatrix --bfile plink --grm agrm --inv --out test 

# 5. 将标记分部分读取、依次计算以节省内存(增加--npart参数值,默认为10)
gmatrix --bfile plink --grm agrm --npart 30 --out test 

# 6. 增加mkl线程数,提高运算效率(默认为10)
gmatrix --bfile plink --grm agrm --mkl_num_threads 30 --out test 

# 7. 矩阵对角线上加值,保证正定(默认为0.001)
gmatrix --bfile plink --grm agrm --val 0.01 --out test 

# 8. 读取Dosage基因型文件,计算加性基因组关系矩阵
gmatrix --dfile plink.traw --skipcols 6 --grm agrm_dosage --out test 

```
#### 显性基因组关系矩阵
参考文献:Zulma G Vitezica, Luis Varona, Andres Legarra, On the Additive and Dominant Variance and Covariance of Individuals Within the Genomic Selection Scope, Genetics, Volume 195, Issue 4, 1 December 2013, Pages 1223–1230. 
```
# 1. 构建生物学意义的显性基因组关系矩阵
gmatrix --bfile plink --grm dgrm_as --out test

# 2. 构建育种学意义的显性基因组关系矩阵
gmatrix --bfile plink --grm dgrm_gs --out test

```
#### 加加上位基因组关系矩阵
依据加性和加性基因组关系矩阵直积计算得到
```
gmatrix --bfile plink --grm aagrm --out test
```
#### 加显上位基因组关系矩阵
依据加性和显性基因组关系矩阵直积计算得到
```
# 1. 显性矩阵为生物学意义的关系矩阵
gmatrix --bfile plink --grm adgrm_as --out test
# 2. 显性矩阵为育种学意义的关系矩阵
gmatrix --bfile plink --grm adgrm_gs --out test

```
#### 显显上位基因组关系矩阵

依据显性和显性基因组关系矩阵直积计算得到
```
# 1. 显性矩阵为生物学意义的关系矩阵
gmatrix --bfile plink --grm ddgrm_as --out test
# 2. 显性矩阵为育种学意义的关系矩阵
gmatrix --bfile plink --grm ddgrm_gs --out test
```
### 3.2 单性状全基因组关联分析
参数详解: 
```
 -d, --data [FILE]                         数据文件,包含表头(列名),第一列必须是个体号,缺失值用NA表示
 -b, --bfile [FILE]                        Plink Binary PED files 前缀
 -o, --out [FILE]                          输出文件前缀
 -g, --grm [FILE]                          加性基因组关系矩阵前缀,
 -t, --trait [str]                         待分析性状列名
 --repeat                              		使用重复力模型(有重复测量值的性状,例如,猪产仔数)
 --covar=str1,str2,...,str(n)            	放入固定效应的连续变量列名
 --class=str1,str2,...,str(n)            	放入固定效应的分类变量列名
 --snp_set=start_snp,end_snp            	起始和终止SNP名字,对此区段的SNP进行关联分析
 --continue                       			估计方差组分时利用上一轮方差组分结果继续迭代
 --condition=SNP1,SNP2,...,SNP(n)        	加入固定效应中的条件SNP名字
 --maxiter=200                    			不加SNP标记时估计方差组分的最大迭代次数
 --maxiter0=30                    			对每个标记逐个检验时的方差组分迭代次数
 --cc_par=1.0e-8                  			方差组分收敛标准,方差组分变化差值的平方和的平方根
 --cc_gra=1.0e-3                  			方差组分收敛标准,一阶偏导的平方和的平方根
 --npart=20                     			将全基因组标记等分,分部分进行关联分析,以节省内存。
 --exact                       				利用精确的方法进行全基因组关联分析,每个标记均进行方差组分估计。
 --p_cut=1.0e-3                   			利用近似方式进行关联分析时,低于此P值的SNP利用精确方法重新计算P值。
 --parallel ith,total               		全基因组标记等分为total份,仅对第i份进行分析。
 --mkl_num_threads=10                		MKL线程数
 -h, --help                       			打印此帮助信息并退出。
```
```
# 计算加性基因组关系矩阵
gmatrix --bfile plink --grm agrm --out test

# 仅估计方差组分
uvlmm --data pheno --grm test --trait trait --class sex,treat --covar age --out test

# 全基因组关联分析
uvlmm --bfile plink --data pheno --grm test --trait trait --class sex,treat --covar age --out test

# 将标记等分为3份进行全基因组关联分析
uvlmm --bfile plink --data pheno --grm test --trait trait --class sex,treat --covar age --out test --parallel 1,3
uvlmm --bfile plink --data pheno --grm test --trait trait --class sex,treat --covar age --out test --parallel 2,3
uvlmm --bfile plink --data pheno --grm test --trait trait --class sex,treat --covar age --out test --parallel 3,3

# 对一区段内标记进行关联分析
uvlmm --bfile plink --data pheno --grm test --trait trait --class sex,treat --covar age --out test --snp_set mCV25266528,rs3663706

#加入Top SNP进行条件GWAS分析
uvlmm --bfile plink --data pheno --grm test --trait trait --class sex,treat --covar age --out test --condition rs13477831

# 利用重复力模型进行全基因组关联分析(yeast数据为例)
gmatrix --bfile yeast --grm agrm --out test
uvlmm --data pheno --grm test --trait trait --class batch --out test --repeat
uvlmm --bfile yeast --data pheno --grm test --trait trait --class batch --out test --repeat
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
