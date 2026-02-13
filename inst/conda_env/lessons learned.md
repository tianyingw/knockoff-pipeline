SAIGE你害得我好苦啊……

历时n天终于装好SAIGE，现在来回忆并整理下完整的步骤（中间踩了太多坑）：

1. 创建 conda 环境

创建环境，然后激活（用我的yml文件！）
```bash
conda env create -f /dssg/home/acct-mashiyang1991/mashiyang1991/yujie/pipeline/conda_env/environment.yml
conda activate RSAIGE
FLAGPATH=`which python | sed 's|/bin/python$||'`
export LDFLAGS="-L${FLAGPATH}/lib"
export CPPFLAGS="-I${FLAGPATH}/include"
```

2. 安装额外依赖

（用我的R文件！）
```bash
Rscript /dssg/home/acct-mashiyang1991/mashiyang1991/yujie/pipeline/conda_env/install_packages.R
```

3. 从源码安装SAIGE

下载库
```bash
src_branch=main
repo_src_url=https://github.com/saigegit/SAIGE
git clone --depth 1 -b $src_branch $repo_src_url
```

接下来哐哐填补库中漏洞（也可能是我“打开方式”不对）：
首先进入SAIGE文件夹目录

(1) 下载plink2源码（include .h文件时会用到）
```bash
git clone https://github.com/chrchang/plink-ng.git
```
(2) 修改src/Makevars:
修改后内容见/dssg/home/acct-mashiyang1991/mashiyang1991/yujie/pipeline/conda_env/Makevars

(3) 获取plink2的静态链接库文件
```r
> library(pgenlibr)
> lib_path <- system.file("libs", package = "pgenlibr")
> files <- list.files(lib_path, full.names = TRUE, recursive = TRUE)
> a_files <- list.files(lib_path, pattern = "\\.a$", full.names = TRUE, recursive = TRUE)
> a_files
[1] "/dssg/home/acct-mashiyang1991/mashiyang1991/.conda/envs/RSAIGE/lib/R/library/pgenlibr/libs/libPGZSTD.a"
[2] "/dssg/home/acct-mashiyang1991/mashiyang1991/.conda/envs/RSAIGE/lib/R/library/pgenlibr/libs/libPLINK2.a"
```
复制libPLINK2.a文件到SAIGE目录下，并改名plink2_includes.a 

4. 编译和安装SAIGE

激动人心的时刻终于到来了，不出意外的话以上流程走完不会报错。祝好运！
```bash
R CMD INSTALL SAIGE
```

5. 安装完毕！开始使用吧
```r
library(SAIGE)
```