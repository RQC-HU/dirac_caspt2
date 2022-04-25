# DIRAC-CASPT2

- DIRAC19or21の1,2電子積分ファイルを用いて、CASCI/CASPT2法またはDMRG/CASPT2法で2次の多配置摂動計算を行います
a
## How to Install

以下のコマンドでmainブランチのソースコードをビルドできます

```sh
    git clone https://github.com/kohei-noda-qcrg/dirac_caspt2
    cd dirac_caspt2
    ./configure
    make
```

- デフォルトでは gfortranを使用して実行ファイルがビルドされます
- 現状GNU Fortranはビルドは成功しますが実行時エラーが発生する可能性があるため非推奨です
- したがって./configure FC=ifort もしくは ./confiugre FC=mpiifort を使用してください

### Requirements

以下のコンパイラおよびライブラリと依存性があり、ビルドを行う計算機でこれらがセットアップされている必要があります

- [GNU Fortran](https://gcc.gnu.org/fortran/) or [Intel Fortran](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html) compiler (並列計算をするために並列コンパイラを使うこともできます)
  - 現時点では並列コンパイラについて、Intel Fortranはmpiifort,GNU Fortranはmpifort以外は./configureで手動設定が必要です
- [Intel MKL(Math Kernel Library)](https://www.intel.com/content/www/us/en/develop/documentation/get-started-with-mkl-for-dpcpp/top.html)
  - 現時点ではMKLのBlas,LapackではなくBlas及びLapack単体でビルドする場合./configureのLDFLAGSを手動設定する必要があります
  - gfortranを使用する場合、MKLをリンクするため環境変数\$MKLROOTが設定されている必要があります
    \$MKLROOTが設定されているか確認するには、使用する計算機にログインして以下のコマンドを実行してMKLにパスが通っているかを確認してください

    ```sh
        echo $MKLROOT
    ```

### ビルド例

各種コンパイラは\$PATHに追加されているか、もしくはフルパスを指定する必要があります

- Intel Fortran

```sh
    ./configure FC=ifort
    make
```

- Intel Fortran(parallel)

```sh
    ./configure FC=mpiifort
    make
```

- GNU Fortran

```sh
    ./configure FC=gfortran
    make
```

- GNU Fortran(parallel)

```sh
    ./configure FC=mpifort
    make
```

### HELP

ヘルプを表示するには以下を実行します

```sh
    ./configure --help
```
