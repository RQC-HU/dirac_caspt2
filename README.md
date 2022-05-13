# DIRAC-CASPT2

#### [DIRAC](http://diracprogram.org/doku.php)の計算結果のうち1,2電子積分ファイルを用いて、CASCI/CASPT2法またはDMRG/CASPT2法で2次の多配置摂動計算を行います

## 目次
- [Requirements](https://github.com/kohei-noda-qcrg/dirac_caspt2#requirements)
- [How to Install](https://github.com/kohei-noda-qcrg/dirac_caspt2#how-to-install)
  - [ビルドオプション](https://github.com/kohei-noda-qcrg/dirac_caspt2#ビルドオプション)
  - [ビルド例](https://github.com/kohei-noda-qcrg/dirac_caspt2#ビルド例)
- [How to use](https://github.com/kohei-noda-qcrg/dirac_caspt2#how-to-use)
- [開発者のかたへ](https://github.com/kohei-noda-qcrg/dirac_caspt2#開発者のかたへ)
  - [環境構築について](https://github.com/kohei-noda-qcrg/dirac_caspt2#環境構築について)
  - [ビルドについて](https://github.com/kohei-noda-qcrg/dirac_caspt2#ビルドについて)
  - [テストについて](https://github.com/kohei-noda-qcrg/dirac_caspt2#テストについて)
## Requirements

以下のコンパイラおよびツール、ライブラリと依存性があり、ビルドを行う計算機でこれらがセットアップされている必要があります

- [GNU Fortran](https://gcc.gnu.org/fortran/) or [Intel Fortran](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html) compiler (並列計算をするために並列コンパイラを使うこともできます)
- [CMake](https://cmake.org/)(version>=3.7 が必要です)
    - cmakeが計算機に入っていないか、バージョンが古い場合[CMakeのGithub](https://github.com/Kitware/CMake/releases)からビルドするもしくはビルド済みのファイルを解凍して使用してください
- [Intel MKL(Math Kernel Library)](https://www.intel.com/content/www/us/en/develop/documentation/get-started-with-mkl-for-dpcpp/top.html)
  - 現時点ではMKLのBlas,LapackではなくBlas及びLapack単体でビルドする場合LDFLAGSを手動設定する必要があります
  - gfortranを使用する場合、MKLをリンクするため環境変数\$MKLROOTが設定されている必要があります
    \$MKLROOTが設定されているか確認するには、使用する計算機にログインして以下のコマンドを実行してMKLにパスが通っているかを確認してください

    ```sh
    echo $MKLROOT
    ```

## How to Install

以下のコマンドでmainブランチのソースコードをビルドできます

```sh
git clone https://github.com/kohei-noda-qcrg/dirac_caspt2
cd dirac_caspt2
mkdir -p build && cd build
FC=ifort cmake ..
make
```

- 現状GNU Fortranはビルドは成功しますが実行時エラーが発生する可能性があるため**非推奨**です
- したがってFC=ifort もしくは FC=mpiifort もしくは FC=mpifort (OpenMPI,ifort) を使用することを推奨します
- CMake version >= 3.13 を使っているなら以下のようなコマンドでもビルドができます

```sh
git clone https://github.com/kohei-noda-qcrg/dirac_caspt2
cd dirac_caspt2
FC=ifort cmake -B build
cmake --build build --clean-first
```

### ビルドオプション

現時点でサポートしているビルドオプションは以下のとおりです

ビルドオプションはcmake -DBUILDOPTION1=on -DBUILDOPTION2=off ,,,のように使います

- MPI
    - MPIを使用するなら必須です.マルチプロセス対応ビルドのためのプリプロセッサの設定を行います
    
        (例)
        ```sh
        mkdir -p build && cd build
        FC=mpiifort cmake -DMPI=on ..
        make
        ```
- OPENMP
    - OpenMPを使用するなら必須です.OpenMP用のビルドオプションを追加します

        (例)
        ```sh
        mkdir -p build && cd build
        FC=ifort cmake -DOPENMP=on ..
        make
        ```
### ビルド例

各種コンパイラは\$PATHに追加されているか、もしくはフルパスを指定する必要があります

- Intel Fortran

    ```sh
    mkdir -p build && cd build
    FC=ifort cmake ..
    make
    ```

- Intel Fortran (with OpenMP)

    ```sh
    mkdir -p build && cd build
    FC=ifort cmake -DOPENMP=on ..
    make
    ```

- Intel Fortran(MPI only, Intel MPI)

    ```sh
    mkdir -p build && cd build
    FC=mpiifort cmake -DMPI=on ..
    make
    ```

- Intel Fortran(MPI/OpenMP hybrid, Intel MPI)

    ```sh
    mkdir -p build && cd build
    FC=mpiifort cmake -DMPI=on -DOPENMP=on ..
    make
    ```

- GNU Fortran

    ```sh
    mkdir -p build && cd build
    FC=gfortran cmake ..
    make
    ```

- GNU Fortran (with OpenMP)

    ```sh
    mkdir -p build && cd build
    FC=gfortran cmake -DOPENMP=on ..
    make
    ```

- OpenMPI Fortran(MPI only)

    ```sh
    mkdir -p build && cd build
    FC=mpifort cmake -DMPI=on ..
    make
    ```

- OpenMPI Fortran(MPI/OpenMP hybrid)

    ```sh
    mkdir -p build && cd build
    FC=mpifort cmake -DMPI=on -DOPENMP=on ..
    make
    ```
## How to use

- active.inpという名前のファイルが必要です(ファイル名は必ずactive.inpとしてください)
- 計算の実行はactive.inpとDIRACの積分ファイル(MDCINT,MRCONEE)があるディレクトリで、ビルドした実行可能ファイル(r4divocoexe, r4dcascicoexe, r4dcaspt2ocoexe)を指定して実行します
- 以下のようなシェルスクリプトを用意して実行すると簡単に実行できます
  - スクリプト(非並列)
  ```sh
  #!/bin/sh

    PGMIVO=/path/to/dirac_caspt2/bin/r4divocoexe
    PGMCASCI=/path/to/dirac_caspt2/bin/r4dcascicoexe
    PGMCASPT2O=/path/to/dirac_caspt2/bin/r4dcaspt2ocoexe

  #-- Execution Sequence ------------------------------------------------------

        $PGMCASCI   &> H2O.caspt2.out
        $PGMCASPT2O &>> H2O.caspt2.out
  ```
  - スクリプト(並列) 
  ```sh
  #!/bin/sh

    PGMIVO=/path/to/dirac_caspt2/bin/r4divocoexe
    PGMCASCI=/path/to/dirac_caspt2/bin/r4dcascicoexe
    PGMCASPT2O=/path/to/dirac_caspt2/bin/r4dcaspt2ocoexe

  #-- Execution Sequence ------------------------------------------------------
  NPROCS=8 # 並列数
        mpiexec -n $NPROCS $PGMCASCI   &> H2O.caspt2.out
        mpiexec -n $NPROCS $PGMCASPT2O &>> H2O.caspt2.out
  ```
- active.inpは以下のような内容を記述してください
```in
8
6
142
2
3
2
3
0
156
0.0
C2
21

Input for CASCI and CASPT2

        read(5,'(I4)')ninact     # of inactive spinors
        read(5,'(I4)')nact       # of active spinors
        read(5,'(I4)')nsec       # of secondary spinors = nbas-ncore-nact-ninact
        read(5,'(I4)')nelec      # of active electrons in active space
        read(5,'(I4)')nroot      # of roots
        read(5,'(I4)')selectroot # which root do you want to obtain
        read(5,'(I4)')totsym     # total symmetry ex. 5 for Ag in C2h closed shell
        read(5,'(I4)')ncore      # of core orbital
        read(5,'(I4)')nbas       # of basis set
        read(5,'(E8.2)')eshift   # for real shift (if you don't write, it will be 0)
        read(5,'(A6)')ptgrp      # point group symmtery
        read(5,'(I4)')dirac_version # DIRAC version
```

  
## 開発者のかたへ

### 環境構築について

#### relqc01のマシンにおいては[野田](https://github.com/kohei-noda-qcrg)がcmake、gitおよびDIRAC(19.0,21.1,22.0)の環境を用意しています
#### 以下の記述を\$HOME/.bashrc に追記するとマシンログイン時に新しいバージョンのcmake,gitが使えます

\$HOME/.bashrc
```bash
module use --append "/home/noda/modulefiles" # Add Noda's modules
module purge            # deactivate all modules
module load cmake       # Load default cmake
module load git         # Load git
source "/home/noda/.config/git/.git-completion.bash" # Activate completions of the git command
##############################
# Git prompt
##############################
source "/home/noda/.config/git/git-prompt.sh" # This script allows you to see repository status in your prompt
export GIT_PS1_SHOWDIRTYSTATE=1 # cf. https://github.com/git/git/blob/e8005e4871f130c4e402ddca2032c111252f070a/contrib/completion/git-prompt.sh#L38-L42
export PS1='\[\033[01;32m\]\u@\h\[\033[01;34m\] \w\[\033[01;33m\]$(__git_ps1)\[\033[01;34m\] \$\[\033[00m\] ' # Change the prompt of your shell
```

#### 用意したDIRACの使い方
- \$HOME/.bashrcにmodule use --append "/home/noda/modulefiles"を記述します
- module load DIRAC/19.0 などと入力するとpam-diracコマンドが使えるようになります
  - DIRACのmoduleはDIRACを使うときだけ一時的にmodule loadすることをお勧めします
  - 従ってDIRACを実行する際は実行用のシェルスクリプト内でmodule loadすることを推奨します
  ```sh
  #!/bin/sh

  module load dirac/21.1 # Load DIRAC 21.1

  MOLECULE=H2O
  INPFILE=${MOLECULE}.inp
  MOLFILE=${MOLECULE}.xyz
  LOGFILE=${MOLECULE}.log
  NPROCS=8
  $PAM --mpi=$NPROCS --get="MRCONEE MDCIN*" '--keep_scratch' --mol=${MOLFILE} --inp=${INPFILE} --noarch &> $LOGFILE
  ```
- 一旦モジュールの読み込みを解除したいときは module unload 解除したいモジュールの名前 を実行します



### ビルドについて 
- デバッグ、リファクタリング時のビルドについて、--clean-first オプションを用いて前のビルド結果を消してから再ビルドすることをお勧めします

  ```sh
  cmake --build build --clean-first
  ```
- ビルドには[CMake](https://cmake.org/)を用います
  - ビルドの設定はCMakeLists.txtに書きます
  - 設定を追加したい場合は[公式ドキュメント](https://cmake.org/cmake/help/v3.7/)が正確でかなりわかりやすいので、"cmake やりたいこと"で検索してオプション名を見つけてから公式ドキュメントをみて追加することをお勧めします

### テストについて
- 今後[安全にコードを変更](https://ja.wikipedia.org/wiki/%E3%82%BD%E3%83%95%E3%83%88%E3%82%A6%E3%82%A7%E3%82%A2%E3%83%86%E3%82%B9%E3%83%88#%E5%A4%89%E6%9B%B4%E3%81%B8%E3%81%AE%E4%BF%A1%E9%A0%BC)できる(変更後と変更前の振る舞いが変わっていないことを確認する)ようにするために[テスト](https://ja.wikipedia.org/wiki/%E3%82%BD%E3%83%95%E3%83%88%E3%82%A6%E3%82%A7%E3%82%A2%E3%83%86%E3%82%B9%E3%83%88)を追加する予定です
- 本来は[単体テスト](https://ja.wikipedia.org/wiki/%E5%8D%98%E4%BD%93%E3%83%86%E3%82%B9%E3%83%88)を用いてプログラムの部品レベルでテストを書くべきですが、本プログラムはテストを前提として書かれておらず[密結合](https://e-words.jp/w/%E5%AF%86%E7%B5%90%E5%90%88.html)のため[単体テストが書きづらい](https://qiita.com/yutachaos/items/857472c7d3c65d3cf316#%E5%8D%98%E4%BD%93%E3%83%86%E3%82%B9%E3%83%88-1)です
- 当面は複数の分子系で、できるだけ違うタイプのインプットを用いて、最初に基準と定めたアウトプットから**自動的に**(ここがテストの良い点です)判定する形式にする予定です
  - 例えばCASPT2 energyが一定以上ずれていないかを判定するようにします
  - いわゆる[統合試験](https://ja.wikipedia.org/wiki/%E3%82%BD%E3%83%95%E3%83%88%E3%82%A6%E3%82%A7%E3%82%A2%E3%83%86%E3%82%B9%E3%83%88#%E7%B5%B1%E5%90%88%E8%A9%A6%E9%A8%93_(Integration_Testing))のみを行います
- ツールはFortranのテストツールは機能が貧弱なので、pythonのunittestかpytestを用いる予定です
  - DIRACもpythonを用いてテストを書いています
  - python側からビルドしたプログラムを実行し、アウトプットをリファレンス値と比較することで自動テストを実現します
