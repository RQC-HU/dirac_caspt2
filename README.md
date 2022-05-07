# DIRAC-CASPT2

#### [DIRAC](http://diracprogram.org/doku.php)の計算結果のうち1,2電子積分ファイルを用いて、CASCI/CASPT2法またはDMRG/CASPT2法で2次の多配置摂動計算を行います

## 目次
- [Requirements](https://github.com/kohei-noda-qcrg/dirac_caspt2#requirements)
- [How to Install](https://github.com/kohei-noda-qcrg/dirac_caspt2#how-to-install)
  - [ビルドオプション](https://github.com/kohei-noda-qcrg/dirac_caspt2#ビルドオプション)
  - [ビルド例](https://github.com/kohei-noda-qcrg/dirac_caspt2#ビルド例)
- [開発者のかたへ](https://github.com/kohei-noda-qcrg/dirac_caspt2#開発者のかたへ)
  - [環境構築について](https://github.com/kohei-noda-qcrg/dirac_caspt2#環境構築について)
  - [ビルドについて](https://github.com/kohei-noda-qcrg/dirac_caspt2#ビルドについて)
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
FC=ifort cmake --build build
cmake -B build --clean-first
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

## 開発者のかたへ

### 環境構築について

#### relqc01のマシンにおいては[野田](https://github.com/kohei-noda-qcrg)がcmakeおよびgitの環境を用意しています
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

### ビルドについて 
- デバッグ、リファクタリング時のビルドについて、--clean-first オプションを用いて前のビルド結果を消してから再ビルドすることをお勧めします

  ```sh
  cmake --build build --clean-first
  ```
- ビルドには[CMake](https://cmake.org/)を用います
  - ビルドの設定はCMakeLists.txtに書きます
  - 設定を追加したい場合は[公式ドキュメント](https://cmake.org/cmake/help/v3.7/)が正確でかなりわかりやすいので、"cmake やりたいこと"で検索してオプション名を見つけてから公式ドキュメントをみて追加することをお勧めします
