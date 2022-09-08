# DIRAC-CASPT2

#### [DIRAC](http://diracprogram.org/doku.php)の計算結果のうち1,2電子積分ファイルを用いて、CASCI/CASPT2法またはDMRG/CASPT2法で2次の多配置摂動計算を行います

## 目次

- [Requirements](https://github.com/kohei-noda-qcrg/dirac_caspt2#requirements)
- [How to Install](https://github.com/kohei-noda-qcrg/dirac_caspt2#how-to-install)
  - [ソフトウェアのテスト](https://github.com/kohei-noda-qcrg/dirac_caspt2#ソフトウェアのテスト)
  - [ビルドオプション](https://github.com/kohei-noda-qcrg/dirac_caspt2#ビルドオプション)
  - [ビルド例](https://github.com/kohei-noda-qcrg/dirac_caspt2#ビルド例)
- [How to use](https://github.com/kohei-noda-qcrg/dirac_caspt2#how-to-use)
- [開発者のかたへ](https://github.com/kohei-noda-qcrg/dirac_caspt2#開発者のかたへ)

## Requirements

以下のコンパイラおよびツール、ライブラリと依存性があり、ビルドを行う計算機でこれらがセットアップされている必要があります

- [GNU Fortran](https://gcc.gnu.org/fortran/) or [Intel Fortran](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html) compiler (並列計算をするために並列コンパイラを使うこともできます)
- [CMake](https://cmake.org/)(version ≧ 3.7 が必要です)
    - cmakeが計算機に入っていないか、バージョンが古い場合[CMakeのGithub](https://github.com/Kitware/CMake/releases)からビルドするもしくはビルド済みのファイルを解凍して使用してください
- [Intel MKL(Math Kernel Library)](https://www.intel.com/content/www/us/en/develop/documentation/get-started-with-mkl-for-dpcpp/top.html)
  - MKLをリンクするため環境変数\$MKLROOTが設定されている必要があります
    \$MKLROOTが設定されているか確認するには、使用する計算機にログインして以下のコマンドを実行してMKLにパスが通っているかを確認してください

    ```sh
    echo $MKLROOT
    ```
  - 現時点ではMKLのBlas,Lapack以外のBlas,Lapackの実装を用いてビルドする場合、-DMKL=offオプションを指定し、かつLDFLAGSを手動設定する必要があります
  - また、MKLのBlas,Lapack以外での動作は現在保障しておりませんのでご了承ください

    ビルド例
    ```sh
    mkdir build
    cd build
    LDFLAGS="Replace this by Your BLAS and LAPACK Library link path" FC=gfortran cmake -DMKL=off ..
    make
    ```


- [Python(version ≧ 3.6)](https://www.python.org/)
  - テストを実行するために使用します
  - Python (version ≧ 3.6)がインストールされておらず、かつルート権限がない場合[pyenv](https://github.com/pyenv/pyenv)などのPythonバージョンマネジメントツールを使用して非ルートユーザーでPythonをインストール、セットアップすることをおすすめします
- [pytest](https://docs.pytest.org/)
  - テストを実行するために使用します
  - python (version ≧ 3.6)をインストールしていれば以下のコマンドで入手できます
  ```sh
  python -m pip install pytest
  ```
## How to Install

以下のコマンドでmainブランチのソースコードをビルドできます

```sh
git clone https://github.com/kohei-noda-qcrg/dirac_caspt2
cd dirac_caspt2
mkdir -p build && cd build
FC=ifort cmake .. --clean-first
make
```

- CMake version ≧ 3.13 を使っているなら以下のようなコマンドでもビルドができます

```sh
git clone https://github.com/kohei-noda-qcrg/dirac_caspt2
cd dirac_caspt2
FC=ifort cmake -B build
cmake --build build --clean-first
```

- 現在Intel Fortranであれば並列ビルドが可能です。並列ビルドは-j並列数のオプションを付ければ実行可能です

```sh
git clone https://github.com/kohei-noda-qcrg/dirac_caspt2
cd dirac_caspt2
FC=ifort cmake -B build
cmake --build build -j4 --clean-first
```

### ソフトウェアのテスト

ビルド後はテストを行うことを推奨します
テストを行うには[Python(version ≧ 3.6)](https://www.python.org/)と[pytest](https://docs.pytest.org/)が必要です
testディレクトリより上位のディレクトリでpytestコマンドを実行することでテストが実行されます

```sh
pytest
```

並列コンパイラでビルドオプション-DMPI=onをつけてMPI並列用のビルドを行った場合
pytestコマンドに--paralles=並列数を付け加え、並列用テストを行うことを推奨します

```sh
pytest --parallel=4
```


### ビルドオプション

現時点でサポートしているビルドオプションは以下のとおりです

ビルドオプションはcmake -DBUILDOPTION1=on -DBUILDOPTION2=off ,,,のように使います

- MPI
    - MPIを使用するなら必須です.マルチプロセス対応ビルドのためのプリプロセッサの設定を行います(default:OFF)

        (例)

        ```sh
        mkdir -p build && cd build
        FC=mpiifort cmake -DMPI=on ..
        make
        ```

- OPENMP

    - OpenMPを使用するなら必須です.OpenMP用のビルドオプションを追加します(default:OFF)

        (例)

        ```sh
        mkdir -p build && cd build
        FC=ifort cmake -DOPENMP=on ..
        make
        ```
- MKL

    - MKLを使わないときはこのビルドオプションをOFFにする必要があります.デフォルトがONなので指定しなければMKLを使う前提でビルドを行います(default:ON)

        (例)

        ```sh
        mkdir -p build && cd build
        LDFLAGS="/your/blas/link/path /your/lapack/link/path" FC=ifort cmake -DMKL=off ..
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
ninact
8
nact
6
nsec
142
nelec
2
nroot
3
selectroot
2
totsym
3
ncore
0
nbas
156
eshift
0.0
ptgrp
C2
diracver
21
ras1
1..4,9,10
1
ras2
5 6
ras3
20..30
3
end
```

各パラメータの意味と必須パラメータかどうかについては以下を参照してください(requiredとあるものは必須パラメータです

```in
Input for CASCI and CASPT2

ninact      : the number of inactive spinors (required)
nact        : the number of active spinors (required)
nsec        : the number of secondary spinors = nbas-ncore-nact-ninact  (required)
nelec       : the number of active electrons in active space (required)
nroot       : the number of roots (required)
selectroot  : which root do you want to obtain (required)
totsym      : total symmetry ex. 5 for Ag in C2h closed shell (required)
ncore       : the number of core orbital (required)
nbas        : the number of basis set (required)
eshift      : for real shift (if you don't write, it will be 0)
ptgrp       : point group symmtery (required)
diracver    : DIRAC version (required)
ras1        : RAS1 spinor list (row 1)and the maximum number of hole allowed in ras1(row 2)
ras2        : RAS2 spinor list
ras3        : RAS3 spinor list (row 1) and the maximum number of electrons in ras3(row2)
minholeras1 : The minimum number of hole in ras1 (If you don't write, it will be 0)
calctype    : The type of calculation. only CASCI or DMRG are currently supported. (if you don't write, it will be CASCI(default))
end         : The identifier at the end of active.inp (required)
```

### active.inpの仕様

- 1行あたり100文字を読み取ります
- endがある行をインプットの終わりと認識します
- end及びrequiredな変数についての指定がないと不正なインプットとしてプログラムを終了します
- \!か\#を書くとそれ以降の文字はコメントと認識します

```in
 nact ! The number of nact
 ↓
 nact
```

- RASについて.(ドット)が2つ以上連続していると範囲指定をしているとみなします(左に小さい数値、右に大きい数値を書く必要があります)

```in
  1..4
  ↓
  1,2,3,4
```

- RASについて,(セミコロン)もしくは半角スペースを数値の区切りであると認識します

```in
  1..4, 5   7 10..13
  ↓
  1,2,3,4,5,7,10,11,12,13
```

## 開発者のかたへ
- [開発者Wikiを参考にして開発を行ってください](https://github.com/kohei-noda-qcrg/dirac_caspt2/wiki/developers-wiki)
