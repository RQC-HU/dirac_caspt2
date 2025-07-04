# DIRAC-CASPT2: 相対論的多配置2次摂動論プログラム


[![DIRAC-CASPT2-CI-test](https://github.com/RQC-HU/dirac_caspt2/actions/workflows/ci.yml/badge.svg)](https://github.com/RQC-HU/dirac_caspt2/actions/workflows/ci.yml) [English README](README.md)

- [DIRAC](http://diracprogram.org/doku.php)の計算結果のうち1,2電子積分ファイルを用いて、CASCI/CASPT2法またはRASCI/RASPT2法で2次の多配置摂動計算を行います

## Contribution

このプロジェクトに貢献(バグレポート、機能追加など)する方法については[CONTRIBUTING.md](CONTRIBUTING.md)を参照してください

## 目次

- [DIRAC-CASPT2: 相対論的多配置2次摂動論プログラム](#dirac-caspt2-相対論的多配置2次摂動論プログラム)
  - [Contribution](#contribution)
  - [目次](#目次)
  - [ダウンロード](#ダウンロード)
  - [Prerequisites for build](#prerequisites-for-build)
  - [How to build](#how-to-build)
    - [Basic build](#basic-build)
    - [MPI Support](#mpi-support)
    - [Install](#install)
    - [CMakeビルドオプション](#cmakeビルドオプション)
  - [How to use](#how-to-use)
    - [User manual](#user-manual)
    - [Prerequisites for execution](#prerequisites-for-execution)
    - [Calculation](#calculation)
    - [input file](#input-file)
    - [インプットファイルの仕様](#インプットファイルの仕様)
  - [License](#license)
  - [Citation](#citation)

## ダウンロード

- GitHubからソースコードをダウンロードします

```sh
git clone --depth=1 https://github.com/RQC-HU/dirac_caspt2.git
```

## Prerequisites for build

以下のコンパイラおよびツール、ライブラリと依存性があり、ビルドを行う計算機でこれらがセットアップされている必要があります

- [GNU Fortran](https://gcc.gnu.org/fortran/) or [Intel Fortran](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html) compiler (並列計算をするために並列コンパイラを使うこともできます)
- [CMake(version ≧ 3.14)](https://cmake.org/)
  - CMakeが計算機に入っていないか、バージョンが古い場合[CMakeのGithub](https://github.com/Kitware/CMake/releases)からビルドするもしくはビルド済みのファイルを解凍して使用してください
  - ifx又はmpiifxをFortranコンパイラとして使用する場合、[CMakeが3.20.2からifx,mpiifxをサポートしたので](https://cmake.org/cmake/help/latest/release/3.20.html#id3:~:text=The%20Intel%20oneAPI%20Fortran%20compiler%20is%20now%20identified%20as%20IntelLLVM)バージョン3.20.2以上を使用してください
- [Intel MKL(Math Kernel Library)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html)
  - MKLをリンクするため環境変数\$MKLROOTが設定されている必要があります
    \$MKLROOTが設定されているか確認するには、以下のコマンドを実行して環境変数\$MKLROOTが設定されているか確認してください

    ```sh
    echo $MKLROOT
    ```

  - 現時点ではMKLのBlas,Lapack以外のBlas,Lapackの実装を用いてビルドする場合、--no-mklオプションを指定し、かつBLAS/LAPACKのリンクを--flagsに指定する必要があります
  - また、MKLのBlas,Lapack以外での動作は現在保障しておりませんのでご了承ください

    ビルド例

    ```sh
    ./setup --no-mkl --flags "Replace this by Your BLAS and LAPACK Library link path" --fc gfortran --build
    ```

- [Python(version ≧ 3.6)](https://www.python.org/)
  - setupスクリプト(ビルド用スクリプト),dcaspt2スクリプト(プログラム実行用スクリプト)およびテストを実行するのに使用します
  - Python (version ≧ 3.6)がインストールされておらず、かつルート権限がない場合[pyenv](https://github.com/pyenv/pyenv?tab=readme-ov-file#installation)などのPythonバージョンマネジメントツールを使用して非ルートユーザーでPythonをインストール、セットアップすることをおすすめします
- [pytest](https://docs.pytest.org/)
  - テストを実行するために使用します
  - Python (version ≧ 3.6)をインストールしていれば以下のコマンドで入手できます

  ```sh
  python -m pip install pytest
  ```

## How to build

- このプログラムはsetupスクリプトを使用してビルドできます
  - CMakeコマンドを直接使用してビルドすることもできますが、setupスクリプトを使用することをおすすめします
  - CMakeを直接使用してビルドしたい場合は、[CMakeビルドオプション](#cmakeビルドオプション)を参照してください
- デフォルトでは、ビルドした結果のバイナリとdcaspt2スクリプトはbuildディレクトリ直下に配置されます

### Basic build

- ソースコードのディレクトリに移動します

```sh
cd /path/to/dirac_caspt2
```

- セットアップスクリプトを実行してビルドします(コンパイラが指定されていない場合、CMakeが最初に見つけたコンパイラが使われます)


  ```sh
  ./setup --build
  ```

  - コンパイラを明示的に指定する場合は--fcオプションを使用します

    ```sh
    ./setup --build --fc ifort
    ```

  - セットアップスクリプトのオプションについては以下のコマンドで確認できます

    ```sh
    ./setup --help
    ```

  - 整数値のデフォルトサイズを64bitでビルドするには--int64オプションを使用します

    ```sh
    ./setup --build --int64
    ```

  - 複数コアを用いた並列ビルドも可能です。並列ビルドは-j 並列数のオプションを付ければ実行できます

    ```sh
    ./setup --build -j 4
    ```

  - スレッド並列実行用にOpenMPを使用する場合は--ompオプションを付けてください

    ```sh
    ./setup --fc=ifort --omp --build
    ```

- ビルドが完了したら問題なくビルドできたか確かめるため、テストを実行することを推奨します

```sh
pytest --all
```

### MPI Support

- プログラムをプロセス並列実行するためにMPIを有効にする場合、--mpiオプションを付けてビルドします(デフォルトで使用するコンパイラはmpiifortです)

  ```sh
  ./setup --mpi --build
  ```

  - コンパイラを指定する場合は--fcオプションを使用します

    ```sh
    ./setup --mpi --fc mpif90 --build -j 4
    ```

  - OpenMPとのハイブリッド並列のビルドも可能です

    ```sh
    ./setup --mpi --omp --fc mpiifort --build -j 4
    ```

- ビルドが完了したら問題なくビルドできたか確かめるため、テストを実行することを推奨します

```sh
# pytest --all --mpi=<number of MPI processes>
pytest --all --mpi=4
```

### Install

- 以下のいずれかのコマンドで--prefixで指定したインストール先にプログラムをインストールできます

```sh
# Use CMake to install the program
cmake --install build
# or use make to install the program
make -C build install
```

### CMakeビルドオプション

(This section is for advanced users.)

CMakeを直接使ってビルドする場合以下のようなコマンドを実行するとビルドできます

```sh
# DCMAKE_Fortran_COMPILER: Fortran compiler, (e.g.) ifort, gfortran, mpiifort
cmake -B build -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_BUILD_TYPE=Release -DOPENMP=ON && cmake --build build
pytest --all
```

ビルドオプションはcmake -DBUILDOPTION1=on -DBUILDOPTION2=off ,,,のように使います

現時点でサポートしているカスタムCMakeビルドオプションは以下のとおりです

- MPI
  - MPIを使用するなら必須です.マルチプロセス対応ビルドのためのプリプロセッサの設定を行います(default:OFF)

      (例)

      ```sh
      cmake -DCMAKE_Fortran_COMPILER=mpiifort -DMPI=on -B build && cmake --build build
      ```

- OPENMP

  - OpenMPを使用するなら必須です.OpenMP用のビルドオプションを追加します(default:OFF)

      (例)

      ```sh
      cmake -DCMAKE_Fortran_COMPILER=ifort -DOPENMP=on -B build && cmake --build build
      ```

- MKL

  - MKLを使わないときはこのビルドオプションをOFFにする必要があります.(default:ON)

      (例)

      ```sh
      LDFLAGS="/your/blas/link/path /your/lapack/link/path" cmake -DCMAKE_Fortran_COMPILER=ifort -DMKL=off -B build && cmake --build build
      ```

- INT64

  - デフォルトのintegerサイズを64bitにします.(default: OFF)

    (例)

    ```sh
    cmake -B build -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_BUILD_TYPE=Release -DINT64=ON && cmake --build build
    ```


## How to use

### User manual

- 以下のセクションでは、最低限の使用方法について説明します
- 本プログラムの使用方法の詳細については[GitHub Wiki](https://github.com/RQC-HU/dirac_caspt2/wiki/%E3%83%A6%E3%83%BC%E3%82%B6%E5%90%91%E3%81%91wiki)を参照してください

### Prerequisites for execution

- [DIRAC](http://diracprogram.org/)の計算で1,2電子積分ファイル(MRCONEE, MDCINT, MDCINXXXX1...)が得られていることを前提としています
  - 1,2電子積分ファイルを得るには[DIRACの**MOLTRAの項](http://www.diracprogram.org/doc/master/manual/moltra.html)を参照してください
  - 1,2電子積分ファイルは同一のディレクトリ上に存在する必要があります
- 任意のファイル名の[インプットファイル](#input-file)が必要です

### Calculation

- ビルド後に作られるbuildディレクトリ直下またはprefixを指定した場合はインストール先のディレクトリ直下のdcaspt2スクリプトを用いて計算を行います
  - dcaspt2スクリプトで使用可能なオプションはdcaspt2 -hで確認できます
  - 例えば以下のように使用します

  ```sh
  dcaspt2 -i h2.caspt2.inp
  ```

### input file

- インプットファイルには以下のような内容を記述してください
  - CASCI/CASPT2 input
    ```in
    .ninact
    8
    .nact
    6
    .nsec
    142
    .nelec
    2
    .caspt2_ciroots
    3 1 2 3  ! calculate CASCI energy for total symmetry 3 (row 1) and CASPT2 energies for the 1st, 2nd, and 3rd roots of total symmetry 3 (row 2, 3 and 4)
    4 1      ! calculate CASCI energy for total symmetry 4 (row 1) and CASPT2 energies for the 1st root of total symmetry 4 (row 2)
    .eshift
    0.0
    .diracver
    21
    .subprograms
    CASCI
    CASPT2
    .end
    ```

  - RASCI/RASPT2 input
    ```in
    .ninact
    30
    .nact ! sum of ras1, ras2 and ras3
    28
    .nsec
    574
    .nelec
    12
    .caspt2_ciroots
    33 1
    .eshift
    0.0
    .diracver
    22
    .ras1
    31..42
    2
    .ras2
    43..48
    .ras3
    49..58
    2
    .subprograms
    CASCI
    CASPT2
    .end
    ```

各パラメータの意味と必須パラメータかどうかについては以下を参照してください

```in
Input for CASCI and CASPT2

[required parameters]
.ninact         : the number of inactive spinors
.nact           : the number of active spinors
.nsec           : the number of secondary spinors
.totsym         : total symmetry (ex. 5 for Ag in C2h closed shell)
.diracver       : DIRAC version
.end            : The identifier at the end of the input file

[required parameters (IVO)]
.nocc           : The number of occupied MO (This option is for molecules without inversion center symmetry)
.noccg          : The number of occupied spinors (gerade)
.noccu          : The number of occupied spinors (ungerade)

[required parameters (CASCI and CASPT2)]
.caspt2_ciroots : Multiple line input. total symmery to calculate CASCI/CASPT2 energy (ex. 5 for Ag in C2h closed shell) (row 1),
                  number of roots that you want to calculate CASPT2 energy (row 2 and later)

[optional parameters]
.eshift         : for real shift (default: 0)
.ras1           : RAS1 spinor list (row 1)and the maximum number of hole allowed in ras1(row 2)
.ras2           : RAS2 spinor list
.ras3           : RAS3 spinor list (row 1) and the maximum number of electrons in ras3(row2)
.minholeras1    : The minimum number of hole in ras1 (default: 0)
.minelecras3    : The minimum number of electrons in ras3 (default: 0)
.scheme         : MOLTRA SCHEME, if you explicitly set the non-default .SCHEME value in **MOLTRA, you must set the same value for this option. (ref .SCHEME: https://diracprogram.org/doc/master/manual/moltra.html#scheme)
.debugprint     : This keyword invokes printing of additional information in the output file
.restart        : Restart calculation from the previous calculation. You need to generate the caspt2_restart file by running gen_dcaspt2_restart [previous_calclation_output] and put it in the same directory as the input file. (default: .false.)
.countndet      : Count and print the number of determinants for each total symmetry of this input file, and skip the any other calculations. If you set this option, you don't need to set .subprograms parameter. (default: false)

[optional parameters (IVO)]
.nhomo          : The number of HOMO-like spinors (default: 0)
.nvcut          : The number of virtual cut MO (default: 0, This option is for molecules without inversion center symmetry)
.nvcutg         : The number of virtual cut MO (default: 0, gerade)
.nvcutu         : The number of virtual cut MO (default: 0, ungerade)
```

### インプットファイルの仕様

- 1行あたり500文字を読み取ります
- endがある行をインプットの終わりと認識します
- end及びrequiredな変数についての指定がないと不正なインプットとしてプログラムを終了します
- \!か\#を書くとそれ以降の文字はコメントと認識します

```in
 .nact ! The number of nact
 ↓
 .nact
```

- .(ドット)が2つ以上連続していると範囲指定をしているとみなします(左に小さい数値、右に大きい数値を書く必要があります)

```in
  1..4
  ↓
  1,2,3,4
```

- ,(セミコロン)もしくは半角スペースを数値の区切りであると認識します

```in
  1..4, 7   8 11..14
  ↓
  1,2,3,4,7,8,11,12,14
```

## License

- LGPL version 2.1 or later
  - 詳細については [LICENSE](LICENSE) ファイルを参照ください
- ただし以下のファイルは別のライセンスの元で配布されています
  - MIT license
    - [src/module_dict.F90](src/module_dict.F90)
    - [setup](setup)
    - [tools/dcaspt2_input](tools/dcaspt2_input)

## Citation

- 本プログラムを使って得られたデータを公開する場合は以下の論文を引用してください
  - Y. Masuda, K. Noda, S. Iwamuro, M. Hada, N. Nakatani, M. Abe. Relativistic CASPT2/RASPT2 Program along with DIRAC software. J. Chem. Theory Comput. 2025, 21, 3, 1249–1258; [https://doi.org/10.1021/acs.jctc.4c01589](https://doi.org/10.1021/acs.jctc.4c01589)
  - M. Abe, T. Nakajima, K. Hirao. The relativistic complete active-space second-order perturbation theory with the four-component Dirac Hamiltonian. J. Chem. Phys. 2006, 125, 234110; [https://doi.org/10.1063/1.2404666](https://doi.org/10.1063/1.2404666)
  - S. Huzinaga, and C. Arnau. Virtual Orbitals in Hartree-Fock Theory. Phys. Rev. A. 1970, 1, 1285-1288; [https://doi.org/10.1103/PhysRevA.1.1285](https://doi.org/10.1103/PhysRevA.1.1285)
  - K. Andersson, P. A. Malmqvist, B. O. Roos, A. J. Sadlej, K. Wolinski. Second-order perturbation theory with a CASSCF reference function. J. Phys. Chem. 1990, 94, 14, 5483–5488; [https://doi.org/10.1021/j100377a012](https://doi.org/10.1021/j100377a012)
