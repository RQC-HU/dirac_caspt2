# DIRAC-CASPT2

- [DIRAC](http://diracprogram.org/doku.php)の計算結果のうち1,2電子積分ファイルを用いて、CASCI/CASPT2法またはDMRG/CASPT2法で2次の多配置摂動計算を行います

## お知らせ

- [開発者の方はGithub Wikiを参考に開発を行ってください](https://github.com/kohei-noda-qcrg/dirac_caspt2/wiki/developers-wiki)

## 目次

- [DIRAC-CASPT2](#dirac-caspt2)
  - [お知らせ](#お知らせ)
  - [目次](#目次)
  - [Requirements](#requirements)
  - [How to build](#how-to-build)
    - [Basic build](#basic-build)
    - [MPI Support](#mpi-support)
    - [Install](#install)
    - [CMakeビルドオプション](#cmakeビルドオプション)
  - [How to use](#how-to-use)
    - [Prerequisites](#prerequisites)
    - [Calculation](#calculation)
    - [input file](#input-file)
    - [インプットファイルの仕様](#インプットファイルの仕様)

## Requirements

以下のコンパイラおよびツール、ライブラリと依存性があり、ビルドを行う計算機でこれらがセットアップされている必要があります

- [GNU Fortran](https://gcc.gnu.org/fortran/) or [Intel Fortran](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html) compiler (並列計算をするために並列コンパイラを使うこともできます)
- [CMake(version ≧ 3.14)](https://cmake.org/)
  - CMakeが計算機に入っていないか、バージョンが古い場合[CMakeのGithub](https://github.com/Kitware/CMake/releases)からビルドするもしくはビルド済みのファイルを解凍して使用してください
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
  - Python (version ≧ 3.6)がインストールされておらず、かつルート権限がない場合[pyenv](https://github.com/pyenv/pyenv)などのPythonバージョンマネジメントツールを使用して非ルートユーザーでPythonをインストール、セットアップすることをおすすめします

    (e.g.) pyenv setup instruction for Bash users

    ```bash
    # Download pyenv
    git clone https://github.com/pyenv/pyenv.git ~/.pyenv

    # Write the environmental variable and setup script for pyenv to the ~/.bashrc file
    echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bashrc
    echo 'command -v pyenv >/dev/null || export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bashrc
    echo 'eval "$(pyenv init -)"' >> ~/.bashrc

    # Reload ~/.bashrc
    source ~/.bashrc

    # Install Python (version ≧ 3.6)
    pyenv install 3.9.9

    # Set default Python version to the one installed with pyenv
    pyenv global 3.9.9
    ```

- [pytest](https://docs.pytest.org/)
  - テストを実行するために使用します
  - Python (version ≧ 3.6)をインストールしていれば以下のコマンドで入手できます

  ```sh
  python -m pip install pytest
  ```

## How to build

- このプログラムはCMakeを使用してビルドを行います
  - CMakeコマンドを直接使用してビルドすることもできますが、setupスクリプトを使用することをおすすめします
  - CMakeを直接使用してビルドしたい場合は、[CMakeビルドオプション](#cmakeビルドオプション)を参照してください

### Basic build

- GitHubからソースコードをダウンロードします(初回のみ)

```sh
git clone --depth=1 https://github.com/kohei-noda-qcrg/dirac_caspt2.git
```

- ソースコードのディレクトリに移動します

```sh
# Change directory to the source code directory
# ( cd /path/to/dirac_caspt2 )
cd dirac_caspt2
```

- セットアップスクリプトを実行します。(--buildオプションをつけるとビルドまで行います。--fcオプションをつけてコンパイラを明示的に指定することを推奨します)

  ```sh
  ./setup --build --fc=ifort
  ```

  - セットアップスクリプトのオプションについては以下のコマンドで確認できます

    ```sh
    ./setup --help
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

## How to use

### Prerequisites

- [DIRAC](http://diracprogram.org/)の計算で1,2電子積分ファイル(MRCONEE, MDCINT, MDCINXXXX1...)が得られていることを前提としています
  - 1,2電子積分ファイルを得るには[DIRACの**MOLTRAの項](http://www.diracprogram.org/doc/master/manual/moltra.html)を参照してください
  - 1,2電子積分ファイルは同一のディレクトリ上に存在する必要があります
- 任意のファイル名の[インプットファイル](#input-file)が必要です

### Calculation

- ビルド後に作られるbinディレクトリ直下またはprefixを指定した場合はインストール先のディレクトリ直下のdcaspt2スクリプトを用いて計算を行います
  - dcaspt2スクリプトで使用可能なオプションはdcaspt2 -hで確認できます
  - 例えば以下のように使用します

  ```sh
  dcaspt2 -i h2.caspt2.inp
  ```


### input file

- インプットファイルには以下のような内容を記述してください

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

各パラメータの意味と必須パラメータかどうかについては以下を参照してください

```in
Input for CASCI and CASPT2

[required parameters]
ninact      : the number of inactive spinors
nact        : the number of active spinors
nsec        : the number of secondary spinors = nbas-ncore-nact-ninact
nelec       : the number of active electrons in active space
nroot       : the number of roots
selectroot  : which root do you want to obtain
totsym      : total symmetry (ex. 5 for Ag in C2h closed shell)
ncore       : the number of core orbital
nbas        : the number of basis set
diracver    : DIRAC version
end         : The identifier at the end of the input file

[optional parameters]
eshift      : for real shift (if you don't write, it will be 0)
ras1        : RAS1 spinor list (row 1)and the maximum number of hole allowed in ras1(row 2)
ras2        : RAS2 spinor list
ras3        : RAS3 spinor list (row 1) and the maximum number of electrons in ras3(row2)
minholeras1 : The minimum number of hole in ras1 (If you don't write, it will be 0)
calctype    : The type of calculation. only CASCI or DMRG are currently supported. (if you don't write, it will be CASCI(default))
```

### インプットファイルの仕様

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
  1..4, 7   8 11..14
  ↓
  1,2,3,4,7,8,11,12,14
```
