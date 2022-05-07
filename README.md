# DIRAC-CASPT2

- DIRAC19or21の1,2電子積分ファイルを用いて、CASCI/CASPT2法またはDMRG/CASPT2法で2次の多配置摂動計算を行います

### Requirements

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

- 現状GNU Fortranはビルドは成功しますが実行時エラーが発生する可能性があるため非推奨です
- したがってFC=ifort もしくは FC=mpiifort もしくは FC=mpifort (OpenMPI,ifort) を使用することを推奨します

### ビルドオプション

現時点でサポートしているビルドオプションは以下のとおりです

ビルドオプションはcmake -DBUILDOPTION1=on -DBUILDOPTION2=off ,,,のように使います

- MPI
    - MPIを使用するなら必須です.プリプロセッサを追加します
    
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
