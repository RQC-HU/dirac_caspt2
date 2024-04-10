# DIRAC-CASPT2: A relativistic second order multi-configuration perturbation calculation program
[![DIRAC-CASPT2-CI-test](https://github.com/RQC-HU/dirac_caspt2/actions/workflows/ci.yml/badge.svg)](https://github.com/RQC-HU/dirac_caspt2/actions/workflows/ci.yml)  [日本語README](README.ja.md)

- This program performs the second order multi-configuration perturbation calculation using the IVO/CASCI/CASPT2 or IVO/RASCI/RASPT2 method with 1- and 2-electron integrals obtained from the [DIRAC](http://diracprogram.org/doku.php) calculation.

## Contribution

If you want to contribute to this project (bug report, feature request, pull request, etc.), please read the [CONTRIBUTING.md](CONTRIBUTING.md) file before you start contributing.

## Table of Contents

- [DIRAC-CASPT2: A relativistic second order multi-configuration perturbation calculation program](#dirac-caspt2-a-relativistic-second-order-multi-configuration-perturbation-calculation-program)
  - [Contribution](#contribution)
  - [Table of Contents](#table-of-contents)
  - [Download](#download)
  - [Prerequisites for build](#prerequisites-for-build)
  - [How to build](#how-to-build)
    - [Basic build](#basic-build)
    - [MPI Support](#mpi-support)
    - [Installation](#installation)
    - [CMake build options](#cmake-build-options)
  - [How to use](#how-to-use)
    - [User manual](#user-manual)
    - [Prerequisites for execution](#prerequisites-for-execution)
    - [Calculation](#calculation)
    - [Input file](#input-file)
    - [Input file specification](#input-file-specification)
  - [License](#license)

## Download

- Download the source code from GitHub.

```sh
git clone --depth=1 https://github.com/RQC-HU/dirac_caspt2.git
```

## Prerequisites for build

If you want to build this program, you need to have the following compilers, tools and libraries installed on your machine.

- [GNU Fortran](https://gcc.gnu.org/fortran/) or [Intel Fortran](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html) compiler (You can use the MPI compiler for parallel calculation)
- [CMake(version >= 3.14)](https://cmake.org/)
  - If CMake is not installed on your machine or the version is too old, please build CMake or use the pre-built CMake binary from [CMake Github](https://github.com/Kitware/CMake/releases).
  - If you use ifx or mpiifx as the Fortran compiler, CMake version >= 3.20.2 is required because [CMake supports ifx and mpiifx from version 3.20.2.](https://cmake.org/cmake/help/latest/release/3.20.html#id3:~:text=The%20Intel%20oneAPI%20Fortran%20compiler%20is%20now%20identified%20as%20IntelLLVM)
- [Intel MKL(Math Kernel Library)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html)
  - You need to configure the environment variable \$MKLROOT to link MKL.
    To verify that \$MKLROOT is configured, run the following command

    ```sh
    echo $MKLROOT
    ```

  - If you want to build with BLAS/LAPACK implementation other than MKL, you need to specify --no-mkl option and specify the link path of BLAS/LAPACK in --flags at this moment.

    Example

    ```sh
    ./setup --no-mkl --flags "Replace this by Your BLAS and LAPACK Library link path" --fc gfortran --build
    ```
- [Python(version >= 3.6)](https://www.python.org/)
  - The setup script(build script), the dcaspt2 script(program execution script) and pytest use Python.
  - If Python (version >= 3.6) is not installed on your machine and you don't have root privileges, it is recommended to install and setup Python with Python version management tool such as [pyenv](https://github.com/pyenv/pyenv)

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

    # Install Python (version >= 3.6)
    pyenv install 3.9.9

    # Set default Python version to the one installed with pyenv
    pyenv global 3.9.9
    ```

- [pytest](https://docs.pytest.org/en/latest/)
  - This program uses pytest to run tests.
  - If pytest is not installed on your machine, please install pytest with pip.

    ```sh
    pip install pytest
    ```

## How to build

- You can build this program with the setup script.
  - You can also build this program directly with CMake command, but we recommend you to use the setup script.
  - If you want to build directly with CMake, please see [CMake build options](#cmake-build-options).

### Basic build

- Change the directory to the source code directory.

```sh
cd /path/to/dirac_caspt2
```

- Build the program with the setup script (If you don't specify the compiler, use the compiler that CMake finds first.)

  ```sh
  ./setup --build
  ```

  - If you want to build with the specific compiler, please specify the compiler with the --fc option

    ```sh
    ./setup --build --fc ifort
    ```

  - You can check the options of setup script with the following command

    ```sh
    ./setup --help
    ```

  - You can also do a parallel build with multiple cores. You can specify the number of cores with the -j option

    ```sh
    ./setup --build -j 4
    ```

  - If you want to use OpenMP for thread parallel execution, please specify the --omp option

    ```sh
    ./setup --build --omp
    ```

- We recommend you to run tests after build

```sh
pytest --all
```

### MPI Support

- You can enable MPI support with the --mpi option (default compiler is mpiifort)

  ```sh
  ./setup --mpi --build
  ```

  - You can specify the compiler with the --fc option

    ```sh
    ./setup --mpi --fc mpif90 --build -j 4
    ```

  - You can also build with hybrid parallel execution with OpenMP

    ```sh
    ./setup --mpi --omp --fc mpiifort --build -j 4
    ```

- We recommend you to run tests after build

```sh
# pytest --all --mpi=<number of MPI processes>
pytest --all --mpi=4
```

### Installation

- You can install the program into the directory specified by the --prefix option with the following command


```sh
# Use CMake to install the program
cmake --install build
# or use make to install the program
make -C build install
```

### CMake build options

(This section is for advanced users.)

If you want to build directly with CMake, you can build with the following command

```sh
# DCMAKE_Fortran_COMPILER: Fortran compiler, (e.g.) ifort, gfortran, mpiifort
cmake -B build -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_BUILD_TYPE=Release -DOPENMP=ON && cmake --build build
pytest --all
```

You can use build options with cmake -DBUILDOPTION1=on -DBUILDOPTION2=off ,,,
The following custom CMake build options are currently supported

- MPI
  - Required if you want to use MPI. You need to build with MPI support to perform parallel execution.(default: OFF)
    (e.g.) Build with MPI support

    ```sh
    cmake -B build -DCMAKE_Fortran_COMPILER=mpiifort -DCMAKE_BUILD_TYPE=Release -DMPI=ON && cmake --build build
    ```

- OpenMP
  - Required if you want to use OpenMP for thread parallel execution.(default: OFF)
    (e.g.) Build with OpenMP support

    ```sh
    cmake -B build -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_BUILD_TYPE=Release -DOPENMP=ON && cmake --build build
    ```

- MKL
  - You need to disable this option if you want to use BLAS and LAPACK libraries other than MKL.(default: ON)
    (e.g.) Build without MKL support

    ```sh
    cmake -B build -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_BUILD_TYPE=Release -DMKL=OFF && cmake --build build
    ```

## How to use

### User manual

- This section explains the minimum usage of the program
- For more detailed usage of this program, please refer to the [GitHub Wiki](https://github.com/RQC-HU/dirac_caspt2/wiki/Users-wiki).

### Prerequisites for execution

- This program assumes that 1 and 2 electron integral files (MRCONEE, MDCINT, MDCINXXXX1...) are obtained by [DIRAC](http://diracprogram.org/) calculation
  - Please refer to the [**MOLTRA section of DIRAC manual](http://www.diracprogram.org/doc/master/manual/moltra.html) to obtain 1 and 2 electron integral files
  - 1 and 2 electron integral files must be in the same directory
- You need an [input file](#input-file) with any file name

### Calculation

- You can use this program with the dcaspt2 script
  - The dcaspt2 script is located in the bin directory or in the destination directory if installed by specifying with a --prefix option.

    ```sh
    # If you did not install with --prefix
    /path/to/dirac_caspt2/bin/dcaspt2 -i input_file
    # If you installed with --prefix
    $PREFIX/dcaspt2 -i input_file
    ```

    (e.g.)

    ```sh
    dcaspt2 -i h2.caspt2.inp
    ```

### Input file

- The input file is a text file with the following format

  - CASCI/CASPT2 input
    ```in
    ninact
    8
    nact
    6
    nsec
    142
    nelec
    2
    totsym
    3
    eshift
    0.0
    diracver
    21
    end
    ```

  - RASCI/RASPT2 input
    ```in
    ninact
    30
    nact ! sum of ras1, ras2 and ras3
    28
    nsec
    574
    nelec
    12
    totsym
    33
    eshift
    0.0
    diracver
    22
    ras1
    31..42
    2
    ras2
    43..48
    ras3
    49..58
    2
    end
    ```

- Please refer to the following for the meaning of each parameter and whether it is a required parameter or not

```in
Input for CASCI and CASPT2

[required parameters]
ninact      : the number of inactive spinors
nact        : the number of active spinors
nsec        : the number of secondary spinors
nelec       : the number of active electrons in active space
totsym      : total symmetry (ex. 5 for Ag in C2h closed shell)
diracver    : DIRAC version
end         : The identifier at the end of the input file

[required parameters (IVO)]
nocc        : The number of occupied MO (This option is for molecules without inversion center symmetry)
noccg       : The number of occupied MO (gerade)
noccu       : The number of occupied MO (ungerade)

[optional parameters]
nroot       : the number of roots (default: 10, max: 500, if the number of CASCI/RASCI configuration is less than nroot, nroot will be replaced by the number of CASCI/RASCI configuration.)
selectroot  : which root do you want to calculate RASPT2/CASPT2 energy? (default: 1, max: 500, the lowest root)
eshift      : for real shift (default: 0)
ras1        : RAS1 spinor list (row 1)and the maximum number of hole allowed in ras1(row 2)
ras2        : RAS2 spinor list
ras3        : RAS3 spinor list (row 1) and the maximum number of electrons in ras3(row2)
minholeras1 : The minimum number of hole in ras1 (default: 0)
calctype    : The type of calculation. only CASCI or DMRG are currently supported. (if you don't write, it will be CASCI(default))

[optional parameters (IVO)]
nhomo       : The number of HOMO-like spinors (default: 0)
nvcut       : The number of virtual cut MO (default: 0, This option is for molecules without inversion center symmetry)
nvcutg      : The number of virtual cut MO (default: 0, gerade)
nvcutu      : The number of virtual cut MO (default: 0, ungerade)
```

### Input file specification

- Reads 500 characters per a line
- Recognizes lines with an end parameter as the end of the input file
- Exits the program as invalid input if the required variable specifications are not filled.
- If you write \! or \#, the rest of the characters are recognized as comments.

  ```in
   nact ! The number of nact
   ↓
   nact
  ```

- Two or more consecutive dots are considered to be a range specification
  - You need to write the smaller number to the left of the dot and the larger number to the right

```in
  RAS2
  2..5
  ↓
  2,3,4,5
```

- Recognizes a , (semicolon) or half-width space as a numeric delimiter

```in
  RAS1
  1..4, 7   8 11..14
  ↓
  1,2,3,4,7,8,11,12,14
```

## License

- This program is mainly licensed under LGPL version 2.1 or later
  - Please see [LICENSE](LICENSE) for details.
- But some files are licensed under other licenses.
  - MIT license
    - [src/module_dict.f90](src/module_dict.f90)
    - [setup](setup)
    - [tools/dcaspt2_input](tools/dcaspt2_input)
