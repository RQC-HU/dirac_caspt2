#!/usr/bin/env bash
abspath="$(cd "$( dirname "$0" )" && pwd )"
echo "runtest abspath $abspath"
# Did you install python on your machine?
if ! (type "python" >/dev/null 2>&1); then
    echo "=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+"
    printf "ERROR:\tPython is not installed."
    echo "You need to install python to run tests."
    echo "You can easily install python using"
    echo "pyenv : https://github.com/pyenv/pyenv"
    echo " even if you are not root user."
    echo "=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+"
    echo "Exit."
    exit
fi

# Built CASCI program already?
r4dcasci="$abspath/bin/r4dcascicoexe"
if [ ! -f  "$r4dcasci" ]; then
    echo "=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+"
    printf "Error:\t$r4dcasci is missing."
    echo "You must build rel-caspt2 program before running tests."
    echo "Read $abspath/README.md and try to build this program."
    echo "=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+"
    echo "Exit."
    exit
fi

# Built CASPT2 program already?
r4dcaspt2="$abspath/bin/r4dcaspt2ocoexe"
if [ ! -f  "$r4dcaspt2" ]; then
    echo "=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+"
    printf "Error:\t$r4dcaspt2 is missing."
    echo "You must build rel-caspt2 program before running tests."
    echo "Read $abspath/README.md and try to build this program."
    echo "=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+"
    echo "Exit."
    exit
fi

# Did you install pytest on your machine?
if type "pytest" >/dev/null 2>&1; then
    cd "$abspath/test" || exit
    pytest --tb=short -vv # Run tests!!
else
    echo "=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+"
    printf "Error:\tpytest is not installed."
    echo "You must install pytest to run tests."
    echo "python3 -m pip install pytest"
    echo "=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+"
    echo "Exit."
fi