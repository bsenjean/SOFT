# SOFT
Site Occupation Functional Theory on the Hubbard model

# Installation

export SOFT_DIR=`pwd`

cd bin

make

cd ..

pip install -e .

export "PATH=$PATH:$SOFT_DIR/bin/"

# execution

cd examples

python3 SOFT_BALDA.py
