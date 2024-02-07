#!/bin/bash
# Fetch the Orekit data file and write/run shell script to define paths
#
# Without current data file, run this script with (optionally) the
# directory where orekit-data.zip should be placed.
# To locate in this directory
#  . get-orekit-data.bash
# To locate elsewhere, for example, the parent directory
#  . get-orekit-data.bash ..
#
# For future Orekit usage with the same Orekit data, run the shell
# script env.sh:
#  . env.sh
# This sets up the correct environment variables without a new download

python -c "import orekit.pyhelpers
orekit.pyhelpers.download_orekit_data_curdir()"

if [ "$#" -eq 1 ] && [ -d "$1" ]; then
    dest=$(realpath $1)
    mv orekit-data.zip $dest
else
    dest=$PWD
fi

echo "export PYTHONPATH=$(dirname $(realpath ${BASH_SOURCE[0]}))" > env.sh
echo "export OREKITDATA=$dest/orekit-data.zip" >> env.sh
. env.sh
