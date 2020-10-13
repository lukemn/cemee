#!/bin/bash
#
# Script to run Choreography, given the path to a folder containing the data from an experiment
# This script exports a 'phenotype_names.txt' file, as a proper way to keep track of the 
# phenotypes that were considered (hence avoiding data format problems as in the data that 
# Bruno processed)
#
# *** Note that the path to the Choreography jar file must be set (see variable 
# 'CHORE_JAR' below)
#
# Modified to include the "BodyOrientation" phenotype, for backwards compatibility
# (argument 'orient')
#
# Created on June, 2016, based on code from script 'ChoreographyHulk.sh'
# in Bruno's home folder (in Hulk: /home/bafonso/Scripts)
#

# Exit immediately if a simple command exits with a nonzero exit value
set -e

# Path to the Choreography jar
CHORE_JAR="/users/gev/mallard/Software/MWT/MWT_1.3.0_r1035/analysis/Chore.jar"

if [ $# -ne 1 ]; then
    echo Usage: ./run_choreograph.sh INPUT_FOLDER_PATH
    exit -1
fi
input_dir_path="$1"

#
# Exit on error, after displaying a message
#
function exitOnError {
    echo "ERROR: $1"
    exit -1
}

# Make sure the paths are absolute
# Keep in mind that using readlink in Mac OS X is not straightforward
if [ "$(uname)" == "Darwin" ]; then # Mac OX X
    input_dir_path=$(greadlink -e "${input_dir_path}")
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    input_dir_path=$(readlink -e "${input_dir_path}")
elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW32_NT" ]; then
    exitOnError "TODO"
fi

if [[ ! -d "${input_dir_path}" ]]; then
    exitOnError "The input must be a directory"
fi

phenotype_str="id,persistence,area,speed,angular,length,width,aspect,midline,morphwidth,bias,pathlen,curve,loc_x,loc_y,orient,area:jitter,midline:jitter,morphwidth:jitter,loc_x:jitter,loc_y:jitter,kink"

java -jar -Xmx6G -Xms6G "${CHORE_JAR}" \
    -S --shadowless -q --plugin Reoutline --plugin Respine \
    -o "${phenotype_str}" \
    -N all \
    "$input_dir_path"
echo ${phenotype_str} > "${input_dir_path}/phenotype_names.txt"
