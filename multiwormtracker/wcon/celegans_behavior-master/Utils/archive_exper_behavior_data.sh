#!/bin/bash
#
# Script to archive experimental data that has been acquired. 
# The input is a set of folders, each folder corresponding to a single "experiment" (acquisition session).
# Each such folder is named based on the format 
# Archive the behavior data: for each experiment folder "[date]_[time]".
# This script browses through all the experiment folders present in the path from which 
# the script is invoked, and for each such folder compresses the data as a single bz2 file in a target folder.
#
# Typical usage: UTILS_PATH/./archive_exper_behavior_data.R "/behavior_data/New_data/Block_10"
# which will process all the experiment folders in the working folder, and store the 
# generated files in the destination folder "/behavior_data/New_data/Block_10"
#
# Created on September 2016
# Added to repository on June 2017
#

set -e

if [ $# -ne 1 ]; then
    echo Usage: ./archive_behavior_data.sh DEST_FOLDER
    exit -1
fi
dest_folder="$1"

#
# Exit on error, after displaying a message
#
function exitOnError {
    echo $1
    exit -1
}

function showWarning {
    echo $1 >&2
}

# List all the folders
readarray -t entry_list < <(exec ls)

echo "Entries: ${entry_list[@]}"

# arr=($(grep -o '[[:digit:]]\{8\}_[[:digit:]]\{6\}' ${exper_list[@]}))
# echo "Filtered: ${arr[@]}"

declare -a exper_list
for entry in "${entry_list[@]}"; do
    #echo "${entry}"
    if [[ "${entry}" =~ ^[0-9]{8}_[0-9]{6}$ ]]; then
        #echo "Match: ${entry}"
        exper_list+=("${entry}")
    else
        echo "No match: ${entry}"
    fi
done

echo "Experiments: ${exper_list[@]}"

echo "${#entry_list[@]} -> ${#exper_list[@]}"

if [[ "${#entry_list[@]}" -ne "${#exper_list[@]}" ]]; then
    showWarning "Some entries do not match the pattern for experiment folders, and will be discarded"
fi

# ls | grep -o '[[:digit:]]\{8\}_[[:digit:]]\{6\}'

# Initialize the output folder
mkdir -p "${dest_folder}"

# Sanitize the folders that (could) correspond to experiment; show warnings for the filtered cases
for exper_id in "${exper_list[@]}"; do

    echo "* Processing experiment ${exper_id}"
    dest_file="${dest_folder}/${exper_id}.tar.bz2"
    echo "  Output file = ${dest_file}"

    # Check if this experiment has already been processed
    if [[ -f "${dest_file}" || -d "${dest_file}" ]]; then
        showWarning "  File corresponding to experiment id ${i} (${dest_file}) already exists. Skipping..."
    else
        echo "  Compressing..."
        tar cjf "${dest_file}" "${exper_id}"

        echo "  Checking..."
        tar -tf "${dest_file}" > /dev/null
    fi
done
