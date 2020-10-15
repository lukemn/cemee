#!/bin/bash

# Run from the parent directory cemee/multiwormtracker (containing conversion scripts in `wcon` and male/herm classification in `male`),
# with a single argument - a file listing the names of experiment directories to process within MWTdata/Raw_behavior/.
# MWT data to be processed should be in MWTdata/Raw_behavior.
# A tar/bzip2 of raw data w made in MWTdata/Raw_behavior,
# Choreography is run, saving raw and processed tars in MWTdata/Processed_behavior
# then the output is parsed and saved as RData in MWTdata/Parsed_behaviour.
# set Google drive rclone remote to also backup here.

#### OPTIONAL
# convert blobs to wcon and also save in Parsed_behaviour
MAKE_WCON=false
NP=3
CONV=`pwd`/wcon/blobs2wcon.R
# requires custom config file, path to which should be listed 2nd (tab delimited) column of the input experiment list. See conversion script CONV for format details,
# and example in wcon/mwt_config
# predict male numbers and frequencies from xgb trained model
PREDICT_MALES=true
# save Processed_behaviour Parsed_behaviour on Gdrive
MALE=`pwd`/male/splitTracksBySex.R
XGB=`pwd`/male/xgb_preds.rda
# backup to google drive via rclone, saving Processed_behaviour, Parsed_behaviour under each experiment dir
BACKUP_GOOGLE=false
GDRIVE=remote:/MWT
####


# RAWD must exist
DD=`pwd`/MWTdata
RAWD=$DD/Raw_behavior/
# pass file listing parent experiments within RAWD, containing MWT directories (YYYYMMDD_HHMMSS format) with blobs to process
if [ "$#" -lt 1 ]; then echo -e "Pass a file listing directories in\n$RAWD"; exit 1; fi
EXP=$(grep -v "#" $1)
if [ $(echo $EXP | wc -w) = 0 ]; then echo -e "No experiments to process.\nPass a file listing directories in $RAWD"; exit 1; fi
# these will be created if they don't exist
PROD=$DD/Processed_behavior
PARD=$DD/Parsed_behavior
# Thiago's scripts
# tar
ARCHIVE=`pwd`/wcon/celegans_behavior-master/Utils/archive_exper_behavior_data.sh 
# call run_choreography.sh (the full call, including traits to extract)
# then parse_raw_exper_data.R to parse, rename traits, save .RData in Parsed_behaviour
PROCESS=`pwd`/wcon/celegans_behavior-master/Utils/process_archived_exper_behavior_data.sh

if [ $(echo $EXP | wc -w) = 0 ]; then echo "No experiments to process in $1"; exit 1; fi
if [ ! -s $RAWD ]; then echo "$RAWD does not exist"; exit 1; fi 
if [ ! -s $ARCHIVE ]; then echo -e "Can't find\n$ARCHIVE"; exit 1; fi
if [ ! -s $PROCESS ]; then echo -e "Can't find\n$PROCESS"; exit 1; fi
if [[ $PREDICT_MALES && ! -s $MALE ]]; then echo -e "Can't find\n$MALE"; exit 1; fi
if [[ $PREDICT_MALES && ! -s $XGB ]]; then echo -e "Can't find\n$XGB"; exit 1; fi
if [[ $MAKE_WCON && ! -s $CONV ]]; then echo -e "Can't find\n$CONV"; exit 1; fi
if $BACKUP_GOOGLE; then
    which rclone
    OK=$?
    if [ OK = 0 ]; then
      REM=$(echo $GDRIVE | cut -f1 -d":")
      REMOTE=$(rclone listremotes | grep $REM | wc -l)
      if [ $REMOTE = 0 ]; then echo "Supplied remote $REM not found by rclone"; exit 1; fi
    else
      echo "Backup to Google drive requested but rclone is not in PATH"
      exit 1
    fi
fi

# echo $MAKE_WCON
# echo $PREDICT_MALES
# echo $BACKUP_GOOGLE

for L in $EXP; do

    E=$(echo $L | cut -f1)
    RAWE=$RAWD/$E
    PROE=$PROD/$E
    PARE=$PARD/$E

    mkdir -p $PROE $PARE
    cd $RAWD/$E && $ARCHIVE $PROE &> $E.archive.log
    cd $PROE && $PROCESS $PARE &> $E.process.log
    
    # convert to wcon (running on blobs in Raw_behaviour)
    if $MAKE_WCON; then
        CFG=$(echo $L | cut -f2)
        $CONV $CFG $NP
        cp $RAWE/*zip $PROE
    fi
    
    # split tracks by sex, predict freqs, numbers
    if $PREDICT_MALES; then
        $MALE $PARE $XGB $E
    fi
    
    # backup
    if $BACKUP_GOOGLE; then
        rclone copy $PARE $GDRIVE/$E
        rclone copy $PROE $GDRIVE/$E
    fi
done

