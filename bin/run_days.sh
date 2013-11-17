#!/bin/bash
# run_days.sh by jinhyok

# -------------------------------------------------------------------------- #
if (( $# < 5 )); then
  echo "Needed to specify a month to run!"
  echo "Usage: $0 PMCAMx_binary Month Is_Remote [# of days] [# of tags]"
  echo "     PMCAMx_binary : PMCAMx binary name"
  echo "     Month        : [January|April|July|October]"
  echo "     Is_Remote    : [remote|local] Remote run through sshfs"
  echo "     [# of days]  : limit the number of days"
  echo "     [# of tags]  : the number of tags"
  exit 1
fi
# -------------------------------------------------------------------------- #

PMCAMx_binary=$1
RUN_MONTH=$2
IS_REMOTE=$3
DAYS=$4
TAGS=$5

CAMx_bin="$HOME/PMCAMx_bin"
if [ $IS_REMOTE == "remote" ]; then
  CAMx_bin="$HOME/surly/PMCAMx_bin"
fi

export LD_LIBRARY_PATH=$CAMx_bin

CAMx_daily="$CAMx_bin/CAMx_daily.PSAT_JH"
OUTPUT_RUNNING="$CAMx_bin/output/${PMCAMx_binary}/${RUN_MONTH}_running"
OUTPUT_FINISHED="$CAMx_bin/output/${PMCAMx_binary}/${RUN_MONTH}_$(date +%F_%H%M)"

# POST="$HOME/PMCAMx_Utils/POST_PROCESSING/after_CB4"

# $CAMx_daily $RUN_MONTH $START_DAY_OF_YEAR

if [ $RUN_MONTH == "January" ]; then
  START_DAY_OF_YEAR="1"
  END_DAY_OF_YEAR="31"
elif [ $RUN_MONTH == "April" ]; then
  START_DAY_OF_YEAR="91"
  END_DAY_OF_YEAR="120"
elif [ $RUN_MONTH == "July" ]; then
  START_DAY_OF_YEAR="193"
  END_DAY_OF_YEAR="209"
elif [ $RUN_MONTH == "October" ]; then
  START_DAY_OF_YEAR="274"
  END_DAY_OF_YEAR="304"
else
  echo "Error: Wrong Month!"
  echo "Usage: $0 [January|April|July|October]"
fi

# set END_DAY_OF_YEAR according to the specified number of days to run
if ((DAYS > 0)); then
  let END_DAY_OF_YEAR=$((START_DAY_OF_YEAR + DAYS - 1))
fi

if [ ! -e $OUTPUT_RUNNING ]; then
  echo "mkdir -p $OUTPUT_RUNNING"
  mkdir -p $OUTPUT_RUNNING
fi

if [ $IS_REMOTE == "remote" ]; then
  touch ${OUTPUT_RUNNING}_`hostname -s`
fi

# if [ -e $OUTPUT_RUNNING/.running ]; then
#  echo "Simulation is runnig ..."
#  echo "Delete $OUTPUT_RUNNING/.running to restart."
#  exit 1
# fi
# touch $OUTPUT_RUNNING/.running

# pass days already done
while ((START_DAY_OF_YEAR <= END_DAY_OF_YEAR)); do
  DAY_PADDED=$(printf "%03d" $START_DAY_OF_YEAR)
  out_file="$OUTPUT_RUNNING/4rpos.baseE.$DAY_PADDED.$PMCAMx_binary.out"
  echo "Passing $START_DAY_OF_YEAR ..."
  if ! tail -n 1 $out_file | grep "TOTAL COARSE GRID TIME STEPS" > /dev/null 2>&1; then
    echo "Starting from $START_DAY_OF_YEAR ..."
    break
  fi
  let START_DAY_OF_YEAR++
done

TRIAL="0"
while ((START_DAY_OF_YEAR <= END_DAY_OF_YEAR)); do
  echo "$CAMx_daily $PMCAMx_binary $RUN_MONTH $START_DAY_OF_YEAR $IS_REMOTE $TAGS"
  while ! $CAMx_daily $PMCAMx_binary $RUN_MONTH $START_DAY_OF_YEAR $IS_REMOTE $TAGS; do
    DAY_PADDED=$(printf "%03d" $START_DAY_OF_YEAR)
    out_file="${OUTPUT_RUNNING}/4rpos.baseE.${DAY_PADDED}.${PMCAMx_binary}.out"
    #if ((TRIAL > 5)); then
    if ((TRIAL > 1)); then
      echo "Error! Look at $out_file"
      exit 1
    fi
    sleep 10
    let TRIAL++
  done
  let START_DAY_OF_YEAR++
done

# mv $OUTPUT_RUNNING $OUTPUT_FINISHED

# mkdir -p ${OUTPUT_FINISHED}_processed
# $POST ${OUTPUT_FINISHED}
# gzip ${OUTPUT_FINISHED}_processed/*

# gzip ${OUTPUT_FINISHED}/*.mass
# gzip ${OUTPUT_FINISHED}/*.out

# done
echo "ALL DONE!"
