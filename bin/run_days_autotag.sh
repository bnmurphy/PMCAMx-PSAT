#!/bin/bash
# run_days.sh by jinhyok

# -------------------------------------------------------------------------- #
if (( $# < 3 )); then
  echo "Needed to specify a month to run!"
  echo "Usage: $0 PMCAMx_binary Month Server [# of days]"
  echo "     PMCAMx_binary: PMCAMx binary name"
  echo "     Month        : [January|April|July|October]"
  echo "     Server       : [local|server] local or remote server name"
  echo "     [# of days]  : limit the number of days, optional"
  exit 1
fi
# -------------------------------------------------------------------------- #
PMCAMx_binary=$1
RUN_MONTH=$2
SERVER=$3
DAYS=$4

CASE=`echo $PMCAMx_binary | cut -d "_" -f 2`
TAGS=${CASE: -2:2}

CAMx_bin="$HOME/PMCAMx_bin"
if [ $SERVER != "local" ]; then
  CAMx_bin="$HOME/$SERVER/PMCAMx_bin"
fi

export LD_LIBRARY_PATH=$CAMx_bin

CAMx_daily="$CAMx_bin/CAMx_daily.PSAT_JH"
OUTPUT_RUNNING="$CAMx_bin/output/${PMCAMx_binary}/${RUN_MONTH}_running"
#OUTPUT_FINISHED="$CAMx_bin/output/${PMCAMx_binary}/${RUN_MONTH}_$(date +%F_%H%M)"
OUTPUT_FINISHED="$CAMx_bin/output/${PMCAMx_binary}/${RUN_MONTH}_completed"
POSTPROCESSOR="$HOME/PMCAMx-Scripts/bin/post"

# $CAMx_daily $RUN_MONTH $START_DAY

if [ $RUN_MONTH == "January" ]; then
  START_DAY="1"
  END_DAY="31"
elif [ $RUN_MONTH == "April" ]; then
  START_DAY="91"
  END_DAY="120"
elif [ $RUN_MONTH == "July" ]; then
  START_DAY="193"
  END_DAY="209"
elif [ $RUN_MONTH == "October" ]; then
  START_DAY="274"
  END_DAY="304"
else
  echo "Error: Wrong Month!"
  echo "Usage: $0 [January|April|July|October]"
fi

# set END_DAY according to the specified number of days to run
if ((DAYS > 0)); then
  let END_DAY=$((START_DAY + DAYS - 1))
fi

if [ ! -e $OUTPUT_RUNNING ]; then
  echo "mkdir -p $OUTPUT_RUNNING"
  mkdir -p $OUTPUT_RUNNING
fi

if [ $SERVER != "local" ]; then
  touch ${OUTPUT_RUNNING}_`hostname -s`
fi

# pass days already done
while ((START_DAY <= END_DAY)); do
  DAY_PADDED=$(printf "%03d" $START_DAY)
  out_file="${OUTPUT_RUNNING}/4rpos.baseE.${DAY_PADDED}.${PMCAMx_binary}.out"
  echo "Passing $START_DAY ..."
  if ! tail -n 1 $out_file | grep "TOTAL COARSE GRID TIME STEPS" > /dev/null 2>&1; then
    echo "Starting from $START_DAY ..."
    break
  fi
  let START_DAY++
done

TRIAL="0"
while ((START_DAY <= END_DAY)); do
  echo "$CAMx_daily $PMCAMx_binary $RUN_MONTH $START_DAY $SERVER $TAGS"

  #$CAMx_daily $PMCAMx_binary $RUN_MONTH $START_DAY $SERVER $TAGS
  #DAY_PADDED=$(printf "%03d" $START_DAY)
  #out_file="${OUTPUT_RUNNING}/4rpos.baseE.${DAY_PADDED}.${PMCAMx_binary}.out"

  DAY_PADDED=$(printf "%03d" $START_DAY)
  out_file="${OUTPUT_RUNNING}/4rpos.baseE.${DAY_PADDED}.${PMCAMx_binary}.out"

  while ! $CAMx_daily $PMCAMx_binary $RUN_MONTH $START_DAY $SERVER $TAGS; do
    if ((TRIAL > 3)); then
      echo "Error! Look at $out_file"
      exit 1
    fi
    echo "Error, but trying again ... $TRIAL"
    sleep 10
    let TRIAL++
  done

  # Stop if it didn't finish w/o error
  if ! tail -n 1 $out_file | grep "TOTAL COARSE GRID TIME STEPS" > /dev/null 2>&1; then
    echo "Error! Look at $out_file:"
    tail -n 10 $out_file
    exit 1
  fi

  let START_DAY++
done

# Rename output folder when it's all done
mv $OUTPUT_RUNNING $OUTPUT_FINISHED
# Post-process output
$POSTPROCESSOR PSAT $OUTPUT_FINISHED

# Done
echo "ALL DONE!"
