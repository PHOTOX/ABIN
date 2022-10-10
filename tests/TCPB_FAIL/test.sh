#/bin/bash

set -euo pipefail
export ABINEXE=$1

if [[ $1 = "clean" ]];then
  rm -f abin.out ERROR
  exit 0
fi

TCPB_PORT=$((1025+RANDOM))
TCPB_HOST=localhost
TCPB_IN="tc.inp"

$ABINEXE -i input.in -x mini.xyz \
  --tcpb-host $TCPB_HOST \
  --tcpb-port $TCPB_PORT \
  --tcpb-input-file $TCPB_IN &> abin.out
