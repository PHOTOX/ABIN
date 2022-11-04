#/bin/bash

set -uo pipefail
export ABINEXE=$1

if [[ $1 = "clean" ]];then
  rm -f abin.out ERROR
  exit 0
fi

TCPB_PORT=$((1025+RANDOM))
TCPB_HOST=localhost
TCPB_IN="tc.inp"

# This one will fail simply because TC server is not available
$ABINEXE -i input.in -x mini.xyz \
  --tcpb-host $TCPB_HOST \
  --tcpb-port $TCPB_PORT \
  --tcpb-input-file $TCPB_IN &> abin.out
mv ERROR ERROR1

# TC input file does not exist
$ABINEXE -i input.in -x mini.xyz \
  --tcpb-port $TCPB_PORT \
  --tcpb-input-file invalid &>> abin.out
mv ERROR ERROR2

# Negative port
$ABINEXE -i input.in -x mini.xyz \
  --tcpb-host $TCPB_HOST \
  --tcpb-port -1 \
  --tcpb-input-file $TCPB_IN &>> abin.out
mv ERROR ERROR3

# Port conflicting with system-reserved ports
$ABINEXE -i input.in -x mini.xyz \
  --tcpb-host $TCPB_HOST \
  --tcpb-port 80 \
  --tcpb-input-file $TCPB_IN &>> abin.out
mv ERROR ERROR4

# Port out of range
$ABINEXE -i input.in -x mini.xyz \
  --tcpb-host $TCPB_HOST \
  --tcpb-port 66000 \
  --tcpb-input-file $TCPB_IN &>> abin.out
mv ERROR ERROR5
