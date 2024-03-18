#!/bin/sh

if [ -z "${VALOR_FMWK}" ] ; then
  if [ -z "${2}" ] ; then
    echo "VALOR framework must already be set up or specified in argument \$2 ($2)"
    return
  fi
  export VALOR_FMWK="$(readlink -m $2)" # Requires GNU readlink - should work on SL builds (but not Macs - try greadlink)
  
  . "${VALOR_FMWK}/setup/setup_at_cvmfs_root6.sh" "${VALOR_FMWK}"
fi

export VALOR_EXPT="$(readlink -m "$1")"

if [ -z "${VALOR_EXPT}" ]; then
  sdir="$( readlink -m "${BASH_SOURCE[0]}" )"
  export VALOR_EXPT="$( dirname "$( dirname "${sdir}" )" )" # Directory above the one containing this file
fi

export LD_LIBRARY_PATH="${VALOR_EXPT}/lib:${LD_LIBRARY_PATH}"
export PATH="${VALOR_EXPT}/bin:${VALOR_EXPT}/prod:${PATH}"

export VALOR_DATA="${VALOR_EXPT}/data"
echo "Set \$VALOR_DATA to ${VALOR_DATA} - you may want to change this."
