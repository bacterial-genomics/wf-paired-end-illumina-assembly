#!/usr/bin/env bash

send_mail=true

# Get path of this directory
DIR="$(dirname "${BASH_SOURCE[0]}")"
DIR="$(realpath "${DIR}")"

msg() {
  # A function, designed for logging, that prints current date and timestamp
  #  as a prefix to the rest of the message. Default sends to stdout but can
  #  redirected during execution.
  echo "[$(date '+%Y-%b-%d %a %H:%M:%S')] $@"
}

check_if_file_exists_allow_seconds() {
  # Parameters:
  #              $1 = file
  #              $2 = maximum seconds to wait for file to appear
  elapsed=0
  while [ ! -f "${1}" ]; do
    sleep 1
    ((elapsed++))
    if [ "${elapsed}" -eq "${2}" ]; then
      msg "ERROR: ${1} cannot be found after waiting ${2} seconds" >&2
      return 1
    fi
  done
  return 0
}

verify_file_minimum_size()
{
  # $1=filename
  # $2=file description
  # $3=size in Bytes (requires c, k, M, or G prefix)
  if [ -f  "${1}" ]; then
    if [ -s  "${1}" ]; then
      if [[ $(find -L "${1}" -type f -size +"${3}") ]]; then
        return
      else
        size=$(echo ${3} | sed 's/c//g')
        msg "ERROR: ${2} file ${1} present but too small (less than ${size} bytes)" >&2
        false
      fi
    else
      msg "ERROR: ${2} file ${1} present but empty" >&2
      false
    fi
  else
    msg "ERROR: ${2} file ${1} absent" >&2
    false
  fi
}

