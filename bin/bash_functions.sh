#!/usr/bin/env bash

send_mail=true

# Get path of this directory
DIR="$(dirname "${BASH_SOURCE[0]}")"
DIR="$(realpath "${DIR}")"

msg() {
  # A function, designed for logging, that prints current date and timestamp
  # as a prefix to the rest of the message. Default sends to stdout but can
  # redirected during execution.
  echo "[$(date '+%Y-%b-%d %a %H:%M:%S')] $@"
}

check_if_file_exists_allow_seconds() {
  # Boolean test to confirm a file exists within a specified maximum time to
  #  wait.

  # Returns:
  #  0: (true) file exists
  #  1: (false) file missing

  # Parameters:
  #  $1 = file
  #  $2 = maximum seconds to wait for file to appear
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

expandsi() {
  #https://stackoverflow.com/a/52799497
  # Yzmir Ramirez: CC BY-SA 4.0
  ## expanded by WAO for case insensitivity
  ## Does not work if b|B is included (e.g. K works, but KB|kb|Kb... will return 0)
  ## Could just remove all B's from input string prior to conversion, if necessary?
  echo ${1} | awk 'function pp(p){printf "%u\n",$0*1024^p} /[0-9]$/{print $0}/[k|K]$/{pp(1)}/[m|M]$/{pp(2)}/[g|G]$/{pp(3)}/[t|T]$/{pp(4)}/[^0-9kmgtKMGT]$/{print 0}'
}

verify_minimum_file_size() {
  # Boolean test to ensure the filepath is a file, is non-zero size, and
  #  is at least the minimum specified size (in Bytes).

  # Returns:
  #  0: (true) file is at least the minimum specified size
  #  1: (false) file is smaller than the specified size

  # Parameters:
  #  $1=filename
  #  $2=file description
  #  $3=minimum size in Bytes
  #   (optionally can specify k, M, or G suffix after a number for big numbers)

  #updated by WAO to remove 'find -size' that is not posix-standard
  min_size=$(expandsi ${3})
  file_size=$(wc -c <"${1}")
  if [ -f "${1}" ]; then
    if [ -s "${1}" ]; then
      if [[ "${file_size}" -gt "${min_size}" ]]; then
        return 0
      else
        msg "ERROR: ${2} file ${1} present but too small (less than ${3})" >&2
        return 1
      fi
    else
      msg "ERROR: ${2} file ${1} present but empty" >&2
      return 1
    fi
  else
    msg "ERROR: ${2} file ${1} absent" >&2
    return 1
  fi
}
