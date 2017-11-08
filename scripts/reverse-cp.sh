#!/bin/sh
# A cp that you simply put the destination FIRST rather than last in the args list, so you can use it with xargs.
dest="$1"
shift
cp -p "$@" "$dest"
