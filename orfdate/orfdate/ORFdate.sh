#!/bin/bash

names=${names:-false}
tree=${tree:-false}

while [ $# -gt 0 ]; do
   if [[ $1 == *"--"* ]]; then
        param="${1/--/}"
        declare $param="$2"
        # echo $1 $2 // Optional to see the parameter:value result
   fi
  shift
done

if ! $names; then
	echo "As no list of names was provided, looking for the provided tree to use as source of names."
