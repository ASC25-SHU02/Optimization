#!/bin/bash
enable_timing="yes"
run_with_timing() {
	local cmd="$@"
    echo "Starting: $cmd"
    if [ "$enable_timing" == "yes" ]; then
        time eval "$@"
    else
        eval "$@"
    fi
	echo "##################...Ending...##################"
	echo
}
