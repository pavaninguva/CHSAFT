#!/bin/sh

while getopts "simtype:" "mode:" option; do
    case "${option}" in
    simtype) SIM=${OPTARG};;
    mode) MODE=${OPTARG};;
    esac
done

echo "Hello $SIM"
echo "Hello $MODE"
