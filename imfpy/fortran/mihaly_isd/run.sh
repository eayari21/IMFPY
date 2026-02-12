#!/usr/bin/env bash
set -e

FC=/opt/homebrew/bin/gfortran
SRC=isd-2022.f
EXE=isd

echo "Compiling $SRC..."
$FC -O3 -Wall $SRC -o $EXE

echo "Running $EXE..."
./$EXE
