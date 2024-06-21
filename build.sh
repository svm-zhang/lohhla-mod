#!/bin/bash

set -e

mkdir -p "$PREFIX/bin"
mkdir -p "$PREFIX/lib/R/library/lohhlamod"

# copy scripts to bin folder
cp "$RECIPE_DIR/R/lohhlamod.R" "$PREFIX/bin/lohhlamod"
cp "$RECIPE_DIR/R/lohhlaplot.R" "$PREFIX/bin/lohhlaplot"

# make it executable
chmod +x "$PREFIX"/bin/*

# copy perl module to lib folder
cp "$RECIPE_DIR/R/bamer.R" "$PREFIX/lib/R/library/lohhlamod"
cp "$RECIPE_DIR/R/cli.R" "$PREFIX/lib/R/library/lohhlamod"
cp "$RECIPE_DIR/R/loh.R" "$PREFIX/lib/R/library/lohhlamod"
cp "$RECIPE_DIR/R/pairwise_aln.R" "$PREFIX/lib/R/library/lohhlamod"
cp "$RECIPE_DIR/R/pathio.R" "$PREFIX/lib/R/library/lohhlamod"
