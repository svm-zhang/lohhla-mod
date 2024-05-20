#!/bin/bash

set -e

mkdir -p "$PREFIX/bin"
mkdir -p "$PREFIX/lib/R/library/hlalohReforged"

# copy scripts to bin folder
cp "$RECIPE_DIR/R/hlalohReforged.R" "$PREFIX/bin/hlalohReforged"
cp "$RECIPE_DIR/R/plot.R" "$PREFIX/bin/hlalohplot"

# make it executable
chmod +x "$PREFIX"/bin/*

# copy perl module to lib folder
cp "$RECIPE_DIR/R/bamer.R" "$PREFIX/lib/R/library/hlalohReforged"
cp "$RECIPE_DIR/R/cli.R" "$PREFIX/lib/R/library/hlalohReforged"
cp "$RECIPE_DIR/R/loh.R" "$PREFIX/lib/R/library/hlalohReforged"
cp "$RECIPE_DIR/R/pairwise_aln.R" "$PREFIX/lib/R/library/hlalohReforged"
cp "$RECIPE_DIR/R/pathio.R" "$PREFIX/lib/R/library/hlalohReforged"
