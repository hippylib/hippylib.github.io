#!/bin/sh

set -ex

conda-env create --file mkdocs-env.yml
source activate mkdocs
mkdocs build
