#!/bin/bash

set -ex

/opt/conda/bin/conda env create --file mkdocs-env.yml
source /opt/conda/bin/activate mkdocs
mkdocs build
