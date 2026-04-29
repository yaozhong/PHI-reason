#!/bin/bash
# PHIreason Pipeline — User Configuration
#
# Set PHI_PROJECT_ROOT to the directory containing ws/, data/, experiments/.
# Clone PHIreason_pipeline/ into that directory (recommended), or override
# this variable explicitly before running any pipeline script.
#
# Quickstart:
#   export PHI_PROJECT_ROOT=/path/to/your/project
#   bash PHIreason_pipeline/03_experiments/run_experiments.sh
#
# Tool overrides (optional — default to PATH lookup):
#   DIAMOND_BIN=/path/to/diamond
#   PRODIGAL_BIN=/path/to/prodigal

PHI_PROJECT_ROOT="${PHI_PROJECT_ROOT:-}"

if [[ -z "$PHI_PROJECT_ROOT" ]]; then
    _PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PHI_PROJECT_ROOT="$(dirname "$_PIPELINE_DIR")"
fi

export PHI_PROJECT_ROOT
