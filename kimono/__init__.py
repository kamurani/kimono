# Kinase Motif Notation (KIMONO) is a Python package for kinase motif analysis.
# Author: Cam İmran  <github.com/kamurani>
# License: Apache License 2.0
# Code Repository: https://github.com/kamurani/kimono

from pathlib import Path


PROJECT_ROOT_DIR = Path(__file__).parent.parent

DATA_DIR = PROJECT_ROOT_DIR / "data"

SEQUENCE_SAVE_DIR = DATA_DIR / "sequences"
STRUCTURE_SAVE_DIR = DATA_DIR / "structures"
RESULTS_SAVE_DIR = DATA_DIR / "results"

if __name__ == "__main__":
    print(PROJECT_ROOT_DIR)
    print(DATA_DIR)
    print(SEQUENCE_SAVE_DIR)
    print(STRUCTURE_SAVE_DIR)