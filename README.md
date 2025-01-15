# Cancer Progression Simulator

[![DOI](https://zenodo.org/badge/836887984.svg)](https://doi.org/10.5281/zenodo.14651855)


This project simulates the progression of cancer in patients over time, including diagnosis, remission, and death. 
The simulation can be run with or without the use of blood tests for cancer detection.

## Features

- Simulates cancer progression through stages I-IV
- Models monthly survival probabilities and remission chances
- Allows for periodic blood tests to improve early detection
- Generates patient trajectories and datasets for analysis

## Installation

1. Clone the repository:
    ```sh
    git clone https://github.com/belapyc/cancer-simulator.git
    cd cancer-simulator
    ```

2. Create a conda environment using the `env.txt` file:
    ```sh
    conda create --name cancer-simulator --file env.txt
    conda activate cancer-simulator
    ```

## Usage

For examples of usage refer to `simulation.ipynb`.
