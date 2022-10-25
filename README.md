# OLKMC

## Intro

This branch contains the code used in the paper: 

[Accelerating off-lattice kinetic Monte Carlo simulations to predict hydrogen vacancy-cluster interactions in alpha–Fe](https://doi.org/10.1016/j.actamat.2022.118452)

For the latest code see https://github.com/ConorWilliams/openFLY

## Building

```bash
git clone https://github.com/ConorWilliams/olkmc
cd olkmc
mkdir build 
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
./olkmc ../data/test.toml 
