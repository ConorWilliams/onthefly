# OLKMC

## Intro

This branch contains the code used in the paper ...

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
