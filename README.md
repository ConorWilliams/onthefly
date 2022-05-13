<br />
<p align="center">
  <img src="./images/dimer.png" height="200" />
</p>
<br />

# OLKMC

Welcome to the ...

Local environment, tolerant, off-lattice, kinetic, Monte Carlo, Simulation/Simulator, Framework, Atomic

L, Le, T, Ol, O, K, Mc, M, S, F, A

## Building

```bash
git clone https://github.com/ConorWilliams/olkmc
cd olkmc
mkdir build 
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
./olkmc ../data/test.toml 

```

- Superbasin dynamic barrier height
- Add resolved mode; throws error for now.
- Hessian
- Explore 3% ?bug

## Minor improvements (QoL/Tidy)
- Pre-prep for local ff, HTST calc, generalised ff
- Decouple ff, sp search, minimiser, topoclassify, reconstruct, superbasin
- Detect too large mechanisms 
- Proper atom-type reading in force/xyz files
- .xyz read/write (dump system object)
- Binary format: [cereal](https://github.com/USCiLab/cereal)
- Meta-data in catalog (force field, etc)
- Argparsing + config file [toml++](https://github.com/marzer/tomlplusplus/) + xyz load: [structopt](https://github.com/p-ranav/structopt)
- Efficient ghost production (edge cell only iteration) 
- Generalise supercell to tri-clinic (abstract interface to supercell)
- Abstract/encapsulate parallelisation
- State hashing for superbasin pre-conditioning


## Majour improvements/extensions

Feature | How much better for experiment | Interest | Score | Rank

    - Error/confidence rate catalogue                   5	5	25	1
    - HTST rate constant + cache SP                     5	5	25	1
    - Think about better topo classification            4	5	20	2
    - Explore alternative force-fields                  5	4	20	3
    - Local forces                                      3	5	15	4
    - Think about simultanious dimer method (swarm)     3	4	12	5
    - Quantum veification of frequent topos             5	2	10	
    - MPI parallelisation master/slave                  3	3	9	
    - LCR sorting (paper + sorting networks)            2	4	8	
    - Two level parallelisation (parallel force field)  2	3	6	
    - Initial displacement (hypersphere)                2	3	6	
    - Preconditioner                                    2	3	6	
    - GPU force field                                   2	2	4	


Since literature:

    - Error estimator on catalogue completeness
    - Degeneracy coefficient in KMC N-fold way, link LE -> atom list
    - Dynamically localised force-fields [[bland2015 p.131]]
    - Dimer method alternatives (k-art, Lancoz)
    - Thermal expansion of lattice (I think in the SP recycling paper)
