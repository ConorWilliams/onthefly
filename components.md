# Components

## Data components

**AtomArray** The default adaptive container type used in the libatom; an AtomArray models a vector of "atom" types but decomposes the atom type and stores each member in a separate vector. This enables efficient cache use. The members of the "atom" are described through a series of template parameters which should inherit from otf::AtomArrayMem. A selection of canonical members are provided in libatom/atomvector.hpp. Example of use:

```c++

#include "libatom/atomvector.hpp"

using namespace otf;

AtomArray<Pos, AtomicNum> atoms;

atoms.push_back({0,0,0}, 1) // Add a hydrogen atom to the origin

Vec3 xyz = atoms(Pos{}, 0) // Get the position of the zeroth atom

std::size_t n = = atoms(AtomicNum{}, 0) // Get the atomic number of the zeroth atom

atoms(Pos{}) += 1. // Add 1 to each of every atoms coordinates

```

**System** : list of **Atoms**, **Simbox** 
**Neighbour List** : list of list of pointers to **Atoms** + data (has a broadcast feature)
**Geometry** : **Key**, **Fingerprint**, list of **Atoms**
**LE** : reference **Geo**, list of **Mechanisms** 
**Catalogue** : collection of **LE**

## Complex transformations



NeighList atom_vec{x_cannon, image_id, neigh idx list}, markers {a, f, g1, g2, g2}, simbox, head, tail



NeighReduce :
 
NeighLister :



NeighList :  **x**, **Box** -> periodically resolved list of neighbour idxs + xyz for each atom

Bases : **x**, **z**, **Box** -> energy, force, hessian 

Minimise :  **x**, **z**, **Box** -> minimum xyz

**System** -> closest MEC
**System** -> closest SP

**System** -> system hash

**System** -> list of **Geometry**


grep for  @code{.cpp} and fixup :)


