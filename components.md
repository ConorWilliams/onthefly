# Components

## Data components

**System** : list of **Atoms**, **Simbox** 
**Neighbour List** : list of list of pointers to **Atoms** + data (has a broadcast feature)
**Geometry** : **Key**, **Fingerprint**, list of **Atoms**
**LE** : reference **Geo**, list of **Mechanisms** 
**Catalogue** : collection of **LE**

## Complex transformations



NeighReduce :
 
NeighLister :



NeighList :  **x**, **Box** -> periodically resolved list of neighbour idxs + xyz for each atom

Potentials : **x**, **z**, **Box** -> energy, force, hessian 

Minimise :  **x**, **z**, **Box** -> minimum xyz

**System** -> closest MEC
**System** -> closest SP

**System** -> system hash

**System** -> list of **Geometry**


