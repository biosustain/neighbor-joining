Neighbor-joining is an NPM package for creating phylogenetic trees. It bases on [Rapid Neighbour-Joining](http://pure.au.dk/ws/files/19821675/rapidNJ.pdf) algorithm.

# Installation

To install Ancestry package with NPM use: `npm install ancestry`

# Usage

```
var RNJ = new RapidNeighborJoining(D, taxa, copyDistanceMatrix, taxonIdAccessor);
```
Description of arguments used in initialization:
* D - distance matrix (two dimensional array of size NxN)
* taxa - array of size N with taxa data (such as strings or objects). Element (taxon) at first index corresponds to first row/column of the distance matrix and so on.
* copyDistanceMatrix - flag specifying whether D can be modified or must be copied. Copying might minimally increase the initialization time. Default: false.
* taxonIdAccessor - function for retrieving the name/identifier of a taxon. It is only called during Newick tree creation. Default: function(t){ return t.name }.

# Example
```

```