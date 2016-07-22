Neighbor-joining is an NPM package for creating phylogenetic trees. It bases on [Rapid Neighbour-Joining](http://pure.au.dk/ws/files/19821675/rapidNJ.pdf) algorithm.

# Installation

To install Neighbor-joining package with NPM use: `npm install neighbor-joining`

# Usage

```javascript
var RNJ = new RapidNeighborJoining(D, taxa, copyDistanceMatrix, taxonIdAccessor);
```
Description of arguments used in initialization:
* **D** - distance matrix (two dimensional array of size NxN)
* **taxa** - array of size N with taxa data (such as strings or objects). Element (taxon) at first index corresponds to first row/column of the distance matrix and so on.
* **copyDistanceMatrix** - flag specifying whether D can be modified or must be copied. Copying might minimally increase the initialization time. Default: `false`.
* **taxonIdAccessor** - function for retrieving the name/identifier of a taxon. It is only called during Newick tree creation. Default: `function(t){ return t.name }`.

# Example
```javascript
var D = [
    [0,  5,  9,  9, 8],
    [5,  0, 10, 10, 9],
    [9, 10,  0,  8, 7],
    [9, 10,  8,  0, 3],
    [8,  9,  7,  3, 0]
];
var taxa = [
    {
        name: "A",
        genotype: "g1"
    },
    {
        name: "B",
        genotype: "g2"
    },
    {
        name: "C",
        genotype: "g3"
    },
    {
        name: "D",
        genotype: "g4"
    },
    {
        name: "E",
        genotype: "g5"
    }
];
var RNJ = new RapidNeighborJoining(D, taxa);
RNJ.run();
var treeObject = RNJ.getAsObject();
var treeNewick = RNJ.getAsNewick();
```

Then, `treeObject` will contain the following object:
```javascript
{
    "taxon": null,
    "length": null,
    "children": [{
        "taxon": {
            "name": "C",
            "genotype": "g3"
        },
        "length": 2,
        "children": []
    }, {
        "taxon": null,
        "length": 2,
        "children": [{
            "taxon": null,
            "length": 3,
            "children": [{
                "taxon": {
                    "name": "A",
                    "genotype": "g1"
                },
                "length": 2,
                "children": []
            }, {
                "taxon": {
                    "name": "B",
                    "genotype": "g2"
                },
                "length": 3,
                "children": []
            }]
        }, {
            "taxon": null,
            "length": 2,
            "children": [{
                "taxon": {
                    "name": "D",
                    "genotype": "g4"
                },
                "length": 2,
                "children": []
            }, {
                "taxon": {
                    "name": "E",
                    "genotype": "g5"
                },
                "length": 1,
                "children": []
            }]
        }]
    }]
}
```
`treeNewick` will keep the following string:
```javascript
"(C:2,((A:2,B:3):3,(D:2,E:1):2):2);"
```