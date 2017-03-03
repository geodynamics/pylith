from collections import OrderedDict
examples = OrderedDict({
    "Step01": {
        "General": {
            "Dimension": 2,
            "Coordinate system": "Cart",
            "Mesh generator": "CUBIT",
            "Cells": "Tri",
            "Problem type": "TD",
            "Time dependence": "S",
        },
        "Boundary Conditions": {
            "Dirichlet": 2,
            "Neumann": 2,
        },
    },
    "Step02": {
        "General": {
            "Dimension": 3,
            "Coordinate system": "Cart",
            "Mesh generator": "CUBIT",
            "Cells": "Tet",
            "Problem type": "TD",
            "Time dependenc": "D",
        },
        "Boundary Conditions": {
            "Absorbing": 5,
        },
    },
})
