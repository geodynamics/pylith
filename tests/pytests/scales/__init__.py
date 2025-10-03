from . import (
    TestGeneral,
    TestQuasistaticElasticity,
    #TestDynamicElasticity,
)

def test_modules():
    return [
        TestGeneral,
        TestQuasistaticElasticity,
        #TestDynamicElasticity,
    ]
