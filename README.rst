Reference Elemental Ground State Test
=====================================

This test finds the lowest energy structure for each single-element species supported by the relevant model.
In addition to simplifying the queries needed to find the lowest energy structures predicted by a model, the results from this test may also be useful when a reference structure is required in some downstream task, such as vacancy tests (used as a reservoir).
This test driver works by taking all supported single-element prototypes for a given species, iterating through each of them, relaxing the structures according to the methodology of the `EquilibriumCrystalStructure <https://openkim.org/id/EquilibriumCrystalStructure__TD_457028483760_002>` test driver, and recording the lowest energy structure.
This test driver was created so that other tests which may depend on reference structures can use tests derived from this driver as a dependency, rather than resorting to depending on several tests (one for each single) element prototype-ehich may become problematic due to failed runs, etc. 
