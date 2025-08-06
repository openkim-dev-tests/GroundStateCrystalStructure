Reference Elemental Ground State Test
=====================================

This test returns reference ground state structures and energies for each element.
The results from this test are useful when a reference structure is required in some downstream test, such as vacancy tests (used as a reservoir).
This test driver works by querying results from the `EquilibriumCrystalStructure <https://openkim.org/id/EquilibriumCrystalStructure__TD_457028483760_002>`_ test driver using element specific reference structures following `CHIPS-FF <https://github.com/usnistgov/chipsff/blob/main/chipsff/chemical_potentials.json>`_.
Although the reference prototypes are independent of model, the resulting structure and energy of the prototypes are model-dependent.
