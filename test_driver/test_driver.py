from kim_tools import SingleCrystalTestDriver
from ase import Atoms
from typing import Any, Optional, List, Union, Dict, IO
from reference_map import reference_map
import kimvv
class TestDriver(SingleCrystalTestDriver):
    # keep same signature so __call__ works
    def _setup(
        self,
        material: str,
        cell_cauchy_stress_eV_angstrom3: Optional[List[float]] = None,
        temperature_K: float = 0,
        **kwargs,
    ) -> None:
        """
        Args:
            material:
                Material in this case is simply a str containing a single element type,
                e.g., 'Au'.
            cell_cauchy_stress_eV_angstrom3:
                Cauchy stress on the cell in eV/angstrom^3 (ASE units) in
                [xx, yy, zz, yz, xz, xy] format. This is a nominal variable, and this
                class simply provides recordkeeping of it. It is up to derived classes
                to implement actually imposing this stress on the system.
            temperature_K:
                The temperature in Kelvin. This is a nominal variable, and this class
                simply provides recordkeeping of it. It is up to derived classes to
                implement actually setting the temperature of the system.
        """
        if cell_cauchy_stress_eV_angstrom3 is None:
            cell_cauchy_stress_eV_angstrom3 = [0, 0, 0, 0, 0, 0]
        self.reference_element =  material

    def _calculate(self, reference_info, temperature_K=None, cell_cauchy_stress_eV_angstrom3=None, **kwargs):
        '''
        Args:
            reference_info:
                Result from a KIM query or EquilibriumCrystalStructure test
                Result will be in form {'<element symbol>': [binding-potential-energy-crystal instance]}
            cell_cauchy_stress_eV_angstrom3:
                Cauchy stress on the cell in eV/angstrom^3 (ASE units) in
                [xx, yy, zz, yz, xz, xy] format. This is a nominal variable, and this
                class simply provides recordkeeping of it. It is up to derived classes
                to implement actually imposing this stress on the system.
            temperature_K:
                The temperature in Kelvin. This is a nominal variable, and this class
                simply provides recordkeeping of it. It is up to derived classes to
                implement actually setting the temperature of the system.
        '''
        # Set reference structure and property instance
        if temperature_K is None:
            temperature_K = 0
        if cell_cauchy_stress_eV_angstrom3 is None:
            cell_cauchy_stress_eV_angstrom3 = [0, 0, 0, 0, 0, 0]

        result = list(reference_info.values())[0][0]
        self._SingleCrystalTestDriver__nominal_crystal_structure_npt = result
        self._SingleCrystalTestDriver__nominal_crystal_structure_npt["temperature"] = {
            "source-value": temperature_K,
            "source-unit": "K",
        }
        self._SingleCrystalTestDriver__nominal_crystal_structure_npt["cell-cauchy-stress"] = {
            "source-value": cell_cauchy_stress_eV_angstrom3,
            "source-unit": "eV/angstrom^3",
        }
        self._add_property_instance_and_common_crystal_genome_keys(
                "reference-elemental-energy",
                write_stress=True, write_temp=True
                )
        self._add_key_to_current_property_instance(
                "binding-potential-energy-per-atom",
                result["binding-potential-energy-per-atom"]["source-value"],
                result["binding-potential-energy-per-atom"]["source-unit"]
                )
        self._add_key_to_current_property_instance(
                "binding-potential-energy-per-formula",
                result["binding-potential-energy-per-formula"]["source-value"],
                result["binding-potential-energy-per-formula"]["source-unit"]
                )

    def _resolve_dependencies(self, **kwargs):
        '''
        Take reference element and map to reference structure
        Compute EquilibriumCrystalStructure and format binding-energy-crystal result
        '''
        print ("Resolving dependencies...")
        atoms = Atoms.fromdict(reference_map[self.reference_element])
        ecs_test = kimvv.EquilibriumCrystalStructure(self._calc)
        ecs_test(atoms)
        kwargs['reference_info'] = {self.reference_element: [ecs_test.property_instances[1]]}
        return kwargs


if __name__ == "__main__":
    #test = TestDriver('EAM_Dynamo_ErcolessiAdams_1994_Al__MO_123629422045_005')
    test = TestDriver('MEAM_LAMMPS_LeeShimBaskes_2003_Al__MO_353977746962_001')
    test._setup('Al')
    kwargs = test._resolve_dependencies()
    test._calculate(**kwargs)
    test.write_property_instances_to_file()
