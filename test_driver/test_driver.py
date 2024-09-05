# change runner to accept different input values
# change template input

from kim_tools import KIMTestDriver, CrystalGenomeTestDriver, aflow_util, KIMTestDriverError
from ase import Atoms
from typing import Any, Optional, List, Union, Dict, IO
from warnings import warn
from ase.constraints import FixSymmetry, UnitCellFilter
from ase.optimize import BFGS
import numpy as np
# TODO: check if shortname is being generated

class TestDriver(CrystalGenomeTestDriver):

    # Overload definition to accomodate multiple atoms or prototypes
    def _setup(self,
               atoms: Optional[List[Atoms]] = None,
               stoichiometric_species: Optional[List[str]] = None,
               prototype_labels: Optional[List[str]] = None,
               parameter_names: Optional[List[List[str]]] = None,
               parameter_values_angstrom: Optional[List[List[float]]] = None,
               library_prototype_label: Optional[str] = None,
               short_name: Optional[Union[List[str],str]] = None,
               cell_cauchy_stress_eV_angstrom3: List[float] = [0,0,0,0,0,0],
               temperature_K: float = 0,
               crystal_genome_source_structure_id: Optional[List[str]] = None,
               rebuild_atoms: bool = True,
               **kwargs
               ):
        # TODO: Modify below
        """
        Args:
            atoms:
                List of ASE Atoms objects to use as the initial configurations or to build supercells. 
                If this is provided, none of the arguments that are part of the Crystal Genome 
                designation should be provided, and vice versa.
            stoichiometric_species:
                List of unique species in the crystal. Required part of the Crystal Genome designation. 
            prototype_labels:
                List of AFLOW prototype label for the crystal. Required part of the Crystal Genome designation. 
            parameter_names:
                Names of free parameters of the crystal besides 'a'. May be None if the crystal is cubic with no internal DOF.
                Should have length one less than `parameter_values_angstrom`. Part of the Crystal Genome designation.
                May be omitted for debugging, required metadata for OpenKIM Pipeline operation.
            parameter_values_angstrom:
                List of AFLOW prototype parameters for the crystal. Required part of the Crystal Genome designation. 
                a (first element, always present) is in angstroms, while the other parameters 
                (present for crystals with more than 1 DOF) are in degrees or unitless. 
            library_prototype_label: 
                AFLOW library prototype label, may be `None`. Optional part of the Crystal Genome designation.  
            short_name: 
                List of any human-readable short names (e.g. "Face-Centered Cubic") associated with the crystal. 
                Optional part of the Crystal Genome designation.
            cell_cauchy_stress_eV_angstrom3:
                Cauchy stress on the cell in eV/angstrom^3 (ASE units) in [xx,yy,zz,yz,xz,xy] format
            temperature_K:
                The temperature in Kelvin
            crystal_genome_source_structure_id:
                Provenance identifiers of the format '[KIM test result uuid]:[instance-id]'. 
                The chains of dependencies of this test that produced this structure ends in the test listed in the test result, 
                and started with the structure computed in the specific test result and instance-id referenced. May be
                None if this test has no dependencies.
            rebuild_atoms:
                Normally, if you provide an Atoms object, it will be analyzed for its symmetry-reduced AFLOW description,
                and then rebuilt so that the orientation is always consistent. This can rarely cause an error due to
                rounding resulting in a higher-symmetry crystal being created during the rebuild. If you do not care
                about having your Atoms in the standard AFLOW orientation, you can turn the rebuild off.
        """ 
        
        self.atoms_list = atoms
        self.poscar = None
        self.stoichiometric_species = stoichiometric_species        
        self.prototype_labels = prototype_labels
        self.prototype_label = None
        self.parameter_names_list = parameter_names
        self.parameter_values_angstrom_list = parameter_values_angstrom
        self.library_prototype_label_list = library_prototype_label
        self.crystal_genome_source_structure_id = crystal_genome_source_structure_id
        if isinstance(short_name,str):
            self.short_name_list = [short_name]
        else:
            self.short_name_list = short_name
        self.cell_cauchy_stress_eV_angstrom3 = cell_cauchy_stress_eV_angstrom3
        self.temperature_K = temperature_K
        self.rebuild_atoms = rebuild_atoms

        if (
            (self.stoichiometric_species is None) !=
            (self.prototype_labels is None) !=
            (self.parameter_values_angstrom_list is None)
        ):
            print (self._setup.__doc__)
            raise KIMTestDriverError ("\n\nYou have provided some but not all of the required parts of the Crystal Genome designation specified in the docstring above.")

        if self.atoms_list is not None:
            if (
                (self.stoichiometric_species is not None) or # only need to check one required part
                ((self.short_name_list is not None)) or
                ((self.library_prototype_label_list is not None)) or
                ((self.parameter_names_list is not None))
            ):
                print (self._setup.__doc__)
                raise KIMTestDriverError ("\n\nYou have provided an Atoms object as well as at least one part of the Crystal Genome designation specified in the docstring above.\n"
                                          "Please provide only an Atoms object or a Crystal Genome designation, not both")  

        elif self.stoichiometric_species is not None: # we've already checked that if this is not None, other required parts exist as well
            # some checks and cleanup
            print (self.prototype_labels, self.parameter_names_list, self.parameter_values_angstrom_list)
            if (
                (len(self.prototype_labels) != len(self.parameter_names_list)) or
                (len(self.prototype_labels) != len(self.parameter_values_angstrom_list))
            ):
                raise KIMTestDriverError ("\n\nYou have provided a mismatched number of AFLOW prototype information. Please check input lists.")

            for i in range (len(self.prototype_labels)):
                if (len(self.parameter_values_angstrom_list[i]) > 1) and (self.parameter_names_list[i] is None):
                    warn("You've provided parameter values besides `a`, but no parameter names.\n"
                         "Placeholders will be inserted for debugging.")
                    # TODO: Modify below
                    self.parameter_names = ["dummy"]*(len(self.parameter_values_angstrom)-1)

            # Move below elsewhere->move to individual atoms setup function                      
        else:
            warn("You've provided neither a Crystal Genome designation nor an Atoms object.\n"
                     "I won't stop you, but you better know what you're doing!") 

    def _atoms_setup(self, atoms):
        self.atoms = atoms
        self.stoichiometric_species = None # hack to avoid verifying symmetry as it *should* change
        self.prototype_label = None
        self._update_crystal_genome_designation_from_atoms()
        if self.rebuild_atoms:
            # rebuild atoms for consistent orientation
            aflow = aflow_util.AFLOW()
            self.atoms = aflow.build_atoms_from_prototype(self.stoichiometric_species, self.prototype_label, self.parameter_values_angstrom)

    def _aflow_setup(self, prototype_label, parameter_values_angstrom):
        aflow = aflow_util.AFLOW()
        self.atoms = aflow.build_atoms_from_prototype(self.stoichiometric_species, prototype_label, parameter_values_angstrom)
        self._update_poscar() 


    def _relax(self, MAXSTEP = 0.05, FTOL = 1e-8, ITER = 10000, **kwargs):
        self.atoms.calc = self._calc
        symmetry = FixSymmetry(self.atoms)
        self.atoms.set_constraint(symmetry)
        atoms_wrapped = UnitCellFilter(self.atoms)

        # Optimize
        opt = BFGS(atoms_wrapped,maxstep=MAXSTEP)
        try:
            converged = opt.run(fmax=FTOL,steps=ITER)
            iteration_limits_reached = not converged
            minimization_stalled = False
        except:
            minimization_stalled = True
            iteration_limits_reached = False               
        forces = self.atoms.get_forces()
        stress = self.atoms.get_stress()        
        # Compute the average energy per atom after subtracting out the energies of the
        # isolated atoms
        energy_isolated = sum(
            [self.get_isolated_energy_per_atom(sym) for sym in self.atoms.get_chemical_symbols()]
        )
        energy_per_atom = (self.atoms.get_potential_energy() - energy_isolated) / self.atoms.get_global_number_of_atoms()

        print("Minimization "+
              ("converged" if not minimization_stalled else "stalled")+
              " after "+
              (("hitting the maximum of "+str(ITER)) if iteration_limits_reached else str(opt.nsteps))+
              " steps.")
        print("Maximum force component: "+str(np.max(np.abs(forces)))+" eV/Angstrom")
        print("Maximum stress component: "+str(np.max(np.abs(stress)))+" eV/Angstrom^3")
        print("==== Minimized structure obtained from ASE ====")
        print("symbols = ", self.atoms.get_chemical_symbols())
        print("basis = ", self.atoms.get_scaled_positions())
        print("cellpar = ", self.atoms.get_cell())
        print("forces = ", forces)
        print("stress = ", stress)
        print("energy per atom = ", energy_per_atom)
        print("===============================================")
        results = {}
        results['energy'] = energy_per_atom
        results['stress'] = stress
        results['forces'] = forces
        results['atoms'] = self.atoms
        return results 

        
    def _calculate(self, **kwargs):
        #TODO: Modify below
        """
        Relaxes multiple prototype structures and records the lowest energy structure 
        
        Args:
            max_volume_scale:
                Maximum fractional change in volume to investigate
            num_steps:
                Number of steps to take in each direction
        """

        lowest_energy = 1e10
        if self.atoms_list is not None:
            for atoms in self.atoms_list:
                self._atoms_setup(atoms)
                results = self._relax(**kwargs)
                if results['energy'] < lowest_energy:
                    print (results)
                    lowest_energy = results['energy']
                    lowest_energy_results = results
        else:
            for i in range (len(self.prototype_labels)):
                self._aflow_setup(self.prototype_labels[i], self.parameter_values_angstrom_list[i])
                results = self._relax(**kwargs)
                if results['energy'] < lowest_energy:
                    lowest_energy = results['energy']
                    lowest_energy_results = results


        self.stoichiometric_species = None # hack to avoid verifying symmetry as it *should* change
        self.prototype_label = None

        self._update_crystal_genome_designation_from_atoms(lowest_energy_results['atoms'])
        
        self._add_property_instance_and_common_crystal_genome_keys("binding-energy-crystal",
                                                                   write_stress=False, write_temp=False)

        self._add_key_to_current_property_instance("binding-potential-energy-per-atom",lowest_energy,
                                                   units="eV")
        # Will always be same for these prototypes
        self._add_key_to_current_property_instance("binding-potential-energy-per-formula",lowest_energy,
                                                   units="eV")

        self._add_property_instance_and_common_crystal_genome_keys("reference-elemental-ground-state-npt",
                                                                   write_stress=True, write_temp=True)

    def __call__(self, atoms: Optional[Atoms] = None,  **kwargs):
        """
        runs test
        """
        self._setup(atoms, **kwargs)
        self._calculate(**kwargs)

    def get_isolated_energy_per_atom(self, symbol):
        """
        Construct a non-periodic cell containing a single atom and compute its energy.
        """
        single_atom = Atoms(
            symbol,
            positions=[(0.1, 0.1, 0.1)],
            cell=(20, 20, 20),
            pbc=(False, False, False),
        )
        single_atom.calc=self._calc
        energy_per_atom = single_atom.get_potential_energy()
        return energy_per_atom

if __name__ == "__main__":
    from ase.build import bulk
    test = TestDriver('EAM_Dynamo_ZhouWadleyJohnson_2001_Al__MO_049243498555_000')
    atoms = bulk('Al','fcc',a=4.04)
    atoms2 = bulk('Al','bcc',a=4.04)
    #test([atoms,atoms2])
    d = {"stoichiometric_species": ["Al"], "prototype_labels": ["A_cF4_225_a", "A_cI2_229_a"], "parameter_names": [["a"], ["a"]], "parameter_values_angstrom": [[4.039], [3.233]]}
    test(**d) 
    test.write_property_instances_to_file()
