# Molecule specification for Hydrogen molecule

from openfermion.hamiltonians import MolecularData

geometry = [['H', [0, 0, 0]],
            ['H', [0, 0, 0.74]]]

basis = 'sto-3g'
multiplicity = 1
charge = 0

h2_molecule = MolecularData(geometry, basis, multiplicity, charge)

print("Current simulation state:\n")
print("Geometry: H [0, 0, 0] - H[0, 0, 0.74]")
print("Basis: Slater-type orbitals, with three gaussians")
print("Spin multiplicity = 1")
print("Charge = 0")

# Second quantization in creation/anhilation operators

print("\nThe molecular Hamiltonian in the second quantization is: H_1 + H_2\n")
print("H = H_1 + H_2 = \u03A3 (h_ij(R)a\u2020_ia_j) + \u03A3 (h_ijkl(R)a\u2020_ia\u2020_ja_k_al)")

integrals = h2_molecule.get_integrals()

one_electron_hamiltonian = integrals[0]
two_electron_hamiltonian = integrals[1]

print("\nH_1 = \n")
h1 = ""
for i in range(2):
    for j in range(2):
        #str_coefficient =  '{:22}'.format( str(one_electron_hamiltonian[i][j]) )
        str_coefficient =  str(one_electron_hamiltonian[i][j])
        #literal_coef = '{:8}'.format("a\u2020" + str(i) + "a"+ str(j))
        literal_coef = "[a\u2020" + str(i) + "a"+ str(j) + "]"
        h1 += str_coefficient  + literal_coef 
        if not (i==1 and j==1):
            h1 += "  + " 

print(h1)

print("\nH_2 = \n")
h2 = ""

for i in range(2):
    for j in range(2):
        for k in range(2):
            for l in range(2):
                str_coefficient =  str(two_electron_hamiltonian[i][j][k][l])
                literal_coef = "[a\u2020" + str(i) + "a\u2020" + str(j) + "a" + str(k) + "a" + str(l) + "]"
                h2 += str_coefficient  + literal_coef 
                if not (i == 1 and j== 1 and k == 1 and l == 1):
                    h2 += " + " 
        h2+= "\n"

print(h2)

print("\nThe latter second quantization has the integral coeficients for the electron integrals in the molecular orbital basis")
print("\n**The molecular hamiltonian's interaction constant, has been manually set to 0 in order to do the following mapping")

# The molecular hamiltonian computed by openfermion sets the interaction constant as "None" type by default when the occupied_indices and active_indices (get_molecular_hamiltonian() arguments) are not specified, which is why is manually set to 0
# If this constant is "none" type instead of "numeric", the next transformation raises an exception
H2_molecular_hamiltonian = h2_molecule.get_molecular_hamiltonian()
H2_molecular_hamiltonian.constant = 0


# Mapping to qubits

from openfermion.transforms import get_fermion_operator, bravyi_kitaev

h2_qubit_hamiltonian = bravyi_kitaev(get_fermion_operator(H2_molecular_hamiltonian))

print("\nThe latter Hamiltoian mapped to a qubit Hamiltonian in through the Bravyi-kitaev transformation is:\n ")
print(h2_qubit_hamiltonian)


# Circuit generation:

import cirq
import openfermioncirq as ofc

qubits = cirq.LineQubit.range(4)
circuit = cirq.Circuit(
        ofc.simulate_trotter(
            qubits,
            H2_molecular_hamiltonian,
            time = 1.0,
            n_steps = 1,
            order = 0,
            algorithm = ofc.trotter.LOW_RANK,
            omit_final_swaps = True
        )
)

cirq.merge_single_qubit_gates_into_phased_x_z(circuit)

print("\nThe corresponding quantum circuit, obtained with cirq:\n")
circuit_cirq_string = circuit.to_text_diagram(use_unicode_characters = True)
circuit_array = circuit_cirq_string.split('\n')

for x in range(4):
    inf_lim = 0 + x*125
    sup_lim = inf_lim + 125
    if sup_lim > 514:
        sup_lim = 514
    for y in range(len(circuit_array)): 
        if x != 0:
            if 0 == y%2:
                print( str(y) + ": ", end = '' )
            else:
                print(3*' ', end = '')
        print(circuit_array[y][inf_lim:sup_lim])
    print("\n\n")




# Optional, it may be done with a psi4 over the h2 molecule:

#Let's reset the original H2_molecular_hamiltonian to it's original form and run psi_4:
H2_molecular_hamiltonian.constant = None

from openfermionpsi4 import run_psi4

h2_molecule_psi4 = run_psi4(h2_molecule,
                            run_mp2 = True,
                            run_cisd = True,
                            run_ccsd = True,
                            run_fci = True)


# Mapping to qubits by Bravyi-Kitaev for the new form of the hamiltonian

print("\n\nNow let's make a psi_4 integral generation and re-map to qubits:\n")
h2_qubit_hamiltonian_psi4 = bravyi_kitaev(get_fermion_operator(h2_molecule_psi4.get_molecular_hamiltonian() ))

print("\nThe latter Hamiltoian mapped to a qubit Hamiltonian in through the Bravyi-kitaev transformation is:\n ")
print(h2_qubit_hamiltonian_psi4)

# circuit generation with psi4

qubits_psi4 = cirq.LineQubit.range(4)
circuit_psi4 = cirq.Circuit(
        ofc.simulate_trotter(
            qubits,
            H2_molecular_hamiltonian,
            time = 1.0,
            n_steps = 1,
            order = 0,
            algorithm = ofc.trotter.LOW_RANK,
            omit_final_swaps = True
        )
)

cirq.merge_single_qubit_gates_into_phased_x_z(circuit)

print("\nThe corresponding quantum circuit with psi4, obtained with cirq:\n")
circuit_cirq_string_psi4 = circuit_psi4.to_text_diagram(use_unicode_characters = True)
circuit_array_psi4 = circuit_cirq_string.split('\n')

for x in range(4):
    inf_lim = 0 + x*125
    sup_lim = inf_lim + 125
    if sup_lim > 514:
        sup_lim = 514
    for y in range(len(circuit_array_psi4)): 
        if x != 0:
            if 0 == y%2:
                print( str(y) + ": ", end = '' )
            else:
                print(3*' ', end = '')
        print(circuit_array_psi4[y][inf_lim:sup_lim])
    print("\n\n")


