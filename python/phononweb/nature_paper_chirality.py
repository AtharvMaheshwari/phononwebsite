"""
Phonon Band Structure Calculation with Phonon Angular Momentum(AM), Pseudo-angular Momentum(PAM) and helicity

Python Version: 3.12.2
Phonopy Version: 2.38.2

Input Files:
    1. POSCAR-unitcell: structure file
    2. phonopy.yaml: Phonopy configuration
    3. FORCE_SETS: Displacement-force dataset
    4. BORN: Born effective charges (optional, mandatory for phonon magnetic moment calculation)

Output Files:
    1. band.dat: Phonon bands with properties (energy, angular momentum, PAM, helicity, group velocities)
    2. kpoints.dat: k-point coordinates
    3. band.png: Visualization of band structure with phonon angular momentum, helicity and PAM
    4. kpaths_nonzero_AM_PAM.json: Angular momentum and PAM allowed path
    5. Magband.dat: Phonon bands with phonon magnetic moment
    6. Magband.png: Visualization of band structure with phonon magnetic momentum

Calculations:
    1. Phonon band structures
    2. Phonon angular momentum components(Jx/Jy/Jz)
    3. Phonon pseudo-angular momentum (PAM)
    4. Group velocities
    5. Helicity 
    6. Visualizes results with specialized plotting
    7. Phonon magnetic moment
    8. gyromagnetic ratios of each atom
"""

###################################################################
# Developers: Yue Yang, Yu Mao, Zhanghuan Li, and Yuanfeng Xu
# Contact: Yuanfeng Xu, Zhejiang University
# Email: y.xu@zju.edu.cn
# Reference: Yue Yang et al., Catalogue of Chiral Phonon Materials.
###################################################################

import os
import numpy as np
import phonopy
import matplotlib.pyplot as plt
import cmath
import json
from phonopy import Phonopy
from phonopy.interface.vasp import read_vasp
from phonopy.file_IO import parse_FORCE_SETS
from phonopy.structure.cells import get_primitive
from library import * 
from multiprocessing import Pool, cpu_count
from matplotlib.colors import Normalize
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import SymLogNorm


# ----------------------
# Global Constants
# ----------------------
nk_p = 40                   # Points per path segment in k-space
sqrt2 = np.sqrt(2)          # Constant for chirality projection operators
ctol = 0.0                  # Tolerance threshold
pi2 = np.pi * 2.0           # 2π constant
pi = np.pi                  # π constant
factor = 64.6541380         # Unit conversion factor (VASP to meV)
H_BAR = 1.0545718e-34       # reduced Planck constant (J·s)
E_CHARGE = 1.60217662e-19   # Coulomb
U_TO_KG = 1.66053906660e-27 # atom weight (u → kg)
MU_NUCLEAR = 5.0507837e-27  # Nuclear magneton (J/T)
MU_BOHR = 9.274009994e-24   # Bohr magneton (J/T) 


def compute_AM(ipath: int, kpath: list, 
                     Jz_path: list, Jx_path: list, Jy_path: list, 
                     eigs: np.ndarray, natom: int, ctol: float, 
                     chirality_band: np.ndarray, chirality_atom_band: np.ndarray,) -> None:
    """
    Compute AM components for specific k-point
    
    Args:
        ipath: Current path segment index
        kpath: Path labels (e.g., ['Γ-X', 'X-M'])
        Jx_path: Paths where Jx is active
        Jy_path: Paths where Jy is active
        Jz_path: Paths where Jz is active
        eigs: Eigenvectors at current k-point
        natom: Number of atoms
        ctol: AM components tolerance
        chirality_band: AM components storage array (updated in-place)
        chirality_atom_band: AM components of each atom
    """
    # XY chirality (Jz component)
    if kpath[ipath] in Jz_path:
        cmat = np.zeros(3 * natom, dtype=np.complex128)
        for ibnd in range(3 * natom):
            for iat in range(natom):
                opr_l, opr_r = chirality_opr_xy(natom, iat)
                tmpl = np.dot(np.conj(eigs[:, ibnd]).T, opr_l)
                tmpr = np.dot(np.conj(eigs[:, ibnd]).T, opr_r)
                amat = tmpr[0] * np.conj(tmpr[0]) - tmpl[0] * np.conj(tmpl[0])
                if abs(np.real(amat)) > ctol:
                    chirality_atom_band[iat, ibnd, 2] = np.real(amat)
                cmat[ibnd] += amat
        for ib in range(3 * natom):
            if abs(np.real(cmat[ib])) > ctol:
                chirality_band[ib, 0] = np.real(cmat[ib])
    else:
        chirality_band[:, 0] = 0.0
        chirality_atom_band[:, :, 2] = 0.0
    
    # YZ chirality (Jx component)
    if kpath[ipath] in Jx_path:
        cmat = np.zeros(3 * natom, dtype=np.complex128)
        for ibnd in range(3 * natom):
            for iat in range(natom):
                opr_l, opr_r = chirality_opr_yz(natom, iat)
                tmpl = np.dot(np.conj(eigs[:, ibnd]).T, opr_l)
                tmpr = np.dot(np.conj(eigs[:, ibnd]).T, opr_r)
                amat = tmpr[0] * np.conj(tmpr[0]) - tmpl[0] * np.conj(tmpl[0])
                if abs(np.real(amat)) > ctol:
                    chirality_atom_band[iat, ibnd, 0] = np.real(amat)
                cmat[ibnd] += amat
        for ib in range(3 * natom):
            if abs(np.real(cmat[ib])) > ctol:
                chirality_band[ib, 1] = np.real(cmat[ib])
    else:
        chirality_band[:, 1] = 0.0
        chirality_atom_band[:, :, 0] = 0.0
    
    # ZX chirality (Jy component)
    if kpath[ipath] in Jy_path:
        cmat = np.zeros(3 * natom, dtype=np.complex128)
        for ibnd in range(3 * natom):
            for iat in range(natom):
                opr_l, opr_r = chirality_opr_zx(natom, iat)
                tmpl = np.dot(np.conj(eigs[:, ibnd]).T, opr_l)
                tmpr = np.dot(np.conj(eigs[:, ibnd]).T, opr_r)
                amat = tmpr[0] * np.conj(tmpr[0]) - tmpl[0] * np.conj(tmpl[0])
                if abs(np.real(amat)) > ctol:
                    chirality_atom_band[iat, ibnd, 1] = np.real(amat)
                cmat[ibnd] += amat
        for ib in range(3 * natom):
            if abs(np.real(cmat[ib])) > ctol:
                chirality_band[ib, 2] = np.real(cmat[ib])
    else:
        chirality_band[:, 2] = 0.0
        chirality_atom_band[:, :, 1] = 0.0

def compute_PAM(qpoint: np.ndarray, natom: int, symm_no: int, 
               rotation: list, translation: list, rprim: np.ndarray, 
               positions: np.ndarray, eigs: np.ndarray, rot_type: str) -> np.ndarray:
    """
    Compute Phonon Pseudo-angular Momentum (PAM)
    
    Args:
        qpoint: Current k-point (fractional)
        natom: Number of atoms
        symm_no: Symmetry operation index
        rotation: List of rotation matrices
        translation: List of translation vectors
        rprim: Primitive cell vectors
        positions: Atomic positions
        eigs: Eigenvectors
        rot_type: Rotation type ('C3' or 'C4')
    
    Returns:
        PAM tensor: [spin, orbital, total] x [mode]
    """
    # Get modified permutation matrices
    DR, DR_s, DR_o = new_get_modified_permutation_matrix(
        qpoint, natom, symm_no, rotation, translation, rprim, positions
    ) 
    
    # Compute trace operations (matrix projections)
    trace_tot = np.dot(np.conj(eigs).T, np.dot(DR, eigs))
    trace_s = np.dot(np.conj(eigs).T, np.dot(DR_s, eigs))
    trace_o = np.dot(np.conj(eigs).T, np.dot(DR_o, eigs))
    
    # Initialize PAM storage
    PAM = np.zeros((3, 3 * natom), dtype=np.float64)
    
    # Set rotation-specific parameters
    if rot_type == 'C3':
        factor_val = -3 / pi2
        mod_val = 3
    elif rot_type == 'C4':
        factor_val = -4 / pi2
        mod_val = 4
    else:
        raise ValueError(f"Invalid rotation type: {rot_type}")
    
    # Compute PAM for each mode
    for ib in range(3 * natom):
        PAM[0, ib] = factor_val * cmath.phase(trace_s[ib, ib])  # Spin PAM
        PAM[1, ib] = factor_val * cmath.phase(trace_o[ib, ib])   # Orbital PAM
        total = factor_val * cmath.phase(trace_tot[ib, ib])      # Total PAM
        PAM[2, ib] = total % mod_val                         # Wrap to [0, mod)
        if PAM[2, ib] < 0:
            PAM[2, ib] += mod_val
        if abs(PAM[2, ib] - mod_val) < 1e-5:  # Near-modulus threshold
            PAM[2, ib] = 0.0
    
    return PAM

# ----------------------
# Parallel K-point Processing
# ----------------------
def process_kpoint(args: tuple) -> tuple:
    """
    Process single k-point: Compute phonon properties
    
    Args:
        args: Packed arguments (see main program)
    
    Returns:
        (index, energy, group velocity, AM, PAM)
    """
    # Unpack arguments
    ipath, ik, nk_p, kpath, Jz_path, Jx_path, Jy_path, symk, bvec, phonon, natom, ctol, \
    rotation_type, rotation_path_num, symm_no_r3, symm_no_r4, TRS, tag_m, tag_r3, tag_r4, \
    rotation, translation, rprim, positions, gamma, option_nac = args
    
    ikk = ipath * nk_p + ik  # Global k-index
    
    # Interpolate k-point along segment
    kpts1 = symk[ipath, :, 0]
    kpts2 = symk[ipath, :, 1]
    qpoint = kpts1 + ik * (kpts2 - kpts1) / (nk_p - 1)
    
    # Compute phonon properties at qpoint
    dmat = phonon.get_dynamical_matrix_at_q(qpoint)
    group_velocity = phonon.get_group_velocity_at_q(qpoint)
    ek, eigs = np.linalg.eigh(dmat)
    enk = (np.sqrt(np.abs(ek)) * np.sign(ek)) * factor  # Convert to meV
    
    # Compute chirality components
    chirality = np.zeros((3 * natom, 3), dtype=np.float64)
    chirality_atom = np.zeros((natom, 3 * natom, 3), dtype=np.float64)
    compute_AM(ipath, kpath, Jz_path, Jx_path, Jy_path, eigs, natom, ctol, chirality, chirality_atom)
    
    # Compute PAM if applicable
    PAM = np.zeros((3, 3 * natom), dtype=np.float64)
    if ipath in rotation_path_num:
        rot_idx = rotation_path_num.index(ipath)
        rot_type = rotation_type[rot_idx]
        if TRS == 0 and tag_m != 0:  # Time-reversal broken
            if rot_type == 'C3' and tag_r3 == 1:
                PAM = compute_PAM(qpoint, natom, symm_no_r3, rotation, translation, rprim, positions, eigs, rot_type)
            elif rot_type == 'C4' and tag_r4 == 1:
                PAM = compute_PAM(qpoint, natom, symm_no_r4, rotation, translation, rprim, positions, eigs, rot_type)
    
    # Compute phonon magnetic moments 
    M = np.zeros((3, 3 * natom), dtype=float)
    if option_nac:
        for i in range(3 * natom):
            for j in range(natom):
                gamma_contribution = gamma[j].diagonal(0) * chirality_atom[j, i, :].T
                M[:, i] += gamma_contribution * H_BAR / MU_NUCLEAR

    return (ikk, enk, group_velocity, chirality, PAM, M)

def get_color_and_size(pam_value):
    """Determine color and marker size based on PAM value
    
    Args:
        PAM_value: Phonon pseudo-angular momentum value
    
    Returns:
        (color, size): Color category and scaled size
    """
    if 0 <= pam_value <= 1:
        return 'orange', pam_value
    elif 1 < pam_value <= 2:
        return 'blue', pam_value - 1
    elif 2 < pam_value <= 3:
        return 'red', pam_value - 2
    elif 3 < pam_value < 4:
        return 'green', pam_value - 3
    else:
        return 'black', 0

def plot_component(ax, data, title,ylab):
    for jj in range(3 * natom):
        ax.plot(lenk, enk_band[jj], color='black', linewidth=0.4)
        point_colors = [cmap(norm(data[jj, ikk])) for ikk in range(nk_total)]   
        ax.scatter(lenk, enk_band[jj], s=fixed_size, c=point_colors, 
                marker='o', edgecolors='none')
    for ii in range(nkp_path-1):
        ax.plot( [lenk[nk_p*(ii+1)-1], lenk[nk_p*(ii+1)-1]], [emin, emax], color='black')  
    ax.set_xlim(lenk[0], lenk[nk_p * nkp_path - 1])
    ax.set_ylim(emin, emax)
    ind=[0.0]
    for ii in range(nkp_path):
        ind.append(lenk[(ii+1)*nk_p-1])
    ax.set_xticks(ind)
    ax.set_xticklabels(symbol)
    if ylab ==1:
        ax.set_ylabel('Frequency (meV)')
    ax.set_title(title)

# ----------------------
# Main Program
# ----------------------
if __name__ == '__main__':
    # ----------------------
    # Read Input Files & Setup
    # ----------------------
    os.chdir("/root/capsule/code/ChiralPY/input")
    option_nac = os.path.isfile("BORN")  # Non-analytical correction flag
    cell = read_vasp("POSCAR-unitcell")  # Unit cell
    
    # Parse configuration files
    phononunit = Phonopy(
        unitcell=cell,
        supercell_matrix=np.eye(3),  # Using unit supercell matrix
        primitive_matrix='auto'      # Auto-generate primitive cell matrix
    )
    primitive_matrix = phononunit.primitive_matrix
    
    # Initialize phonopy object
    phonon = phonopy.load(
        phonopy_yaml='phonopy.yaml',
        primitive_matrix=primitive_matrix,
        is_nac=option_nac,
        log_level=0
    )
    phonon.save()
    # Get primitive cell information
    cell1 = get_primitive(cell, primitive_matrix, symprec=1e-4)
    rprim = cell1.cell
    positions = cell1.scaled_positions
    natom, _ = positions.shape
    mass = cell1.masses
    
    # Get NAC information and calculate gyromagnetic ratios of each atom
    gamma=[[] for i in range(natom)]
    if option_nac:
        born_charges, dielectric_matrix = parse_nac_data("phonopy_params.yaml")
        for i in range(natom):
            scalar = E_CHARGE / (2 * mass[i])/U_TO_KG
            gamma[i] = np.multiply(born_charges[i], scalar)

    # ----------------------
    # Symmetry Analysis
    # ----------------------
    symmetry = phonon.symmetry
    prim_symm = phonon.primitive_symmetry
    rotation = prim_symm.symmetry_operations['rotations']
    translation = prim_symm.symmetry_operations['translations']
    nsym = len(rotation)
    
    # Get space group information
    spg_info = symmetry.get_international_table().split()
    spg_symbol = spg_info[0]
    spg_number = int(spg_info[1].strip('()'))
    print(f"Space group: {spg_symbol} ({spg_number})")
    
    # Determine chirality and PAM paths
    Jz_path, Jx_path, Jy_path, ls_path, PAM_rotation_type, PAM_rotation_axis, PAM_path_num = cal_sp(spg_number, option_nac)
    data_paths_AM_PAM = {
    "Jx": Jx_path,
    "Jy": Jy_path,
    "Jz": Jz_path,
    "PAM": ls_path
    }
    
    # ----------------------
    # Force Constants & Path Setup
    # ----------------------
    force_sets = parse_FORCE_SETS()
    phonon.dataset = force_sets
    phonon.produce_force_constants()
    
    # Generate k-path in reciprocal space
    bvec = rec_latt(rprim)  # Reciprocal lattice vectors
    nkp_path, symk, symbol, lenk, kpath = hsym_path(spg_number, bvec)
    nk_total = nk_p * nkp_path  # Total k-points
    
    # ----------------------
    # Initialize Storage Arrays
    # ----------------------
    qpts_band = np.zeros((nk_total, 3))  # Fractional k-points
    kpts_band = np.zeros((nk_total, 3))  # Cartesian k-points
    enk_band = np.zeros((3 * natom, nk_total))  # Phonon energies
    chirality_band = np.zeros((3 * natom, nk_total, 3))  # AM components
    helicity_band = np.zeros((3 * natom, nk_total))  # Helicity
    group_velocity_band = np.zeros((3 * natom, nk_total, 3))  # Group velocities
    PAM_band = np.zeros((3, 3 * natom, nk_total))  # PAM
    Mag_band =  np.zeros((3, 3 * natom, nk_total)) #Phonon magnetic moments
    Mag_abs_band =  np.zeros((3 * natom, nk_total)) #Module of phonon magnetic moments
    P_vec_band = np.zeros((4, 3 * natom, nk_total)) # Cycloidicity

    # ----------------------
    # Symmetry Precomputation
    # ----------------------
    symm_no_r3_list = []
    symm_no_r4_list = []
    TRS_list = []
    tag_m_list = []
    tag_r3_list = []
    tag_r4_list = []
    high_symm_pt = np.zeros((2 * nkp_path, 3))  # High-symmetry points
    
    # ----------------------
    # Parallel Task Preparation
    # ----------------------
    tasks = []
    for ipath in range(nkp_path):
        kpts1 = symk[ipath, :, 0]
        kpts2 = symk[ipath, :, 1]
        kpts3 = symk[ipath, :, 2]
        # Compute symmetry properties at path center
        tag_m, _, tag_r3, tag_r4, _, _, symm_no_r3, symm_no_r4, _, _, _, _, _, _, _, _, _ = kpt_pg(kpts3, nsym, rotation, rprim)
        TRS = time_rev(kpts3)
        
        # Store symmetry data
        high_symm_pt[2*ipath] = kpts1
        high_symm_pt[2*ipath+1] = kpts2
        TRS_list.append(TRS)
        tag_m_list.append(tag_m)
        tag_r3_list.append(tag_r3)
        tag_r4_list.append(tag_r4)
        symm_no_r3_list.append(symm_no_r3)
        symm_no_r4_list.append(symm_no_r4)

        # Prepare tasks for each k-point in segment
        for ik in range(nk_p):
            qpoint = kpts1 + ik * (kpts2 - kpts1) / (nk_p - 1)
            qpts_band[ipath*nk_p + ik] = qpoint
            kpts_band[ipath*nk_p + ik] = qpoint @ bvec.T  # Convert to Cartesian
            
            args = (
                ipath, ik, nk_p, kpath, Jz_path, Jx_path, Jy_path,
                symk, bvec, phonon, natom, ctol,
                PAM_rotation_type, PAM_path_num, 
                symm_no_r3_list[ipath], symm_no_r4_list[ipath],
                TRS_list[ipath], tag_m_list[ipath], 
                tag_r3_list[ipath], tag_r4_list[ipath],
                rotation, translation, rprim, positions, gamma, option_nac
            )
            tasks.append(args)

    # ----------------------
    # Parallel Computation
    # ----------------------
    num_cores = cpu_count()
    print(f"Processing {nk_total} k-points using {num_cores} cores")
    with Pool(num_cores) as pool:
        results = pool.map(process_kpoint, tasks)
    
    # Store results
    for ikk, enk, group_velocity, chirality, PAM, M in results:
        enk_band[:, ikk] = enk
        group_velocity_band[:, ikk] = group_velocity
        chirality_band[:, ikk] = chirality
        PAM_band[:, :, ikk] = PAM
        Mag_band[:, :, ikk] = M

    # ----------------------
    # Post-Processing: Degeneracy Handling
    # ----------------------
    # For Jz-active paths
    for ipath in range(nkp_path):
        if kpath[ipath] in Jz_path:
            for jj in range(3 * natom - 1):
                e_diff = np.mean(np.abs(enk_band[jj+1, ipath*nk_p:(ipath+1)*nk_p] - enk_band[jj, ipath*nk_p:(ipath+1)*nk_p]))
                if e_diff < 1e-3:  # Degeneracy threshold
                    # Average value for degenerate modes
                    avg = (chirality_band[jj+1, ipath*nk_p:(ipath+1)*nk_p, 0] + 
                           chirality_band[jj, ipath*nk_p:(ipath+1)*nk_p, 0]) / 2
                    chirality_band[jj+1, ipath*nk_p:(ipath+1)*nk_p, 0] = avg
                    chirality_band[jj, ipath*nk_p:(ipath+1)*nk_p, 0] = avg
    # For Jx-active paths
    for ipath in range(nkp_path):
        if kpath[ipath] in Jx_path:
            for jj in range(3 * natom - 1):
                e_diff = np.mean(np.abs(enk_band[jj+1, ipath*nk_p:(ipath+1)*nk_p] - enk_band[jj, ipath*nk_p:(ipath+1)*nk_p]))
                if e_diff < 1e-3:  # Degeneracy threshold
                    # Average value for degenerate modes
                    avg = (chirality_band[jj+1, ipath*nk_p:(ipath+1)*nk_p, 1] + 
                           chirality_band[jj, ipath*nk_p:(ipath+1)*nk_p, 1]) / 2
                    chirality_band[jj+1, ipath*nk_p:(ipath+1)*nk_p, 1] = avg
                    chirality_band[jj, ipath*nk_p:(ipath+1)*nk_p, 1] = avg
    # For Jy-active paths
    for ipath in range(nkp_path):
        if kpath[ipath] in Jy_path:
            for jj in range(3 * natom - 1):
                e_diff = np.mean(np.abs(enk_band[jj+1, ipath*nk_p:(ipath+1)*nk_p] - enk_band[jj, ipath*nk_p:(ipath+1)*nk_p]))
                if e_diff < 1e-3:  # Degeneracy threshold
                    # Average value for degenerate modes
                    avg = (chirality_band[jj+1, ipath*nk_p:(ipath+1)*nk_p, 2] + 
                           chirality_band[jj, ipath*nk_p:(ipath+1)*nk_p, 2]) / 2
                    chirality_band[jj+1, ipath*nk_p:(ipath+1)*nk_p, 2] = avg
                    chirality_band[jj, ipath*nk_p:(ipath+1)*nk_p, 2] = avg
    
    # For PAM-active paths
    for ipath in PAM_path_num:
        for jj in range(3 * natom - 1):
            e_diff = np.mean(np.abs(enk_band[jj+1, ipath*nk_p:(ipath+1)*nk_p] - enk_band[jj, ipath*nk_p:(ipath+1)*nk_p]))
            if e_diff < 1e-3:
                avg = (PAM_band[2, jj+1, ipath*nk_p:(ipath+1)*nk_p] + 
                       PAM_band[2, jj, ipath*nk_p:(ipath+1)*nk_p]) / 2
                PAM_band[2, jj+1, ipath*nk_p:(ipath+1)*nk_p] = avg
                PAM_band[2, jj, ipath*nk_p:(ipath+1)*nk_p] = avg

    # For phonon magnetic moment
    if option_nac:
        for ipath in range(nkp_path):
            for jj in range(3 * natom - 1):
                e_diff = np.mean(np.abs(enk_band[jj+1, ipath*nk_p:(ipath+1)*nk_p] - enk_band[jj, ipath*nk_p:(ipath+1)*nk_p]))
                if e_diff < 1e-3:  # Degeneracy threshold
                    # Average value for degenerate modes
                    avgx = (Mag_band[0, jj+1, ipath*nk_p:(ipath+1)*nk_p] +  Mag_band[0, jj, ipath*nk_p:(ipath+1)*nk_p]) / 2
                    avgy = (Mag_band[1, jj+1, ipath*nk_p:(ipath+1)*nk_p] +  Mag_band[1, jj, ipath*nk_p:(ipath+1)*nk_p]) / 2
                    avgz = (Mag_band[2, jj+1, ipath*nk_p:(ipath+1)*nk_p] +  Mag_band[2, jj, ipath*nk_p:(ipath+1)*nk_p]) / 2
                    Mag_band[0, jj+1, ipath*nk_p:(ipath+1)*nk_p] = Mag_band[0, jj, ipath*nk_p:(ipath+1)*nk_p] = avgx
                    Mag_band[1, jj+1, ipath*nk_p:(ipath+1)*nk_p] = Mag_band[1, jj, ipath*nk_p:(ipath+1)*nk_p] = avgy
                    Mag_band[2, jj+1, ipath*nk_p:(ipath+1)*nk_p] = Mag_band[2, jj, ipath*nk_p:(ipath+1)*nk_p] = avgz

    # Module of phonon magnetic moment
    for ik in range(nk_total):
        for i in range(3 * natom):
                Mag_abs_band[i,ik]=np.sqrt(Mag_band[0,i,ik]**2+Mag_band[1,i,ik]**2+Mag_band[2,i,ik]**2)

    # ----------------------
    # Process High-Symmetry Points
    # ----------------------
    for ik, hp in enumerate(high_symm_pt):
        ipath = ik // 2
        ik_local = 0 if (ik % 2 == 0) else (nk_p - 1)
        ikk_global = ipath * nk_p + ik_local
        
        # Recompute at high-symmetry points
        dmat = phonon.get_dynamical_matrix_at_q(hp)
        ek, eigs = np.linalg.eigh(dmat)
        
        # Compute PAM at HSPs
        tag_m, _, tag_r3, tag_r4, _, _, symm_no_r3, symm_no_r4, _, _, _, _, _, _, _, _, _ = kpt_pg(hp, nsym, rotation, rprim)
        TRS = time_rev(hp)
        
        temp_PAM = np.zeros((3, 3 * natom), dtype=np.float64)
        if TRS == 0 and tag_m != 0:  # Time-reversal broken
            if tag_r3 == 1:  # C3 symmetry
                temp_PAM = compute_PAM(hp, natom, symm_no_r3, rotation, translation, rprim, positions, eigs, 'C3')
            elif tag_r4 == 1:  # C4 symmetry
                temp_PAM = compute_PAM(hp, natom, symm_no_r4, rotation, translation, rprim, positions, eigs, 'C4')
        PAM_band[:, :, ikk_global] = temp_PAM
        
        # Compute AM at HSPs
        cmat = np.zeros((3 * natom, 3), dtype=np.complex128)
        for ibnd in range(3 * natom):
            for iat in range(natom):
                # XY plane
                opr_lxy, opr_rxy = chirality_opr_xy(natom, iat)
                tmplxy = np.conj(eigs[:, ibnd]).T @ opr_lxy
                tmprxy = np.conj(eigs[:, ibnd]).T @ opr_rxy
                # YZ plane
                opr_lyz, opr_ryz = chirality_opr_yz(natom, iat)
                tmplyz = np.conj(eigs[:, ibnd]).T @ opr_lyz
                tmpryz = np.conj(eigs[:, ibnd]).T @ opr_ryz
                # ZX plane
                opr_lzx, opr_rzx = chirality_opr_zx(natom, iat)
                tmplzx = np.conj(eigs[:, ibnd]).T @ opr_lzx
                tmprzx = np.conj(eigs[:, ibnd]).T @ opr_rzx
                
                # Compute differences
                cmat[ibnd, 0] += tmprxy[0] * tmprxy[0].conj() - tmplxy[0] * tmplxy[0].conj()
                cmat[ibnd, 1] += tmpryz[0] * tmpryz[0].conj() - tmplyz[0] * tmplyz[0].conj()
                cmat[ibnd, 2] += tmprzx[0] * tmprzx[0].conj() - tmplzx[0] * tmplzx[0].conj()
        
        # Update AM array
        idkk = (ik+1)//2 * nk_p - (ik % 2)
        for ib in range(3 * natom):
            for comp in range(3):
                val = np.real(cmat[ib, comp])
                if abs(val) > ctol:
                    chirality_band[ib, idkk, comp] = val

    # ----------------------
    # Helicity Calculation
    # ----------------------
    list_path_zero_helicity = []  # Paths with zero helicity
    for ipath in range(nkp_path):
        kpts3 = symk[ipath, :, 2]
        lg_nsym, lg_symm_list = little_group(kpts3, nsym, rotation)
        for iii in lg_symm_list:
            lg_matrix = similarity_transformation(np.transpose(rprim), rotation[int(iii)])
            if abs(np.linalg.det(lg_matrix) + 1) < 1e-5:  # Mirror symmetry
                list_path_zero_helicity.append(ipath)
                
    # Compute helicity for non-mirror paths
    for ipath in range(nkp_path):
        if ipath not in list_path_zero_helicity:
            for ik in range(nk_p):
                ikk = ipath * nk_p + ik
                qpoint_cart = kpts_band[ikk]
                if np.linalg.norm(qpoint_cart) != 0:
                    for ib in range(3 * natom):
                        J = chirality_band[ib, ikk, [1, 2, 0]]  # AM vector
                        helicity_band[ib, ikk] = np.dot(qpoint_cart, J) / np.linalg.norm(qpoint_cart)
    
    # ----------------------
    # Cycloidicity Calculation
    # ----------------------
    P_vec_band = np.zeros((4, 3 * natom, nk_total)) 
    for ib in range(3 * natom):
        for ik in range(nk_total):
            tempJ = np.array([chirality_band[ib, ik, 1],chirality_band[ib, ik, 2],chirality_band[ib, ik, 0]])
            tempk = np.array([kpts_band[ik,0],kpts_band[ik,1],kpts_band[ik,2]])
            k_norm = np.linalg.norm(tempk)
            if k_norm > 1e-12:
                k_unit = tempk / k_norm
            else:
                k_unit = np.array([0, 0, 0])
            P_vec_band[:3,ib,ik]=np.cross(k_unit,tempJ)
            P_vec_band[3,ib,ik] = np.linalg.norm(P_vec_band[:3,ib,ik])

    # ----------------------
    # Write Output Files
    # ----------------------
    # Band data file
    with open('/root/capsule/results/band.dat', 'w') as fenk:
        # Header with PAM information
        first_line = '# PAM:'
        for i in range(len(PAM_path_num)):
            pam_path = kpath[PAM_path_num[i]]
            rot_type = PAM_rotation_type[i]
            axis = PAM_rotation_axis[i]
            first_line += f" {pam_path}({PAM_path_num[i]+1}) {rot_type} [{int(axis[0])},{int(axis[1])},{int(axis[2])}],"
        first_line = first_line.rstrip(',')
        fenk.write(f"{first_line}\n")
        
        # Column headers
        line_enk='{:18s}{:14s}{:11s}{:10s}{:16s}{:16s}{:16s}{:16s}{:16s}{:16s}{:14s}{:6s}'.format('#   lenk','enk(meV)','l_s','l_o','l_ph','v_x','v_y','v_z','J_x','J_y','J_z','Helicity')
        fenk.write(line_enk+'\n')
        
        # Write data for each mode and k-point
        for ib in range(3 * natom):
            for ik in range(nk_total):
                line = f"{lenk[ik]:10.7f} {enk_band[ib, ik]:15.7f} {PAM_band[0, ib, ik]:10.5f} {PAM_band[1, ib, ik]:10.5f} {PAM_band[2, ib, ik]:10.5f} "
                line += f"{group_velocity_band[ib, ik, 0]:15.7f} {group_velocity_band[ib, ik, 1]:15.7f} {group_velocity_band[ib, ik, 2]:15.7f} "
                line += f"{chirality_band[ib, ik, 1]:15.7f} {chirality_band[ib, ik, 2]:15.7f} {chirality_band[ib, ik, 0]:15.7f} "
                line += f"{helicity_band[ib, ik]:15.7f}\n"
                fenk.write(line)
            fenk.write("\n\n")
    
    # Write k-point coordinates (kpoints.dat)
    with open('/root/capsule/results/kpoints.dat', 'w') as fk:
        fk.write("#(k1,k2,k3) and (kx,ky,kz) represent direct and Cartesian coordinates of k points.\n")
        line_k='{:10s}{:15s}{:15s}{:15s}{:15s}{:15s}{:15s}'.format('#','k1','k2','k3','kx','ky','kz')
        fk.write(line_k+'\n')
        for ikk in range(nk_total):
            line_k='{:15.7f}{:15.7f}{:15.7f}{:15.7f}{:15.7f}{:15.7f}'.format(qpts_band[ikk,0],qpts_band[ikk,1],qpts_band[ikk,2],kpts_band[ikk,0],kpts_band[ikk,1],kpts_band[ikk,2])
            fk.write(line_k)
            fk.write('\n')

    # Write Cycloidicity data file (Cycloidicity.dat)
    with open('/root/capsule/results/Cycloidicity.dat', 'w') as f:
        # File header with column descriptions
        f.write("# P:phonon cycloidicity vector\n")
        line_enk='{:19s}{:16s}{:14s}{:14s}{:14s}{:16s}'.format('#   lenk','enk(meV)','P_x','P_y','P_z','|P|')
        f.write(line_enk+'\n')
        
        for ib in range(3*natom):
            for ik in range(nk_total):
                line = f" {lenk[ik]:10.7f} {enk_band[ib, ik]:15.7f} {P_vec_band[0,ib,ik]:13.7f} {P_vec_band[1,ib,ik]:13.7f} {P_vec_band[2,ib,ik]:13.7f} {P_vec_band[3,ib,ik]:13.7f}\n"
                f.write(line)
            f.write("\n\n")

    # Write k-paths with nonzero AM and PAM
    with open('/root/capsule/results/kpaths_nonzero_AM_PAM.json', 'w', encoding='utf-8') as f:
        json.dump(data_paths_AM_PAM, f, indent=4, ensure_ascii=False)
    
    if option_nac:
        with open('/root/capsule/results/MagBand.dat', 'w') as fenk:
            first_line = '# M:phonon magnetic moment'
            fenk.write(f"{first_line}\n")
            # Column headers
            line_enk='{:18s}{:16s}{:17s}{:16s}{:16s}{:16s}'.format('#   lenk','enk(meV)','M_x(μₙ)','M_y','M_z','M_abs')
            fenk.write(line_enk+'\n')
            for ib in range(3 * natom):
                for ik in range(nk_total):
                    line = f"{lenk[ik]:10.7f} {enk_band[ib, ik]:15.7f} "
                    line += f"{Mag_band[0, ib, ik]:15.9f} {Mag_band[1, ib, ik]:15.9f} {Mag_band[2, ib, ik]:15.9f} {Mag_abs_band[ib, ik]:15.9f}\n"
                    fenk.write(line)
                fenk.write("\n\n")
        print("Calculation completed. Output files: band.dat, kpoints.dat, Cycloidicity.dat, kpaths_nonzero_AM_PAM.json, MagBand.dat")

    if not option_nac:
            print("Calculation completed. Output files: band.dat, kpoints.dat, Cycloidicity.dat, kpaths_nonzero_AM_PAM.json")

    # =====================
    # VISUALIZATION SECTION
    # =====================
    # Determine energy range for plotting
    emin = min(enk_band[0])-1
    emax = max(enk_band[3*natom-1])+1

    if len(PAM_rotation_type) != 0:
        first_line = '# PAM:'
        pam_info = ''
        for i in range(len(PAM_path_num)):
            pam_path = kpath[PAM_path_num[i]]
            rot_type = PAM_rotation_type[i]
            axis_f = 1.0000001*PAM_rotation_axis[i]
            first_line += ' '+pam_path+'('+str(PAM_path_num[i]+1)+') '+rot_type+ ' ['+str(int(axis_f[0]))+','+str(int(axis_f[1]))+','+str(int(axis_f[2]))+'] ,'
            pam_info += pam_path + ' ' + rot_type + ','
        first_line = first_line.strip(',')
        # print(first_line)

    # ----------------------
    # Modify PAM due to degenerate state
    # ----------------------
        for iii in PAM_path_num:
            # if kpath[iii] in ls_path:
            if True:
                for jj_ls in range(3*natom-1):
                    e_diff=0
                    for ikk in range(nk_p):
                        e_diff+=abs(enk_band[jj_ls+1,iii*nk_p+ikk]-enk_band[jj_ls,iii*nk_p+ikk])
                    e_gap=e_diff/nk_p
                    if abs(e_gap)<1e-3:
                        for ikk in range(nk_p):
                            lph_average = (PAM_band[2,jj_ls,iii*nk_p+ikk] + PAM_band[2,jj_ls+1,iii*nk_p+ikk]) / 2
                            PAM_band[2,jj_ls,iii*nk_p+ikk] = lph_average
                            PAM_band[2,jj_ls+1,iii*nk_p+ikk] = lph_average
                    else:
                        continue

        # Create figure
        fig = plt.figure(figsize=(10, 6), dpi=300)

        # Define layout, "." represents empty space
        layout = [
            ['a', 'b', '.'],
            ['a', 'b', 'c'],
            ['d', 'e', 'c'],
            ['d', 'e', '.']
        ]

        # Create subplots
        ax_dict = fig.subplot_mosaic(layout)

        # subfig (a) for Jx (top-left)
        ax = ax_dict['a']
        for jj in range(3*natom):
            ax.plot(lenk, enk_band[jj], color='black', linewidth='0.4')
            scal = np.array([0.0]*nk_p*nkp_path)
            colors = np.array(['black']*nk_p*nkp_path)
            for ikk in range(nk_p*nkp_path):
                scal[ikk] = abs(chirality_band[jj,ikk,1])*3
                color_ = np.sign(chirality_band[jj,ikk,1])
                if color_ > 0:
                    colors[ikk] = 'red'
                elif color_ < 0:
                    colors[ikk] = 'blue'
                else:
                    colors[ikk] = 'black'
            ax.scatter( np.array(lenk), np.array(enk_band[jj]), s=scal, c=colors, marker='o')
        for ii in range(nkp_path-1):
            ax.plot( [lenk[nk_p*(ii+1)-1], lenk[nk_p*(ii+1)-1]], [emin, emax],color='black')  
        ax.set_xlim(lenk[0], lenk[nk_p * nkp_path - 1])
        ax.set_ylim(emin, emax)
        ind=[0.0]
        for ii in range(nkp_path):
            ind.append(lenk[(ii+1)*nk_p-1])
        ax.set_xticks(ind)
        ax.set_xticklabels(symbol)
        ax.set_ylabel('Frequency (meV)')
        ax.set_title('J$_{x}$')

        # subfig (b) for Jy (top-right)
        ax = ax_dict['b']
        for jj in range(3*natom):
            ax.plot(lenk, enk_band[jj], color='black', linewidth='0.4')
            scal = np.array([0.0]*nk_p*nkp_path)
            colors = np.array(['black']*nk_p*nkp_path)
            for ikk in range(nk_p*nkp_path):
                scal[ikk] = abs(chirality_band[jj,ikk,2])*3
                color_ = np.sign(chirality_band[jj,ikk,2])
                if color_ > 0:
                    colors[ikk] = 'red'
                elif color_ < 0:
                    colors[ikk] = 'blue'
                else:
                    colors[ikk] = 'black'
            ax.scatter( np.array(lenk), np.array(enk_band[jj]), s=scal, c=colors, marker='o')
        for ii in range(nkp_path-1):
            ax.plot( [lenk[nk_p*(ii+1)-1], lenk[nk_p*(ii+1)-1]], [emin, emax],color='black')  
        ax.set_xlim(lenk[0], lenk[nk_p * nkp_path - 1])
        ax.set_ylim(emin, emax)
        ind=[0.0]
        for ii in range(nkp_path):
            ind.append(lenk[(ii+1)*nk_p-1])
        ax.set_xticks(ind)
        ax.set_xticklabels(symbol)
        # ax.ylabel('Frequency (meV)')
        ax.set_title('J$_{y}$')

        # subfig (d) for Jz (bottom-left)
        ax = ax_dict['d']
        for jj in range(3*natom):
            ax.plot(lenk, enk_band[jj], color='black', linewidth='0.4')
            scal = np.array([0.0]*nk_p*nkp_path)
            colors = np.array(['black']*nk_p*nkp_path)
            for ikk in range(nk_p*nkp_path):
                scal[ikk] = abs(chirality_band[jj,ikk,0])*3
                color_ = np.sign(chirality_band[jj,ikk,0])
                if color_ > 0:
                    colors[ikk] = 'red'
                elif color_ < 0:
                    colors[ikk] = 'blue'
                else:
                    colors[ikk] = 'black'
            ax.scatter( np.array(lenk), np.array(enk_band[jj]), s=scal, c=colors,alpha=0.5, marker='o')
        for ii in range(nkp_path-1):
            ax.plot( [lenk[nk_p*(ii+1)-1], lenk[nk_p*(ii+1)-1]], [emin, emax],color='black')  
        ax.set_xlim(lenk[0], lenk[nk_p * nkp_path - 1])
        ax.set_ylim(emin, emax)
        ind=[0.0]
        for ii in range(nkp_path):
            ind.append(lenk[(ii+1)*nk_p-1])
        ax.set_xticks(ind)
        ax.set_xticklabels(symbol)
        ax.set_ylabel('Frequency (meV)')
        ax.set_title('J$_{z}$')

        # subfig (e) for helicity (bottom-left)
        ax = ax_dict['e']
        for jj in range(3*natom):
            ax.plot(lenk, enk_band[jj], color='black', linewidth='0.4')
            scal = np.array([0.0]*nk_p*nkp_path)
            colors = np.array(['black']*nk_p*nkp_path)
            for ikk in range(nk_p*nkp_path):
                scal[ikk] = abs(helicity_band[jj,ikk])*3
                color_ = np.sign(helicity_band[jj,ikk])
                if color_ > 0:
                    colors[ikk] = 'red'
                elif color_ < 0:
                    colors[ikk] = 'blue'
                else:
                    colors[ikk] = 'black'
            ax.scatter( np.array(lenk), np.array(enk_band[jj]), s=scal, c=colors,alpha=0.5, marker='o')
        for ii in range(nkp_path-1):
            ax.plot( [lenk[nk_p*(ii+1)-1], lenk[nk_p*(ii+1)-1]], [emin, emax],color='black')  
        ax.set_xlim(lenk[0], lenk[nk_p * nkp_path - 1])
        ax.set_ylim(emin, emax)
        ind=[0.0]
        for ii in range(nkp_path):
            ind.append(lenk[(ii+1)*nk_p-1])
        ax.set_xticks(ind)
        ax.set_xticklabels(symbol)
        # ax.ylabel('Frequency (meV)')
        ax.set_title('Helicity')

        # subfig (c) for PAM (right-center)
        ax = ax_dict['c']
        for jj in range(3*natom):
            ax.plot( lenk, enk_band[jj],color='black',linewidth='0.4')
        
            sizes = []
            colors = []
            for ikk in range(nk_p*nkp_path):
                pam_val = PAM_band[2,jj,ikk]
                color, size = get_color_and_size(pam_val)
                colors.append(color)
                sizes.append(5 * size)

            ax.scatter(lenk, enk_band[jj], c=colors, s=sizes, marker='o')

        for ii in range(nkp_path-1):
            ax.plot( [lenk[nk_p*(ii+1)-1], lenk[nk_p*(ii+1)-1]], [emin, emax],color='black')
        ax.set_xlim(lenk[0], lenk[nk_p * nkp_path - 1])
        ax.set_ylim(emin, emax)
        ind=[0.0]
        for ii in range(nkp_path):
            ind.append(lenk[(ii+1)*nk_p-1])
        ax.set_xticks(ind)
        ax.set_xticklabels(symbol)
        # ax.ylabel('Frequency (meV)')
        ax.set_title('PAM')

        plt.tight_layout()
        plt.subplots_adjust(wspace=0.2,hspace=0.5)
        plt.savefig('/root/capsule/results/band.png',bbox_inches='tight',dpi=300)

    else:
        first_line = '# No PAM.'

        # Create figure
        fig = plt.figure(figsize=(8, 8), dpi=300)

        # Define layout, "." represents empty space
        layout = [
            ['a', 'b'],
            ['c', 'd'],
        ]

        # Create subplots
        ax_dict = fig.subplot_mosaic(layout)

        # subfig (a) for Jx (top-left)
        ax = ax_dict['a']
        for jj in range(3*natom):
            ax.plot(lenk, enk_band[jj], color='black', linewidth='0.4')
            scal = np.array([0.0]*nk_p*nkp_path)
            colors = np.array(['black']*nk_p*nkp_path)
            for ikk in range(nk_p*nkp_path):
                scal[ikk] = abs(chirality_band[jj,ikk,1])*3
                color_ = np.sign(chirality_band[jj,ikk,1])
                if color_ > 0:
                    colors[ikk] = 'red'
                elif color_ < 0:
                    colors[ikk] = 'blue'
                else:
                    colors[ikk] = 'black'
            ax.scatter( np.array(lenk), np.array(enk_band[jj]), s=scal, c=colors, marker='o')
        for ii in range(nkp_path-1):
            ax.plot( [lenk[nk_p*(ii+1)-1], lenk[nk_p*(ii+1)-1]], [emin, emax],color='black')  
        ax.set_xlim(lenk[0], lenk[nk_p * nkp_path - 1])
        ax.set_ylim(emin, emax)
        ind=[0.0]
        for ii in range(nkp_path):
            ind.append(lenk[(ii+1)*nk_p-1])
        ax.set_xticks(ind)
        ax.set_xticklabels(symbol)
        ax.set_ylabel('Frequency (meV)')
        ax.set_title('J$_{x}$')

        # subfig (b) for Jy (top-right)
        ax = ax_dict['b']
        for jj in range(3*natom):
            ax.plot(lenk, enk_band[jj], color='black', linewidth='0.4')
            scal = np.array([0.0]*nk_p*nkp_path)
            colors = np.array(['black']*nk_p*nkp_path)
            for ikk in range(nk_p*nkp_path):
                scal[ikk] = abs(chirality_band[jj,ikk,2])*3
                color_ = np.sign(chirality_band[jj,ikk,2])
                if color_ > 0:
                    colors[ikk] = 'red'
                elif color_ < 0:
                    colors[ikk] = 'blue'
                else:
                    colors[ikk] = 'black'
            ax.scatter( np.array(lenk), np.array(enk_band[jj]), s=scal, c=colors, marker='o')
        for ii in range(nkp_path-1):
            ax.plot( [lenk[nk_p*(ii+1)-1], lenk[nk_p*(ii+1)-1]], [emin, emax],color='black')  
        ax.set_xlim(lenk[0], lenk[nk_p * nkp_path - 1])
        ax.set_ylim(emin, emax)
        ind=[0.0]
        for ii in range(nkp_path):
            ind.append(lenk[(ii+1)*nk_p-1])
        ax.set_xticks(ind)
        ax.set_xticklabels(symbol)
        # ax.ylabel('Frequency (meV)')
        ax.set_title('J$_{y}$')

        # subfig (c) for Jz (bottom-left)
        ax = ax_dict['c']
        for jj in range(3*natom):
            ax.plot(lenk, enk_band[jj], color='black', linewidth='0.4')
            scal = np.array([0.0]*nk_p*nkp_path)
            colors = np.array(['black']*nk_p*nkp_path)
            for ikk in range(nk_p*nkp_path):
                scal[ikk] = abs(chirality_band[jj,ikk,0])*3
                color_ = np.sign(chirality_band[jj,ikk,0])
                if color_ > 0:
                    colors[ikk] = 'red'
                elif color_ < 0:
                    colors[ikk] = 'blue'
                else:
                    colors[ikk] = 'black'
            ax.scatter( np.array(lenk), np.array(enk_band[jj]), s=scal, c=colors,alpha=0.5, marker='o')
        for ii in range(nkp_path-1):
            ax.plot( [lenk[nk_p*(ii+1)-1], lenk[nk_p*(ii+1)-1]], [emin, emax],color='black')  
        ax.set_xlim(lenk[0], lenk[nk_p * nkp_path - 1])
        ax.set_ylim(emin, emax)
        ind=[0.0]
        for ii in range(nkp_path):
            ind.append(lenk[(ii+1)*nk_p-1])
        ax.set_xticks(ind)
        ax.set_xticklabels(symbol)
        ax.set_ylabel('Frequency (meV)')
        ax.set_title('J$_{z}$')

        # subfig (d) for helicity (bottom-left)
        ax = ax_dict['d']
        for jj in range(3*natom):
            ax.plot(lenk, enk_band[jj], color='black', linewidth='0.4')
            scal = np.array([0.0]*nk_p*nkp_path)
            colors = np.array(['black']*nk_p*nkp_path)
            for ikk in range(nk_p*nkp_path):
                scal[ikk] = abs(helicity_band[jj,ikk])*3
                color_ = np.sign(helicity_band[jj,ikk])
                if color_ > 0:
                    colors[ikk] = 'red'
                elif color_ < 0:
                    colors[ikk] = 'blue'
                else:
                    colors[ikk] = 'black'
            ax.scatter( np.array(lenk), np.array(enk_band[jj]), s=scal, c=colors,alpha=0.5, marker='o')
        for ii in range(nkp_path-1):
            ax.plot( [lenk[nk_p*(ii+1)-1], lenk[nk_p*(ii+1)-1]], [emin, emax],color='black')  
        ax.set_xlim(lenk[0], lenk[nk_p * nkp_path - 1])
        ax.set_ylim(emin, emax)
        ind=[0.0]
        for ii in range(nkp_path):
            ind.append(lenk[(ii+1)*nk_p-1])
        ax.set_xticks(ind)
        ax.set_xticklabels(symbol)
        # ax.ylabel('Frequency (meV)')
        ax.set_title('Helicity')

        plt.tight_layout()
        plt.subplots_adjust(wspace=0.2,hspace=0.2)
        plt.savefig('/root/capsule/results/band.png',bbox_inches='tight',dpi=300)

    #For phonon magnetic moment
    if option_nac:
        fig = plt.figure(figsize=(14, 12), dpi=300)
        layout = [
            ['a', 'b'],
            ['a', 'b'],
            ['c', 'd'],
            ['c', 'd'],
        ]

        ax_dict = fig.subplot_mosaic(layout)
        all_M_data = np.concatenate([Mag_band.flatten()])
        max_abs = max(abs(np.min(all_M_data)), abs(np.max(all_M_data)))
        vmin, vmax = -max_abs, max_abs

        absmax = Mag_abs_band.max()
        normabs = Normalize(0, absmax)
        cmapabs = plt.cm.Oranges

        colors = ['blue','white','red']
        cmap = LinearSegmentedColormap.from_list('enhanced_colors', colors, N=256)
        linthresh = max_abs * 0.1
        norm = SymLogNorm(linthresh=linthresh, linscale=1, vmin=-max_abs, vmax=max_abs)
        fixed_size = 6

        plot_component(ax_dict['a'], Mag_band[0,:,:], 'M$_{x}$', 1)
        plot_component(ax_dict['b'], Mag_band[1,:,:], 'M$_{y}$', 0)
        plot_component(ax_dict['c'], Mag_band[2,:,:], 'M$_{z}$', 1)

        ax = ax_dict['d']
        for jj in range(3 * natom):
            ax.plot(lenk, enk_band[jj], color='black', linewidth=0.4)
            point_colors = [cmapabs(normabs(Mag_abs_band[jj, ikk])) for ikk in range(nk_total)]   
            ax.scatter(lenk, enk_band[jj], s=fixed_size, c=point_colors, 
                    marker='o', edgecolors='none')
        for ii in range(nkp_path-1):
            ax.plot( [lenk[nk_p*(ii+1)-1], lenk[nk_p*(ii+1)-1]], [emin, emax], color='black')  
        ax.set_xlim(lenk[0], lenk[nk_p * nkp_path - 1])
        ax.set_ylim(emin, emax)
        ind=[0.0]
        for ii in range(nkp_path):
            ind.append(lenk[(ii+1)*nk_p-1])
        ax.set_xticks(ind)
        ax.set_xticklabels(symbol)
        ax.set_title('M$_{abs}$')

        cbar_ax = fig.add_axes([0.92, 0.52, 0.015, 0.35])
        cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), 
                        cax=cbar_ax, 
                        orientation='vertical')
        cbar.set_ticks([vmin, 0, vmax])
        cbar.set_ticklabels([f'{vmin:.3f}', '0', f'{vmax:.3f}'])
        cbar.ax.yaxis.set_label_position('left') 
        cbar.ax.tick_params(labelsize=9)

        cbar_ax = fig.add_axes([0.92, 0.12, 0.015, 0.35])
        cbar = plt.colorbar(plt.cm.ScalarMappable(norm=normabs, cmap=cmapabs), 
                        cax=cbar_ax, 
                        orientation='vertical')
        cbar.set_ticks([0, absmax])
        cbar.set_ticklabels([f'0', f'{absmax:.3f}'])
        cbar.ax.yaxis.set_label_position('left') 
        cbar.ax.tick_params(labelsize=9)

        plt.savefig('/root/capsule/results/MAGBAND.png', bbox_inches='tight', dpi=300)
        plt.close()

    #For Cycloidicity
    fig = plt.figure(figsize=(14, 12), dpi=300)
    layout = [
        ['a', 'b'],
        ['a', 'b'],
        ['c', 'd'],
        ['c', 'd'],
    ]

    ax_dict = fig.subplot_mosaic(layout)
    all_P_data = np.concatenate([P_vec_band[:3,:,:].flatten()])
    max_abs = max(abs(np.min(all_P_data)), abs(np.max(all_P_data)))
    vmin, vmax = -max_abs, max_abs

    absmax = P_vec_band[3,:,:].max()
    normabs = Normalize(0, absmax)
    cmapabs = plt.cm.Oranges

    colors = ['blue','white','red']
    cmap = LinearSegmentedColormap.from_list('enhanced_colors', colors, N=256)
    linthresh = max_abs * 0.1
    norm = SymLogNorm(linthresh=linthresh, linscale=1, vmin=-max_abs, vmax=max_abs)
    fixed_size = 6

    plot_component(ax_dict['a'], P_vec_band[0,:,:], 'P$_{x}$', 1)
    plot_component(ax_dict['b'], P_vec_band[1,:,:], 'P$_{y}$', 0)
    plot_component(ax_dict['c'], P_vec_band[2,:,:], 'P$_{z}$', 1)

    ax = ax_dict['d']
    for jj in range(3 * natom):
        ax.plot(lenk, enk_band[jj], color='black', linewidth=0.4)
        point_colors = [cmapabs(normabs(P_vec_band[3, jj, ikk])) for ikk in range(nk_total)]   
        ax.scatter(lenk, enk_band[jj], s=fixed_size, c=point_colors, 
                marker='o', edgecolors='none')
    for ii in range(nkp_path-1):
        ax.plot( [lenk[nk_p*(ii+1)-1], lenk[nk_p*(ii+1)-1]], [emin, emax], color='black')  
    ax.set_xlim(lenk[0], lenk[nk_p * nkp_path - 1])
    ax.set_ylim(emin, emax)
    ind=[0.0]
    for ii in range(nkp_path):
        ind.append(lenk[(ii+1)*nk_p-1])
    ax.set_xticks(ind)
    ax.set_xticklabels(symbol)
    ax.set_title('P$_{abs}$')

    cbar_ax = fig.add_axes([0.92, 0.52, 0.015, 0.35])
    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), 
                    cax=cbar_ax, 
                    orientation='vertical')
    cbar.set_ticks([vmin, 0, vmax])
    cbar.set_ticklabels([f'{vmin:.3f}', '0', f'{vmax:.3f}'])
    cbar.ax.yaxis.set_label_position('left') 
    cbar.ax.tick_params(labelsize=9)

    cbar_ax = fig.add_axes([0.92, 0.12, 0.015, 0.35])
    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=normabs, cmap=cmapabs), 
                    cax=cbar_ax, 
                    orientation='vertical')
    cbar.set_ticks([0, absmax])
    cbar.set_ticklabels([f'0', f'{absmax:.3f}'])
    cbar.ax.yaxis.set_label_position('left') 
    cbar.ax.tick_params(labelsize=9)

    plt.savefig('/root/capsule/results/Cycloidicity.png', bbox_inches='tight', dpi=300)
    plt.close()
