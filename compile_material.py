import sys
import json
import numpy as np

# Add phononweb to path
sys.path.insert(0, 'python/')
from phononweb.phonon import Phonon
from phononweb.runtime_dynamical_matrix import build_dynamical_matrix_from_payload

def compile_material(json_path, output_path=None, on_the_fly=False):
    if output_path is None:
        output_path = json_path
        
    print(f"Reading {json_path}...")
    with open(json_path, 'r') as f:
        data = json.load(f)
        
    qpoints = np.array(data['qpoints'], dtype=float)
    nqpoints = len(qpoints)
    natoms = data['natoms']
    
    # We will determine nphons dynamically from the data, but default to natoms * 3
    nphons = natoms * 3

    if 'dynamical_matrix' in data:
        print("Evaluating eigenvectors statically from dynamical matrix...")
        eigs_val = np.zeros((nqpoints, nphons), dtype=float)
        eigs_real_imag = np.zeros((nqpoints, nphons, natoms, 3, 2), dtype=float)
        
        def get_nac_q_direction(q_idx):
            q = qpoints[q_idx]
            if np.linalg.norm(q) > 1e-8:
                return None
            line_breaks = data.get('line_breaks', [])
            seg_start, seg_end = 0, len(qpoints) - 1
            for b_start, b_end in line_breaks:
                if b_start <= q_idx < b_end:
                    seg_start, seg_end = b_start, b_end - 1
                    break
            d = qpoints[seg_start] - qpoints[seg_end]
            if np.linalg.norm(d) > 1e-8: return d.tolist()
            if q_idx > seg_start:
                d = qpoints[q_idx-1] - q
                if np.linalg.norm(d) > 1e-8: return d.tolist()
            if q_idx < seg_end:
                d = q - qpoints[q_idx+1]
                if np.linalg.norm(d) > 1e-8: return d.tolist()
            return None

        # We first calculate pure energy-sorted eigenvectors
        for i, q in enumerate(qpoints):
            q_dir = get_nac_q_direction(i)
            dm = build_dynamical_matrix_from_payload(data['dynamical_matrix'], q, q_direction=q_dir)
            w, v = np.linalg.eigh(dm)
            
            # Sort by eigenvalue
            idx = np.argsort(w)
            w = w[idx]
            v = v[:, idx] # eigs as columns
            
            # Phonopy's default VASP conversion factor to THz is 15.633302
            # 1 THz = 33.35641 cm^-1
            freqs = np.sqrt(np.abs(w)) * np.sign(w) * 15.633302 * 33.35641
            eigs_val[i, :] = freqs
            
            # v[:, j] is the j-th eigenvector. We must shape it to (natoms, 3)
            for ibnd in range(nphons):
                e = v[:, ibnd]
                e_comp = e.reshape(natoms, 3)
                eigs_real_imag[i, ibnd, ..., 0] = np.real(e_comp)
                eigs_real_imag[i, ibnd, ..., 1] = np.imag(e_comp)
    else:
        print("Loading pre-computed eigenvectors from JSON...")
        eigs_val = np.array(data['eigenvalues'], dtype=float)
        eigs_real_imag = np.array(data['vectors'], dtype=float)
        nphons = eigs_val.shape[1]

    print("Building Phonon compilation object...")
    p = Phonon()
    p.name = data.get('name', 'Unknown')
    p.nqpoints = nqpoints
    p.natoms = natoms
    p.nphons = nphons
    p.qpoints = qpoints.tolist()
    p.cell = data.get('lattice')
    p.pos = data.get('atom_pos_red')
    p.atom_numbers = data.get('atom_numbers')
    p.atom_types = data.get('atom_types', [])
    p.chemical_symbols = p.atom_types # fallback
    p.atomic_numbers = list(set(p.atom_numbers)) # simplified
    p.chemical_formula = data.get('formula', '')
    p.reps = data.get('repetitions', [3, 3, 3])
    p.distances = data.get('distances', [])
    p.highsym_qpts = data.get('highsym_qpts', [])
    p.line_breaks = data.get('line_breaks', [])
    
    p.eigenvalues = eigs_val
    p.eigenvectors = eigs_real_imag
    
    print("Reordering bands continuously...")
    p.reorder_eigenvalues()
    
    if on_the_fly:
        print("Skipping PAM and Chirality computation (--on-the-fly enabled).")
        # Ensure we delete the heavy arrays if they exist in the payload
        for key in ['eigenvalues', 'vectors', 'pam_total_compensated', 'pam_total_uncompensated', 'angular_momentum_x', 'angular_momentum_y', 'angular_momentum_z', 'cycloidicity_x', 'cycloidicity_y', 'cycloidicity_z', 'helicity', 'magnetic_moment_x', 'magnetic_moment_y', 'magnetic_moment_z', 'pam_total', 'pam_rotation_only']:
            data.pop(key, None)
    else:
        print("Evaluating physical properties (PAM, Symmetry)...")
        from phononweb.pam import compute_pam_properties
        pam_props = compute_pam_properties(p)
        data.update(pam_props)
        
        from phononweb.chirality import compute_chiral_properties
        chiral_props = compute_chiral_properties(p)
        data.update(chiral_props)
    
    from phononweb.sym_using_spglib import get_crystal_symmetry, get_symmetry_labels_for_path
    try:
        sym_cell = (np.array(p.cell), np.array(p.pos), np.array(p.atom_numbers))
        data["crystal_symmetry"] = get_crystal_symmetry(sym_cell)
        sym_labels = get_symmetry_labels_for_path(sym_cell, p.qpoints, p.highsym_qpts)
        data["highsym_point_groups"] = sym_labels["highsym_point_groups"]
        data["segment_point_groups"] = sym_labels["segment_point_groups"]
    except Exception as e:
        print(f"Warning: could not compute symmetry labels: {e}")
        
    if not on_the_fly:
        data['eigenvalues'] = p.eigenvalues.tolist()
        data['vectors'] = p.eigenvectors.tolist()
    
    print(f"Exporting compiled material to {output_path}...")
    with open(output_path, 'w') as f:
        json.dump(data, f)
    print("Done!")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 compile_material.py <path_to_data.json> [--on-the-fly]")
        sys.exit(1)
        
    on_the_fly = "--on-the-fly" in sys.argv
    path = sys.argv[1] if sys.argv[1] != "--on-the-fly" else sys.argv[2]
    compile_material(path, on_the_fly=on_the_fly)

