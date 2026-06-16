import numpy as np
from phononweb.sym_using_spglib import get_little_group_ops

def get_rotation_axis_angle(R_car):
    """
    Extract the rotation axis and angle from a 3x3 Cartesian rotation matrix.
    
    Parameters:
        R_car: 3x3 array-like, rotation matrix in Cartesian coordinates.
        
    Returns:
        axis: 3D unit vector representing rotation axis.
        angle: float, rotation angle in radians.
    """
    R_car = np.array(R_car, dtype=float)
    tr = np.trace(R_car)
    cos_val = (tr - 1.0) / 2.0
    cos_val = np.clip(cos_val, -1.0, 1.0)
    angle = np.arccos(cos_val)
    
    if np.isclose(angle, 0.0, atol=1e-5):
        # Identity matrix: rotation axis is arbitrary, angle is 0
        return np.array([0.0, 0.0, 1.0]), 0.0
        
    if np.isclose(angle, np.pi, atol=1e-5):
        # 180 degrees rotation: R_car is symmetric.
        # Find the eigenvector corresponding to eigenvalue 1.0
        w, v = np.linalg.eig(R_car)
        idx = np.argmin(np.abs(w - 1.0))
        axis = np.real(v[:, idx])
        axis = axis / np.linalg.norm(axis)
        return axis, np.pi
        
    # Standard rotation: axis is antisymmetric component
    axis = np.array([
        R_car[2, 1] - R_car[1, 2],
        R_car[0, 2] - R_car[2, 0],
        R_car[1, 0] - R_car[0, 1]
    ], dtype=float)
    
    norm = np.linalg.norm(axis)
    if norm > 1e-5:
        axis = axis / norm
    else:
        # Fallback to eigenvalue method
        w, v = np.linalg.eig(R_car)
        idx = np.argmin(np.abs(w - 1.0))
        axis = np.real(v[:, idx])
        axis = axis / np.linalg.norm(axis)
        
    return axis, angle

def evaluate_pam_expectations(cell, q, e, lg_ops=None, symprec=1e-5):
    """
    Evaluate the Spin, Orbital, and Total Pseudo-Angular Momentum expectation values
    and phases for a phonon mode eigenvector under all little group symmetry operations.
    
    Parameters:
        cell: tuple (lattice, positions, numbers)
        q: array-like of shape (3,), wavevector in fractional coordinates
        e: array-like of shape (3*Natoms,) or (Natoms, 3), complex eigenvector.
           Will be normalized internally.
        lg_ops: dict, pre-computed little group operations. If None, evaluated internally.
        symprec: float, symmetry precision
        
    Returns:
        list of dicts, one for each little group operation containing expectation values.
    """
    lattice = np.array(cell[0], dtype=float)
    pos_red = np.array(cell[1], dtype=float)
    natoms = len(pos_red)
    
    q = np.array(q, dtype=float)
    e = np.array(e, dtype=complex).reshape(natoms, 3)
    
    # Normalize eigenvector
    norm = np.linalg.norm(e)
    if norm > 1e-12:
        e = e / norm
        
    # Get direct lattice transformation matrix (columns are lattice vectors)
    # x_car = A * x_red
    A = lattice.T
    A_inv = np.linalg.inv(A)
    
    # Retrieve little group operations
    if lg_ops is None:
        lg_ops = get_little_group_ops(cell, q, symprec=symprec)
    results = []
    
    for idx, R in enumerate(lg_ops["rotations"]):
        tau = lg_ops["translations"][idx]
        original_idx = lg_ops["indices"][idx]
        
        # Transform rotation matrix to Cartesian coordinates
        R_car = np.dot(A, np.dot(R, A_inv))
        
        # Get rotation axis and angle in Cartesian
        axis, angle = get_rotation_axis_angle(R_car)
        
        # Initialize the transformed eigenvectors
        De = np.zeros_like(e, dtype=complex) # Total
        Se = np.zeros_like(e, dtype=complex) # Spin only
        Le = np.zeros_like(e, dtype=complex) # Orbital only
        
        for i in range(natoms):
            # Rotated fractional position
            r_prime = np.dot(R, pos_red[i]) + tau
            
            # Find matching atom index j in the primitive cell (modulo 1)
            diffs = pos_red - r_prime
            diffs_mod1 = diffs - np.round(diffs)
            dists = np.linalg.norm(diffs_mod1, axis=1)
            j = np.argmin(dists)
            
            if dists[j] > 1e-4:
                raise ValueError(f"Symmetry mapping failed: no matching atom for rotated pos of atom {i}")
                
            # Integer cell shift vector R_L
            R_L = np.round(r_prime - pos_red[j])
            
            # Bloch phase factor from integer translation R_L
            phase_RL = np.exp(-1j * 2.0 * np.pi * np.dot(q, R_L))
            
            # Apply operators
            De[j] = np.dot(R_car, e[i]) * phase_RL
            Le[j] = e[i] * phase_RL
            Se[i] = np.dot(R_car, e[i])
            
        # Compute expectation values: <e | Op | e>
        exp_D = np.vdot(e, De) 
        exp_S = np.vdot(e, Se)
        exp_L = np.vdot(e, Le)
        
        # Phases
        phase_tot = np.angle(exp_D)
        phase_translation = -2.0 * np.pi * np.dot(q, tau)
        
        # Subtract fractional translation phase shift to get rotation-only phase
        exp_D_rot_only = exp_D * np.exp(-1j * phase_translation)
        phase_rotation_only = np.angle(exp_D_rot_only)
        
        phase_spin = np.angle(exp_S)
        phase_orbital = np.angle(exp_L)
        
        results.append({
            "index": int(original_idx),
            "rotation_axis": axis.tolist(),
            "rotation_angle": float(angle),
            "det": int(np.round(np.linalg.det(R))),
            "translation": tau.tolist(),
            "exp_D": complex(exp_D),
            "exp_S": complex(exp_S),
            "exp_L": complex(exp_L),
            "phase_tot": float(phase_tot),
            "phase_translation": float(phase_translation),
            "phase_rotation_only": float(phase_rotation_only),
            "phase_spin": float(phase_spin),
            "phase_orbital": float(phase_orbital)
        })
        
    return results

def compute_pam_properties(phonon_obj):
    """
    Compute Phonon Pseudo-Angular Momentum (PAM) values for all k-points and all bands.
    Returns PAM values both with and without the fractional translation phase.
    """
    if not hasattr(phonon_obj, 'eigenvectors') or phonon_obj.eigenvectors is None:
        return {}
        
    eigs_real_imag = phonon_obj.eigenvectors
    nqpoints, nphons, natoms, _, _ = eigs_real_imag.shape
    
    eigs = eigs_real_imag[..., 0] + 1j * eigs_real_imag[..., 1]
    eigs = eigs.reshape(nqpoints, nphons, natoms * 3)
    
    qpoints = np.array(phonon_obj.qpoints, dtype=float)
    cell = (phonon_obj.cell, phonon_obj.pos, phonon_obj.atom_numbers)
    
    pam_tot = np.zeros((nqpoints, nphons))
    pam_rot = np.zeros((nqpoints, nphons))
    
    line_breaks = getattr(phonon_obj, 'line_breaks', None)
    if line_breaks and len(line_breaks) > 0:
        segments = line_breaks
    else:
        segments = [[0, nqpoints]]
        
    A = np.array(cell[0], dtype=float).T
    A_inv = np.linalg.inv(A)
    
    from phononweb.lattice import rec_lat, red_car
    B_rec = rec_lat(cell[0])

    # use local get_rotation_axis_angle defined above

    for segment in segments:
        start_idx, end_idx = segment[0], segment[1]
        
        # We determine the symmetry for the entire path segment using its midpoint
        mid_idx = (start_idx + end_idx - 1) // 2
        q_mid = qpoints[mid_idx]
        lg_ops_segment = get_little_group_ops(cell, q_mid)
        
        # Determine wave propagation direction in Cartesian coordinates
        q_start = qpoints[start_idx]
        q_end_actual = qpoints[min(end_idx, nqpoints - 1)]
        delta_q_red = q_end_actual - q_start
        delta_q_car = red_car([delta_q_red], B_rec)[0]
        
        norm_dq = np.linalg.norm(delta_q_car)
        if norm_dq > 1e-5:
            dir_q_car = delta_q_car / norm_dq
        else:
            dir_q_car = np.array([0.0, 0.0, 1.0]) # Fallback for a 0-length segment
            
        # Find the primary generator Cn operation ONCE for the segment
        best_op_index = -1
        min_angle = 10.0 # larger than 2*pi
        
        for idx_op, R in enumerate(lg_ops_segment["rotations"]):
            det = int(np.round(np.linalg.det(R)))
            if det > 0: # proper rotations only
                R_car = np.dot(A, np.dot(R, A_inv))
                axis, angle = get_rotation_axis_angle(R_car)
                
                if angle > 1e-3:
                    # Check if rotation axis is parallel to wave propagation direction
                    dot_val = np.abs(np.dot(axis, dir_q_car))
                    if np.isclose(dot_val, 1.0, atol=1e-2):
                        # We want the generator, which has the smallest fundamental angle
                        if angle < min_angle:
                            min_angle = angle
                            best_op_index = idx_op
                            
        if best_op_index == -1:
            filtered_lg_ops = None
        else:
            filtered_lg_ops = {
                "rotations": [lg_ops_segment["rotations"][best_op_index]],
                "translations": [lg_ops_segment["translations"][best_op_index]],
                "rec_rotations": [lg_ops_segment["rec_rotations"][best_op_index]],
                "indices": [lg_ops_segment["indices"][best_op_index]]
            }
        
        for iq in range(start_idx, end_idx):
            if iq >= nqpoints:
                continue
            q = qpoints[iq]
            
            for ibnd in range(nphons):
                if filtered_lg_ops is None:
                    pam_tot[iq, ibnd] = 0.0
                    pam_rot[iq, ibnd] = 0.0
                    continue
                    
                e = eigs[iq, ibnd, :]
                
                # Evaluate expectation ONLY for the primary Cn operation
                results = evaluate_pam_expectations(cell, q, e, lg_ops=filtered_lg_ops)
                best_res = results[0]
                
                angle = best_res['rotation_angle']
                m_tot = best_res['phase_tot'] / angle
                m_rot = best_res['phase_rotation_only'] / angle
                
                if abs(m_tot - round(m_tot)) < 1e-2:
                    m_tot = round(m_tot)
                if abs(m_rot - round(m_rot)) < 1e-2:
                    m_rot = round(m_rot)
                    
                pam_tot[iq, ibnd] = m_tot
                pam_rot[iq, ibnd] = m_rot

    return {
        'pam_total': pam_tot.tolist(),
        'pam_rotation_only': pam_rot.tolist()
    }
