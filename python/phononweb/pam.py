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

def build_representation_matrices(cell, q, R, tau):
    lattice = np.array(cell[0], dtype=float)
    pos_red = np.array(cell[1], dtype=float)
    natoms = len(pos_red)
    q = np.array(q, dtype=float)
    
    A = lattice.T
    A_inv = np.linalg.inv(A)
    R_car = np.dot(A, np.dot(R, A_inv))
    
    D = np.zeros((3*natoms, 3*natoms), dtype=complex)
    
    for i in range(natoms):
        r_prime = np.dot(R, pos_red[i]) + tau
        diffs = pos_red - r_prime
        diffs_mod1 = diffs - np.round(diffs)
        diffs_car = np.dot(diffs_mod1, A.T)
        dists = np.linalg.norm(diffs_car, axis=1)
        j = np.argmin(dists)
        if dists[j] > 1e-3:
            raise ValueError(f"Symmetry mapping failed for atom {i}")
        R_L = np.round(r_prime - pos_red[j])
        
        # Phonopy Type II convention:
        # D_ji = R_car * exp(-i 2pi G . r_j)
        G = np.round(np.dot(R, q) - q)
        phase_RL = np.exp(-2.0j * np.pi * np.dot(G, pos_red[j]))
            
        D[3*j:3*j+3, 3*i:3*i+3] = phase_RL * R_car
    return D

def compute_pam_properties(phonon_obj):
    if not hasattr(phonon_obj, 'eigenvectors') or phonon_obj.eigenvectors is None:
        return {}
    
    #initialization of arrays:

    eigs_real_imag = phonon_obj.eigenvectors
    nqpoints, nphons, natoms, _, _ = eigs_real_imag.shape
    
    eigs = eigs_real_imag[..., 0] + 1j * eigs_real_imag[..., 1]
    eigs = eigs.reshape(nqpoints, nphons, natoms * 3)
    
    freqs = np.array(phonon_obj.eigenvalues)
    qpoints = np.array(phonon_obj.qpoints, dtype=float)
    cell = (phonon_obj.cell, phonon_obj.pos, phonon_obj.atom_numbers)
    
    pam_tot_uncomp = np.zeros((nqpoints, nphons))
    pam_tot_comp = np.zeros((nqpoints, nphons))
    
    # Store symmetry-adapted eigenvectors back to the array
    eigs_adapted = np.zeros_like(eigs)

    line_breaks = getattr(phonon_obj, 'line_breaks', None)
    if line_breaks and len(line_breaks) > 0:
        segments = line_breaks
    else:
        segments = [[0, nqpoints]]
        
    A = np.array(cell[0], dtype=float).T
    A_inv = np.linalg.inv(A)
    
    from phononweb.lattice import rec_lat, red_car
    B_rec = rec_lat(cell[0])

    for segment in segments:
        start_idx, end_idx = segment[0], segment[1]
        
        mid_idx = (start_idx + end_idx - 1) // 2
        q_mid = qpoints[mid_idx]
        lg_ops_segment = get_little_group_ops(cell, q_mid)
        
        q_start = qpoints[start_idx]
        q_end_actual = qpoints[min(end_idx, nqpoints - 1)]
        delta_q_red = q_end_actual - q_start
        delta_q_car = red_car([delta_q_red], B_rec)[0]
        
        norm_dq = np.linalg.norm(delta_q_car)
        if norm_dq > 1e-5:
            dir_q_car = delta_q_car / norm_dq
        else:
            dir_q_car = np.array([0.0, 0.0, 1.0])
            
        best_op_index = -1
        min_angle = 10.0
        
        for idx_op, R in enumerate(lg_ops_segment["rotations"]):
            det = int(np.round(np.linalg.det(R)))
            if det > 0:
                R_car = np.dot(A, np.dot(R, A_inv))
                axis, angle = get_rotation_axis_angle(R_car)
                if angle > 1e-3:
                    dot_val = np.abs(np.dot(axis, dir_q_car))
                    if np.isclose(dot_val, 1.0, atol=1e-2):
                        if angle < min_angle:
                            min_angle = angle
                            best_op_index = idx_op
                            
        if best_op_index == -1:
            # Fallback to no rotation processing for this segment
            # Just copy eigs over un-adapted.
            for iq in range(start_idx, end_idx):
                if iq < nqpoints:
                    eigs_adapted[iq] = eigs[iq]
            continue
            
        R_gen = lg_ops_segment["rotations"][best_op_index]
        tau_gen = lg_ops_segment["translations"][best_op_index]
        
        R_car_gen = np.dot(A, np.dot(R_gen, A_inv))
        _, angle_gen = get_rotation_axis_angle(R_car_gen)
        n_axis = int(np.round(2.0 * np.pi / angle_gen)) if angle_gen > 1e-4 else 1

        # Calculate projective representation compensation translation vector t
        # t = 1/n sum_{i=0}^{n-1} (R_gen)^i tau_gen
        t_vec = np.zeros(3)
        R_pow = np.eye(3)
        for _ in range(n_axis):
            t_vec += np.dot(R_pow, tau_gen)
            R_pow = np.dot(R_gen, R_pow)
        t_vec /= n_axis

        for iq in range(start_idx, end_idx):
            if iq >= nqpoints:
                continue
                
            # --- GLOBAL TRACKING FIX ---
            # Reorder eigs[iq] to match eigs_adapted[iq-1] across the whole segment!
            if iq > start_idx:
                from scipy.optimize import linear_sum_assignment
                overlap_mat = np.abs(np.dot(eigs_adapted[iq-1].conj(), eigs[iq].T))
                row_ind, col_ind = linear_sum_assignment(-overlap_mat)
                eigs[iq] = eigs[iq, col_ind]
                freqs[iq] = freqs[iq, col_ind]
            # ---------------------------

            q = qpoints[iq]
            w = freqs[iq]
            
            # Construct full representation matrices for this q and operation
            try:
                D = build_representation_matrices(cell, q, R_gen, tau_gen)
            except ValueError:
                # If symmetry mapping fails, just fallback
                eigs_adapted[iq] = eigs[iq]
                continue
            
            # Group bands into degenerate manifolds
            # Use a tolerance of 0.1 cm⁻¹ to reliably catch DFT degeneracies
            tol = 1e-1
            visited = np.zeros(nphons, dtype=bool)
            
            for ibnd in range(nphons):
                if visited[ibnd]:
                    continue
                    
                # Find all bands degenerate with ibnd
                manifold = [j for j in range(nphons) if abs(w[j] - w[ibnd]) < tol]
                for j in manifold:
                    visited[j] = True
                    
                m_size = len(manifold)
                basis = eigs[iq, manifold, :] # shape (m_size, 3N)
                
                if m_size == 1:
                    # Non-degenerate
                    e_vec = basis[0]
                    # <e | D | e>
                    lambda_D = np.vdot(e_vec, np.dot(D, e_vec))
                    
                    # 1. Uncompensated (drifting)
                    phase_uncomp = np.angle(lambda_D)
                    pam_tot_uncomp[iq, ibnd] = (phase_uncomp / angle_gen) % n_axis
                    
                    # 2. Compensated (quantized for non-symmorphic)
                    # D_comp = D * exp(i 2pi q . t)
                    phase_comp_factor = np.exp(2.0j * np.pi * np.dot(q, t_vec))
                    lambda_D_comp = lambda_D * phase_comp_factor
                    phase_comp = np.angle(lambda_D_comp)
                    pam_comp = (phase_comp / angle_gen) % n_axis
                    pam_comp_err = pam_comp - np.round(pam_comp)
                    if np.abs(pam_comp_err) < 1e-2:
                        pam_comp = np.round(pam_comp) % n_axis
                    pam_tot_comp[iq, ibnd] = pam_comp
                    
                    eigs_adapted[iq, ibnd] = e_vec
                else:
                    # Degenerate Subspace: Project D matrix into the basis
                    M_D = np.zeros((m_size, m_size), dtype=complex)
                    for a in range(m_size):
                        for b in range(m_size):
                            M_D[a, b] = np.vdot(basis[a], np.dot(D, basis[b]))
                            
                    # Diagonalize M_D to get symmetry-adapted (circularly polarized) eigenstates.
                    # M_D is the projection of a unitary operator, so its eigenvectors are the
                    # correct basis that diagonalizes the symmetry operation.
                    eigvals_D, eigvecs_D = np.linalg.eig(M_D)
                    
                    # Compute adapted modes
                    adapted_modes = []
                    for m_idx in range(m_size):
                        V = eigvecs_D[:, m_idx]
                        mode = np.zeros(3*natoms, dtype=complex)
                        for a in range(m_size):
                            mode += V[a] * basis[a]
                        norm = np.linalg.norm(mode)
                        if norm > 1e-12:
                            mode /= norm
                        adapted_modes.append(mode)

                    # Match with previous q-point to preserve band tracking
                    ref_modes = eigs_adapted[iq - 1, manifold] if iq > start_idx else basis
                    
                    # Compute overlap matrix: overlap[a, b] = |<adapted_modes[a] | ref_modes[b]>|
                    overlap = np.zeros((m_size, m_size))
                    for a in range(m_size):
                        for b in range(m_size):
                            overlap[a, b] = np.abs(np.vdot(adapted_modes[a], ref_modes[b]))
                            
                    # Strict maximum weight bipartite matching using the Hungarian algorithm
                    from scipy.optimize import linear_sum_assignment
                    # linear_sum_assignment minimizes cost, so we supply -overlap to maximize
                    row_ind, col_ind = linear_sum_assignment(-overlap)
                    
                    # col_ind[a] gives the index in ref_modes (i.e. 'b_idx') that adapted_modes[a] maps to
                    best_match = {col_ind[a]: a for a in range(m_size)}

                    # Store adapted eigenvectors and their phases
                    for b_idx, j in enumerate(manifold):
                        m_idx = best_match[b_idx]
                        adapted_mode = adapted_modes[m_idx]
                            
                        # Original logic (uncompensated D matrix)
                        # phase = np.angle(eigvals_D[m_idx])
                        # pam = (phase / angle_gen) % n_axis
                        # pam_tot[iq, j] = pam
                        
                        # 1. Uncompensated (drifting)
                        phase_uncomp = np.angle(eigvals_D[m_idx])
                        pam_tot_uncomp[iq, j] = (phase_uncomp / angle_gen) % n_axis
                        
                        # 2. Compensated (quantized for non-symmorphic)
                        phase_comp_factor = np.exp(2.0j * np.pi * np.dot(q, t_vec))
                        lambda_D_comp = eigvals_D[m_idx] * phase_comp_factor
                        phase_comp = np.angle(lambda_D_comp)
                        pam_comp = (phase_comp / angle_gen) % n_axis
                        pam_comp_err = pam_comp - np.round(pam_comp)
                        if np.abs(pam_comp_err) < 1e-2:
                            pam_comp = np.round(pam_comp) % n_axis
                        pam_tot_comp[iq, j] = pam_comp
                        
                        eigs_adapted[iq, j] = adapted_mode

    # Push the symmetry adapted vectors back to the original payload structure
    eigs_real_imag_adapted = np.zeros_like(eigs_real_imag)
    eigs_adapted = eigs_adapted.reshape(nqpoints, nphons, natoms, 3)
    eigs_real_imag_adapted[..., 0] = np.real(eigs_adapted)
    eigs_real_imag_adapted[..., 1] = np.imag(eigs_adapted)
    
    phonon_obj.eigenvectors = eigs_real_imag_adapted

    return {
        'pam_total_uncompensated': pam_tot_uncomp.tolist(),
        'pam_total_compensated': pam_tot_comp.tolist(),
    }
