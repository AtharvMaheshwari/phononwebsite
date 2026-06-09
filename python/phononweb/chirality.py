import numpy as np

def chirality_opr_xy(natom, iatom):
    """
    Generates Jxy(Jz) chirality projection operators (left/right-handed).
    """
    Left = np.zeros((natom * 3, 1), dtype=np.complex128)
    Right = np.zeros((natom * 3, 1), dtype=np.complex128)
    for i in range(natom):
        if i == iatom:
            Left[i * 3, 0] = np.sqrt(2) / 2.
            Left[i * 3 + 1, 0] = -1j * np.sqrt(2) / 2.
            Right[i * 3, 0] = np.sqrt(2) / 2.
            Right[i * 3 + 1, 0] = 1j * np.sqrt(2) / 2.
    return Left, Right

def chirality_opr_yz(natom, iatom):
    """
    Generates Jyz(Jx) chirality projection operators (left/right-handed).
    """
    Left = np.zeros((natom * 3, 1), dtype=np.complex128)
    Right = np.zeros((natom * 3, 1), dtype=np.complex128)
    for i in range(natom):
        if i == iatom:
            Left[i * 3 + 1, 0] = np.sqrt(2) / 2.
            Left[i * 3 + 2, 0] = -1j * np.sqrt(2) / 2.
            Right[i * 3 + 1, 0] = np.sqrt(2) / 2.
            Right[i * 3 + 2, 0] = 1j * np.sqrt(2) / 2.
    return Left, Right

def chirality_opr_zx(natom, iatom):
    """
    Generates Jzx(Jy) chirality projection operators (left/right-handed).
    """
    Left = np.zeros((natom * 3, 1), dtype=np.complex128)
    Right = np.zeros((natom * 3, 1), dtype=np.complex128)
    for i in range(natom):
        if i == iatom:
            Left[i * 3 + 2, 0] = np.sqrt(2) / 2.
            Left[i * 3 + 0, 0] = -1j * np.sqrt(2) / 2.
            Right[i * 3 + 2, 0] = np.sqrt(2) / 2.
            Right[i * 3 + 0, 0] = 1j * np.sqrt(2) / 2.
    return Left, Right

def compute_chiral_properties(phonon_obj):
    """
    Compute Phonon Angular Momentum (Jx, Jy, Jz), Helicity, and Cycloidicity
    for all k-points and all bands.
    
    Args:
        phonon_obj: instance of Phonon containing .eigenvectors, .qpoints
    
    Returns:
        dict with angular_momentum, helicity, cycloidicity
    """
    if not hasattr(phonon_obj, 'eigenvectors') or phonon_obj.eigenvectors is None:
        return {}
        
    eigs_real_imag = phonon_obj.eigenvectors
    nqpoints, nphons, natoms, _, _ = eigs_real_imag.shape
    
    # Reshape and convert to complex
    # eigs_real_imag is (nqpoints, nphons, natoms, 3, 2)
    # Target shape: (nqpoints, nphons, natoms*3)
    eigs = eigs_real_imag[..., 0] + 1j * eigs_real_imag[..., 1]
    eigs = eigs.reshape(nqpoints, nphons, natoms * 3)
    
    # Initialize output arrays
    # Jx, Jy, Jz
    angular_momentum = np.zeros((nqpoints, nphons, 3))
    
    # Pre-generate operators for all atoms
    opr_xy = [chirality_opr_xy(natoms, i) for i in range(natoms)] # gives Jz
    opr_yz = [chirality_opr_yz(natoms, i) for i in range(natoms)] # gives Jx
    opr_zx = [chirality_opr_zx(natoms, i) for i in range(natoms)] # gives Jy
    
    for iq in range(nqpoints):
        for ibnd in range(nphons):
            eig_vec = eigs[iq, ibnd, :] # shape (natoms*3,)
            
            cmat_x, cmat_y, cmat_z = 0, 0, 0
            
            for iat in range(natoms):
                # Jz
                L_z, R_z = opr_xy[iat]
                tmpl_z = np.dot(np.conj(eig_vec).T, L_z)
                tmpr_z = np.dot(np.conj(eig_vec).T, R_z)
                cmat_z += (tmpr_z[0] * np.conj(tmpr_z[0]) - tmpl_z[0] * np.conj(tmpl_z[0]))
                
                # Jx
                L_x, R_x = opr_yz[iat]
                tmpl_x = np.dot(np.conj(eig_vec).T, L_x)
                tmpr_x = np.dot(np.conj(eig_vec).T, R_x)
                cmat_x += (tmpr_x[0] * np.conj(tmpr_x[0]) - tmpl_x[0] * np.conj(tmpl_x[0]))
                
                # Jy
                L_y, R_y = opr_zx[iat]
                tmpl_y = np.dot(np.conj(eig_vec).T, L_y)
                tmpr_y = np.dot(np.conj(eig_vec).T, R_y)
                cmat_y += (tmpr_y[0] * np.conj(tmpr_y[0]) - tmpl_y[0] * np.conj(tmpl_y[0]))
                
            angular_momentum[iq, ibnd, 0] = np.real(cmat_x)
            angular_momentum[iq, ibnd, 1] = np.real(cmat_y)
            angular_momentum[iq, ibnd, 2] = np.real(cmat_z)
            
    # Helicity = J dot k_unit
    # Cycloidicity = k_unit x J
    helicity = np.zeros((nqpoints, nphons))
    cycloidicity = np.zeros((nqpoints, nphons, 4)) # Cx, Cy, Cz, |C|
    
    qpoints = phonon_obj.qpoints
    
    for iq in range(nqpoints):
        q = qpoints[iq]
        q_norm = np.linalg.norm(q)
        if q_norm > 1e-6:
            q_unit = q / q_norm
        else:
            q_unit = np.array([0.0, 0.0, 0.0])
            
        for ibnd in range(nphons):
            J = angular_momentum[iq, ibnd, :]
            
            # Helicity
            helicity[iq, ibnd] = np.dot(J, q_unit)
            
            # Cycloidicity
            C = np.cross(q_unit, J)
            cycloidicity[iq, ibnd, 0] = C[0]
            cycloidicity[iq, ibnd, 1] = C[1]
            cycloidicity[iq, ibnd, 2] = C[2]
            cycloidicity[iq, ibnd, 3] = np.linalg.norm(C)
            
    return {
        'angular_momentum_x': angular_momentum[:, :, 0].tolist(),
        'angular_momentum_y': angular_momentum[:, :, 1].tolist(),
        'angular_momentum_z': angular_momentum[:, :, 2].tolist(),
        'helicity': helicity.tolist(),
        'cycloidicity_x': cycloidicity[:, :, 0].tolist(),
        'cycloidicity_y': cycloidicity[:, :, 1].tolist(),
        'cycloidicity_z': cycloidicity[:, :, 2].tolist()
    }
