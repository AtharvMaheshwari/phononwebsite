import numpy as np
import spglib

# List of the 73 symmorphic space group numbers out of 230
SYMMORPHIC_SPACE_GROUPS = {
    1, 3, 6, 10, 16, 21, 22, 23, 25, 35, 38, 42, 47, 65, 69, 71, 75, 79, 81, 83,
    87, 89, 97, 99, 107, 111, 115, 119, 123, 139, 143, 147, 149, 150, 155, 156,
    157, 160, 162, 164, 166, 168, 174, 175, 177, 183, 187, 189, 191, 195, 197,
    200, 202, 204, 207, 209, 211, 215, 216, 217, 219, 221, 225, 229
}

def get_crystal_symmetry(cell, symprec=1e-4):
    """
    Get the spacegroup name, spacegroup number, and whether the crystal is symmorphic.
    
    Parameters:
        cell: tuple (lattice, positions, numbers)
        symprec: float, symmetry precision
        
    Returns:
        dict containing spacegroup symbol, number, and symmorphic boolean
    """
    spacegroup_str = spglib.get_spacegroup(cell, symprec=symprec)
    if not spacegroup_str:
        return {
            "spacegroup": "Unknown",
            "spacegroup_number": 0,
            "symmorphic": False
        }
        
    sg_name, sg_num_str = spacegroup_str.split(' (')
    sg_num = int(sg_num_str.rstrip(')'))
    is_symmorphic = sg_num in SYMMORPHIC_SPACE_GROUPS
    
    return {
        "spacegroup": sg_name,
        "spacegroup_number": sg_num,
        "symmorphic": is_symmorphic
    }

def get_little_group_ops(cell, q, symprec=1e-4):
    """
    Find the symmetry operations in the space group that leave wavevector q invariant
    (modulo reciprocal lattice vectors).
    
    Parameters:
        cell: tuple (lattice, positions, numbers)
        q: array-like of shape (3,), q-point in reciprocal fractional coordinates
        symprec: float, symmetry precision
        
    Returns:
        dict containing:
            'rotations': list of 3x3 integer matrices in direct space
            'translations': list of shape (3,) fractional translation vectors
            'rec_rotations': list of 3x3 float matrices representing reciprocal rotations R_q = (R^-1)^T
            'indices': list of indices matching original space group operations
    """
    symmetry = spglib.get_symmetry(cell, symprec=symprec)
    rotations = symmetry['rotations']
    translations = symmetry['translations']
    
    lg_rotations = []
    lg_translations = []
    lg_rec_rotations = []
    lg_indices = []
    
    q = np.array(q, dtype=float)
    
    for idx, R in enumerate(rotations):
        # In reciprocal space, fractional coordinates transform as R_q = (R^-1)^T
        try:
            R_inv = np.linalg.inv(R)
            R_q = R_inv.T
        except np.linalg.LinAlgError:
            continue
            
        # Rotate wavevector
        q_rot = np.dot(R_q, q)
        
        # Check modulo 1 equivalence (R_q * q - q is integer vector)
        diff = q_rot - q
        if np.all(np.isclose(diff, np.round(diff), atol=symprec)):
            lg_rotations.append(R)
            lg_translations.append(translations[idx])
            lg_rec_rotations.append(R_q)
            lg_indices.append(idx)
            
    return {
        "rotations": lg_rotations,
        "translations": lg_translations,
        "rec_rotations": lg_rec_rotations,
        "indices": lg_indices
    }

def get_little_group_point_group(cell, q, symprec=1e-4):
    """
    Determine the Hermann-Mauguin point group symbol of the little group at wavevector q.
    
    Parameters:
        cell: tuple (lattice, positions, numbers)
        q: array-like of shape (3,), wavevector
        symprec: float
        
    Returns:
        str, Hermann-Mauguin point group symbol (e.g., '3m', '4mm', '2/m')
    """
    ops = get_little_group_ops(cell, q, symprec=symprec)
    lg_rotations = np.array(ops["rotations"], dtype=int)
    
    if len(lg_rotations) == 0:
        return "1"
        
    symbol, pg_number, transformation = spglib.get_pointgroup(lg_rotations)
    hm_symbol = symbol.strip()
    
    # Map from Hermann-Mauguin to Schoenflies
    hm_to_schoenflies = {
        "1": "C1", "-1": "Ci",
        "2": "C2", "m": "Cs", "2/m": "C2h",
        "222": "D2", "mm2": "C2v", "mmm": "D2h",
        "4": "C4", "-4": "S4", "4/m": "C4h",
        "422": "D4", "4mm": "C4v", "-42m": "D2d", "4/mmm": "D4h",
        "3": "C3", "-3": "S6", "32": "D3", "3m": "C3v", "-3m": "D3d",
        "6": "C6", "-6": "C3h", "6/m": "C6h",
        "622": "D6", "6mm": "C6v", "-6m2": "D3h", "6/mmm": "D6h",
        "23": "T", "m-3": "Th", "432": "O", "-43m": "Td", "m-3m": "Oh"
    }
    
    return hm_to_schoenflies.get(hm_symbol, hm_symbol)


def get_symmetry_labels_for_path(cell, qpoints, highsym_qpts, symprec=1e-4):
    """
    Compute the point group label at every high-symmetry q-point, and at
    the midpoint of every path segment between consecutive high-symmetry points.
    
    Parameters:
        cell: tuple (lattice, positions, numbers)
        qpoints: array-like of shape (nqpoints, 3), fractional q-point coordinates
        highsym_qpts: list of [index, label] pairs from the phonon data,
                      e.g. [[0, 'GAMMA'], [20, 'X'], [41, 'M'], ...]
        symprec: float
        
    Returns:
        dict with two keys:
            'highsym_point_groups': list of dicts, one per high-symmetry point:
                { "index": int, "label": str, "point_group": str }
            'segment_point_groups': list of dicts, one per path segment:
                { "start_index": int, "end_index": int,
                  "start_label": str, "end_label": str,
                  "midpoint_q": [float, float, float],
                  "point_group": str }
    """
    qpoints = np.array(qpoints, dtype=float)
    
    # 1. Point group at each high-symmetry point
    highsym_point_groups = []
    for item in highsym_qpts:
        idx = int(item[0])
        label = str(item[1])
        q = qpoints[idx]
        pg = get_little_group_point_group(cell, q, symprec=symprec)
        highsym_point_groups.append({
            "index": idx,
            "label": label,
            "point_group": pg
        })
    
    # 2. Point group at the midpoint of each segment
    segment_point_groups = []
    for i in range(len(highsym_qpts) - 1):
        idx_start = int(highsym_qpts[i][0])
        idx_end = int(highsym_qpts[i + 1][0])
        label_start = str(highsym_qpts[i][1])
        label_end = str(highsym_qpts[i + 1][1])
        
        # Use the actual midpoint q-vector of the segment
        mid_idx = (idx_start + idx_end) // 2
        q_mid = qpoints[mid_idx]
        
        lg_ops = get_little_group_ops(cell, q_mid, symprec=symprec)
        pg_mid = get_little_group_point_group(cell, q_mid, symprec=symprec)
        segment_point_groups.append({
            "start_index": idx_start,
            "end_index": idx_end,
            "start_label": label_start,
            "end_label": label_end,
            "midpoint_q": q_mid.tolist(),
            "point_group": pg_mid,
            "rotations": [r.tolist() for r in lg_ops["rotations"]],
            "translations": [t.tolist() for t in lg_ops["translations"]]
        })
    
    return {
        "highsym_point_groups": highsym_point_groups,
        "segment_point_groups": segment_point_groups
    }
