import spglibFactory from './spglib.js';
import * as mat from './mat.js';

let spglibInstance = null;
let spglibPromise = null;

export async function initSpglib() {
    if (spglibInstance) return spglibInstance;
    if (!spglibPromise) {
        spglibPromise = spglibFactory();
    }
    spglibInstance = await spglibPromise;
    return spglibInstance;
}

export async function computeSymmetry(phonon) {
    let spglib = await initSpglib();
    
    let lat = phonon.lat; // 3x3 array [ [a1,a2,a3], [b1,b2,b3], [c1,c2,c3] ]
    let pos = phonon.atom_pos_red; // Nx3 array
    let types = phonon.atomic_numbers || phonon.atom_numbers; // array of N ints
    let num_atom = pos.length;
    let symprec = 1e-4; // Increased tolerance for detecting symmetries (useful in case of relaxed lattice structures)
    
    let lat_ptr = spglib._malloc(9 * 8);
    let pos_ptr = spglib._malloc(3 * num_atom * 8);
    let typ_ptr = spglib._malloc(num_atom * 4);
    
    let rot_ptr = spglib._malloc(192 * 9 * 4);
    let trans_ptr = spglib._malloc(192 * 3 * 8);
    
    let flat_lat = new Float64Array(9);
    for (let i=0; i<3; i++) {
        for (let j=0; j<3; j++) {
            flat_lat[i*3 + j] = lat[i][j];
        }
    }
    
    let flat_pos = new Float64Array(num_atom * 3);
    for (let i=0; i<num_atom; i++) {
        for (let j=0; j<3; j++) {
            flat_pos[i*3 + j] = pos[i][j];
        }
    }
    
    let flat_types = new Int32Array(num_atom);
    for (let i=0; i<num_atom; i++) {
        flat_types[i] = types[i];
    }
    
    spglib.HEAPF64.set(flat_lat, lat_ptr / 8);
    spglib.HEAPF64.set(flat_pos, pos_ptr / 8);
    spglib.HEAP32.set(flat_types, typ_ptr / 4);
    
    let num_sym = spglib.ccall(
        'get_symmetry',
        'number',
        ['number', 'number', 'number', 'number', 'number', 'number', 'number'],
        [lat_ptr, pos_ptr, typ_ptr, num_atom, symprec, rot_ptr, trans_ptr]
    );
    
    let rotations = [];
    let translations = [];
    
    let rot_out = new Int32Array(spglib.HEAP32.buffer, rot_ptr, num_sym * 9);
    let trans_out = new Float64Array(spglib.HEAPF64.buffer, trans_ptr, num_sym * 3);
    
    for (let i=0; i<num_sym; i++) {
        let R = [
            [rot_out[i*9 + 0], rot_out[i*9 + 1], rot_out[i*9 + 2]],
            [rot_out[i*9 + 3], rot_out[i*9 + 4], rot_out[i*9 + 5]],
            [rot_out[i*9 + 6], rot_out[i*9 + 7], rot_out[i*9 + 8]]
        ];
        let t = [trans_out[i*3 + 0], trans_out[i*3 + 1], trans_out[i*3 + 2]];
        rotations.push(R);
        translations.push(t);
    }
    
    spglib._free(lat_ptr);
    spglib._free(pos_ptr);
    spglib._free(typ_ptr);
    spglib._free(rot_ptr);
    spglib._free(trans_ptr);
    
    // Now we must find the little group for each segment!
    let segment_point_groups = [];
    let line_breaks = phonon.line_breaks || [[0, phonon.kpoints.length]];
    
    for (let i=0; i<line_breaks.length; i++) {
        let start = line_breaks[i][0];
        let end = line_breaks[i][1];
        if (start === undefined || end === undefined) continue;
        
        // Find segment little group from midpoint
        let mid_idx = Math.floor((start + end) / 2);
        if (mid_idx === start && end > start) {
            mid_idx = start + 1; // force it to be inside the segment
        }
        if (mid_idx >= phonon.kpoints.length) mid_idx = phonon.kpoints.length - 1;
        
        let q_mid = phonon.kpoints[mid_idx];
        
        let lg_rotations = [];
        let lg_translations = [];
        
        for (let op=0; op<rotations.length; op++) {
            let R = rotations[op];
            
            // R_q = (R^-1)^T
            let R_inv = mat.matrix_inverse(R);
            if (!R_inv) continue;
            let R_q = mat.matrix_transpose(R_inv);
            
            // q_rot = R_q * q
            let q_rot = [
                R_q[0][0]*q_mid[0] + R_q[0][1]*q_mid[1] + R_q[0][2]*q_mid[2],
                R_q[1][0]*q_mid[0] + R_q[1][1]*q_mid[1] + R_q[1][2]*q_mid[2],
                R_q[2][0]*q_mid[0] + R_q[2][1]*q_mid[1] + R_q[2][2]*q_mid[2]
            ];
            
            // diff = q_rot - q
            let diff = [
                q_rot[0] - q_mid[0],
                q_rot[1] - q_mid[1],
                q_rot[2] - q_mid[2]
            ];
            
            // is integer?
            let is_int = (Math.abs(diff[0] - Math.round(diff[0])) < symprec) &&
                         (Math.abs(diff[1] - Math.round(diff[1])) < symprec) &&
                         (Math.abs(diff[2] - Math.round(diff[2])) < symprec);
                         
            if (is_int) {
                lg_rotations.push(R);
                lg_translations.push(translations[op]);
            }
        }
        
        segment_point_groups.push({
            start: phonon.distances[start],
            end: phonon.distances[end-1],
            rotations: lg_rotations,
            translations: lg_translations,
            point_group: identifyPointGroupSymbol(lg_rotations)
        });
    }
    
    phonon.segment_point_group_list = segment_point_groups;
    phonon.segment_point_groups = segment_point_groups; // for PhononPropertyCalculator
    
    // Now compute the point group for each high-symmetry point
    let highsym_point_group_map = {};
    if (phonon.highsym_qpts && phonon.qindex) {
        for (let dist in phonon.highsym_qpts) {
            let q_idx = phonon.qindex[dist];
            if (q_idx === undefined) continue;
            let q = phonon.kpoints[q_idx];
            
            let lg_rotations = [];
            for (let op=0; op<rotations.length; op++) {
                let R = rotations[op];
                let R_inv = mat.matrix_inverse(R);
                if (!R_inv) continue;
                let R_q = mat.matrix_transpose(R_inv);
                let q_rot = [
                    R_q[0][0]*q[0] + R_q[0][1]*q[1] + R_q[0][2]*q[2],
                    R_q[1][0]*q[0] + R_q[1][1]*q[1] + R_q[1][2]*q[2],
                    R_q[2][0]*q[0] + R_q[2][1]*q[1] + R_q[2][2]*q[2]
                ];
                let diff = [ q_rot[0] - q[0], q_rot[1] - q[1], q_rot[2] - q[2] ];
                let is_int = (Math.abs(diff[0] - Math.round(diff[0])) < symprec) &&
                             (Math.abs(diff[1] - Math.round(diff[1])) < symprec) &&
                             (Math.abs(diff[2] - Math.round(diff[2])) < symprec);
                if (is_int) lg_rotations.push(R);
            }
            highsym_point_group_map[dist] = identifyPointGroupSymbol(lg_rotations);
        }
    }
    phonon.highsym_point_group_map = highsym_point_group_map;
    
    // Store full crystal operations for future 3D animations
    phonon.crystal_symmetries = {
        rotations: rotations,
        translations: translations
    };
}

export function identifyPointGroupSymbol(rotations) {
    let counts = { E:0, C2:0, C3:0, C4:0, C6:0, i:0, m:0, S3:0, S4:0, S6:0 };
    for (let R of rotations) {
        let det = Math.round(R[0][0]*(R[1][1]*R[2][2]-R[2][1]*R[1][2]) - R[0][1]*(R[1][0]*R[2][2]-R[1][2]*R[2][0]) + R[0][2]*(R[1][0]*R[2][1]-R[1][1]*R[2][0]));
        let tr = Math.round(R[0][0] + R[1][1] + R[2][2]);
        
        if (det === 1) {
            if (tr === 3) counts.E++;
            else if (tr === -1) counts.C2++;
            else if (tr === 0) counts.C3++;
            else if (tr === 1) counts.C4++;
            else if (tr === 2) counts.C6++;
        } else if (det === -1) {
            if (tr === -3) counts.i++;
            else if (tr === 1) counts.m++;
            else if (tr === -2) counts.S3++;
            else if (tr === -1) counts.S4++;
            else if (tr === 0) counts.S6++;
        }
    }
    
    const sig = `${counts.E},${counts.C2},${counts.C3},${counts.C4},${counts.C6},${counts.i},${counts.m},${counts.S3},${counts.S4},${counts.S6}`;
    
    const groups = {
        "1,0,0,0,0,0,0,0,0,0": "C1",
        "1,0,0,0,0,1,0,0,0,0": "Ci",
        "1,0,0,0,0,0,1,0,0,0": "Cs",
        "1,1,0,0,0,0,0,0,0,0": "C2",
        "1,0,2,0,0,0,0,0,0,0": "C3",
        "1,1,0,2,0,0,0,0,0,0": "C4",
        "1,1,0,0,0,0,0,0,2,0": "S4",
        "1,1,2,0,2,0,0,0,0,0": "C6",
        "1,0,2,0,0,0,1,2,0,0": "C3h",
        "1,0,2,0,0,1,0,0,0,2": "S6",
        "1,1,0,0,0,1,1,0,0,0": "C2h",
        "1,3,0,0,0,0,0,0,0,0": "D2",
        "1,1,0,0,0,0,2,0,0,0": "C2v",
        "1,1,0,2,0,1,1,0,2,0": "C4h",
        "1,5,0,2,0,0,0,0,0,0": "D4",
        "1,1,0,2,0,0,4,0,0,0": "C4v",
        "1,3,0,0,0,0,2,0,2,0": "D2d",
        "1,3,0,0,0,1,3,0,0,0": "D2h",
        "1,1,2,0,2,1,1,2,0,2": "C6h",
        "1,7,2,0,2,0,0,0,0,0": "D6",
        "1,1,2,0,2,0,6,0,0,0": "C6v",
        "1,3,2,0,0,0,4,2,0,0": "D3h",
        "1,3,2,0,0,1,3,0,0,2": "D3d",
        "1,3,2,0,0,0,0,0,0,0": "D3",
        "1,0,2,0,0,0,3,0,0,0": "C3v",
        "1,3,8,0,0,0,0,0,0,0": "T",
        "1,3,8,0,0,1,3,0,0,8": "Th",
        "1,9,8,6,0,0,0,0,0,0": "O",
        "1,3,8,0,0,0,6,0,6,0": "Td",
        "1,5,0,2,0,1,5,0,2,0": "D4h",
        "1,7,2,0,2,1,7,2,0,2": "D6h",
        "1,9,8,6,0,1,9,0,6,8": "Oh"
    };
    
    return groups[sig] || "C1"; // fallback to C1 if somehow signature doesn't match
}
