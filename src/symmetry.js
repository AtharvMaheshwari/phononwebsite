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
    let symprec = 1e-5;
    
    let lat_ptr = spglib._malloc(9 * 8);
    let pos_ptr = spglib._malloc(3 * num_atom * 8);
    let typ_ptr = spglib._malloc(num_atom * 4);
    
    let rot_ptr = spglib._malloc(48 * 9 * 4);
    let trans_ptr = spglib._malloc(48 * 3 * 8);
    
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
        
        // Find segment little group from midpoint
        let mid_idx = Math.floor((start + end) / 2);
        if (mid_idx === start && end > start) {
            mid_idx = start + 1; // force it to be inside the segment
        }
        if (mid_idx >= phonon.kpoints.length) mid_idx = phonon.kpoints.length - 1;
        
        let q_mid = phonon.kpoints[mid_idx];
        
        let lg_rotations = [];
        let lg_translations = [];
        
        for (let op=0; op<num_sym; op++) {
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
            translations: lg_translations
        });
    }
    
    phonon.segment_point_group_list = segment_point_groups;
    phonon.segment_point_groups = segment_point_groups; // for PhononPropertyCalculator
}
