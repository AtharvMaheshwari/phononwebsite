export class PhononPropertyCalculator {
    constructor(phononJson) {
        this.phonon = phononJson;
    }

    // --- Helper Functions ---

    complexDot(a, b) {
        let real = 0, imag = 0;
        for (let i = 0; i < a.length; i++) {
            let ar = a[i][0], ai = -a[i][1]; // conjugate
            let br = b[i][0], bi = b[i][1];
            real += ar * br - ai * bi;
            imag += ar * bi + ai * br;
        }
        return [real, imag];
    }

    chiralityOprXy(natom, iatom) {
        let L = Array(natom * 3).fill(0).map(() => [0, 0]);
        let R = Array(natom * 3).fill(0).map(() => [0, 0]);
        let val = Math.SQRT2 / 2;
        L[iatom * 3] = [val, 0];
        L[iatom * 3 + 1] = [0, -val];
        R[iatom * 3] = [val, 0];
        R[iatom * 3 + 1] = [0, val];
        return { L, R };
    }

    chiralityOprYz(natom, iatom) {
        let L = Array(natom * 3).fill(0).map(() => [0, 0]);
        let R = Array(natom * 3).fill(0).map(() => [0, 0]);
        let val = Math.SQRT2 / 2;
        L[iatom * 3 + 1] = [val, 0];
        L[iatom * 3 + 2] = [0, -val];
        R[iatom * 3 + 1] = [val, 0];
        R[iatom * 3 + 2] = [0, val];
        return { L, R };
    }

    chiralityOprZx(natom, iatom) {
        let L = Array(natom * 3).fill(0).map(() => [0, 0]);
        let R = Array(natom * 3).fill(0).map(() => [0, 0]);
        let val = Math.SQRT2 / 2;
        L[iatom * 3 + 2] = [val, 0];
        L[iatom * 3 + 0] = [0, -val];
        R[iatom * 3 + 2] = [val, 0];
        R[iatom * 3 + 0] = [0, val];
        return { L, R };
    }

    getRecLat(lat) {
        let cross = (a, b) => [
            a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0]
        ];
        let dot = (a, b) => a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
        let a1 = lat[0], a2 = lat[1], a3 = lat[2];
        let b1 = cross(a2, a3);
        let b2 = cross(a3, a1);
        let b3 = cross(a1, a2);
        let v = dot(a1, b1);
        return [ b1.map(x => x / v), b2.map(x => x / v), b3.map(x => x / v) ];
    }

    redCar(vec, lat) {
        let r = [0, 0, 0];
        for (let i = 0; i < 3; i++) {
            r[i] = lat[0][i]*vec[0] + lat[1][i]*vec[1] + lat[2][i]*vec[2];
        }
        return r;
    }

    matMul3x3(A, B) {
        let C = [[0,0,0], [0,0,0], [0,0,0]];
        for(let i=0; i<3; i++)
            for(let j=0; j<3; j++)
                for(let k=0; k<3; k++)
                    C[i][j] += A[i][k] * B[k][j];
        return C;
    }

    matInv3x3(m) {
        let det = m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) -
                  m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                  m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        return [
            [(m[1][1] * m[2][2] - m[2][1] * m[1][2]) / det, (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / det, (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / det],
            [(m[1][2] * m[2][0] - m[1][0] * m[2][2]) / det, (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / det, (m[1][0] * m[0][2] - m[0][0] * m[1][2]) / det],
            [(m[1][0] * m[2][1] - m[2][0] * m[1][1]) / det, (m[2][0] * m[0][1] - m[0][0] * m[2][1]) / det, (m[0][0] * m[1][1] - m[1][0] * m[0][1]) / det]
        ];
    }

    getRotationAxisAngle(R_car) {
        let tr = R_car[0][0] + R_car[1][1] + R_car[2][2];
        let cos_val = Math.max(-1.0, Math.min(1.0, (tr - 1.0) / 2.0));
        let angle = Math.acos(cos_val);
        
        if (Math.abs(angle) < 1e-5) return { axis: [0, 0, 1], angle: 0.0 };
        
        if (Math.abs(angle - Math.PI) < 1e-5) {
            let R_plus_I = [
                [R_car[0][0] + 1, R_car[0][1], R_car[0][2]],
                [R_car[1][0], R_car[1][1] + 1, R_car[1][2]],
                [R_car[2][0], R_car[2][1], R_car[2][2] + 1]
            ];
            for (let i = 0; i < 3; i++) {
                let col = [R_plus_I[0][i], R_plus_I[1][i], R_plus_I[2][i]];
                let norm = Math.hypot(...col);
                if (norm > 1e-5) return { axis: col.map(x => x / norm), angle: Math.PI };
            }
        }
        
        let axis = [
            R_car[2][1] - R_car[1][2],
            R_car[0][2] - R_car[2][0],
            R_car[1][0] - R_car[0][1]
        ];
        let norm = Math.hypot(...axis);
        return norm > 1e-5 ? { axis: axis.map(x => x / norm), angle } : { axis: [0, 0, 1], angle };
    }

    buildRepresentationMatrices(cell, q, R, tau) {
        let lattice = cell[0];
        let pos_red = cell[1];
        let natoms = pos_red.length;
        
        let A = [
            [lattice[0][0], lattice[1][0], lattice[2][0]],
            [lattice[0][1], lattice[1][1], lattice[2][1]],
            [lattice[0][2], lattice[1][2], lattice[2][2]]
        ];
        let A_inv = this.matInv3x3(A);
        let R_car = this.matMul3x3(A, this.matMul3x3(R, A_inv));
        
        let D = Array(3 * natoms).fill(0).map(() => Array(3 * natoms).fill(0).map(() => [0, 0]));
        
        for (let i = 0; i < natoms; i++) {
            let r_prime = [0, 0, 0];
            for(let x = 0; x < 3; x++) r_prime[x] = R[x][0]*pos_red[i][0] + R[x][1]*pos_red[i][1] + R[x][2]*pos_red[i][2] + tau[x];
            
            let min_dist = 1e10, j_best = -1;
            for (let j = 0; j < natoms; j++) {
                let diff_mod1 = [
                    (pos_red[j][0] - r_prime[0]) - Math.round(pos_red[j][0] - r_prime[0]),
                    (pos_red[j][1] - r_prime[1]) - Math.round(pos_red[j][1] - r_prime[1]),
                    (pos_red[j][2] - r_prime[2]) - Math.round(pos_red[j][2] - r_prime[2])
                ];
                let diff_car = [0, 0, 0];
                for(let x=0; x<3; x++) diff_car[x] = diff_mod1[0]*lattice[0][x] + diff_mod1[1]*lattice[1][x] + diff_mod1[2]*lattice[2][x];
                let dist = Math.hypot(...diff_car);
                if (dist < min_dist) { min_dist = dist; j_best = j; }
            }
            
            if (min_dist > 1e-3) throw new Error("Symmetry mapping failed for atom " + i);
            
            let Rq = [
                R[0][0]*q[0] + R[0][1]*q[1] + R[0][2]*q[2],
                R[1][0]*q[0] + R[1][1]*q[1] + R[1][2]*q[2],
                R[2][0]*q[0] + R[2][1]*q[1] + R[2][2]*q[2]
            ];
            let G = [Math.round(Rq[0] - q[0]), Math.round(Rq[1] - q[1]), Math.round(Rq[2] - q[2])];
            
            let dotG = G[0]*pos_red[j_best][0] + G[1]*pos_red[j_best][1] + G[2]*pos_red[j_best][2];
            let phase_val = -2.0 * Math.PI * dotG;
            let phase = [Math.cos(phase_val), Math.sin(phase_val)];
            
            for (let r = 0; r < 3; r++) {
                for (let c = 0; c < 3; c++) {
                    D[3 * j_best + r][3 * i + c] = [phase[0] * R_car[r][c], phase[1] * R_car[r][c]];
                }
            }
        }
        return D;
    }

    greedyAssignment(matrix, maximize = false) {
        let n = matrix.length;
        if (n === 0) return [];
        let m = matrix[0].length;
        let row_ind = [], col_ind = [];
        let used_cols = Array(m).fill(false), used_rows = Array(n).fill(false);
        
        let elements = [];
        for (let i = 0; i < n; i++) for (let j = 0; j < m; j++) elements.push({ r: i, c: j, v: matrix[i][j] });
        elements.sort((a, b) => maximize ? b.v - a.v : a.v - b.v);
        
        for (let el of elements) {
            if (!used_rows[el.r] && !used_cols[el.c]) {
                row_ind.push(el.r); col_ind.push(el.c);
                used_rows[el.r] = true; used_cols[el.c] = true;
            }
            if (row_ind.length === n) break;
        }
        
        for(let i = 0; i < n; i++) {
            if (!used_rows[i]) {
                for(let j = 0; j < m; j++) {
                    if (!used_cols[j]) {
                        row_ind.push(i); col_ind.push(j);
                        used_rows[i] = true; used_cols[j] = true;
                        break;
                    }
                }
            }
        }
        
        return row_ind.map((r, i) => [r, col_ind[i]]).sort((a, b) => a[0] - b[0]).map(pair => pair[1]);
    }

    realJacobi(M) {
        let n = M.length;
        let V = Array(n).fill(0).map((_, i) => Array(n).fill(0).map((_, j) => i === j ? 1 : 0));
        let A = JSON.parse(JSON.stringify(M));
        
        for (let iter = 0; iter < 100 * n * n; iter++) {
            let maxVal = 0, p = 0, q = 1;
            for (let i = 0; i < n; i++) {
                for (let j = i + 1; j < n; j++) {
                    if (Math.abs(A[i][j]) > maxVal) { maxVal = Math.abs(A[i][j]); p = i; q = j; }
                }
            }
            if (maxVal < 1e-10) break;
            
            let app = A[p][p], aqq = A[q][q], apq = A[p][q];
            let phi = 0.5 * Math.atan2(2 * apq, app - aqq);
            let c = Math.cos(phi), s = Math.sin(phi);
            
            for (let i = 0; i < n; i++) {
                if (i !== p && i !== q) {
                    let aip = A[i][p], aiq = A[i][q];
                    A[i][p] = A[p][i] = c * aip - s * aiq;
                    A[i][q] = A[q][i] = s * aip + c * aiq;
                }
                let vip = V[i][p], viq = V[i][q];
                V[i][p] = c * vip - s * viq;
                V[i][q] = s * vip + c * viq;
            }
            A[p][p] = c * c * app + 2 * s * c * apq + s * s * aqq;
            A[q][q] = s * s * app - 2 * s * c * apq + c * c * aqq;
            A[p][q] = A[q][p] = 0;
        }
        return A.map((row, i) => ({ val: row[i], vec: V.map(vRow => vRow[i]) })).sort((a, b) => b.val - a.val);
    }

    diagonalizeComplexUnitary(M) {
        let n = M.length;
        let R = Array(2 * n).fill(0).map(() => Array(2 * n).fill(0));
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < n; j++) {
                let r = M[i][j][0] + M[j][i][0] - 0.5 * (M[i][j][1] + M[j][i][1]);
                let im = M[i][j][1] - M[j][i][1] + 0.5 * (M[i][j][0] - M[j][i][0]);
                R[i][j] = R[i + n][j + n] = r;
                R[i][j + n] = -im; R[i + n][j] = im;
            }
        }
        let evals = this.realJacobi(R);
        let evecs = [];
        for (let k = 0; k < 2 * n; k += 2) {
            evecs.push(Array.from({length: n}, (_, i) => [evals[k].vec[i], evals[k].vec[i + n]]));
        }
        
        let M_evals = evecs.map(v => {
            let dotR = 0, dotI = 0;
            for (let i = 0; i < n; i++) {
                let mvr = 0, mvi = 0;
                for (let j = 0; j < n; j++) {
                    mvr += M[i][j][0]*v[j][0] - M[i][j][1]*v[j][1];
                    mvi += M[i][j][0]*v[j][1] + M[i][j][1]*v[j][0];
                }
                dotR += v[i][0]*mvr + v[i][1]*mvi;
                dotI += v[i][0]*mvi - v[i][1]*mvr;
            }
            return [dotR, dotI];
        });
        return { evals: M_evals, evecs };
    }

    // --- Core API ---

    computeChiralProperties() {
        if (!this.phonon.vec) return {};
        
        let eigs = this.phonon.vec;
        let nqpoints = eigs.length, nphons = eigs[0].length, natoms = eigs[0][0].length;
        
        let angular_momentum = Array(nqpoints).fill(0).map(() => Array(nphons).fill(0).map(() => [0,0,0]));
        let helicity = Array(nqpoints).fill(0).map(() => Array(nphons).fill(0));
        let cycloidicity = Array(nqpoints).fill(0).map(() => Array(nphons).fill(0).map(() => [0,0,0,0]));
        
        let opr_xy = Array.from({length: natoms}, (_, i) => this.chiralityOprXy(natoms, i));
        let opr_yz = Array.from({length: natoms}, (_, i) => this.chiralityOprYz(natoms, i));
        let opr_zx = Array.from({length: natoms}, (_, i) => this.chiralityOprZx(natoms, i));
        
        for (let iq = 0; iq < nqpoints; iq++) {
            for (let ibnd = 0; ibnd < nphons; ibnd++) {
                let eig_vec = [];
                for (let iat = 0; iat < natoms; iat++) {
                    for (let x = 0; x < 3; x++) eig_vec.push(eigs[iq][ibnd][iat][x]);
                }
                
                let cmat = [0, 0, 0];
                for (let iat = 0; iat < natoms; iat++) {
                    let Lz = this.complexDot(eig_vec, opr_xy[iat].L), Rz = this.complexDot(eig_vec, opr_xy[iat].R);
                    cmat[2] += (Rz[0]**2 + Rz[1]**2) - (Lz[0]**2 + Lz[1]**2);
                    
                    let Lx = this.complexDot(eig_vec, opr_yz[iat].L), Rx = this.complexDot(eig_vec, opr_yz[iat].R);
                    cmat[0] += (Rx[0]**2 + Rx[1]**2) - (Lx[0]**2 + Lx[1]**2);
                    
                    let Ly = this.complexDot(eig_vec, opr_zx[iat].L), Ry = this.complexDot(eig_vec, opr_zx[iat].R);
                    cmat[1] += (Ry[0]**2 + Ry[1]**2) - (Ly[0]**2 + Ly[1]**2);
                }
                angular_momentum[iq][ibnd] = cmat;
            }
        }
        
        let rec = this.getRecLat(this.phonon.cell || this.phonon.lat);
        for (let iq = 0; iq < nqpoints; iq++) {
            let q_car = this.redCar(this.phonon.kpoints[iq], rec);
            let q_norm = Math.hypot(...q_car);
            let q_unit = q_norm > 1e-6 ? q_car.map(x => x / q_norm) : [0,0,0];
            
            for (let ibnd = 0; ibnd < nphons; ibnd++) {
                let J = angular_momentum[iq][ibnd];
                helicity[iq][ibnd] = J[0]*q_unit[0] + J[1]*q_unit[1] + J[2]*q_unit[2];
                
                let C = [
                    q_unit[1]*J[2] - q_unit[2]*J[1],
                    q_unit[2]*J[0] - q_unit[0]*J[2],
                    q_unit[0]*J[1] - q_unit[1]*J[0]
                ];
                cycloidicity[iq][ibnd] = [C[0], C[1], C[2], Math.hypot(...C)];
            }
        }
        
        return {
            angular_momentum_x: angular_momentum.map(q => q.map(b => b[0])),
            angular_momentum_y: angular_momentum.map(q => q.map(b => b[1])),
            angular_momentum_z: angular_momentum.map(q => q.map(b => b[2])),
            helicity,
            cycloidicity_x: cycloidicity.map(q => q.map(b => b[0])),
            cycloidicity_y: cycloidicity.map(q => q.map(b => b[1])),
            cycloidicity_z: cycloidicity.map(q => q.map(b => b[2]))
        };
    }

    computePamProperties() {
        if (!this.phonon.vec) return {};
        
        let eigs = this.phonon.vec;
        let nqpoints = eigs.length, nphons = eigs[0].length, natoms = eigs[0][0].length;
        
        let eigs_flat = eigs.map(q_eigs => q_eigs.map(bnd => {
            let flat = [];
            for (let iat = 0; iat < natoms; iat++) for (let x = 0; x < 3; x++) flat.push([...bnd[iat][x]]);
            return flat;
        }));
        
        let freqs = JSON.parse(JSON.stringify(this.phonon.eigenvalues));
        let qpoints = this.phonon.kpoints;
        let lattice = this.phonon.cell || this.phonon.lat;
        let pos_red = this.phonon.atom_pos_red || this.phonon.atom_pos_car; 
        let cell = [lattice, pos_red, null];
        
        let pam_uncomp = Array(nqpoints).fill(0).map(() => Array(nphons).fill(0));
        let pam_comp = Array(nqpoints).fill(0).map(() => Array(nphons).fill(0));
        let eigs_adapted = JSON.parse(JSON.stringify(eigs_flat));
        
        let actual_segments = this.phonon.line_breaks || [[0, nqpoints]];
        let B_rec = this.getRecLat(lattice);
        let A = [[lattice[0][0], lattice[1][0], lattice[2][0]], [lattice[0][1], lattice[1][1], lattice[2][1]], [lattice[0][2], lattice[1][2], lattice[2][2]]];
        let A_inv = this.matInv3x3(A);
        
        for (let seg_idx = 0; seg_idx < actual_segments.length; seg_idx++) {
            let start_idx = actual_segments[seg_idx][0], end_idx = actual_segments[seg_idx][1];
            let lg_ops_segment = (this.phonon.segment_point_groups || [])[seg_idx];
            if (!lg_ops_segment) continue;
            
            let q_start = qpoints[start_idx], q_end = qpoints[Math.min(end_idx, nqpoints - 1)];
            let dq_car = this.redCar([q_end[0]-q_start[0], q_end[1]-q_start[1], q_end[2]-q_start[2]], B_rec);
            let norm_dq = Math.hypot(...dq_car);
            let dir_q_car = norm_dq > 1e-5 ? dq_car.map(x => x / norm_dq) : [0, 0, 1];
            
            let best_op = -1, min_angle = 10.0;
            for (let idx = 0; idx < lg_ops_segment.rotations.length; idx++) {
                let R = lg_ops_segment.rotations[idx];
                let det = Math.round(R[0][0]*(R[1][1]*R[2][2]-R[2][1]*R[1][2]) - R[0][1]*(R[1][0]*R[2][2]-R[1][2]*R[2][0]) + R[0][2]*(R[1][0]*R[2][1]-R[1][1]*R[2][0]));
                if (det > 0) {
                    let R_car = this.matMul3x3(A, this.matMul3x3(R, A_inv));
                    let { axis, angle } = this.getRotationAxisAngle(R_car);
                    if (angle > 1e-3 && Math.abs(Math.abs(axis[0]*dir_q_car[0] + axis[1]*dir_q_car[1] + axis[2]*dir_q_car[2]) - 1.0) < 1e-2) {
                        if (angle < min_angle) { min_angle = angle; best_op = idx; }
                    }
                }
            }
            if (best_op === -1) continue;
            
            let R_gen = lg_ops_segment.rotations[best_op], tau_gen = lg_ops_segment.translations[best_op];
            let angle_gen = this.getRotationAxisAngle(this.matMul3x3(A, this.matMul3x3(R_gen, A_inv))).angle;
            let n_axis = angle_gen > 1e-4 ? Math.round(2.0 * Math.PI / angle_gen) : 1;
            
            let t_vec = [0, 0, 0], R_pow = [[1,0,0],[0,1,0],[0,0,1]];
            for (let i = 0; i < n_axis; i++) {
                for(let x=0; x<3; x++) t_vec[x] += R_pow[x][0]*tau_gen[0] + R_pow[x][1]*tau_gen[1] + R_pow[x][2]*tau_gen[2];
                R_pow = this.matMul3x3(R_gen, R_pow);
            }
            t_vec = t_vec.map(x => x / n_axis);
            
            for (let iq = start_idx; iq < end_idx; iq++) {
                if (iq >= nqpoints) continue;
                
                if (iq > start_idx) {
                    let overlap = Array(nphons).fill(0).map((_, i) => Array(nphons).fill(0).map((_, j) => {
                        let dot = this.complexDot(eigs_adapted[iq-1][i], eigs_flat[iq][j]);
                        return Math.hypot(dot[0], dot[1]);
                    }));
                    let col_ind = this.greedyAssignment(overlap, true);
                    eigs_flat[iq] = col_ind.map(j => eigs_flat[iq][j]);
                    freqs[iq] = col_ind.map(j => freqs[iq][j]);
                }
                
                let q = qpoints[iq], w = freqs[iq];
                let D;
                try { D = this.buildRepresentationMatrices(cell, q, R_gen, tau_gen); } catch(e) { eigs_adapted[iq] = eigs_flat[iq].slice(); continue; }
                
                let visited = Array(nphons).fill(false);
                for (let ibnd = 0; ibnd < nphons; ibnd++) {
                    if (visited[ibnd]) continue;
                    let manifold = [];
                    for (let j = 0; j < nphons; j++) if (Math.abs(w[j] - w[ibnd]) < 1e-5) { manifold.push(j); visited[j] = true; }
                    
                    let basis = manifold.map(j => eigs_flat[iq][j]), m_size = manifold.length;
                    
                    if (m_size === 1) {
                        let De = Array(3*natoms).fill(0).map(() => [0, 0]);
                        for(let a=0; a<3*natoms; a++) for(let b=0; b<3*natoms; b++) {
                            De[a][0] += D[a][b][0]*basis[0][b][0] - D[a][b][1]*basis[0][b][1];
                            De[a][1] += D[a][b][0]*basis[0][b][1] + D[a][b][1]*basis[0][b][0];
                        }
                        let lambda_D = this.complexDot(basis[0], De);
                        
                        pam_uncomp[iq][ibnd] = ((Math.atan2(lambda_D[1], lambda_D[0]) / angle_gen) % n_axis + n_axis) % n_axis;
                        let dot_qt = q[0]*t_vec[0] + q[1]*t_vec[1] + q[2]*t_vec[2];
                        let pc = [Math.cos(2*Math.PI*dot_qt), Math.sin(2*Math.PI*dot_qt)];
                        let ld_comp = [lambda_D[0]*pc[0] - lambda_D[1]*pc[1], lambda_D[0]*pc[1] + lambda_D[1]*pc[0]];
                        let p_comp = ((Math.atan2(ld_comp[1], ld_comp[0]) / angle_gen) % n_axis + n_axis) % n_axis;
                        pam_comp[iq][ibnd] = Math.abs(p_comp - Math.round(p_comp)) < 1e-2 ? (Math.round(p_comp) % n_axis + n_axis) % n_axis : p_comp;
                        eigs_adapted[iq][ibnd] = basis[0];
                    } else {
                        let M_D = Array(m_size).fill(0).map((_, a) => Array(m_size).fill(0).map((_, b) => {
                            let Db = Array(3*natoms).fill(0).map(() => [0, 0]);
                            for(let r=0; r<3*natoms; r++) for(let c=0; c<3*natoms; c++) {
                                Db[r][0] += D[r][c][0]*basis[b][c][0] - D[r][c][1]*basis[b][c][1];
                                Db[r][1] += D[r][c][0]*basis[b][c][1] + D[r][c][1]*basis[b][c][0];
                            }
                            return this.complexDot(basis[a], Db);
                        }));
                        let { evals: eigvals_D, evecs: eigvecs_D } = this.diagonalizeComplexUnitary(M_D);
                        
                        let adapted = eigvecs_D.map(V => {
                            let mode = Array(3*natoms).fill(0).map(() => [0, 0]);
                            for(let a=0; a<m_size; a++) for(let c=0; c<3*natoms; c++) {
                                mode[c][0] += V[a][0]*basis[a][c][0] - V[a][1]*basis[a][c][1];
                                mode[c][1] += V[a][0]*basis[a][c][1] + V[a][1]*basis[a][c][0];
                            }
                            let norm = Math.sqrt(mode.reduce((sum, v) => sum + v[0]**2 + v[1]**2, 0));
                            return norm > 1e-12 ? mode.map(v => [v[0]/norm, v[1]/norm]) : mode;
                        });
                        
                        let ref_modes = iq > start_idx ? manifold.map(j => eigs_adapted[iq - 1][j]) : basis;
                        let overlap = adapted.map(a_mode => ref_modes.map(r_mode => {
                            let dot = this.complexDot(a_mode, r_mode);
                            return Math.hypot(dot[0], dot[1]);
                        }));
                        let match_ind = this.greedyAssignment(overlap, true);
                        
                        for (let b_idx = 0; b_idx < m_size; b_idx++) {
                            let j = manifold[b_idx], m_idx = match_ind.indexOf(b_idx);
                            pam_uncomp[iq][j] = ((Math.atan2(eigvals_D[m_idx][1], eigvals_D[m_idx][0]) / angle_gen) % n_axis + n_axis) % n_axis;
                            let dot_qt = q[0]*t_vec[0] + q[1]*t_vec[1] + q[2]*t_vec[2];
                            let pc = [Math.cos(2*Math.PI*dot_qt), Math.sin(2*Math.PI*dot_qt)];
                            let ld = eigvals_D[m_idx];
                            let ld_comp = [ld[0]*pc[0] - ld[1]*pc[1], ld[0]*pc[1] + ld[1]*pc[0]];
                            let p_comp = ((Math.atan2(ld_comp[1], ld_comp[0]) / angle_gen) % n_axis + n_axis) % n_axis;
                            pam_comp[iq][j] = Math.abs(p_comp - Math.round(p_comp)) < 1e-2 ? (Math.round(p_comp) % n_axis + n_axis) % n_axis : p_comp;
                            eigs_adapted[iq][j] = adapted[m_idx];
                        }
                    }
                }
            }
        }
        
        // Push the tracked eigenvalues and vectors back!
        this.phonon.eigenvalues = freqs;
        let eigs_real_imag = Array(nqpoints).fill(0).map((_, iq) => Array(nphons).fill(0).map((_, ibnd) => {
            let atom_vec = Array(natoms).fill(0).map(() => Array(3).fill(0).map(() => [0, 0]));
            for (let iat = 0; iat < natoms; iat++) {
                for (let x = 0; x < 3; x++) {
                    atom_vec[iat][x] = [...eigs_adapted[iq][ibnd][iat*3 + x]];
                }
            }
            return atom_vec;
        }));
        this.phonon.vec = eigs_real_imag;
        
        return { pam_total_uncompensated: pam_uncomp, pam_total_compensated: pam_comp };
    }
}
