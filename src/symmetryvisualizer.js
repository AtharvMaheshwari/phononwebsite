/**
 * SymmetryVisualizer — Animates crystallographic symmetry operations
 * on the main VibCrystal viewer.
 *
 * Architecture:
 *   - Reads raw fractional symmetry operations from phonon.crystal_symmetries
 *   - Converts them to Cartesian space using the lattice vectors
 *   - Applies Gram-Schmidt to guarantee orthonormality (avoids NaN crashes)
 *   - Extracts axis + angle, uses SLERP for smooth rotation animation
 *   - Creates a semi-transparent "ghost" lattice as a visual reference
 *   - Drives the animation via a 0→1 slider value
 */

import * as THREE from 'three';
import * as mat from './mat.js';

export class SymmetryVisualizer {

    constructor(vibcrystal) {
        /** @type {VibCrystal} */
        this.crystal = vibcrystal;
        this.phonon = null;

        // State
        this.active = false;
        this.currentOpIndex = -1;
        this.sliderValue = 0.0; // 0 → 1

        // Ghost lattice meshes (THREE objects added to scene)
        this.ghostMeshes = [];

        // Saved state to restore on deactivation
        this.savedAmplitude = 0;
        this.savedPaused = false;

        // Precomputed Cartesian operations
        this.cartesianOps = []; // [{ R_cart, t_cart, axis, angle, det, label }]

        // DOM references (set externally)
        this.dropdownEl = null;
        this.sliderEl = null;
        this.panelEl = null;
        this.labelEl = null;
        this.toggleBtn = null;
    }

    // ─────────────────────────────────────────────
    // DOM BINDING
    // ─────────────────────────────────────────────

    bindDOM(panelEl, dropdownEl, sliderEl, labelEl, toggleBtn) {
        this.panelEl = panelEl;
        this.dropdownEl = dropdownEl;
        this.sliderEl = sliderEl;
        this.labelEl = labelEl;
        this.toggleBtn = toggleBtn;

        // Toggle button shows/hides the panel
        if (this.toggleBtn) {
            this.toggleBtn.on('click', () => {
                if (this.active) {
                    this.deactivate();
                } else {
                    this.activate();
                }
            });
        }

        // Dropdown: select a symmetry operation
        if (this.dropdownEl) {
            this.dropdownEl.on('change', () => {
                let idx = parseInt(this.dropdownEl.val(), 10);
                if (Number.isFinite(idx) && idx >= 0) {
                    this.selectOperation(idx);
                }
            });
        }

        // Slider: animate the operation
        if (this.sliderEl) {
            this.sliderEl.on('input', () => {
                let t = parseFloat(this.sliderEl.val()) / 100.0;
                this.setSliderValue(t);
            });
        }
    }

    // ─────────────────────────────────────────────
    // ACTIVATION / DEACTIVATION
    // ─────────────────────────────────────────────

    activate() {
        if (!this.crystal || !this.crystal.phonon) return;
        this.phonon = this.crystal.phonon;

        if (!this.phonon.crystal_symmetries) {
            console.warn('SymmetryVisualizer: No crystal_symmetries available.');
            return;
        }

        this.active = true;
        this.crystal.symmetryAnimationActive = true;

        // Freeze phonon vibrations
        this.savedAmplitude = this.crystal.amplitude;
        this.savedPaused = this.crystal.paused;
        this.crystal.amplitude = 0;
        this.crystal.paused = false; // Keep rendering but with 0 amplitude

        // Precompute Cartesian operations
        this.precomputeOperations();

        // Populate the dropdown
        this.populateDropdown();

        // Show UI panel
        if (this.panelEl) this.panelEl.show();
        if (this.toggleBtn) {
            this.toggleBtn.text('✕ Close');
            this.toggleBtn.css('background', '#dc2626');
        }

        // Auto-select the first non-identity operation
        let firstNonIdentity = this.cartesianOps.findIndex(op => op.label !== 'E (Identity)');
        if (firstNonIdentity < 0) firstNonIdentity = 0;
        if (this.dropdownEl) this.dropdownEl.val(firstNonIdentity);
        this.selectOperation(firstNonIdentity);
    }

    deactivate() {
        this.active = false;
        this.currentOpIndex = -1;
        this.sliderValue = 0;

        // Restore phonon vibrations
        this.crystal.symmetryAnimationActive = false;
        this.crystal.amplitude = this.savedAmplitude;
        this.crystal.paused = this.savedPaused;

        // Remove ghost lattice
        this.removeGhostLattice();

        // Reset atom positions to equilibrium
        this.resetAtomPositions();

        // Hide UI panel
        if (this.panelEl) this.panelEl.hide();
        if (this.toggleBtn) {
            this.toggleBtn.text('⚛ Symmetry');
            this.toggleBtn.css('background', 'rgba(2, 132, 199, 0.9)');
        }

        // Trigger a render
        this.crystal.needsRender = true;
        this.crystal.startAnimationLoop();
    }

    // ─────────────────────────────────────────────
    // MATH ENGINE
    // ─────────────────────────────────────────────

    /**
     * Convert a fractional 3×3 rotation matrix R_frac to Cartesian space.
     * R_cart = L · R_frac · L^{-1}
     * where L is the lattice matrix (rows are lattice vectors a, b, c).
     */
    fractionalToCartesianRotation(R_frac, lat) {
        let L = lat;             // 3×3 lattice vectors
        let L_inv = mat.matrix_inverse(L);
        if (!L_inv) return null;
        let temp = mat.matrix_multiply(L, R_frac);
        return mat.matrix_multiply(temp, L_inv);
    }

    /**
     * Convert a fractional translation vector to Cartesian.
     * t_cart = L^T · t_frac  (since positions in Cartesian = frac · L)
     * Actually: r_cart = frac_coords · lat_matrix for row-vector convention.
     * So t_cart[i] = sum_j t_frac[j] * lat[j][i]
     */
    fractionalToCartesianTranslation(t_frac, lat) {
        let t_cart = [0, 0, 0];
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
                t_cart[i] += t_frac[j] * lat[j][i];
            }
        }
        return t_cart;
    }

    /**
     * Gram-Schmidt orthonormalization for a 3×3 matrix.
     * Guarantees R is a proper rotation (det=+1 orthogonal) or
     * improper rotation (det=-1) suitable for Three.js.
     */
    orthonormalize(M) {
        // Extract column vectors
        let u0 = [M[0][0], M[1][0], M[2][0]];
        let u1 = [M[0][1], M[1][1], M[2][1]];
        let u2 = [M[0][2], M[1][2], M[2][2]];

        // Normalize u0
        let n0 = Math.sqrt(u0[0]*u0[0] + u0[1]*u0[1] + u0[2]*u0[2]);
        if (n0 < 1e-12) return M; // degenerate
        u0 = [u0[0]/n0, u0[1]/n0, u0[2]/n0];

        // u1 = u1 - proj(u1, u0)
        let d10 = u1[0]*u0[0] + u1[1]*u0[1] + u1[2]*u0[2];
        u1 = [u1[0] - d10*u0[0], u1[1] - d10*u0[1], u1[2] - d10*u0[2]];
        let n1 = Math.sqrt(u1[0]*u1[0] + u1[1]*u1[1] + u1[2]*u1[2]);
        if (n1 < 1e-12) return M;
        u1 = [u1[0]/n1, u1[1]/n1, u1[2]/n1];

        // u2 = u0 × u1 (guarantees right-handed orthonormal frame)
        u2 = [
            u0[1]*u1[2] - u0[2]*u1[1],
            u0[2]*u1[0] - u0[0]*u1[2],
            u0[0]*u1[1] - u0[1]*u1[0]
        ];

        // Check if original determinant was negative (improper rotation)
        let detOrig = mat.matrix_determinant(M);
        if (detOrig < 0) {
            u2 = [-u2[0], -u2[1], -u2[2]];
        }

        // Reconstruct matrix from column vectors
        return [
            [u0[0], u1[0], u2[0]],
            [u0[1], u1[1], u2[1]],
            [u0[2], u1[2], u2[2]]
        ];
    }

    /**
     * Extract rotation axis and angle from a 3×3 proper rotation matrix.
     * For improper rotations (det=-1), extracts from -R.
     * Returns { axis: THREE.Vector3, angle: Number (radians), isImproper: Boolean }
     */
    extractAxisAngle(R_cart) {
        let detR = mat.matrix_determinant(R_cart);
        let isImproper = (detR < 0);

        // Work with proper part
        let M = isImproper ? mat.matrix_scale(R_cart, -1) : R_cart;

        let trace = M[0][0] + M[1][1] + M[2][2];
        let cosAngle = (trace - 1.0) / 2.0;
        cosAngle = Math.max(-1, Math.min(1, cosAngle)); // clamp

        let angle = Math.acos(cosAngle);

        let axis;
        if (Math.abs(angle) < 1e-8) {
            // Identity: no rotation. Pick arbitrary axis.
            axis = new THREE.Vector3(0, 0, 1);
            angle = 0;
        } else if (Math.abs(angle - Math.PI) < 1e-8) {
            // 180° rotation. Find axis from eigenvector of eigenvalue +1.
            // (M + I) has rank ≤ 2; find the largest diagonal and solve.
            let Mpi = [
                [M[0][0] + 1, M[0][1], M[0][2]],
                [M[1][0], M[1][1] + 1, M[1][2]],
                [M[2][0], M[2][1], M[2][2] + 1]
            ];
            // Use the row with the largest norm
            let best = 0;
            let bestNorm = 0;
            for (let i = 0; i < 3; i++) {
                let n = Mpi[i][0]*Mpi[i][0] + Mpi[i][1]*Mpi[i][1] + Mpi[i][2]*Mpi[i][2];
                if (n > bestNorm) { bestNorm = n; best = i; }
            }
            let v = Mpi[best];
            let vn = Math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
            axis = new THREE.Vector3(v[0]/vn, v[1]/vn, v[2]/vn);
        } else {
            // General case: axis from skew-symmetric part
            let ax = M[2][1] - M[1][2];
            let ay = M[0][2] - M[2][0];
            let az = M[1][0] - M[0][1];
            let an = Math.sqrt(ax*ax + ay*ay + az*az);
            if (an < 1e-12) {
                axis = new THREE.Vector3(0, 0, 1);
            } else {
                axis = new THREE.Vector3(ax/an, ay/an, az/an);
            }
        }

        return { axis, angle, isImproper };
    }

    /**
     * Classify a symmetry operation and return a human-readable label.
     */
    classifyOperation(R_frac, t_frac) {
        let detR = Math.round(
            R_frac[0][0]*(R_frac[1][1]*R_frac[2][2] - R_frac[2][1]*R_frac[1][2]) -
            R_frac[0][1]*(R_frac[1][0]*R_frac[2][2] - R_frac[1][2]*R_frac[2][0]) +
            R_frac[0][2]*(R_frac[1][0]*R_frac[2][1] - R_frac[1][1]*R_frac[2][0])
        );
        let tr = Math.round(R_frac[0][0] + R_frac[1][1] + R_frac[2][2]);

        // Check if translation is non-zero (screw/glide)
        let hasTranslation = (Math.abs(t_frac[0]) > 1e-6 ||
                              Math.abs(t_frac[1]) > 1e-6 ||
                              Math.abs(t_frac[2]) > 1e-6);
        let tLabel = hasTranslation ? ' + τ' : '';

        if (detR === 1) {
            if (tr === 3)  return 'E (Identity)';
            if (tr === -1) return 'C₂ (180° rotation)' + tLabel;
            if (tr === 0)  return 'C₃ (120° rotation)' + tLabel;
            if (tr === 1)  return 'C₄ (90° rotation)' + tLabel;
            if (tr === 2)  return 'C₆ (60° rotation)' + tLabel;
        } else if (detR === -1) {
            if (tr === -3) return 'i (Inversion)' + tLabel;
            if (tr === 1)  return 'σ (Mirror)' + tLabel;
            if (tr === -2) return 'S₃ (Rotoreflection 120°)' + tLabel;
            if (tr === -1) return 'S₄ (Rotoreflection 90°)' + tLabel;
            if (tr === 0)  return 'S₆ (Rotoreflection 60°)' + tLabel;
        }
        return 'Unknown operation';
    }

    // ─────────────────────────────────────────────
    // PRECOMPUTATION
    // ─────────────────────────────────────────────

    precomputeOperations() {
        let rots = this.phonon.crystal_symmetries.rotations;
        let trans = this.phonon.crystal_symmetries.translations;
        let lat = this.phonon.lat;

        this.cartesianOps = [];

        for (let i = 0; i < rots.length; i++) {
            let R_frac = rots[i];
            let t_frac = trans[i];

            // Convert to Cartesian
            let R_cart = this.fractionalToCartesianRotation(R_frac, lat);
            if (!R_cart) continue;

            // Orthonormalize to prevent numerical drift
            R_cart = this.orthonormalize(R_cart);

            let t_cart = this.fractionalToCartesianTranslation(t_frac, lat);

            // Extract axis and angle
            let { axis, angle, isImproper } = this.extractAxisAngle(R_cart);

            // Human-readable label
            let label = this.classifyOperation(R_frac, t_frac);

            this.cartesianOps.push({
                R_frac,
                t_frac,
                R_cart,
                t_cart,
                axis,
                angle,
                isImproper,
                label,
                index: i
            });
        }
    }

    // ─────────────────────────────────────────────
    // DROPDOWN POPULATION
    // ─────────────────────────────────────────────

    populateDropdown() {
        if (!this.dropdownEl) return;
        this.dropdownEl.empty();

        // Group operations by type for a cleaner UI
        for (let i = 0; i < this.cartesianOps.length; i++) {
            let op = this.cartesianOps[i];
            let angleDeg = Math.round(op.angle * 180 / Math.PI);
            let displayLabel = `#${i}: ${op.label}`;
            if (op.label !== 'E (Identity)' && !op.isImproper && angleDeg > 0) {
                displayLabel += ` [${angleDeg}°]`;
            }
            this.dropdownEl.append(`<option value="${i}">${displayLabel}</option>`);
        }
    }

    // ─────────────────────────────────────────────
    // OPERATION SELECTION & ANIMATION
    // ─────────────────────────────────────────────

    selectOperation(idx) {
        if (idx < 0 || idx >= this.cartesianOps.length) return;
        this.currentOpIndex = idx;

        // Reset slider
        this.sliderValue = 0;
        if (this.sliderEl) this.sliderEl.val(0);

        // Update label
        let op = this.cartesianOps[idx];
        if (this.labelEl) {
            let angleDeg = Math.round(op.angle * 180 / Math.PI);
            let info = op.label;
            if (!op.isImproper && op.angle > 0.01) {
                info += ` — Axis: (${op.axis.x.toFixed(2)}, ${op.axis.y.toFixed(2)}, ${op.axis.z.toFixed(2)}), Angle: ${angleDeg}°`;
            }
            this.labelEl.text(info);
        }

        // Create ghost lattice
        this.removeGhostLattice();
        this.createGhostLattice();

        // Set atoms to equilibrium
        this.resetAtomPositions();

        this.crystal.needsRender = true;
        this.crystal.startAnimationLoop();
    }

    setSliderValue(t) {
        this.sliderValue = Math.max(0, Math.min(1, t));
        if (this.currentOpIndex < 0) return;

        let op = this.cartesianOps[this.currentOpIndex];
        this.applyInterpolatedOperation(op, this.sliderValue);

        this.crystal.needsRender = true;
        this.crystal.startAnimationLoop();
    }

    /**
     * Interpolate atom positions from equilibrium (t=0) to
     * the symmetry-operated position (t=1).
     *
     * For proper rotations: use quaternion SLERP for smooth curved paths.
     * For improper rotations (inversion, mirrors): use linear interpolation.
     */
    applyInterpolatedOperation(op, t) {
        if (!this.crystal.atomobjects || !this.crystal.atompos) return;

        let nAtoms = this.crystal.atomobjects.length;
        let center = this.crystal.geometricCenter;

        // Build Three.js matrix for the full operation
        let R_mat4 = new THREE.Matrix4();
        let quat_full = new THREE.Quaternion();

        if (!op.isImproper) {
            // Proper rotation: SLERP from identity to target quaternion
            let quat_target = new THREE.Quaternion();
            quat_target.setFromAxisAngle(op.axis, op.angle);

            // SLERP: identity → target
            let quat_identity = new THREE.Quaternion(); // identity
            quat_full.slerpQuaternions(quat_identity, quat_target, t);
        }

        for (let i = 0; i < nAtoms; i++) {
            let eqPos = this.crystal.atompos[i]; // THREE.Vector3 (centered)

            let newPos = new THREE.Vector3();

            if (op.isImproper) {
                // For inversion/mirror: linear interpolation from r to R·r + t
                // Compute target position: R_cart · r + t_cart (centered)
                let rx = eqPos.x, ry = eqPos.y, rz = eqPos.z;
                let R = op.R_cart;
                let tx = op.t_cart[0], ty = op.t_cart[1], tz = op.t_cart[2];

                let targetX = R[0][0]*rx + R[0][1]*ry + R[0][2]*rz + tx;
                let targetY = R[1][0]*rx + R[1][1]*ry + R[1][2]*rz + ty;
                let targetZ = R[2][0]*rx + R[2][1]*ry + R[2][2]*rz + tz;

                // Linear interpolation
                newPos.set(
                    rx + t * (targetX - rx),
                    ry + t * (targetY - ry),
                    rz + t * (targetZ - rz)
                );
            } else {
                // Proper rotation: apply interpolated quaternion
                newPos.copy(eqPos);

                // Also interpolate translation
                let tVec = new THREE.Vector3(
                    t * op.t_cart[0],
                    t * op.t_cart[1],
                    t * op.t_cart[2]
                );

                newPos.applyQuaternion(quat_full);
                newPos.add(tVec);
            }

            // Update atom position
            this.crystal.atomobjects[i].position.copy(newPos);

            // Update InstancedMesh matrix
            let atomInstance = this.crystal.atomInstanceRefs[i];
            if (atomInstance) {
                this.crystal.instanceDummy.position.copy(newPos);
                this.crystal.instanceDummy.quaternion.set(0, 0, 0, 1);
                this.crystal.instanceDummy.scale.set(1, 1, 1);
                this.crystal.instanceDummy.updateMatrix();
                atomInstance.mesh.setMatrixAt(atomInstance.instanceId, this.crystal.instanceDummy.matrix);
            }
        }

        // Mark all instanced meshes as needing update
        for (let i = 0; i < this.crystal.atommeshes.length; i++) {
            this.crystal.atommeshes[i].instanceMatrix.needsUpdate = true;
        }

        // Update bonds
        this.updateBondsForCurrentPositions();
    }

    /**
     * Reset all atom positions back to their equilibrium positions.
     */
    resetAtomPositions() {
        if (!this.crystal.atomobjects || !this.crystal.atompos) return;

        for (let i = 0; i < this.crystal.atomobjects.length; i++) {
            let eqPos = this.crystal.atompos[i];
            this.crystal.atomobjects[i].position.copy(eqPos);

            let atomInstance = this.crystal.atomInstanceRefs[i];
            if (atomInstance) {
                this.crystal.instanceDummy.position.copy(eqPos);
                this.crystal.instanceDummy.quaternion.set(0, 0, 0, 1);
                this.crystal.instanceDummy.scale.set(1, 1, 1);
                this.crystal.instanceDummy.updateMatrix();
                atomInstance.mesh.setMatrixAt(atomInstance.instanceId, this.crystal.instanceDummy.matrix);
            }
        }

        for (let i = 0; i < this.crystal.atommeshes.length; i++) {
            this.crystal.atommeshes[i].instanceMatrix.needsUpdate = true;
        }

        this.updateBondsForCurrentPositions();
    }

    /**
     * Recalculate bond positions based on current atom positions.
     * Uses the same bond-update logic as VibCrystal.render().
     */
    updateBondsForCurrentPositions() {
        if (!this.crystal.bonds) return;

        for (let i = 0; i < this.crystal.bonds.length; i++) {
            let bond = this.crystal.bonds[i];
            let a = bond.a; // THREE.Vector3 reference (points to atom.position)
            let b = bond.b;

            // Compute midpoint and direction
            let midpoint = new THREE.Vector3().addVectors(a, b).multiplyScalar(0.5);
            let direction = new THREE.Vector3().subVectors(b, a);
            let lengthNow = direction.length();
            direction.normalize();

            // Quaternion to orient the cylinder
            let quaternion = new THREE.Quaternion();
            let yAxis = new THREE.Vector3(0, 1, 0);
            quaternion.setFromUnitVectors(yAxis, direction);

            if (this.crystal.bondColorByAtom && this.crystal.splitBondObjects && this.crystal.splitBondObjects.length) {
                let offset = direction.clone().multiplyScalar(lengthNow * 0.25);
                let bondPair = this.crystal.splitBondObjects[i];
                if (bondPair && bondPair.meshA && bondPair.meshB) {
                    bondPair.meshA.position.copy(midpoint).sub(offset);
                    bondPair.meshB.position.copy(midpoint).add(offset);
                    bondPair.meshA.setRotationFromQuaternion(quaternion);
                    bondPair.meshB.setRotationFromQuaternion(quaternion);
                    bondPair.meshA.scale.set(1, lengthNow * 0.5, 1);
                    bondPair.meshB.scale.set(1, lengthNow * 0.5, 1);
                }
            } else if (this.crystal.bondColorByAtom && this.crystal.bondmeshes && this.crystal.bondmeshes.length >= 2) {
                let offset = direction.clone().multiplyScalar(lengthNow * 0.25);
                this.crystal.instanceDummy.quaternion.copy(quaternion);
                this.crystal.instanceDummy.scale.set(1, lengthNow * 0.5, 1);

                this.crystal.instanceDummy.position.copy(midpoint).sub(offset);
                this.crystal.instanceDummy.updateMatrix();
                this.crystal.bondmeshes[0].setMatrixAt(i, this.crystal.instanceDummy.matrix);

                this.crystal.instanceDummy.position.copy(midpoint).add(offset);
                this.crystal.instanceDummy.updateMatrix();
                this.crystal.bondmeshes[1].setMatrixAt(i, this.crystal.instanceDummy.matrix);
            } else if (this.crystal.bondmesh) {
                this.crystal.instanceDummy.position.copy(midpoint);
                this.crystal.instanceDummy.quaternion.copy(quaternion);
                this.crystal.instanceDummy.scale.set(1, lengthNow, 1);
                this.crystal.instanceDummy.updateMatrix();
                this.crystal.bondmesh.setMatrixAt(i, this.crystal.instanceDummy.matrix);
            }
        }

        // Mark bond meshes as needing update
        if (this.crystal.bondmeshes && this.crystal.bondmeshes.length) {
            for (let i = 0; i < this.crystal.bondmeshes.length; i++) {
                this.crystal.bondmeshes[i].instanceMatrix.needsUpdate = true;
            }
        } else if (this.crystal.bondmesh) {
            this.crystal.bondmesh.instanceMatrix.needsUpdate = true;
        }
    }

    // ─────────────────────────────────────────────
    // GHOST LATTICE
    // ─────────────────────────────────────────────

    /**
     * Create semi-transparent "ghost" copies of all atoms at their
     * equilibrium positions as a visual reference.
     */
    createGhostLattice() {
        if (!this.crystal.atomobjects || !this.crystal.atompos) return;
        this.removeGhostLattice();

        let ghostMaterial = new THREE.MeshLambertMaterial({
            color: 0x888888,
            transparent: true,
            opacity: 0.25,
            depthWrite: false
        });

        let sphereGeom = new THREE.SphereGeometry(0.3, 16, 12);

        for (let i = 0; i < this.crystal.atompos.length; i++) {
            let pos = this.crystal.atompos[i];
            let mesh = new THREE.Mesh(sphereGeom, ghostMaterial);
            mesh.position.copy(pos);
            mesh.name = 'symmetry-ghost';
            this.crystal.scene.add(mesh);
            this.ghostMeshes.push(mesh);
        }

        // Ghost bonds
        if (this.crystal.bonds && this.crystal.bonds.length > 0) {
            let bondMat = new THREE.MeshLambertMaterial({
                color: 0x666666,
                transparent: true,
                opacity: 0.15,
                depthWrite: false
            });

            let bondGeom = new THREE.CylinderGeometry(0.06, 0.06, 1, 6);

            for (let i = 0; i < this.crystal.bonds.length; i++) {
                let bond = this.crystal.bonds[i];
                let a = this.crystal.atompos[bond.atomIndexA !== undefined ? bond.atomIndexA : 0];
                let b = this.crystal.atompos[bond.atomIndexB !== undefined ? bond.atomIndexB : 0];

                // Use stored atom positions for bond endpoints
                if (!a || !b) continue;

                let midpoint = new THREE.Vector3().addVectors(a, b).multiplyScalar(0.5);
                let dir = new THREE.Vector3().subVectors(b, a);
                let len = dir.length();
                dir.normalize();

                let bondMesh = new THREE.Mesh(bondGeom, bondMat);
                bondMesh.position.copy(midpoint);
                bondMesh.scale.set(1, len, 1);
                let yAxis = new THREE.Vector3(0, 1, 0);
                bondMesh.quaternion.setFromUnitVectors(yAxis, dir);
                bondMesh.name = 'symmetry-ghost';
                this.crystal.scene.add(bondMesh);
                this.ghostMeshes.push(bondMesh);
            }
        }
    }

    removeGhostLattice() {
        for (let mesh of this.ghostMeshes) {
            this.crystal.scene.remove(mesh);
            if (mesh.geometry) mesh.geometry.dispose();
            if (mesh.material) mesh.material.dispose();
        }
        this.ghostMeshes = [];
    }

    // ─────────────────────────────────────────────
    // CLEANUP (called when material changes)
    // ─────────────────────────────────────────────

    onMaterialChanged() {
        if (this.active) {
            this.deactivate();
        }
    }
}
