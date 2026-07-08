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
import { createCellLineObject } from './viewergeometry.js';

export class SymmetryVisualizer {

    constructor(vibcrystal) {
        /** @type {VibCrystal} */
        this.crystal = vibcrystal;
        this.phonon = null;

        // State
        this.active = false;
        this.currentOpIndex = -1;
        this.sliderRotValue = 0.0;
        this.sliderTransValue = 0.0;

        // Ghost lattice meshes (THREE objects added to scene)
        this.ghostMeshes = [];

        // Saved state to restore on deactivation
        this.savedAmplitude = 0;
        this.savedPaused = false;

        // Precomputed Cartesian operations
        this.cartesianOps = []; // [{ R_cart, t_cart, axis, angle, det, label }]

        // DOM references (set externally)
        this.dropdownEl = null;
        this.sliderRotEl = null;
        this.sliderTransEl = null;
        this.sliderRotContainer = null;
        this.sliderTransContainer = null;
        this.sliderRotLabel = null;
        this.panelEl = null;
        this.labelEl = null;
        this.toggleBtn = null;
        this.ghostAtomsCheckboxEl = null;
        this.bondsCheckboxEl = null;

        // Hook into structure updates (e.g. changing cell repetitions)
        this.crystal.onStructureRebuilt = () => this.refreshGhostLattice();
    }

    // ─────────────────────────────────────────────
    // DOM BINDING
    // ─────────────────────────────────────────────

    bindDOM(panelEl, dropdownEl, sliderRotEl, sliderTransEl, rotContainer, transContainer, rotLabel, labelEl, toggleBtn, ghostAtomsCheckboxEl, bondsCheckboxEl, planesCheckboxEl, autoplayCheckboxEl) {
        this.panelEl = panelEl;
        this.dropdownEl = dropdownEl;
        this.sliderRotEl = sliderRotEl;
        this.sliderTransEl = sliderTransEl;
        this.sliderRotContainer = rotContainer;
        this.sliderTransContainer = transContainer;
        this.sliderRotLabel = rotLabel;
        this.labelEl = labelEl;
        this.toggleBtn = toggleBtn;
        this.ghostAtomsCheckboxEl = ghostAtomsCheckboxEl;
        this.bondsCheckboxEl = bondsCheckboxEl;
        this.planesCheckboxEl = planesCheckboxEl;
        this.autoplayCheckboxEl = autoplayCheckboxEl;

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

        // Sliders: animate the operation
        if (this.sliderRotEl) {
            this.sliderRotEl.on('input', () => {
                let t = parseFloat(this.sliderRotEl.val()) / 100.0;
                this.setSliderValues(t, this.sliderTransValue);
            });
        }
        if (this.sliderTransEl) {
            this.sliderTransEl.on('input', () => {
                let t = parseFloat(this.sliderTransEl.val()) / 100.0;
                this.setSliderValues(this.sliderRotValue, t);
            });
        }
        
        // Checkbox: toggle ghost atoms
        if (this.ghostAtomsCheckboxEl) {
            this.ghostAtomsCheckboxEl.on('change', () => {
                let isChecked = this.ghostAtomsCheckboxEl.is(':checked');
                if (!isChecked && this.bondsCheckboxEl) {
                    this.bondsCheckboxEl.prop('checked', false);
                    this.bondsCheckboxEl.prop('disabled', true);
                } else if (this.bondsCheckboxEl) {
                    this.bondsCheckboxEl.prop('disabled', false);
                }
                this.refreshGhostLattice();
            });
        }

        // Checkbox: toggle ghost bonds
        if (this.bondsCheckboxEl) {
            this.bondsCheckboxEl.on('change', () => {
                this.refreshGhostBondsVisibility();
            });
        }

        // Checkbox: toggle planes and axes
        if (this.planesCheckboxEl) {
            this.planesCheckboxEl.on('change', () => {
                let isChecked = this.planesCheckboxEl.is(':checked');
                if (this.axesAndPlanesGroup) {
                    this.axesAndPlanesGroup.visible = isChecked;
                    this.crystal.needsRender = true;
                    this.crystal.startAnimationLoop();
                }
            });
        }

        // Checkbox: Autoplay
        if (this.autoplayCheckboxEl) {
            this.autoplayCheckboxEl.on('change', () => {
                if (this.autoplayCheckboxEl.is(':checked')) {
                    this.startAutoplay();
                } else {
                    this.stopAutoplay();
                }
            });
        }
    }

    bindReferenceDOM(refCrystal, showBtn, popupEl, closeBtn, syncBtn, headerEl) {
        this.refCrystal = refCrystal;
        this.refPopupEl = popupEl;

        if (showBtn) {
            showBtn.on('click', () => {
                if (!this.refPopupEl) return;
                
                if (this.refPopupEl.is(':visible')) {
                    this.refPopupEl.hide();
                    showBtn.text('Show Reference Crystal');
                } else {
                    this.refPopupEl.show();
                    showBtn.text('Hide Reference Crystal');
                    // trigger resize for the ref crystal to ensure canvas sizes correctly
                    if (this.refCrystal && this.refCrystal.onWindowResize) {
                        setTimeout(() => this.refCrystal.onWindowResize(), 10);
                    }
                }
            });
        }

        if (closeBtn) {
            closeBtn.on('click', () => {
                if (this.refPopupEl) {
                    this.refPopupEl.hide();
                    if (showBtn) showBtn.text('Show Reference Crystal');
                }
            });
        }

        if (syncBtn) {
            syncBtn.on('click', () => {
                if (this.crystal && this.refCrystal) {
                    this.refCrystal.camera.position.copy(this.crystal.camera.position);
                    this.refCrystal.camera.quaternion.copy(this.crystal.camera.quaternion);
                    if (this.crystal.controls && this.refCrystal.controls) {
                        this.refCrystal.controls.target.copy(this.crystal.controls.target);
                    }
                    this.refCrystal.needsRender = true;
                    this.refCrystal.startAnimationLoop();
                }
            });
        }

        // Custom drag logic for the popup window
        if (headerEl && popupEl) {
            let isDragging = false;
            let startX, startY, initialTop, initialLeft;

            headerEl.on('mousedown', (e) => {
                isDragging = true;
                startX = e.clientX;
                startY = e.clientY;
                let pos = popupEl.position();
                initialTop = pos.top;
                initialLeft = pos.left;
                // prevent selection
                e.preventDefault();
            });

            $(document).on('mousemove', (e) => {
                if (isDragging) {
                    let dx = e.clientX - startX;
                    let dy = e.clientY - startY;
                    popupEl.css({
                        top: initialTop + dy + 'px',
                        left: initialLeft + dx + 'px',
                        right: 'auto' // clear 'right' to respect 'left'
                    });
                }
            });

            $(document).on('mouseup', () => {
                isDragging = false;
            });
        }

        // Add ResizeObserver to update the Three.js canvas when the popup is resized
        if (popupEl && popupEl.length > 0) {
            try {
                const resizeObserver = new ResizeObserver(() => {
                    if (this.refCrystal && this.refCrystal.onWindowResize && popupEl.is(':visible')) {
                        this.refCrystal.onWindowResize();
                    }
                });
                resizeObserver.observe(popupEl[0]);
            } catch (e) {
                console.warn('ResizeObserver not supported in this browser.', e);
            }
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

        // Freeze phonon vibrations — snap the current phase so arrows show
        // the eigenvector direction at the exact moment symmetry was opened.
        this.savedAmplitude = this.crystal.amplitude;
        this.savedPaused = this.crystal.paused;
        this.crystal.symmetryPhaseSnap = this.crystal.time; // capture current phase
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

        // Add axes and planes
        this.createAxesAndPlanes();

        if (this.autoplayCheckboxEl && this.autoplayCheckboxEl.is(':checked')) {
            this.startAutoplay();
        }
    }

    deactivate() {
        this.active = false;
        this.currentOpIndex = -1;
        this.sliderRotValue = 0;
        this.sliderTransValue = 0;

        this.stopAutoplay();

        // Restore phonon vibrations
        this.crystal.symmetryAnimationActive = false;
        this.crystal.amplitude = this.savedAmplitude;
        this.crystal.paused = this.savedPaused;

        if (this.symElementMesh) {
            this.crystal.scene.remove(this.symElementMesh);
            if (this.symElementMesh.geometry) this.symElementMesh.geometry.dispose();
            if (this.symElementMesh.material) this.symElementMesh.material.dispose();
            this.symElementMesh = null;
        }

        // Remove ghost lattice
        this.removeGhostLattice();

        // Remove axes and planes
        this.removeAxesAndPlanes();

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
    // AXES AND PLANES (Feature 1)
    // ─────────────────────────────────────────────

    createAxesAndPlanes() {
        this.removeAxesAndPlanes();

        let maxDist = 5;
        if (this.crystal.atompos && this.crystal.atompos.length > 0) {
            for (let pos of this.crystal.atompos) {
                maxDist = Math.max(maxDist, pos.length());
            }
        }
        let size = maxDist * 2.5;

        // The physical crystallographic origin is at -geometricCenter in viewer space
        let origin = new THREE.Vector3().copy(this.crystal.geometricCenter).multiplyScalar(-1);

        this.axesAndPlanesGroup = new THREE.Group();
        this.axesAndPlanesGroup.position.copy(origin);

        // Axes (Red=X, Green=Y, Blue=Z)
        let axesHelper = new THREE.AxesHelper(size * 0.8);
        this.axesAndPlanesGroup.add(axesHelper);

        // Planes (XY, YZ, ZX)
        let createPlane = (color, rotX, rotY, rotZ) => {
            // Add grid helper on the plane
            let grid = new THREE.GridHelper(size, 50, color, color);
            grid.material.opacity = 0.12;
            grid.material.transparent = true;
            grid.rotation.set(rotX, rotY, rotZ);

            // GridHelper in Three.js is created on the XZ plane by default.
            // PlaneGeometry is created on the XY plane.
            // So we must rotate the grid by 90 degrees on X to match the plane coordinates.
            grid.rotateX(Math.PI / 2);

            this.axesAndPlanesGroup.add(grid);
        };

        // XY plane (blueish) - normal is Z
        createPlane(0x0000ff, 0, 0, 0);
        // YZ plane (reddish) - normal is X
        createPlane(0xff0000, 0, Math.PI/2, 0);
        // ZX plane (greenish) - normal is Y
        createPlane(0x00ff00, Math.PI/2, 0, 0);

        if (this.planesCheckboxEl && !this.planesCheckboxEl.is(':checked')) {
            this.axesAndPlanesGroup.visible = false;
        }

        this.crystal.scene.add(this.axesAndPlanesGroup);
    }

    removeAxesAndPlanes() {
        if (this.axesAndPlanesGroup) {
            this.crystal.scene.remove(this.axesAndPlanesGroup);
            this.axesAndPlanesGroup.traverse((child) => {
                if (child.geometry) child.geometry.dispose();
                if (child.material) child.material.dispose();
            });
            this.axesAndPlanesGroup = null;
        }
    }

    createReferenceAxesAndPlanes() {
        if (!this.refCrystal) return;

        // Remove old if any
        if (this.refAxesAndPlanesGroup) {
            this.refCrystal.scene.remove(this.refAxesAndPlanesGroup);
            this.refAxesAndPlanesGroup.traverse((child) => {
                if (child.geometry) child.geometry.dispose();
                if (child.material) child.material.dispose();
            });
        }

        let maxDist = 5;
        if (this.refCrystal.atompos && this.refCrystal.atompos.length > 0) {
            for (let pos of this.refCrystal.atompos) {
                maxDist = Math.max(maxDist, pos.length());
            }
        }
        let size = maxDist * 2.5;

        // Physical crystallographic origin
        let origin = new THREE.Vector3().copy(this.refCrystal.geometricCenter).multiplyScalar(-1);

        this.refAxesAndPlanesGroup = new THREE.Group();
        this.refAxesAndPlanesGroup.position.copy(origin);

        // Axes (Red=X, Green=Y, Blue=Z)
        let axesHelper = new THREE.AxesHelper(size * 0.8);
        this.refAxesAndPlanesGroup.add(axesHelper);

        let createPlane = (color, rotX, rotY, rotZ) => {
            let grid = new THREE.GridHelper(size, 50, color, color);
            grid.material.opacity = 0.12;
            grid.material.transparent = true;
            grid.rotation.set(rotX, rotY, rotZ);
            grid.rotateX(Math.PI / 2);
            this.refAxesAndPlanesGroup.add(grid);
        };

        // XY plane (blueish) - normal is Z
        createPlane(0x0000ff, 0, 0, 0);
        // YZ plane (reddish) - normal is X
        createPlane(0xff0000, 0, Math.PI/2, 0);
        // ZX plane (greenish) - normal is Y
        createPlane(0x00ff00, Math.PI/2, 0, 0);

        this.refCrystal.scene.add(this.refAxesAndPlanesGroup);
        this.refCrystal.needsRender = true;
        this.refCrystal.startAnimationLoop();
    }

    // ─────────────────────────────────────────────
    // AUTOPLAY ENGINE
    // ─────────────────────────────────────────────

    startAutoplay() {
        if (this.isAutoplaying) return;
        this.isAutoplaying = true;
        this.lastTime = performance.now();
        this.autoplayPhase = 0;
        this.autoplayLoop();
    }

    stopAutoplay() {
        this.isAutoplaying = false;
        if (this.autoplayFrameId) {
            cancelAnimationFrame(this.autoplayFrameId);
            this.autoplayFrameId = null;
        }
    }

    autoplayLoop() {
        if (!this.isAutoplaying || !this.active) return;

        let now = performance.now();
        let dt = (now - this.lastTime) / 1000.0;
        this.lastTime = now;

        // Speed: 0.5 per second (2s per animation phase)
        let speed = 0.5;
        this.autoplayPhase += speed * dt;

        let hasTrans = this.sliderTransContainer && this.sliderTransContainer.is(':visible');
        let maxPhase = 1.0; // Both animate simultaneously

        // Add a 1 second pause (0.5 phase at speed 0.5) at the end of the animation
        if (this.autoplayPhase > maxPhase + 0.5) {
            this.autoplayPhase = 0;
        }

        let rotVal = this.autoplayPhase;
        let transVal = hasTrans ? this.autoplayPhase : 0;

        rotVal = Math.min(Math.max(rotVal, 0), 1);
        transVal = Math.min(Math.max(transVal, 0), 1);

        if (this.sliderRotEl) this.sliderRotEl.val(Math.round(rotVal * 100));
        if (this.sliderTransEl) this.sliderTransEl.val(Math.round(transVal * 100));
        
        this.setSliderValues(rotVal, transVal);

        this.autoplayFrameId = requestAnimationFrame(() => this.autoplayLoop());
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
     * Convert a Cartesian translation vector to fractional.
     */
    cartesianToFractionalTranslation(t_cart, lat) {
        let L_inv = mat.matrix_inverse(lat);
        if (!L_inv) return [0, 0, 0];
        let t_frac = [0, 0, 0];
        // Fractional = Cartesian * L^-1
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
                t_frac[i] += t_cart[j] * L_inv[j][i];
            }
        }
        return t_frac;
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

        // Calculate the exact mathematical center of the drawn supercell
        let centerCart = new THREE.Vector3(0, 0, 0);
        if (this.crystal.atoms && this.crystal.atoms.length > 0) {
            let box = new THREE.Box3();
            for (let atom of this.crystal.atoms) {
                box.expandByPoint(new THREE.Vector3(atom[1], atom[2], atom[3]));
            }
            box.getCenter(centerCart);
        }
        let centerFrac = this.cartesianToFractionalTranslation([centerCart.x, centerCart.y, centerCart.z], lat);

        this.cartesianOps = [];

        for (let i = 0; i < rots.length; i++) {
            let R_frac = rots[i];
            let t_frac = trans[i];

            // Auto-shift the operation pivot to the center of the drawn supercell
            // T_frac = round(center_frac - R_frac * center_frac - t_frac)
            let R_center = [0, 0, 0];
            for (let r = 0; r < 3; r++) {
                for (let c = 0; c < 3; c++) {
                    R_center[r] += R_frac[r][c] * centerFrac[c];
                }
            }
            
            let T_frac = [
                Math.round(centerFrac[0] - R_center[0] - t_frac[0]),
                Math.round(centerFrac[1] - R_center[1] - t_frac[1]),
                Math.round(centerFrac[2] - R_center[2] - t_frac[2])
            ];
            
            let t_frac_shifted = [
                t_frac[0] + T_frac[0],
                t_frac[1] + T_frac[1],
                t_frac[2] + T_frac[2]
            ];

            // Convert to Cartesian
            let R_cart = this.fractionalToCartesianRotation(R_frac, lat);
            if (!R_cart) continue;

            // Orthonormalize to prevent numerical drift
            R_cart = this.orthonormalize(R_cart);

            let t_cart = this.fractionalToCartesianTranslation(t_frac_shifted, lat);

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

        let q_cart = null;
        let lat = null;
        if (this.crystal && this.crystal.phononweb && this.crystal.phononweb.phonon) {
            let k_idx = this.crystal.phononweb.k;
            if (this.crystal.phononweb.phonon.kpoints && this.crystal.phononweb.phonon.kpoints[k_idx]) {
                let q = this.crystal.phononweb.phonon.kpoints[k_idx];
                lat = this.crystal.phononweb.phonon.lat;
                
                // Calculate reciprocal lattice vectors b1, b2, b3 without 2pi factor
                // b_i dot a_j = delta_ij
                let a1 = lat[0], a2 = lat[1], a3 = lat[2];
                let b1 = mat.vec_cross(a2, a3);
                let b2 = mat.vec_cross(a3, a1);
                let b3 = mat.vec_cross(a1, a2);
                let v = mat.vec_dot(a1, b1);
                b1 = mat.vec_scale(b1, 1/v);
                b2 = mat.vec_scale(b2, 1/v);
                b3 = mat.vec_scale(b3, 1/v);
                
                // Cartesian q
                q_cart = [
                    q[0]*b1[0] + q[1]*b2[0] + q[2]*b3[0],
                    q[0]*b1[1] + q[1]*b2[1] + q[2]*b3[1],
                    q[0]*b1[2] + q[1]*b2[2] + q[2]*b3[2]
                ];
            }
        }

        let validOpsCount = 0;

        // Group operations by type for a cleaner UI
        for (let i = 0; i < this.cartesianOps.length; i++) {
            let op = this.cartesianOps[i];
            if (op.label === 'E (Identity)') continue; // Skip identity operation
            
            if (q_cart && lat) {
                // R_cart * q_cart
                let Rq_cart = [
                    op.R_cart[0][0]*q_cart[0] + op.R_cart[0][1]*q_cart[1] + op.R_cart[0][2]*q_cart[2],
                    op.R_cart[1][0]*q_cart[0] + op.R_cart[1][1]*q_cart[1] + op.R_cart[1][2]*q_cart[2],
                    op.R_cart[2][0]*q_cart[0] + op.R_cart[2][1]*q_cart[1] + op.R_cart[2][2]*q_cart[2]
                ];
                let dq_cart = [
                    Rq_cart[0] - q_cart[0],
                    Rq_cart[1] - q_cart[1],
                    Rq_cart[2] - q_cart[2]
                ];
                
                // Check if dq_cart is a reciprocal lattice vector
                let G_frac = [
                    mat.vec_dot(dq_cart, lat[0]),
                    mat.vec_dot(dq_cart, lat[1]),
                    mat.vec_dot(dq_cart, lat[2])
                ];
                
                let isInteger = (val) => Math.abs(val - Math.round(val)) < 1e-3;
                if (!isInteger(G_frac[0]) || !isInteger(G_frac[1]) || !isInteger(G_frac[2])) {
                    continue; // Skip operations not in the little group
                }
            }
            
            validOpsCount++;

            let angleDeg = Math.round(op.angle * 180 / Math.PI);
            let displayLabel = `#${i}: ${op.label}`;
            if (!op.isImproper && angleDeg > 0) {
                displayLabel += ` [${angleDeg}°]`;
            }
            this.dropdownEl.append(`<option value="${i}">${displayLabel}</option>`);
        }
        
        if (validOpsCount === 0) {
            this.dropdownEl.append(`<option value="-1">No applicable operations for this mode</option>`);
        }
    }

    // ─────────────────────────────────────────────
    // OPERATION SELECTION & ANIMATION
    // ─────────────────────────────────────────────

    selectOperation(idx) {
        if (idx < 0 || idx >= this.cartesianOps.length) return;
        this.currentOpIndex = idx;
        let op = this.cartesianOps[idx];
        if (!op) return;

        // Reset autoplay phase if active
        if (this.isAutoplaying) {
            this.autoplayPhase = 0;
        }

        // Reset sliders
        this.sliderRotValue = 0;
        this.sliderTransValue = 0;
        if (this.sliderRotEl) this.sliderRotEl.val(0);
        if (this.sliderTransEl) this.sliderTransEl.val(0);

        // Check if there is a translation component
        let hasTrans = (Math.abs(op.t_cart[0]) > 1e-4 || Math.abs(op.t_cart[1]) > 1e-4 || Math.abs(op.t_cart[2]) > 1e-4);
        let isIdentityRot = (op.label === 'E (Identity)' || (op.angle < 1e-4 && !op.isImproper));

        // Show/hide sliders
        if (this.sliderTransContainer) {
            if (hasTrans) this.sliderTransContainer.show();
            else this.sliderTransContainer.hide();
        }
        
        if (this.sliderRotContainer) {
            if (isIdentityRot && hasTrans) {
                this.sliderRotContainer.hide(); // Pure translation
            } else {
                this.sliderRotContainer.show();
                if (this.sliderRotLabel) {
                    this.sliderRotLabel.text(op.isImproper ? "Reflection / Inversion:" : "Rotation:");
                }
            }
        }

        // Update label
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

        // Set atoms and main arrows to equilibrium (t=0)
        this.applyInterpolatedOperation(op, 0, 0);

        this.refreshGhostBondsVisibility();
        this.drawSymmetryElement(op);

        this.crystal.needsRender = true;
        this.crystal.startAnimationLoop();
    }

    refreshGhostLattice() {
        if (!this.active) return;

        // Always ensure axes/planes are created/restored in case scene was cleared by vibcrystal.updatelocal
        this.createAxesAndPlanes();

        this.precomputeOperations();
        this.populateDropdown();
        
        let validOps = [];
        if (this.dropdownEl) {
            // Need to handle standard DOM element or jQuery depending on what dropdownEl is
            let options = this.dropdownEl[0].options;
            for (let i = 0; i < options.length; i++) {
                let val = parseInt(options[i].value, 10);
                if (val >= 0) validOps.push(val);
            }
        }

        if (validOps.length > 0 && !validOps.includes(this.currentOpIndex)) {
            this.currentOpIndex = validOps[0];
        } else if (validOps.length === 0) {
            this.currentOpIndex = -1;
        }

        if (this.dropdownEl && this.currentOpIndex >= 0) {
            this.dropdownEl.val(this.currentOpIndex);
        }

        if (this.currentOpIndex >= 0) {
            this.createGhostLattice();
            let op = this.cartesianOps[this.currentOpIndex];
            if (op) {
                this.applyInterpolatedOperation(op, this.sliderRotValue, this.sliderTransValue);
                this.drawSymmetryElement(op);
            }
        } else {
            this.removeGhostLattice();
        }

        this.refreshGhostBondsVisibility();
        this.crystal.needsRender = true;
        this.crystal.startAnimationLoop();
    }

    refreshGhostBondsVisibility() {
        if (!this.ghostMeshes) return;
        let show = this.bondsCheckboxEl ? this.bondsCheckboxEl.is(':checked') : true;
        for (let mesh of this.ghostMeshes) {
            if (mesh.name === 'symmetry-ghost-bond') {
                mesh.visible = show;
            }
        }
        this.crystal.needsRender = true;
        this.crystal.startAnimationLoop();
    }

    setSliderValues(tRot, tTrans) {
        this.sliderRotValue = Math.max(0, Math.min(1, tRot));
        this.sliderTransValue = Math.max(0, Math.min(1, tTrans));
        if (this.currentOpIndex < 0) return;

        let op = this.cartesianOps[this.currentOpIndex];
        this.applyInterpolatedOperation(op, this.sliderRotValue, this.sliderTransValue);

        this.crystal.needsRender = true;
        this.crystal.startAnimationLoop();
    }

    applyInterpolatedOperation(op, tRot, tTrans) {
        if (!this.crystal.atomobjects || !this.crystal.atompos) return;

        let nAtoms = this.crystal.atomobjects.length;
        let center = this.crystal.geometricCenter;

        let quat_full = new THREE.Quaternion();

        if (!op.isImproper) {
            let quat_target = new THREE.Quaternion();
            quat_target.setFromAxisAngle(op.axis, op.angle);
            let quat_identity = new THREE.Quaternion();
            quat_full.slerpQuaternions(quat_identity, quat_target, tRot);
        }

        for (let i = 0; i < nAtoms; i++) {
            let eqPos = this.crystal.atompos[i];
            let newPos = new THREE.Vector3();

            // 1. Shift to true crystallographic space (where origin is 0,0,0)
            let truePos = new THREE.Vector3().copy(eqPos).add(center);

            if (op.isImproper) {
                // Linear interpolation for improper operations
                let rx = truePos.x, ry = truePos.y, rz = truePos.z;
                let R = op.R_cart;
                
                let targetX = R[0][0]*rx + R[0][1]*ry + R[0][2]*rz;
                let targetY = R[1][0]*rx + R[1][1]*ry + R[1][2]*rz;
                let targetZ = R[2][0]*rx + R[2][1]*ry + R[2][2]*rz;

                truePos.set(
                    rx + tRot * (targetX - rx),
                    ry + tRot * (targetY - ry),
                    rz + tRot * (targetZ - rz)
                );
            } else {
                truePos.applyQuaternion(quat_full);
            }

            // 2. Shift back to screen space
            newPos.copy(truePos).sub(center);

            // 3. Apply translation phase
            let tx = op.t_cart[0], ty = op.t_cart[1], tz = op.t_cart[2];
            newPos.x += tTrans * tx;
            newPos.y += tTrans * ty;
            newPos.z += tTrans * tz;

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

            // Update main moving arrow if exists
            if (this.mainArrows && this.mainArrows[i] && this.crystal.arrows) {
                let vibrations = this.crystal.vibrationComponents[i];
                let snapTime = (typeof this.crystal.symmetryPhaseSnap === 'number') ? this.crystal.symmetryPhaseSnap : 0;
                let snapAngle = snapTime * 2.0 * Math.PI;
                let snapRe = Math.cos(snapAngle);
                let snapIm = Math.sin(snapAngle);

                let vx = snapRe * vibrations[0][0] - snapIm * vibrations[0][1];
                let vy = snapRe * vibrations[1][0] - snapIm * vibrations[1][1];
                let vz = snapRe * vibrations[2][0] - snapIm * vibrations[2][1];

                let v_orig = new THREE.Vector3(vx, vy, vz);
                let v_interp = new THREE.Vector3();

                if (op.isImproper) {
                    let R = op.R_cart;
                    let targetVx = R[0][0]*vx + R[0][1]*vy + R[0][2]*vz;
                    let targetVy = R[1][0]*vx + R[1][1]*vy + R[1][2]*vz;
                    let targetVz = R[2][0]*vx + R[2][1]*vy + R[2][2]*vz;
                    v_interp.set(
                        vx + tRot * (targetVx - vx),
                        vy + tRot * (targetVy - vy),
                        vz + tRot * (targetVz - vz)
                    );
                } else {
                    v_interp.copy(v_orig).applyQuaternion(quat_full);
                }

                let vlength = v_interp.length();
                if (vlength > 1e-10) {
                    let halfVisual = vlength * this.crystal.arrowScale * 0.5;
                    let nx = v_interp.x / vlength;
                    let ny = v_interp.y / vlength;
                    let nz = v_interp.z / vlength;
                    
                    this.mainArrows[i].position.set(
                        newPos.x + nx * halfVisual,
                        newPos.y + ny * halfVisual,
                        newPos.z + nz * halfVisual
                    );
                    this.mainArrows[i].scale.y = vlength * this.crystal.arrowScale;
                    this.mainArrows[i].quaternion.setFromUnitVectors(new THREE.Vector3(0,1,0), v_interp.normalize());
                } else {
                    this.mainArrows[i].scale.y = 0;
                }

                // Update ghost arrows position offset and scale dynamically if arrow scale changes
                let origLength = v_orig.length();
                let origHalfVisual = origLength * this.crystal.arrowScale * 0.5;
                if (origLength > 1e-10) {
                    let nxO = v_orig.x / origLength;
                    let nyO = v_orig.y / origLength;
                    let nzO = v_orig.z / origLength;
                    
                    if (this.initialGhostArrows && this.initialGhostArrows[i]) {
                        this.initialGhostArrows[i].scale.y = origLength * this.crystal.arrowScale;
                        this.initialGhostArrows[i].position.set(
                            eqPos.x + nxO * origHalfVisual,
                            eqPos.y + nyO * origHalfVisual,
                            eqPos.z + nzO * origHalfVisual
                        );
                    }
                } else {
                    if (this.initialGhostArrows && this.initialGhostArrows[i]) this.initialGhostArrows[i].scale.y = 0;
                }

                if (this.finalGhostArrows && this.finalGhostArrows[i]) {
                    let j = (this.currentMapping && this.currentMapping[i] !== undefined) ? this.currentMapping[i] : i;
                    let vibrations_j = this.crystal.vibrationComponents[j];
                    let vx_j = 0, vy_j = 0, vz_j = 0;
                    if (vibrations_j) {
                        vx_j = snapRe * vibrations_j[0][0] - snapIm * vibrations_j[0][1];
                        vy_j = snapRe * vibrations_j[1][0] - snapIm * vibrations_j[1][1];
                        vz_j = snapRe * vibrations_j[2][0] - snapIm * vibrations_j[2][1];
                    }
                    let v_j = new THREE.Vector3(vx_j, vy_j, vz_j);
                    let vlength_j = v_j.length();
                    
                    if (vlength_j > 1e-10) {
                        let nx_j = vx_j / vlength_j;
                        let ny_j = vy_j / vlength_j;
                        let nz_j = vz_j / vlength_j;
                        let halfVisual_j = vlength_j * this.crystal.arrowScale * 0.5;

                        this.finalGhostArrows[i].scale.y = vlength_j * this.crystal.arrowScale;
                        // finalPos was already calculated, but we must re-calculate it to set position correctly
                        let truePos = new THREE.Vector3().copy(eqPos).add(center);
                        let R = op.R_cart;
                        let t = op.t_cart;
                        let targetX = R[0][0]*truePos.x + R[0][1]*truePos.y + R[0][2]*truePos.z + t[0];
                        let targetY = R[1][0]*truePos.x + R[1][1]*truePos.y + R[1][2]*truePos.z + t[1];
                        let targetZ = R[2][0]*truePos.x + R[2][1]*truePos.y + R[2][2]*truePos.z + t[2];
                        let finalPos = new THREE.Vector3(targetX, targetY, targetZ).sub(center);
                        
                        this.finalGhostArrows[i].position.set(
                            finalPos.x + nx_j * halfVisual_j,
                            finalPos.y + ny_j * halfVisual_j,
                            finalPos.z + nz_j * halfVisual_j
                        );
                    } else {
                        this.finalGhostArrows[i].scale.y = 0;
                    }
                }
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
    // VISUALIZE SYMMETRY ELEMENTS (Bonus Idea)
    // ─────────────────────────────────────────────

    drawSymmetryElement(op) {
        if (this.symElementMesh) {
            this.crystal.scene.remove(this.symElementMesh);
            if (this.symElementMesh.geometry) this.symElementMesh.geometry.dispose();
            if (this.symElementMesh.material) this.symElementMesh.material.dispose();
            this.symElementMesh = null;
        }

        if (!op || op.label === 'E (Identity)') return;

        // The physical crystallographic origin is at -geometricCenter in viewer space
        let center = new THREE.Vector3().copy(this.crystal.geometricCenter).multiplyScalar(-1);

        if (op.label.startsWith('i')) {
            // Draw Inversion center as a glowing amber dot
            let geom = new THREE.SphereGeometry(0.2, 16, 16);
            let mat = new THREE.MeshPhongMaterial({ color: 0xffea00, emissive: 0xaa8800, shininess: 100, transparent: true, opacity: 0.9, depthWrite: false });
            let mesh = new THREE.Mesh(geom, mat);
            mesh.position.copy(center);
            this.symElementMesh = mesh;
            this.crystal.scene.add(mesh);
        } else if (op.label.startsWith('σ')) {
            // Draw Mirror plane as a semi-transparent cyan surface
            let geom = new THREE.PlaneGeometry(25, 25);
            let mat = new THREE.MeshBasicMaterial({ color: 0x00ffff, transparent: true, opacity: 0.35, side: THREE.DoubleSide, depthWrite: false });
            let mesh = new THREE.Mesh(geom, mat);
            
            // Align plane with the normal (op.axis)
            let zAxis = new THREE.Vector3(0, 0, 1);
            let quaternion = new THREE.Quaternion().setFromUnitVectors(zAxis, op.axis);
            mesh.quaternion.copy(quaternion);
            mesh.position.copy(center);

            this.symElementMesh = mesh;
            this.crystal.scene.add(mesh);
        } else if (op.label.startsWith('C') || op.label.startsWith('S')) {
            // Draw Rotation axis as a glowing magenta line/cylinder
            let geom = new THREE.CylinderGeometry(0.02, 0.02, 50, 12);
            let mat = new THREE.MeshBasicMaterial({ color: 0xff00ff, transparent: true, opacity: 0.5, depthWrite: false });
            let mesh = new THREE.Mesh(geom, mat);
            
            // Align cylinder with axis
            let yAxis = new THREE.Vector3(0, 1, 0);
            let quaternion = new THREE.Quaternion().setFromUnitVectors(yAxis, op.axis);
            mesh.quaternion.copy(quaternion);
            mesh.position.copy(center);

            this.symElementMesh = mesh;
            this.crystal.scene.add(mesh);
        }
    }

    // ─────────────────────────────────────────────
    // GHOST LATTICE
    // ─────────────────────────────────────────────

    /**
     * Helper to create an arrow mesh (mirroring vibcrystal's style)
     */
    createArrowMesh(color, opacity = 1.0) {
        let arrowGeometry = new THREE.CylinderGeometry(
            0,
            this.crystal.arrowHeadRadiusRatio * this.crystal.arrowRadius,
            this.crystal.arrowLength * this.crystal.arrowHeadLengthRatio,
            16, 1, true // openEnded to prevent overlapping caps in transparency
        );

        let axisGeometry = new THREE.CylinderGeometry(
            this.crystal.arrowRadius,
            this.crystal.arrowRadius,
            this.crystal.arrowLength,
            16, 1, true // openEnded to prevent overlapping caps in transparency
        );

        let AxisMaterial = new THREE.MeshLambertMaterial({
            color: color,
            transparent: opacity < 1.0,
            opacity: opacity,
            depthWrite: opacity < 1.0 ? false : true,
            blending: THREE.NormalBlending
        });

        let object = new THREE.Group();
        let axisMesh = new THREE.Mesh(axisGeometry, AxisMaterial);
        let arrowMesh = new THREE.Mesh(arrowGeometry, AxisMaterial);
        let length = (this.crystal.arrowLength + this.crystal.arrowLength * this.crystal.arrowHeadLengthRatio) / 2;

        arrowMesh.position.y = length;
        object.add(axisMesh);
        object.add(arrowMesh);
        return object;
    }

    /**
     * Compute the atom permutation map for a given symmetry operation.
     * Returns an array `mapping` where `mapping[i] = j`, meaning atom `i`
     * maps to atom `j` under the operation `R * r_i + t`.
     */
    getAtomMapping(op) {
        let phonon = this.crystal.phonon;
        if (!phonon || !phonon.atom_pos_red) return null;
        
        let mapping = [];
        let R = op.R_frac;
        let t = op.t_frac;
        let pos = phonon.atom_pos_red;
        let types = phonon.atomic_numbers || phonon.atom_numbers;
        
        for (let i = 0; i < pos.length; i++) {
            let r_i = pos[i];
            
            // R * r_i + t
            let rx = R[0][0]*r_i[0] + R[0][1]*r_i[1] + R[0][2]*r_i[2] + t[0];
            let ry = R[1][0]*r_i[0] + R[1][1]*r_i[1] + R[1][2]*r_i[2] + t[1];
            let rz = R[2][0]*r_i[0] + R[2][1]*r_i[1] + R[2][2]*r_i[2] + t[2];
            
            let match_j = i;
            let minDist = 1e9;
            
            for (let j = 0; j < pos.length; j++) {
                if (types && types[i] !== types[j]) continue;
                
                let r_j = pos[j];
                let dx = rx - r_j[0];
                let dy = ry - r_j[1];
                let dz = rz - r_j[2];
                
                // Wrap to primitive cell [-0.5, 0.5)
                dx = dx - Math.round(dx);
                dy = dy - Math.round(dy);
                dz = dz - Math.round(dz);
                
                let dist = dx*dx + dy*dy + dz*dz;
                if (dist < minDist) {
                    minDist = dist;
                    match_j = j;
                }
            }
            
            mapping[i] = match_j;
        }
        return mapping;
    }

    /**
     * Create semi-transparent "ghost" copies of all atoms at their
     * equilibrium positions as a visual reference.
     */
    createGhostLattice() {
        if (!this.crystal.atomobjects || !this.crystal.atompos) return;
        
        this.removeGhostLattice();

        if (this.ghostAtomsCheckboxEl && !this.ghostAtomsCheckboxEl.is(':checked')) {
            return;
        }

        if (this.currentOpIndex < 0) return;
        let op = this.cartesianOps[this.currentOpIndex];
        let R = op.R_cart;
        let t = op.t_cart;
        let center = this.crystal.geometricCenter;

        let sphereGeom = new THREE.SphereGeometry(0.3, 16, 12);
        let materialCache = {};

        // 1. Ghost Atoms at final positions
        let snapTime = (typeof this.crystal.symmetryPhaseSnap === 'number') ? this.crystal.symmetryPhaseSnap : 0;
        let snapAngle = snapTime * 2.0 * Math.PI;
        let snapRe = Math.cos(snapAngle);
        let snapIm = Math.sin(snapAngle);
        let vec_y = new THREE.Vector3(0, 1, 0);

        this.mainArrows = [];
        this.initialGhostArrows = [];
        this.finalGhostArrows = [];

        this.currentMapping = this.getAtomMapping(op);

        for (let i = 0; i < this.crystal.atompos.length; i++) {
            let eqPos = this.crystal.atompos[i];
            let atomNumber = this.crystal.atomobjects[i].atom_number;
            
            if (!materialCache[atomNumber]) {
                let colorHex = this.crystal.getAtomColorHex(atomNumber);
                materialCache[atomNumber] = new THREE.MeshLambertMaterial({
                    color: colorHex,
                    transparent: true,
                    opacity: 0.15,
                    depthWrite: false
                });
            }

            // --- INITIAL GHOST ATOM ---
            let initialMesh = new THREE.Mesh(sphereGeom, materialCache[atomNumber]);
            initialMesh.position.copy(eqPos);
            initialMesh.name = 'symmetry-ghost-initial';
            this.crystal.scene.add(initialMesh);
            this.ghostMeshes.push(initialMesh);

            // --- FINAL GHOST ATOM ---
            let truePos = new THREE.Vector3().copy(eqPos).add(center);
            let targetX = R[0][0]*truePos.x + R[0][1]*truePos.y + R[0][2]*truePos.z + t[0];
            let targetY = R[1][0]*truePos.x + R[1][1]*truePos.y + R[1][2]*truePos.z + t[1];
            let targetZ = R[2][0]*truePos.x + R[2][1]*truePos.y + R[2][2]*truePos.z + t[2];
            let finalPos = new THREE.Vector3(targetX, targetY, targetZ).sub(center);

            let finalMesh = new THREE.Mesh(sphereGeom, materialCache[atomNumber]);
            finalMesh.position.copy(finalPos);
            finalMesh.name = 'symmetry-ghost-final';
            this.crystal.scene.add(finalMesh);
            this.ghostMeshes.push(finalMesh);

            // Add arrows if enabled
            if (this.crystal.arrows && this.crystal.vibrationComponents && this.crystal.vibrationComponents[i]) {
                let vibrations = this.crystal.vibrationComponents[i];
                let vx = snapRe * vibrations[0][0] - snapIm * vibrations[0][1];
                let vy = snapRe * vibrations[1][0] - snapIm * vibrations[1][1];
                let vz = snapRe * vibrations[2][0] - snapIm * vibrations[2][1];

                let v = new THREE.Vector3(vx, vy, vz);
                let vlength = v.length(); // normalized (effAmp=1)
                
                // MAIN moving arrow (same color as original native arrows)
                let mainArrow = this.createArrowMesh(this.crystal.arrowcolor, 1.0);
                this.crystal.scene.add(mainArrow);
                this.mainArrows.push(mainArrow);
                // (position and rotation for mainArrow are set in applyInterpolatedOperation)

                let halfVisual = vlength * this.crystal.arrowScale * 0.5;

                // INITIAL GHOST ARROW (black, 0.15 opacity)
                let initialGhostArrow = this.createArrowMesh(0x000000, 0.15);
                if (vlength > 1e-10) {
                    let nx = vx / vlength;
                    let ny = vy / vlength;
                    let nz = vz / vlength;
                    initialGhostArrow.position.set(
                        eqPos.x + nx * halfVisual,
                        eqPos.y + ny * halfVisual,
                        eqPos.z + nz * halfVisual
                    );
                    initialGhostArrow.scale.y = vlength * this.crystal.arrowScale;
                    initialGhostArrow.quaternion.setFromUnitVectors(vec_y, v.normalize());
                } else {
                    initialGhostArrow.scale.y = 0;
                }
                initialGhostArrow.name = 'symmetry-ghost-arrow-initial';
                this.crystal.scene.add(initialGhostArrow);
                this.ghostMeshes.push(initialGhostArrow);
                if (!this.initialGhostArrows) this.initialGhostArrows = [];
                this.initialGhostArrows.push(initialGhostArrow);

                // FINAL GHOST ARROW (black, 0.15 opacity)
                let finalGhostArrow = this.createArrowMesh(0x000000, 0.15);
                
                // Get atom j that atom i maps to
                let j = this.currentMapping ? this.currentMapping[i] : i;
                let vibrations_j = this.crystal.vibrationComponents[j];
                let vx_j = vx, vy_j = vy, vz_j = vz;
                if (vibrations_j) {
                    vx_j = snapRe * vibrations_j[0][0] - snapIm * vibrations_j[0][1];
                    vy_j = snapRe * vibrations_j[1][0] - snapIm * vibrations_j[1][1];
                    vz_j = snapRe * vibrations_j[2][0] - snapIm * vibrations_j[2][1];
                }
                let v_j = new THREE.Vector3(vx_j, vy_j, vz_j);
                let vlength_j = v_j.length();
                let halfVisual_j = vlength_j * this.crystal.arrowScale * 0.5;

                if (vlength_j > 1e-10) {
                    let nx_j = vx_j / vlength_j;
                    let ny_j = vy_j / vlength_j;
                    let nz_j = vz_j / vlength_j;
                    finalGhostArrow.position.set(
                        finalPos.x + nx_j * halfVisual_j,
                        finalPos.y + ny_j * halfVisual_j,
                        finalPos.z + nz_j * halfVisual_j
                    );
                    finalGhostArrow.scale.y = vlength_j * this.crystal.arrowScale;
                    finalGhostArrow.quaternion.setFromUnitVectors(vec_y, v_j.normalize());
                } else {
                    finalGhostArrow.scale.y = 0;
                }
                finalGhostArrow.name = 'symmetry-ghost-arrow-final';
                this.crystal.scene.add(finalGhostArrow);
                this.ghostMeshes.push(finalGhostArrow);
                if (!this.finalGhostArrows) this.finalGhostArrows = [];
                this.finalGhostArrows.push(finalGhostArrow);
            }
        }

        // 2. Ghost Bonds at initial and final positions
        if (this.crystal.bonds && this.crystal.bonds.length > 0) {
            let bondMat = new THREE.MeshLambertMaterial({
                color: 0x666666,
                transparent: true,
                opacity: 0.15,
                depthWrite: false
            });
            let bondGeom = new THREE.CylinderGeometry(0.06, 0.06, 1, 6);
            let yAxis = new THREE.Vector3(0, 1, 0);

            for (let i = 0; i < this.crystal.bonds.length; i++) {
                let bond = this.crystal.bonds[i];
                if (!bond.a || !bond.b) continue;

                // Find equilibrium positions for this bond by matching the position object reference
                let aIndex = this.crystal.atomobjects.findIndex(atom => atom.position === bond.a);
                let bIndex = this.crystal.atomobjects.findIndex(atom => atom.position === bond.b);
                if (aIndex < 0 || bIndex < 0) continue;

                let eqA = this.crystal.atompos[aIndex];
                let eqB = this.crystal.atompos[bIndex];

                // INITIAL GHOST BOND
                let initialMidpoint = new THREE.Vector3().addVectors(eqA, eqB).multiplyScalar(0.5);
                let initialDir = new THREE.Vector3().subVectors(eqB, eqA);
                let initialLen = initialDir.length();
                initialDir.normalize();

                let initialBondMesh = new THREE.Mesh(bondGeom, bondMat);
                initialBondMesh.position.copy(initialMidpoint);
                initialBondMesh.scale.set(1, initialLen, 1);
                initialBondMesh.quaternion.setFromUnitVectors(yAxis, initialDir);
                initialBondMesh.name = 'symmetry-ghost-bond';
                this.crystal.scene.add(initialBondMesh);
                this.ghostMeshes.push(initialBondMesh);

                // Transform endpoints for FINAL GHOST BOND
                let trueA = new THREE.Vector3().copy(eqA).add(center);
                let targetAX = R[0][0]*trueA.x + R[0][1]*trueA.y + R[0][2]*trueA.z + t[0];
                let targetAY = R[1][0]*trueA.x + R[1][1]*trueA.y + R[1][2]*trueA.z + t[1];
                let targetAZ = R[2][0]*trueA.x + R[2][1]*trueA.y + R[2][2]*trueA.z + t[2];
                let finalA = new THREE.Vector3(targetAX, targetAY, targetAZ).sub(center);

                let trueB = new THREE.Vector3().copy(eqB).add(center);
                let targetBX = R[0][0]*trueB.x + R[0][1]*trueB.y + R[0][2]*trueB.z + t[0];
                let targetBY = R[1][0]*trueB.x + R[1][1]*trueB.y + R[1][2]*trueB.z + t[1];
                let targetBZ = R[2][0]*trueB.x + R[2][1]*trueB.y + R[2][2]*trueB.z + t[2];
                let finalB = new THREE.Vector3(targetBX, targetBY, targetBZ).sub(center);

                let midpoint = new THREE.Vector3().addVectors(finalA, finalB).multiplyScalar(0.5);
                let dir = new THREE.Vector3().subVectors(finalB, finalA);
                let len = dir.length();
                dir.normalize();

                let finalBondMesh = new THREE.Mesh(bondGeom, bondMat);
                finalBondMesh.position.copy(midpoint);
                finalBondMesh.scale.set(1, len, 1);
                finalBondMesh.quaternion.setFromUnitVectors(yAxis, dir);
                finalBondMesh.name = 'symmetry-ghost-bond';
                this.crystal.scene.add(finalBondMesh);
                this.ghostMeshes.push(finalBondMesh);
            }
        }

        // Add emphasis to the primitive unit cell
        if (this.crystal.phonon && this.crystal.phonon.lat) {
            let cellObj = createCellLineObject(this.crystal.phonon.lat, this.crystal.geometricCenter, 0x0284c7);
            cellObj.name = 'symmetry-ghost-cell';
            
            // Make the lines slightly thicker and transparent if possible
            cellObj.material.transparent = true;
            cellObj.material.opacity = 0.8;
            
            this.crystal.scene.add(cellObj);
            this.ghostMeshes.push(cellObj);
        }
    }

    removeGhostLattice() {
        for (let mesh of this.ghostMeshes) {
            this.crystal.scene.remove(mesh);
            mesh.traverse((child) => {
                if (child.geometry) child.geometry.dispose();
                if (child.material) child.material.dispose();
            });
        }
        this.ghostMeshes = [];

        if (this.mainArrows) {
            for (let mesh of this.mainArrows) {
                this.crystal.scene.remove(mesh);
                mesh.traverse((child) => {
                    if (child.geometry) child.geometry.dispose();
                    if (child.material) child.material.dispose();
                });
            }
            this.mainArrows = [];
        }
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
