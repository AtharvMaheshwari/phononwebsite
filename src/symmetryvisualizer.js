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
            this.isCameraSynced = false;
            
            // Handlers for two-way sync to avoid infinite loop
            let isSyncing = false;
            
            let syncMainToRef = () => {
                if (!this.isCameraSynced || isSyncing) return;
                if (this.crystal && this.refCrystal && this.crystal.camera && this.refCrystal.camera) {
                    isSyncing = true;
                    this.refCrystal.camera.position.copy(this.crystal.camera.position);
                    this.refCrystal.camera.quaternion.copy(this.crystal.camera.quaternion);
                    this.refCrystal.camera.up.copy(this.crystal.camera.up);
                    this.refCrystal.camera.zoom = this.crystal.camera.zoom;
                    this.refCrystal.camera.updateProjectionMatrix();
                    if (this.crystal.controls && this.refCrystal.controls) {
                        this.refCrystal.controls.target.copy(this.crystal.controls.target);
                        this.refCrystal.controls.update(); // Important for TrackballControls
                    }
                    this.refCrystal.needsRender = true;
                    this.refCrystal.startAnimationLoop();
                    isSyncing = false;
                }
            };
            
            let syncRefToMain = () => {
                if (!this.isCameraSynced || isSyncing) return;
                if (this.crystal && this.refCrystal && this.crystal.camera && this.refCrystal.camera) {
                    isSyncing = true;
                    this.crystal.camera.position.copy(this.refCrystal.camera.position);
                    this.crystal.camera.quaternion.copy(this.refCrystal.camera.quaternion);
                    this.crystal.camera.up.copy(this.refCrystal.camera.up);
                    this.crystal.camera.zoom = this.refCrystal.camera.zoom;
                    this.crystal.camera.updateProjectionMatrix();
                    if (this.crystal.controls && this.refCrystal.controls) {
                        this.crystal.controls.target.copy(this.refCrystal.controls.target);
                        this.crystal.controls.update(); // Important for TrackballControls
                    }
                    this.crystal.needsRender = true;
                    this.crystal.startAnimationLoop();
                    isSyncing = false;
                }
            };

            syncBtn.on('click', () => {
                this.isCameraSynced = !this.isCameraSynced;
                if (this.isCameraSynced) {
                    syncBtn.css({ 'background': '#f59e0b', 'color': 'white' }); // Amber to show active
                    syncBtn.text('Unsync Camera');
                    
                    // Do an initial sync
                    syncMainToRef();
                    
                    // Attach listeners if controls exist
                    if (this.crystal && this.crystal.controls) {
                        this.crystal.controls.addEventListener('change', syncMainToRef);
                    }
                    if (this.refCrystal && this.refCrystal.controls) {
                        this.refCrystal.controls.addEventListener('change', syncRefToMain);
                    }
                } else {
                    syncBtn.css({ 'background': '#10b981', 'color': 'white' }); // Green for inactive
                    syncBtn.text('Sync Camera');
                    
                    // Remove listeners
                    if (this.crystal && this.crystal.controls) {
                        this.crystal.controls.removeEventListener('change', syncMainToRef);
                    }
                    if (this.refCrystal && this.refCrystal.controls) {
                        this.refCrystal.controls.removeEventListener('change', syncRefToMain);
                    }
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
    classifyOperation(R_frac, t_frac, axis, hasResidualTranslation) {
        let detR = Math.round(
            R_frac[0][0]*(R_frac[1][1]*R_frac[2][2] - R_frac[2][1]*R_frac[1][2]) -
            R_frac[0][1]*(R_frac[1][0]*R_frac[2][2] - R_frac[1][2]*R_frac[2][0]) +
            R_frac[0][2]*(R_frac[1][0]*R_frac[2][1] - R_frac[1][1]*R_frac[2][0])
        );
        let tr = Math.round(R_frac[0][0] + R_frac[1][1] + R_frac[2][2]);

        // Only show + τ if there is a genuine residual translation (screw/glide)
        // that could not be absorbed by choosing a better pivot point
        let tLabel = hasResidualTranslation ? ' + τ' : '';
        
        let formatAxis = (ax) => {
            if (!ax) return "";
            
            // Try to convert to fractional Miller indices
            let phonon = this.crystal ? this.crystal.phonon : null;
            if (phonon && phonon.lat) {
                let A = mat.matrix_transpose(phonon.lat);
                let A_inv = mat.matrix_inverse(A);
                if (A_inv) {
                    let u = A_inv[0][0]*ax.x + A_inv[0][1]*ax.y + A_inv[0][2]*ax.z;
                    let v = A_inv[1][0]*ax.x + A_inv[1][1]*ax.y + A_inv[1][2]*ax.z;
                    let w = A_inv[2][0]*ax.x + A_inv[2][1]*ax.y + A_inv[2][2]*ax.z;
                    let u_orig = u;
                    let v_orig = v;
                    let w_orig = w;
                    
                    let maxVal = Math.max(Math.abs(u), Math.abs(v), Math.abs(w));
                    if (maxVal > 1e-6) {
                        u /= maxVal;
                        v /= maxVal;
                        w /= maxVal;
                    }
                    
                    let bestMult = 1;
                    for (let m = 1; m <= 12; m++) {
                        let diff = Math.abs(u*m - Math.round(u*m)) + Math.abs(v*m - Math.round(v*m)) + Math.abs(w*m - Math.round(w*m));
                        if (diff < 1e-4) {
                            bestMult = m;
                            break;
                        }
                    }
                    
                    let iu = Math.round(u * bestMult);
                    let iv = Math.round(v * bestMult);
                    let iw = Math.round(w * bestMult);
                    
                    // Check if the lattice is hexagonal/trigonal
                    let isHex = false;
                    let lat = phonon.lat;
                    if (lat) {
                        let a = Math.sqrt(lat[0][0]*lat[0][0] + lat[0][1]*lat[0][1] + lat[0][2]*lat[0][2]);
                        let b = Math.sqrt(lat[1][0]*lat[1][0] + lat[1][1]*lat[1][1] + lat[1][2]*lat[1][2]);
                        let dot = lat[0][0]*lat[1][0] + lat[0][1]*lat[1][1] + lat[0][2]*lat[1][2];
                        let gamma = Math.acos(dot / (a*b)) * 180 / Math.PI;
                        // Hexagonal: a=b, gamma=120 or 60
                        if (Math.abs(a - b) < 1e-3 && (Math.abs(gamma - 120) < 1e-3 || Math.abs(gamma - 60) < 1e-3)) {
                            isHex = true;
                        }
                    }
                    
                    if (isHex) {
                        let U = (2*u_orig - v_orig)/3;
                        let V = (2*v_orig - u_orig)/3;
                        let T = -(U + V);
                        let W = w_orig;
                        
                        let maxValHex = Math.max(Math.abs(U), Math.abs(V), Math.abs(T), Math.abs(W));
                        if (maxValHex > 1e-6) {
                            U /= maxValHex; V /= maxValHex; T /= maxValHex; W /= maxValHex;
                        }
                        
                        let bestMultHex = 1;
                        for (let m = 1; m <= 12; m++) {
                            let diff = Math.abs(U*m - Math.round(U*m)) + Math.abs(V*m - Math.round(V*m)) + Math.abs(T*m - Math.round(T*m)) + Math.abs(W*m - Math.round(W*m));
                            if (diff < 1e-4) { bestMultHex = m; break; }
                        }
                        let iU = Math.round(U * bestMultHex);
                        let iV = Math.round(V * bestMultHex);
                        let iT = Math.round(T * bestMultHex);
                        let iW = Math.round(W * bestMultHex);
                        
                        if (iU < 0 || (iU === 0 && iV < 0) || (iU === 0 && iV === 0 && iT < 0) || (iU === 0 && iV === 0 && iT === 0 && iW < 0)) {
                            iU = -iU; iV = -iV; iT = -iT; iW = -iW;
                        }
                        return `[${iU} ${iV} ${iT} ${iW}]`;
                    }
                    
                    // Standardize sign: first non-zero should be positive
                    if (iu < 0 || (iu === 0 && iv < 0) || (iu === 0 && iv === 0 && iw < 0)) {
                        iu = -iu; iv = -iv; iw = -iw;
                    }
                    
                    // Formatting negative numbers with an overline would be nice, but standard minus is fine
                    return `[${iu} ${iv} ${iw}]`;
                }
            }

            // Fallback to Cartesian
            let x = Math.abs(ax.x) < 1e-4 ? 0 : ax.x;
            let y = Math.abs(ax.y) < 1e-4 ? 0 : ax.y;
            let z = Math.abs(ax.z) < 1e-4 ? 0 : ax.z;
            return `[${x.toFixed(2)}, ${y.toFixed(2)}, ${z.toFixed(2)}]`;
        };
        let axLabel = axis ? ` ∥ ${formatAxis(axis)}` : '';
        let perpLabel = axis ? ` ⟂ ${formatAxis(axis)}` : '';

        if (detR === 1) {
            if (tr === 3)  return 'E (Identity)';
            if (tr === -1) return 'C₂ (180°)' + axLabel + tLabel;
            if (tr === 0)  return 'C₃ (120°)' + axLabel + tLabel;
            if (tr === 1)  return 'C₄ (90°)' + axLabel + tLabel;
            if (tr === 2)  return 'C₆ (60°)' + axLabel + tLabel;
        } else if (detR === -1) {
            if (tr === -3) return 'i (Inversion)' + tLabel;
            if (tr === 1)  return 'σ (Mirror)' + perpLabel + tLabel;
            if (tr === -2) return 'S₃ (120°)' + axLabel + tLabel;
            if (tr === -1) return 'S₄ (90°)' + axLabel + tLabel;
            if (tr === 0)  return 'S₆ (60°)' + axLabel + tLabel;
        }
        return 'Unknown operation';
    }

    // ─────────────────────────────────────────────
    // PRECOMPUTATION
    // ─────────────────────────────────────────────

    /**
     * Find the optimal pivot point p (in Cartesian) such that applying
     * the operation {R|t} about p eliminates as much translation as possible.
     *
     * The operation r → R·r + t about the origin is equivalent to
     * r → R·(r−p) + p + t_residual about pivot p, where:
     *   (I − R)·p = t − t_residual
     *
     * For inversion (R = −I): p = t/2, residual = 0
     * For reflection: absorb the component of t along the mirror normal,
     *   residual = component of t in the mirror plane (glide component)
     * For rotation: absorb the component of t perpendicular to the axis,
     *   residual = component of t along the axis (screw component)
     * For improper rotation (S_n): solve via pseudo-inverse
     */
    computePivotAndResidual(R_cart, t_cart, axis, angle, isImproper) {
        let tx = t_cart[0], ty = t_cart[1], tz = t_cart[2];
        let tMag = Math.sqrt(tx*tx + ty*ty + tz*tz);
        if (tMag < 1e-8) {
            // No translation — pivot at origin, no residual
            return { pivot: [0, 0, 0], t_residual: [0, 0, 0] };
        }

        let R = R_cart;
        let det = Math.round(
            R[0][0]*(R[1][1]*R[2][2] - R[2][1]*R[1][2]) -
            R[0][1]*(R[1][0]*R[2][2] - R[1][2]*R[2][0]) +
            R[0][2]*(R[1][0]*R[2][1] - R[1][1]*R[2][0])
        );
        let tr = R[0][0] + R[1][1] + R[2][2];
        let trRound = Math.round(tr);

        // Inversion: R = -I, so (I - R) = 2I → p = t/2
        if (det === -1 && trRound === -3) {
            return {
                pivot: [tx/2, ty/2, tz/2],
                t_residual: [0, 0, 0]
            };
        }

        // Reflection (det = -1, trace = 1): axis is the mirror normal
        // Component of t along normal can be absorbed (shift the plane)
        // Component of t in the plane is the glide and cannot be absorbed
        if (det === -1 && trRound === 1 && axis) {
            let nx = axis.x, ny = axis.y, nz = axis.z;
            // t_perp = (t · n) n  — absorbable by shifting the plane
            let t_dot_n = tx*nx + ty*ny + tz*nz;
            let t_perp = [t_dot_n * nx, t_dot_n * ny, t_dot_n * nz];
            // t_parallel = t - t_perp — glide component (residual)
            let t_par = [tx - t_perp[0], ty - t_perp[1], tz - t_perp[2]];
            // Pivot: shift along normal by t_perp/2
            let pivot = [t_perp[0]/2, t_perp[1]/2, t_perp[2]/2];
            return { pivot: pivot, t_residual: t_par };
        }

        // Proper rotation (det = 1, not identity): axis is the rotation axis
        // Component of t along axis is the screw (residual)
        // Component of t perpendicular to axis can be absorbed (shift the axis)
        if (det === 1 && trRound !== 3 && axis && angle > 1e-4) {
            let nx = axis.x, ny = axis.y, nz = axis.z;
            // t_along = (t · n) n — screw component (residual)
            let t_dot_n = tx*nx + ty*ny + tz*nz;
            let t_along = [t_dot_n * nx, t_dot_n * ny, t_dot_n * nz];
            // t_perp = t - t_along — can be absorbed
            let t_perp = [tx - t_along[0], ty - t_along[1], tz - t_along[2]];
            // Solve (I - R) * p = t_perp for the perpendicular components
            // Use the formula: p = [(I - R)^T (I - R)]^{-1} (I - R)^T t_perp
            // For a rotation by angle θ about axis n, the pseudo-inverse gives:
            // p = (1/2) t_perp + (1/2) cot(θ/2) (n × t_perp)
            let halfAngle = angle / 2;
            let cotHalf = Math.cos(halfAngle) / Math.sin(halfAngle);
            // n × t_perp
            let cross = [
                ny * t_perp[2] - nz * t_perp[1],
                nz * t_perp[0] - nx * t_perp[2],
                nx * t_perp[1] - ny * t_perp[0]
            ];
            let pivot = [
                0.5 * t_perp[0] + 0.5 * cotHalf * cross[0],
                0.5 * t_perp[1] + 0.5 * cotHalf * cross[1],
                0.5 * t_perp[2] + 0.5 * cotHalf * cross[2]
            ];
            return { pivot: pivot, t_residual: t_along };
        }

        // Improper rotation S_n (det = -1, not mirror, not inversion)
        // For S_n, the proper part is C_n, and the reflection is across the plane ⊥ axis.
        // We can solve (I - R) p = t using the pseudo-inverse.
        // (I - R) for S_n is always invertible (no eigenvalue 1), so we can solve directly.
        if (det === -1 && trRound !== -3 && trRound !== 1) {
            let ImR = [
                [1 - R[0][0], -R[0][1], -R[0][2]],
                [-R[1][0], 1 - R[1][1], -R[1][2]],
                [-R[2][0], -R[2][1], 1 - R[2][2]]
            ];
            let ImR_inv = mat.matrix_inverse(ImR);
            if (ImR_inv) {
                let pivot = [
                    ImR_inv[0][0]*tx + ImR_inv[0][1]*ty + ImR_inv[0][2]*tz,
                    ImR_inv[1][0]*tx + ImR_inv[1][1]*ty + ImR_inv[1][2]*tz,
                    ImR_inv[2][0]*tx + ImR_inv[2][1]*ty + ImR_inv[2][2]*tz
                ];
                return { pivot: pivot, t_residual: [0, 0, 0] };
            }
        }

        // Fallback: no optimization possible
        return { pivot: [0, 0, 0], t_residual: [tx, ty, tz] };
    }

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

            // The full Cartesian translation (for atom mapping / final position)
            let t_cart_full = this.fractionalToCartesianTranslation(t_frac_shifted, lat);

            // The INTRINSIC crystallographic translation (from spglib, no centering)
            // This is what determines screw/glide character
            let t_cart_intrinsic = this.fractionalToCartesianTranslation(t_frac, lat);

            // The centering shift (a pure lattice vector, not part of the symmetry)
            let T_cart = this.fractionalToCartesianTranslation(T_frac, lat);

            // Extract axis and angle
            let { axis, angle, isImproper } = this.extractAxisAngle(R_cart);

            // Compute optimal pivot point from INTRINSIC translation only
            // This way pure rotations/reflections (t_frac=0) get pivot=0, residual=0
            let { pivot, t_residual } = this.computePivotAndResidual(R_cart, t_cart_intrinsic, axis, angle, isImproper);

            // Check if residual translation is significant
            let residualMag = Math.sqrt(t_residual[0]*t_residual[0] + t_residual[1]*t_residual[1] + t_residual[2]*t_residual[2]);
            let hasResidualTranslation = residualMag > 1e-4;

            // Human-readable label (only shows + τ for genuine screws/glides)
            let label = this.classifyOperation(R_frac, t_frac, axis, hasResidualTranslation);

            this.cartesianOps.push({
                R_frac,
                t_frac,
                R_cart,
                t_cart: t_cart_full,
                pivot_cart: pivot,
                t_residual: t_residual,
                T_centering: T_cart,
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
                // The correct little group condition in fractional space:
                // A symmetry operation R (integer matrix in direct-lattice basis) acts on
                // a q-vector in fractional reciprocal coordinates as q → R^{-T} q.
                // The operation is in the little group iff R^{-T} q - q is an integer vector
                // (i.e., a reciprocal lattice vector G).
                let q = this.crystal.phononweb.phonon.kpoints[this.crystal.phononweb.k];
                let R = op.R_frac;
                let R_inv = mat.matrix_inverse(R);
                if (!R_inv) continue;
                let Rinv_T = mat.matrix_transpose(R_inv);

                let Rq = [
                    Rinv_T[0][0]*q[0] + Rinv_T[0][1]*q[1] + Rinv_T[0][2]*q[2],
                    Rinv_T[1][0]*q[0] + Rinv_T[1][1]*q[1] + Rinv_T[1][2]*q[2],
                    Rinv_T[2][0]*q[0] + Rinv_T[2][1]*q[1] + Rinv_T[2][2]*q[2]
                ];
                let dq = [Rq[0] - q[0], Rq[1] - q[1], Rq[2] - q[2]];

                let isInteger = (val) => Math.abs(val - Math.round(val)) < 1e-2;
                if (!isInteger(dq[0]) || !isInteger(dq[1]) || !isInteger(dq[2])) {
                    continue; // Skip operations not in the little group
                }
            }

            
            validOpsCount++;

            let angleDeg = Math.round(op.angle * 180 / Math.PI);
            let displayLabel = `${op.label}`;
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

        // Check if there is a residual translation component (genuine screw/glide)
        let tr = op.t_residual || [0,0,0];
        let hasTrans = (Math.abs(tr[0]) > 1e-4 || Math.abs(tr[1]) > 1e-4 || Math.abs(tr[2]) > 1e-4);
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
                    let rotText = "Rotation:";
                    if (op.isImproper) {
                        if (op.label.startsWith('i')) rotText = "Inversion:";
                        else if (op.label.startsWith('σ')) rotText = "Reflection:";
                        else rotText = "Improper Rotation:";
                    }
                    this.sliderRotLabel.text(rotText);
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

        // Pivot point in true crystallographic space
        let pivot = op.pivot_cart || [0, 0, 0];
        let pivotVec = new THREE.Vector3(pivot[0], pivot[1], pivot[2]);

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

            // 2. Translate to pivot frame
            truePos.sub(pivotVec);

            if (op.isImproper) {
                // Linear interpolation for improper operations (about pivot)
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

            // 3. Translate back from pivot frame
            truePos.add(pivotVec);

            // 4. Shift back to screen space
            newPos.copy(truePos).sub(center);

            // 5. Apply residual translation phase (only genuine screws/glides)
            let t_res = op.t_residual || [0, 0, 0];
            newPos.x += tTrans * t_res[0];
            newPos.y += tTrans * t_res[1];
            newPos.z += tTrans * t_res[2];

            // 6. Apply centering lattice shift (smooth, proportional to rotation)
            // This is a pure lattice vector that ensures atoms map to the
            // nearest periodic image, not part of the symmetry operation itself
            let T = op.T_centering || [0, 0, 0];
            newPos.x += tRot * T[0];
            newPos.y += tRot * T[1];
            newPos.z += tRot * T[2];

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
                    let j = (this.currentMapping && this.currentMapping[i] !== undefined) ? this.currentMapping[i] : undefined;
                    
                    // Use ADAPTED vibration components so degenerate bands scale cleanly
                    let vibrations_j = j !== undefined ? this.currentAdaptedComponents[j] : null;
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

        // The physical crystallographic origin is at -geometricCenter in viewer space.
        // The pivot in viewer space = pivot_cart - geometricCenter
        let geomCenter = this.crystal.geometricCenter;
        let pivot = op.pivot_cart || [0, 0, 0];
        let pivotInViewerSpace = new THREE.Vector3(
            pivot[0] - geomCenter.x,
            pivot[1] - geomCenter.y,
            pivot[2] - geomCenter.z
        );

        if (op.label.startsWith('i')) {
            // Draw Inversion center as a glowing amber dot at the pivot
            let geom = new THREE.SphereGeometry(0.2, 16, 16);
            let mat = new THREE.MeshPhongMaterial({ color: 0xffea00, emissive: 0xaa8800, shininess: 100, transparent: true, opacity: 0.9, depthWrite: false });
            let mesh = new THREE.Mesh(geom, mat);
            mesh.position.copy(pivotInViewerSpace);
            this.symElementMesh = mesh;
            this.crystal.scene.add(mesh);
        } else if (op.label.startsWith('σ')) {
            // Draw Mirror plane as a semi-transparent cyan surface at the pivot
            let geom = new THREE.PlaneGeometry(25, 25);
            let mat = new THREE.MeshBasicMaterial({ color: 0x00ffff, transparent: true, opacity: 0.35, side: THREE.DoubleSide, depthWrite: false });
            let mesh = new THREE.Mesh(geom, mat);
            
            // Align plane with the normal (op.axis)
            let zAxis = new THREE.Vector3(0, 0, 1);
            let quaternion = new THREE.Quaternion().setFromUnitVectors(zAxis, op.axis);
            mesh.quaternion.copy(quaternion);
            mesh.position.copy(pivotInViewerSpace);

            this.symElementMesh = mesh;
            this.crystal.scene.add(mesh);
        } else if (op.label.startsWith('C') || op.label.startsWith('S')) {
            // Draw Rotation axis as a glowing magenta line/cylinder, passing through the pivot
            let geom = new THREE.CylinderGeometry(0.02, 0.02, 50, 12);
            let mat = new THREE.MeshBasicMaterial({ color: 0xff00ff, transparent: true, opacity: 0.5, depthWrite: false });
            let mesh = new THREE.Mesh(geom, mat);
            
            // Align cylinder with axis, centered on the pivot
            let yAxis = new THREE.Vector3(0, 1, 0);
            let quaternion = new THREE.Quaternion().setFromUnitVectors(yAxis, op.axis);
            mesh.quaternion.copy(quaternion);
            mesh.position.copy(pivotInViewerSpace);

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

    // ─────────────────────────────────────────────
    // DEGENERATE BAND ADAPTATION (D-matrix)
    // ─────────────────────────────────────────────

    /**
     * For degenerate bands, build the representation matrix D for the selected
     * symmetry operation, project it into the degenerate subspace, diagonalize,
     * and return symmetry-adapted (circularly polarized) vibration components.
     *
     * Returns adapted vibrationComponents array, or the original if no
     * adaptation is needed.
     */
    getAdaptedVibrationComponents(op) {
        let phonon = this.crystal.phonon;
        let phononweb = this.crystal.phononweb;
        if (!phonon || !phononweb || !phonon.vec || !phonon.eigenvalues) {
            return this.crystal.vibrationComponents;
        }

        let k = phononweb.k;
        let n = phononweb.n;
        if (!phonon.eigenvalues[k]) return this.crystal.vibrationComponents;

        let freqs = phonon.eigenvalues[k];
        let currentFreq = freqs[n];

        // Find degenerate manifold (all modes within tolerance of the current mode)
        let tol = 0.1; // cm⁻¹
        let manifold = [];
        for (let m = 0; m < freqs.length; m++) {
            if (Math.abs(freqs[m] - currentFreq) < tol) {
                manifold.push(m);
            }
        }

        // If non-degenerate, no adaptation needed
        if (manifold.length <= 1) {
            return this.crystal.vibrationComponents;
        }

        // Build the full representation matrix D for this symmetry operation.
        // D is a (3*natoms x 3*natoms) complex matrix.
        let natoms = phonon.natoms;
        let D = this._buildRepresentationMatrix(op, k);
        if (!D) return this.crystal.vibrationComponents;

        // Get eigenvectors for all modes in the manifold.
        // phonon.vec[k][m][atom][xyz] = [re, im]
        let basis = []; // array of complex flat vectors, each length 3*natoms
        for (let mi = 0; mi < manifold.length; mi++) {
            let m = manifold[mi];
            let vec_m = phonon.vec[k][m]; // [natoms][3][2]
            let flat = new Array(3 * natoms * 2); // interleaved [re, im, re, im, ...]
            for (let a = 0; a < natoms; a++) {
                for (let d = 0; d < 3; d++) {
                    flat[(a * 3 + d) * 2] = vec_m[a][d][0];     // re
                    flat[(a * 3 + d) * 2 + 1] = vec_m[a][d][1]; // im
                }
            }
            basis.push(flat);
        }

        // Project D into the degenerate subspace: M_D[a][b] = <basis[a] | D | basis[b]>
        let mSize = manifold.length;
        let M_D = []; // mSize x mSize complex (stored as [re, im])
        for (let a = 0; a < mSize; a++) {
            M_D[a] = [];
            for (let b = 0; b < mSize; b++) {
                // D * basis[b]
                let Db = this._complexMatVec(D, basis[b], 3 * natoms);
                // <basis[a] | D*basis[b]>
                let dot = this._complexDot(basis[a], Db, 3 * natoms);
                M_D[a][b] = dot; // [re, im]
            }
        }

        // Diagonalize M_D to find symmetry-adapted eigenvectors
        let { eigenvalues, eigenvectors } = this._complexEig(M_D, mSize);

        // Find which adapted mode corresponds to the currently selected mode n.
        // The adapted mode that has the largest overlap with original mode n.
        let nIdx = manifold.indexOf(n);
        if (nIdx < 0) return this.crystal.vibrationComponents;

        // Construct the adapted eigenvector for mode n:
        // adapted_vec = sum_a eigenvectors[bestIdx][a] * basis[a]
        // We pick the eigenvector with max overlap with the original mode n.
        let bestIdx = 0;
        let bestOverlap = 0;
        for (let ei = 0; ei < mSize; ei++) {
            // Overlap = |sum_a conj(eigvec[ei][a]) * delta(a, nIdx)|
            // = |eigvec[ei][nIdx]|
            let re = eigenvectors[ei][nIdx][0];
            let im = eigenvectors[ei][nIdx][1];
            let overlap = re * re + im * im;
            if (overlap > bestOverlap) {
                bestOverlap = overlap;
                bestIdx = ei;
            }
        }

        // Now reconstruct vibrationComponents using the adapted eigenvector.
        // The adapted primitive-cell eigenvector:
        let adaptedPrim = new Array(natoms);
        for (let a = 0; a < natoms; a++) {
            adaptedPrim[a] = [[0, 0], [0, 0], [0, 0]];
            for (let d = 0; d < 3; d++) {
                let re_sum = 0, im_sum = 0;
                for (let mi = 0; mi < mSize; mi++) {
                    let coeff_re = eigenvectors[bestIdx][mi][0];
                    let coeff_im = eigenvectors[bestIdx][mi][1];
                    let vec_re = basis[mi][(a * 3 + d) * 2];
                    let vec_im = basis[mi][(a * 3 + d) * 2 + 1];
                    // complex multiply: (coeff) * (vec)
                    re_sum += coeff_re * vec_re - coeff_im * vec_im;
                    im_sum += coeff_re * vec_im + coeff_im * vec_re;
                }
                adaptedPrim[a][d] = [re_sum, im_sum];
            }
        }

        // Rebuild full vibrationComponents for the supercell using the adapted vector
        // (replicate the same logic as getVibrations in phononwebpage.js)
        let kpt = phonon.kpoints[k];
        let nx = phononweb.nx || 1;
        let ny = phononweb.ny || 1;
        let nz = phononweb.nz || 1;
        let nx_int = parseInt(nx);
        let ny_int = parseInt(ny);
        let nz_int = parseInt(nz);
        let ix_start = -Math.floor(nx_int / 2);
        let iy_start = -Math.floor(ny_int / 2);
        let iz_start = -Math.floor(nz_int / 2);

        let adaptedComponents = [];

        for (let ix = ix_start; ix < ix_start + nx_int; ix++) {
            for (let iy = iy_start; iy < iy_start + ny_int; iy++) {
                for (let iz = iz_start; iz < iz_start + nz_int; iz++) {
                    for (let a = 0; a < natoms; a++) {
                        let atom_phase = 0;
                        if (phonon.addatomphase) {
                            atom_phase = mat.vec_dot(kpt, phonon.atom_pos_red[a]);
                        }
                        let sprod = mat.vec_dot(kpt, [ix, iy, iz]) + atom_phase;
                        let angle = sprod * 2.0 * Math.PI;
                        let phase_re = Math.cos(angle);
                        let phase_im = Math.sin(angle);

                        let comp = [];
                        for (let d = 0; d < 3; d++) {
                            let v_re = adaptedPrim[a][d][0];
                            let v_im = adaptedPrim[a][d][1];
                            // (v) * phase
                            let out_re = v_re * phase_re - v_im * phase_im;
                            let out_im = v_re * phase_im + v_im * phase_re;
                            comp.push([out_re, out_im]);
                        }
                        adaptedComponents.push(comp);
                    }
                }
            }
        }

        return adaptedComponents;
    }

    /**
     * Build the (3N x 3N) representation matrix D for a symmetry operation at q-point k.
     * D is stored as a flat array of interleaved [re, im] pairs, row-major.
     * D[(i*dim + j)*2] = re, D[(i*dim + j)*2 + 1] = im
     */
    _buildRepresentationMatrix(op, k) {
        let phonon = this.crystal.phonon;
        if (!phonon || !phonon.atom_pos_red || !phonon.lat) return null;

        let natoms = phonon.natoms;
        let dim = 3 * natoms;
        let lat = phonon.lat; // 3x3 row-major
        let pos = phonon.atom_pos_red;
        let q = phonon.kpoints[k];

        let R_frac = op.R_frac;
        let tau = op.t_frac;

        // A = lat^T (column vectors of lattice)
        let A = mat.matrix_transpose(lat);
        let A_inv = mat.matrix_inverse(A);
        if (!A_inv) return null;

        // R_cart = A * R_frac * A_inv
        let R_cart = mat.matrix_multiply(A, mat.matrix_multiply(R_frac, A_inv));

        // D is dim x dim complex, stored interleaved
        let D = new Float64Array(dim * dim * 2);

        // G = round(R_frac^T_inv * q - q) ... actually G = round(R_q * q - q)
        // where R_q = (R^-1)^T
        let R_inv = mat.matrix_inverse(R_frac);
        if (!R_inv) return null;
        let R_q = mat.matrix_transpose(R_inv);
        let Gx = R_q[0][0]*q[0] + R_q[0][1]*q[1] + R_q[0][2]*q[2] - q[0];
        let Gy = R_q[1][0]*q[0] + R_q[1][1]*q[1] + R_q[1][2]*q[2] - q[1];
        let Gz = R_q[2][0]*q[0] + R_q[2][1]*q[1] + R_q[2][2]*q[2] - q[2];
        let G = [Math.round(Gx), Math.round(Gy), Math.round(Gz)];

        for (let i = 0; i < natoms; i++) {
            // r_prime = R * pos[i] + tau
            let rp0 = R_frac[0][0]*pos[i][0] + R_frac[0][1]*pos[i][1] + R_frac[0][2]*pos[i][2] + tau[0];
            let rp1 = R_frac[1][0]*pos[i][0] + R_frac[1][1]*pos[i][1] + R_frac[1][2]*pos[i][2] + tau[1];
            let rp2 = R_frac[2][0]*pos[i][0] + R_frac[2][1]*pos[i][1] + R_frac[2][2]*pos[i][2] + tau[2];

            // Find j such that pos[j] ≡ r_prime (mod 1)
            let best_j = 0;
            let best_dist = 1e9;
            for (let j = 0; j < natoms; j++) {
                let dx = rp0 - pos[j][0];
                let dy = rp1 - pos[j][1];
                let dz = rp2 - pos[j][2];
                dx -= Math.round(dx);
                dy -= Math.round(dy);
                dz -= Math.round(dz);
                let dist = dx*dx + dy*dy + dz*dz;
                if (dist < best_dist) {
                    best_dist = dist;
                    best_j = j;
                }
            }

            // phase = exp(-2πi G · pos[j])
            let Gdot = G[0]*pos[best_j][0] + G[1]*pos[best_j][1] + G[2]*pos[best_j][2];
            let phase_angle = -2.0 * Math.PI * Gdot;
            let phase_re = Math.cos(phase_angle);
            let phase_im = Math.sin(phase_angle);

            // D[3j:3j+3, 3i:3i+3] = phase * R_cart
            for (let r = 0; r < 3; r++) {
                for (let c = 0; c < 3; c++) {
                    let row = 3 * best_j + r;
                    let col = 3 * i + c;
                    D[(row * dim + col) * 2]     = phase_re * R_cart[r][c];
                    D[(row * dim + col) * 2 + 1] = phase_im * R_cart[r][c];
                }
            }
        }

        return D;
    }

    /**
     * Complex matrix-vector multiply: result = M * v
     * M is dim x dim stored interleaved, v is dim stored interleaved.
     */
    _complexMatVec(M, v, dim) {
        let result = new Array(dim * 2);
        for (let i = 0; i < dim; i++) {
            let re = 0, im = 0;
            for (let j = 0; j < dim; j++) {
                let m_re = M[(i * dim + j) * 2];
                let m_im = M[(i * dim + j) * 2 + 1];
                let v_re = v[j * 2];
                let v_im = v[j * 2 + 1];
                re += m_re * v_re - m_im * v_im;
                im += m_re * v_im + m_im * v_re;
            }
            result[i * 2] = re;
            result[i * 2 + 1] = im;
        }
        return result;
    }

    /**
     * Complex dot product: <a|b> = sum conj(a_i) * b_i
     */
    _complexDot(a, b, dim) {
        let re = 0, im = 0;
        for (let i = 0; i < dim; i++) {
            let a_re = a[i * 2], a_im = a[i * 2 + 1];
            let b_re = b[i * 2], b_im = b[i * 2 + 1];
            // conj(a) * b = (a_re - i*a_im)(b_re + i*b_im)
            re += a_re * b_re + a_im * b_im;
            im += a_re * b_im - a_im * b_re;
        }
        return [re, im];
    }

    /**
     * Eigendecomposition of a small NxN complex matrix via QR iteration.
     * For typical degenerate multiplets N = 2 or 3.
     * Returns { eigenvalues: [[re,im],...], eigenvectors: [[[re,im],...],...] }
     * eigenvectors[i] is the i-th eigenvector as array of [re,im] components.
     */
    _complexEig(M, n) {
        if (n === 1) {
            return {
                eigenvalues: [M[0][0]],
                eigenvectors: [[[1, 0]]]
            };
        }
        if (n === 2) {
            return this._complexEig2x2(M);
        }
        // General case: power iteration / direct for small matrices
        return this._complexEigGeneral(M, n);
    }

    /**
     * Analytic eigendecomposition of a 2x2 complex matrix.
     */
    _complexEig2x2(M) {
        let a = M[0][0], b = M[0][1], c = M[1][0], d = M[1][1];
        // trace = a + d
        let tr_re = a[0] + d[0], tr_im = a[1] + d[1];
        // det = a*d - b*c
        let det_re = (a[0]*d[0] - a[1]*d[1]) - (b[0]*c[0] - b[1]*c[1]);
        let det_im = (a[0]*d[1] + a[1]*d[0]) - (b[0]*c[1] + b[1]*c[0]);
        // discriminant = tr^2 - 4*det
        let disc_re = (tr_re*tr_re - tr_im*tr_im) - 4*det_re;
        let disc_im = 2*tr_re*tr_im - 4*det_im;
        // sqrt(disc)
        let [sqrt_re, sqrt_im] = this._complexSqrt(disc_re, disc_im);

        let eigvals = [
            [(tr_re + sqrt_re) / 2, (tr_im + sqrt_im) / 2],
            [(tr_re - sqrt_re) / 2, (tr_im - sqrt_im) / 2]
        ];

        let eigvecs = [];
        for (let i = 0; i < 2; i++) {
            let lambda = eigvals[i];
            // (A - lambda*I) * v = 0
            // Use first row: (a - lambda)*v0 + b*v1 = 0
            // v = [b, lambda - a] (unnormalized)
            let v0_re = b[0], v0_im = b[1];
            let v1_re = lambda[0] - a[0], v1_im = lambda[1] - a[1];
            let norm = Math.sqrt(v0_re*v0_re + v0_im*v0_im + v1_re*v1_re + v1_im*v1_im);
            if (norm < 1e-15) {
                // Fallback: use second row
                v0_re = lambda[0] - d[0]; v0_im = lambda[1] - d[1];
                v1_re = c[0]; v1_im = c[1];
                norm = Math.sqrt(v0_re*v0_re + v0_im*v0_im + v1_re*v1_re + v1_im*v1_im);
            }
            if (norm < 1e-15) {
                eigvecs.push([[1, 0], [0, 0]]);
            } else {
                eigvecs.push([
                    [v0_re / norm, v0_im / norm],
                    [v1_re / norm, v1_im / norm]
                ]);
            }
        }

        return { eigenvalues: eigvals, eigenvectors: eigvecs };
    }

    /**
     * General eigendecomposition for NxN complex matrix.
     * Uses inverse iteration to find eigenvectors after computing eigenvalues
     * via the characteristic polynomial for small N, or direct iteration.
     */
    _complexEigGeneral(M, n) {
        // For small degenerate multiplets (typically n=3), we use a direct approach:
        // Find eigenvalues by iterative QR, then eigenvectors by inverse iteration.

        // Since M_D is (approximately) unitary, its eigenvalues lie on the unit circle.
        // We can find them by diagonalizing M_D using Jacobi-like rotations,
        // but for simplicity, use a direct power-method approach for each eigenvector.

        let eigvals = [];
        let eigvecs = [];

        // Work with a copy
        let A = [];
        for (let i = 0; i < n; i++) {
            A[i] = [];
            for (let j = 0; j < n; j++) {
                A[i][j] = [M[i][j][0], M[i][j][1]];
            }
        }

        // Simple Schur-like iteration: repeatedly find eigenvalue/eigenvector and deflate
        for (let found = 0; found < n; found++) {
            let dim = n - found;
            // Power iteration on current matrix to find dominant eigenvector
            let v = [];
            for (let i = 0; i < dim; i++) v[i] = [i === 0 ? 1 : 0, 0];

            let lambda = [0, 0];
            for (let iter = 0; iter < 200; iter++) {
                // w = A * v
                let w = [];
                for (let i = 0; i < dim; i++) {
                    let re = 0, im = 0;
                    for (let j = 0; j < dim; j++) {
                        re += A[i][j][0] * v[j][0] - A[i][j][1] * v[j][1];
                        im += A[i][j][0] * v[j][1] + A[i][j][1] * v[j][0];
                    }
                    w[i] = [re, im];
                }
                // Rayleigh quotient: lambda = <v|w> / <v|v>
                let num_re = 0, num_im = 0;
                for (let i = 0; i < dim; i++) {
                    num_re += v[i][0] * w[i][0] + v[i][1] * w[i][1];
                    num_im += v[i][0] * w[i][1] - v[i][1] * w[i][0];
                }
                lambda = [num_re, num_im];
                // Normalize w
                let norm = 0;
                for (let i = 0; i < dim; i++) norm += w[i][0]*w[i][0] + w[i][1]*w[i][1];
                norm = Math.sqrt(norm);
                if (norm < 1e-15) break;
                for (let i = 0; i < dim; i++) { w[i][0] /= norm; w[i][1] /= norm; }
                v = w;
            }

            eigvals.push(lambda);

            // Map back to full space
            let fullVec = [];
            // The deflated matrix indices map to the remaining original indices
            // For simplicity, we track a permutation
            if (found === 0) {
                for (let i = 0; i < n; i++) fullVec[i] = i < v.length ? v[i] : [0, 0];
            } else {
                for (let i = 0; i < n; i++) fullVec[i] = [0, 0];
                for (let i = 0; i < dim; i++) {
                    fullVec[found + i] = v[i];
                }
            }
            eigvecs.push(fullVec);

            if (dim <= 1) break;

            // Deflate: A' = A - lambda * v * v^H (Hotelling deflation)
            // Actually, for better stability, do a Householder deflation
            // But for small matrices (n<=3), direct deflation is fine
            let newA = [];
            let newDim = dim - 1;
            // Build projector P = I - v*v^H, then A' = P*A*P
            // But simpler: use similarity transform to put v as first column
            // For now, use simple Wielandt deflation:
            // A_new[i][j] = A[i][j] - lambda * v[i] * conj(v[j])
            for (let i = 0; i < dim; i++) {
                newA[i] = [];
                for (let j = 0; j < dim; j++) {
                    // lambda * v[i] * conj(v[j])
                    let vi_re = v[i][0], vi_im = v[i][1];
                    let vj_re = v[j][0], vj_im = -v[j][1]; // conj
                    let prod_re = vi_re * vj_re - vi_im * vj_im;
                    let prod_im = vi_re * vj_im + vi_im * vj_re;
                    let lp_re = lambda[0] * prod_re - lambda[1] * prod_im;
                    let lp_im = lambda[0] * prod_im + lambda[1] * prod_re;
                    newA[i][j] = [A[i][j][0] - lp_re, A[i][j][1] - lp_im];
                }
            }
            // Reduce dimension by removing the row/col most aligned with v
            // For simplicity, remove row 0, col 0 after rotating v to e_0
            // Actually just use the deflated matrix directly
            A = [];
            for (let i = 1; i < dim; i++) {
                A[i - 1] = [];
                for (let j = 1; j < dim; j++) {
                    A[i - 1][j - 1] = newA[i][j];
                }
            }
        }

        return { eigenvalues: eigvals, eigenvectors: eigvecs };
    }

    /**
     * Complex square root of (a + bi).
     */
    _complexSqrt(a, b) {
        let r = Math.sqrt(a * a + b * b);
        let re = Math.sqrt((r + a) / 2);
        let im = (b >= 0 ? 1 : -1) * Math.sqrt((r - a) / 2);
        return [re, im];
    }

    /**
     * Compute the atom permutation map for a given symmetry operation.
     * Returns an array `mapping` where `mapping[i] = j`, meaning atom `i`
     * maps to atom `j` under the operation `R * r_i + t`.
     * If the atom maps outside the visible supercell bounds, it returns undefined.
     */
    getAtomMapping(op) {
        if (!this.crystal.atompos || !this.crystal.atomobjects) return null;
        
        let mapping = [];
        let R = op.R_cart;
        let t = op.t_cart;
        let center = this.crystal.geometricCenter;
        let pos = this.crystal.atompos;
        
        for (let i = 0; i < pos.length; i++) {
            let truePos = new THREE.Vector3().copy(pos[i]).add(center);
            
            // R * r_i + t
            let rx = R[0][0]*truePos.x + R[0][1]*truePos.y + R[0][2]*truePos.z + t[0];
            let ry = R[1][0]*truePos.x + R[1][1]*truePos.y + R[1][2]*truePos.z + t[1];
            let rz = R[2][0]*truePos.x + R[2][1]*truePos.y + R[2][2]*truePos.z + t[2];
            
            let match_j = undefined;
            let minDist = 1e9;
            
            for (let j = 0; j < pos.length; j++) {
                if (this.crystal.atomobjects[i].atom_number !== this.crystal.atomobjects[j].atom_number) continue;
                
                let targetPos = new THREE.Vector3().copy(pos[j]).add(center);
                let dx = rx - targetPos.x;
                let dy = ry - targetPos.y;
                let dz = rz - targetPos.z;
                
                let dist = dx*dx + dy*dy + dz*dz;
                if (dist < minDist) {
                    minDist = dist;
                    match_j = j;
                }
            }
            
            // If the atom lands outside the drawn supercell bounds, minDist will be large.
            // We use a tolerance of 1e-2 Angstroms squared.
            if (minDist < 1e-2) {
                mapping[i] = match_j;
            } else {
                mapping[i] = undefined;
            }
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
        
        // SYMMETRY ADAPTATION
        // Compute D-matrix adapted vibrations for the target positions.
        this.currentAdaptedComponents = this.getAdaptedVibrationComponents(op);

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
                let j = (this.currentMapping && this.currentMapping[i] !== undefined) ? this.currentMapping[i] : undefined;
                
                // Use ADAPTED vibration components so degenerate bands show cleanly
                let vibrations_j = j !== undefined ? this.currentAdaptedComponents[j] : null;
                let vx_j = 0, vy_j = 0, vz_j = 0;
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
