
export class PhononHighcharts {

    constructor(container) {
        this.container = container;
        this.phonon = { highsym_qpts: [] };
        this.highcharts = [];
        this.showModeWeights = false;
        this.atomTypeLegend = [];
        this.getAtomColorHex = null;
        this.getAtomLabel = null;
        this.weightLineWidthMin = 1.5;
        this.weightLineWidthScale = 8.0;
        this.legendVisibility = {};
        this.currentOptions = {};
        this.selectedPoint = null;
        this.selectedBandIndex = null;
        this.selectedX = null;

        let phonon = this.phonon;

        this.labels_formatter = function(phonon) {
            return function() {
                if ( phonon.highsym_qpts[this.value] ) {
                    let label = phonon.highsym_qpts[this.value];
                    label = label.replace("$","").replace("$","");
                    label = label.replace("\\Gamma","Γ");
                    label = label.replace("GAMMA","Γ");
                    label = label.replace("DELTA","Δ");
                    label = label.replace("\\Sigma","Σ");
                    label = label.replace("_","");
                    return label;
                }
                return ''
            }
        }

        this.top_labels_formatter = function(phonon) {
            return function() {
                let dist = this.value;
                // Check if this tick is at a high-symmetry point
                if (phonon.highsym_point_group_map && phonon.highsym_point_group_map[dist]) {
                    let pg = phonon.highsym_point_group_map[dist];
                    return '<span style="color:#0066cc; font-size:11px; font-weight:600;">' + pg + '</span>';
                }
                // Check if this tick is a segment midpoint
                if (phonon.segment_point_group_list) {
                    let tol = 1e-6;
                    for (let i = 0; i < phonon.segment_point_group_list.length; i++) {
                        let seg = phonon.segment_point_group_list[i];
                        let mid = (seg.start + seg.end) / 2.0;
                        if (Math.abs(dist - mid) < tol) {
                            let pg = seg.point_group;
                            return '<span style="color:#888888; font-size:10px; font-style:italic;">' + pg + '</span>';
                        }
                    }
                }
                return '';
            }
        }

        this.HighchartsOptions = {
            chart: { type: 'line',
                     zoomType: 'xy' },
            accessibility: { enabled: false },
            title: { text: null },
            xAxis: [
                {
                    plotLines: [],
                    lineWidth: 0,
                    minorGridLineWidth: 0,
                    lineColor: 'transparent',
                    minorTickLength: 0,
                    tickLength: 0,
                    crosshair: { width: 1, color: '#999999', dashStyle: 'solid' },
                    labels: {
                        style: { fontSize:'20px' },
                        allowOverlap: true
                    }
                },
                {
                    linkedTo: 0,
                    opposite: true,
                    lineWidth: 0,
                    minorGridLineWidth: 0,
                    lineColor: 'transparent',
                    minorTickLength: 0,
                    tickLength: 0,
                    labels: {
                        useHTML: true,
                        allowOverlap: true,
                        formatter: this.top_labels_formatter(phonon)
                    }
                }
            ],
            yAxis: { title: { text: 'Frequency (cm<sup>-1</sup>)' },
                     crosshair: { width: 1, color: '#999999', dashStyle: 'solid' },
                     plotLines: [ {value: 0, color: '#000000', width: 2} ]
                   },
            tooltip: { 
                enabled: false
            },
            legend: {
                enabled: false,
                floating: true,
                layout: 'horizontal',
                align: 'center',
                verticalAlign: 'top',
                y: 8,
                backgroundColor: 'rgba(255,255,255,0.85)',
                borderWidth: 0,
                itemStyle: { fontWeight: 'normal' },
                symbolRadius: 0
            },
            series: [],
            plotOptions: { line:   { animation: false },
                           series: { 
                                findNearestPointBy: 'xy',
                                turboThreshold: 0,
                                allowPointSelect: true,
                                stickyTracking: false,      // ADD: prevents hover state bleeding across series
                                marker: { 
                                    states: { 
                                        select: { 
                                            enabled: true,
                                            fillColor: '#000000',
                                            radius: 4,
                                            lineWidth: 0,
                                            symbol: 'circle'
                                        },
                                        hover: { enabled: true, radius: 4 }
                                    } 
                                },
                                states: {
                                    hover: {
                                        halo: false,
                                        lineWidthPlus: 0
                                    },
                                    inactive: { opacity: 1 }   // ADD: prevents non-hovered series from fading
                                },
                                cursor: 'pointer',
                                point: { events: { } }
                            }
                         }
        };
        this.HighchartsOptions.__phononHighcharts = this;
    }

    setClickEvent(phononweb) {
        let click_event = function () {
            if (this.series.options.isLegendSeries || this.series.options.isWeightSeries) { return; }
            let k = Number.isFinite(this.options && this.options.kIndex)
                ? this.options.kIndex
                : phononweb.phonon.qindex[this.x];
            let n = this.series.options.bandIndex;
            if (!Number.isFinite(n)) {
                n = Number(this.series.name);
            }
            
            // If the user clicks on the already selected point, 
            // cycle through other bands that overlap at this frequency
            if (phononweb.k === k && phononweb.n === n) {
                let freqs = phononweb.phonon.eigenvalues[k];
                let overlapping_bands = [];
                for (let b = 0; b < freqs.length; b++) {
                    // Check if frequency is within 3 cm^-1
                    if (Math.abs(freqs[b] - this.y) < 3.0) {
                        overlapping_bands.push(b);
                    }
                }
                if (overlapping_bands.length > 1) {
                    let idx = overlapping_bands.indexOf(n);
                    if (idx !== -1) {
                        n = overlapping_bands[(idx + 1) % overlapping_bands.length];
                    }
                }
            }
            
            phononweb.selectModeByBandIndex(k, n, true);
            return false;
        };

        // ADD: re-assert selected point appearance on any point hover
        let mouseOver_event = function() {
            let phHC = this.series.chart.userOptions.__phononHighcharts;
            if (!phHC || !phHC.selectedPoint || phHC.selectedPoint === this) { return; }
            // Defer so Highcharts finishes its own hover state logic first
            setTimeout(() => {
                if (phHC.selectedPoint && typeof phHC.selectedPoint.select === 'function') {
                    try { phHC.selectedPoint.select(true, false); } catch(e) {}
                }
            }, 0);
        };

        this.HighchartsOptions.plotOptions.series.point.events.click = click_event;
        this.HighchartsOptions.plotOptions.series.point.events.mouseOver = mouseOver_event;
    }

    selectModePoint(phonon, k, n) {
        if (!this.chart || !this.chart.series || !phonon || !phonon.distances) { return; }
        let targetX = phonon.distances[k];
        if (!Number.isFinite(targetX)) { return; }

        if (this.selectedPoint &&
            this.selectedPoint.selected &&
            this.selectedBandIndex === n &&
            this.selectedK === k &&
            this.selectedX === targetX) {
            return;
        }

        // Remove patch from old selected point before deselecting
        if (this.selectedPoint) {
            if (this.selectedPoint.__origSetState) {
                this.selectedPoint.setState = this.selectedPoint.__origSetState;
                delete this.selectedPoint.__origSetState;
            }
            try { this.selectedPoint.select(false, false); } catch(e) {}
            this.selectedPoint = null;
        }

        let bestPoint = null;
        let bestDist = Infinity;

        for (let i = 0; i < this.chart.series.length; i++) {
            let series = this.chart.series[i];
            if (series.options.isLegendSeries || series.options.isWeightSeries) { continue; }
            if (Number(series.options.bandIndex) !== Number(n)) { continue; }

            for (let j = 0; j < series.points.length; j++) {
                let point = series.points[j];
                let pointK = point.options ? point.options.kIndex : undefined;
                if (Number.isFinite(pointK) && Number(pointK) === Number(k)) {
                    bestPoint = point;
                    bestDist = -1;
                    break;
                }
                let d = Math.abs(point.x - targetX);
                if (d < bestDist) {
                    bestDist = d;
                    bestPoint = point;
                }
            }
            if (bestDist === -1) { break; }
        }

        if (bestPoint) {
            bestPoint.select(true, false);

            // Patch setState so Highcharts can never visually reset this point
            bestPoint.__origSetState = bestPoint.setState;
            bestPoint.setState = function(state) {
                // Ignore any attempt to clear or change our select state
                if (state === '' || state === 'normal' || state === 'inactive') { return; }
                bestPoint.__origSetState.call(this, state);
            };

            this.selectedPoint = bestPoint;
            this.selectedBandIndex = n;
            this.selectedK = k;
            this.selectedX = targetX;
        } else {
            this.selectedBandIndex = null;
            this.selectedK = null;
            this.selectedX = null;
        }
    }

    setModeWeightsOptions(options = {}) {
        this.currentOptions = options;
        this.showModeWeights = !!options.enabled;
        this.getAtomColorHex = typeof options.getAtomColorHex === 'function' ? options.getAtomColorHex : null;
        this.getAtomLabel = typeof options.getAtomLabel === 'function' ? options.getAtomLabel : null;
        this.heatmapProperty = options.heatmapProperty || "none";
    }

    getHeatmapValue(phonon, property, k, n) {
        if (!phonon) return null;
        let arr = phonon[property];

        if (arr && arr[k] && arr[k][n] !== undefined) {
            return arr[k][n];
        }
        return null;
    }

    getHeatmapRange(phonon, property) {
        let min = Infinity, max = -Infinity;
        if (!phonon || !phonon.eigenvalues) return {min: -1, max: 1};
        let nbands = phonon.eigenvalues[0].length;
        for (let k = 0; k < phonon.kpoints.length; k++) {
            for (let n = 0; n < nbands; n++) {
                let val = this.getHeatmapValue(phonon, property, k, n);
                if (val !== null) {
                    if (val < min) min = val;
                    if (val > max) max = val;
                }
            }
        }
        if (min === Infinity) return {min: -1, max: 1};
        
        if (property && property.includes('pam')) {
            if (max === min) max = min + 1;
            return {min: 0, max: max};
        }
        
        // Zero-center the heatmap
        let maxAbs = Math.max(Math.abs(min), Math.abs(max));
        if (maxAbs === 0) maxAbs = 1;
        return {min: -maxAbs, max: maxAbs};
    }

    valueToColor(value, min, max, propertyName = null) {
        if (!Number.isFinite(value)) return '#888888';
        
        if (propertyName && propertyName.includes('pam')) {
            let pamColors = [
                [144, 144, 144], // 0.0: Gray
                [207, 135, 72],  // 0.5: Soft Orange
                [198, 171, 36],  // 1.0: Mustard Yellow
                [90, 171, 90],   // 1.5: Soft Green
                [54, 117, 207],  // 2.0: Soft Blue
                [135, 81, 180],  // 2.5: Soft Purple
                [198, 63, 63],   // 3.0: Soft Red
                [135, 108, 135], // 3.5: Soft Mauve
                [36, 162, 162],  // 4.0: Soft Cyan
                [90, 135, 180],  // 4.5: Soft Blue
                [180, 72, 180],  // 5.0: Soft Magenta
                [117, 108, 198], // 5.5: Soft Periwinkle
                [36, 162, 162],  // 6.0: Soft Cyan
            ];
            
            let n_fold = Math.max(1, Math.round(max));
            let v = Math.abs(value) % n_fold;
            if (v < 0) v += n_fold;
            
            let scaled_v = v * 2.0;
            let idx1 = Math.floor(scaled_v);
            let idx2 = idx1 + 1;
            let frac = scaled_v - idx1;
            
            let numSteps = n_fold * 2;
            let c1 = pamColors[idx1 % pamColors.length];
            let c2 = pamColors[(idx1 + 1 === numSteps) ? 0 : (idx2 % pamColors.length)];
            
            let r = Math.round(c1[0] + (c2[0] - c1[0]) * frac);
            let g = Math.round(c1[1] + (c2[1] - c1[1]) * frac);
            let b = Math.round(c1[2] + (c2[2] - c1[2]) * frac);
            
            return '#' + r.toString(16).padStart(2, '0') + g.toString(16).padStart(2, '0') + b.toString(16).padStart(2, '0');
        }

        let norm = (max === min) ? 0.5 : (value - min) / (max - min);
        norm = Math.max(0, Math.min(1, norm));
        
        let r, g, b;
        if (norm < 0.5) {
            let t = norm * 2.0; 
            r = Math.round(t * 200);
            g = Math.round(t * 200);
            b = Math.round(255 - t * 55);
        } else {
            let t = (norm - 0.5) * 2.0; 
            r = Math.round(200 + t * 55);
            g = Math.round((1 - t) * 200);
            b = Math.round((1 - t) * 200);
        }
        return '#' + r.toString(16).padStart(2, '0') + g.toString(16).padStart(2, '0') + b.toString(16).padStart(2, '0');
    }

    syncLegendVisibility(reset = false) {
        let nextVisibility = {};
        for (let i = 0; i < this.atomTypeLegend.length; i++) {
            let atomNumber = this.atomTypeLegend[i].atomNumber;
            nextVisibility[atomNumber] = reset
                ? true
                : this.legendVisibility[atomNumber] !== false;
        }
        this.legendVisibility = nextVisibility;
    }

    isAtomNumberVisible(atomNumber) {
        return this.legendVisibility[atomNumber] !== false;
    }

    refreshLegendAndWeights() {
        this.updateWeightedSeriesStyles();
        this.applyLegendStyles();
    }

    ensureAtomTypeWeights(phonon) {
        let drawHeatmap = this.heatmapProperty && this.heatmapProperty !== 'none';
        if ((this.showModeWeights || drawHeatmap) && phonon && typeof phonon.ensureAllEigenvectors === 'function') {
            phonon.ensureAllEigenvectors();
        }
        if (!phonon || !phonon.vec || !phonon.atom_numbers) {
            return { atomNumbers: [], weights: [] };
        }
        if (phonon.atomTypeWeightsCache) {
            return phonon.atomTypeWeightsCache;
        }

        let uniqueAtomNumbers = [];
        let atomTypeIndex = {};
        for (let i = 0; i < phonon.atom_numbers.length; i++) {
            let atomNumber = phonon.atom_numbers[i];
            if (!(atomNumber in atomTypeIndex)) {
                atomTypeIndex[atomNumber] = uniqueAtomNumbers.length;
                uniqueAtomNumbers.push(atomNumber);
            }
        }

        let weights = [];
        for (let k = 0; k < phonon.vec.length; k++) {
            if (!phonon.vec[k] || !phonon.vec[k].length) {
                weights.push([]);
                continue;
            }
            let qpointWeights = [];
            for (let n = 0; n < phonon.vec[k].length; n++) {
                let mode = phonon.vec[k][n];
                let totals = new Array(uniqueAtomNumbers.length).fill(0);
                let totalWeight = 0;

                for (let atomIndex = 0; atomIndex < mode.length; atomIndex++) {
                    let components = mode[atomIndex];
                    let atomWeight = 0;
                    for (let axis = 0; axis < components.length; axis++) {
                        let component = components[axis];
                        let re = component[0];
                        let im = component[1];
                        atomWeight += re * re + im * im;
                    }
                    let typeIndex = atomTypeIndex[phonon.atom_numbers[atomIndex]];
                    totals[typeIndex] += atomWeight;
                    totalWeight += atomWeight;
                }

                if (totalWeight > 0) {
                    for (let i = 0; i < totals.length; i++) {
                        totals[i] /= totalWeight;
                    }
                }
                qpointWeights.push(totals);
            }
            weights.push(qpointWeights);
        }

        phonon.atomTypeWeightsCache = {
            atomNumbers: uniqueAtomNumbers,
            weights: weights
        };
        return phonon.atomTypeWeightsCache;
    }

    getAtomTypeLegend(phonon) {
        let cache = this.ensureAtomTypeWeights(phonon);
        return cache.atomNumbers.map((atomNumber) => ({
            atomNumber: atomNumber,
            label: this.getAtomLabel ? this.getAtomLabel(atomNumber) : String(atomNumber),
            color: this.getAtomColorHex ? ('#' + Number(this.getAtomColorHex(atomNumber)).toString(16).padStart(6, '0')) : '#0066FF'
        }));
    }

    reflow() {
        if (this.chart && this.chart.reflow) {
            this.chart.reflow();
        }
    }

    handleLegendToggle(atomNumber) {
        if (!this.phonon || !this.chart) {
            return false;
        }

        if (!(atomNumber in this.legendVisibility)) {
            return false;
        }

        this.legendVisibility[atomNumber] = !this.isAtomNumberVisible(atomNumber);
        this.refreshLegendAndWeights();
        return false;
    }

    applyLegendStyles() {
        if (!this.chart || !this.chart.legend || !this.chart.legend.allItems) {
            return;
        }

        for (let i = 0; i < this.chart.legend.allItems.length; i++) {
            let item = this.chart.legend.allItems[i];
            if (!item || !item.options || !item.options.isLegendSeries || !item.legendItem) {
                continue;
            }

            let atomNumber = item.options.atomNumber;
            let visible = this.isAtomNumberVisible(atomNumber);
            let textDecoration = visible ? 'none' : 'line-through';
            let opacity = visible ? 1 : 0.65;
            let legendItem = item.legendItem;

            if (legendItem && typeof legendItem.css === 'function') {
                legendItem.css({
                    textDecoration: textDecoration,
                    opacity: opacity
                });
            } else if (legendItem && typeof legendItem.attr === 'function') {
                legendItem.attr({
                    opacity: opacity
                });
                if (legendItem.element && legendItem.element.style) {
                    legendItem.element.style.textDecoration = textDecoration;
                }
            } else if (legendItem && legendItem.style) {
                legendItem.style.textDecoration = textDecoration;
                legendItem.style.opacity = String(opacity);
            }

            if (item.legendGroup && typeof item.legendGroup.attr === 'function') {
                item.legendGroup.attr({ opacity: opacity });
            }
        }
    }

    getVisibleAtomTypeIndices() {
        let visible = [];
        for (let i = 0; i < this.atomTypeLegend.length; i++) {
            let atomNumber = this.atomTypeLegend[i].atomNumber;
            if (this.isAtomNumberVisible(atomNumber)) {
                visible.push(i);
            }
        }
        return visible;
    }

    getWeightedColorForIndices(atomIndices, weights) {
        let total = 0;
        let red = 0;
        let green = 0;
        let blue = 0;

        for (let i = 0; i < atomIndices.length; i++) {
            let index = atomIndices[i];
            let weight = weights[index];
            let color = this.getAtomColorHex
                ? this.getAtomColorHex(this.atomTypeLegend[index].atomNumber)
                : 0x0066ff;
            let r = (color >> 16) & 255;
            let g = (color >> 8) & 255;
            let b = color & 255;
            red += r * weight;
            green += g * weight;
            blue += b * weight;
            total += weight;
        }

        if (total <= 0) {
            let fallback = this.getAtomColorHex
                ? this.getAtomColorHex(this.atomTypeLegend[atomIndices[0]].atomNumber)
                : 0x0066ff;
            return '#' + Number(fallback).toString(16).padStart(6, '0');
        }

        red = Math.round(red / total);
        green = Math.round(green / total);
        blue = Math.round(blue / total);
        return '#' + ((red << 16) | (green << 8) | blue).toString(16).padStart(6, '0');
    }

    updateWeightedSeriesStyles() {
        if (!this.chart || !this.phonon || !this.showModeWeights) {
            return;
        }

        let visibleTypeIndices = this.getVisibleAtomTypeIndices();
        let singleVisibleType = visibleTypeIndices.length === 1 ? visibleTypeIndices[0] : null;

        for (let i = 0; i < this.chart.series.length; i++) {
            let series = this.chart.series[i];
            if (!series.options.isWeightSeries) {
                continue;
            }

            let avgWeights = series.options.avgWeights;
            if (!avgWeights) {
                continue;
            }

            let visible = visibleTypeIndices.length > 0;
            let color = '#0066ff';
            let lineWidth = this.weightLineWidthMin + this.weightLineWidthScale;

            if (visible) {
                color = this.getWeightedColorForIndices(visibleTypeIndices, avgWeights);

                if (singleVisibleType !== null) {
                    lineWidth = this.weightLineWidthMin + avgWeights[singleVisibleType] * this.weightLineWidthScale;
                    color = '#' + Number(this.getAtomColorHex(this.atomTypeLegend[singleVisibleType].atomNumber)).toString(16).padStart(6, '0');
                }
            }

            series.visible = visible;
            series.options.visible = visible;
            series.color = color;
            series.options.color = color;
            series.options.lineWidth = lineWidth;

            if (series.group) {
                if (visible) {
                    series.group.show();
                } else {
                    series.group.hide();
                }
            }

            if (series.graph) {
                series.graph.attr({
                    stroke: color,
                    'stroke-width': lineWidth
                });
            }
        }
    }

    refreshWeightedSeriesColors() {
        if (!this.chart || !this.phonon || !this.showModeWeights) {
            return;
        }

        let visibleTypeIndices = this.getVisibleAtomTypeIndices();
        let singleVisibleType = visibleTypeIndices.length === 1 ? visibleTypeIndices[0] : null;

        for (let i = 0; i < this.chart.series.length; i++) {
            let series = this.chart.series[i];
            if (!series.options.isWeightSeries) {
                continue;
            }

            let avgWeights = series.options.avgWeights;
            if (!avgWeights) {
                continue;
            }

            let color = this.getWeightedColorForIndices(visibleTypeIndices, avgWeights);
            if (singleVisibleType !== null) {
                color = '#' + Number(this.getAtomColorHex(this.atomTypeLegend[singleVisibleType].atomNumber)).toString(16).padStart(6, '0');
            }

            series.color = color;
            series.options.color = color;
            if (series.graph) {
                series.graph.attr({ stroke: color });
            }
        }
    }

    refreshAppearance(options = {}) {
        if (!this.chart || !this.phonon) {
            return;
        }

        let previousShowModeWeights = this.showModeWeights;
        let previousHeatmapProperty = this.heatmapProperty;
        this.setModeWeightsOptions(options);

        if (previousShowModeWeights !== this.showModeWeights || previousHeatmapProperty !== this.heatmapProperty) {
            this.update(this.phonon, options);
            return;
        }

        if (!this.showModeWeights) {
            return;
        }

        this.atomTypeLegend = this.getAtomTypeLegend(this.phonon);

        for (let i = 0; i < this.chart.series.length; i++) {
            let series = this.chart.series[i];
            if (!series.options.isLegendSeries) {
                continue;
            }

            let atomNumber = series.options.atomNumber;
            let legendEntry = null;
            for (let j = 0; j < this.atomTypeLegend.length; j++) {
                if (this.atomTypeLegend[j].atomNumber === atomNumber) {
                    legendEntry = this.atomTypeLegend[j];
                    break;
                }
            }
            if (!legendEntry) {
                continue;
            }

            series.update({
                name: legendEntry.label,
                color: legendEntry.color
            }, false);
        }

        this.refreshWeightedSeriesColors();
        this.applyLegendStyles();
        this.chart.redraw(false);
    }

    update(phonon, options = {}) {
        /*
        update phonon dispersion plot
        */

        this.phonon = phonon;
        this.setModeWeightsOptions(options);
        this.atomTypeLegend = this.getAtomTypeLegend(phonon);
        this.syncLegendVisibility(!!options.resetLegendVisibility);

        if (this.selectedPoint && this.selectedPoint.__origSetState) {
            this.selectedPoint.setState = this.selectedPoint.__origSetState;
            delete this.selectedPoint.__origSetState;
        }

        this.selectedPoint = null;
        this.selectedBandIndex = null;
        this.selectedK = null;
        this.selectedX = null;

    
        //set the minimum of the plot with the smallest phonon frequency
        let minVal = 0;
        for (let i=0; i<phonon.eigenvalues.length; i++) {
            let min = Math.min.apply(null, phonon.eigenvalues[i])
            if ( minVal > min ) {
                minVal = min;
            }
        }
        if (minVal > -1) minVal = 0;


        //get positions of high symmetry qpoints
        let ticks = Object.keys(phonon.highsym_qpts)
            .map((k) => Number(k))
            .filter((k) => Number.isFinite(k))
            .sort((a, b) => a - b);

        //get the high symmetry qpoints for highcharts
        let plotLines = []
        for (let i=0; i<ticks.length ; i++ ) {
            plotLines.push({ value: ticks[i],
                             color: '#555555',
                             width: 1,
                             zIndex: 6 })
        }

        //actually set the eigenvalues
        this.getGraph(phonon);

        this.HighchartsOptions.series = this.highcharts;
        this.HighchartsOptions.legend.enabled = this.showModeWeights && this.atomTypeLegend.length > 0;
        //calculate midpoints for line symmetry labels
        let topTicks = [...ticks];
        for (let i = 0; i < ticks.length - 1; i++) {
            topTicks.push((ticks[i] + ticks[i+1]) / 2.0);
        }
        topTicks.sort((a, b) => a - b);

        this.HighchartsOptions.xAxis[0].tickPositions = ticks;
        this.HighchartsOptions.xAxis[0].plotLines = plotLines;
        this.HighchartsOptions.xAxis[0].labels.formatter = this.labels_formatter(phonon);
        this.HighchartsOptions.xAxis[1].tickPositions = topTicks;
        this.HighchartsOptions.xAxis[1].labels.formatter = this.top_labels_formatter(phonon);
        this.HighchartsOptions.yAxis.min = minVal;
        if (this.chart) {
            this.chart.destroy();
            this.chart = null;
        }
        this.chart = globalThis.Highcharts.chart(this.container[0], this.HighchartsOptions);
        this.refreshLegendAndWeights();
    }

    isGammaLabel(label) {
        if (typeof label !== 'string') {
            return false;
        }
        let normalized = label.replace(/\$/g, '').replace(/\s+/g, '').toUpperCase();
        return normalized === 'G' || normalized === 'GAMMA' || normalized === '\\GAMMA' || normalized === 'Γ';
    }

    hasGammaDiscontinuity(phonon, previousIndex, nextIndex) {
        if (!phonon || !phonon.kpoints || !phonon.distances) {
            return false;
        }

        let prevQ = phonon.kpoints[previousIndex];
        let nextQ = phonon.kpoints[nextIndex];
        if (!prevQ || !nextQ) {
            return false;
        }

        let sameQpoint = true;
        for (let axis = 0; axis < 3; axis++) {
            if (Math.abs(Number(prevQ[axis]) - Number(nextQ[axis])) > 1e-8) {
                sameQpoint = false;
                break;
            }
        }
        if (!sameQpoint) {
            return false;
        }

        let sameDistance = Math.abs(Number(phonon.distances[previousIndex]) - Number(phonon.distances[nextIndex])) < 1e-10;
        if (!sameDistance) {
            return false;
        }

        let label = phonon.highsym_qpts ? phonon.highsym_qpts[phonon.distances[nextIndex]] : null;
        return this.isGammaLabel(label);
    }

    getPlotSegments(phonon, startk, endk) {
        let segments = [];
        let segmentStart = startk;

        for (let k = startk + 1; k < endk; k++) {
            if (this.hasGammaDiscontinuity(phonon, k - 1, k)) {
                if (k - segmentStart > 0) {
                    segments.push([segmentStart, k]);
                }
                segmentStart = k;
            }
        }

        if (endk - segmentStart > 0) {
            segments.push([segmentStart, endk]);
        }

        return segments;
    }

    getGraph(phonon) {
        /*
        From a phonon object containing:
            distances : distance between the k-points
            eigenvalues : eigenvalues
        put the data in the highcharts format
        */

        let eival = phonon.eigenvalues;
        let dists = phonon.distances;
        let line_breaks = phonon.line_breaks;

        let nbands = eival[0].length;
        this.highcharts = [];
        let weightCache = this.ensureAtomTypeWeights(phonon);
        let drawHeatmap = this.heatmapProperty && this.heatmapProperty !== 'none';
        let baseColor = (this.showModeWeights || drawHeatmap) ? '#94a3b8' : '#0066FF';
        let visibleTypeIndices = this.getVisibleAtomTypeIndices();
        let singleVisibleType = visibleTypeIndices.length === 1 ? visibleTypeIndices[0] : null;
        
        if (drawHeatmap) {
            this.heatmapRange = this.getHeatmapRange(phonon, this.heatmapProperty);
            if (typeof this.updateColorbar === 'function') {
                this.updateColorbar(this.heatmapRange.min, this.heatmapRange.max, this.heatmapProperty);
            }
        }

        let bucketCache = {};

        //go through the eigenvalues and create eival list
        for (let n=0; n<nbands; n++) {
            //iterate over the line breaks
            for (let i=0; i<line_breaks.length; i++) {
                let startk = line_breaks[i][0];
                let endk = line_breaks[i][1];
                let plotSegments = this.getPlotSegments(phonon, startk, endk);

                for (let segmentIndex = 0; segmentIndex < plotSegments.length; segmentIndex++) {
                    let segmentStart = plotSegments[segmentIndex][0];
                    let segmentEnd = plotSegments[segmentIndex][1];
                    let data = [];
                    for (let k=segmentStart; k<segmentEnd; k++) {
                        let propVal = null;
                        if (drawHeatmap) {
                            propVal = this.getHeatmapValue(phonon, this.heatmapProperty, k, n);
                        }
                        data.push({
                            x: dists[k],
                            y: eival[k][n],
                            kIndex: k,
                            propValue: propVal
                        });
                    }

                    this.highcharts.push({
                        name:  n+"",
                        bandIndex: n,
                        color: drawHeatmap ? 'transparent' : baseColor,
                        lineWidth: this.showModeWeights ? 0.8 : 2,
                        zIndex: 5,
                        showInLegend: false,
                        marker: { radius: 1, symbol: "circle"},
                        data: data
                    });

                    if (this.showModeWeights && !drawHeatmap) {
                        for (let k=segmentStart; k<segmentEnd - 1; k++) {
                            if (!visibleTypeIndices.length) {
                                continue;
                            }

                            let weights0 = weightCache.weights[k][n];
                            let weights1 = weightCache.weights[k + 1][n];
                            let avgWeights = weights0.map((value, index) => (value + weights1[index]) / 2);
                            let lineWidth = this.weightLineWidthMin + this.weightLineWidthScale;
                            let color = this.getWeightedColorForIndices(visibleTypeIndices, avgWeights);

                            if (singleVisibleType !== null) {
                                lineWidth = this.weightLineWidthMin + avgWeights[singleVisibleType] * this.weightLineWidthScale;
                                color = '#' + Number(this.getAtomColorHex(this.atomTypeLegend[singleVisibleType].atomNumber)).toString(16).padStart(6, '0');
                            }

                            let bucketKey = 'weights_' + color + '_' + lineWidth;
                            if (!bucketCache[bucketKey]) {
                                bucketCache[bucketKey] = { color: color, lineWidth: lineWidth, isHeatmapSeries: false, segments: [] };
                            }
                            bucketCache[bucketKey].segments.push([
                                [dists[k], eival[k][n]],
                                [dists[k + 1], eival[k + 1][n]]
                            ]);
                        }
                    } else if (drawHeatmap) {
                        for (let k=segmentStart; k<segmentEnd - 1; k++) {
                            let val0 = this.getHeatmapValue(phonon, this.heatmapProperty, k, n);
                            let val1 = this.getHeatmapValue(phonon, this.heatmapProperty, k + 1, n);
                            let color = '#aaaaaa';
                            if (val0 !== null && val1 !== null) {
                                let avgVal = (val0 + val1) / 2;
                                color = this.valueToColor(avgVal, this.heatmapRange.min, this.heatmapRange.max, this.heatmapProperty);
                            }
                            
                            let bucketKey = 'heatmap_' + color;
                            if (!bucketCache[bucketKey]) {
                                bucketCache[bucketKey] = { color: color, lineWidth: 4, isHeatmapSeries: true, segments: [] };
                            }
                            bucketCache[bucketKey].segments.push([
                                [dists[k], eival[k][n]],
                                [dists[k + 1], eival[k + 1][n]]
                            ]);
                        }
                    }
                }
            }
        }

        for (let key in bucketCache) {
            let bucket = bucketCache[key];
            
            // Highcharts requires data to be sorted by X for binary search zooming to work!
            bucket.segments.sort((a, b) => a[0][0] - b[0][0]);
            
            let data = [];
            for (let i = 0; i < bucket.segments.length; i++) {
                data.push(bucket.segments[i][0]);
                data.push(bucket.segments[i][1]);
                data.push([bucket.segments[i][1][0], null]);
            }
            
            this.highcharts.push({
                name: bucket.isHeatmapSeries ? 'heatmap' : 'weights',
                isHeatmapSeries: bucket.isHeatmapSeries,
                isWeightSeries: !bucket.isHeatmapSeries,
                color: bucket.color,
                lineWidth: bucket.lineWidth,
                zIndex: 3,
                enableMouseTracking: false,
                showInLegend: false,
                states: { inactive: { opacity: 1 } },
                marker: { enabled: false },
                data: data
            });
        }

        if (this.showModeWeights) {
            for (let i = 0; i < this.atomTypeLegend.length; i++) {
                let atomType = this.atomTypeLegend[i];
                this.highcharts.push({
                    id: 'legend-' + atomType.atomNumber,
                    name: atomType.label,
                    atomNumber: atomType.atomNumber,
                    isLegendSeries: true,
                    color: atomType.color,
                    data: [],
                    enableMouseTracking: false,
                    showInLegend: true,
                    lineWidth: 8,
                    marker: { enabled: false },
                    states: { inactive: { opacity: 1 } },
                    events: {
                        legendItemClick: function() {
                            return this.chart.userOptions.__phononHighcharts.handleLegendToggle(atomType.atomNumber);
                        }
                    }
                });
            }
        }
    }
}
