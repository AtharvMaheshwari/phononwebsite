const fs = require('fs');
let data = JSON.parse(fs.readFileSync('data/localdb/al4li4o8/data.json'));

let kpoints = data.qpoints;
let nqpoints = kpoints.length;
// We don't have eigenvectors in data.json, so let's mock the structure just to test array allocation
let nbands = 48;
let natoms = data.natoms;

let vec = [];
for (let iq=0; iq<nqpoints; iq++) {
    let qvec = [];
    for (let nb=0; nb<nbands; nb++) {
        let bvec = [];
        for (let na=0; na<natoms; na++) {
            bvec.push([ [1,2], [3,4], [5,6] ]);
        }
        qvec.push(bvec);
    }
    vec.push(qvec);
}

let Jx = new Array(nqpoints).fill(0).map(() => new Array(nbands).fill(0));

for (let iq=0; iq<nqpoints; iq++) {
    let k_x = kpoints[iq][0], k_y = kpoints[iq][1], k_z = kpoints[iq][2];
    let norm_k = Math.sqrt(k_x*k_x + k_y*k_y + k_z*k_z);
    if (norm_k > 1e-8) {
        k_x /= norm_k; k_y /= norm_k; k_z /= norm_k;
    } else {
        k_x = 0; k_y = 0; k_z = 0;
    }

    for (let ibnd=0; ibnd<nbands; ibnd++) {
        let lx = 0, ly = 0, lz = 0;
        for (let iat=0; iat<natoms; iat++) {
            let v = vec[iq][ibnd][iat];
            let rx = v[0][0], ix = v[0][1];
            let ry = v[1][0], iy = v[1][1];
            let rz = v[2][0], iz = v[2][1];

            lx += 2 * (ry * iz - iy * rz);
        }
        Jx[iq][ibnd] = lx;
    }
}
console.log("No crash!");
