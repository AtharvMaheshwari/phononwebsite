const fs = require('fs');

let filepath = 'data/phonondb2017/mp-3427.json.gz';
let buffer = fs.readFileSync(filepath);
let data;
try {
    const zlib = require('zlib');
    let unzipped = zlib.gunzipSync(buffer);
    data = JSON.parse(unzipped.toString('utf8'));
} catch (e) {
    data = JSON.parse(buffer.toString('utf8'));
}

// data.eigenvalues is shape [nqpoints][nbands]
let nq = data.eigenvalues.length;
let nb = data.eigenvalues[0].length;

let am = [];
let c = [];
let h = [];

for(let k=0; k<nq; k++) {
    am.push([]);
    c.push([]);
    h.push([]);
    for(let b=0; b<nb; b++) {
        // dummy values
        let val_x = Math.sin(k * 0.1 + b);
        let val_y = Math.cos(k * 0.1 + b);
        let val_z = Math.sin(k * 0.2 - b);
        am[k].push([val_x, val_y, val_z]);
        
        c[k].push([val_x * 0.5, val_y * -0.5, Math.cos(k * 0.1) * 0.5, 1.0]);
        h[k].push(Math.sin(b) * Math.cos(k));
    }
}

data.angular_momentum = am;
data.cycloidicity = c;
data.helicity = h;

let newData = Buffer.from(JSON.stringify(data), 'utf8');

fs.writeFileSync(filepath, newData);
fs.writeFileSync('build/' + filepath, newData);
console.log('Done adding mock properties to Al4Li4O8 (mp-3427).');
