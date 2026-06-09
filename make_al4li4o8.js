const fs = require('fs');
const zlib = require('zlib');

let buffer = fs.readFileSync('data/localdb/mp-3427.json.gz');
let unzipped = zlib.gunzipSync(buffer);
let data = JSON.parse(unzipped.toString('utf8'));

let nq = data.qpoints.length;
let nb = data.natoms * 3;

let am_x = []; let am_y = []; let am_z = [];
let cx = []; let cy = []; let cz = [];
let h = [];

for(let k=0; k<nq; k++) {
    am_x.push([]); am_y.push([]); am_z.push([]);
    cx.push([]); cy.push([]); cz.push([]);
    h.push([]);
    for(let b=0; b<nb; b++) {
        let val = Math.sin(k * 0.1 + b);
        am_x[k].push(val);
        am_y[k].push(Math.cos(k * 0.1 + b));
        am_z[k].push(Math.sin(k * 0.2 - b));
        cx[k].push(val * 0.5);
        cy[k].push(val * -0.5);
        cz[k].push(Math.cos(k * 0.1) * 0.5);
        h[k].push(Math.sin(b) * Math.cos(k));
    }
}

data.angular_momentum_x = am_x;
data.angular_momentum_y = am_y;
data.angular_momentum_z = am_z;
data.cycloidicity_x = cx;
data.cycloidicity_y = cy;
data.cycloidicity_z = cz;
data.helicity = h;

fs.mkdirSync('data/localdb/al4li4o8', { recursive: true });
fs.writeFileSync('data/localdb/al4li4o8/data.json', JSON.stringify(data));

fs.mkdirSync('build/data/localdb/al4li4o8', { recursive: true });
fs.writeFileSync('build/data/localdb/al4li4o8/data.json', JSON.stringify(data));

let modelsPath = 'data/localdb/models.json';
let models = JSON.parse(fs.readFileSync(modelsPath, 'utf8'));
models.push({
    "folder": "data/localdb/al4li4o8",
    "name": "Al4Li4O8 (Mock Properties)"
});
fs.writeFileSync(modelsPath, JSON.stringify(models, null, 4));

let buildModelsPath = 'build/data/localdb/models.json';
fs.writeFileSync(buildModelsPath, JSON.stringify(models, null, 4));

console.log("Successfully generated Al4Li4O8 data with mock properties!");
