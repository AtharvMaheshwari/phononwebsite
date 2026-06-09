const fs = require('fs');
const zlib = require('zlib');
let buffer = fs.readFileSync('data/localdb/mp-3427.json.gz');
let unzipped = zlib.gunzipSync(buffer);
let data = JSON.parse(unzipped.toString('utf8'));
console.log(Object.keys(data));
console.log("eigenvalues shape:", data.eigenvalues.length, data.eigenvalues[0].length);
