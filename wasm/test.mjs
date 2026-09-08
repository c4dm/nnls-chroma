// Node smoke test for the WebAssembly build: synthesizes C, F, G, Am (two
// seconds each) and checks the labels come back in order.
//   node wasm/test.mjs            (run from the build's wasm/ directory)
import { loadChordino } from './chordino-wasm.js';

const rate = 44100;
const prog = [
  ['C', [36, 48, 52, 55, 60]],
  ['F', [41, 53, 57, 60, 65]],
  ['G', [43, 55, 59, 62, 67]],
  ['Am', [45, 57, 60, 64, 69]],
];
const perChord = rate * 2;
const mono = new Float32Array(perChord * prog.length);
prog.forEach(([, notes], k) => {
  for (const m of notes) {
    const f0 = 440 * Math.pow(2, (m - 69) / 12);
    for (let i = 0; i < perChord; i++) {
      const t = i / rate;
      const env = Math.exp(-t * 0.8) * (1 - Math.exp(-t * 200));
      let s = 0;
      for (let h = 1; h <= 6; h++) {
        if (f0 * h > rate * 0.45) break;
        s += Math.sin(2 * Math.PI * f0 * h * t) * Math.pow(0.6, h - 1);
      }
      mono[k * perChord + i] += s * env * 0.12;
    }
  }
});

const chordino = await loadChordino();
const t0 = Date.now();
const r = chordino.analyze(mono, rate);
console.log(`analysed ${mono.length} samples in ${Date.now() - t0} ms`);
for (const c of r.chords) console.log(c.time.toFixed(3).padStart(8), c.label);
const labels = r.chords.map((c) => c.label).join(' ');
if (!/^N C F G Am N$/.test(labels)) {
  console.log('MISMATCH:', labels);
  process.exit(1);
}
const harte = chordino.analyze(mono, rate, { harte: true });
if (!harte.chords.some((c) => c.label === 'A:min')) {
  console.log('MISMATCH (harte):', harte.chords.map((c) => c.label).join(' '));
  process.exit(1);
}
console.log('wasm OK');
