// Small wrapper around the Emscripten module produced by the CMake
// build under emcmake.  Usage (browser or Node):
//
//   import { loadChordino } from './chordino-wasm.js';
//   const chordino = await loadChordino();       // loads chordino.js + .wasm
//   const result = chordino.analyze(float32Samples, sampleRate, { harte: false });
//   // result.chords = [{time, duration, label}, ...]
//
// The samples must be mono; mix stereo down first (average the channels).

export async function loadChordino(options = {}) {
  const factory = (await import('./chordino.js')).default;
  const Module = await factory(options);
  const analyzeJson = Module.cwrap('chordino_analyze_json', 'number',
    ['number', 'number', 'number', 'number', 'number']);
  const freeResult = Module.cwrap('chordino_free', null, ['number']);

  return {
    module: Module,
    analyze(samples, sampleRate, opts = {}) {
      const n = samples.length;
      const ptr = Module._malloc(n * 4);
      try {
        Module.HEAPF32.set(samples, ptr >> 2);
        let flags = 0;
        if (opts.nnls === false) flags |= 1;
        if (opts.localTuning) flags |= 2;
        if (opts.harte) flags |= 4;
        if (opts.harmonicChange) flags |= 8;
        const boostN = typeof opts.boostN === 'number' ? opts.boostN : -1;
        const out = analyzeJson(ptr, n, sampleRate, flags, boostN);
        try {
          return JSON.parse(Module.UTF8ToString(out));
        } finally {
          freeResult(out);
        }
      } finally {
        Module._free(ptr);
      }
    },
  };
}
