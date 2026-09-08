/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
  NNLS-Chroma / Chordino

  This file copyright 2026 Keith Adler.  GPL-2.0-or-later, see COPYING.
*/

/*
  Minimal audio file reading for the chordino tool: WAV, AIFF, FLAC,
  MP3 and Ogg Vorbis via the single-header decoders in third_party/.
*/

#ifndef CHORDINO_AUDIOFILE_H
#define CHORDINO_AUDIOFILE_H

#include <string>
#include <vector>

namespace chordino {

struct AudioData {
    float sampleRate;
    int channels;
    std::vector<float> mono;   // mixed down to one channel
};

// Decode a file, mixing to mono. Returns false and sets error on failure.
bool readAudioFile(const std::string &path, AudioData &out, std::string *error);

// Write a 16-bit PCM mono WAV (used by the tests and demos).
bool writeWavMono(const std::string &path, const float *samples, size_t frames,
                  float sampleRate, std::string *error);

} // namespace chordino

#endif
