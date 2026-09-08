/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
  NNLS-Chroma / Chordino

  This file copyright 2026 Keith Adler.  GPL-2.0-or-later, see COPYING.
*/

#include "cli/audiofile.h"

#define DR_WAV_IMPLEMENTATION
#define DR_FLAC_IMPLEMENTATION
#define DR_MP3_IMPLEMENTATION
#include "dr_libs/dr_wav.h"
#include "dr_libs/dr_flac.h"
#include "dr_libs/dr_mp3.h"

#define STB_VORBIS_NO_PUSHDATA_API
#define STB_VORBIS_NO_INTEGER_CONVERSION
#include "stb/stb_vorbis.c"

#include <cstdio>
#include <cstring>
#include <cstdint>
#include <cctype>
#include <algorithm>

namespace chordino {

namespace {

std::string lowerExt(const std::string &path)
{
    size_t dot = path.rfind('.');
    size_t slash = path.find_last_of("/\\");
    if (dot == std::string::npos || (slash != std::string::npos && dot < slash)) return "";
    std::string e = path.substr(dot + 1);
    for (size_t i = 0; i < e.size(); ++i) e[i] = (char)tolower((unsigned char)e[i]);
    return e;
}

void mixDown(const float *interleaved, uint64_t frames, int channels, std::vector<float> &mono)
{
    mono.resize((size_t)frames);
    if (channels == 1) {
        std::copy(interleaved, interleaved + frames, mono.begin());
        return;
    }
    const float scale = 1.f / (float)channels;
    for (uint64_t i = 0; i < frames; ++i) {
        float s = 0.f;
        for (int c = 0; c < channels; ++c) s += interleaved[i * channels + c];
        mono[(size_t)i] = s * scale;
    }
}

bool readWav(const std::string &path, AudioData &out)
{
    unsigned int channels = 0, rate = 0;
    drwav_uint64 frames = 0;
    float *data = drwav_open_file_and_read_pcm_frames_f32(path.c_str(), &channels, &rate, &frames, NULL);
    if (!data) return false;
    out.sampleRate = (float)rate;
    out.channels = (int)channels;
    mixDown(data, frames, (int)channels, out.mono);
    drwav_free(data, NULL);
    return true;
}

bool readFlac(const std::string &path, AudioData &out)
{
    unsigned int channels = 0, rate = 0;
    drflac_uint64 frames = 0;
    float *data = drflac_open_file_and_read_pcm_frames_f32(path.c_str(), &channels, &rate, &frames, NULL);
    if (!data) return false;
    out.sampleRate = (float)rate;
    out.channels = (int)channels;
    mixDown(data, frames, (int)channels, out.mono);
    drflac_free(data, NULL);
    return true;
}

bool readMp3(const std::string &path, AudioData &out)
{
    drmp3_config config;
    drmp3_uint64 frames = 0;
    float *data = drmp3_open_file_and_read_pcm_frames_f32(path.c_str(), &config, &frames, NULL);
    if (!data) return false;
    out.sampleRate = (float)config.sampleRate;
    out.channels = (int)config.channels;
    mixDown(data, frames, (int)config.channels, out.mono);
    drmp3_free(data, NULL);
    return true;
}

bool readOgg(const std::string &path, AudioData &out)
{
    int err = 0;
    stb_vorbis *v = stb_vorbis_open_filename(path.c_str(), &err, NULL);
    if (!v) return false;
    stb_vorbis_info info = stb_vorbis_get_info(v);
    out.sampleRate = (float)info.sample_rate;
    out.channels = info.channels;
    out.mono.clear();
    std::vector<float> buf((size_t)info.channels * 4096);
    for (;;) {
        int got = stb_vorbis_get_samples_float_interleaved(v, info.channels, &buf[0], (int)buf.size());
        if (got <= 0) break;
        const float scale = 1.f / (float)info.channels;
        for (int i = 0; i < got; ++i) {
            float s = 0.f;
            for (int c = 0; c < info.channels; ++c) s += buf[(size_t)i * info.channels + c];
            out.mono.push_back(s * scale);
        }
    }
    stb_vorbis_close(v);
    return true;
}

// AIFF/AIFC (16/24/32-bit PCM, big-endian, uncompressed) reader; the
// dr_libs do not cover it and it is the other common Mac format.
uint32_t be32(const unsigned char *p) { return ((uint32_t)p[0] << 24) | ((uint32_t)p[1] << 16) | ((uint32_t)p[2] << 8) | p[3]; }
uint16_t be16(const unsigned char *p) { return (uint16_t)(((uint16_t)p[0] << 8) | p[1]); }

double ieee80ToDouble(const unsigned char *p)
{
    int sign = (p[0] & 0x80) ? -1 : 1;
    int exponent = ((p[0] & 0x7f) << 8) | p[1];
    uint64_t mantissa = ((uint64_t)be32(p + 2) << 32) | be32(p + 6);
    if (exponent == 0 && mantissa == 0) return 0.0;
    return sign * ldexp((double)mantissa, exponent - 16383 - 63);
}

bool readAiff(const std::string &path, AudioData &out)
{
    FILE *f = fopen(path.c_str(), "rb");
    if (!f) return false;
    std::vector<unsigned char> bytes;
    {
        unsigned char chunk[65536];
        size_t n;
        while ((n = fread(chunk, 1, sizeof chunk, f)) > 0) bytes.insert(bytes.end(), chunk, chunk + n);
        fclose(f);
    }
    if (bytes.size() < 12 || memcmp(&bytes[0], "FORM", 4) != 0) return false;
    bool aifc = memcmp(&bytes[8], "AIFC", 4) == 0;
    if (!aifc && memcmp(&bytes[8], "AIFF", 4) != 0) return false;

    int channels = 0, bits = 0;
    uint32_t frames = 0;
    double rate = 0;
    const unsigned char *ssnd = NULL;
    size_t ssndSize = 0;
    size_t pos = 12;
    while (pos + 8 <= bytes.size()) {
        const unsigned char *id = &bytes[pos];
        uint32_t size = be32(&bytes[pos + 4]);
        size_t body = pos + 8;
        if (body + size > bytes.size()) size = (uint32_t)(bytes.size() - body);
        if (memcmp(id, "COMM", 4) == 0 && size >= 18) {
            channels = be16(&bytes[body]);
            frames = be32(&bytes[body + 2]);
            bits = be16(&bytes[body + 6]);
            rate = ieee80ToDouble(&bytes[body + 8]);
            if (aifc && size >= 22) {
                const unsigned char *comp = &bytes[body + 18];
                if (memcmp(comp, "NONE", 4) != 0 && memcmp(comp, "twos", 4) != 0) return false;
            }
        } else if (memcmp(id, "SSND", 4) == 0 && size >= 8) {
            uint32_t offset = be32(&bytes[body]);
            ssnd = &bytes[body + 8 + offset];
            ssndSize = size - 8 - offset;
        }
        pos = body + size + (size & 1);
    }
    if (!ssnd || channels <= 0 || rate <= 0 || (bits != 8 && bits != 16 && bits != 24 && bits != 32)) return false;
    size_t bytesPer = bits / 8;
    size_t avail = ssndSize / (bytesPer * channels);
    if (frames > avail) frames = (uint32_t)avail;

    out.sampleRate = (float)rate;
    out.channels = channels;
    out.mono.assign(frames, 0.f);
    const float scale = 1.f / (float)channels;
    const unsigned char *p = ssnd;
    for (uint32_t i = 0; i < frames; ++i) {
        float s = 0.f;
        for (int c = 0; c < channels; ++c) {
            int32_t v = 0;
            switch (bits) {
            case 8:  v = (int8_t)p[0]; s += v / 128.f; break;
            case 16: v = (int16_t)be16(p); s += v / 32768.f; break;
            case 24: v = ((int32_t)((p[0] << 24) | (p[1] << 16) | (p[2] << 8))) >> 8; s += v / 8388608.f; break;
            case 32: v = (int32_t)be32(p); s += v / 2147483648.f; break;
            }
            p += bytesPer;
        }
        out.mono[i] = s * scale;
    }
    return true;
}

} // namespace

bool readAudioFile(const std::string &path, AudioData &out, std::string *error)
{
    std::string ext = lowerExt(path);
    typedef bool (*Reader)(const std::string &, AudioData &);
    Reader readers[5];
    int n = 0;
    // Try the reader matching the extension first, then the rest.
    if (ext == "wav" || ext == "wave") readers[n++] = readWav;
    else if (ext == "flac") readers[n++] = readFlac;
    else if (ext == "mp3") readers[n++] = readMp3;
    else if (ext == "ogg" || ext == "oga") readers[n++] = readOgg;
    else if (ext == "aif" || ext == "aiff" || ext == "aifc") readers[n++] = readAiff;
    Reader all[] = { readWav, readFlac, readOgg, readAiff, readMp3 };
    for (int i = 0; i < 5; ++i) {
        bool dup = false;
        for (int j = 0; j < n; ++j) if (readers[j] == all[i]) dup = true;
        if (!dup) readers[n++] = all[i];
    }
    FILE *probe = fopen(path.c_str(), "rb");
    if (!probe) {
        if (error) *error = "cannot open " + path;
        return false;
    }
    fclose(probe);
    for (int i = 0; i < n; ++i) {
        out = AudioData();
        if (readers[i](path, out) && !out.mono.empty()) return true;
    }
    if (error) *error = "unsupported or corrupt audio file: " + path + " (WAV, AIFF, FLAC, MP3 and Ogg Vorbis are supported)";
    return false;
}

bool writeWavMono(const std::string &path, const float *samples, size_t frames,
                  float sampleRate, std::string *error)
{
    drwav_data_format fmt;
    fmt.container = drwav_container_riff;
    fmt.format = DR_WAVE_FORMAT_PCM;
    fmt.channels = 1;
    fmt.sampleRate = (drwav_uint32)sampleRate;
    fmt.bitsPerSample = 16;
    drwav wav;
    if (!drwav_init_file_write(&wav, path.c_str(), &fmt, NULL)) {
        if (error) *error = "cannot write " + path;
        return false;
    }
    std::vector<drwav_int16> pcm(frames);
    for (size_t i = 0; i < frames; ++i) {
        float v = samples[i];
        if (v > 1.f) v = 1.f;
        if (v < -1.f) v = -1.f;
        pcm[i] = (drwav_int16)(v * 32767.f);
    }
    drwav_write_pcm_frames(&wav, frames, &pcm[0]);
    drwav_uninit(&wav);
    return true;
}

} // namespace chordino
