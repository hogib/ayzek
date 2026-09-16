#include "weights.hpp"

#include <bit>
#include <cstring>
#include <format>
#include <fstream>
#include <stdexcept>

namespace ayzek {

namespace {

std::size_t item_size(DType t) {
    switch (t) {
        case DType::F32: return 4;
        case DType::F64: return 8;
        case DType::I64: return 8;
        case DType::U8: return 1;
    }
    throw std::runtime_error("unknown dtype");
}

template <typename T>
T read_le(std::ifstream& f, const std::string& path) {
    T v{};
    if (!f.read(reinterpret_cast<char*>(&v), sizeof v)) throw std::runtime_error(path + ": truncated");
    // The format is little-endian; every target this runs on is too.
    if constexpr (std::endian::native != std::endian::little) v = std::byteswap(v);
    return v;
}

template <typename T>
std::span<const T> typed(const Tensor& t, DType want) {
    if (t.dtype != want) throw std::runtime_error("tensor dtype mismatch");
    return {reinterpret_cast<const T*>(static_cast<const void*>(t.bytes.data())), t.bytes.size() / sizeof(T)};
}

std::string shape_str(const std::vector<std::size_t>& s) {
    std::string out = "(";
    for (std::size_t i = 0; i < s.size(); ++i) out += (i ? ", " : "") + std::to_string(s[i]);
    return out + ")";
}

}  // namespace

std::size_t Tensor::numel() const noexcept {
    std::size_t n = 1;
    for (auto d : shape) n *= d;
    return n;
}

std::span<const float> Tensor::f32() const { return typed<float>(*this, DType::F32); }
std::span<const double> Tensor::f64() const { return typed<double>(*this, DType::F64); }
std::span<const std::int64_t> Tensor::i64() const { return typed<std::int64_t>(*this, DType::I64); }

Weights Weights::load(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + path);
    char magic[4];
    if (!f.read(magic, 4) || std::memcmp(magic, "AYZW", 4) != 0) throw std::runtime_error(path + ": not an AYZW file");
    if (read_le<std::uint32_t>(f, path) != 1) throw std::runtime_error(path + ": unsupported version");
    const auto count = read_le<std::uint32_t>(f, path);

    Weights w;
    w.path_ = path;
    for (std::uint32_t i = 0; i < count; ++i) {
        std::string name(read_le<std::uint16_t>(f, path), '\0');
        f.read(name.data(), static_cast<std::streamsize>(name.size()));
        Tensor t;
        t.dtype = static_cast<DType>(read_le<std::uint8_t>(f, path));
        const auto ndim = read_le<std::uint8_t>(f, path);
        for (int d = 0; d < ndim; ++d) t.shape.push_back(read_le<std::uint32_t>(f, path));
        const auto nbytes = read_le<std::uint64_t>(f, path);
        if (nbytes != t.numel() * item_size(t.dtype)) throw std::runtime_error(path + ": size mismatch in " + name);
        // operator new returns storage aligned for every fundamental type, so the
        // reinterpret in `typed` lands on correctly aligned float and double.
        t.bytes.resize(nbytes);
        if (!f.read(reinterpret_cast<char*>(t.bytes.data()), static_cast<std::streamsize>(nbytes)))
            throw std::runtime_error(path + ": truncated in " + name);
        if (name == "__meta__") w.meta_.assign(reinterpret_cast<const char*>(t.bytes.data()), t.bytes.size());
        else w.tensors_.emplace(std::move(name), std::move(t));
    }
    return w;
}

const Tensor& Weights::at(const std::string& name, std::vector<std::size_t> shape) const {
    auto it = tensors_.find(name);
    if (it == tensors_.end()) throw std::runtime_error(path_ + ": missing tensor " + name);
    if (!shape.empty() && it->second.shape != shape)
        throw std::runtime_error(std::format("{}: {} has shape {}, expected {}", path_, name,
                                             shape_str(it->second.shape), shape_str(shape)));
    return it->second;
}

std::vector<float> Weights::vec(const std::string& name, std::vector<std::size_t> shape) const {
    auto s = at(name, std::move(shape)).f32();
    return {s.begin(), s.end()};
}

}  // namespace ayzek
