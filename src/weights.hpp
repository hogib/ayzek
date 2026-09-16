#pragma once

// AYZW: the flat tensor file tools/export_models.py writes. See
// docs/impl/02-weights.md for the layout.

#include <cstddef>
#include <cstdint>
#include <map>
#include <span>
#include <string>
#include <vector>

namespace ayzek {

enum class DType : std::uint8_t { F32 = 0, F64 = 1, I64 = 2, U8 = 3 };

struct Tensor {
    DType dtype = DType::F32;
    std::vector<std::size_t> shape;
    std::vector<std::byte> bytes;

    [[nodiscard]] std::size_t numel() const noexcept;
    [[nodiscard]] std::span<const float> f32() const;
    [[nodiscard]] std::span<const double> f64() const;
    [[nodiscard]] std::span<const std::int64_t> i64() const;
};

class Weights {
public:
    // Throws std::runtime_error on a missing or malformed file. Loading happens
    // once at startup, so an exception here costs nothing on the hot path.
    static Weights load(const std::string& path);

    [[nodiscard]] bool has(const std::string& name) const { return tensors_.contains(name); }
    // Throws if absent, or if `shape` is given and does not match.
    [[nodiscard]] const Tensor& at(const std::string& name, std::vector<std::size_t> shape = {}) const;
    [[nodiscard]] std::vector<float> vec(const std::string& name, std::vector<std::size_t> shape = {}) const;
    [[nodiscard]] const std::string& meta() const noexcept { return meta_; }
    [[nodiscard]] const std::string& path() const noexcept { return path_; }

private:
    std::map<std::string, Tensor> tensors_;
    std::string meta_;
    std::string path_;
};

}  // namespace ayzek
