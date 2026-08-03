#pragma once

#include <atomic>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>

// Small RAII helpers shared by tests that need files with names expected by
// the production readers.  Keep these helpers independent of the C++ std
// module because this header is included before `import std;` in the tests.

class TemporaryFile
{
 public:
  TemporaryFile(std::filesystem::path path, std::string_view contents) : path_(std::move(path))
  {
    std::error_code error;
    existed_ = std::filesystem::exists(path_, error);
    if (error)
    {
      throw std::runtime_error("Unable to inspect temporary test file: " + path_.string());
    }

    if (existed_)
    {
      std::ifstream input(path_, std::ios::binary);
      if (!input)
      {
        throw std::runtime_error("Unable to preserve existing test file: " + path_.string());
      }
      previousContents_.assign(std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>());
    }

    write(contents);
  }

  TemporaryFile(const TemporaryFile &) = delete;
  TemporaryFile &operator=(const TemporaryFile &) = delete;
  TemporaryFile(TemporaryFile &&) = delete;
  TemporaryFile &operator=(TemporaryFile &&) = delete;

  ~TemporaryFile()
  {
    if (existed_)
    {
      std::ofstream output(path_, std::ios::binary | std::ios::trunc);
      if (output) output << previousContents_;
    }
    else
    {
      std::error_code error;
      std::filesystem::remove(path_, error);
    }
  }

  [[nodiscard]] const std::filesystem::path &path() const noexcept { return path_; }

 private:
  void write(std::string_view contents)
  {
    std::ofstream output(path_, std::ios::binary | std::ios::trunc);
    if (!output)
    {
      throw std::runtime_error("Unable to create temporary test file: " + path_.string());
    }
    output.write(contents.data(), static_cast<std::streamsize>(contents.size()));
    if (!output)
    {
      throw std::runtime_error("Unable to write temporary test file: " + path_.string());
    }
  }

  std::filesystem::path path_;
  bool existed_{false};
  std::string previousContents_;
};

class TemporaryDirectory
{
 public:
  TemporaryDirectory() : path_(makeUniqueDirectory()) {}

  TemporaryDirectory(const TemporaryDirectory &) = delete;
  TemporaryDirectory &operator=(const TemporaryDirectory &) = delete;

  TemporaryDirectory(TemporaryDirectory &&other) noexcept : path_(std::move(other.path_)) { other.path_.clear(); }

  TemporaryDirectory &operator=(TemporaryDirectory &&other) noexcept
  {
    if (this != &other)
    {
      remove();
      path_ = std::move(other.path_);
      other.path_.clear();
    }
    return *this;
  }

  ~TemporaryDirectory() { remove(); }

  [[nodiscard]] const std::filesystem::path &path() const noexcept { return path_; }

  void write(const std::filesystem::path &relativePath, std::string_view contents) const
  {
    const std::filesystem::path outputPath = path_ / relativePath;
    std::error_code error;
    std::filesystem::create_directories(outputPath.parent_path(), error);
    if (error)
    {
      throw std::runtime_error("Unable to create temporary test directory: " + outputPath.parent_path().string());
    }

    std::ofstream output(outputPath, std::ios::binary | std::ios::trunc);
    if (!output)
    {
      throw std::runtime_error("Unable to create temporary test file: " + outputPath.string());
    }
    output.write(contents.data(), static_cast<std::streamsize>(contents.size()));
    if (!output)
    {
      throw std::runtime_error("Unable to write temporary test file: " + outputPath.string());
    }
  }

 private:
  static std::filesystem::path makeUniqueDirectory()
  {
    static std::atomic<std::uint64_t> counter{0};
    const auto timestamp = std::chrono::steady_clock::now().time_since_epoch().count();
    const std::filesystem::path base = std::filesystem::temp_directory_path();

    for (std::uint64_t attempt = 0; attempt < 128; ++attempt)
    {
      const std::uint64_t sequence = counter.fetch_add(1, std::memory_order_relaxed);
      const std::filesystem::path candidate =
          base / ("raspa3-test-" + std::to_string(timestamp) + "-" + std::to_string(sequence));
      std::error_code error;
      if (std::filesystem::create_directory(candidate, error)) return candidate;
      if (error && error != std::errc::file_exists)
      {
        throw std::runtime_error("Unable to create temporary test directory: " + error.message());
      }
    }

    throw std::runtime_error("Unable to allocate a unique temporary test directory");
  }

  void remove() noexcept
  {
    if (path_.empty()) return;
    std::error_code error;
    std::filesystem::remove_all(path_, error);
  }

  std::filesystem::path path_;
};
