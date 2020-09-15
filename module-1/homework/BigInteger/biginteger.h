#ifndef BIGINTEGER_H
#define BIGINTEGER_H

#include <cstdint>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

namespace big_integer {

class BigInteger {
 public:
  BigInteger();
  BigInteger(std::string_view string);  // NOLINT(google-explicit-constructor)
  BigInteger(int value);                // NOLINT(google-explicit-constructor)

  BigInteger(const BigInteger& other);
  BigInteger(BigInteger&& other) noexcept;

  BigInteger& operator=(const BigInteger& other) = default;
  BigInteger& operator=(BigInteger&& other) noexcept = default;

  template <typename T>
  BigInteger operator+(const T& other) const {
    return *this + BigInteger(other);
  }
  BigInteger operator+(const BigInteger& other) const;

  BigInteger operator-(const BigInteger& other) const;
  template <typename T>
  BigInteger operator-(const T& other) const {
    return *this - BigInteger(other);
  }

  BigInteger operator*(const BigInteger& other) const;
  template <typename T>
  BigInteger operator*(const T& other) const {
    return *this * BigInteger(other);
  }

  BigInteger operator/(const BigInteger& other) const;
  template <typename T>
  BigInteger operator/(const T& other) const {
    return *this / BigInteger(other);
  }

  BigInteger operator%(const BigInteger& other) const;
  template <typename T>
  BigInteger operator%(const T& other) const {
    return *this % BigInteger(other);
  }

  BigInteger& operator++();
  BigInteger operator++(int);  // NOLINT(cert-dcl21-cpp)
  BigInteger& operator--();
  BigInteger operator--(int);  // NOLINT(cert-dcl21-cpp)

  BigInteger operator-() const;
  BigInteger operator+() const;

  BigInteger& operator+=(const BigInteger& other);
  BigInteger& operator-=(const BigInteger& other);
  BigInteger& operator*=(const BigInteger& other);
  BigInteger& operator/=(const BigInteger& other);
  BigInteger& operator%=(const BigInteger& other);

  template <typename T>
  bool operator<(const T& other) const {
    return *this < BigInteger(other);
  }
  bool operator<(const BigInteger& other) const noexcept;

  template <typename T>
  bool operator<=(const T& other) const {
    return *this <= BigInteger(other);
  }
  bool operator<=(const BigInteger& other) const noexcept;

  bool operator>(const BigInteger& other) const noexcept;
  template <typename T>
  bool operator>(const T& other) const {
    return *this > BigInteger(other);
  }

  bool operator>=(const BigInteger& other) const noexcept;
  template <typename T>
  bool operator>=(const T& other) const {
    return *this >= BigInteger(other);
  }

  bool operator==(const BigInteger& other) const noexcept;
  template <typename T>
  bool operator==(const T& other) const {
    return *this == BigInteger(other);
  }

  bool operator!=(const BigInteger& other) const noexcept;
  template <typename T>
  bool operator!=(const T& other) const {
    return *this != BigInteger(other);
  }

  operator bool() const noexcept;  // NOLINT(google-explicit-constructor)

  [[nodiscard]] std::string toString() const;

  ~BigInteger() = default;

 private:
  static uint32_t GetHighPart(uint64_t value) noexcept;
  static uint32_t GetLowPart(uint64_t value) noexcept;

  enum CompareResult {
    LESS = -1,
    EQUAL = 0,
    GREATER = 1,
  };
  static CompareResult CompareAbsoluteValues(
      const std::vector<uint32_t>& left, const std::vector<uint32_t>& right);
  static void SummarizeAbsoluteValuesInplace(
      std::vector<uint32_t>& left, const std::vector<uint32_t>& right);
  static std::vector<uint32_t> SummarizeAbsoluteValues(
      const std::vector<uint32_t>& left, const std::vector<uint32_t>& right);
  static std::vector<uint32_t> SubtractAbsoluteValues(
      const std::vector<uint32_t>& left, const std::vector<uint32_t>& right);
  static std::vector<uint32_t> MultiplyAbsoluteValues(
      const std::vector<uint32_t>& left, const std::vector<uint32_t>& right);
  static std::vector<uint32_t> MultiplyAbsoluteValuesWithKaratsubaAlgo(
      const std::vector<uint32_t>& left, const std::vector<uint32_t>& right);
  static std::vector<uint32_t> MultiplyAbsoluteValuesNaive(
      const std::vector<uint32_t>& shorter,
      const std::vector<uint32_t>& longer);
  static std::pair<std::vector<uint32_t>, std::vector<uint32_t>>
  DivideAbsoluteValues(const std::vector<uint32_t>& left,
                       const std::vector<uint32_t>& right);
  static std::pair<std::vector<uint32_t>, uint32_t> DivideAbsoluteValueByLimb(
      const std::vector<uint32_t>& left, uint32_t limb);
  static void ShiftLimbsLeft(std::vector<uint32_t>& limbs, size_t count);
  static void ShiftLimbsRight(std::vector<uint32_t>& limbs, long count);
  static void ShiftBitsLeft(std::vector<uint32_t>& limbs, size_t count);
  static void ShiftBitsRight(std::vector<uint32_t>& limbs, long count);
  static constexpr uint32_t GetLeastBitsMask(size_t bits_count);
  static void Normalize(std::vector<uint32_t>& limbs) noexcept;
  static size_t GetMostSignificantBitPosition(uint32_t value) noexcept;
  static size_t GetMostSignificantBitPosition(
      const std::vector<uint32_t>& value) noexcept;

  /**
   * This constructor should be used from internal arithmetic functions
   * It provides necessary sign checks
   * @param is_negative is the number negative
   * @param digits the digits of the number
   */
  BigInteger(bool is_negative, std::vector<uint32_t> digits);
  [[nodiscard]] bool IsZeroed() const;

 private:
  bool is_negative_{};

  /**
   * I use uint32_t because I want to multiply numbers in uint64_t and
   * handle carried bits with bitwise operations
   */
  std::vector<uint32_t> limbs_;

  static constexpr size_t KARATSUBA_STOP_RECURSION_SIZE = 70;
};

std::istream& operator>>(std::istream& in, BigInteger& value);
std::ostream& operator<<(std::ostream& out, const BigInteger& value);

}  // namespace big_integer

#endif
