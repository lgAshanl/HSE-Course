//
// Created by aminimulin on 08.09.2020.
//

#include "biginteger.h"

#include <exception>
#include <istream>
#include <stdexcept>

big_integer::BigInteger::BigInteger() : is_negative_(false), limbs_{0} {}

big_integer::BigInteger::BigInteger(int value) {
  is_negative_ = value < 0;
  limbs_.push_back(uint32_t(std::abs(value)));
}

big_integer::BigInteger::BigInteger(const big_integer::BigInteger& other) =
    default;

big_integer::BigInteger::BigInteger(big_integer::BigInteger&& other) noexcept
    : is_negative_(other.is_negative_), limbs_(std::move(other.limbs_)) {}

big_integer::BigInteger big_integer::BigInteger::operator+(
    const big_integer::BigInteger& other) const {
  if (is_negative_ == other.is_negative_) {
    return BigInteger(is_negative_,
                      SummarizeAbsoluteValues(limbs_, other.limbs_));
  } else {
    auto compare_result = CompareAbsoluteValues(limbs_, other.limbs_);
    if (compare_result == CompareResult::LESS) {
      return BigInteger(!is_negative_,
                        SubtractAbsoluteValues(other.limbs_, limbs_));
    } else if (compare_result == CompareResult::GREATER) {
      return BigInteger(is_negative_,
                        SubtractAbsoluteValues(limbs_, other.limbs_));
    } else if (compare_result == CompareResult::EQUAL) {
      return 0;
    } else {
      throw std::runtime_error("Comparing result has unexpected value.");
    }
  }
}

uint32_t big_integer::BigInteger::GetHighPart(uint64_t value) noexcept {
  return static_cast<uint32_t>(value >> (sizeof(uint32_t) * 8U));
}

uint32_t big_integer::BigInteger::GetLowPart(uint64_t value) noexcept {
  return static_cast<uint32_t>(value & (~uint32_t(0)));
}

big_integer::BigInteger big_integer::BigInteger::operator-(
    const big_integer::BigInteger& other) const {
  if (is_negative_ != other.is_negative_) {
    return BigInteger(is_negative_,
                      SummarizeAbsoluteValues(limbs_, other.limbs_));
  } else {
    auto compare_result = CompareAbsoluteValues(limbs_, other.limbs_);
    if (compare_result == CompareResult::LESS) {
      return BigInteger(!is_negative_,
                        SubtractAbsoluteValues(other.limbs_, limbs_));
    } else if (compare_result == CompareResult::GREATER) {
      return BigInteger(is_negative_,
                        SubtractAbsoluteValues(limbs_, other.limbs_));
    } else if (compare_result == CompareResult::EQUAL) {
      return 0;
    } else {
      throw std::runtime_error("Comparing result has unexpected value");
    }
  }
}

big_integer::BigInteger::CompareResult
big_integer::BigInteger::CompareAbsoluteValues(
    const std::vector<uint32_t>& left, const std::vector<uint32_t>& right) {
  if (left.size() < right.size()) {
    return CompareResult::LESS;
  } else if (left.size() > right.size()) {
    return CompareResult::GREATER;
  } else {
    int pos = int(left.size() - 1U);
    while (pos >= 0) {
      auto index = size_t(pos);
      if (left[index] < right[index]) {
        return CompareResult::LESS;
      } else if (left[index] > right[index]) {
        return CompareResult::GREATER;
      } else {
        --pos;
      }
    }
    return CompareResult::EQUAL;
  }
}

std::vector<uint32_t> big_integer::BigInteger::SummarizeAbsoluteValues(
    const std::vector<uint32_t>& left, const std::vector<uint32_t>& right) {
  std::vector<uint32_t> result(std::max(left.size(), right.size()));

  std::size_t pos = 0;
  uint32_t carry = 0;
  while (pos < result.size()) {
    uint64_t limbs_sum;
    if (pos < right.size()) {
      limbs_sum = uint64_t(carry) + uint64_t(left[pos]) + uint64_t(right[pos]);
    } else {
      limbs_sum = uint64_t(left[pos]) + carry;
    }
    result[pos] = GetLowPart(limbs_sum);
    carry = GetHighPart(limbs_sum);
    pos++;
  }

  if (carry != 0) {
    result.push_back(carry);
  }

  return result;
}

std::vector<uint32_t> big_integer::BigInteger::SubtractAbsoluteValues(
    const std::vector<uint32_t>& left, const std::vector<uint32_t>& right) {
  // It is private function. We assume that left > right.
  std::vector<uint32_t> result(left.size(), 0);
  size_t pos = 0;
  uint64_t borrow = 0;
  while (pos < right.size()) {
    // Use uint64_t to calculate because uint32_t overflow occurs when
    // (right[pos] == std::numeric_limits<uint32_t>::max() && borrow == 1)
    if (uint64_t(left[pos]) >= uint64_t(right[pos]) + borrow) {
      result[pos] = uint32_t(uint64_t(left[pos] - right[pos]) - borrow);
      borrow = 0;
    } else {
      result[pos] = uint32_t(uint64_t(left[pos]) +
                             (uint64_t(1) << (sizeof(uint32_t) * 8U)) -
                             uint64_t(right[pos]) - borrow);
      borrow = 1;
    }
    pos++;
  }
  while (pos < left.size() && borrow == 1) {
    if (left[pos] >= borrow) {
      result[pos] = uint32_t(left[pos] - borrow);
      borrow = 0;
    } else {
      result[pos] = uint32_t(std::numeric_limits<uint32_t>::max() - borrow);
      // here borrow still must be equal to 1
    }
    pos++;
  }
  while (pos < left.size()) {
    result[pos] = left[pos];
    ++pos;
  }
  Normalize(result);
  if (result.empty()) {
    throw std::runtime_error("Unexpected empty result of subtraction");
  }
  return result;
}

big_integer::BigInteger big_integer::BigInteger::operator-() const {
  if (IsZeroed()) return *this;
  return big_integer::BigInteger(!is_negative_, limbs_);
}

big_integer::BigInteger big_integer::BigInteger::operator+() const {
  return *this;
}

big_integer::BigInteger& big_integer::BigInteger::operator+=(
    const big_integer::BigInteger& other) {
  auto sum = *this + other;
  this->is_negative_ = sum.is_negative_;
  this->limbs_ = std::move(sum.limbs_);
  return *this;
}

big_integer::BigInteger& big_integer::BigInteger::operator-=(
    const big_integer::BigInteger& other) {
  auto subtraction_result = *this - other;
  this->is_negative_ = subtraction_result.is_negative_;
  this->limbs_ = std::move(subtraction_result.limbs_);
  return *this;
}

big_integer::BigInteger& big_integer::BigInteger::operator*=(
    const big_integer::BigInteger& other) {
  auto multiplication = *this * other;
  this->is_negative_ = multiplication.is_negative_;
  this->limbs_ = std::move(multiplication.limbs_);
  return *this;
}

big_integer::BigInteger& big_integer::BigInteger::operator/=(
    const big_integer::BigInteger& other) {
  auto division = *this / other;
  this->limbs_ = std::move(division.limbs_);
  this->is_negative_ = division.is_negative_;
  return *this;
}

big_integer::BigInteger& big_integer::BigInteger::operator%=(
    const big_integer::BigInteger& other) {
  auto modulo = *this % other;
  this->is_negative_ = modulo.is_negative_;
  this->limbs_ = std::move(modulo.limbs_);
  return *this;
}

bool big_integer::BigInteger::operator<(
    const big_integer::BigInteger& other) const noexcept {
  auto compare_result = CompareAbsoluteValues(this->limbs_, other.limbs_);
  return compare_result == CompareResult::LESS;
}

bool big_integer::BigInteger::operator<=(
    const big_integer::BigInteger& other) const noexcept {
  auto compare_result = CompareAbsoluteValues(this->limbs_, other.limbs_);
  return compare_result == CompareResult::LESS ||
         compare_result == CompareResult::EQUAL;
}

bool big_integer::BigInteger::operator>(
    const big_integer::BigInteger& other) const noexcept {
  auto compare_result = CompareAbsoluteValues(this->limbs_, other.limbs_);
  return compare_result == CompareResult::GREATER;
}

bool big_integer::BigInteger::operator>=(
    const big_integer::BigInteger& other) const noexcept {
  auto compare_result = CompareAbsoluteValues(this->limbs_, other.limbs_);
  return compare_result == CompareResult::GREATER ||
         compare_result == CompareResult::EQUAL;
}

bool big_integer::BigInteger::operator==(
    const big_integer::BigInteger& other) const noexcept {
  if (is_negative_ != other.is_negative_) return false;
  return limbs_ == other.limbs_;
}

bool big_integer::BigInteger::operator!=(
    const big_integer::BigInteger& other) const noexcept {
  auto compare_result = CompareAbsoluteValues(this->limbs_, other.limbs_);
  return compare_result != CompareResult::EQUAL;
}

big_integer::BigInteger::operator bool() const noexcept {
  return !(limbs_.size() == 1 && limbs_[0] == 0);
}

big_integer::BigInteger  // NOLINT(cert-dcl21-cpp)
big_integer::BigInteger::operator--(int) {
  auto result = *this;
  *this = *this - BigInteger(1);
  return result;
}

big_integer::BigInteger  // NOLINT(cert-dcl21-cpp)
big_integer::BigInteger::operator++(int) {
  auto result = *this;
  *this = *this + BigInteger(1);
  return result;
}

big_integer::BigInteger& big_integer::BigInteger::operator++() {
  *this = *this + BigInteger(1);
  return *this;
}

big_integer::BigInteger& big_integer::BigInteger::operator--() {
  *this = *this - BigInteger(1);
  return *this;
}

std::vector<uint32_t> big_integer::BigInteger::MultiplyAbsoluteValues(
    const std::vector<uint32_t>& left, const std::vector<uint32_t>& right) {
  return MultiplyAbsoluteValuesWithKaratsubaAlgo(left, right);
}

std::vector<uint32_t>
big_integer::BigInteger::MultiplyAbsoluteValuesWithKaratsubaAlgo(
    const std::vector<uint32_t>& left, const std::vector<uint32_t>& right) {
  if (left.size() < KARATSUBA_STOP_RECURSION_SIZE) {
    return MultiplyAbsoluteValuesNaive(right, left);
  } else if (right.size() < KARATSUBA_STOP_RECURSION_SIZE) {
    return MultiplyAbsoluteValuesNaive(left, right);
  } else {
    auto half_size = std::min(std::min(left.size(), right.size()),
                              std::max(left.size() / 2, right.size() / 2));
    std::vector a(left.begin() + half_size, left.end());
    std::vector b(left.begin(), left.begin() + half_size);
    std::vector c(right.begin() + half_size, right.end());
    std::vector d(right.begin(), right.begin() + half_size);

    auto mul_ac = MultiplyAbsoluteValuesWithKaratsubaAlgo(a, c);
    auto mul_bd = MultiplyAbsoluteValuesWithKaratsubaAlgo(b, d);
    auto sum_ab = SummarizeAbsoluteValues(a, b);
    auto sum_cd = SummarizeAbsoluteValues(c, d);

    auto mul_sum_ab_sum_cd =
        MultiplyAbsoluteValuesWithKaratsubaAlgo(sum_ab, sum_cd);
    auto center = SubtractAbsoluteValues(
        mul_sum_ab_sum_cd, SummarizeAbsoluteValues(mul_ac, mul_bd));

    // I am trying to avoid unnecessary copying and memory allocations,
    // so I can use mul_bd instance to store result
    // mul_bd is part of the final answer, so I need to add remaining values
    auto& result = mul_bd;

    ShiftLimbsLeft(center, half_size);
    ShiftLimbsLeft(mul_ac, half_size * 2);

    result = SummarizeAbsoluteValues(result,
                                     SummarizeAbsoluteValues(mul_ac, center));
    return result;
  }
}

std::vector<uint32_t> big_integer::BigInteger::MultiplyAbsoluteValuesNaive(
    const std::vector<uint32_t>& shorter, const std::vector<uint32_t>& longer) {
  std::vector<uint32_t> result(shorter.size() + longer.size(), 0);
  for (size_t short_pos = 0; short_pos < shorter.size(); ++short_pos) {
    uint64_t carry = 0;
    for (size_t long_pos = 0; long_pos < longer.size(); ++long_pos) {
      uint64_t value =
          carry + uint64_t(shorter[short_pos]) * uint64_t(longer[long_pos]) +
          uint64_t(result[long_pos + short_pos]);
      result[long_pos + short_pos] = GetLowPart(value);
      carry = GetHighPart(value);
    }
    if (carry) {
      result[short_pos + longer.size()] = uint32_t(carry);
    }
  }
  Normalize(result);
  return result;
}

std::pair<std::vector<uint32_t>, uint32_t>
big_integer::BigInteger::DivideAbsoluteValueByLimb(
    const std::vector<uint32_t>& left, uint32_t right) {
  std::vector<uint32_t> result;
  result.resize(left.size());
  uint64_t carry = 0;
  int pos = int(left.size()) - 1;
  while (pos >= 0) {
    uint64_t cur = left[pos] + (uint64_t(carry) << (sizeof(uint32_t) * 8U));
    result[pos] = cur / right;
    carry = cur % right;
    pos--;
  }
  Normalize(result);
  return std::make_pair(result, carry);
}

std::pair<std::vector<uint32_t>, std::vector<uint32_t>>
big_integer::BigInteger::DivideAbsoluteValues(
    const std::vector<uint32_t>& left, const std::vector<uint32_t>& right) {
  if (CompareAbsoluteValues(left, right) == CompareResult::LESS) {
    return {{0}, left};
  }
  std::vector<uint32_t> dividend = left;
  std::vector<uint32_t> result = {0U};
  std::vector<uint32_t> current_value = right;
  {
    auto left_msb_position = GetMostSignificantBitPosition(left);
    auto right_msb_position = GetMostSignificantBitPosition(right);
    ShiftBitsLeft(current_value, left_msb_position - right_msb_position);
  }
  while (CompareAbsoluteValues(current_value, right) != LESS) {
    ShiftBitsLeft(result, 1);
    if (CompareAbsoluteValues(dividend, current_value) != CompareResult::LESS) {
      // set least bit
      result[0] |= 1U;
      dividend = SubtractAbsoluteValues(dividend, current_value);
    }
    ShiftBitsRight(current_value, 1);
  }
  return {result, dividend};
}

big_integer::BigInteger::BigInteger(bool is_negative,
                                    std::vector<uint32_t> limbs)
    : is_negative_(is_negative), limbs_(std::move(limbs)) {
  if (IsZeroed()) is_negative_ = false;
}

big_integer::BigInteger big_integer::BigInteger::operator*(
    const big_integer::BigInteger& other) const {
  if (this->IsZeroed() || other.IsZeroed()) {
    return 0;
  }
  auto result = big_integer::BigInteger(
      this->is_negative_ != other.is_negative_,
      MultiplyAbsoluteValuesWithKaratsubaAlgo(this->limbs_, other.limbs_));
  return result;
}

big_integer::BigInteger big_integer::BigInteger::operator/(
    const big_integer::BigInteger& other) const {
  if (other.IsZeroed()) {
    throw std::logic_error("Division by zero");
  }
  auto result = big_integer::BigInteger(
      this->is_negative_ != other.is_negative_,
      DivideAbsoluteValues(this->limbs_, other.limbs_).first);
  if (result.IsZeroed()) {
    result.is_negative_ = false;
  }
  return result;
}

big_integer::BigInteger big_integer::BigInteger::operator%(
    const big_integer::BigInteger& other) const {
  if (other.IsZeroed()) {
    throw std::logic_error("Division by zero");
  }
  auto result = big_integer::BigInteger(
      this->is_negative_,
      DivideAbsoluteValues(this->limbs_, other.limbs_).second);
  if (result.IsZeroed()) {
    result.is_negative_ = false;
  }
  return result;
}

std::string big_integer::BigInteger::toString() const {
  if (IsZeroed()) return "0";
  std::string result;
  if (is_negative_) {
    result += '-';
  }
  std::vector<uint32_t> limbs = limbs_;
  for (;;) {
    if (limbs.size() == 1U && limbs[0] == 0) break;
    auto division_result = DivideAbsoluteValueByLimb(limbs, 10U);
    result += char(division_result.second + '0');
    std::swap(limbs, division_result.first);
  }
  if (is_negative_) {
    std::reverse(result.begin() + 1, result.end());
  } else {
    std::reverse(result.begin(), result.end());
  }
  return result;
}

void big_integer::BigInteger::ShiftLimbsLeft(std::vector<uint32_t>& limbs,
                                             size_t count) {
  limbs.insert(limbs.begin(), count, 0);
}

void big_integer::BigInteger::ShiftLimbsRight(std::vector<uint32_t>& limbs,
                                              long count) {
  count = std::min(count, long(limbs.size()));
  limbs.erase(limbs.begin(), limbs.begin() + count);
}

void big_integer::BigInteger::ShiftBitsLeft(std::vector<uint32_t>& limbs,
                                            size_t count) {
  if (count == 0) return;
  auto limbs_count = count / (sizeof(uint32_t) * 8U);
  ShiftLimbsLeft(limbs, limbs_count);

  count %= sizeof(uint32_t) * 8U;
  uint32_t carry = 0;
  for (uint32_t& limb_value : limbs) {
    uint64_t new_limb_value = (uint64_t(limb_value) << count) + carry;
    limb_value = GetLowPart(new_limb_value);
    carry = GetHighPart(new_limb_value);
  }
  if (carry != 0) {
    limbs.push_back(carry);
  }
}

void big_integer::BigInteger::ShiftBitsRight(std::vector<uint32_t>& limbs,
                                             long count) {
  if (count == 0) return;
  auto limbs_count = count / long(sizeof(uint32_t) * 8U);
  ShiftLimbsRight(limbs, limbs_count);

  count %= sizeof(uint32_t) * 8U;
  uint32_t borrow = 0;
  for (auto i = long(limbs.size()) - 1; i >= 0; --i) {
    uint32_t new_limb_value = (limbs[size_t(i)] >> size_t(count)) | borrow;
    borrow = (GetLeastBitsMask(size_t(count)) & limbs[size_t(i)])
             << (sizeof(uint32_t) * 8U - size_t(count));
    limbs[size_t(i)] = new_limb_value;
  }
  Normalize(limbs);
}

constexpr uint32_t big_integer::BigInteger::GetLeastBitsMask(
    size_t bits_count) {
  uint32_t result = 0;
  for (size_t i = 0; i < bits_count; ++i) {
    result |= (uint32_t(1) << i);
  }
  return result;
}

void big_integer::BigInteger::Normalize(std::vector<uint32_t>& limbs) noexcept {
  if (limbs.empty()) {
    limbs = {0};
    return;
  }
  auto count_lead_zeros =
      size_t(std::find_if(limbs.rbegin(), limbs.rend(),
                          [](uint32_t value) { return value != 0; }) -
             limbs.rbegin());
  limbs.erase(limbs.end() - (long)std::min(limbs.size() - 1U, count_lead_zeros),
              limbs.end());
}

size_t big_integer::BigInteger::GetMostSignificantBitPosition(
    uint32_t value) noexcept {
  if (!value) {
    return ~size_t(0);
  }
  size_t pos = 0;
  while (value) {
    value >>= 1U;
    pos++;
  }
  return pos - 1;
}

size_t big_integer::BigInteger::GetMostSignificantBitPosition(
    const std::vector<uint32_t>& value) noexcept {
  return GetMostSignificantBitPosition(value.back()) +
         (value.size() - 1) * sizeof(uint32_t) * 8U;
}

big_integer::BigInteger::BigInteger(std::string_view string) {
  limbs_ = {0};
  is_negative_ = false;
  if (string.empty()) {
    throw std::runtime_error(
        "Cannot construct BigInteger from provided string");
  }
  size_t pos = 0;
  if (string[0] == '-') {
    is_negative_ = true;
    pos += 1;
  }
  for (; pos < string.size(); ++pos) {
    if (!isdigit(string[pos])) {
      throw std::runtime_error(
          "Cannot construct BigInteger from provided string");
    }
    *this = *this * 10 + (string[pos] - '0');
  }
}

std::ostream& big_integer::operator<<(std::ostream& out,
                                      const big_integer::BigInteger& value) {
  out << value.toString();
  return out;
}

std::istream& big_integer::operator>>(std::istream& in,
                                      big_integer::BigInteger& value) {
  char c;
  bool first = true;
  bool is_negative = false;
  value = 0;
  in >> c;
  while (!in.eof()) {
    if (first && c == '-') {
      is_negative = true;
    } else if (isdigit(c)) {
      value = value * 10 + c - '0';
    } else {
      in.putback(c);
      break;
    }
    first = false;

    // c can become -1 when we read EOF symbol from the input stream
    c = char(in.get());
  }
  if (is_negative && value != 0) {
    value = -value;
  }
  return in;
}

bool big_integer::BigInteger::IsZeroed() const {
  return limbs_.size() == 1 && limbs_[0] == 0;
}
