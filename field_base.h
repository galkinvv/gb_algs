#pragma once
#include <exception>
enum class ExactSubtractionResultInfo
{
	Zero,
	NonZero
};
enum class ApproxSubtractionResultInfo
{
	MaybeZero,
	Nonzero
};
struct unexact_divisor_exception: std::exception{};
struct cant_detect_zero_equality_exception: std::exception{};