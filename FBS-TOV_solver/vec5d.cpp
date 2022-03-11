#include "vec5d.hpp"

// Constructor:
vec5d::vec5d()
{
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
    v[3] = 0.0;
    v[4] = 0.0;
}
vec5d::vec5d(const double& a, const double& b, const double& c, const double& d, const double& e)
{
    v[0] = a;
    v[1] = b;
    v[2] = c;
    v[3] = d;
    v[4] = e;
}

// []-Operator:
double& vec5d::operator[] (const int& idx)
{
    return v[idx];
}
const double& vec5d::operator[] (const int& idx) const
{
    return v[idx];
}

// Vector x Scalar:
vec5d operator* (const vec5d& v, const double& c)
{
    return vec5d(c * v[0], c * v[1], c * v[2], c * v[3], c * v[4]);
}
vec5d operator* (const double& c, const vec5d& v)
{
    return vec5d(c * v[0], c * v[1], c * v[2], c * v[3], c * v[4]);
}

// Vector x Vector:
vec5d vec5d::operator+ (const vec5d& rhs)
{
    return vec5d(v[0] + rhs[0], v[1] + rhs[1], v[2] + rhs[2], v[3] + rhs[3], v[4] + rhs[4]);
}
vec5d vec5d::operator- (const vec5d& rhs)
{
    return vec5d(v[0] - rhs[0], v[1] - rhs[1], v[2] - rhs[2], v[3] - rhs[3], v[4] - rhs[4]);
}
const vec5d vec5d::operator+ (const vec5d& rhs) const
{
    return vec5d(v[0] + rhs[0], v[1] + rhs[1], v[2] + rhs[2], v[3] + rhs[3], v[4] + rhs[4]);
}
const vec5d vec5d::operator- (const vec5d& rhs) const
{
    return vec5d(v[0] - rhs[0], v[1] - rhs[1], v[2] - rhs[2], v[3] - rhs[3], v[4] - rhs[4]);
}
vec5d& vec5d::operator+= (const vec5d& rhs)
{
    v[0] += rhs[0];
	v[1] += rhs[1];
	v[2] += rhs[2];
	v[3] += rhs[3];
    v[4] += rhs[4];
	return *this;
}
vec5d& vec5d::operator-= (const vec5d& rhs)
{
    v[0] -= rhs[0];
	v[1] -= rhs[1];
	v[2] -= rhs[2];
	v[3] -= rhs[3];
    v[4] -= rhs[4];
	return *this;
}

// =-Operator:
vec5d& vec5d::operator= (const vec5d& rhs)
{
    v[0] = rhs[0];
	v[1] = rhs[1];
	v[2] = rhs[2];
	v[3] = rhs[3];
    v[4] = rhs[4];
	return *this;
}

//Static functions:
double vec5d::max(const vec5d& v)
{
    double max = std::abs(v[0]);
	if (std::abs(v[1]) >= max) max = std::abs(v[1]);
	if (std::abs(v[2]) >= max) max = std::abs(v[2]);
	if (std::abs(v[3]) >= max) max = std::abs(v[3]);
    if (std::abs(v[4]) >= max) max = std::abs(v[4]);
	return max;
}
double vec5d::dot(const vec5d& lhs, const vec5d& rhs)
{
    return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2] + lhs[3] * rhs[3] + lhs[4] * rhs[4];
}