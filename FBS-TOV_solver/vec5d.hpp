#ifndef VEC5D_HPP
#define VEC5D_HPP

#include <cmath>

// 5-dim vector class with double precision
class vec5d
{
public:
	// Constructor:
	vec5d();
    vec5d(const double& a, const double& b, const double& c, const double& d, const double& e);

	// []-Operator:
	double& operator[] (const int& idx);
	const double& operator[] (const int& idx) const;

	// Vector x Scalar:
	friend vec5d operator* (const vec5d& v, const double& c);
	friend vec5d operator* (const double& c, const vec5d& v);

	// Vector x Vector:
	vec5d operator+ (const vec5d& rhs);
	vec5d operator- (const vec5d& rhs);
	const vec5d operator+ (const vec5d& rhs) const;
	const vec5d operator- (const vec5d& rhs) const;
	vec5d& operator+= (const vec5d& rhs);
	vec5d& operator-= (const vec5d& rhs);

	// = Operator:
	vec5d& operator= (const vec5d& rhs);

	//Static functions:
	static double max(const vec5d& v);  // obtain maximum (infinity norm)
	static double dot(const vec5d& lhs, const vec5d& rhs);  // dor product of two vectors

private:
	double v[5];
};

#endif