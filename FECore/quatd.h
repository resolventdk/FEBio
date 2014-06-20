#ifndef _QUATD_H_02082007_
#define _QUATD_H_02082007_

#include "vec3d.h"
#include "mat3d.h"

//-----------------------------------------------------------------------------
//! This class implements a quaternion. 

class quatd
{
public:
	// constructors
	quatd () { x = y = z = w = 0.f; }

	quatd( const double angle, vec3d v)
	{
		w = (double) cos(angle * 0.5);

		double sina = (double) sin(angle * 0.5);

		v.unit();
		
		x = v.x*sina;
		y = v.y*sina;
		z = v.z*sina;
	}

	quatd(vec3d v)
	{
        double angle = v.unit();
        
		w = (double) cos(angle * 0.5);
        
		double sina = (double) sin(angle * 0.5);
        
		x = v.x*sina;
		y = v.y*sina;
		z = v.z*sina;
	}
    
	quatd (vec3d& v1, vec3d& v2)
	{
		vec3d n = v1^v2;
		n.unit();

		double d = v1*v2;

		double sina = (double) sqrt((1.0-d)*0.5);
		double cosa = (double) sqrt((1.0+d)*0.5);

		w = cosa;

		x = n.x*sina;
		y = n.y*sina;
		z = n.z*sina;

	}

	quatd(const double qx, const double qy, const double qz, const double qw = 0.0)
	{
		w = qw;
		x = qx;
		y = qy;
		z = qz;
	}

	bool operator != (const quatd& q) { return ((x!=q.x) || (y!=q.y) || (z!=q.z) || (w!=q.w)); }

	// addition and substraction

	quatd operator + (const quatd& q) const
	{
		return quatd(x + q.x, y + q.y, z + q.z, w + q.w);
	}

	quatd operator - (const quatd& q) const
	{
		return quatd(x - q.x, y - q.y, z - q.z, w - q.w);
	}

	quatd& operator += (const quatd& q)
	{
		x += q.x;
		y += q.y;
		z += q.z;
		w += q.w;

		return *this;
	}

	quatd& operator -= (const quatd& q)
	{
		x -= q.x;
		y -= q.y;
		z -= q.z;
		w -= q.w;

		return *this;
	}


	// multiplication

	quatd operator * (const quatd& q) const
	{
		double qw = w*q.w - x*q.x - y*q.y - z*q.z;
		double qx = w*q.x + x*q.w + y*q.z - z*q.y;
		double qy = w*q.y + y*q.w + z*q.x - x*q.z;
		double qz = w*q.z + z*q.w + x*q.y - y*q.x;

		return quatd(qx, qy, qz, qw);
	}

	quatd& operator *= (const quatd& q)
	{
		double qw = w*q.w - x*q.x - y*q.y - z*q.z;
		double qx = w*q.x + x*q.w + y*q.z - z*q.y;
		double qy = w*q.y + y*q.w + z*q.x - x*q.z;
		double qz = w*q.z + z*q.w + x*q.y - y*q.x;

		x = qx;
		y = qy;
		z = qz;
		w = qw;

		return *this;
	}

	quatd operator*(const double a) const
	{
		return quatd(x*a, y*a, z*a, w*a);
	}

	// division

	quatd operator / (const double a) const
	{
		return quatd(x/a, y/a, z/a, w/a);
	}

	quatd& operator /= (const double a)
	{
		x /= a;
		y /= a;
		z /= a;
		w /= a;

		return *this;
	}

	// Special ops

	quatd Conjugate() const { return quatd(-x, -y, -z, w); }

	double Norm() const { return w*w + x*x + y*y + z*z; } 

	void MakeUnit() 
	{
		double N = (double) sqrt(w*w + x*x + y*y + z*z);

		if (N != 0)
		{
			x /= N;
			y /= N;
			z /= N;
			w /= N;
		}
		else w = 1.f;
	}

	quatd Inverse() const
	{
		double N = w*w + x*x + y*y + z*z;

		return quatd(-x/N, -y/N, -z/N, w/N);
	}

	quatd Cayley() const
    {
        double a = GetAngle();
        vec3d v = GetVector();
        return quatd(2*atan(a/2),v);
    }
    
	double DotProduct(const quatd& q) const
	{
		return w*q.w + x*q.x + y*q.y + z*q.z;
	}

	vec3d GetVector() const
	{
		vec3d r(x,y,z);
		r.unit();
		return r;
	}

	double GetAngle() const
	{
		return (double)(acos(w)*2.0);
	}

	// use only when *this is unit vector
	void RotateVector(vec3d& v) const
	{
		if ((w == 0) || ((x==0) && (y==0) && (z==0))) return;

		// v*q^-1
		double qw = v.x*x + v.y*y + v.z*z;
		double qx = v.x*w - v.y*z + v.z*y;
		double qy = v.y*w - v.z*x + v.x*z;
		double qz = v.z*w - v.x*y + v.y*x;

		// q* (v* q^-1)
		v.x = w*qx + x*qw + y*qz - z*qy;
		v.y = w*qy + y*qw + z*qx - x*qz;
		v.z = w*qz + z*qw + x*qy - y*qx;
	}

	// use only when *this is unit vector
	vec3d operator * (const vec3d& r)
	{
		vec3d n = r;

		// v*q^-1
		double qw = n.x*x + n.y*y + n.z*z;
		double qx = n.x*w - n.y*z + n.z*y;
		double qy = n.y*w - n.z*x + n.x*z;
		double qz = n.z*w - n.x*y + n.y*x;

		// q* (v* q^-1)
		n.x = w*qx + x*qw + y*qz - z*qy;
		n.y = w*qy + y*qw + z*qx - x*qz;
		n.z = w*qz + z*qw + x*qy - y*qx;

		return n;
	}

	void RotateVectorP(double* v, double* r) const
	{
		double qw = v[0]*x + v[1]*y + v[2]*z;
		double qx = v[0]*w - v[1]*z + v[2]*y;
		double qy = v[1]*w - v[2]*x + v[0]*z;
		double qz = v[2]*w - v[0]*y + v[1]*x;

		r[0] = w*qx + x*qw + y*qz - z*qy;
		r[1] = w*qy + y*qw + z*qx - x*qz;
		r[2] = w*qz + z*qw + x*qy - y*qx;
	}

	//! Convert a quaternion to a matrix
	mat3d RotationMatrix()
	{
		return mat3d(
			w*w + x*x - y*y - z*z,
			2.0*(x*y - w*z),
			2.0*(x*z + w*y),
			2.0*(x*y + w*z),
			w*w - x*x + y*y - z*z,
			2.0*(y*z - w*x),
			2.0*(x*z - w*y),
			2.0*(y*z + w*x),
			w*w - x*x - y*y + z*z);
	}

public:
	double w;
	double x, y, z;
};

inline quatd operator * (const double a, const quatd& q)
{
	return q*a;
}


#endif //_QUATD_H_02082007_
