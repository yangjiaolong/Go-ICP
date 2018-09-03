/********************************************************************
An ICP Implementation for 3D Pointset Registration
Last modified: Feb 13, 2014

"Go-ICP: Solving 3D Registration Efficiently and Globally Optimally"
Jiaolong Yang, Hongdong Li, Yunde Jia
International Conference on Computer Vision (ICCV), 2013

Copyright (C) 2013 Jiaolong Yang (BIT and ANU)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*********************************************************************/

#ifndef JLY_ICP3D_HPP
#define JLY_ICP3D_HPP

#include "matrix.h"
#include "nanoflann.hpp"
using namespace nanoflann;


// A custom data set class to use nanoflann
template <typename T>
struct PointCloud
{
	struct Point
	{
		T  x,y,z;
	};

	std::vector<Point>  pts;

	// Must return the number of data points
	inline size_t kdtree_get_point_count() const { return pts.size(); }

	// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
	inline T kdtree_distance(const T *p1, const size_t idx_p2,size_t size) const
	{
		const T d0=p1[0]-pts[idx_p2].x;
		const T d1=p1[1]-pts[idx_p2].y;
		const T d2=p1[2]-pts[idx_p2].z;
		return d0*d0+d1*d1+d2*d2;
	}

	// Returns the dim'th component of the idx'th point in the class:
	// Since this is inlined and the "dim" argument is typically an immediate value, the
	//  "if/else's" are actually solved at compile time.
	inline T kdtree_get_pt(const size_t idx, int dim) const
	{
		if (dim==0) return pts[idx].x;
		else if (dim==1) return pts[idx].y;
		else return pts[idx].z;
	}

	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX &bb) const { return false; }

};

struct POINTREF
{
	double dis;
	int id_data;
	int id_model;
};

template <typename T>
class ICP3D
{
public:
	
	ICP3D();
	~ICP3D();
	size_t max_iter_def;
	T err_diff_def;
	T trim_fraction;
	bool do_trim;
	void Build(T * model, size_t n);
	T Run(T * data, size_t n, Matrix & R, Matrix & t);
	T Run(T * data, size_t n, Matrix & R, Matrix & t, size_t max_iter);
	T Run(T * data, size_t n, Matrix & R, Matrix & t, T err_diff);
	T Run(T * data, size_t n, Matrix & R, Matrix & t, size_t max_iter, T err_diff);

private:

	static int cmp(const void * a, const void * b);

	PointCloud<T> model_;

	KDTreeSingleIndexAdaptor<
		L2_Simple_Adaptor<T, PointCloud<T> > ,
		PointCloud<T>,
		3 /* dim */
	> * kdtree;
};

template <typename T>
ICP3D<T>::ICP3D()
{
	max_iter_def = 10000;
	err_diff_def = 0.000001;
	trim_fraction = 0;
	do_trim = true;

	kdtree = NULL;
}

template <typename T>
ICP3D<T>::~ICP3D()
{
	if(kdtree != NULL)
		delete(kdtree);
}

template <typename T>
void ICP3D<T>::Build(T * model, size_t n)
{

	if(kdtree != NULL)
		delete(kdtree);

	model_.pts.resize(n);

	size_t i, idx;
	for(i = 0; i < n; i++)
	{
		idx = i*3;
		model_.pts[i].x = model[idx];
		model_.pts[i].y = model[idx+1];
		model_.pts[i].z = model[idx+2];
	}

	kdtree = new KDTreeSingleIndexAdaptor<
		L2_Simple_Adaptor<T, PointCloud<T> > ,
		PointCloud<T>,
		3 /* dim */
	>(3 /*dim*/, model_, KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );

	kdtree->buildIndex();
}

template <typename T>
int ICP3D<T>::cmp(const void * a, const void * b)
{
	return ((struct POINTREF*)a)->dis > ((struct POINTREF*)b)->dis ? 1:-1;
}

template <typename T>
T ICP3D<T>::Run(T * data, size_t n, Matrix & R, Matrix & t)
{
	return Run(data, n, R, t, max_iter_def, err_diff_def);
}

template <typename T>
T ICP3D<T>::Run(T * data, size_t n, Matrix & R, Matrix & t, size_t max_iter)
{
	return Run(data, n, R, t, max_iter, err_diff_def);
}

template <typename T>
T ICP3D<T>::Run(T * data, size_t n, Matrix & R, Matrix & t, T err_diff)
{
	return Run(data, n, R,  t, max_iter_def, err_diff);
}

template <typename T>
T ICP3D<T>::Run(T * data, size_t n, Matrix & R, Matrix & t, size_t max_iter, T err_diff)
{
  size_t num;

	T query[3];
	std::vector<size_t> ret_index(1);
	std::vector<T> out_dist_sqr(1);

	if(do_trim)
	{
	  num = (int)(n*(1-trim_fraction));
	}
	else
	{
	  num = n;
	}

	struct POINTREF * points = (struct POINTREF *)malloc(sizeof(struct POINTREF)*n);

	// init matrix for point correspondences
	Matrix p_m(num,3); // model
	Matrix p_d(num,3); // data

	// init mean
	Matrix mu_m(1,3);
	Matrix mu_d(1,3);

	size_t iter, idx, i;
	T err = -1, err_new;
	for(iter = 0; iter < max_iter; iter++)
	{
		T r00 = R.val[0][0]; T r01 = R.val[0][1]; T r02 = R.val[0][2];
		T r10 = R.val[1][0]; T r11 = R.val[1][1]; T r12 = R.val[1][2];
		T r20 = R.val[2][0]; T r21 = R.val[2][1]; T r22 = R.val[2][2];
		T t0  = t.val[0][0]; T t1  = t.val[1][0]; T t2  = t.val[2][0];

		err_new = 0;
		for(i = 0; i < n; i++)
		{
			idx = i*3;

			//transform point according to R and T
			query[0] = r00*data[idx+0] + r01*data[idx+1] + r02*data[idx+2] + t0;
			query[1] = r10*data[idx+0] + r11*data[idx+1] + r12*data[idx+2] + t1;
			query[2] = r20*data[idx+0] + r21*data[idx+1] + r22*data[idx+2] + t2;


			//search nearest neighbor
			kdtree->knnSearch(&query[0], 1, &ret_index[0], &out_dist_sqr[0]);

			points[i].dis = out_dist_sqr[0];
			points[i].id_data = i;
			points[i].id_model = ret_index[0];
		}

		if(do_trim)
		{
		  qsort(points, n, sizeof(struct POINTREF), cmp);
		}

		for(i = 0; i < num; i++)
		{
			// set model point
			p_m.val[i][0] = model_.pts[points[i].id_model].x; mu_m.val[0][0] += p_m.val[i][0];
			p_m.val[i][1] = model_.pts[points[i].id_model].y; mu_m.val[0][1] += p_m.val[i][1];
			p_m.val[i][2] = model_.pts[points[i].id_model].z; mu_m.val[0][2] += p_m.val[i][2];

			idx = points[i].id_data*3;
			// set query point
			p_d.val[i][0] = r00*data[idx+0] + r01*data[idx+1] + r02*data[idx+2] + t0; mu_d.val[0][0] += p_d.val[i][0];
			p_d.val[i][1] = r10*data[idx+0] + r11*data[idx+1] + r12*data[idx+2] + t1; mu_d.val[0][1] += p_d.val[i][1];
			p_d.val[i][2] = r20*data[idx+0] + r21*data[idx+1] + r22*data[idx+2] + t2; mu_d.val[0][2] += p_d.val[i][2];

			err_new += points[i].dis;
		}

		if(err > 0 && err - err_new < err_diff*num)
			break;
		err = err_new;

		// subtract mean
		mu_m = mu_m/(T)n;
		mu_d = mu_d/(T)n;
		Matrix q_m = p_m - Matrix::ones(num,1)*mu_m;
		Matrix q_t = p_d - Matrix::ones(num,1)*mu_d;

		// compute relative rotation matrix R and translation vector T
		Matrix H = ~q_t*q_m;
		Matrix U,W,V;
		H.svd(U,W,V);
		Matrix R_ = V*~U;

		//There are some problems with Matrix::det(), so it is not used
		//R11*(R22*R33-R23*R32)
		T a = R_.val[0][0]*(R_.val[1][1]* R_.val[2][2] - R_.val[1][2]*R_.val[2][1]);
		//R12*(R21*R33-R23*R31)
		T b = -R_.val[0][1]*(R_.val[1][0]* R_.val[2][2] - R_.val[1][2]*R_.val[2][0]);
		//R13*(R21*R32-R22*R31)
		T c = R_.val[0][2]*(R_.val[1][0]* R_.val[2][1] - R_.val[1][1]*R_.val[2][0]);
		T det = a+b+c;

		Matrix tmp = Matrix::eye(3);
		tmp.val[2][2] = det;

		R_ = V*tmp*~U;

		Matrix t_ = ~mu_m - R_*~mu_d;

		// compose transformation
		R = R_*R;
		t = R_*t + t_;
	}

	return err_new;
}


#endif
