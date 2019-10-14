#pragma once
#include "BlackMesh.h"
#include "Primal_Dual_graph.h"
#include "MeshSampling.h"
#include "Eigen/Dense"
#include "SparseTensor.h"
#include "VisualProcesser.h"
#include<math.h>

#include <boost/iterator.hpp>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>

struct Unary_Element {
	// ccid key
	double weight;

	double pos_score;
	double neg_score;

	std::set<std::pair<int, int>> confid_pos;// curvid-confid
	std::set<std::pair<int, int>> confid_neg;

	Unary_Element() {
		pos_score = 1.0;
		neg_score = 1.0;
		weight = 1.0;
	}
	Unary_Element(int ccsiz) {
		pos_score = 1.0;
		neg_score = 1.0;
		weight = 1.0/(double)ccsiz;
	}
};

struct Binary_Element{
	//std::pair<int,int> curvid_dualid; key
	int inidi;
	int inidj;
	int diri;
	int dirj;
	double score;
	double weight;
	std::set<int> confid;

	Binary_Element() {
		inidi = -1;
		inidj = -1;
		diri = -1;
		dirj = -1;
		score = 1.0;
		weight = 1.0;
	}

};

class Tribunal
{
private:
	typedef CGAL::Simple_cartesian<double> K;
	// My own point type
	struct My_point {
		double x;
		double y;
		double z;
		My_point(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
		My_point(const std::vector<double> &pos) { x = pos[0]; y = pos[1]; z = pos[2]; }
	};
	// The triangles are stored in a flat vector of points (a triangle soup): 
	// three consecutive points represent a triangle
	typedef std::vector<My_point>::const_iterator Point_iterator;
	// defines the iterator over triangles needed by the tree:
	class Triangle_iterator
		: public boost::iterator_adaptor<
		Triangle_iterator               // Derived
		, Point_iterator                  // Base
		, boost::use_default              // Value
		, boost::forward_traversal_tag    // CategoryOrTraversal
		>
	{
	public:
		Triangle_iterator()
			: Triangle_iterator::iterator_adaptor_() {}
		explicit Triangle_iterator(Point_iterator p)
			: Triangle_iterator::iterator_adaptor_(p) {}
	private:
		friend class boost::iterator_core_access;
		void increment() { this->base_reference() += 3; }
	};
	// The following primitive provides the conversion facilities between
	// my own triangle and point types and the CGAL ones
	struct My_triangle_primitive {
	public:
		typedef Triangle_iterator Id;
		// the CGAL types returned
		typedef K::Point_3    Point;
		typedef K::Triangle_3 Datum;
	private:
		Id m_it; // this is what the AABB tree will store internally
	public:
		My_triangle_primitive() {} // default constructor needed
								   // the following constructor is the one that receives the iterators from the 
								   // iterator range given as input to the AABB_tree
		My_triangle_primitive(Triangle_iterator a)
			: m_it(a) {}
		Id id() const { return m_it; }
		// on the fly conversion from the internal data
		// to the CGAL types
		Datum datum() const
		{
			Point_iterator p_it = m_it.base();
			Point p(p_it->x, p_it->y, p_it->z);
			++p_it;
			Point q(p_it->x, p_it->y, p_it->z);
			++p_it;
			Point r(p_it->x, p_it->y, p_it->z);
			return Datum(p, q, r); // assembles a triangle from three points
		}
		// returns one point which must be on the primitive
		Point reference_point() const
		{
			return Point(m_it->x, m_it->y, m_it->z);
		}
	};
	// types
	typedef CGAL::AABB_traits<K, My_triangle_primitive> My_AABB_traits;
	typedef CGAL::AABB_tree<My_AABB_traits> Tree;


	BlackMesh::BlackMesh<double> *_msh;
	Primal_Dual_graph *_gh;
	MeshSampling *_sampler;

	std::vector<My_point> _triangles;
	Tree _tree;

	//std::map<std::pair<int, int>, vector<vector<Binary_element>>> *_binary_score;// ccidi, ccidj -> vec(00/01/10/11)(curvid, confid, val)
	//std::vector<std::pair<double, std::vector<std::pair<int, int>>>> *_unary_score;// 2*ccid+0/1 -> score val, vec(edgeid, confid)//->set(edgeid, conid)
	//std::map<int, vector<Binary_element>> *_unary_score_from_binary;//2*ccid+0/1 -> vec(edgeid, confid, val)


	std::map<std::pair<int, int>, Binary_Element> *_binary_score;
	std::map<std::pair<int, int>, Binary_Element> *_unary_score_from_binary;
	std::vector<Unary_Element> *_unary_score;

	void build_cgal_mesh();
	void build_AABB_Tree();
	bool if_intersects(array<double, 3>& pos, std::vector<double> &nm, double CUTTING_EPS);

	//useful func
	inline double dot_prod(double* v1, double *v2) {
		return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
	}
	inline double norm_vec(double *v) {
		return sqrt(dot_prod(v, v));
	}
	inline double cos_angle(double *v1, double *v2) {
		return dot_prod(v1, v2);
	}
	inline double edge_length(double *v1, double *v2) {
		double tmp[3] = { v1[0] - v2[0],v1[1] - v2[1] ,v1[2] - v2[2]};
		return norm_vec(&(tmp[0]));
	}

	std::vector<int> cc_incident_edge_num;
	void assemble_number_of_edge();


	//for visualize socre

	void Tribunal::GroundColorMix(const double x_in, double & r, double & g, double & b);
	
	//dual edge weight
	double total_dual_edge_weight;

public:
	Tribunal();
	~Tribunal();
	void setMesh(BlackMesh::BlackMesh<double> *msh) { _msh = msh; }
	void setSamplig(MeshSampling *sampler) { _sampler = sampler; }
	void setPrimalDualGraph(Primal_Dual_graph *gh) { _gh = gh; }
	void setScoreList(std::map<std::pair<int, int>, Binary_Element> *binary_score,
	std::map<std::pair<int, int>, Binary_Element> *unary_score_from_binary,
	std::vector<Unary_Element> *unary_score) {
		_unary_score = unary_score;
		_binary_score = binary_score;
		_unary_score_from_binary = unary_score_from_binary;
	}

	inline void init() {
		assemble_number_of_edge();
	}

	void update_unary_score(double CUTTING_EPS,double IGNOR_TRI_UNDER=1e-7, 
		std::vector<double>* pt_posSco=NULL, std::vector<double>* pt_negSco=NULL );

	void upate_binary_score(double IGNOR_TRI_UNDER = 1e-7);

	void update_weight(std::vector<double> &curve_len, double uw, double bw);

	void visualize_unary_score(const char* filename);

	inline double get_dual_edge_total_weight() { return total_dual_edge_weight; }
};

