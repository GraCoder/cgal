#ifndef CGAL_POLYGON_2_CLIPPER_H
#define CGAL_POLYGON_2_CLIPPER_H

#include <CGAL/config.h>
#include <CGAL/Simple_cartesian.h>

#include <vector>
#include <list>
#include <iterator>
#include <queue>

namespace CGAL {
namespace Polygon2Clipper {

class ClipperException : public std::exception
{
public:
	ClipperException(const char* description) : _descrip(description) {}
	virtual ~ClipperException() {}
	virtual const char* what() const throw() { return _descrip.c_str(); }
private:
	std::string		_descrip;
};

template<typename K>
struct Rectangle2 {
	typedef typename K::FT RT;

	RT left;	RT top;
	RT right;	RT bottom;
};

enum class PolyType		{ PT_Subject, PT_Clip };
enum class ClipType		{ CT_Intersect, CT_Union, CT_Difference, CT_Xor };
//http://glprogramming.com/red/chapter11.html
enum class PolyFillType { PFT_EvenOdd, PFT_NonZero, PFT_Positive, PFT_Negative };
enum class InitOptions	{ IO_ReverseSolution = 1, IO_StrictSimple = 2, IO_PreserveCollinear = 4 };
enum class JointType	{ JT_Square, JT_Round, JT_Miter };
enum class EndType		{ ET_ClosedPolygon, ET_ClosedLine, ET_OpenButt, ET_OpenSquare, ET_OpenRound };
enum class NodeType		{ NT_Any, NT_Open, NT_Closed };

template<typename>
class Clipper;
template<typename>
class ClipperOffset;

template<typename K>
class PolyNode
{
	friend class Clipper<K>;
	friend class ClipperOffset<K>;

	typedef typename Clipper<K>::Point2 Point2;
	//typedef typename K::Point_2 Point2;
public:
	PolyNode()
		: _parent(nullptr)
		, _index(0)
		, _is_open(false)
	{ };
	virtual ~PolyNode() {};

	bool is_hole() const
	{
		bool ret = true;
		PolyNode* node = _parent;
		while (node)
		{
			ret = !ret;
			node = node->_parent;
		}
		return ret;
	}
	bool is_open() const { return _is_open; }
	int	 child_count() const { return _childs.size(); }
	PolyNode* getNext() const
	{
		if (_childs.empty())
			return get_next_sibling_up();
		else
			return _childs[0];
	}
protected:
	void add_child(PolyNode *child)
	{
		_childs.push_back(child);
		child->_parent = this;
		child->_index = _childs.size() - 1;
	}
	PolyNode* get_next_sibling_up() const
	{
		if (!_parent)
			return nullptr;
		else if (_index == _parent->_childs.size() - 1)
			return _parent->get_next_sibling_up();
		else
			return _parent->_childs[_index + 1];
	}

	bool						_is_open;
	unsigned				_index; //node index in _parent._childs
	JointType				_join_type;
	EndType					_end_type;
	std::vector<Point2>		_contour;

	PolyNode* _parent;
	std::vector<PolyNode*> _childs;
};

template<typename K>
class PolyTree : public PolyNode<K>
{
	friend class Clipper<K>;
	friend class ClipperOffset<K>;
public:
	PolyTree()
	{
	}

	~PolyTree() { clear(); };

	int	total() const
	{
		int ret = _all_nodes.size();
		if (ret > 0 && _childs[0] != _all_nodes[0])
			ret--;
		return ret;
	}
	void clear()
	{
		_childs.clear();
		for (int i = 0; i < _all_nodes.size(); i++)
			delete _all_nodes[i];
		_all_nodes.clear();
	}
	PolyNode<K>* get_first() const
	{
		if (_childs.empty())
			return nullptr;

		return _childs.front();
	}

private:
	std::vector<PolyNode<K>*> _all_nodes;
};

template<typename K>
using Path = std::vector<typename K::Point_2>;

template<typename K>
using Paths = std::vector<Path<K> >;

//------------------------------------------------------------------------------

template<typename K>
void reverse_path(Path<K> & p)
{
	std::reverse(p.begin(), p.end());
}

template<typename K>
void reverse_paths(Paths<K>& p)
{
	for (typename Paths<K>::size_type i = 0; i < p.size(); ++i)
		reverse_path(p[i]);
}

template<typename K>
double area(const Path<K> & path)
{
	int size = path.size();
	if (size < 3)
		return 0;

	double a = 0;
	for (int i = 0, j = size - 1; i < size; ++i)
	{
		a += ((double)path[j].x() + path[i].x()) * ((double)path[j].y() - path[i].y());
		j = i;
	}
	return -a * 0.5;
}

template<typename K>
bool orientation(const Path<K>& p)
{
    return area<K>(p) < 0;
}

template<typename K>
void simplify_polygon(const Path<K>& in_poly, Paths<K>& out_polys, PolyFillType fillType)
{
	Clipper<K> c;
	c.setStrictlySimple(true);
	c.add_path(in_poly, PolyType::PT_Subject, true);
	c.execute(ClipType::CT_Union, out_polys, fillType, fillType);
}

template<typename K>
void simplify_polygons(const Paths<K>& in_polys, Paths<K>& out_polys, PolyFillType fillType)
{
	Clipper<K> c;
	c.setStrictlySimple(true);
	c.add_paths(in_polys, PolyType::PT_Subject, true);
	c.execute(ClipType::CT_Union, out_polys, fillType, fillType);
}

template<typename K>
void simplify_polygons(Paths<K>& polys, PolyFillType fillType)
{
	simplify_polygons(polys, polys, fillType);
}

template<typename K>
bool slopes_near_collinear(const typename K::Point_2& pt1, const typename K::Point_2& pt2,
	const typename K::Point_2& pt3, double distSqrd)
{
	//this function is more accurate when the point that's geometrically
	//between the other 2 points is the one that's tested for distance.
	//ie makes it more likely to pick up 'spikes' ...
	if (std::abs(pt1.x() - pt2.x()) > std::abs(pt1.y() - pt2.y()))
	{
		if ((pt1.x() > pt2.x()) == (pt1.x() < pt3.x()))
			return distance_from_linesqrd(pt1, pt2, pt3) < distSqrd;
		else if ((pt2.x() > pt1.x()) == (pt2.x() < pt3.x()))
			return distance_from_linesqrd(pt2, pt1, pt3) < distSqrd;
		else
			return distance_from_linesqrd(pt3, pt1, pt2) < distSqrd;
	}
	else
	{
		if ((pt1.y() > pt2.y()) == (pt1.y() < pt3.y()))
			return distance_from_linesqrd(pt1, pt2, pt3) < distSqrd;
		else if ((pt2.y() > pt1.y()) == (pt2.y() < pt3.y()))
			return distance_from_linesqrd(pt2, pt1, pt3) < distSqrd;
		else
			return distance_from_linesqrd(pt3, pt1, pt2) < distSqrd;
	}
}

template<typename K>
bool points_are_close(typename K::Point_2 pt1, typename K::Point_2 pt2, double distSqrd)
{
	double dx = (double)pt1.x() - pt2.x();
	double dy = (double)pt1.y() - pt2.y();
	return ((dx * dx) + (dy * dy) <= distSqrd);
}

template<typename K>
void clean_polygon(const Path<K>& in_poly, Path<K>& out_poly, double distance)
{
	using OutPt = Clipper<K>::OutPt;
	//distance = proximity in units/pixels below which vertices
	//will be stripped. Default ~= sqrt(2).

	size_t size = in_poly.size();

	if (size == 0)
	{
		out_poly.clear();
		return;
	}

	OutPt* outPts = new OutPt[size];
	for (size_t i = 0; i < size; ++i)
	{
		outPts[i].pt = in_poly[i];
		outPts[i].next = &outPts[(i + 1) % size];
		outPts[i].next->prev = &outPts[i];
		outPts[i].idx = 0;
	}

	double distSqrd = distance * distance;
	OutPt* op = &outPts[0];
	while (op->idx == 0 && op->next != op->prev)
	{
		if (points_are_close(op->pt, op->prev->pt, distSqrd))
		{
			op = exclude_op(op);
			size--;
		}
		else if (points_are_close(op->prev->pt, op->next->pt, distSqrd))
		{
			exclude_op(op->next);
			op = exclude_op(op);
			size -= 2;
		}
		else if (slopes_near_collinear(op->prev->pt, op->pt, op->next->pt, distSqrd))
		{
			op = exclude_op(op);
			size--;
		}
		else
		{
			op->idx = 1;
			op = op->next;
		}
	}

	if (size < 3) size = 0;
	out_poly.resize(size);
	for (size_t i = 0; i < size; ++i)
	{
		out_poly[i] = op->pt;
		op = op->next;
	}
	delete[] outPts;
}

template<typename K>
void clean_polygon(Path<K>& poly, double distance)
{
	clean_polygon(poly, poly, distance);
}

template<typename K>
void clean_polygons(const Paths<K>& in_polys, Paths<K>& out_polys, double distance)
{
	out_polys.resize(in_polys.size());
	for (typename Paths<K>::size_type i = 0; i < in_polys.size(); ++i)
		clean_polygon(in_polys[i], out_polys[i], distance);
}

template<typename K>
void clean_polygons(Paths<K>& polys, double distance)
{
	clean_polygons(polys, polys, distance);
}

template<typename K>
void minkowski(const Path<K>& poly, const Path<K>& path, Paths<K>& solution, bool isSum, bool isClosed)
{
	int delta = (isClosed ? 1 : 0);
	size_t polyCnt = poly.size();
	size_t pathCnt = path.size();
	Paths<K> pp;
	pp.reserve(pathCnt);
	if (isSum)
		for (size_t i = 0; i < pathCnt; ++i)
		{
			Path<K> p;
			p.reserve(polyCnt);
			for (size_t j = 0; j < poly.size(); ++j)
				p.push_back(Point2(path[i].x() + poly[j].x(), path[i].y() + poly[j].y()));
			pp.push_back(p);
		}
	else
		for (size_t i = 0; i < pathCnt; ++i)
		{
			Path<K> p;
			p.reserve(polyCnt);
			for (size_t j = 0; j < poly.size(); ++j)
				p.push_back(Point2(path[i].x() - poly[j].x(), path[i].y() - poly[j].y()));
			pp.push_back(p);
		}

	solution.clear();
	solution.reserve((pathCnt + delta) * (polyCnt + 1));
	for (size_t i = 0; i < pathCnt - 1 + delta; ++i)
		for (size_t j = 0; j < polyCnt; ++j)
		{
			Path<K> quad;
			quad.reserve(4);
			quad.push_back(pp[i % pathCnt][j % polyCnt]);
			quad.push_back(pp[(i + 1) % pathCnt][j % polyCnt]);
			quad.push_back(pp[(i + 1) % pathCnt][(j + 1) % polyCnt]);
			quad.push_back(pp[i % pathCnt][(j + 1) % polyCnt]);
			if (!orientation(quad)) 
				reverse_path(quad);
			solution.push_back(quad);
		}
}

template<typename K>
void minkowski_sum(const Path<K>& pattern, const Path<K>& path, Paths<K>& solution, bool closed)
{
	minkowski(pattern, path, solution, true, closed);
	Clipper<K> c;
	c.add_paths(solution, PolyType::PT_Subject, true);
	c.execute(ClipType::CT_Union, solution, PolyFillType::PFT_NonZero, PolyFillType::PFT_NonZero);
}

template<typename K>
void translate_path(const Path<K>& input, Path<K>& output, const typename K::Point_2 delta)
{
	//precondition: input != output
	output.resize(input.size());
	for (size_t i = 0; i < input.size(); ++i)
		output[i] = Point2(input[i].x() + delta.x(), input[i].y() + delta.y());
}

template<typename K>
void minkowski_sum(const Path<K>& pattern, const Paths<K>& paths, Paths<K>& solution, bool closed)
{
	Clipper<K> c;
	for (size_t i = 0; i < paths.size(); ++i)
	{
		Paths<K> tmp;
		minkowski(pattern, paths[i], tmp, true, closed);
		c.add_paths(tmp, PolyType::PT_Subject, true);
		if (closed)
		{
			Path<K> tmp2;
			translate_path(paths[i], tmp2, pattern[0]);
			c.add_path(tmp2, PolyType::PT_Clip, true);
		}
	}
	c.execute(ClipType::CT_Union, solution, PolyFillType::PFT_NonZero, PolyFillType::PFT_NonZero);
}

template<typename K>
void minkowski_diff(const Path<K>& poly1, const Path<K>& poly2, Paths<K>& solution)
{
	minkowski(poly1, poly2, solution, false, true);
	Clipper<K> c;
	c.add_paths(solution, PolyType::PT_Subject, true);
	c.execute(ClipType::CT_Union, solution, PolyFillType::PFT_NonZero, PolyFillType::PFT_NonZero);
}

template<typename K>
void add_polynode_to_paths(const PolyNode<K>& polynode, NodeType nodetype, Paths<K>& paths)
{
	bool match = true;
	if (nodetype == NodeType::NT_Closed) match = !polynode.is_open();
	else if (nodetype == NodeType::NT_Open) return;

	if (!polynode._contour.empty() && match)
		paths.push_back(polynode._contour);
	for (int i = 0; i < polynode.child_count(); ++i)
		add_polynode_to_paths(*polynode._childs[i], nodetype, paths);
}

template<typename K>
void polytree_to_paths(const PolyTree<K>& polytree, Paths<K>& paths)
{
	paths.resize(0);
	paths.reserve(polytree.total());
	add_polynode_to_paths(polytree, NodeType::NT_Any, paths);
}

template<typename K>
void closed_paths_from_polytree(const PolyTree<K>& polytree, Paths<K>& paths)
{
	paths.resize(0);
	paths.reserve(polytree.total());
	add_polynode_to_paths(polytree, NodeType::NT_Closed, paths);
}

template<typename K>
void open_paths_from_polytree(PolyTree<K>& polytree, Paths<K>& paths)
{
	paths.resize(0);
	paths.reserve(polytree.total());
	//Open paths are top level only, so ...
	for (int i = 0; i < polytree.child_count(); ++i)
		if (polytree._childs[i]->is_open())
			paths.push_back(polytree._childs[i]->_contour);
}

//------------------------------------------------------------------------------

template<typename K = Simple_cartesian<double> >
class Clipper {
	friend class ClipperOffset<K>;

	enum class EdgeSide { esLeft = 1, esRight = 2};
public:

	typedef K							KT;
	typedef typename K::FT				FT;
	typedef typename K::FT				CT;
	typedef typename K::Point_2			Point2;
	typedef typename K::Segment_2		Segment2;

	typedef std::pair<double, double> DebugPt;
	typedef Path<K>			Path;
	typedef Paths<K>		Paths;
public:
	Clipper();
	virtual ~Clipper();

	bool	add_path(const Path& pg, PolyType type = PolyType::PT_Subject, bool closed = true);
	bool	add_paths(const Paths& pg, PolyType PolyTyp, bool closed);
	bool	execute(ClipType clipType, Paths& solution, PolyFillType fillType = PolyFillType::PFT_EvenOdd);
	bool	execute(ClipType clipType, Paths& solution, PolyFillType subjFillType, PolyFillType clipFillType);
	bool	execute(ClipType clipType, PolyTree<K>& polytree, PolyFillType fillType = PolyFillType::PFT_EvenOdd);
	bool	execute(ClipType clipType, PolyTree<K>& polytree, PolyFillType subjFillType, PolyFillType clipFillType);
	void	clear();

	Rectangle2<K>	get_bounds();

	inline bool		preserve_collinear() { return _preserve_collinear; }
	inline void		set_preserve_collinear(bool value) { _preserve_collinear = value; };
	inline bool		reverse_solution() { return _reverse_output; };
	inline void		set_reverse_solution(bool value) { _reverse_output = value; };
	inline bool		strictly_simple() { return _strict_simple; };
	inline void		set_strictly_simple(bool value) { _strict_simple = value; };
protected:	
	struct TEdge : Segment2 {
		Point2		cur;
#ifdef DEBUG
		DebugPt     dbot;
		DebugPt     dtop;
#endif
		double		dx;
		PolyType	poly_type;
		EdgeSide	side;
		int			wind_delta;
		int			wind_cnt;
		int			wind_cnt2;
		int			out_idx;
		TEdge		*next = nullptr;
		TEdge		*prev = nullptr;
		TEdge		*next_in_lml = nullptr;
		TEdge		*next_in_ael = nullptr;
		TEdge		*prev_in_ael = nullptr;
		TEdge		*next_in_sel = nullptr;
		TEdge		*prev_in_sel = nullptr;

		inline Point2 &bot() {
			return ((Point2 *)&rep())[0];
		}

		inline const Point2 &bot() const {
			return ((Point2 *)&rep())[0];
		}
		inline Point2 &top() {
			return ((Point2 *)&rep())[1];
		}

		inline const Point2 &top() const {
			return ((Point2 *)&rep())[1];
		}
	};

	struct IntersectNode {
		TEdge		*edge1 = nullptr;
		TEdge		*edge2 = nullptr;
		Point2		  pt;
	};

	struct LocalMinimum {
		TEdge		*left_bound = nullptr;
		TEdge		*right_bound = nullptr;
		FT				y;
	};

	struct OutPt;
	struct OutRec {
		int			idx;
		bool		is_hole;
		bool		is_open;
		OutPt		*pts;
		OutPt		*bot_pt;
		OutRec		*first_left;  //see comments in clipper.pas
		PolyNode<K>	*poly_node;
	};

	struct OutPt {
		int			 idx;
		OutPt		*next;
		OutPt		*prev;
		Point2		 pt;
#ifdef DEBUG
		std::pair<double, double> debug_pt;
#endif
	};

	struct Joint {
		OutPt		*out_pt1;
		OutPt		*out_pt2;
		Point2		 offPt;
	};
	
	typedef std::vector<std::vector<TEdge> > EdgeList;
	typedef std::vector<Joint> JoinList;
	typedef std::vector<IntersectNode> IntersectList;
	typedef std::vector<LocalMinimum>		MinimaList;
	typedef std::vector<OutRec *> PolyOutList;
	typedef std::vector<PolyNode<K> *> PolyNodes;
	typedef std::list<CT> MaximaList;
	typedef std::priority_queue<CT, std::vector<CT>, std::greater<CT>> ScanbeamList;

	static inline bool		is_horizontal(TEdge &e);
	static inline void		set_delta(TEdge &e);
	static inline double	get_delta(const Point2& pt1, const Point2& pt2);
	static inline void		init_edge(TEdge *e, TEdge *enext, TEdge *eprev, const Point2 &pt);
	static inline void		init_edge2(TEdge &e, PolyType type);
	static inline bool		slopes_equal(const Point2 &pt1, const Point2 &pt2, const Point2 &pt3);
	static inline bool		slopes_equal(const Point2 &pt1, const Point2 &pt2, const Point2 &pt3, const Point2 &pt4);
	static inline bool		slopes_equal(const TEdge& e1, const TEdge& e2);
	static inline bool		pt2_between_pt1_pt3(const Point2 pt1, const Point2 pt2, const Point2 pt3);
	static inline void		reverse_horizontal(TEdge &e);
	static inline void		reverse_poly_pt_links(OutPt *pp);
	static inline void		dispose_out_pts(OutPt *&pp); 
	static inline bool		horz_segments_overlap(FT, FT, FT, FT);
	static inline void		swap_sides(TEdge &edge1, TEdge &edge2);
	static inline void		swap_poly_indexs(TEdge &edge1, TEdge &edge2);
	static inline bool		e2_insert_before_e1(TEdge &e1, TEdge &e2);
	static inline bool		is_minima(TEdge *e, const FT y);
	static inline bool		is_maxima(TEdge *e, const FT y);
	static inline bool		is_intermediate(TEdge *e, const FT y);
	static inline bool		edges_adjacent(const IntersectNode &inode);
	static inline FT		topx(TEdge &edge, const FT currentY);
	static inline TEdge*	remove_edge(TEdge *e);
	static inline TEdge*	get_maxima_pair(TEdge *e);
	static inline TEdge*	get_maxima_pair_ex(TEdge *e);
	static inline TEdge*	get_next_in_ael(TEdge *e, Orientation dir);
	static inline double	area(const OutPt *op);

	//--------------------------------------------------------------------------------------------------------------

	TEdge*		find_next_local_min(TEdge *e);
	TEdge*		process_bound(TEdge *e, bool nextIsForward);

	//--------------------------------------------------------------------------------------------------------------

	bool		execute_internal();

	//--------------------------------------------------------------------------------------------------------------

	void		reset();
	void		disposeLocalMinimaList();
	bool		pop_local_minima(CT y, const LocalMinimum *&loc_min);
	void		insert_scan_beam(const CT y);
	bool		pop_scan_beam(CT &y);
	void		dispose_all_outrecs();
	void		dispose_out_rec(typename PolyOutList::size_type index);
	void		swap_pos_in_ael(TEdge *edge1, TEdge *edge2);
	void		delete_from_ael(TEdge *e);
	OutRec*		create_out_rec();
	void		update_edge_into_ael(TEdge *&e);
	bool		local_mimima_pending();

	//--------------------------------------------------------------------------------------------------------------

	void		set_winding_count(TEdge &edge);
	bool		is_even_odd_fill_type(const TEdge &edge) const;
	bool		is_even_odd_alt_fill_type(const TEdge &edge) const;
	void		get_horz_direction(TEdge& horz_edge, Orientation& dir, CT& left, CT& right);
	bool		get_overlap(const FT a1, const FT a2, const FT b1, const FT b2, FT& left, FT& right);
	void		insert_local_minimal_into_ael(const FT bot_y);
	void		insert_edge_into_ael(TEdge *edge, TEdge *start_edge);
	void		add_edge_to_sel(TEdge *edge);
	bool		pop_edge_from_sel(TEdge *&edge);
	void		copy_ael_to_sel();
	void		delete_from_sel(TEdge *e);
	void		swap_pos_in_sel(TEdge *edge1, TEdge *edge2);
	bool		is_contributing(const TEdge &edge) const;
	void		do_maxima(TEdge *e);
	void		process_horizontals();
	void		process_horizontal(TEdge *horz_edge);
	OutPt*		add_out_pt(TEdge *e, const Point2 &pt);
	OutPt*		get_last_outpt(TEdge *e);
	OutPt*		add_local_min_poly(TEdge *e1, TEdge *e2, const Point2 &pt);
	OutPt*		dup_outpt(OutPt* outPt, bool InsertAfter);
	OutPt*		get_bottom_pt(OutPt* pp);
	void		add_local_max_poly(TEdge *e1, TEdge *e2, const Point2 &pt);
	OutRec*		get_outrec(int idx);
	OutRec*		get_lowermost_rec(OutRec *outRec1, OutRec *outRec2);
	void		skip_outrec(TEdge *e);
	void		append_polygon(TEdge *e1, TEdge *e2);
	void		intersect_edges(TEdge *e1, TEdge *e2, Point2 &pt);
	bool		process_intersections(const CT top_y);
	void		build_intersect_list(const CT top_y);
	void		process_intersect_list();
	void		process_top_edges_scanbeam(const CT top_y);
	void		build_result(Paths &polys);
	void		build_result2(PolyTree<K> &polytree);
	void		set_hole_state(TEdge *e, OutRec *outrec);
	void		dispose_intersect_nodes();
	bool		fixup_intersection_order();
	void		fixup_out_polygon(OutRec &outrec);
	void		fixup_out_polyline(OutRec &outrec);
	bool		outrec1_right_of_outrec2(OutRec *outRec1, OutRec *outRec2);
	//bool is_hole(TEdge *e);
	//bool FindOwnerFromSplitRecs(OutRec &outRec, OutRec *&currOrfl);
	void		fix_hole_linkage(OutRec &outrec);
	void		add_join(OutPt *op1, OutPt *op2, const Point2 offPt);
	void		clear_joints();
	void		clear_ghost_joints();
	void		add_ghost_join(OutPt *op, const Point2 offPt);
	bool		join_horz(OutPt* op1, OutPt* op1b, OutPt* op2, OutPt* op2b, const Point2 pt, bool discardLeft);
	bool		join_points(Joint *j, OutRec *outRec1, OutRec *outRec2);
	void		join_common_edges();
	int			point_count(OutPt *pts);
	void		update_outpt_idxs(OutRec& outrec);
	void		do_simple_polygons();
	bool		first_is_bottom_pt(const OutPt* btmPt1, const OutPt* btmPt2);
	OutRec*		parse_first_left(OutRec* first_left);
	bool		poly2_contains_poly1(OutPt* out_pt1, OutPt* out_pt2);
	void		fixup_first_lefts1(OutRec *oldOutRec, OutRec *newOutRec);
	void		fixup_first_lefts2(OutRec *inner_outrec, OutRec *outer_outrec);
	void		fixup_first_lefts3(OutRec *oldOutRec, OutRec *newOutRec);
	Point2		intersect_point(TEdge &edge1, TEdge &edge2);

private:

	const static	int _unassigned = -1;  //edge not currently 'owning' a solution
	const static	int _skip = -2;        //edge that would otherwise close a path

	bool			_open_path;
	bool			_use_full_range;
	bool			_preserve_collinear;
	bool			_excute_locked;
	bool			_reverse_output;
	bool			_using_polytree;
	bool			_strict_simple;
	TEdge			*_active_edges;
	TEdge			*_sorted_edges;
#ifdef DEBUG
	int				_active_count;
	int				_sorted_count;
#endif

	typename MinimaList::iterator	_curr_lm;

	EdgeList		_edges;
	MinimaList		_minimal_list;
	PolyOutList		_poly_outs;
	ScanbeamList	_scanbeam;
	JoinList		_joints;
	JoinList		_ghost_joints;
	IntersectList	_intersect_list;
	ClipType		_clip_type;
	MaximaList		_maxima;
	PolyFillType	_clip_fill_type;
	PolyFillType	_sub_fill_type;
};

//------------------------------------------------------------------------------

template<typename K = Simple_cartesian<double> >
class ClipperOffset
{
	using PolyNode = PolyNode<K>;
	using PolyTree = PolyTree<K>;

	using OutPt    = typename Clipper<K>::OutPt;

	typedef typename Clipper<K>::Point2 Point2;
	typedef typename K::Vector_2		Vector2;

public:
	typedef Path<K>				Path;
	typedef Paths<K>			Paths;
public:
	ClipperOffset(double miterLimit = 2.0, double roundPrecision = 0.25);
	~ClipperOffset();
	void	add_path(const Path& path, JointType joinType, EndType endType);
	void	add_paths(const Paths& paths, JointType joinType, EndType endType);
	void	execute(Paths& solution, double delta);
	void	execute(PolyTree& solution, double delta);
	void	clear();
private:
	double distance_sqrt(const Point2& pt1, const Point2& pt2);
	double distance_from_linesqrd(const Point2& pt, const Point2& ln1, const Point2& ln2);

	OutPt*		exclude_op(OutPt* op);
	Vector2		get_unit_normal(const Point2& pt1, const Point2& pt2);

	void	fix_orientations();
	void	do_offset(double delta);
	void	offset_point(int j, int& k, JointType jointype);
	void	do_square(int j, int k);
	void	do_miter(int j, int k, double r);
	void	do_round(int j, int k);

	double		_miter_limit;
	double		_arc_tolerance;
	Paths		_dest_polys;
	Path		_src_poly;
	Path		_dest_poly;
	double		_delta, _sinA, _sin, _cos;
	double		_miterlimit, _steps_per_rad;
	int			_lowest_x;
	int			_lowest_y;
	PolyNode	_poly_nodes;

	std::vector<Vector2>	_normals;

	static constexpr double pi = 3.141592653589793238;
	static constexpr double pi2 = pi *2;
	static constexpr double arc_tolerance = 0.25;
};

//------------------------------------------------------------------------------

template<typename K>
Clipper<K>::Clipper()
	: _open_path(false)
	, _strict_simple(false)
	, _excute_locked(false)
	, _sorted_edges(nullptr)
	, _active_edges(nullptr)
{
}

template<typename K>
Clipper<K>::~Clipper() {}

template<typename K>
bool Clipper<K>::add_path(const Path &pg, PolyType type, bool closed)
{
    if (!closed && type == PolyType::PT_Clip)
        throw ClipperException("add_path: open paths must be subject.");

    using std::vector;
    int high_idx = pg.size() - 1;
    if (closed)
        while (high_idx > 0 && (pg[high_idx] == pg[0]))
            --high_idx;
    while (high_idx > 0 && (pg[high_idx] == pg[high_idx - 1]))
        --high_idx;
    if (closed && high_idx < 2 || !closed && high_idx < 1)
        return false;

    std::vector<TEdge> edges(high_idx + 1);
    bool is_flat = true;

    edges[1].cur = pg[1];
    init_edge(&edges[0], &edges[1], &edges[high_idx], pg[0]);
    init_edge(&edges[high_idx], &edges[0], &edges[high_idx - 1], pg[high_idx]);
    for (int i = high_idx - 1; i >= 1; --i)
    {
        init_edge(&edges[i], &edges[i + 1], &edges[i - 1], pg[i]);
    }

    TEdge *edge_start = &edges[0];

    TEdge *e = edge_start, *eLoopStop = edge_start;
    for (;;)
    {
        if (e->cur == e->next->cur && (closed || e->next != edge_start))
        {
            if (e == e->next) 
				break;
            if (e == edge_start) 
				edge_start = e->next;
            e = remove_edge(e);
            eLoopStop = e;
            continue;
        }
        if (e->prev == e->next)
            break; //only two vertices

        else if (closed && slopes_equal(e->prev->cur, e->cur, e->next->cur) && !pt2_between_pt1_pt3(e->prev->cur, e->cur, e->next->cur))
        {
            if (e == edge_start)
                edge_start = e->next;
            e = remove_edge(e);
            e = e->prev;
            eLoopStop = e;
            continue;
        }
        e = e->next;
        if (e == eLoopStop || !closed && e->next == edge_start)
            break;
    }

    if (!closed && e == e->next || closed && e->prev == e->next)
        return false;

    if (!closed) {
        _open_path = true;
        edge_start->prev->out_idx = _skip;
    }

    //3. Do second stage of edge initialization ...
    e = edge_start;
    do
    {
        init_edge2(*e, type);
        e = e->next;
        if (is_flat && e->cur.y() != edge_start->cur.y())
            is_flat = false;
    } while (e != edge_start);

    //4. Finally, add edge bounds to LocalMinima list ...

    //Totally flat paths must be handled differently when adding them
    //to LocalMinima list to avoid endless loops etc ...
	if (is_flat)
	{
		e->prev->out_idx = _skip;
		typename MinimaList::value_type loc_min;
		loc_min.y = e->bot().y();
		//loc_min.left_bound = 0;
		loc_min.left_bound = e->prev;
		loc_min.right_bound = e;
		loc_min.right_bound->side = EdgeSide::esRight;
		loc_min.right_bound->wind_delta = 0;
		for (;;)
		{
			if (e->bot().x() != e->prev->top().x())
				reverse_horizontal(*e);
			if (e->next->out_idx == _skip)
				break;
			e->next_in_lml = e->next;
			e = e->next;
		}
		_minimal_list.push_back(loc_min);
		_edges.push_back(edges);
		return true;
	}

    _edges.push_back(std::move(edges));
    bool left_bound_forward;
    TEdge *eMin = 0;

    //workaround to avoid an endless loop in the while loop below when
    //open paths have matching start and end points ...
    if (e->prev->bot() == e->prev->top())
        e = e->next;

    for (;;)
    {
        e = find_next_local_min(e);
        if (e == eMin)
            break;
        else if (!eMin)
            eMin = e;

        //E and E.prev now share a local minima (left aligned if horizontal).
        //Compare their slopes to find which starts which bound ...
        typename MinimaList::value_type loc_min;
        loc_min.y = e->bot().y();
        if (e->dx > e->prev->dx)
        {
            loc_min.left_bound = e->prev;
            loc_min.right_bound = e;
            left_bound_forward = false;
        }
        else
        {
            loc_min.left_bound = e;
            loc_min.right_bound = e->prev;
            left_bound_forward = true;
        }

        if(!closed)
            loc_min.left_bound->wind_delta = 0;
        else if (loc_min.left_bound->next == loc_min.right_bound)
            loc_min.left_bound->wind_delta =  1;
        else
            loc_min.left_bound->wind_delta = -1;
        loc_min.right_bound->wind_delta = -loc_min.left_bound->wind_delta;

        e = process_bound(loc_min.left_bound, left_bound_forward);
        if (e->out_idx == _skip)
            e = process_bound(e, left_bound_forward);

        TEdge *e2 = process_bound(loc_min.right_bound, !left_bound_forward);
        if (e2->out_idx == _skip)
            e2 = process_bound(e2, !left_bound_forward);

        //if (loc_min.left_bound->out_idx == _skip)
        //    loc_min.left_bound = 0;
        //else if (loc_min.right_bound->out_idx == _skip)
        //    loc_min.right_bound = 0;
        _minimal_list.push_back(loc_min);
        if (!left_bound_forward)
            e = e2;
    }
    return true;
}

template<typename K>
bool Clipper<K>::add_paths(const Paths& pg, PolyType PolyTyp, bool closed)
{
	bool result = false;
	for (typename Paths::size_type i = 0; i < pg.size(); ++i)
		if (add_path(pg[i], PolyTyp, closed)) 
			result = true;
	return result;
}

template<typename K>
bool Clipper<K>::execute(ClipType clipType, Paths &solution, PolyFillType fillType)
{
    return execute(clipType, solution, fillType, fillType);
}

template<typename K>
bool Clipper<K>::execute(ClipType clipType, Paths &solution, PolyFillType subjFillType, PolyFillType clipFillType)
{
    if (_excute_locked)
        return false;
    if (_open_path)
        throw ClipperException("Error: PolyTree struct is needed for open path clipping.");
    _excute_locked = true;
    _sub_fill_type = subjFillType;
    _clip_fill_type = clipFillType;
    _clip_type = clipType;
    _using_polytree = false;
    solution.resize(0);
    bool succeeded = execute_internal();
    if (succeeded)
        build_result(solution);
    dispose_all_outrecs();
    _excute_locked = false;
    return succeeded;
}

template<typename K>
bool Clipper<K>::execute(ClipType clipType, PolyTree<K>& polytree, PolyFillType fillType)
{
	return execute(clipType, polytree, fillType, fillType);
}

template<typename K>
bool Clipper<K>::execute(ClipType clipType, PolyTree<K>& polytree, PolyFillType subjFillType, PolyFillType clipFillType)
{
	if (_excute_locked)
		return false;
	_excute_locked = true;
	_sub_fill_type = subjFillType;
	_clip_fill_type = clipFillType;
	_clip_type = clipType;
	_using_polytree = true;

	bool succeeded = execute_internal();
	if (succeeded) 
		build_result2(polytree);
	dispose_all_outrecs();
	_excute_locked = false;
	return succeeded;
}

template<typename K>
void Clipper<K>::clear()
{
    disposeLocalMinimaList();
    _edges.clear();
    _open_path = false;
    _use_full_range = false;
}

template<typename K>
Rectangle2<K> Clipper<K>::get_bounds()
{
    Rectangle2<K> result;
    MinimaList::iterator lm = _minimal_list.begin();
    if (lm == _minimal_list.end())
    {
        result.left = result.top = result.right = result.bottom = 0;
        return result;
    }
    result.left = lm->left_bound->bot().x();
    result.top = lm->left_bound->bot().y();
    result.right = lm->left_bound->bot().x();
    result.bottom = lm->left_bound->bot().y();
    while (lm != _minimal_list.end())
    {
        //todo - needs fixing for open paths
        result.bottom = std::max<FT>(result.bottom, lm->left_bound->bot().y());
        TEdge *e = lm->left_bound;
        for (;;) {
            TEdge *bottomE = e;
            while (e->next_in_lml)
            {
                if (e->bot().x() < result.left)
					result.left = e->bot().x();
                if (e->bot().x() > result.right) 
					result.right = e->bot().x();
                e = e->next_in_lml;
            }
            result.left = std::min<FT>(result.left, e->bot().x());
            result.right = std::max<FT>(result.right, e->bot().x());
            result.left = std::min<FT>(result.left, e->top().x());
            result.right = std::max<FT>(result.right, e->top().x());
            result.top = std::min<FT>(result.top, e->top().y());
            if (bottomE == lm->left_bound) 
				e = lm->right_bound;
            else break;
        }
        ++lm;
    }
    return result;
}

template<typename K>
bool Clipper<K>::is_horizontal(TEdge &e)
{
    return e.dx == -DBL_MAX;
}

template<typename K>
void Clipper<K>::set_delta(TEdge &e)
{
    if (e.is_horizontal())
        e.dx = -DBL_MAX;
    else
        e.dx = (e.top().x() - e.bot().x()) / (e.top().y() - e.bot().y());
}

template<typename K>
double Clipper<K>::get_delta(const Point2& pt1, const Point2& pt2)
{
  return (pt1.y() == pt2.y()) ?  -DBL_MAX : (double)(pt2.x() - pt1.x()) / (pt2.y() - pt1.y());
}

template<typename K>
void Clipper<K>::init_edge(TEdge *e, TEdge *enext, TEdge *eprev, const Point2 &pt)
{
    std::memset(e, 0, sizeof(TEdge));
    e->next = enext;
    e->prev = eprev;
    e->cur = pt;
    e->out_idx = _unassigned;
}

template<typename K>
void Clipper<K>::init_edge2(TEdge &e, PolyType type)
{
    if (e.cur.y() <= e.next->cur.y())
    {
        e.bot() = e.cur;
        e.top() = e.next->cur;
    }
    else
    {
        e.top() = e.cur;
        e.bot() = e.next->cur;
    }
#ifdef DEBUG
    e.dbot = std::make_pair(e.bot().x(), e.bot().y());
    e.dtop = std::make_pair(e.top().x(), e.top().y());
#endif
    set_delta(e);
    e.poly_type = type;
}

template<typename K>
bool Clipper<K>::slopes_equal(const Point2 &pt1, const Point2 &pt2, const Point2 &pt3)
{
    return compare((pt1.y() - pt2.y()) * (pt2.x() - pt3.x()), (pt1.x() - pt2.x()) * (pt2.y() - pt3.y())) == EQUAL;
}

template<typename K>
bool Clipper<K>::slopes_equal(const Point2 &pt1, const Point2 &pt2, const Point2 &pt3, const Point2 &pt4)
{
    return (pt1.y() - pt2.y()) * (pt3.x() - pt4.x()) == (pt1.x() - pt2.x()) * (pt3.y() - pt4.y());
}

template<typename K>
bool Clipper<K>::slopes_equal(const TEdge& e1, const TEdge& e2)
{
    return (e1.top().y() - e1.bot().y()) * (e2.top().x() - e2.bot().x()) == 
		(e1.top().x() - e1.bot().x()) * (e2.top().y() - e2.bot().y());
}

template<typename K>
bool Clipper<K>::pt2_between_pt1_pt3(const Point2 pt1, const Point2 pt2, const Point2 pt3)
{
    if ((pt1 == pt3) || (pt1 == pt2) || (pt3 == pt2))
        return false;
    else if (pt1.x() != pt3.x())
        return (pt2.x() > pt1.x()) == (pt2.x() < pt3.x());
    else
        return (pt2.y() > pt1.y()) == (pt2.y() < pt3.y());
}

template<typename K>
void Clipper<K>::reverse_horizontal(TEdge &e)
{
    //swap horizontal edges' top() and Bottom x's so they follow the natural
    //progression of the bounds - ie so their xbots will align with the
    //adjoining lower edge. [Helpful in the process_horizontal() method.]
    std::swap(e.top(), e.bot());
}

template<typename K>
void Clipper<K>::reverse_poly_pt_links(OutPt *pp)
{
    if (!pp)
        return;
    OutPt *pp1, *pp2;
    pp1 = pp;
    do {
        pp2 = pp1->next;
        pp1->next = pp1->prev;
        pp1->prev = pp2;
        pp1 = pp2;
    } while (pp1 != pp);
}


template<typename K>
void Clipper<K>::dispose_out_pts(OutPt *&pp)
{
    if (pp == 0)
        return;
    pp->prev->next = 0;
    while (pp)
    {
        OutPt *tmpPp = pp;
        pp = pp->next;
        delete tmpPp;
    }
}

template<typename K>
bool Clipper<K>::horz_segments_overlap(FT seg1a, FT seg1b, FT seg2a, FT seg2b)
{
	if (seg1a > seg1b)
		std::swap(seg1a, seg1b);
	if (seg2a > seg2b)
		std::swap(seg2a, seg2b);
	return (seg1a < seg2b) && (seg2a < seg1b);
}

template<typename K>
void Clipper<K>::swap_sides(TEdge &e1, TEdge &e2)
{
    EdgeSide side = e1.side;
    e1.side = e2.side;
    e2.side = side;
}

template<typename K>
void Clipper<K>::swap_poly_indexs(TEdge &e1, TEdge &e2)
{
    int out_idx = e1.out_idx;
    e1.out_idx = e2.out_idx;
    e2.out_idx = out_idx;
}

template<typename K>
bool Clipper<K>::e2_insert_before_e1(TEdge &e1, TEdge &e2)
{
    if (e2.cur.x() == e1.cur.x())
    {
        if (e2.top().y() > e1.top().y())
            return e2.top().x() < topx(e1, e2.top().y());
        else
            return e1.top().x() > topx(e2, e1.top().y());
    }
    else
        return e2.cur.x() < e1.cur.x();
}

template<typename K>
bool Clipper<K>::is_minima(TEdge *e, const FT y)
{
    return e && (e->prev->next_in_lml != e) && (e->next->next_in_lml != e);
}

template<typename K>
bool Clipper<K>::is_maxima(TEdge *e, const FT y)
{
    return e && e->top().y() == y && !e->next_in_lml;
}

template<typename K>
bool Clipper<K>::is_intermediate(TEdge *e, const FT y)
{
    return e->top().y() == y && e->next_in_lml;
}

template<typename K>
bool Clipper<K>::edges_adjacent(const IntersectNode &inode)
{
    return (inode.edge1->next_in_sel == inode.edge2) ||
        (inode.edge1->prev_in_sel == inode.edge2);
}

template<typename K>
typename Clipper<K>::FT
Clipper<K>::topx(TEdge &edge, const FT currentY)
{
    return (currentY == edge.top().y()) ?
        edge.top().x() : edge.bot().x() + (edge.dx * (currentY - edge.bot().y()));
}

template<typename K>
typename Clipper<K>::TEdge *
Clipper<K>::remove_edge(TEdge *e)
{
    e->prev->next = e->next;
    e->next->prev = e->prev;
    TEdge *result = e->next;
    e->prev = 0;
    return result;
}

template<typename K>
typename Clipper<K>::TEdge *
Clipper<K>::get_maxima_pair(TEdge *e)
{
    if ((e->next->top() == e->top()) && !e->next->next_in_lml)
        return e->next;
    else if ((e->prev->top() == e->top()) && !e->prev->next_in_lml)
        return e->prev;
    else 
		return 0;
}

template<typename K>
typename Clipper<K>::TEdge *
Clipper<K>::get_maxima_pair_ex(TEdge *e)
{
    //as get_maxima_pair() but returns 0 if MaxPair isn't in AEL (unless it's horizontal)
    TEdge *result = get_maxima_pair(e);
    if (result && (result->out_idx == _skip ||
        (result->next_in_ael == result->prev_in_ael && !is_horizontal(*result))))
        return 0;
    return result;
}

template<typename K>
typename Clipper<K>::TEdge * 
Clipper<K>::get_next_in_ael(TEdge* e, Orientation dir)
{
	return dir == RIGHT_TURN ? e->next_in_ael : e->prev_in_ael;
}

template<typename K>
double Clipper<K>::area(const OutPt *op)
{
    if (!op)
        return 0;
    const OutPt *sop = op;
    double ret = 0;
    do {
        ret += (double)(op->prev->pt.x() + op->pt.x()) * (op->prev->pt.y() - op->pt.y());
        op = op->next;
    }
	while (op != sop);
    return ret * 0.5;
}

template<typename K>
typename Clipper<K>::TEdge *
Clipper<K>::find_next_local_min(TEdge *e)
{
    for (;;)
    {
        while (e->bot() != e->prev->bot() || e->cur == e->top())
            e = e->next;
        if (!is_horizontal(*e) && !is_horizontal(*e->prev))
            break;
        while (is_horizontal(*e->prev))
            e = e->prev;
        TEdge *e2 = e;
        while (is_horizontal(*e))
            e = e->next;
        if (e->top().y() == e->prev->bot().y())
            continue; //ie just an intermediate horz.
        if (e2->prev->bot().x() < e->bot().x())
            e = e2;
        break;
    }
    return e;
}

template<typename K>
typename Clipper<K>::TEdge *
Clipper<K>::process_bound(TEdge *e, bool next_forward)
{
    TEdge *result = e, *horz = 0;

    if (e->out_idx == _skip)
    {
    	//if edges still remain in the current bound beyond the skip edge then
    	//create another LocMin and call process_bound once more
    	if (next_forward)
    	{
    		while (e->top().y() == e->next->bot().y())
				e = e->next;
    		//don't include top() horizontals when parsing a bound a second time,
    		//they will be contained in the opposite bound ...
    		while (e != result && is_horizontal(*e))
				e = e->prev;
    	}
    	else
    	{
    		while (e->top().y() == e->prev->bot().y())
				e = e->prev;
    		while (e != result && is_horizontal(*e))
				e = e->next;
    	}

    	if (e == result)
    	{
    		if (next_forward) 
				result = e->next;
    		else 
				result = e->prev;
    	}
    	else
    	{
    		//there are more edges in the bound beyond result starting with e
    		if (next_forward)
    			e = result->next;
    		else
    			e = result->prev;
    		MinimaList::value_type loc_min;
    		loc_min.y = e->bot().y();
    		loc_min.left_bound = result;
    		loc_min.right_bound = e;
    		e->wind_delta = 0;
    		result = process_bound(e, next_forward);
    		_minimal_list.push_back(loc_min);
    	}
    	return result;
    }

    TEdge *edge_start;

    if (is_horizontal(*e)) {
        //We need to be careful with open paths because this may not be a
        //true local minima (ie e may be following a skip edge).
        //Also, consecutive horz. edges may start heading left before going right.
        if (next_forward)
            edge_start = e->prev;
        else
            edge_start = e->next;
        if (is_horizontal(*edge_start)) //ie an adjoining horizontal skip edge
        {
            if (edge_start->bot().x() != e->bot().x() && edge_start->top().x() != e->bot().x())
                reverse_horizontal(*e);
        }
        else if (edge_start->bot().x() != e->bot().x())
            reverse_horizontal(*e);
    }

    edge_start = e;
    if (next_forward)
    {
        while (result->top().y() == result->next->bot().y() && result->next->out_idx != _skip)
            result = result->next;
        if (is_horizontal(*result) && result->next->out_idx != _skip)
        {
            //nb: at the top() of a bound, horizontals are added to the bound
            //only when the preceding edge attaches to the horizontal's left vertex
            //unless a _skip edge is encountered when that becomes the top() divide
            horz = result;
            while (is_horizontal(*horz->prev))
				horz = horz->prev;
            if (horz->prev->top().x() > result->next->top().x()) 
				result = horz->prev;
        }
        while (e != result)
        {
            e->next_in_lml = e->next;
            if (is_horizontal(*e) && e != edge_start && e->bot().x() != e->prev->top().x())
				reverse_horizontal(*e);
            e = e->next;
        }
        if (is_horizontal(*e) && e != edge_start && e->bot().x() != e->prev->top().x())
            reverse_horizontal(*e);
        result = result->next; //move to the edge just beyond current bound
    }
    else {
        while (result->top().y() == result->prev->bot().y() && result->prev->out_idx != _skip)
            result = result->prev;
        if (is_horizontal(*result) && result->prev->out_idx != _skip)
        {
            horz = result;
            while (is_horizontal(*horz->next)) 
				horz = horz->next;
            if (horz->next->top().x() == result->prev->top().x() || horz->next->top().x() > result->prev->top().x())
				result = horz->next;
        }

        while (e != result)
        {
            e->next_in_lml = e->prev;
            if (is_horizontal(*e) && e != edge_start && e->bot().x() != e->next->top().x())
                reverse_horizontal(*e);
            e = e->prev;
        }
        if (is_horizontal(*e) && e != edge_start && e->bot().x() != e->next->top().x())
            reverse_horizontal(*e);
        result = result->prev; //move to the edge just beyond current bound
    }

    return result;
}

template<typename K>
bool Clipper<K>::execute_internal()
{
    bool succeeded = true;
    try {
        reset();
        _maxima = MaximaList();
        _sorted_edges = 0;

        succeeded = true;
        CT bot_y, top_y;
        if (!pop_scan_beam(bot_y))
            return false;
        insert_local_minimal_into_ael(bot_y);
        while (pop_scan_beam(top_y) || local_mimima_pending())
        {
            process_horizontals();
            clear_ghost_joints();
            if (!process_intersections(top_y))
            {
                succeeded = false;
                break;
            }
            process_top_edges_scanbeam(top_y);
            bot_y = top_y;
            insert_local_minimal_into_ael(bot_y);
        }
    }
    catch (...)
    {
        succeeded = false;
    }

    if (succeeded)
    {
        //fix orientations ...
        for (PolyOutList::size_type i = 0; i < _poly_outs.size(); ++i)
        {
            OutRec *outRec = _poly_outs[i];
            if (!outRec->pts || outRec->is_open)
                continue;
            if ((outRec->is_hole ^ _reverse_output) == (area(outRec->pts) > 0))
                reverse_poly_pt_links(outRec->pts);
        }

        if (!_joints.empty())
            join_common_edges();

        //unfortunately fixup_out_polygon() must be done after join_common_edges()
        for (PolyOutList::size_type i = 0; i < _poly_outs.size(); ++i) {
            OutRec *outRec = _poly_outs[i];
            if (!outRec->pts)
                continue;

            if (outRec->is_open)
                fixup_out_polyline(*outRec);
            else
                fixup_out_polygon(*outRec);
        }

        if (_strict_simple)
            do_simple_polygons();
    }

    clear_joints();
    clear_ghost_joints();
    return succeeded;
}

template<typename K>
void Clipper<K>::reset()
{
    if (_minimal_list.empty())
        return; //ie nothing to process
    std::sort(_minimal_list.begin(), _minimal_list.end(),
        [](const LocalMinimum &locMin1, const LocalMinimum &locMin2) {
            return locMin1.y < locMin2.y;
        });

    _scanbeam = ScanbeamList(); //clears/resets priority_queue
    //reset all edges ...
    for (MinimaList::iterator lm = _minimal_list.begin(); lm != _minimal_list.end(); ++lm)
    {
        insert_scan_beam(lm->y);
        TEdge *e = lm->left_bound;
        if (e)
        {
            e->cur = e->bot();
            e->side = EdgeSide::esLeft;
            e->out_idx = _unassigned;
        }

        e = lm->right_bound;
        if (e)
        {
            e->cur = e->bot();
            e->side = EdgeSide::esRight;
            e->out_idx = _unassigned;
        }
    }
    _active_edges = 0;
    _curr_lm = _minimal_list.begin();

#ifdef DEBUG
	_active_count = 0;
	_sorted_count = 0;
#endif
}

template<typename K>
void Clipper<K>::disposeLocalMinimaList()
{
    _minimal_list.clear();
    _curr_lm = _minimal_list.begin();
}

template<typename K>
bool Clipper<K>::pop_local_minima(CT y, const LocalMinimum *&loc_min)
{
    if (_curr_lm == _minimal_list.end() || (*_curr_lm).y != y)
        return false;
    loc_min = &(*_curr_lm);
    ++_curr_lm;
    return true;
}

template<typename K>
void Clipper<K>::insert_scan_beam(const CT y)
{
    _scanbeam.push(y);
}

template<typename K>
bool Clipper<K>::pop_scan_beam(CT &y)
{
    if (_scanbeam.empty())
        return false;
    y = _scanbeam.top();
    _scanbeam.pop();
    while (!_scanbeam.empty() && y == _scanbeam.top())
        _scanbeam.pop(); // Pop duplicates.
    return true;
}

template<typename K>
void Clipper<K>::dispose_all_outrecs() 
{
    for (PolyOutList::size_type i = 0; i < _poly_outs.size(); ++i)
        dispose_out_rec(i);
    _poly_outs.clear();
}

template<typename K>
void Clipper<K>::dispose_out_rec(typename PolyOutList::size_type index)
{
    OutRec *outRec = _poly_outs[index];
    if (outRec->pts) 
		dispose_out_pts(outRec->pts);
    delete outRec;
    _poly_outs[index] = 0;
}

template<typename K>
void Clipper<K>::swap_pos_in_ael(TEdge *edge1, TEdge *edge2)
{
    //check that one or other edge hasn't already been removed from AEL ...
    if (edge1->next_in_ael == edge1->prev_in_ael ||
        edge2->next_in_ael == edge2->prev_in_ael)
        return;

    if (edge1->next_in_ael == edge2)
    {
        TEdge *next = edge2->next_in_ael;
        if (next)
            next->prev_in_ael = edge1;
        TEdge *prev = edge1->prev_in_ael;
        if (prev)
            prev->next_in_ael = edge2;
        edge2->prev_in_ael = prev;
        edge2->next_in_ael = edge1;
        edge1->prev_in_ael = edge2;
        edge1->next_in_ael = next;
    }
    else if (edge2->next_in_ael == edge1)
    {
        TEdge *next = edge1->next_in_ael;
        if (next)
            next->prev_in_ael = edge2;
        TEdge *prev = edge2->prev_in_ael;
        if (prev)
            prev->next_in_ael = edge1;
        edge1->prev_in_ael = prev;
        edge1->next_in_ael = edge2;
        edge2->prev_in_ael = edge1;
        edge2->next_in_ael = next;
    }
    else
    {
        TEdge *next = edge1->next_in_ael;
        TEdge *prev = edge1->prev_in_ael;
        edge1->next_in_ael = edge2->next_in_ael;
        if (edge1->next_in_ael)
            edge1->next_in_ael->prev_in_ael = edge1;
        edge1->prev_in_ael = edge2->prev_in_ael;
        if (edge1->prev_in_ael)
            edge1->prev_in_ael->next_in_ael = edge1;
        edge2->next_in_ael = next;
        if (edge2->next_in_ael)
            edge2->next_in_ael->prev_in_ael = edge2;
        edge2->prev_in_ael = prev;
        if (edge2->prev_in_ael)
            edge2->prev_in_ael->next_in_ael = edge2;
    }

    if (!edge1->prev_in_ael)
        _active_edges = edge1;
    else if (!edge2->prev_in_ael)
        _active_edges = edge2;
}

template<typename K>
void Clipper<K>::delete_from_ael(TEdge *e)
{
    TEdge *ael_prev = e->prev_in_ael;
    TEdge *ael_next = e->next_in_ael;
    if (!ael_prev && !ael_next && (e != _active_edges))
		return; //already deleted
    if (ael_prev) 
		ael_prev->next_in_ael = ael_next;
    else 
		_active_edges = ael_next;
    if (ael_next) 
		ael_next->prev_in_ael = ael_prev;
    e->next_in_ael = 0;
    e->prev_in_ael = 0;
}

template<typename K>
typename Clipper<K>::OutRec *
Clipper<K>::create_out_rec()
{
    OutRec *result = new OutRec;
    result->is_hole = false;
    result->is_open = false;
    result->first_left = 0;
    result->pts = 0;
    result->bot_pt = 0;
    result->poly_node = 0;
    _poly_outs.push_back(result);
    result->idx = (int)_poly_outs.size() - 1;
    return result;
}

template<typename K>
void Clipper<K>::update_edge_into_ael(TEdge *&e)
{
    if (!e->next_in_lml)
        throw ClipperException("update_edge_into_ael: invalid call");

    e->next_in_lml->out_idx = e->out_idx;
    TEdge *ael_prev = e->prev_in_ael;
    TEdge *ael_next = e->next_in_ael;
    if (ael_prev)
        ael_prev->next_in_ael = e->next_in_lml;
    else
        _active_edges = e->next_in_lml;
    if (ael_next)
        ael_next->prev_in_ael = e->next_in_lml;
    e->next_in_lml->side = e->side;
    e->next_in_lml->wind_delta = e->wind_delta;
    e->next_in_lml->wind_cnt = e->wind_cnt;
    e->next_in_lml->wind_cnt2 = e->wind_cnt2;
    e = e->next_in_lml;
    e->cur = e->bot();
    e->prev_in_ael = ael_prev;
    e->next_in_ael = ael_next;
    if (!is_horizontal(*e))
        insert_scan_beam(e->top().y());
}

template<typename K>
bool Clipper<K>::local_mimima_pending()
{
    return (_curr_lm != _minimal_list.end());
}

template<typename K>
void Clipper<K>::set_winding_count(TEdge &edge)
{
    TEdge *e = edge.prev_in_ael;
    //find the edge of the same polytype that immediately preceeds 'edge' in AEL
    while (e && ((e->poly_type != edge.poly_type) || (e->wind_delta == 0)))
        e = e->prev_in_ael;
    if (!e)
    {
        if (edge.wind_delta == 0)
        {
            PolyFillType pft = edge.poly_type == PolyType::PT_Subject ? _sub_fill_type : _clip_fill_type;
            edge.wind_cnt = (pft == PolyFillType::PFT_Positive ? 1 : -1);
        }
        else
            edge.wind_cnt = edge.wind_delta;

        edge.wind_cnt2 = 0;
        e = _active_edges; //ie get ready to calc wind_cnt2
    }
    else if (edge.wind_delta == 0 && _clip_type != ClipType::CT_Union)
    {
        edge.wind_cnt = 1;
        edge.wind_cnt2 = e->wind_cnt2;
        e = e->next_in_ael; //ie get ready to calc wind_cnt2
    }
    else if (is_even_odd_fill_type(edge))
    {
        //EvenOdd filling ...
        if (edge.wind_delta == 0)
        {
            //are we inside a subj polygon ...
            bool inside = true;
            TEdge *e2 = e->prev_in_ael;
            while (e2) {
                if (e2->poly_type == e->poly_type && e2->wind_delta != 0)
                    inside = !inside;
                e2 = e2->prev_in_ael;
            }
            edge.wind_cnt = (inside ? 0 : 1);
        }
        else
        {
            edge.wind_cnt = edge.wind_delta;
        }
        edge.wind_cnt2 = e->wind_cnt2;
        e = e->next_in_ael; //ie get ready to calc wind_cnt2
    }
    else
    {
        //nonZero, Positive or Negative filling ...
        if (e->wind_cnt * e->wind_delta < 0)
        {
            //prev edge is 'decreasing' WindCount (WC) toward zero
            //so we're outside the previous polygon ...
            if (std::abs(e->wind_cnt) > 1)
            {
                //outside prev poly but still inside another.
                //when reversing direction of prev poly use the same WC 
                if (e->wind_delta * edge.wind_delta < 0)
                    edge.wind_cnt = e->wind_cnt;
                //otherwise continue to 'decrease' WC ...
                else
                    edge.wind_cnt = e->wind_cnt + edge.wind_delta;
            }
            else
                //now outside all polys of same polytype so set own WC ...
                edge.wind_cnt = (edge.wind_delta == 0 ? 1 : edge.wind_delta);
        }
        else
        {
            //prev edge is 'increasing' WindCount (WC) away from zero
            //so we're inside the previous polygon ...
            if (edge.wind_delta == 0)
                edge.wind_cnt = (e->wind_cnt < 0 ? e->wind_cnt - 1 : e->wind_cnt + 1);
            //if wind direction is reversing prev then use same WC
            else if (e->wind_delta * edge.wind_delta < 0)
                edge.wind_cnt = e->wind_cnt;
            //otherwise add to WC ...
            else
                edge.wind_cnt = e->wind_cnt + edge.wind_delta;
        }
        edge.wind_cnt2 = e->wind_cnt2;
        e = e->next_in_ael; //ie get ready to calc wind_cnt2
    }

    //update wind_cnt2 ...
    if (is_even_odd_alt_fill_type(edge))
    {
        //EvenOdd filling ...
        while (e != &edge)
        {
            if (e->wind_delta != 0)
                edge.wind_cnt2 = (edge.wind_cnt2 == 0 ? 1 : 0);
            e = e->next_in_ael;
        }
    }
    else
    {
        //nonZero, Positive or Negative filling ...
        while (e != &edge)
        {
            edge.wind_cnt2 += e->wind_delta;
            e = e->next_in_ael;
        }
    }
}

template<typename K>
bool Clipper<K>::is_even_odd_fill_type(const TEdge &edge) const
{
    if (edge.poly_type == PolyType::PT_Subject)
        return _sub_fill_type == PolyFillType::PFT_EvenOdd;
    else
        return _clip_fill_type == PolyFillType::PFT_EvenOdd;
}

template<typename K>
bool Clipper<K>::is_even_odd_alt_fill_type(const TEdge &edge) const
{
    if (edge.poly_type == PolyType::PT_Subject)
        return _clip_fill_type == PolyFillType::PFT_EvenOdd;
    else
        return _sub_fill_type == PolyFillType::PFT_EvenOdd;
}

template<typename K>
void Clipper<K>::get_horz_direction(TEdge& horz_edge, Orientation& dir, CT& left, CT& right)
{
	if (horz_edge.bot().x() < horz_edge.top().x())
	{
		left = horz_edge.bot().x();
		right = horz_edge.top().x();
		dir = RIGHT_TURN;
	}
	else
	{
		left = horz_edge.top().x();
		right = horz_edge.bot().x();
		dir = LEFT_TURN;
	}
}

template<typename K>
bool Clipper<K>::get_overlap(const FT a1, const FT a2, const FT b1, const FT b2, FT& left, FT& right)
{
	if (a1 < a2)
	{
		if (b1 < b2) { 
			left = std::max<double>(a1, b1); 
			right = std::min<double>(a2, b2); 
		} else { 
			left = std::max<double>(a1, b2); 
			right = std::min<double>(a2, b1); 
		}
	}
	else
	{
		if (b1 < b2) { 
			left = std::max<double>(a2, b1); 
			right = std::min<double>(a1, b2); 
		}
		else { 
			left = std::max<double>(a2, b2); 
			right = std::min<double>(a1, b1); 
		}
	}
	return left < right;
}

template<typename K>
void Clipper<K>::insert_local_minimal_into_ael(const FT bot_y)
{
    const LocalMinimum *lm;
    while (pop_local_minima(bot_y, lm))
    {
        TEdge *lb = lm->left_bound;
        TEdge *rb = lm->right_bound;

        OutPt *op1 = 0;
        if (!lb)
        {
            //nb: don't insert LB into either AEL or SEL
            insert_edge_into_ael(rb, 0);
            set_winding_count(*rb);
            if (is_contributing(*rb))
                op1 = add_out_pt(rb, rb->bot());
		}
        else if (!rb)
        {
            insert_edge_into_ael(lb, 0);
            set_winding_count(*lb);
            if (is_contributing(*lb))
                op1 = add_out_pt(lb, lb->bot());
            insert_scan_beam(lb->top().y());
        }
        else
        {
            insert_edge_into_ael(lb, 0);
            insert_edge_into_ael(rb, lb);
            set_winding_count(*lb);
            rb->wind_cnt = lb->wind_cnt;
            rb->wind_cnt2 = lb->wind_cnt2;
            if (is_contributing(*lb))
                op1 = add_local_min_poly(lb, rb, lb->bot());
            insert_scan_beam(lb->top().y());
        }

        if (rb)
        {
            if (is_horizontal(*rb))
            {
                add_edge_to_sel(rb);
                if (rb->next_in_lml)
                    insert_scan_beam(rb->next_in_lml->bot().y());
            }
            else
                insert_scan_beam(rb->top().y());
        }

        if (!lb || !rb)
            continue;

        //if any output polygons share an edge, they'll need joining later ...
        if (op1 && is_horizontal(*rb) && _ghost_joints.size() > 0 && (rb->wind_delta != 0))
        {
            for (JoinList::size_type i = 0; i < _ghost_joints.size(); ++i)
            {
                Joint *jr = &_ghost_joints[i];
                //if the horizontal Rb and a 'ghost' horizontal overlap, then convert
                //the 'ghost' join to a real join ready for later ...
                if (horz_segments_overlap(jr->out_pt1->pt.x(), jr->offPt.x(), rb->bot().x(), rb->top().x()))
                    add_join(jr->out_pt1, op1, jr->offPt);
            }
        }

        if (lb->out_idx >= 0 && lb->prev_in_ael && lb->prev_in_ael->cur.x() == lb->bot().x() &&
            lb->prev_in_ael->out_idx >= 0 && slopes_equal(lb->prev_in_ael->bot(), lb->prev_in_ael->top(), lb->cur, lb->top()) &&
            (lb->wind_delta != 0) && (lb->prev_in_ael->wind_delta != 0))
        {
            OutPt *Op2 = add_out_pt(lb->prev_in_ael, lb->bot());
            add_join(op1, Op2, lb->top());
        }

        if (lb->next_in_ael != rb)
        {

            if (rb->out_idx >= 0 && rb->prev_in_ael->out_idx >= 0 &&
                slopes_equal(rb->prev_in_ael->cur, rb->prev_in_ael->top(), rb->cur, rb->top()) &&
                (rb->wind_delta != 0) && (rb->prev_in_ael->wind_delta != 0))
            {
                OutPt *Op2 = add_out_pt(rb->prev_in_ael, rb->bot());
                add_join(op1, Op2, rb->top());
            }

            TEdge *e = lb->next_in_ael;
            if (e)
            {
                while (e != rb)
                {
                    //nb: For calculating winding counts etc, intersect_edges() assumes
                    //that param1 will be to the right of param2 ABOVE the intersection ...
                    intersect_edges(rb, e, lb->cur); //order important here
                    e = e->next_in_ael;
                }
            }
        }
    }
}

template<typename K>
void Clipper<K>::insert_edge_into_ael(TEdge *edge, TEdge *start_edge)
{
    if (!_active_edges)
    {
        edge->prev_in_ael = 0;
        edge->next_in_ael = 0;
        _active_edges = edge;
    }
    else if (!start_edge && e2_insert_before_e1(*_active_edges, *edge))
    {
        edge->prev_in_ael = 0;
        edge->next_in_ael = _active_edges;
        _active_edges->prev_in_ael = edge;
        _active_edges = edge;
    }
    else
    {
        if (!start_edge)
            start_edge = _active_edges;
        while (start_edge->next_in_ael && !e2_insert_before_e1(*start_edge->next_in_ael, *edge))
            start_edge = start_edge->next_in_ael;
        edge->next_in_ael = start_edge->next_in_ael;
        if (start_edge->next_in_ael)
            start_edge->next_in_ael->prev_in_ael = edge;
        edge->prev_in_ael = start_edge;
        start_edge->next_in_ael = edge;
    }
#ifdef DEBUG
	_active_count++;
#endif
}

template<typename K>
void Clipper<K>::add_edge_to_sel(TEdge *edge)
{
    //SEL pointers in PEdge are reused to build a list of horizontal edges.
    //However, we don't need to worry about order with horizontal edge processing.
    if (!_sorted_edges)
    {
        _sorted_edges = edge;
        edge->prev_in_sel = 0;
        edge->next_in_sel = 0;
    }
    else
    {
        edge->next_in_sel = _sorted_edges;
        edge->prev_in_sel = 0;
        _sorted_edges->prev_in_sel = edge;
        _sorted_edges = edge;
    }

#ifdef DEBUG
	_sorted_count++;
#endif
}

template<typename K>
bool Clipper<K>::pop_edge_from_sel(TEdge *&edge)
{
    if (!_sorted_edges)
        return false;
    edge = _sorted_edges;
    delete_from_sel(_sorted_edges);
    return true;
}

template<typename K>
void Clipper<K>::copy_ael_to_sel()
{
    TEdge *e = _active_edges;
    _sorted_edges = e;
    while (e)
    {
        e->prev_in_sel = e->prev_in_ael;
        e->next_in_sel = e->next_in_ael;
        e = e->next_in_ael;
    }
}

template<typename K>
void Clipper<K>::delete_from_sel(TEdge *e)
{
    TEdge *selPrev = e->prev_in_sel;
    TEdge *selNext = e->next_in_sel;
    if (!selPrev && !selNext && (e != _sorted_edges)) 
		return; //already deleted
    if (selPrev) 
		selPrev->next_in_sel = selNext;
    else 
		_sorted_edges = selNext;
    if (selNext) 
		selNext->prev_in_sel = selPrev;
    e->next_in_sel = 0;
    e->prev_in_sel = 0;

#ifdef DEBUG
	_sorted_count--;
#endif
}

#ifdef use_xyz
void Clipper<K>::SetZ(Point2 &pt, TEdge &e1, TEdge &e2)
{
    if (pt.Z != 0 || !m_ZFill) return;
    else if (pt == e1.bot()) pt.Z = e1.bot().Z;
    else if (pt == e1.top()) pt.Z = e1.top().Z;
    else if (pt == e2.bot()) pt.Z = e2.bot().Z;
    else if (pt == e2.top()) pt.Z = e2.top().Z;
    else (*m_ZFill)(e1.bot(), e1.top(), e2.bot(), e2.top(), pt);
}
//------------------------------------------------------------------------------
#endif

template<typename K>
void Clipper<K>::swap_pos_in_sel(TEdge *edge1, TEdge *edge2)
{
    if (!(edge1->next_in_sel) && !(edge1->prev_in_sel))
        return;
    if (!(edge2->next_in_sel) && !(edge2->prev_in_sel))
        return;

    if (edge1->next_in_sel == edge2)
    {
        TEdge *next = edge2->next_in_sel;
        if (next)
            next->prev_in_sel = edge1;
        TEdge *prev = edge1->prev_in_sel;
        if (prev)
            prev->next_in_sel = edge2;
        edge2->prev_in_sel = prev;
        edge2->next_in_sel = edge1;
        edge1->prev_in_sel = edge2;
        edge1->next_in_sel = next;
    }
    else if (edge2->next_in_sel == edge1)
    {
        TEdge *next = edge1->next_in_sel;
        if (next)
            next->prev_in_sel = edge2;
        TEdge *prev = edge2->prev_in_sel;
        if (prev)
            prev->next_in_sel = edge1;
        edge1->prev_in_sel = prev;
        edge1->next_in_sel = edge2;
        edge2->prev_in_sel = edge1;
        edge2->next_in_sel = next;
    }
    else
    {
        TEdge *next = edge1->next_in_sel;
        TEdge *prev = edge1->prev_in_sel;
        edge1->next_in_sel = edge2->next_in_sel;
        if (edge1->next_in_sel)
            edge1->next_in_sel->prev_in_sel = edge1;
        edge1->prev_in_sel =
            edge2->prev_in_sel;
        if (edge1->prev_in_sel)
            edge1->prev_in_sel->next_in_sel = edge1;
        edge2->next_in_sel = next;
        if (edge2->next_in_sel)
            edge2->next_in_sel->prev_in_sel = edge2;
        edge2->prev_in_sel = prev;
        if (edge2->prev_in_sel)
            edge2->prev_in_sel->next_in_sel = edge2;
    }

    if (!edge1->prev_in_sel)
        _sorted_edges = edge1;
    else if (!edge2->prev_in_sel)
        _sorted_edges = edge2;
}

template<typename K>
bool Clipper<K>::is_contributing(const TEdge &edge) const
{
    PolyFillType pft, pft2;
    if (edge.poly_type == PolyType::PT_Subject)
    {
        pft = _sub_fill_type;
        pft2 = _clip_fill_type;
    }
    else
    {
        pft = _clip_fill_type;
        pft2 = _sub_fill_type;
    }

    switch (pft)
    {
    case PolyFillType::PFT_EvenOdd:
        //return false if a subj line has been flagged as inside a subj polygon
        if (edge.wind_delta == 0 && edge.wind_cnt != 1)
            return false;
        break;
    case PolyFillType::PFT_NonZero:
        if (std::abs(edge.wind_cnt) != 1)
            return false;
        break;
    case PolyFillType::PFT_Positive:
        if (edge.wind_cnt != 1)
            return false;
        break;
    default: //PolyFillType::PFT_Negative
        if (edge.wind_cnt != -1)
            return false;
    }

    switch (_clip_type)
    {
    case ClipType::CT_Intersect:
        switch (pft2)
        {
        case PolyFillType::PFT_EvenOdd:
        case PolyFillType::PFT_NonZero:
            return (edge.wind_cnt2 != 0);
        case PolyFillType::PFT_Positive:
            return (edge.wind_cnt2 > 0);
        default:
            return (edge.wind_cnt2 < 0);
        }
        break;
    case ClipType::CT_Union:
        switch (pft2)
        {
        case PolyFillType::PFT_EvenOdd:
        case PolyFillType::PFT_NonZero:
            return (edge.wind_cnt2 == 0);
        case PolyFillType::PFT_Positive:
            return (edge.wind_cnt2 <= 0);
        default:
            return (edge.wind_cnt2 >= 0);
        }
        break;
    case ClipType::CT_Difference:
        if (edge.poly_type == PolyType::PT_Subject)
            switch (pft2)
            {
            case PolyFillType::PFT_EvenOdd:
            case PolyFillType::PFT_NonZero:
                return (edge.wind_cnt2 == 0);
            case PolyFillType::PFT_Positive:
                return (edge.wind_cnt2 <= 0);
            default:
                return (edge.wind_cnt2 >= 0);
            }
        else
            switch (pft2)
            {
            case PolyFillType::PFT_EvenOdd:
            case PolyFillType::PFT_NonZero:
                return (edge.wind_cnt2 != 0);
            case PolyFillType::PFT_Positive:
                return (edge.wind_cnt2 > 0);
            default:
                return (edge.wind_cnt2 < 0);
            }
        break;
    case ClipType::CT_Xor:
        if (edge.wind_delta == 0) //XOr always contributing unless open
            switch (pft2)
            {
            case PolyFillType::PFT_EvenOdd:
            case PolyFillType::PFT_NonZero:
                return (edge.wind_cnt2 == 0);
            case PolyFillType::PFT_Positive:
                return (edge.wind_cnt2 <= 0);
            default:
                return (edge.wind_cnt2 >= 0);
            }
        else
            return true;
        break;
    default:
        return true;
    }
}

template<typename K>
void Clipper<K>::do_maxima(TEdge *e)
{
    TEdge *eMaxPair = get_maxima_pair_ex(e);
    if (!eMaxPair)
    {
        if (e->out_idx >= 0)
            add_out_pt(e, e->top());
        delete_from_ael(e);
        return;
    }

    TEdge *edge_next = e->next_in_ael;
    while (edge_next && edge_next != eMaxPair)
    {
        intersect_edges(e, edge_next, e->top());
        swap_pos_in_ael(e, edge_next);
        edge_next = e->next_in_ael;
    }

    if (e->out_idx == _unassigned && eMaxPair->out_idx == _unassigned)
    {
        delete_from_ael(e);
        delete_from_ael(eMaxPair);
    }
    else if (e->out_idx >= 0 && eMaxPair->out_idx >= 0)
    {
        if (e->out_idx >= 0) 
			add_local_max_poly(e, eMaxPair, e->top());
        delete_from_ael(e);
        delete_from_ael(eMaxPair);
    }
#ifdef USE_LINES
    else if (e->wind_delta == 0)
    {
        if (e->out_idx >= 0)
        {
            add_out_pt(e, e->top());
            e->out_idx = _unassigned;
        }
        delete_from_ael(e);

        if (eMaxPair->out_idx >= 0)
        {
            add_out_pt(eMaxPair, e->top());
            eMaxPair->out_idx = _unassigned;
        }
        delete_from_ael(eMaxPair);
    }
#endif
    else 
		throw ClipperException("do_maxima error");
}

template<typename K>
void Clipper<K>::process_horizontals()
{
    TEdge *horz_edge;
    while (pop_edge_from_sel(horz_edge))
        process_horizontal(horz_edge);
}

template<typename K>
void Clipper<K>::process_horizontal(TEdge *horz_edge)
{
	Orientation dir;
	CT horzLeft, horzRight;
	bool is_open = (horz_edge->wind_delta == 0);

	get_horz_direction(*horz_edge, dir, horzLeft, horzRight);

	TEdge* eLastHorz = horz_edge, * eMaxPair = 0;
	while (eLastHorz->next_in_lml && is_horizontal(*eLastHorz->next_in_lml))
		eLastHorz = eLastHorz->next_in_lml;
	if (!eLastHorz->next_in_lml)
		eMaxPair = get_maxima_pair(eLastHorz);

	MaximaList::const_iterator maxIt;
	MaximaList::const_reverse_iterator maxRit;
	if (_maxima.size() > 0)
	{
		//get the first maxima in range (x()) ...
		if (dir == RIGHT_TURN)
		{
			maxIt = _maxima.begin();
			while (maxIt != _maxima.end() && *maxIt <= horz_edge->bot().x()) 
				maxIt++;
			if (maxIt != _maxima.end() && *maxIt >= eLastHorz->top().x())
				maxIt = _maxima.end();
		}
		else
		{
			maxRit = _maxima.rbegin();
			while (maxRit != _maxima.rend() && *maxRit > horz_edge->bot().x()) 
				maxRit++;
			if (maxRit != _maxima.rend() && *maxRit <= eLastHorz->top().x())
				maxRit = _maxima.rend();
		}
	}

	OutPt* op1 = 0;

	for (;;) //loop through consec. horizontal edges
	{

		bool isLastHorz = (horz_edge == eLastHorz);
		TEdge* e = get_next_in_ael(horz_edge, dir);
		while (e)
		{

			//this code block inserts extra coords into horizontal edges (in output
			//polygons) whereever maxima touch these horizontal edges. This helps
			//'simplifying' polygons (ie if the Simplify property is set).
			if (_maxima.size() > 0)
			{
				if (dir == RIGHT_TURN)
				{
					while (maxIt != _maxima.end() && *maxIt < e->cur.x())
					{
						if (horz_edge->out_idx >= 0 && !is_open)
							add_out_pt(horz_edge, Point2(*maxIt, horz_edge->bot().y()));
						maxIt++;
					}
				}
				else
				{
					while (maxRit != _maxima.rend() && *maxRit > e->cur.x())
					{
						if (horz_edge->out_idx >= 0 && !is_open)
							add_out_pt(horz_edge, Point2(*maxRit, horz_edge->bot().y()));
						maxRit++;
					}
				}
			};

			if ((dir == RIGHT_TURN && e->cur.x() > horzRight) || (dir == LEFT_TURN && e->cur.x() < horzLeft)) 
				break;

			//Also break if we've got to the end of an intermediate horizontal edge ...
			//nb: Smaller dx's are to the right of larger dx's ABOVE the horizontal.
			if (e->cur.x() == horz_edge->top().x() && horz_edge->next_in_lml && e->dx < horz_edge->next_in_lml->dx)
				break;

			if (horz_edge->out_idx >= 0 && !is_open)  //note: may be done multiple times
			{
#ifdef use_xyz
				if (dir == RIGHT_TURN) SetZ(e->cur, *horz_edge, *e);
				else SetZ(e->cur, *e, *horz_edge);
#endif      
				op1 = add_out_pt(horz_edge, e->cur);
				TEdge* eNextHorz = _sorted_edges;
				while (eNextHorz)
				{
					if (eNextHorz->out_idx >= 0 &&
						horz_segments_overlap(horz_edge->bot().x(), horz_edge->top().x(), eNextHorz->bot().x(), eNextHorz->top().x()))
					{
						OutPt* op2 = get_last_outpt(eNextHorz);
						add_join(op2, op1, eNextHorz->top());
					}
					eNextHorz = eNextHorz->next_in_sel;
				}
				add_ghost_join(op1, horz_edge->bot());
			}

			//OK, so far we're still in range of the horizontal Edge  but make sure
					//we're at the last of consec. horizontals when matching with eMaxPair
			if (e == eMaxPair && isLastHorz)
			{
				if (horz_edge->out_idx >= 0)
					add_local_max_poly(horz_edge, eMaxPair, horz_edge->top());
				delete_from_ael(horz_edge);
				delete_from_ael(eMaxPair);
				return;
			}

			if (dir == RIGHT_TURN)
			{
				Point2 pt = Point2(e->cur.x(), horz_edge->cur.y());
				intersect_edges(horz_edge, e, pt);
			}
			else
			{
				Point2 pt = Point2(e->cur.x(), horz_edge->cur.y());
				intersect_edges(e, horz_edge, pt);
			}
			TEdge* edge_next = get_next_in_ael(e, dir);
			swap_pos_in_ael(horz_edge, e);
			e = edge_next;
		} //end while(e)

		//Break out of loop if horz_edge.next_in_lml is not also horizontal ...
		if (!horz_edge->next_in_lml || !is_horizontal(*horz_edge->next_in_lml))
			break;

		update_edge_into_ael(horz_edge);
		if (horz_edge->out_idx >= 0) 
			add_out_pt(horz_edge, horz_edge->bot());
		get_horz_direction(*horz_edge, dir, horzLeft, horzRight);

	} //end for (;;)

	if (horz_edge->out_idx >= 0 && !op1)
	{
		op1 = get_last_outpt(horz_edge);
		TEdge* eNextHorz = _sorted_edges;
		while (eNextHorz)
		{
			if (eNextHorz->out_idx >= 0 &&
				horz_segments_overlap(horz_edge->bot().x(), horz_edge->top().x(), eNextHorz->bot().x(), eNextHorz->top().x()))
			{
				OutPt* op2 = get_last_outpt(eNextHorz);
				add_join(op2, op1, eNextHorz->top());
			}
			eNextHorz = eNextHorz->next_in_sel;
		}
		add_ghost_join(op1, horz_edge->top());
	}

	if (horz_edge->next_in_lml)
	{
		if (horz_edge->out_idx >= 0)
		{
			op1 = add_out_pt(horz_edge, horz_edge->top());
			update_edge_into_ael(horz_edge);
			if (horz_edge->wind_delta == 0)
				return;
			//nb: horz_edge is no longer horizontal here
			TEdge* prev_edge = horz_edge->prev_in_ael;
			TEdge* edge_next = horz_edge->next_in_ael;
			if (prev_edge && prev_edge->cur.x() == horz_edge->bot().x() &&
				prev_edge->cur.y() == horz_edge->bot().y() && prev_edge->wind_delta != 0 &&
				(prev_edge->out_idx >= 0 && prev_edge->cur.y() > prev_edge->top().y() &&
					slopes_equal(*horz_edge, *prev_edge)))
			{
				OutPt* op2 = add_out_pt(prev_edge, horz_edge->bot());
				add_join(op1, op2, horz_edge->top());
			}
			else if (edge_next && edge_next->cur.x() == horz_edge->bot().x() &&
				edge_next->cur.y() == horz_edge->bot().y() && edge_next->wind_delta != 0 &&
				edge_next->out_idx >= 0 && edge_next->cur.y() > edge_next->top().y() &&
				slopes_equal(*horz_edge, *edge_next))
			{
				OutPt* op2 = add_out_pt(edge_next, horz_edge->bot());
				add_join(op1, op2, horz_edge->top());
			}
		}
		else
			update_edge_into_ael(horz_edge);
	}
	else
	{
		if (horz_edge->out_idx >= 0)
			add_out_pt(horz_edge, horz_edge->top());
		delete_from_ael(horz_edge);
	}
}

template<typename K>
typename Clipper<K>::OutPt *
Clipper<K>::add_out_pt(TEdge *e, const Point2 &pt)
{
    if (e->out_idx < 0)
    {
        OutRec *outRec = create_out_rec();
        //outRec->is_open = (e->wind_delta == 0);
		outRec->is_open = (e->wind_delta == _skip);
        OutPt *newOp = new OutPt;
        outRec->pts = newOp;
        newOp->idx = outRec->idx;
        newOp->pt = pt;
        newOp->next = newOp;
        newOp->prev = newOp;
#ifdef DEBUG
		newOp->debug_pt = std::make_pair(pt.x(), pt.y());
#endif
        if (!outRec->is_open)
            set_hole_state(e, outRec);
        e->out_idx = outRec->idx;
        return newOp;
    }
    else
    {
        OutRec *outRec = _poly_outs[e->out_idx];
        //OutRec.pts is the 'left-most' point & OutRec.pts.prev is the 'right-most'
        OutPt *op = outRec->pts;

        bool tofront = (e->side == EdgeSide::esLeft);
        if (tofront && (pt == op->pt))
            return op;
        else if (!tofront && (pt == op->prev->pt))
            return op->prev;

        OutPt *newOp = new OutPt;
        newOp->idx = outRec->idx;
        newOp->pt = pt;
        newOp->next = op;
        newOp->prev = op->prev;
        newOp->prev->next = newOp;
#ifdef DEBUG
		newOp->debug_pt = std::make_pair(pt.x(), pt.y());
#endif
        op->prev = newOp;
        if (tofront)
            outRec->pts = newOp;
        return newOp;
    }
    return nullptr;
}

template<typename K>
typename Clipper<K>::OutPt*
Clipper<K>::get_last_outpt(TEdge* e)
{
	OutRec* outRec = _poly_outs[e->out_idx];
	if (e->side == EdgeSide::esLeft)
		return outRec->pts;
	else
		return outRec->pts->prev;
}

template<typename K>
typename Clipper<K>::OutPt *
Clipper<K>::add_local_min_poly(TEdge *e1, TEdge *e2, const Point2 &pt)
{
    OutPt *result = nullptr;
    TEdge *e, *prev_edge;
    if (is_horizontal(*e2) || (e1->dx > e2->dx))
    {
        result = add_out_pt(e1, pt);
        e2->out_idx = e1->out_idx;
        e1->side = EdgeSide::esRight;
        e2->side = EdgeSide::esLeft;
        e = e1;
        if (e->prev_in_ael == e2)
            prev_edge = e2->prev_in_ael;
        else
            prev_edge = e->prev_in_ael;
    }
    else
    {
        result = add_out_pt(e2, pt);
        e1->out_idx = e2->out_idx;
        e1->side = EdgeSide::esLeft;
        e2->side = EdgeSide::esRight;
        e = e2;
        if (e->prev_in_ael == e1)
            prev_edge = e1->prev_in_ael;
        else
            prev_edge = e->prev_in_ael;
    }

    if (prev_edge && prev_edge->out_idx >= 0 && prev_edge->top().y() < pt.y() && e->top().y() < pt.y())
    {
        CT xPrev = topx(*prev_edge, pt.y());
        CT xE = topx(*e, pt.y());
        if (xPrev == xE && (e->wind_delta != 0) && (prev_edge->wind_delta != 0) &&
            slopes_equal(Point2(xPrev, pt.y()), prev_edge->top(), Point2(xE, pt.y()), e->top()))
        {
            OutPt *outPt = add_out_pt(prev_edge, pt);
            add_join(result, outPt, e->top());
        }
    }
    return result;
}

template<typename K>
typename Clipper<K>::OutPt* 
Clipper<K>::dup_outpt(OutPt* outPt, bool InsertAfter)
{
	OutPt* result = new OutPt;
	result->pt = outPt->pt;
	result->idx = outPt->idx;
#ifdef DEBUG
	result->debug_pt = std::make_pair(outPt->pt.x(), outPt->pt.y());
#endif
	if (InsertAfter)
	{
		result->next = outPt->next;
		result->prev = outPt;
		outPt->next->prev = result;
		outPt->next = result;
	}
	else
	{
		result->prev = outPt->prev;
		result->next = outPt;
		outPt->prev->next = result;
		outPt->prev = result;
	}
	return result;
}

template<typename K>
typename Clipper<K>::OutPt*
Clipper<K>::get_bottom_pt(Clipper<K>::OutPt* pp)
{
	OutPt* dups = 0;
	OutPt* p = pp->next;
	while (p != pp)
	{
		if (p->pt.y() < pp->pt.y())
		{
			pp = p;
			dups = 0;
		}
		else if (p->pt.y() == pp->pt.y() && p->pt.x() <= pp->pt.x())
		{
			if (p->pt.x() < pp->pt.x())
			{
				dups = 0;
				pp = p;
			}
			else
			{
				if (p->next != pp && p->prev != pp) 
					dups = p;
			}
		}
		p = p->next;
	}
	if (dups)
	{
		//there appears to be at least 2 vertices at BottomPt so ...
		while (dups != p)
		{
			if (!first_is_bottom_pt(p, dups))
				pp = dups;
			dups = dups->next;
			while (dups->pt != pp->pt) 
				dups = dups->next;
		}
	}
	return pp;
}

template<typename K>
void Clipper<K>::add_local_max_poly(TEdge *e1, TEdge *e2, const Point2 &pt)
{
    add_out_pt(e1, pt);
    if (e2->wind_delta == 0)
        add_out_pt(e2, pt);
    if (e1->out_idx == e2->out_idx)
    {
        e1->out_idx = _unassigned;
        e2->out_idx = _unassigned;
    }
    else if (e1->out_idx < e2->out_idx)
        append_polygon(e1, e2);
    else
        append_polygon(e2, e1);
}

template<typename K>
typename Clipper<K>::OutRec *
Clipper<K>::get_outrec(int idx)
{
    OutRec* outrec = _poly_outs[idx];
    while (outrec != _poly_outs[outrec->idx])
    	outrec = _poly_outs[outrec->idx];
    return outrec;
}

template<typename K>
typename Clipper<K>::OutRec *
Clipper<K>::get_lowermost_rec(OutRec* outRec1, OutRec* outRec2)
{
	//work out which polygon fragment has the correct hole state ...
	if (!outRec1->bot_pt)
		outRec1->bot_pt = get_bottom_pt(outRec1->pts);
	if (!outRec2->bot_pt)
		outRec2->bot_pt = get_bottom_pt(outRec2->pts);
	OutPt* out_pt1 = outRec1->bot_pt;
	OutPt* out_pt2 = outRec2->bot_pt;
	if (out_pt1->pt.y() > out_pt2->pt.y())
		return outRec1;
	else if (out_pt1->pt.y() < out_pt2->pt.y()) 
		return outRec2;
	else if (out_pt1->pt.x() < out_pt2->pt.x())
		return outRec1;
	else if (out_pt1->pt.x() > out_pt2->pt.x())
		return outRec2;
	else if (out_pt1->next == out_pt1) 
		return outRec2;
	else if (out_pt2->next == out_pt2) 
		return outRec1;
	else if (first_is_bottom_pt(out_pt1, out_pt2)) 
		return outRec1;
	else 
		return outRec2;
}

template<typename K>
void Clipper<K>::skip_outrec(TEdge* e)
{
	if (e->out_idx < 0)
		return;

	if (e->out_idx != _skip)
		return;
}

template<typename K>
void Clipper<K>::append_polygon(TEdge *e1, TEdge *e2)
{
    //get the start and ends of both output polygons ...
    OutRec* outRec1 = _poly_outs[e1->out_idx];
    OutRec* outRec2 = _poly_outs[e2->out_idx];

    OutRec* hole_state_rec;
    if (outrec1_right_of_outrec2(outRec1, outRec2))
    	hole_state_rec = outRec2;
    else if (outrec1_right_of_outrec2(outRec2, outRec1))
    	hole_state_rec = outRec1;
    else
    	hole_state_rec = get_lowermost_rec(outRec1, outRec2);

    //get the start and ends of both output polygons and
    //join e2 poly onto e1 poly and delete pointers to e2 ...

    OutPt* p1_lft = outRec1->pts;
    OutPt* p1_rt = p1_lft->prev;
    OutPt* p2_lft = outRec2->pts;
    OutPt* p2_rt = p2_lft->prev;

    //join e2 poly onto e1 poly and delete pointers to e2 ...
    if (e1->side == EdgeSide::esLeft)
    {
    	if (e2->side == EdgeSide::esLeft)
    	{
    		//z y x a b c
    		reverse_poly_pt_links(p2_lft);
    		p2_lft->next = p1_lft;
    		p1_lft->prev = p2_lft;
    		p1_rt->next = p2_rt;
    		p2_rt->prev = p1_rt;
    		outRec1->pts = p2_rt;
    	}
    	else
    	{
    		//x y z a b c
    		p2_rt->next = p1_lft;
    		p1_lft->prev = p2_rt;
    		p2_lft->prev = p1_rt;
    		p1_rt->next = p2_lft;
    		outRec1->pts = p2_lft;
    	}
    }
    else
    {
    	if (e2->side == EdgeSide::esRight)
    	{
    		//a b c z y x
    		reverse_poly_pt_links(p2_lft);
    		p1_rt->next = p2_rt;
    		p2_rt->prev = p1_rt;
    		p2_lft->next = p1_lft;
    		p1_lft->prev = p2_lft;
    	}
    	else
    	{
    		//a b c x y z
    		p1_rt->next = p2_lft;
    		p2_lft->prev = p1_rt;
    		p1_lft->prev = p2_rt;
    		p2_rt->next = p1_lft;
    	}
    }

    outRec1->bot_pt = 0;
    if (hole_state_rec == outRec2)
    {
    	if (outRec2->first_left != outRec1)
    		outRec1->first_left = outRec2->first_left;
    	outRec1->is_hole = outRec2->is_hole;
    }
    outRec2->pts = 0;
    outRec2->bot_pt = 0;
    outRec2->first_left = outRec1;

    int OKIdx = e1->out_idx;
    int obsolete_idx = e2->out_idx;

    e1->out_idx = _unassigned; //nb: safe because we only get here via add_local_max_poly
    e2->out_idx = _unassigned;

    TEdge* e = _active_edges;
    while (e)
    {
    	if (e->out_idx == obsolete_idx)
    	{
    		e->out_idx = OKIdx;
    		e->side = e1->side;
    		break;
    	}
    	e = e->next_in_ael;
    }

    outRec2->idx = outRec1->idx;
}

template<typename K>
void Clipper<K>::intersect_edges(TEdge *e1, TEdge *e2, Point2 &pt)
{
    bool e1_contribute = (e1->out_idx >= 0);
    bool e2_contribute = (e2->out_idx >= 0);

#ifdef use_xyz
    SetZ(pt, *e1, *e2);
#endif

#ifdef USE_LINES
    //if either edge is on an OPEN path ...
    if (e1->wind_delta == 0 || e2->wind_delta == 0)
    {
        //ignore subject-subject open path intersections UNLESS they
        //are both open paths, AND they are both 'contributing maximas' ...
        if (e1->wind_delta == 0 && e2->wind_delta == 0)
            return;

        //if intersecting a subj line with a subj poly ...
        else if (e1->poly_type == e2->poly_type &&
            e1->wind_delta != e2->wind_delta && _clip_type == ClipType::CT_Union)
        {
            if (e1->wind_delta == 0)
            {
                if (e2_contribute)
                {
                    add_out_pt(e1, pt);
                    if (e1_contribute)
                        e1->out_idx = _unassigned;
                }
            }
            else
            {
                if (e1_contribute)
                {
                    add_out_pt(e2, pt);
                    if (e2_contribute)
                        e2->out_idx = _unassigned;
                }
            }
        }
        else if (e1->poly_type != e2->poly_type)
        {
            //toggle subj open path out_idx on/off when std::abs(clip.WndCnt) == 1 ...
            if ((e1->wind_delta == 0) && abs(e2->wind_cnt) == 1 &&
                (_clip_type != ClipType::CT_Union || e2->wind_cnt2 == 0))
            {
                add_out_pt(e1, pt);
                if (e1_contribute)
                    e1->out_idx = _unassigned;
            }
            else if ((e2->wind_delta == 0) && (abs(e1->wind_cnt) == 1) &&
                (_clip_type != ClipType::CT_Union || e1->wind_cnt2 == 0))
            {
                add_out_pt(e2, pt);
                if (e2_contribute)
                    e2->out_idx = _unassigned;
            }
        }
        return;
    }
#endif

    //update winding counts...
    //assumes that e1 will be to the right of e2 ABOVE the intersection
    if (e1->poly_type == e2->poly_type) {
        if (is_even_odd_fill_type(*e1))
        {
            int old_e1_windcnt = e1->wind_cnt;
            e1->wind_cnt = e2->wind_cnt;
            e2->wind_cnt = old_e1_windcnt;
        }
        else
        {
            if (e1->wind_cnt + e2->wind_delta == 0)
                e1->wind_cnt = -e1->wind_cnt;
            else
                e1->wind_cnt += e2->wind_delta;

            if (e2->wind_cnt - e1->wind_delta == 0)
                e2->wind_cnt = -e2->wind_cnt;
            else
                e2->wind_cnt -= e1->wind_delta;
        }
    }
    else {
        if (is_even_odd_fill_type(*e2))
            e1->wind_cnt2 = (e1->wind_cnt2 == 0) ? 1 : 0;
        else
            e1->wind_cnt2 += e2->wind_delta;

        if (is_even_odd_fill_type(*e1))
            e2->wind_cnt2 = (e2->wind_cnt2 == 0) ? 1 : 0;
        else
            e2->wind_cnt2 -= e1->wind_delta;
    }

    PolyFillType e1FillType, e2FillType,
        e1FillType2, e2FillType2;
    if (e1->poly_type == PolyType::PT_Subject)
    {
        e1FillType = _sub_fill_type;
        e1FillType2 = _clip_fill_type;
    }
    else
    {
        e1FillType = _clip_fill_type;
        e1FillType2 = _sub_fill_type;
    }
    if (e2->poly_type == PolyType::PT_Subject)
    {
        e2FillType = _sub_fill_type;
        e2FillType2 = _clip_fill_type;
    }
    else
    {
        e2FillType = _clip_fill_type;
        e2FillType2 = _sub_fill_type;
    }

    int e1_wc, e2_wc;
    switch (e1FillType)
    {
    case PolyFillType::PFT_Positive:
        e1_wc = e1->wind_cnt;
        break;
    case PolyFillType::PFT_Negative:
        e1_wc = -e1->wind_cnt;
        break;
    default:
        e1_wc = std::abs(e1->wind_cnt);
    }
    switch (e2FillType)
    {
    case PolyFillType::PFT_Positive:
        e2_wc = e2->wind_cnt; 
		break;
    case PolyFillType::PFT_Negative:
        e2_wc = -e2->wind_cnt;
		break;
    default:
        e2_wc = std::abs(e2->wind_cnt);
    }

    if (e1_contribute && e2_contribute)
    {
        if ((e1_wc != 0 && e1_wc != 1) || (e2_wc != 0 && e2_wc != 1) || (e1->poly_type != e2->poly_type && _clip_type != ClipType::CT_Xor)) 
        {
            add_local_max_poly(e1, e2, pt);
        }
        else
        {
            add_out_pt(e1, pt);
            add_out_pt(e2, pt);
            swap_sides(*e1, *e2);
            swap_poly_indexs(*e1, *e2);
        }
    }
    else if (e1_contribute)
    {
        if (e2_wc == 0 || e2_wc == 1)
        {
            add_out_pt(e1, pt);
            swap_sides(*e1, *e2);
            swap_poly_indexs(*e1, *e2);
        }
    }
    else if (e2_contribute)
    {
        if (e1_wc == 0 || e1_wc == 1)
        {
            add_out_pt(e2, pt);
            swap_sides(*e1, *e2);
            swap_poly_indexs(*e1, *e2);
        }
    }
    else if ((e1_wc == 0 || e1_wc == 1) && (e2_wc == 0 || e2_wc == 1)) {
        if (e1->poly_type != e2->poly_type)
        {
            add_local_min_poly(e1, e2, pt);
        }
		else if (e1_wc == 1 && e2_wc == 1) {
			CT e1_wc2, e2_wc2;
			switch (e1FillType2)
			{
			case PolyFillType::PFT_Positive:
				e1_wc2 = e1->wind_cnt2;
				break;
			case PolyFillType::PFT_Negative:
				e1_wc2 = -e1->wind_cnt2;
				break;
			default:
				e1_wc2 = std::abs(e1->wind_cnt2);
			}
			switch (e2FillType2)
			{
			case PolyFillType::PFT_Positive:
				e2_wc2 = e2->wind_cnt2;
				break;
			case PolyFillType::PFT_Negative:
				e2_wc2 = -e2->wind_cnt2;
				break;
			default:
				e2_wc2 = std::abs(e2->wind_cnt2);
			}

			switch (_clip_type) {
			case ClipType::CT_Intersect:
				if (e1_wc2 > 0 && e2_wc2 > 0)
					add_local_min_poly(e1, e2, pt);
				break;
			case ClipType::CT_Union:
				if (e1_wc2 <= 0 && e2_wc2 <= 0)
					add_local_min_poly(e1, e2, pt);
				break;
			case ClipType::CT_Difference:
				if (((e1->poly_type == PolyType::PT_Clip) && (e1_wc2 > 0) && (e2_wc2 > 0)) ||
					((e1->poly_type == PolyType::PT_Subject) && (e1_wc2 <= 0) && (e2_wc2 <= 0)))
					add_local_min_poly(e1, e2, pt);
				break;
			case ClipType::CT_Xor:
				add_local_min_poly(e1, e2, pt);
			}
		} else
            swap_sides(*e1, *e2);
    }
}

template<typename K>
bool Clipper<K>::process_intersections(const CT top_y)
{
    if (!_active_edges)
        return true;
    try {
        build_intersect_list(top_y);
        size_t ilSize = _intersect_list.size();
        if (ilSize == 0)
            return true;
        if (ilSize == 1 || fixup_intersection_order())
            process_intersect_list();
        else
            return false;
    }
    catch (...)
    {
        _sorted_edges = 0;
        dispose_intersect_nodes();
        throw ClipperException("process_intersections error");
    }
    _sorted_edges = 0;
    return true;
}

template<typename K>
void Clipper<K>::build_intersect_list(const CT top_y)
{
    if (!_active_edges)
        return;

    //prepare for sorting ...
    TEdge *e = _active_edges;
    _sorted_edges = e;
    while (e)
    {
        e->prev_in_sel = e->prev_in_ael;
        e->next_in_sel = e->next_in_ael;
        e->cur = Point2(topx(*e, top_y), e->cur.y());
        e = e->next_in_ael;
    }

    //bubblesort ...
    bool is_modified;
    do
    {
        is_modified = false;
        e = _sorted_edges;
        while (e->next_in_sel)
        {
            TEdge *edge_next = e->next_in_sel;
            if (e->cur.x() > edge_next->cur.x())
            {
                Point2 pt = intersect_point(*e, *edge_next);
                if (pt.y() > top_y)
                    pt = Point2(topx(*e, top_y), top_y);
                IntersectNode new_node;
                new_node.edge1 = e;
                new_node.edge2 = edge_next;
                new_node.pt = pt;
                _intersect_list.push_back(new_node);

                swap_pos_in_sel(e, edge_next);
                is_modified = true;
            }
            else
                e = edge_next;
        }
        if (e->prev_in_sel)
            e->prev_in_sel->next_in_sel = 0;
        else
            break;
    } while (is_modified);
    _sorted_edges = 0; //important
}

template<typename K>
void Clipper<K>::process_intersect_list()
{
    for (size_t i = 0; i < _intersect_list.size(); ++i)
    {
        IntersectNode *iNode = &_intersect_list[i];
        {
            intersect_edges(iNode->edge1, iNode->edge2, iNode->pt);
            swap_pos_in_ael(iNode->edge1, iNode->edge2);
        }
    }
    _intersect_list.clear();
}

template<typename K>
void Clipper<K>::process_top_edges_scanbeam(const CT top_y)
{
    TEdge *e = _active_edges;
    while (e)
    {
        //1. process maxima, treating them as if they're 'bent' horizontal edges,
        //   but exclude maxima with horizontal edges. nb: e can't be a horizontal.
        bool is_max = is_maxima(e, top_y);

        if (is_max)
        {
            TEdge *eMaxPair = get_maxima_pair_ex(e);
            is_max = (!eMaxPair || !is_horizontal(*eMaxPair));
        }

        if (is_max)
        {
			if (e->out_idx >= 0 && e->next && e->next->out_idx == _skip) 
				_poly_outs[e->out_idx]->is_open = true;

            if (_strict_simple)
                _maxima.push_back(e->top().x());

            TEdge *prev_edge = e->prev_in_ael;
            do_maxima(e);
            if (!prev_edge)
                e = _active_edges;
            else
                e = prev_edge->next_in_ael;
        }
        else
        {
            //2. promote horizontal edges, otherwise update cur.x() and cur.y() ...
            if (is_intermediate(e, top_y) && is_horizontal(*e->next_in_lml))
            {
                update_edge_into_ael(e);
                if (e->out_idx >= 0)
                    add_out_pt(e, e->bot());
                add_edge_to_sel(e);
            }
            else
            {
                e->cur = Point2(topx(*e, top_y), top_y);
#ifdef use_xyz
                e->cur.Z = top_y == e->top().y() ? e->top().Z : (top_y == e->bot().y() ? e->bot().Z : 0);
#endif
            }

            //When strictlySimple and 'e' is being touched by another edge, then
            //make sure both edges have a vertex here ...
            if (_strict_simple)
            {
                TEdge *prev_edge = e->prev_in_ael;
                if ((e->out_idx >= 0) && (e->wind_delta != 0) && prev_edge && (prev_edge->out_idx >= 0) &&
                    (prev_edge->cur.x() == e->cur.x()) && (prev_edge->wind_delta != 0))
                {
                    Point2 pt = e->cur;
#ifdef use_xyz
                    SetZ(pt, *prev_edge, *e);
#endif
                    OutPt *op = add_out_pt(prev_edge, pt);
                    OutPt *op2 = add_out_pt(e, pt);
                    add_join(op, op2, pt); //strictlySimple (type-3) join
                }
            }

            e = e->next_in_ael;
        }
    }

    //3. Process horizontals at the top() of the scanbeam ...
    _maxima.sort();
    process_horizontals();
    _maxima.clear();

    //4. Promote intermediate vertices ...
    e = _active_edges;
    while (e)
    {
        if (is_intermediate(e, top_y))
        {
            OutPt *op = 0;
            if (e->out_idx >= 0)
                op = add_out_pt(e, e->top());
            update_edge_into_ael(e);

            //if output polygons share an edge, they'll need joining later ...
            TEdge *prev_edge = e->prev_in_ael;
            TEdge *edge_next = e->next_in_ael;
            if (prev_edge && prev_edge->cur.x() == e->bot().x() &&
                prev_edge->cur.y() == e->bot().y() && op &&
                prev_edge->out_idx >= 0 && prev_edge->cur.y() > prev_edge->top().y() &&
                slopes_equal(e->cur, e->top(), prev_edge->cur, prev_edge->top()) &&
                (e->wind_delta != 0) && (prev_edge->wind_delta != 0))
            {
                OutPt *op2 = add_out_pt(prev_edge, e->bot());
                add_join(op, op2, e->top());
            }
            else if (edge_next && edge_next->cur.x() == e->bot().x() &&
                edge_next->cur.y() == e->bot().y() && op &&
                edge_next->out_idx >= 0 && edge_next->cur.y() > edge_next->top().y() &&
                slopes_equal(e->cur, e->top(), edge_next->cur, edge_next->top()) &&
                (e->wind_delta != 0) && (edge_next->wind_delta != 0))
            {
                OutPt *op2 = add_out_pt(edge_next, e->bot());
                add_join(op, op2, e->top());
            }
        }
        e = e->next_in_ael;
    }
}

template<typename K>
void Clipper<K>::build_result(Paths &polys)
{
    polys.reserve(_poly_outs.size());
    for (PolyOutList::size_type i = 0; i < _poly_outs.size(); ++i)
    {
        if (!_poly_outs[i]->pts)
            continue;
        OutPt *p = _poly_outs[i]->pts->next;
        int cnt = point_count(p);
        if (cnt < 2)
            continue;
        Path pg;
        pg.reserve(cnt);
        for (int i = 0; i < cnt; ++i)
        {
            pg.push_back(p->pt);
            p = p->next;
        }
        polys.push_back(pg);
    }
}

template<typename K>
void Clipper<K>::build_result2(PolyTree<K>& polytree)
{
	polytree.clear();
	polytree._all_nodes.reserve(_poly_outs.size());
	//add each output polygon/contour to polytree ...
	for (PolyOutList::size_type i = 0; i < _poly_outs.size(); i++)
	{
		OutRec* outRec = _poly_outs[i];
		int cnt = point_count(outRec->pts);
		if ((outRec->is_open && cnt < 2) || (!outRec->is_open && cnt < 3)) 
			continue;
		fix_hole_linkage(*outRec);
		PolyNode<K> * pn = new PolyNode<K>();
		//nb: polytree takes ownership of all the PolyNodes
		polytree._all_nodes.push_back(pn);
		outRec->poly_node = pn;
		pn->_parent = 0;
		pn->_index = 0;
		pn->_contour.reserve(cnt);
		OutPt* op = outRec->pts->next;
		for (int j = 0; j < cnt; j++)
		{
			pn->_contour.push_back(op->pt);
			op = op->next;
		}
	}

	//fixup PolyNode links etc ...
	polytree._childs.reserve(_poly_outs.size());
	for (PolyOutList::size_type i = 0; i < _poly_outs.size(); i++)
	{
		OutRec* outRec = _poly_outs[i];
		if (!outRec->poly_node) 
			continue;
		if (outRec->is_open)
		{
			outRec->poly_node->_is_open = true;
			polytree.add_child(outRec->poly_node);
		}
		else if (outRec->first_left && outRec->first_left->poly_node)
			outRec->first_left->poly_node->add_child(outRec->poly_node);
		else
			polytree.add_child(outRec->poly_node);
	}
}

template<typename K>
void Clipper<K>::set_hole_state(TEdge *e, OutRec *outrec)
{
    TEdge *e2 = e->prev_in_ael;
    TEdge *eTmp = 0;
    while (e2)
    {
        if (e2->out_idx >= 0 && e2->wind_delta != 0)
        {
            if (!eTmp)
                eTmp = e2;
            else if (eTmp->out_idx == e2->out_idx)
                eTmp = 0;
        }
        e2 = e2->prev_in_ael;
    }
    if (!eTmp)
    {
        outrec->first_left = 0;
        outrec->is_hole = false;
    }
    else
    {
        outrec->first_left = _poly_outs[eTmp->out_idx];
        outrec->is_hole = !outrec->first_left->is_hole;
    }
}

template<typename K>
void Clipper<K>::dispose_intersect_nodes()
{
    _intersect_list.clear();
}

template<typename K>
bool Clipper<K>::fixup_intersection_order()
{
    //pre-condition: intersections are sorted Bottom-most first.
    //Now it's crucial that intersections are made only between adjacent edges,
    //so to ensure this the order of intersections may need adjusting ...
    copy_ael_to_sel();

    std::sort(_intersect_list.begin(), _intersect_list.end(),
        [](const IntersectNode &node1, const IntersectNode &node2) {
            return node2.pt.y() < node1.pt.y();
        });

    size_t cnt = _intersect_list.size();
    for (size_t i = 0; i < cnt; ++i)
    {
        if (!edges_adjacent(_intersect_list[i]))
        {
            size_t j = i + 1;
            while (j < cnt && !edges_adjacent(_intersect_list[j]))
                j++;
            if (j == cnt)
                return false;
            std::swap(_intersect_list[i], _intersect_list[j]);
        }
        swap_pos_in_sel(_intersect_list[i].edge1, _intersect_list[i].edge2);
    }
    return true;
}

template<typename K>
void Clipper<K>::fixup_out_polygon(OutRec &outrec)
{
    //fixup_out_polygon() - removes duplicate points and simplifies consecutive
    //parallel edges by removing the middle vertex.
    OutPt *lastOK = 0;
    outrec.bot_pt = 0;
    OutPt *pp = outrec.pts;
    bool preserveCol = _preserve_collinear || _strict_simple;

    for (;;)
    {
        if (pp->prev == pp || pp->prev == pp->next)
        {
            dispose_out_pts(pp);
            outrec.pts = 0;
            return;
        }

        //test for duplicate points and collinear edges ...
        if ((pp->pt == pp->next->pt) || (pp->pt == pp->prev->pt) || (slopes_equal(pp->prev->pt, pp->pt, pp->next->pt) &&
                (!preserveCol || !pt2_between_pt1_pt3(pp->prev->pt, pp->pt, pp->next->pt))))
        {
            lastOK = 0;
            OutPt *tmp = pp;
            pp->prev->next = pp->next;
            pp->next->prev = pp->prev;
            pp = pp->prev;
            delete tmp;
        }
        else if (pp == lastOK)
			break;
        else
        {
            if (!lastOK) 
				lastOK = pp;
            pp = pp->next;
        }
    }
    outrec.pts = pp;
}

template<typename K>
void Clipper<K>::fixup_out_polyline(OutRec &outrec)
{
    OutPt *pp = outrec.pts;
    OutPt *lastPP = pp->prev;
    while (pp != lastPP)
    {
        pp = pp->next;
        if (pp->pt == pp->prev->pt)
        {
            if (pp == lastPP)
                lastPP = pp->prev;
            OutPt *tmpPP = pp->prev;
            tmpPP->next = pp->next;
            pp->next->prev = tmpPP;
            delete pp;
            pp = tmpPP;
        }
    }

    if (pp == pp->prev)
    {
        dispose_out_pts(pp);
        outrec.pts = 0;
        return;
    }
}

template<typename K>
bool Clipper<K>::outrec1_right_of_outrec2(OutRec* outRec1, OutRec* outRec2)
{
	do
	{
		outRec1 = outRec1->first_left;
		if (outRec1 == outRec2)
			return true;
	} while (outRec1);
	return false;
}

template<typename K>
void Clipper<K>::fix_hole_linkage(OutRec &outrec)
{
    //skip OutRecs that (a) contain outermost polygons or
    //(b) already have the correct owner/child linkage ...
    if (!outrec.first_left || (outrec.is_hole != outrec.first_left->is_hole && outrec.first_left->pts)) 
		return;

    OutRec *orfl = outrec.first_left;
    while (orfl && ((orfl->is_hole == outrec.is_hole) || !orfl->pts))
        orfl = orfl->first_left;
    outrec.first_left = orfl;
}

template<typename K>
void Clipper<K>::add_join(OutPt *op1, OutPt *op2, const Point2 offPt)
{
    Joint j;
    j.out_pt1 = op1;
    j.out_pt2 = op2;
    j.offPt = offPt;
    _joints.push_back(j);
}

template<typename K>
void Clipper<K>::clear_joints()
{
    _joints.clear();
}

template<typename K>
void Clipper<K>::clear_ghost_joints()
{
    _ghost_joints.clear();
}

template<typename K>
void Clipper<K>::add_ghost_join(OutPt *op, const Point2 offPt)
{
    Joint j;
    j.out_pt1 = op;
    j.out_pt2 = 0;
    j.offPt = offPt;
    _ghost_joints.push_back(j);
}

template<typename K>
bool Clipper<K>::join_horz(OutPt* op1, OutPt* op1b, OutPt* op2, OutPt* op2b, const Point2 pt, bool discardLeft)
{
	Orientation Dir1 = (op1->pt.x() > op1b->pt.x() ? LEFT_TURN : RIGHT_TURN);
	Orientation Dir2 = (op2->pt.x() > op2b->pt.x() ? LEFT_TURN : RIGHT_TURN);
	if (Dir1 == Dir2) return false;

	//When discardLeft, we want Op1b to be on the left of Op1, otherwise we
	//want Op1b to be on the right. (And likewise with Op2 and Op2b.)
	//So, to facilitate this while inserting Op1b and Op2b ...
	//when discardLeft, make sure we're AT or RIGHT of pt before adding Op1b,
	//otherwise make sure we're AT or LEFT of pt. (Likewise with Op2b.)
	if (Dir1 == RIGHT_TURN)
	{
		while (op1->next->pt.x() <= pt.x() && op1->next->pt.x() >= op1->pt.x() && op1->next->pt.y() == pt.y())
			op1 = op1->next;
		if (discardLeft && (op1->pt.x() != pt.x())) 
			op1 = op1->next;
		op1b = dup_outpt(op1, !discardLeft);
		if (op1b->pt != pt)
		{
			op1 = op1b;
			op1->pt = pt;
			op1b = dup_outpt(op1, !discardLeft);
		}
	}
	else
	{
		while (op1->next->pt.x() >= pt.x() && op1->next->pt.x() <= op1->pt.x() && op1->next->pt.y() == pt.y())
			op1 = op1->next;
		if (!discardLeft && (op1->pt.x() != pt.x()))
			op1 = op1->next;
		op1b = dup_outpt(op1, discardLeft);
		if (op1b->pt != pt)
		{
			op1 = op1b;
			op1->pt = pt;
			op1b = dup_outpt(op1, discardLeft);
		}
	}

	if (Dir2 == RIGHT_TURN)
	{
		while (op2->next->pt.x() <= pt.x() && op2->next->pt.x() >= op2->pt.x() && op2->next->pt.y() == pt.y())
			op2 = op2->next;
		if (discardLeft && (op2->pt.x() != pt.x())) 
			op2 = op2->next;
		op2b = dup_outpt(op2, !discardLeft);
		if (op2b->pt != pt)
		{
			op2 = op2b;
			op2->pt = pt;
			op2b = dup_outpt(op2, !discardLeft);
		};
	}
	else
	{
		while (op2->next->pt.x() >= pt.x() && op2->next->pt.x() <= op2->pt.x() && op2->next->pt.y() == pt.y())
			op2 = op2->next;
		if (!discardLeft && (op2->pt.x() != pt.x())) 
			op2 = op2->next;
		op2b = dup_outpt(op2, discardLeft);
		if (op2b->pt != pt)
		{
			op2 = op2b;
			op2->pt = pt;
			op2b = dup_outpt(op2, discardLeft);
		};
	};

	if ((Dir1 == RIGHT_TURN) == discardLeft)
	{
		op1->prev = op2;
		op2->next = op1;
		op1b->next = op2b;
		op2b->prev = op1b;
	}
	else
	{
		op1->next = op2;
		op2->prev = op1;
		op1b->prev = op2b;
		op2b->next = op1b;
	}
	return true;
}

template<typename K>
bool Clipper<K>::join_points(Joint* j, OutRec* outRec1, OutRec* outRec2)
{
	OutPt* op1 = j->out_pt1, * op1b;
	OutPt* op2 = j->out_pt2, * op2b;

	//There are 3 kinds of joins for output polygons ...
	//1. Horizontal joins where Joint.out_pt1 & Joint.out_pt2 are vertices anywhere
	//along (horizontal) collinear edges (& Joint.offPt is on the same horizontal).
	//2. Non-horizontal joins where Joint.out_pt1 & Joint.out_pt2 are at the same
	//location at the Bottom of the overlapping segment (& Joint.offPt is above).
	//3. StrictSimple joins where edges touch but are not collinear and where
	//Joint.out_pt1, Joint.out_pt2 & Joint.offPt all share the same point.
	bool is_horizontal = (j->out_pt1->pt.y() == j->offPt.y());

	if (is_horizontal && (j->offPt == j->out_pt1->pt) && (j->offPt == j->out_pt2->pt))
	{
		//Strictly Simple join ...
		if (outRec1 != outRec2)
			return false;
		op1b = j->out_pt1->next;
		while (op1b != op1 && (op1b->pt == j->offPt))
			op1b = op1b->next;
		bool reverse1 = (op1b->pt.y() > j->offPt.y());
		op2b = j->out_pt2->next;
		while (op2b != op2 && (op2b->pt == j->offPt))
			op2b = op2b->next;
		bool reverse2 = (op2b->pt.y() > j->offPt.y());
		if (reverse1 == reverse2)
			return false;
		if (reverse1)
		{
			op1b = dup_outpt(op1, false);
			op2b = dup_outpt(op2, true);
			op1->prev = op2;
			op2->next = op1;
			op1b->next = op2b;
			op2b->prev = op1b;
			j->out_pt1 = op1;
			j->out_pt2 = op1b;
			return true;
		}
		else
		{
			op1b = dup_outpt(op1, true);
			op2b = dup_outpt(op2, false);
			op1->next = op2;
			op2->prev = op1;
			op1b->prev = op2b;
			op2b->next = op1b;
			j->out_pt1 = op1;
			j->out_pt2 = op1b;
			return true;
		}
	}
	else if (is_horizontal)
	{
		//treat horizontal joins differently to non-horizontal joins since with
		//them we're not yet sure where the overlapping is. out_pt1.pt & out_pt2.pt
		//may be anywhere along the horizontal edge.
		op1b = op1;
		while (op1->prev->pt.y() == op1->pt.y() && op1->prev != op1b && op1->prev != op2)
			op1 = op1->prev;
		while (op1b->next->pt.y() == op1b->pt.y() && op1b->next != op1 && op1b->next != op2)
			op1b = op1b->next;
		if (op1b->next == op1 || op1b->next == op2) 
			return false; //a flat 'polygon'

		op2b = op2;
		while (op2->prev->pt.y() == op2->pt.y() && op2->prev != op2b && op2->prev != op1b)
			op2 = op2->prev;
		while (op2b->next->pt.y() == op2b->pt.y() && op2b->next != op2 && op2b->next != op1)
			op2b = op2b->next;
		if (op2b->next == op2 || op2b->next == op1) 
			return false; //a flat 'polygon'

		CT left, right;
		//Op1 --> Op1b & Op2 --> Op2b are the extremites of the horizontal edges
		if (!get_overlap(op1->pt.x(), op1b->pt.x(), op2->pt.x(), op2b->pt.x(), left, right))
			return false;

		//discardLeftSide: when overlapping edges are joined, a spike will created
		//which needs to be cleaned up. However, we don't want Op1 or Op2 caught up
		//on the discard side as either may still be needed for other joins ...
		Point2 pt;
		bool discardLeftSide;
		if (op1->pt.x() >= left && op1->pt.x() <= right)
		{
			pt = op1->pt; discardLeftSide = (op1->pt.x() > op1b->pt.x());
		}
		else if (op2->pt.x() >= left && op2->pt.x() <= right)
		{
			pt = op2->pt; discardLeftSide = (op2->pt.x() > op2b->pt.x());
		}
		else if (op1b->pt.x() >= left && op1b->pt.x() <= right)
		{
			pt = op1b->pt; discardLeftSide = op1b->pt.x() > op1->pt.x();
		}
		else
		{
			pt = op2b->pt; discardLeftSide = (op2b->pt.x() > op2->pt.x());
		}
		j->out_pt1 = op1; j->out_pt2 = op2;
		return join_horz(op1, op1b, op2, op2b, pt, discardLeftSide);
	}
	else
	{
		//nb: For non-horizontal joins ...
		//    1. Jr.out_pt1.pt.y() == Jr.out_pt2.pt.y()
		//    2. Jr.out_pt1.pt > Jr.offPt.y()

		//make sure the polygons are correctly oriented ...
		op1b = op1->next;
		while ((op1b->pt == op1->pt) && (op1b != op1)) 
			op1b = op1b->next;
		bool Reverse1 = ((op1b->pt.y() > op1->pt.y()) || !slopes_equal(op1->pt, op1b->pt, j->offPt));
		if (Reverse1)
		{
			op1b = op1->prev;
			while ((op1b->pt == op1->pt) && (op1b != op1)) 
				op1b = op1b->prev;
			if ((op1b->pt.y() > op1->pt.y()) || !slopes_equal(op1->pt, op1b->pt, j->offPt))
				return false;
		};
		op2b = op2->next;
		while ((op2b->pt == op2->pt) && (op2b != op2))
			op2b = op2b->next;
		bool Reverse2 = ((op2b->pt.y() > op2->pt.y()) || !slopes_equal(op2->pt, op2b->pt, j->offPt));
		if (Reverse2)
		{
			op2b = op2->prev;
			while ((op2b->pt == op2->pt) && (op2b != op2))
				op2b = op2b->prev;
			if ((op2b->pt.y() > op2->pt.y()) || !slopes_equal(op2->pt, op2b->pt, j->offPt)) 
				return false;
		}

		if ((op1b == op1) || (op2b == op2) || (op1b == op2b) || ((outRec1 == outRec2) && (Reverse1 == Reverse2)))
			return false;

		if (Reverse1)
		{
			op1b = dup_outpt(op1, false);
			op2b = dup_outpt(op2, true);
			op1->prev = op2;
			op2->next = op1;
			op1b->next = op2b;
			op2b->prev = op1b;
			j->out_pt1 = op1;
			j->out_pt2 = op1b;
			return true;
		}
		else
		{
			op1b = dup_outpt(op1, true);
			op2b = dup_outpt(op2, false);
			op1->next = op2;
			op2->prev = op1;
			op1b->prev = op2b;
			op2b->next = op1b;
			j->out_pt1 = op1;
			j->out_pt2 = op1b;
			return true;
		}
	}
	return false;
}

template<typename K>
void Clipper<K>::join_common_edges()
{
	for (JoinList::size_type i = 0; i < _joints.size(); i++)
	{
		Joint* join = &_joints[i];

		OutRec* outRec1 = get_outrec(join->out_pt1->idx);
		OutRec* outRec2 = get_outrec(join->out_pt2->idx);

		if (!outRec1->pts || !outRec2->pts) 
			continue;
		if (outRec1->is_open || outRec2->is_open)
			continue;

		//get the polygon fragment with the correct hole state (first_left)
		//before calling join_points() ...
		OutRec* hole_state_rec;
		if (outRec1 == outRec2) 
			hole_state_rec = outRec1;
		else if (outrec1_right_of_outrec2(outRec1, outRec2)) 
			hole_state_rec = outRec2;
		else if (outrec1_right_of_outrec2(outRec2, outRec1))
			hole_state_rec = outRec1;
		else 
			hole_state_rec = get_lowermost_rec(outRec1, outRec2);

		if (!join_points(join, outRec1, outRec2))
			continue;

		if (outRec1 == outRec2)
		{
			//instead of joining two polygons, we've just created a new one by
			//splitting one polygon into two.
			outRec1->pts = join->out_pt1;
			outRec1->bot_pt = 0;
			outRec2 = create_out_rec();
			outRec2->pts = join->out_pt2;

			//update all OutRec2.pts idx's ...
			update_outpt_idxs(*outRec2);

			if (poly2_contains_poly1(outRec2->pts, outRec1->pts))
			{
				//outRec1 contains outRec2 ...
				outRec2->is_hole = !outRec1->is_hole;
				outRec2->first_left = outRec1;

				if (_using_polytree)
					fixup_first_lefts2(outRec2, outRec1);

				if ((outRec2->is_hole ^ _reverse_output) == (area(outRec2->pts) > 0))
					reverse_poly_pt_links(outRec2->pts);

			}
			else if (poly2_contains_poly1(outRec1->pts, outRec2->pts))
			{
				//outRec2 contains outRec1 ...
				outRec2->is_hole = outRec1->is_hole;
				outRec1->is_hole = !outRec2->is_hole;
				outRec2->first_left = outRec1->first_left;
				outRec1->first_left = outRec2;

				if (_using_polytree) 
					fixup_first_lefts2(outRec1, outRec2);

				if ((outRec1->is_hole ^ _reverse_output) == (area(outRec1->pts) > 0))
					reverse_poly_pt_links(outRec1->pts);
			}
			else
			{
				//the 2 polygons are completely separate ...
				outRec2->is_hole = outRec1->is_hole;
				outRec2->first_left = outRec1->first_left;

				//fixup first_left pointers that may need reassigning to OutRec2
				if (_using_polytree) 
					fixup_first_lefts1(outRec1, outRec2);
			}

		}
		else
		{
			//joined 2 polygons together ...

			outRec2->pts = 0;
			outRec2->bot_pt = 0;
			outRec2->idx = outRec1->idx;

			outRec1->is_hole = hole_state_rec->is_hole;
			if (hole_state_rec == outRec2)
				outRec1->first_left = outRec2->first_left;
			outRec2->first_left = outRec1;

			if (_using_polytree)
				fixup_first_lefts3(outRec2, outRec1);
		}
	}
}

template<typename K>
int Clipper<K>::point_count(OutPt *pts)
{
    if (!pts) 
		return 0;
    int result = 0;
    OutPt *p = pts;
    do
    {
        result++;
        p = p->next;
    } while (p != pts);

    return result;
}

template<typename K>
void Clipper<K>::update_outpt_idxs(OutRec& outrec)
{
	OutPt* op = outrec.pts;
	do
	{
		op->idx = outrec.idx;
		op = op->prev;
	} while (op != outrec.pts);
}

template<typename K>
inline void Clipper<K>::do_simple_polygons()
{
	PolyOutList::size_type i = 0;
	while (i < _poly_outs.size())
	{
		OutRec* outrec = _poly_outs[i++];
		OutPt* op = outrec->pts;
		if (!op || outrec->is_open) 
			continue;
		do //for each pt in Polygon until duplicate found do ...
		{
			OutPt* op2 = op->next;
			while (op2 != outrec->pts)
			{
				if ((op->pt == op2->pt) && op2->next != op && op2->prev != op)
				{
					//split the polygon into two ...
					OutPt* op3 = op->prev;
					OutPt* op4 = op2->prev;
					op->prev = op4;
					op4->next = op;
					op2->prev = op3;
					op3->next = op2;

					outrec->pts = op;
					OutRec* outrec2 = create_out_rec();
					outrec2->pts = op2;
					update_outpt_idxs(*outrec2);
					if (poly2_contains_poly1(outrec2->pts, outrec->pts))
					{
						//OutRec2 is contained by OutRec1 ...
						outrec2->is_hole = !outrec->is_hole;
						outrec2->first_left = outrec;
						if (_using_polytree) 
							fixup_first_lefts2(outrec2, outrec);
					}
					else
						if (poly2_contains_poly1(outrec->pts, outrec2->pts))
						{
							//OutRec1 is contained by OutRec2 ...
							outrec2->is_hole = outrec->is_hole;
							outrec->is_hole = !outrec2->is_hole;
							outrec2->first_left = outrec->first_left;
							outrec->first_left = outrec2;
							if (_using_polytree) 
								fixup_first_lefts2(outrec, outrec2);
						}
						else
						{
							//the 2 polygons are separate ...
							outrec2->is_hole = outrec->is_hole;
							outrec2->first_left = outrec->first_left;
							if (_using_polytree) 
								fixup_first_lefts1(outrec, outrec2);
						}
					op2 = op; //ie get ready for the next iteration
				}
				op2 = op2->next;
			}
			op = op->next;
		} while (op != outrec->pts);
	}
}

template<typename K>
bool Clipper<K>::first_is_bottom_pt(const OutPt* btmPt1, const OutPt* btmPt2)
{
	OutPt* p = btmPt1->prev;
	while ((p->pt == btmPt1->pt) && (p != btmPt1)) 
		p = p->prev;
	double dx1p = std::fabs(get_delta(btmPt1->pt, p->pt));
	p = btmPt1->next;
	while ((p->pt == btmPt1->pt) && (p != btmPt1)) 
		p = p->next;
	double dx1n = std::fabs(get_delta(btmPt1->pt, p->pt));

	p = btmPt2->prev;
	while ((p->pt == btmPt2->pt) && (p != btmPt2))
		p = p->prev;
	double dx2p = std::fabs(get_delta(btmPt2->pt, p->pt));
	p = btmPt2->next;
	while ((p->pt == btmPt2->pt) && (p != btmPt2))
		p = p->next;
	double dx2n = std::fabs(get_delta(btmPt2->pt, p->pt));

	if (std::max<double>(dx1p, dx1n) == std::max<double>(dx2p, dx2n)
		&& std::min<double>(dx1p, dx1n) == std::min<double>(dx2p, dx2n))
		return area(btmPt1) > 0; //if otherwise identical use orientation
	else
		return (dx1p >= dx2p && dx1p >= dx2n) || (dx1n >= dx2p && dx1n >= dx2n);
}

template<typename K>
typename Clipper<K>:: OutRec* 
Clipper<K>::parse_first_left(OutRec* first_left)
{
	while (first_left && !first_left->pts)
		first_left = first_left->first_left;
	return first_left;
}

template<typename K>
bool Clipper<K>::poly2_contains_poly1(OutPt* out_pt1, OutPt* out_pt2)
{
	class PtIter : public std::iterator<std::forward_iterator_tag, Point2>
	{
	public:
		explicit PtIter(OutPt* op = nullptr)
			: _pt(op)
			, _ori_pt(op)
		{
		}

		PtIter& operator++()
		{
			if (_pt)
				_pt = _pt->next;
			if (_pt == _ori_pt)
				_pt = nullptr;
			return *this;
		}

		PtIter& operator--()
		{
			_pt = _pt->prev;
			return *this;
		}

		bool operator==(const PtIter& other) const
		{
			return _pt == other._pt;
		}

		bool operator!=(const PtIter& other) const
		{
			return !(*this == other);
		}

		reference operator*() const
		{
			return _pt->pt;
		}

	private:
		OutPt* _pt;
		OutPt* _ori_pt;
	};

	PtIter iterb(out_pt2);
	auto orient = CGAL::orientation_2(iterb, PtIter());
	OutPt* op = out_pt1;
	do
	{
		auto ret = CGAL::oriented_side_2(iterb, PtIter(), op->pt);
		if (ret * orient < 0)
			return false;
		op = op->next;
	} while (op != out_pt1);

	return true;
}

template<typename K>
void Clipper<K>::fixup_first_lefts1(OutRec* oldOutRec, OutRec* newOutRec)
{
	//tests if newOutRec contains the polygon before reassigning first_left
	for (PolyOutList::size_type i = 0; i < _poly_outs.size(); ++i)
	{
		OutRec* outRec = _poly_outs[i];
		OutRec* first_left = parse_first_left(outRec->first_left);
		if (outRec->pts && first_left == oldOutRec)
		{
			if (poly2_contains_poly1(outRec->pts, newOutRec->pts))
				outRec->first_left = newOutRec;
		}
	}
}

template<typename K>
void Clipper<K>::fixup_first_lefts2(OutRec* inner_outrec, OutRec* outer_outrec)
{
	//A polygon has split into two such that one is now the inner of the other.
	//It's possible that these polygons now wrap around other polygons, so check
	//every polygon that's also contained by outer_outrec's first_left container
	//(including 0) to see if they've become inner to the new inner polygon ...
	OutRec* orfl = outer_outrec->first_left;
	for (PolyOutList::size_type i = 0; i < _poly_outs.size(); ++i)
	{
		OutRec* outRec = _poly_outs[i];

		if (!outRec->pts || outRec == outer_outrec || outRec == inner_outrec)
			continue;
		OutRec* first_left = parse_first_left(outRec->first_left);
		if (first_left != orfl && first_left != inner_outrec && first_left != outer_outrec)
			continue;
		if (poly2_contains_poly1(outRec->pts, inner_outrec->pts))
			outRec->first_left = inner_outrec;
		else if (poly2_contains_poly1(outRec->pts, outer_outrec->pts))
			outRec->first_left = outer_outrec;
		else if (outRec->first_left == inner_outrec || outRec->first_left == outer_outrec)
			outRec->first_left = orfl;
	}
}

template<typename K>
void Clipper<K>::fixup_first_lefts3(OutRec* oldOutRec, OutRec* newOutRec)
{
	//reassigns first_left WITHOUT testing if newOutRec contains the polygon
	for (PolyOutList::size_type i = 0; i < _poly_outs.size(); ++i)
	{
		OutRec* outRec = _poly_outs[i];
		OutRec* first_left = parse_first_left(outRec->first_left);
		if (outRec->pts && first_left == oldOutRec)
			outRec->first_left = newOutRec;
	}
}

template<typename K>
typename Clipper<K>::Point2 
Clipper<K>::intersect_point(TEdge &edge1, TEdge &edge2)
{
    auto result = intersection(edge1, edge2);
    Point2 ret= boost::get<Point2>(result.get());
    return ret;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

template<typename K>
double ClipperOffset<K>::distance_sqrt(const Point2& pt1, const Point2& pt2)
{
	double dx = ((double)pt1.x() - pt2.x());
	double dy = ((double)pt1.y() - pt2.y());
	return (dx * dx + dy * dy);
}

template<typename K>
double ClipperOffset<K>::distance_from_linesqrd( const Point2& pt, const Point2& ln1, const Point2& ln2)
{
	//The equation of a line in general form (Ax + By + C = 0)
	//given 2 points (x?y? & (x?y? is ...
	//(y?- y?x + (x?- x?y + (y?- y?x?- (x?- x?y?= 0
	//A = (y?- y?; B = (x?- x?; C = (y?- y?x?- (x?- x?y?
	//perpendicular distance of point (x?y? = (Ax?+ By?+ C)/Sqrt(A?+ B?
	//see http://en.wikipedia.org/wiki/Perpendicular_distance
	double A = double(ln1.y() - ln2.y());
	double B = double(ln2.x() - ln1.x());
	double C = A * ln1.x() + B * ln1.y();
	C = A * pt.x() + B * pt.y() - C;
	return (C * C) / (A * A + B * B);
}

template<typename K>
typename Clipper<K>:: OutPt* 
ClipperOffset<K>::exclude_op(typename Clipper<K>::OutPt* op)
{
	Clipper<K>::OutPt* result = op->prev;
	result->next = op->next;
	op->next->prev = result;
	result->idx = 0;
	return result;
}

template<typename K>
typename ClipperOffset<K>:: Vector2 
ClipperOffset<K>::get_unit_normal(const Point2& pt1, const Point2& pt2)
{
	if (pt2.x() == pt1.x() && pt2.y() == pt1.y())
		return Vector2(0, 0);

	double dx = (double)(pt2.x() - pt1.x());
	double dy = (double)(pt2.y() - pt1.y());
	double f = 1 * 1.0 / std::sqrt(dx * dx + dy * dy);
	dx *= f;
	dy *= f;
	return Vector2(dy, -dx);
}

//---------------------------------------------------------------------------

template<typename K>
ClipperOffset<K>::ClipperOffset(double miterLimit, double roundPrecision)
{
	_miter_limit = miterLimit;
	_arc_tolerance = roundPrecision;
	_lowest_x = -1;
}

template<typename K>
ClipperOffset<K>::~ClipperOffset()
{
	clear();
}

template<typename K>
void ClipperOffset<K>::add_path(const Path& path, JointType joinType, EndType endType) 
{
	int highI = (int)path.size() - 1;
	if (highI < 0)
		return;
	PolyNode* new_node = new PolyNode();
	new_node->_join_type = joinType;
	new_node->_end_type = endType;

	//strip duplicate points from path and also get index to the lowest point ...
	if (endType == EndType::ET_ClosedLine || endType == EndType::ET_ClosedPolygon)
	{
		while (highI > 0 && path[0] == path[highI])
			highI--;
	}
	new_node->_contour.reserve(highI + 1);
	new_node->_contour.push_back(path[0]);
	int j = 0, k = 0;
	for (int i = 1; i <= highI; i++)
	{
		if (new_node->_contour[j] != path[i])
		{
			j++;
			new_node->_contour.push_back(path[i]);
			if (path[i].y() < new_node->_contour[k].y() ||
				(path[i].y() == new_node->_contour[k].y() &&
					path[i].x() < new_node->_contour[k].x()))
				k = j;
		}
	}
	if (endType == EndType::ET_ClosedPolygon && j < 2)
	{
		delete new_node;
		return;
	}
	_poly_nodes.add_child(new_node);

	//if this path's lowest pt is lower than all the others then update _lowest_x
	if (endType != EndType::ET_ClosedPolygon)
		return;
	if (_lowest_x < 0)
	{
		_lowest_x = _poly_nodes.child_count() - 1;
		_lowest_y = k;
	}
	else
	{
		Point2 ip = _poly_nodes._childs[_lowest_x]->_contour[_lowest_y];
		if (new_node->_contour[k].y() < ip.y() 
			|| (new_node->_contour[k].y() == ip.y() && new_node->_contour[k].x() < ip.x()))
		{
			_lowest_x = _poly_nodes.child_count() - 1;
			_lowest_y = k;
		}
	}

}

template<typename K>
void ClipperOffset<K>::add_paths(const Paths& paths, JointType joinType, EndType endType)
{
	for (Paths::size_type i = 0; i < paths.size(); ++i)
		add_path(paths[i], joinType, endType);
}

template<typename K>
void ClipperOffset<K>::execute(Paths& solution, double delta)
{
	solution.clear();
	fix_orientations();
	do_offset(delta);

	//now clean up 'corners' ...
	Clipper<K> clpr;
	clpr.add_paths(_dest_polys, PolyType::PT_Subject, true);
	
	if (delta > 0)
	{
		clpr.execute(ClipType::CT_Union, solution, PolyFillType::PFT_Positive, PolyFillType::PFT_Positive);
	}
	else
	{
		Rectangle2<K> r = clpr.get_bounds();
		Path outer(4);
		outer[0] = Point2(r.left - 10, r.bottom + 10);
		outer[1] = Point2(r.right + 10, r.bottom + 10);
		outer[2] = Point2(r.right + 10, r.top - 10);
		outer[3] = Point2(r.left - 10, r.top - 10);

		clpr.add_path(outer, PolyType::PT_Subject, true);
		clpr.set_reverse_solution(true);
		clpr.execute(ClipType::CT_Union, solution, PolyFillType::PFT_Negative, PolyFillType::PFT_Negative);
		if (solution.size() > 0)
			solution.erase(solution.begin());
	}
}

template<typename K>
void ClipperOffset<K>::execute(PolyTree& solution, double delta)
{
	solution.clear();
	fix_orientations();
	do_offset(delta);

	//now clean up 'corners' ...
	Clipper<K> clpr;
	clpr.add_paths(_dest_polys, PolyType::PT_Subject, true);
	if (delta > 0)
	{
		clpr.execute(ClipType::CT_Union, solution, PolyFillType::PFT_Positive, PolyFillType::PFT_Positive);
	}
	else
	{
		Rectangle2<K> r = clpr.get_bounds();
		Path outer(4);
		outer[0] = Point2(r.left - 10, r.bottom + 10);
		outer[1] = Point2(r.right + 10, r.bottom + 10);
		outer[2] = Point2(r.right + 10, r.top - 10);
		outer[3] = Point2(r.left - 10, r.top - 10);

		clpr.add_path(outer, PolyType::PT_Subject, true);
		clpr.set_reverse_solution(true);
		clpr.execute(ClipType::CT_Union, solution, PolyFillType::PFT_Negative, PolyFillType::PFT_Negative);
		//remove the outer PolyNode rectangle ...
		if (solution.child_count() == 1 && solution._childs[0]->child_count() > 0)
		{
			PolyNode* outerNode = solution._childs[0];
			solution._childs.reserve(outerNode->child_count());
			solution._childs[0] = outerNode->_childs[0];
			solution._childs[0]->_parent = outerNode->_parent;
			for (int i = 1; i < outerNode->child_count(); ++i)
				solution.add_child(outerNode->_childs[i]);
		}
		else
			solution.clear();
	}
}

template<typename K>
void ClipperOffset<K>::clear()
{
	for (int i = 0; i < _poly_nodes.child_count(); ++i)
		delete _poly_nodes._childs[i];
	_poly_nodes._childs.clear();
	_lowest_x = -1;
}

template<typename K>
void ClipperOffset<K>::fix_orientations()
{
	//fixup orientations of all closed paths if the orientation of the
	//closed path with the lowermost vertex is wrong ...
	if (orientation<K>(_poly_nodes._childs[_lowest_x]->_contour))
	{
		for (int i = 0; i < _poly_nodes.child_count(); ++i)
		{
			PolyNode& node = *_poly_nodes._childs[i];
			if (node._end_type == EndType::ET_ClosedPolygon ||
				(node._end_type == EndType::ET_ClosedLine && orientation<K>(node._contour)))
				reverse_path<K>(node._contour);
		}
	}
	else
	{
		for (int i = 0; i < _poly_nodes.child_count(); ++i)
		{
			PolyNode& node = *_poly_nodes._childs[i];
			if (node._end_type == EndType::ET_ClosedLine && !orientation<K>(node._contour))
				reverse_path<K>(node._contour);
		}
	}
}

template<typename K>
void ClipperOffset<K>::do_offset(double delta)
{
	_dest_polys.clear();
	_delta = delta;

	//if Zero offset, just copy any CLOSED polygons to m_p and return ...
	if (CGAL::compare(delta, 0) == EQUAL)
	{
		_dest_polys.reserve(_poly_nodes.child_count());
		for (int i = 0; i < _poly_nodes.child_count(); i++)
		{
			PolyNode *node = _poly_nodes._childs[i];
			if (node->_end_type == EndType::ET_ClosedPolygon)
				_dest_polys.push_back(node->_contour);
		}
		return;
	}

	//see offset_triginometry3.svg in the documentation folder ...
	if (_miter_limit > 2)
		_miterlimit = 2 / (_miter_limit * _miter_limit);
	else 
		_miterlimit = 0.5;

	double y;
	if (_arc_tolerance <= 0.0)
		y = arc_tolerance;
	else if (_arc_tolerance > std::fabs(delta) * arc_tolerance)
		y = std::fabs(delta) * arc_tolerance;
	else
		y = _arc_tolerance;

	//see offset_triginometry2.svg in the documentation folder ...
	double steps = pi / std::acos(1 - y / std::fabs(delta));
	if (steps > std::fabs(delta) * pi)
		steps = std::fabs(delta) * pi;  //ie excessive precision check
	_sin = std::sin(pi2 / steps);
	_cos = std::cos(pi2 / steps);
	_steps_per_rad = steps / pi2;
	if (delta < 0.0)
		_sin = -_sin;

	_dest_polys.reserve(_poly_nodes.child_count() * 2);
	for (int i = 0; i < _poly_nodes.child_count(); i++)
	{
		PolyNode& node = *_poly_nodes._childs[i];
		_src_poly = node._contour;

		int len = (int)_src_poly.size();
		if (len == 0 || (delta <= 0 && (len < 3 || node._end_type != EndType::ET_ClosedPolygon)))
			continue;

		_dest_poly.clear();
		if (len == 1)
		{
			if (node._join_type == JointType::JT_Round)
			{
				double x = 1.0, y = 0.0;
				for (int j = 1; j <= steps; j++)
				{
					_dest_poly.push_back(Point2( _src_poly[0].x() + x * delta, _src_poly[0].y() + y * delta));
					double x2 = x;
					x = x * _cos - _sin * y;
					y = x2 * _sin + y * _cos;
				}
			}
			else
			{
				double x = -1.0, y = -1.0;
				for (int j = 0; j < 4; ++j)
				{
					_dest_poly.push_back(Point2( _src_poly[0].x() + x * delta, _src_poly[0].y() + y * delta));
					if (x < 0)
						x = 1;
					else if (y < 0) 
						y = 1;
					else 
						x = -1;
				}
			}
			_dest_polys.push_back(_dest_poly);
			continue;
		}

		//build _normals ...
		_normals.clear();
		_normals.reserve(len);
		for (int j = 0; j < len - 1; ++j)
			_normals.push_back(get_unit_normal(_src_poly[j], _src_poly[j + 1]));
		if (node._end_type == EndType::ET_ClosedLine || node._end_type == EndType::ET_ClosedPolygon)
			_normals.push_back(get_unit_normal(_src_poly[len - 1], _src_poly[0]));
		else
			_normals.push_back(_normals[len - 2]);

		if (node._end_type == EndType::ET_ClosedPolygon)
		{
			int k = len - 1;
			for (int j = 0; j < len; ++j)
				offset_point(j, k, node._join_type);
			_dest_polys.push_back(_dest_poly);
		}
		else if (node._end_type == EndType::ET_ClosedLine)
		{
			int k = len - 1;
			for (int j = 0; j < len; ++j)
				offset_point(j, k, node._join_type);
			_dest_polys.push_back(_dest_poly);
			_dest_poly.clear();
			//re-build _normals ...
			Vector2 n = _normals[len - 1];
			for (int j = len - 1; j > 0; j--)
				_normals[j] = Vector2(-_normals[j - 1].x(), -_normals[j - 1].y());
			_normals[0] = Vector2(-n.x(), -n.y());
			k = 0;
			for (int j = len - 1; j >= 0; j--)
				offset_point(j, k, node._join_type);
			_dest_polys.push_back(_dest_poly);
		}
		else
		{
			int k = 0;
			for (int j = 1; j < len - 1; ++j)
				offset_point(j, k, node._join_type);

			Point2 pt1;
			if (node._end_type == EndType::ET_OpenButt)
			{
				int j = len - 1;
				pt1 = Point2((_src_poly[j].x() + _normals[j].x() * delta), (_src_poly[j].y() + _normals[j].y() * delta));
				_dest_poly.push_back(pt1);
				pt1 = Point2((_src_poly[j].x() - _normals[j].x() * delta), (_src_poly[j].y() - _normals[j].y() * delta));
				_dest_poly.push_back(pt1);
			}
			else
			{
				int j = len - 1;
				k = len - 2;
				_sinA = 0;
				_normals[j] = Vector2(-_normals[j].x(), -_normals[j].y());
				if (node._end_type == EndType::ET_OpenSquare)
					do_square(j, k);
				else
					do_round(j, k);
			}

			//re-build _normals ...
			for (int j = len - 1; j > 0; j--)
				_normals[j] = Vector2(-_normals[j - 1].x(), -_normals[j - 1].y());
			_normals[0] = Vector2(-_normals[1].x(), -_normals[1].y());

			k = len - 1;
			for (int j = k - 1; j > 0; --j) 
				offset_point(j, k, node._join_type);

			if (node._end_type == EndType::ET_OpenButt)
			{
				pt1 = Point2((_src_poly[0].x() - _normals[0].x() * delta), (_src_poly[0].y() - _normals[0].y() * delta));
				_dest_poly.push_back(pt1);
				pt1 = Point2((_src_poly[0].x() + _normals[0].x() * delta), (_src_poly[0].y() + _normals[0].y() * delta));
				_dest_poly.push_back(pt1);
			}
			else
			{
				k = 1;
				_sinA = 0;
				if (node._end_type == EndType::ET_OpenSquare)
					do_square(0, 1);
				else
					do_round(0, 1);
			}
			_dest_polys.push_back(_dest_poly);
		}
	}
}

template<typename K>
void ClipperOffset<K>::offset_point(int j, int& k, JointType jointype)
{
	_sinA = _normals[k].x() * _normals[j].y() - _normals[k].y() * _normals[j].x();

	if (_sinA * _delta < 0)
	{
		_dest_poly.push_back(Point2((_src_poly[j].x() + _normals[k].x() * _delta),
			(_src_poly[j].y() + _normals[k].y() * _delta)));
		_dest_poly.push_back(_src_poly[j]);
		_dest_poly.push_back(Point2((_src_poly[j].x() + _normals[j].x() * _delta),
			(_src_poly[j].y() + _normals[j].y() * _delta)));
	}
	else
	{
		switch (jointype)
		{
		case JointType::JT_Miter:
		{
			double r = 1 + _normals[j].x() * _normals[k].x() + _normals[j].y() * _normals[k].y();
			if (r >= _miterlimit)
				do_miter(j, k, r);
			else
				do_square(j, k);
			break;
		}
		case JointType::JT_Square:
			do_square(j, k);
			break;
		case JointType::JT_Round:
			do_round(j, k);
			break;
		}
	}
	k = j;
}

template<typename K>
void ClipperOffset<K>::do_square(int j, int k)
{
	double dx = std::tan(std::atan2(_sinA,
		_normals[k].x() * _normals[j].x() + _normals[k].y() * _normals[j].y()) / 4);
	_dest_poly.push_back(Point2(
		_src_poly[j].x() + _delta * (_normals[k].x() - _normals[k].y() * dx),
		_src_poly[j].y() + _delta * (_normals[k].y() + _normals[k].x() * dx)));
	_dest_poly.push_back(Point2(
		_src_poly[j].x() + _delta * (_normals[j].x() + _normals[j].y() * dx),
		_src_poly[j].y() + _delta * (_normals[j].y() - _normals[j].x() * dx)));
}

template<typename K>
void ClipperOffset<K>::do_miter(int j, int k, double r)
{
	double q = _delta / r;
	_dest_poly.push_back(Point2(_src_poly[j].x() + (_normals[k].x() + _normals[j].x()) * q,
		_src_poly[j].y() + (_normals[k].y() + _normals[j].y()) * q));
}

template<typename K>
void ClipperOffset<K>::do_round(int j, int k)
{
	double a = std::atan2(_sinA, _normals[k].x() * _normals[j].x() + _normals[k].y() * _normals[j].y());
	int steps = std::max<double>(_steps_per_rad * std::fabs(a), 1);

	double x = _normals[k].x(), y = _normals[k].y(), X2;
	for (int i = 0; i < steps; ++i)
	{
		_dest_poly.push_back(Point2(_src_poly[j].x() + x * _delta, _src_poly[j].y() + y * _delta));
		X2 = x;
		x = x * _cos - _sin * y;
		y = X2 * _sin + y * _cos;
	}
	_dest_poly.push_back(Point2(_src_poly[j].x() + _normals[j].x() * _delta, _src_poly[j].y() + _normals[j].y() * _delta));
}



template<typename K>
void test()
{
	Clipper<K> cl;
	typedef Clipper<K>::Point2 Point2;

	int ex = 5;

	if (ex == 1) {
		Clipper<K>::Path poly, poly1;
		poly.push_back(Point2(0, 0));
		poly.push_back(Point2(10, 10));
		poly.push_back(Point2(0, 20));
		poly.push_back(Point2(0, 0));

		poly1.push_back(Point2(15, 0));
		poly1.push_back(Point2(15, 20));
		poly1.push_back(Point2(5, 10));
		poly1.push_back(Point2(15, 0));

		cl.add_path(poly, PolyType::PT_Subject);
		cl.add_path(poly1, PolyType::PT_Clip);

		Clipper<K>::Paths out;
		cl.execute(ClipType::CT_Intersect, out);
	}
	else if (ex == 2) {
		Clipper<K>::Path poly, poly1;
		poly.push_back(Point2(0, 0));
		poly.push_back(Point2(-2, 2));
		poly.push_back(Point2(5, 10));
		poly.push_back(Point2(0, 0));

		poly1.push_back(Point2(1, 1));
		poly1.push_back(Point2(0, 8));
		poly1.push_back(Point2(3, 1.5));
		poly1.push_back(Point2(1, 1));

		cl.add_path(poly, PolyType::PT_Subject);
		cl.add_path(poly1, PolyType::PT_Clip);

		Clipper<K>::Paths out;
		cl.execute(ClipType::CT_Intersect, out);
	}
	else if (ex == 3) {
		Clipper<K>::Path poly, poly1;
		poly.push_back(Point2(0, 0));
		poly.push_back(Point2(-5, 5));
		poly.push_back(Point2(0, 2));
		poly.push_back(Point2(5, 5));
		poly.push_back(Point2(0, 0));

		poly1.push_back(Point2(-5, 0));
		poly1.push_back(Point2(0, 5));
		poly1.push_back(Point2(5, -5));
		poly1.push_back(Point2(-5, 0));

		cl.add_path(poly, PolyType::PT_Subject);
		cl.add_path(poly1, PolyType::PT_Clip);

		Clipper<K>::Paths out;
		cl.execute(ClipType::CT_Intersect, out);
	}
	else if (ex == 4)
	{
		Clipper<K>::Path poly, poly1;
		poly.push_back(Point2(0, 0));
		poly.push_back(Point2(10, 10));
		poly.push_back(Point2(-10, 10));
		poly.push_back(Point2(0, 0));

		poly1.push_back(Point2(0, 5));
		poly1.push_back(Point2(-10, -10));
		poly1.push_back(Point2(0, -10));
		poly1.push_back(Point2(10, -10));
		poly1.push_back(Point2(0, 5));

		cl.add_path(poly, PolyType::PT_Subject);
		cl.add_path(poly1, PolyType::PT_Clip);

		Clipper<K>::Paths out;
		cl.execute(ClipType::CT_Intersect, out);
	}
	else if (ex == 5)
	{
		Clipper<K>::Path poly;
		poly.push_back(Point2(0, 0));
		poly.push_back(Point2(0, 10));
		poly.push_back(Point2(10, 5));
		cl.add_path(poly, PolyType::PT_Subject, true);

		//cl.setStrictlySimple(true);
		PolyTree<K> out;
		cl.execute(ClipType::CT_Union, out);
		printf("");
	}
	else if (ex == 6)
	{
		Clipper<K>::Path poly;
		poly.push_back(Point2(0, 10));
		poly.push_back(Point2(0, -10));
		poly.push_back(Point2(10, 0));
		poly.push_back(Point2(-10, 5));
		//poly.push_back(Point2(0, 10));
		cl.add_path(poly, PolyType::PT_Subject, false);

		//cl.setStrictlySimple(true);
		PolyTree<K> out;
		cl.execute(ClipType::CT_Union, out, PolyFillType::PFT_Positive);
		printf("");
	}
	else if (ex == 7)
	{
		Clipper<K>::Path poly;
		poly.push_back(Point2(0, 10));
		poly.push_back(Point2(0, -10));
		poly.push_back(Point2(10, 0));
		poly.push_back(Point2(-10, 5));
		cl.add_path(poly, PolyType::PT_Subject, false);

		cl.setStrictlySimple(true);
		PolyTree<K> out;
		cl.execute(ClipType::CT_Union, out, PolyFillType::PFT_Negative);
		printf("");
	}
	else if (ex == 8)
	{
		Clipper<K>::Path poly;
		poly.push_back(Point2(0, 10));
		poly.push_back(Point2(0, -10));
		poly.push_back(Point2(10, 0));
		cl.add_path(poly, PolyType::PT_Subject, false);

		//cl.setStrictlySimple(true);
		PolyTree<K> out;
		cl.execute(ClipType::CT_Union, out, PolyFillType::PFT_Positive);
		printf("");
	}
	else if (ex == 9)
	{
		Clipper<K>::Path poly;
		poly.push_back(Point2(0, -10));
		poly.push_back(Point2(10, 0));
		poly.push_back(Point2(0, 10));
		cl.add_path(poly, PolyType::PT_Subject, false);

		//cl.setStrictlySimple(true);
		PolyTree<K> out;
		cl.execute(ClipType::CT_Union, out, PolyFillType::PFT_Positive);
		printf("");
	}
	else if (ex == 10)
	{
		Clipper<K>::Path poly;
		poly.push_back(Point2(-10, 0));
		poly.push_back(Point2(10, 0));
		poly.push_back(Point2(0, 10));
		cl.add_path(poly, PolyType::PT_Subject, false);

		//cl.setStrictlySimple(true);
		PolyTree<K> out;
		cl.execute(ClipType::CT_Union, out, PolyFillType::PFT_Positive);
		printf("");
	}
	else if (ex == 11)
	{
		Clipper<K>::Path poly;
		poly.push_back(Point2(0, 0));
		poly.push_back(Point2(10, 10));
		poly.push_back(Point2(20, 0));
		poly.push_back(Point2(5, 4));
		poly.push_back(Point2(15, 4));
		cl.add_path(poly, PolyType::PT_Subject, true);

		//cl.setStrictlySimple(true);
		PolyTree<K> out;
		cl.setStrictlySimple(true);
		cl.execute(ClipType::CT_Union, out);
		printf("");

	}
	else if (ex == 12)
	{
		Clipper<K>::Path poly;
		poly.push_back(Point2(0, 0));
		poly.push_back(Point2(10, 10));
		poly.push_back(Point2(0, 10));
		poly.push_back(Point2(10, 0));
		cl.add_path(poly, PolyType::PT_Subject, true);

		PolyTree<K> out;
		cl.execute(ClipType::CT_Union, out);
		printf("");
	}
	else if (ex == 101)
	{
		ClipperOffset<K>::Path poly;
		ClipperOffset<K>::Paths outPath;
		poly.push_back(Point2(0, 0));
		poly.push_back(Point2(10, 0));
		poly.push_back(Point2(0, 10));
		poly.push_back(Point2(10, 10));
		poly.push_back(Point2(0, 0));
		ClipperOffset<K> co(4, 2.5);
		co.add_path(poly, JointType::JT_Miter, EndType::ET_ClosedPolygon);

		//co.setStrictlySimple(true);
		co.execute(outPath, 2);
		printf("");
	}
	else if (ex == 102) {
		ClipperOffset<K>::Path poly;
		ClipperOffset<K>::Paths outPath;
		poly.push_back(Point2(0, 0));
		poly.push_back(Point2(10, 0));
		poly.push_back(Point2(10, 10));
		poly.push_back(Point2(0, 0));
		ClipperOffset<K> co(4, 2.5);
		co.add_path(poly, JointType::JT_Miter, EndType::ET_ClosedPolygon);

		//co.setStrictlySimple(true);
		co.execute(outPath, 2);
		printf("");
	}
}

}
}

#endif
