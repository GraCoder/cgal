#ifndef CGAL_POLYGON_2_CLIPPER_H
#define CGAL_POLYGON_2_CLIPPER_H

#include <CGAL/config.h>
#include <vector>
#include <list>
#include <iterator>

#include <CGAL/Rectangular_p_center_traits_2.h>
#include <CGAL/Boolean_set_operations_2.h>

namespace CGAL {
namespace Polygon2Clipper {

class ClipperException : public std::exception
{
public:
	ClipperException(const char* description) : _descr(description) {}
	virtual ~ClipperException() {}
	virtual const char* what() const throw() { return _descr.c_str(); }
private:
	std::string		_descr;
};

template<typename K>
struct Rectangle2 {
	typedef typename K::FT RT;

	RT left;	RT top;
	RT right;	RT bottom;
};

enum class PolyType { ptSubject, ptClip };
enum class ClipType { ctIntersection, ctUnion, ctDifference, ctXor };
//http://glprogramming.com/red/chapter11.html
enum class PolyFillType { pftEvenOdd, pftNonZero, pftPositive, pftNegative };
enum class InitOptions { ioReverseSolution = 1, ioStrictlySimple = 2, ioPreserveCollinear = 4 };
enum class JointType { jtSquare, jtRound, jtMiter };
enum class EndType	{ etClosedPolygon, etClosedLine, etOpenButt, etOpenSquare, etOpenRound };

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
		, _isOpen(false)
	{ };
	virtual ~PolyNode() {};

	bool isHole() const
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
	bool isOpen() const { return _isOpen; }
	int	 childCount() const { return _childs.size(); }
	PolyNode* getNext() const
	{
		if (_childs.empty())
			return getNextSiblingUp();
		else
			return _childs[0];
	}
protected:
	void addChild(PolyNode *child)
	{
		_childs.push_back(child);
		child->_parent = this;
		child->_index = _childs.size() - 1;
	}
	PolyNode* getNextSiblingUp() const
	{
		if (!_parent)
			return nullptr;
		else if (_index == _parent->_childs.size() - 1)
			return _parent->getNextSiblingUp();
		else
			return _parent->_childs[_index + 1];
	}

	bool					_isOpen;
	unsigned				_index; //node index in _parent._childs
	JointType				_joinType;
	EndType					_endType;
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
		int ret = _allNodes.size();
		if (ret > 0 && _childs[0] != _allNodes[0])
			ret--;
		return ret;
	}
	void clear()
	{
		_childs.clear();
		for (int i = 0; i < _allNodes.size(); i++)
			delete _allNodes[i];
		_allNodes.clear();
	}
	PolyNode<K>* getFirst() const
	{
		if (_childs.empty())
			return nullptr;

		return _childs.front();
	}

private:
	std::vector<PolyNode<K>*> _allNodes;
};


template<typename K = Exact_predicates_inexact_constructions_kernel>
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
	typedef std::vector<typename K::Point_2>  Path;
	typedef std::vector<Path> Paths;
public:
	Clipper();
	virtual ~Clipper();

	bool	addPath(const Path& pg, PolyType type = PolyType::ptSubject, bool closed = true);
	bool	addPaths(const Paths& pg, PolyType PolyTyp, bool closed);
	bool	execute(ClipType clipType, Paths& solution, PolyFillType fillType = PolyFillType::pftEvenOdd);
	bool	execute(ClipType clipType, Paths& solution, PolyFillType subjFillType, PolyFillType clipFillType);
	bool	execute(ClipType clipType, PolyTree<K>& polytree, PolyFillType fillType = PolyFillType::pftEvenOdd);
	bool	execute(ClipType clipType, PolyTree<K>& polytree, PolyFillType subjFillType, PolyFillType clipFillType);
	void	clear();

	Rectangle2<K>	getBounds();

	inline bool		preserveCollinear() { return _preserveCollinear; }
	inline void		setPreserveCollinear(bool value) { _preserveCollinear = value; };
	inline bool		reverseSolution() { return _reverseOutput; };
	inline void		setReverseSolution(bool value) { _reverseOutput = value; };
	inline bool		strictlySimple() { return _strictSimple; };
	inline void		setStrictlySimple(bool value) { _strictSimple = value; };
protected:	
	struct TEdge : Segment2 {
		Point2		cur;
#ifdef DEBUG
		DebugPt     dbot;
		DebugPt     dtop;
#endif
		double		dx;
		PolyType	polyType;
		EdgeSide	side;
		int			windDelta;
		int			windCnt;
		int			windCnt2;
		int			outIdx;
		TEdge		*next = nullptr;
		TEdge		*prev = nullptr;
		TEdge		*nextInLML = nullptr;
		TEdge		*nextInAEL = nullptr;
		TEdge		*prevInAEL = nullptr;
		TEdge		*nextInSEL = nullptr;
		TEdge		*prevInSEL = nullptr;

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
		TEdge		*leftBound = nullptr;
		TEdge		*rightBound = nullptr;
		FT				y;
	};

	struct OutPt;
	struct OutRec {
		int			idx;
		bool		isHole;
		bool		isOpen;
		OutPt		*pts;
		OutPt		*botPt;
		OutRec		*firstLeft;  //see comments in clipper.pas
		PolyNode<K>	*polyNode;
	};

	struct OutPt {
		int			 idx;
		OutPt		*next;
		OutPt		*prev;
		Point2		 pt;
	};

	struct Joint {
		OutPt		*outPt1;
		OutPt		*outPt2;
		Point2		 offPt;
	};
	
	typedef std::vector<std::vector<TEdge> > EdgeList;
	typedef std::vector<Joint> JoinList;
	typedef std::vector<IntersectNode> IntersectList;
	typedef std::vector<LocalMinimum>		MinimaList;
	typedef std::vector<OutRec *> PolyOutList;
	typedef std::vector<PolyNode<K> *> PolyNodes;
	typedef std::list<CT> MaximaList;
	typedef std::priority_queue<CT, std::vector<CT>, std::greater<int>> ScanbeamList;

	static inline bool		isHorizontal(TEdge &e);
	static inline void		setDelta(TEdge &e);
	static inline double	getDelta(const Point2& pt1, const Point2& pt2);
	static inline void		initEdge(TEdge *e, TEdge *enext, TEdge *eprev, const Point2 &pt);
	static inline void		initEdge2(TEdge &e, PolyType type);
	static inline bool		slopesEqual(const Point2 &pt1, const Point2 &pt2, const Point2 &pt3);
	static inline bool		slopesEqual(const Point2 &pt1, const Point2 &pt2, const Point2 &pt3, const Point2 &pt4);
	static inline bool		slopesEqual(const TEdge& e1, const TEdge& e2);
	static inline bool		pt2BetweenPt1AndPt3(const Point2 pt1, const Point2 pt2, const Point2 pt3);
	static inline void		reverseHorizontal(TEdge &e);
	static inline void		reversePolyPtLinks(OutPt *pp);
	static inline void		disposeOutPts(OutPt *&pp); 
	static inline bool		horzSegmentsOverlap(FT, FT, FT, FT);
	static inline void		swapSides(TEdge &edge1, TEdge &edge2);
	static inline void		swapPolyIndexes(TEdge &edge1, TEdge &edge2);
	static inline bool		e2InsertsBeforeE1(TEdge &e1, TEdge &e2);
	static inline bool		isMinima(TEdge *e, const FT y);
	static inline bool		isMaxima(TEdge *e, const FT y);
	static inline bool		isIntermediate(TEdge *e, const FT y);
	static inline bool		edgesAdjacent(const IntersectNode &inode);
	static inline FT		topX(TEdge &edge, const FT currentY);
	static inline TEdge*	removeEdge(TEdge *e);
	static inline TEdge*	getMaximaPair(TEdge *e);
	static inline TEdge*	getMaximaPairEx(TEdge *e);
	static inline TEdge*	getNextInAEL(TEdge *e, Orientation dir);
	static inline double	area(const OutPt *op);
	static inline double	area(const Path& path);

	//--------------------------------------------------------------------------------------------------------------

	TEdge*		findNextLocMin(TEdge *e);
	TEdge*		processBound(TEdge *e, bool nextIsForward);

	//--------------------------------------------------------------------------------------------------------------

	bool		executeInternal();

	//--------------------------------------------------------------------------------------------------------------

	void		reset();
	void		disposeLocalMinimaList();
	bool		popLocalMinima(CT y, const LocalMinimum *&locMin);
	void		insertScanbeam(const CT y);
	bool		popScanbeam(CT &y);
	void		disposeAllOutRecs();
	void		disposeOutRec(typename PolyOutList::size_type index);
	void		swapPositionsInAEL(TEdge *edge1, TEdge *edge2);
	void		deleteFromAEL(TEdge *e);
	OutRec*		createOutRec();
	void		updateEdgeIntoAEL(TEdge *&e);
	bool		localMinimaPending();

	//--------------------------------------------------------------------------------------------------------------

	void		setWindingCount(TEdge &edge);
	bool		isEvenOddFillType(const TEdge &edge) const;
	bool		isEvenOddAltFillType(const TEdge &edge) const;
	void		getHorzDirection(TEdge& horzEdge, Orientation& dir, CT& left, CT& right);
	bool		getOverlap(const FT a1, const FT a2, const FT b1, const FT b2, FT& left, FT& right);
	void		insertLocalMinimaIntoAEL(const FT botY);
	void		insertEdgeIntoAEL(TEdge *edge, TEdge *startEdge);
	void		addEdgeToSEL(TEdge *edge);
	bool		popEdgeFromSEL(TEdge *&edge);
	void		copyAELToSEL();
	void		deleteFromSEL(TEdge *e);
	void		swapPositionInSEL(TEdge *edge1, TEdge *edge2);
	bool		isContributing(const TEdge &edge) const;
	void		doMaxima(TEdge *e);
	void		processHorizontals();
	void		processHorizontal(TEdge *horzEdge);
	OutPt*		addOutPt(TEdge *e, const Point2 &pt);
	OutPt*		getLastOutPt(TEdge *e);
	OutPt*		addLocalMinPoly(TEdge *e1, TEdge *e2, const Point2 &pt);
	OutPt*		dupOutPt(OutPt* outPt, bool InsertAfter);
	OutPt*		getBottomPt(OutPt* pp);
	void		addLocalMaxPoly(TEdge *e1, TEdge *e2, const Point2 &pt);
	OutRec*		getOutRec(int idx);
	OutRec*		getLowermostRec(OutRec *outRec1, OutRec *outRec2);
	void		appendPolygon(TEdge *e1, TEdge *e2);
	void		intersectEdges(TEdge *e1, TEdge *e2, Point2 &pt);
	bool		processIntersections(const CT topY);
	void		buildIntersectList(const CT topY);
	void		processIntersectList();
	void		processEdgesAtTopOfScanbeam(const CT topY);
	void		buildResult(Paths &polys);
	void		buildResult2(PolyTree<K> &polytree);
	void		setHoleState(TEdge *e, OutRec *outrec);
	void		disposeIntersectNodes();
	bool		fixupIntersectionOrder();
	void		fixupOutPolygon(OutRec &outrec);
	void		fixupOutPolyline(OutRec &outrec);
	bool		outRec1RightOfOutRec2(OutRec *outRec1, OutRec *outRec2);
	//bool isHole(TEdge *e);
	//bool FindOwnerFromSplitRecs(OutRec &outRec, OutRec *&currOrfl);
	void		fixHoleLinkage(OutRec &outrec);
	void		addJoin(OutPt *op1, OutPt *op2, const Point2 offPt);
	void		clearJoints();
	void		clearGhostJoints();
	void		addGhostJoin(OutPt *op, const Point2 offPt);
	bool		joinHorz(OutPt* op1, OutPt* op1b, OutPt* op2, OutPt* op2b, const Point2 pt, bool discardLeft);
	bool		joinPoints(Joint *j, OutRec *outRec1, OutRec *outRec2);
	void		joinCommonEdges();
	int			pointCount(OutPt *pts);
	void		updateOutPtIdxs(OutRec& outrec);
	void		doSimplePolygons();
	bool		firstIsBottomPt(const OutPt* btmPt1, const OutPt* btmPt2);
	OutRec*		parseFirstLeft(OutRec* firstLeft);
	bool		poly2ContainsPoly1(OutPt* outPt1, OutPt* outPt2);
	void		fixupFirstLefts1(OutRec *oldOutRec, OutRec *newOutRec);
	void		fixupFirstLefts2(OutRec *innerOutRec, OutRec *outerOutRec);
	void		fixupFirstLefts3(OutRec *oldOutRec, OutRec *newOutRec);
	Point2		intersectPoint(TEdge &edge1, TEdge &edge2);

private:

	const static	int _unassigned = -1;  //edge not currently 'owning' a solution
	const static	int _skip = -2;        //edge that would otherwise close a path

	bool			_openPath;
	bool			_useFullRange;
	bool			_preserveCollinear;
	bool			_excuteLocked;
	bool			_reverseOutput;
	bool			_usingPolyTree;
	bool			_strictSimple;
	TEdge			*_activeEdges;
	TEdge			*_sortedEdges;
#ifdef DEBUG
	int				_activeCount;
	int				_sortedCount;
#endif

	typename MinimaList::iterator	_currentLM;

	EdgeList		_edges;
	MinimaList		_minimaList;
	PolyOutList		_polyOuts;
	ScanbeamList	_scanbeam;
	JoinList		_joints;
	JoinList		_ghostJoints;
	IntersectList	_intersectList;
	ClipType		_clipType;
	MaximaList		_maxima;
	PolyFillType	_clipFillType;
	PolyFillType	_subFillType;
};

//------------------------------------------------------------------------------

template<typename K = Exact_predicates_inexact_constructions_kernel>
class ClipperOffset
{
	enum NodeType { ntAny, ntOpen, ntClosed };
	using OutPt    = typename Clipper<K>::OutPt;
	using PolyNode = typename Clipper<K>::PolyNode;
	using PolyTree = typename Clipper<K>::PolyTree;

	typedef typename Clipper<K>::Point2 Point2;
	typedef typename K::Point_2 DoublePoint;

public:
	typedef std::vector<typename K::Point_2>  Path;
	typedef std::vector<Path> Paths;
public:
	ClipperOffset(double miterLimit = 2.0, double roundPrecision = 0.25);
	~ClipperOffset();
	void	addPath(const Path& path, JointType joinType, EndType endType);
	void	addPaths(const Paths& paths, JointType joinType, EndType endType);
	void	execute(Paths& solution, double delta);
	void	execute(PolyTree& solution, double delta);
	void	clear();
private:

	static inline void reversePath(Path& p);
	static inline void reversePaths(Paths& p);
	static inline bool orientation(const Path &p);
	static inline void simplifyPolygon(const Path& in_poly, Paths& out_polys, PolyFillType fillType);
	static inline void simplifyPolygons(const Paths& in_polys, Paths& out_polys, PolyFillType fillType);
	static inline void simplifyPolygons(Paths& polys, PolyFillType fillType);
	static inline bool slopesNearCollinear(const Point2& pt1, const Point2& pt2, const Point2& pt3, double distSqrd);
	static inline bool pointsAreClose(Point2 pt1, Point2 pt2, double distSqrd);
	static inline void cleanPolygon(Path& poly, double distance);
	static inline void cleanPolygon(const Path& in_poly, Path& out_poly, double distance);
	static inline void cleanPolygons(const Paths& in_polys, Paths& out_polys, double distance);
	static inline void cleanPolygons(Paths& polys, double distance);
	static inline void minkowski(const Path& poly, const Path& path, Paths& solution, bool isSum, bool isClosed);
	static inline void minkowskiSum(const Path& pattern, const Path& path, Paths& solution, bool pathIsClosed);
	static inline void minkowskiSum(const Path& pattern, const Paths& paths, Paths& solution, bool pathIsClosed);
	static inline void minkowskiDiff(const Path& poly1, const Path& poly2, Paths& solution);
	static inline void translatePath(const Path& input, Path& output, const Point2 delta);
	static inline void addPolyNodeToPaths(const PolyNode& polynode, NodeType nodetype, Paths& paths);
	static inline void polyTreeToPaths(const PolyTree& polytree, Paths& paths);
	static inline void closedPathsFromPolyTree(const PolyTree& polytree, Paths& paths);
	static inline void openPathsFromPolyTree(PolyTree& polytree, Paths& paths);

	static inline double distanceSqrt(const Point2& pt1, const Point2& pt2);
	static inline double distanceFromLineSqrd(const Point2& pt, const Point2& ln1, const Point2& ln2);

	static OutPt*		excludeOp(OutPt* op);
	static DoublePoint	getUnitNormal(const Point2& pt1, const Point2& pt2);

	void	fixOrientations();
	void	doOffset(double delta);
	void	offsetPoint(int j, int& k, JointType jointype);
	void	doSquare(int j, int k);
	void	doMiter(int j, int k, double r);
	void	doRound(int j, int k);

	double		_miterLimit;
	double		_arcTolerance;
	Paths		_destPolys;
	Path		_srcPoly;
	Path		_destPoly;
	double		_delta, _sinA, _sin, _cos;
	double		_miterLim, _stepsPerRad;
	Point2		_lowest;
	PolyNode	_polyNodes;

	std::vector<DoublePoint>	_normals;

	static constexpr double pi = 3.141592653589793238;
	static constexpr double pi2 = pi *2;
	static constexpr double arcTolerance = 0.25;
};

//------------------------------------------------------------------------------

template<typename K>
Clipper<K>::Clipper()
	: _openPath(false)
  , _strictSimple(false)
  , _excuteLocked(false)
  , _sortedEdges(nullptr)
  , _activeEdges(nullptr)
{
}

template<typename K>
Clipper<K>::~Clipper() {}


template<typename K>
bool Clipper<K>::addPath(const Path &pg, PolyType type, bool closed)
{
    if (!closed && type == PolyType::ptClip)
        throw ClipperException("addPath: open paths must be subject.");

    using std::vector;
    int highIdx = pg.size() - 1;
    if (closed)
        while (highIdx > 0 && (pg[highIdx] == pg[0]))
            --highIdx;
    while (highIdx > 0 && (pg[highIdx] == pg[highIdx - 1]))
        --highIdx;
    if (closed && highIdx < 2 || !closed && highIdx < 1)
        return false;

    std::vector<TEdge> edges(highIdx + 1);
    bool isFlat = true;

    edges[1].cur = pg[1];
    initEdge(&edges[0], &edges[1], &edges[highIdx], pg[0]);
    initEdge(&edges[highIdx], &edges[0], &edges[highIdx - 1], pg[highIdx]);
    for (int i = highIdx - 1; i >= 1; --i)
    {
        initEdge(&edges[i], &edges[i + 1], &edges[i - 1], pg[i]);
    }

    TEdge *eStart = &edges[0];

    TEdge *e = eStart, *eLoopStop = eStart;
    for (;;)
    {
        if (e->cur == e->next->cur && (closed || e->next != eStart))
        {
            if (e == e->next) 
				break;
            if (e == eStart) 
				eStart = e->next;
            e = removeEdge(e);
            eLoopStop = e;
            continue;
        }
        if (e->prev == e->next)
            break; //only two vertices

        else if (closed && slopesEqual(e->prev->cur, e->cur, e->next->cur) && !pt2BetweenPt1AndPt3(e->prev->cur, e->cur, e->next->cur))
        {
            if (e == eStart)
                eStart = e->next;
            e = removeEdge(e);
            e = e->prev;
            eLoopStop = e;
            continue;
        }
        e = e->next;
        if (e == eLoopStop || !closed && e->next == eStart)
            break;
    }

    if (!closed && e == e->next || closed && e->prev == e->next)
        return false;

    if (!closed) {
        _openPath = true;
        eStart->prev->outIdx = _skip;
    }

    //3. Do second stage of edge initialization ...
    e = eStart;
    do
    {
        initEdge2(*e, type);
        e = e->next;
        if (isFlat && e->cur.y() != eStart->cur.y())
            isFlat = false;
    } while (e != eStart);

    //4. Finally, add edge bounds to LocalMinima list ...

    //Totally flat paths must be handled differently when adding them
    //to LocalMinima list to avoid endless loops etc ...
	if (isFlat)
	{
		e->prev->outIdx = _skip;
		MinimaList::value_type locMin;
		locMin.y = e->bot().y();
		locMin.leftBound = 0;
		locMin.rightBound = e;
		locMin.rightBound->side = EdgeSide::esRight;
		locMin.rightBound->windDelta = 0;
		for (;;)
		{
			if (e->bot().x() != e->prev->top().x())
				reverseHorizontal(*e);
			if (e->next->outIdx == _skip)
				break;
			e->nextInLML = e->next;
			e = e->next;
		}
		_minimaList.push_back(locMin);
		_edges.push_back(edges);
		return true;
	}

    _edges.push_back(std::move(edges));
    bool leftBoundForward;
    TEdge *eMin = 0;

    //workaround to avoid an endless loop in the while loop below when
    //open paths have matching start and end points ...
    if (e->prev->bot() == e->prev->top())
        e = e->next;

    for (;;)
    {
        e = findNextLocMin(e);
        if (e == eMin)
            break;
        else if (!eMin)
            eMin = e;

        //E and E.prev now share a local minima (left aligned if horizontal).
        //Compare their slopes to find which starts which bound ...
        MinimaList::value_type locMin;
        locMin.y = e->bot().y();
        if (e->dx > e->prev->dx)
        {
            locMin.leftBound = e->prev;
            locMin.rightBound = e;
            leftBoundForward = false;
        }
        else
        {
            locMin.leftBound = e;
            locMin.rightBound = e->prev;
            leftBoundForward = true;
        }

        if (!closed)
            locMin.leftBound->windDelta = 0;
        else 
		if (locMin.leftBound->next == locMin.rightBound)
            locMin.leftBound->windDelta =  1;
        else
            locMin.leftBound->windDelta = -1;
        locMin.rightBound->windDelta = -locMin.leftBound->windDelta;

        e = processBound(locMin.leftBound, leftBoundForward);
        if (e->outIdx == _skip)
            e = processBound(e, leftBoundForward);

        TEdge *e2 = processBound(locMin.rightBound, !leftBoundForward);
        if (e2->outIdx == _skip)
            e2 = processBound(e2, !leftBoundForward);

        if (locMin.leftBound->outIdx == _skip)
            locMin.leftBound = 0;
        else if (locMin.rightBound->outIdx == _skip)
            locMin.rightBound = 0;
        _minimaList.push_back(locMin);
        if (!leftBoundForward)
            e = e2;
    }
    return true;
}

template<typename K>
bool Clipper<K>::addPaths(const Paths& pg, PolyType PolyTyp, bool closed)
{
	bool result = false;
	for (Paths::size_type i = 0; i < pg.size(); ++i)
		if (addPath(pg[i], PolyTyp, closed)) 
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
    if (_excuteLocked)
        return false;
    if (_openPath)
        throw ClipperException("Error: PolyTree struct is needed for open path clipping.");
    _excuteLocked = true;
    _subFillType = subjFillType;
    _clipFillType = clipFillType;
    _clipType = clipType;
    _usingPolyTree = false;
    solution.resize(0);
    bool succeeded = executeInternal();
    if (succeeded)
        buildResult(solution);
    disposeAllOutRecs();
    _excuteLocked = false;
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
	if (_excuteLocked)
		return false;
	_excuteLocked = true;
	_subFillType = subjFillType;
	_clipFillType = clipFillType;
	_clipType = clipType;
	_usingPolyTree = true;

	bool succeeded = executeInternal();
	if (succeeded) 
		buildResult2(polytree);
	disposeAllOutRecs();
	_excuteLocked = false;
	return succeeded;
}

template<typename K>
void Clipper<K>::clear()
{
    disposeLocalMinimaList();
    _edges.clear();
    _openPath = false;
    _useFullRange = false;
}

template<typename K>
Rectangle2<K> Clipper<K>::getBounds()
{
    Rectangle2<K> result;
    MinimaList::iterator lm = _minimaList.begin();
    if (lm == _minimaList.end())
    {
        result.left = result.top = result.right = result.bottom = 0;
        return result;
    }
    result.left = lm->leftBound->bot().x();
    result.top = lm->leftBound->bot().y();
    result.right = lm->leftBound->bot().x();
    result.bottom = lm->leftBound->bot().y();
    while (lm != _minimaList.end())
    {
        //todo - needs fixing for open paths
        result.bottom = std::max(result.bottom, lm->leftBound->bot().y());
        TEdge *e = lm->leftBound;
        for (;;) {
            TEdge *bottomE = e;
            while (e->nextInLML)
            {
                if (e->bot().x() < result.left)
					result.left = e->bot().x();
                if (e->bot().x() > result.right) 
					result.right = e->bot().x();
                e = e->nextInLML;
            }
            result.left = std::min(result.left, e->bot().x());
            result.right = std::max(result.right, e->bot().x());
            result.left = std::min(result.left, e->top().x());
            result.right = std::max(result.right, e->top().x());
            result.top = std::min(result.top, e->top().y());
            if (bottomE == lm->leftBound) 
				e = lm->rightBound;
            else break;
        }
        ++lm;
    }
    return result;
}

template<typename K>
bool Clipper<K>::isHorizontal(TEdge &e)
{
    return e.dx == -DBL_MAX;
}

template<typename K>
void Clipper<K>::setDelta(TEdge &e)
{
    if (e.is_horizontal())
        e.dx = -DBL_MAX;
    else
        e.dx = (e.top().x() - e.bot().x()) / (e.top().y() - e.bot().y());
}

template<typename K>
double Clipper<K>::getDelta(const Point2& pt1, const Point2& pt2)
{
  return (pt1.y() == pt2.y()) ?  -DBL_MAX : (double)(pt2.x() - pt1.x()) / (pt2.y() - pt1.y());
}

template<typename K>
void Clipper<K>::initEdge(TEdge *e, TEdge *enext, TEdge *eprev, const Point2 &pt)
{
    std::memset(e, 0, sizeof(TEdge));
    e->next = enext;
    e->prev = eprev;
    e->cur = pt;
    e->outIdx = _unassigned;
}

template<typename K>
void Clipper<K>::initEdge2(TEdge &e, PolyType type)
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
    setDelta(e);
    e.polyType = type;
}

template<typename K>
bool Clipper<K>::slopesEqual(const Point2 &pt1, const Point2 &pt2, const Point2 &pt3)
{
    return compare((pt1.y() - pt2.y()) * (pt2.x() - pt3.x()), (pt1.x() - pt2.x()) * (pt2.y() - pt3.y())) == EQUAL;
}

template<typename K>
bool Clipper<K>::slopesEqual(const Point2 &pt1, const Point2 &pt2, const Point2 &pt3, const Point2 &pt4)
{
    return (pt1.y() - pt2.y()) * (pt3.x() - pt4.x()) == (pt1.x() - pt2.x()) * (pt3.y() - pt4.y());
}

template<typename K>
bool Clipper<K>::slopesEqual(const TEdge& e1, const TEdge& e2)
{
    return (e1.top().y() - e1.bot().y()) * (e2.top().x() - e2.bot().x()) == 
		(e1.top().x() - e1.bot().x()) * (e2.top().y() - e2.bot().y());
}

template<typename K>
bool Clipper<K>::pt2BetweenPt1AndPt3(const Point2 pt1, const Point2 pt2, const Point2 pt3)
{
    if ((pt1 == pt3) || (pt1 == pt2) || (pt3 == pt2))
        return false;
    else if (pt1.x() != pt3.x())
        return (pt2.x() > pt1.x()) == (pt2.x() < pt3.x());
    else
        return (pt2.y() > pt1.y()) == (pt2.y() < pt3.y());
}

template<typename K>
void Clipper<K>::reverseHorizontal(TEdge &e)
{
    //swap horizontal edges' top() and Bottom x's so they follow the natural
    //progression of the bounds - ie so their xbots will align with the
    //adjoining lower edge. [Helpful in the processHorizontal() method.]
    std::swap(e.top(), e.bot());
}

template<typename K>
void Clipper<K>::reversePolyPtLinks(OutPt *pp)
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
void Clipper<K>::disposeOutPts(OutPt *&pp)
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
bool Clipper<K>::horzSegmentsOverlap(FT seg1a, FT seg1b, FT seg2a, FT seg2b)
{
	if (seg1a > seg1b)
		std::swap(seg1a, seg1b);
	if (seg2a > seg2b)
		std::swap(seg2a, seg2b);
	return (seg1a < seg2b) && (seg2a < seg1b);
}

template<typename K>
void Clipper<K>::swapSides(TEdge &e1, TEdge &e2)
{
    EdgeSide side = e1.side;
    e1.side = e2.side;
    e2.side = side;
}

template<typename K>
void Clipper<K>::swapPolyIndexes(TEdge &e1, TEdge &e2)
{
    int outIdx = e1.outIdx;
    e1.outIdx = e2.outIdx;
    e2.outIdx = outIdx;
}

template<typename K>
bool Clipper<K>::e2InsertsBeforeE1(TEdge &e1, TEdge &e2)
{
    if (e2.cur.x() == e1.cur.x())
    {
        if (e2.top().y() > e1.top().y())
            return e2.top().x() < topX(e1, e2.top().y());
        else
            return e1.top().x() > topX(e2, e1.top().y());
    }
    else
        return e2.cur.x() < e1.cur.x();
}

template<typename K>
bool Clipper<K>::isMinima(TEdge *e, const FT y)
{
    return e && (e->prev->nextInLML != e) && (e->next->nextInLML != e);
}

template<typename K>
bool Clipper<K>::isMaxima(TEdge *e, const FT y)
{
    return e && e->top().y() == y && !e->nextInLML;
}

template<typename K>
bool Clipper<K>::isIntermediate(TEdge *e, const FT y)
{
    return e->top().y() == y && e->nextInLML;
}

template<typename K>
bool Clipper<K>::edgesAdjacent(const IntersectNode &inode)
{
    return (inode.edge1->nextInSEL == inode.edge2) ||
        (inode.edge1->prevInSEL == inode.edge2);
}

template<typename K>
typename Clipper<K>::FT
Clipper<K>::topX(TEdge &edge, const FT currentY)
{
    return (currentY == edge.top().y()) ?
        edge.top().x() : edge.bot().x() + (edge.dx * (currentY - edge.bot().y()));
}

template<typename K>
typename Clipper<K>::TEdge *
Clipper<K>::removeEdge(TEdge *e)
{
    e->prev->next = e->next;
    e->next->prev = e->prev;
    TEdge *result = e->next;
    e->prev = 0;
    return result;
}

template<typename K>
typename Clipper<K>::TEdge *
Clipper<K>::getMaximaPair(TEdge *e)
{
    if ((e->next->top() == e->top()) && !e->next->nextInLML)
        return e->next;
    else if ((e->prev->top() == e->top()) && !e->prev->nextInLML)
        return e->prev;
    else 
		return 0;
}

template<typename K>
typename Clipper<K>::TEdge *
Clipper<K>::getMaximaPairEx(TEdge *e)
{
    //as getMaximaPair() but returns 0 if MaxPair isn't in AEL (unless it's horizontal)
    TEdge *result = getMaximaPair(e);
    if (result && (result->outIdx == _skip ||
        (result->nextInAEL == result->prevInAEL && !isHorizontal(*result))))
        return 0;
    return result;
}

template<typename K>
typename Clipper<K>::TEdge * 
Clipper<K>::getNextInAEL(TEdge* e, Orientation dir)
{
	return dir == RIGHT_TURN ? e->nextInAEL : e->prevInAEL;
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
double Clipper<K>::area(const Path& path)
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
typename Clipper<K>::TEdge *
Clipper<K>::findNextLocMin(TEdge *e)
{
    for (;;)
    {
        while (e->bot() != e->prev->bot() || e->cur == e->top())
            e = e->next;
        if (!isHorizontal(*e) && !isHorizontal(*e->prev))
            break;
        while (isHorizontal(*e->prev))
            e = e->prev;
        TEdge *e2 = e;
        while (isHorizontal(*e))
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
Clipper<K>::processBound(TEdge *e, bool nextForward)
{
    TEdge *result = e, *horz = 0;

    if (e->outIdx == _skip)
    {
    	//if edges still remain in the current bound beyond the skip edge then
    	//create another LocMin and call processBound once more
    	if (nextForward)
    	{
    		while (e->top().y() == e->next->bot().y())
				e = e->next;
    		//don't include top() horizontals when parsing a bound a second time,
    		//they will be contained in the opposite bound ...
    		while (e != result && isHorizontal(*e))
				e = e->prev;
    	}
    	else
    	{
    		while (e->top().y() == e->prev->bot().y())
				e = e->prev;
    		while (e != result && isHorizontal(*e))
				e = e->next;
    	}

    	if (e == result)
    	{
    		if (nextForward) 
				result = e->next;
    		else 
				result = e->prev;
    	}
    	else
    	{
    		//there are more edges in the bound beyond result starting with e
    		if (nextForward)
    			e = result->next;
    		else
    			e = result->prev;
    		MinimaList::value_type locMin;
    		locMin.y = e->bot().y();
    		locMin.leftBound = 0;
    		locMin.rightBound = e;
    		e->windDelta = 0;
    		result = processBound(e, nextForward);
    		_minimaList.push_back(locMin);
    	}
    	return result;
    }

    TEdge *eStart;

    if (isHorizontal(*e)) {
        //We need to be careful with open paths because this may not be a
        //true local minima (ie e may be following a skip edge).
        //Also, consecutive horz. edges may start heading left before going right.
        if (nextForward)
            eStart = e->prev;
        else
            eStart = e->next;
        if (isHorizontal(*eStart)) //ie an adjoining horizontal skip edge
        {
            if (eStart->bot().x() != e->bot().x() && eStart->top().x() != e->bot().x())
                reverseHorizontal(*e);
        }
        else if (eStart->bot().x() != e->bot().x())
            reverseHorizontal(*e);
    }

    eStart = e;
    if (nextForward)
    {
        while (result->top().y() == result->next->bot().y() && result->next->outIdx != _skip)
            result = result->next;
        if (isHorizontal(*result) && result->next->outIdx != _skip)
        {
            //nb: at the top() of a bound, horizontals are added to the bound
            //only when the preceding edge attaches to the horizontal's left vertex
            //unless a _skip edge is encountered when that becomes the top() divide
            horz = result;
            while (isHorizontal(*horz->prev))
				horz = horz->prev;
            if (horz->prev->top().x() > result->next->top().x()) 
				result = horz->prev;
        }
        while (e != result)
        {
            e->nextInLML = e->next;
            if (isHorizontal(*e) && e != eStart && e->bot().x() != e->prev->top().x())
				reverseHorizontal(*e);
            e = e->next;
        }
        if (isHorizontal(*e) && e != eStart && e->bot().x() != e->prev->top().x())
            reverseHorizontal(*e);
        result = result->next; //move to the edge just beyond current bound
    }
    else {
        while (result->top().y() == result->prev->bot().y() && result->prev->outIdx != _skip)
            result = result->prev;
        if (isHorizontal(*result) && result->prev->outIdx != _skip)
        {
            horz = result;
            while (isHorizontal(*horz->next)) 
				horz = horz->next;
            if (horz->next->top().x() == result->prev->top().x() || horz->next->top().x() > result->prev->top().x())
				result = horz->next;
        }

        while (e != result)
        {
            e->nextInLML = e->prev;
            if (isHorizontal(*e) && e != eStart && e->bot().x() != e->next->top().x())
                reverseHorizontal(*e);
            e = e->prev;
        }
        if (isHorizontal(*e) && e != eStart && e->bot().x() != e->next->top().x())
            reverseHorizontal(*e);
        result = result->prev; //move to the edge just beyond current bound
    }

    return result;
}

template<typename K>
bool Clipper<K>::executeInternal()
{
    bool succeeded = true;
    try {
        reset();
        _maxima = MaximaList();
        _sortedEdges = 0;

        succeeded = true;
        CT botY, topY;
        if (!popScanbeam(botY))
            return false;
        insertLocalMinimaIntoAEL(botY);
        while (popScanbeam(topY) || localMinimaPending())
        {
            processHorizontals();
            clearGhostJoints();
            if (!processIntersections(topY))
            {
                succeeded = false;
                break;
            }
            processEdgesAtTopOfScanbeam(topY);
            botY = topY;
            insertLocalMinimaIntoAEL(botY);
        }
    }
    catch (...)
    {
        succeeded = false;
    }

    if (succeeded)
    {
        //fix orientations ...
        for (PolyOutList::size_type i = 0; i < _polyOuts.size(); ++i)
        {
            OutRec *outRec = _polyOuts[i];
            if (!outRec->pts || outRec->isOpen)
                continue;
            if ((outRec->isHole ^ _reverseOutput) == (area(outRec->pts) > 0))
                reversePolyPtLinks(outRec->pts);
        }

        if (!_joints.empty())
            joinCommonEdges();

        //unfortunately fixupOutPolygon() must be done after joinCommonEdges()
        for (PolyOutList::size_type i = 0; i < _polyOuts.size(); ++i) {
            OutRec *outRec = _polyOuts[i];
            if (!outRec->pts)
                continue;

            if (outRec->isOpen)
                fixupOutPolyline(*outRec);
            else
                fixupOutPolygon(*outRec);
        }

        if (_strictSimple)
            doSimplePolygons();
    }

    clearJoints();
    clearGhostJoints();
    return succeeded;
}

template<typename K>
void Clipper<K>::reset()
{
    if (_minimaList.empty())
        return; //ie nothing to process
    std::sort(_minimaList.begin(), _minimaList.end(),
        [](const LocalMinimum &locMin1, const LocalMinimum &locMin2) {
            return locMin1.y < locMin2.y;
        });

    _scanbeam = ScanbeamList(); //clears/resets priority_queue
    //reset all edges ...
    for (MinimaList::iterator lm = _minimaList.begin(); lm != _minimaList.end(); ++lm)
    {
        insertScanbeam(lm->y);
        TEdge *e = lm->leftBound;
        if (e)
        {
            e->cur = e->bot();
            e->side = EdgeSide::esLeft;
            e->outIdx = _unassigned;
        }

        e = lm->rightBound;
        if (e)
        {
            e->cur = e->bot();
            e->side = EdgeSide::esRight;
            e->outIdx = _unassigned;
        }
    }
    _activeEdges = 0;
    _currentLM = _minimaList.begin();

#ifdef DEBUG
	_activeCount = 0;
	_sortedCount = 0;
#endif
}

template<typename K>
void Clipper<K>::disposeLocalMinimaList()
{
    _minimaList.clear();
    _currentLM = _minimaList.begin();
}

template<typename K>
bool Clipper<K>::popLocalMinima(CT y, const LocalMinimum *&locMin)
{
    if (_currentLM == _minimaList.end() || (*_currentLM).y != y)
        return false;
    locMin = &(*_currentLM);
    ++_currentLM;
    return true;
}

template<typename K>
void Clipper<K>::insertScanbeam(const CT y)
{
    _scanbeam.push(y);
}

template<typename K>
bool Clipper<K>::popScanbeam(CT &y)
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
void Clipper<K>::disposeAllOutRecs() 
{
    for (PolyOutList::size_type i = 0; i < _polyOuts.size(); ++i)
        disposeOutRec(i);
    _polyOuts.clear();
}

template<typename K>
void Clipper<K>::disposeOutRec(typename PolyOutList::size_type index)
{
    OutRec *outRec = _polyOuts[index];
    if (outRec->pts) 
		disposeOutPts(outRec->pts);
    delete outRec;
    _polyOuts[index] = 0;
}

template<typename K>
void Clipper<K>::swapPositionsInAEL(TEdge *edge1, TEdge *edge2)
{
    //check that one or other edge hasn't already been removed from AEL ...
    if (edge1->nextInAEL == edge1->prevInAEL ||
        edge2->nextInAEL == edge2->prevInAEL)
        return;

    if (edge1->nextInAEL == edge2)
    {
        TEdge *next = edge2->nextInAEL;
        if (next)
            next->prevInAEL = edge1;
        TEdge *prev = edge1->prevInAEL;
        if (prev)
            prev->nextInAEL = edge2;
        edge2->prevInAEL = prev;
        edge2->nextInAEL = edge1;
        edge1->prevInAEL = edge2;
        edge1->nextInAEL = next;
    }
    else if (edge2->nextInAEL == edge1)
    {
        TEdge *next = edge1->nextInAEL;
        if (next)
            next->prevInAEL = edge2;
        TEdge *prev = edge2->prevInAEL;
        if (prev)
            prev->nextInAEL = edge1;
        edge1->prevInAEL = prev;
        edge1->nextInAEL = edge2;
        edge2->prevInAEL = edge1;
        edge2->nextInAEL = next;
    }
    else
    {
        TEdge *next = edge1->nextInAEL;
        TEdge *prev = edge1->prevInAEL;
        edge1->nextInAEL = edge2->nextInAEL;
        if (edge1->nextInAEL)
            edge1->nextInAEL->prevInAEL = edge1;
        edge1->prevInAEL = edge2->prevInAEL;
        if (edge1->prevInAEL)
            edge1->prevInAEL->nextInAEL = edge1;
        edge2->nextInAEL = next;
        if (edge2->nextInAEL)
            edge2->nextInAEL->prevInAEL = edge2;
        edge2->prevInAEL = prev;
        if (edge2->prevInAEL)
            edge2->prevInAEL->nextInAEL = edge2;
    }

    if (!edge1->prevInAEL)
        _activeEdges = edge1;
    else if (!edge2->prevInAEL)
        _activeEdges = edge2;
}

template<typename K>
void Clipper<K>::deleteFromAEL(TEdge *e)
{
    TEdge *aelPrev = e->prevInAEL;
    TEdge *aelNext = e->nextInAEL;
    if (!aelPrev && !aelNext && (e != _activeEdges))
		return; //already deleted
    if (aelPrev) 
		aelPrev->nextInAEL = aelNext;
    else 
		_activeEdges = aelNext;
    if (aelNext) 
		aelNext->prevInAEL = aelPrev;
    e->nextInAEL = 0;
    e->prevInAEL = 0;
}

template<typename K>
typename Clipper<K>::OutRec *
Clipper<K>::createOutRec()
{
    OutRec *result = new OutRec;
    result->isHole = false;
    result->isOpen = false;
    result->firstLeft = 0;
    result->pts = 0;
    result->botPt = 0;
    result->polyNode = 0;
    _polyOuts.push_back(result);
    result->idx = (int)_polyOuts.size() - 1;
    return result;
}

template<typename K>
void Clipper<K>::updateEdgeIntoAEL(TEdge *&e)
{
    if (!e->nextInLML)
        throw ClipperException("updateEdgeIntoAEL: invalid call");

    e->nextInLML->outIdx = e->outIdx;
    TEdge *aelPrev = e->prevInAEL;
    TEdge *aelNext = e->nextInAEL;
    if (aelPrev)
        aelPrev->nextInAEL = e->nextInLML;
    else
        _activeEdges = e->nextInLML;
    if (aelNext)
        aelNext->prevInAEL = e->nextInLML;
    e->nextInLML->side = e->side;
    e->nextInLML->windDelta = e->windDelta;
    e->nextInLML->windCnt = e->windCnt;
    e->nextInLML->windCnt2 = e->windCnt2;
    e = e->nextInLML;
    e->cur = e->bot();
    e->prevInAEL = aelPrev;
    e->nextInAEL = aelNext;
    if (!isHorizontal(*e))
        insertScanbeam(e->top().y());
}

template<typename K>
bool Clipper<K>::localMinimaPending()
{
    return (_currentLM != _minimaList.end());
}

template<typename K>
void Clipper<K>::setWindingCount(TEdge &edge)
{
    TEdge *e = edge.prevInAEL;
    //find the edge of the same polytype that immediately preceeds 'edge' in AEL
    while (e && ((e->polyType != edge.polyType) || (e->windDelta == 0)))
        e = e->prevInAEL;
    if (!e)
    {
        if (edge.windDelta == 0)
        {
            PolyFillType pft = edge.polyType == PolyType::ptSubject ? _subFillType : _clipFillType;
            edge.windCnt = (pft == PolyFillType::pftPositive ? 1 : -1);
        }
        else
            edge.windCnt = edge.windDelta;

        edge.windCnt2 = 0;
        e = _activeEdges; //ie get ready to calc windCnt2
    }
    else if (edge.windDelta == 0 && _clipType != ClipType::ctUnion)
    {
        edge.windCnt = 1;
        edge.windCnt2 = e->windCnt2;
        e = e->nextInAEL; //ie get ready to calc windCnt2
    }
    else if (isEvenOddFillType(edge))
    {
        //EvenOdd filling ...
        if (edge.windDelta == 0)
        {
            //are we inside a subj polygon ...
            bool inside = true;
            TEdge *e2 = e->prevInAEL;
            while (e2) {
                if (e2->polyType == e->polyType && e2->windDelta != 0)
                    inside = !inside;
                e2 = e2->prevInAEL;
            }
            edge.windCnt = (inside ? 0 : 1);
        }
        else
        {
            edge.windCnt = edge.windDelta;
        }
        edge.windCnt2 = e->windCnt2;
        e = e->nextInAEL; //ie get ready to calc windCnt2
    }
    else
    {
        //nonZero, Positive or Negative filling ...
        if (e->windCnt * e->windDelta < 0)
        {
            //prev edge is 'decreasing' WindCount (WC) toward zero
            //so we're outside the previous polygon ...
            if (std::abs(e->windCnt) > 1)
            {
                //outside prev poly but still inside another.
                //when reversing direction of prev poly use the same WC 
                if (e->windDelta * edge.windDelta < 0)
                    edge.windCnt = e->windCnt;
                //otherwise continue to 'decrease' WC ...
                else
                    edge.windCnt = e->windCnt + edge.windDelta;
            }
            else
                //now outside all polys of same polytype so set own WC ...
                edge.windCnt = (edge.windDelta == 0 ? 1 : edge.windDelta);
        }
        else
        {
            //prev edge is 'increasing' WindCount (WC) away from zero
            //so we're inside the previous polygon ...
            if (edge.windDelta == 0)
                edge.windCnt = (e->windCnt < 0 ? e->windCnt - 1 : e->windCnt + 1);
            //if wind direction is reversing prev then use same WC
            else if (e->windDelta * edge.windDelta < 0)
                edge.windCnt = e->windCnt;
            //otherwise add to WC ...
            else
                edge.windCnt = e->windCnt + edge.windDelta;
        }
        edge.windCnt2 = e->windCnt2;
        e = e->nextInAEL; //ie get ready to calc windCnt2
    }

    //update windCnt2 ...
    if (isEvenOddAltFillType(edge))
    {
        //EvenOdd filling ...
        while (e != &edge)
        {
            if (e->windDelta != 0)
                edge.windCnt2 = (edge.windCnt2 == 0 ? 1 : 0);
            e = e->nextInAEL;
        }
    }
    else
    {
        //nonZero, Positive or Negative filling ...
        while (e != &edge)
        {
            edge.windCnt2 += e->windDelta;
            e = e->nextInAEL;
        }
    }
}

template<typename K>
bool Clipper<K>::isEvenOddFillType(const TEdge &edge) const
{
    if (edge.polyType == PolyType::ptSubject)
        return _subFillType == PolyFillType::pftEvenOdd;
    else
        return _clipFillType == PolyFillType::pftEvenOdd;
}

template<typename K>
bool Clipper<K>::isEvenOddAltFillType(const TEdge &edge) const
{
    if (edge.polyType == PolyType::ptSubject)
        return _clipFillType == PolyFillType::pftEvenOdd;
    else
        return _subFillType == PolyFillType::pftEvenOdd;
}

template<typename K>
void Clipper<K>::getHorzDirection(TEdge& horzEdge, Orientation& dir, CT& left, CT& right)
{
	if (horzEdge.bot().x() < horzEdge.top().x())
	{
		left = horzEdge.bot().x();
		right = horzEdge.top().x();
		dir = RIGHT_TURN;
	}
	else
	{
		left = horzEdge.top().x();
		right = horzEdge.bot().x();
		dir = LEFT_TURN;
	}
}

template<typename K>
bool Clipper<K>::getOverlap(const FT a1, const FT a2, const FT b1, const FT b2, FT& left, FT& right)
{
	if (a1 < a2)
	{
		if (b1 < b2) { left = std::max(a1, b1); right = std::min(a2, b2); }
		else { left = std::max(a1, b2); right = std::min(a2, b1); }
	}
	else
	{
		if (b1 < b2) { left = std::max(a2, b1); right = std::min(a1, b2); }
		else { left = std::max(a2, b2); right = std::min(a1, b1); }
	}
	return left < right;
}

template<typename K>
void Clipper<K>::insertLocalMinimaIntoAEL(const FT botY)
{
    const LocalMinimum *lm;
    while (popLocalMinima(botY, lm))
    {
        TEdge *lb = lm->leftBound;
        TEdge *rb = lm->rightBound;

        OutPt *op1 = 0;
        if (!lb)
        {
            //nb: don't insert LB into either AEL or SEL
            insertEdgeIntoAEL(rb, 0);
            setWindingCount(*rb);
            if (isContributing(*rb))
                op1 = addOutPt(rb, rb->bot());
        }
        else if (!rb)
        {
            insertEdgeIntoAEL(lb, 0);
            setWindingCount(*lb);
            if (isContributing(*lb))
                op1 = addOutPt(lb, lb->bot());
            insertScanbeam(lb->top().y());
        }
        else
        {
            insertEdgeIntoAEL(lb, 0);
            insertEdgeIntoAEL(rb, lb);
            setWindingCount(*lb);
            rb->windCnt = lb->windCnt;
            rb->windCnt2 = lb->windCnt2;
            if (isContributing(*lb))
                op1 = addLocalMinPoly(lb, rb, lb->bot());
            insertScanbeam(lb->top().y());
        }

        if (rb)
        {
            if (isHorizontal(*rb))
            {
                addEdgeToSEL(rb);
                if (rb->nextInLML)
                    insertScanbeam(rb->nextInLML->bot().y());
            }
            else
                insertScanbeam(rb->top().y());
        }

        if (!lb || !rb)
            continue;

        //if any output polygons share an edge, they'll need joining later ...
        if (op1 && isHorizontal(*rb) && _ghostJoints.size() > 0 && (rb->windDelta != 0))
        {
            for (JoinList::size_type i = 0; i < _ghostJoints.size(); ++i)
            {
                Joint *jr = &_ghostJoints[i];
                //if the horizontal Rb and a 'ghost' horizontal overlap, then convert
                //the 'ghost' join to a real join ready for later ...
                if (horzSegmentsOverlap(jr->outPt1->pt.x(), jr->offPt.x(), rb->bot().x(), rb->top().x()))
                    addJoin(jr->outPt1, op1, jr->offPt);
            }
        }

        if (lb->outIdx >= 0 && lb->prevInAEL && lb->prevInAEL->cur.x() == lb->bot().x() &&
            lb->prevInAEL->outIdx >= 0 && slopesEqual(lb->prevInAEL->bot(), lb->prevInAEL->top(), lb->cur, lb->top()) &&
            (lb->windDelta != 0) && (lb->prevInAEL->windDelta != 0))
        {
            OutPt *Op2 = addOutPt(lb->prevInAEL, lb->bot());
            addJoin(op1, Op2, lb->top());
        }

        if (lb->nextInAEL != rb)
        {

            if (rb->outIdx >= 0 && rb->prevInAEL->outIdx >= 0 &&
                slopesEqual(rb->prevInAEL->cur, rb->prevInAEL->top(), rb->cur, rb->top()) &&
                (rb->windDelta != 0) && (rb->prevInAEL->windDelta != 0))
            {
                OutPt *Op2 = addOutPt(rb->prevInAEL, rb->bot());
                addJoin(op1, Op2, rb->top());
            }

            TEdge *e = lb->nextInAEL;
            if (e)
            {
                while (e != rb)
                {
                    //nb: For calculating winding counts etc, intersectEdges() assumes
                    //that param1 will be to the right of param2 ABOVE the intersection ...
                    intersectEdges(rb, e, lb->cur); //order important here
                    e = e->nextInAEL;
                }
            }
        }
    }
}

template<typename K>
void Clipper<K>::insertEdgeIntoAEL(TEdge *edge, TEdge *startEdge)
{
    if (!_activeEdges)
    {
        edge->prevInAEL = 0;
        edge->nextInAEL = 0;
        _activeEdges = edge;
    }
    else if (!startEdge && e2InsertsBeforeE1(*_activeEdges, *edge))
    {
        edge->prevInAEL = 0;
        edge->nextInAEL = _activeEdges;
        _activeEdges->prevInAEL = edge;
        _activeEdges = edge;
    }
    else
    {
        if (!startEdge)
            startEdge = _activeEdges;
        while (startEdge->nextInAEL && !e2InsertsBeforeE1(*startEdge->nextInAEL, *edge))
            startEdge = startEdge->nextInAEL;
        edge->nextInAEL = startEdge->nextInAEL;
        if (startEdge->nextInAEL)
            startEdge->nextInAEL->prevInAEL = edge;
        edge->prevInAEL = startEdge;
        startEdge->nextInAEL = edge;
    }
#ifdef DEBUG
	_activeCount++;
#endif
}

template<typename K>
void Clipper<K>::addEdgeToSEL(TEdge *edge)
{
    //SEL pointers in PEdge are reused to build a list of horizontal edges.
    //However, we don't need to worry about order with horizontal edge processing.
    if (!_sortedEdges)
    {
        _sortedEdges = edge;
        edge->prevInSEL = 0;
        edge->nextInSEL = 0;
    }
    else
    {
        edge->nextInSEL = _sortedEdges;
        edge->prevInSEL = 0;
        _sortedEdges->prevInSEL = edge;
        _sortedEdges = edge;
    }

#ifdef DEBUG
	_sortedCount++;
#endif
}

template<typename K>
bool Clipper<K>::popEdgeFromSEL(TEdge *&edge)
{
    if (!_sortedEdges)
        return false;
    edge = _sortedEdges;
    deleteFromSEL(_sortedEdges);
    return true;
}

template<typename K>
void Clipper<K>::copyAELToSEL()
{
    TEdge *e = _activeEdges;
    _sortedEdges = e;
    while (e)
    {
        e->prevInSEL = e->prevInAEL;
        e->nextInSEL = e->nextInAEL;
        e = e->nextInAEL;
    }
}

template<typename K>
void Clipper<K>::deleteFromSEL(TEdge *e)
{
    TEdge *selPrev = e->prevInSEL;
    TEdge *selNext = e->nextInSEL;
    if (!selPrev && !selNext && (e != _sortedEdges)) 
		return; //already deleted
    if (selPrev) 
		selPrev->nextInSEL = selNext;
    else 
		_sortedEdges = selNext;
    if (selNext) 
		selNext->prevInSEL = selPrev;
    e->nextInSEL = 0;
    e->prevInSEL = 0;

#ifdef DEBUG
	_sortedCount--;
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
void Clipper<K>::swapPositionInSEL(TEdge *edge1, TEdge *edge2)
{
    if (!(edge1->nextInSEL) && !(edge1->prevInSEL))
        return;
    if (!(edge2->nextInSEL) && !(edge2->prevInSEL))
        return;

    if (edge1->nextInSEL == edge2)
    {
        TEdge *next = edge2->nextInSEL;
        if (next)
            next->prevInSEL = edge1;
        TEdge *prev = edge1->prevInSEL;
        if (prev)
            prev->nextInSEL = edge2;
        edge2->prevInSEL = prev;
        edge2->nextInSEL = edge1;
        edge1->prevInSEL = edge2;
        edge1->nextInSEL = next;
    }
    else if (edge2->nextInSEL == edge1)
    {
        TEdge *next = edge1->nextInSEL;
        if (next)
            next->prevInSEL = edge2;
        TEdge *prev = edge2->prevInSEL;
        if (prev)
            prev->nextInSEL = edge1;
        edge1->prevInSEL = prev;
        edge1->nextInSEL = edge2;
        edge2->prevInSEL = edge1;
        edge2->nextInSEL = next;
    }
    else
    {
        TEdge *next = edge1->nextInSEL;
        TEdge *prev = edge1->prevInSEL;
        edge1->nextInSEL = edge2->nextInSEL;
        if (edge1->nextInSEL)
            edge1->nextInSEL->prevInSEL = edge1;
        edge1->prevInSEL =
            edge2->prevInSEL;
        if (edge1->prevInSEL)
            edge1->prevInSEL->nextInSEL = edge1;
        edge2->nextInSEL = next;
        if (edge2->nextInSEL)
            edge2->nextInSEL->prevInSEL = edge2;
        edge2->prevInSEL = prev;
        if (edge2->prevInSEL)
            edge2->prevInSEL->nextInSEL = edge2;
    }

    if (!edge1->prevInSEL)
        _sortedEdges = edge1;
    else if (!edge2->prevInSEL)
        _sortedEdges = edge2;
}

template<typename K>
bool Clipper<K>::isContributing(const TEdge &edge) const
{
    PolyFillType pft, pft2;
    if (edge.polyType == PolyType::ptSubject)
    {
        pft = _subFillType;
        pft2 = _clipFillType;
    }
    else
    {
        pft = _clipFillType;
        pft2 = _subFillType;
    }

    switch (pft)
    {
    case PolyFillType::pftEvenOdd:
        //return false if a subj line has been flagged as inside a subj polygon
        if (edge.windDelta == 0 && edge.windCnt != 1)
            return false;
        break;
    case PolyFillType::pftNonZero:
        if (std::abs(edge.windCnt) != 1)
            return false;
        break;
    case PolyFillType::pftPositive:
        if (edge.windCnt != 1)
            return false;
        break;
    default: //PolyFillType::pftNegative
        if (edge.windCnt != -1)
            return false;
    }

    switch (_clipType)
    {
    case ClipType::ctIntersection:
        switch (pft2)
        {
        case PolyFillType::pftEvenOdd:
        case PolyFillType::pftNonZero:
            return (edge.windCnt2 != 0);
        case PolyFillType::pftPositive:
            return (edge.windCnt2 > 0);
        default:
            return (edge.windCnt2 < 0);
        }
        break;
    case ClipType::ctUnion:
        switch (pft2)
        {
        case PolyFillType::pftEvenOdd:
        case PolyFillType::pftNonZero:
            return (edge.windCnt2 == 0);
        case PolyFillType::pftPositive:
            return (edge.windCnt2 <= 0);
        default:
            return (edge.windCnt2 >= 0);
        }
        break;
    case ClipType::ctDifference:
        if (edge.polyType == PolyType::ptSubject)
            switch (pft2)
            {
            case PolyFillType::pftEvenOdd:
            case PolyFillType::pftNonZero:
                return (edge.windCnt2 == 0);
            case PolyFillType::pftPositive:
                return (edge.windCnt2 <= 0);
            default:
                return (edge.windCnt2 >= 0);
            }
        else
            switch (pft2)
            {
            case PolyFillType::pftEvenOdd:
            case PolyFillType::pftNonZero:
                return (edge.windCnt2 != 0);
            case PolyFillType::pftPositive:
                return (edge.windCnt2 > 0);
            default:
                return (edge.windCnt2 < 0);
            }
        break;
    case ClipType::ctXor:
        if (edge.windDelta == 0) //XOr always contributing unless open
            switch (pft2)
            {
            case PolyFillType::pftEvenOdd:
            case PolyFillType::pftNonZero:
                return (edge.windCnt2 == 0);
            case PolyFillType::pftPositive:
                return (edge.windCnt2 <= 0);
            default:
                return (edge.windCnt2 >= 0);
            }
        else
            return true;
        break;
    default:
        return true;
    }
}

template<typename K>
void Clipper<K>::doMaxima(TEdge *e)
{
    TEdge *eMaxPair = getMaximaPairEx(e);
    if (!eMaxPair)
    {
        if (e->outIdx >= 0)
            addOutPt(e, e->top());
        deleteFromAEL(e);
        return;
    }

    TEdge *eNext = e->nextInAEL;
    while (eNext && eNext != eMaxPair)
    {
        intersectEdges(e, eNext, e->top());
        swapPositionsInAEL(e, eNext);
        eNext = e->nextInAEL;
    }

    if (e->outIdx == _unassigned && eMaxPair->outIdx == _unassigned)
    {
        deleteFromAEL(e);
        deleteFromAEL(eMaxPair);
    }
    else if (e->outIdx >= 0 && eMaxPair->outIdx >= 0)
    {
        if (e->outIdx >= 0) 
			addLocalMaxPoly(e, eMaxPair, e->top());
        deleteFromAEL(e);
        deleteFromAEL(eMaxPair);
    }
#ifdef USE_LINES
    else if (e->windDelta == 0)
    {
        if (e->outIdx >= 0)
        {
            addOutPt(e, e->top());
            e->outIdx = _unassigned;
        }
        deleteFromAEL(e);

        if (eMaxPair->outIdx >= 0)
        {
            addOutPt(eMaxPair, e->top());
            eMaxPair->outIdx = _unassigned;
        }
        deleteFromAEL(eMaxPair);
    }
#endif
    else 
		throw ClipperException("doMaxima error");
}

template<typename K>
void Clipper<K>::processHorizontals()
{
    TEdge *horzEdge;
    while (popEdgeFromSEL(horzEdge))
        processHorizontal(horzEdge);
}

template<typename K>
void Clipper<K>::processHorizontal(TEdge *horzEdge)
{
	Orientation dir;
	CT horzLeft, horzRight;
	bool isOpen = (horzEdge->windDelta == 0);

	getHorzDirection(*horzEdge, dir, horzLeft, horzRight);

	TEdge* eLastHorz = horzEdge, * eMaxPair = 0;
	while (eLastHorz->nextInLML && isHorizontal(*eLastHorz->nextInLML))
		eLastHorz = eLastHorz->nextInLML;
	if (!eLastHorz->nextInLML)
		eMaxPair = getMaximaPair(eLastHorz);

	MaximaList::const_iterator maxIt;
	MaximaList::const_reverse_iterator maxRit;
	if (_maxima.size() > 0)
	{
		//get the first maxima in range (x()) ...
		if (dir == RIGHT_TURN)
		{
			maxIt = _maxima.begin();
			while (maxIt != _maxima.end() && *maxIt <= horzEdge->bot().x()) 
				maxIt++;
			if (maxIt != _maxima.end() && *maxIt >= eLastHorz->top().x())
				maxIt = _maxima.end();
		}
		else
		{
			maxRit = _maxima.rbegin();
			while (maxRit != _maxima.rend() && *maxRit > horzEdge->bot().x()) 
				maxRit++;
			if (maxRit != _maxima.rend() && *maxRit <= eLastHorz->top().x())
				maxRit = _maxima.rend();
		}
	}

	OutPt* op1 = 0;

	for (;;) //loop through consec. horizontal edges
	{

		bool isLastHorz = (horzEdge == eLastHorz);
		TEdge* e = getNextInAEL(horzEdge, dir);
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
						if (horzEdge->outIdx >= 0 && !isOpen)
							addOutPt(horzEdge, Point2(*maxIt, horzEdge->bot().y()));
						maxIt++;
					}
				}
				else
				{
					while (maxRit != _maxima.rend() && *maxRit > e->cur.x())
					{
						if (horzEdge->outIdx >= 0 && !isOpen)
							addOutPt(horzEdge, Point2(*maxRit, horzEdge->bot().y()));
						maxRit++;
					}
				}
			};

			if ((dir == RIGHT_TURN && e->cur.x() > horzRight) || (dir == LEFT_TURN && e->cur.x() < horzLeft)) 
				break;

			//Also break if we've got to the end of an intermediate horizontal edge ...
			//nb: Smaller Dx's are to the right of larger Dx's ABOVE the horizontal.
			if (e->cur.x() == horzEdge->top().x() && horzEdge->nextInLML && e->dx < horzEdge->nextInLML->dx)
				break;

			if (horzEdge->outIdx >= 0 && !isOpen)  //note: may be done multiple times
			{
#ifdef use_xyz
				if (dir == RIGHT_TURN) SetZ(e->cur, *horzEdge, *e);
				else SetZ(e->cur, *e, *horzEdge);
#endif      
				op1 = addOutPt(horzEdge, e->cur);
				TEdge* eNextHorz = _sortedEdges;
				while (eNextHorz)
				{
					if (eNextHorz->outIdx >= 0 &&
						horzSegmentsOverlap(horzEdge->bot().x(), horzEdge->top().x(), eNextHorz->bot().x(), eNextHorz->top().x()))
					{
						OutPt* op2 = getLastOutPt(eNextHorz);
						addJoin(op2, op1, eNextHorz->top());
					}
					eNextHorz = eNextHorz->nextInSEL;
				}
				addGhostJoin(op1, horzEdge->bot());
			}

			//OK, so far we're still in range of the horizontal Edge  but make sure
					//we're at the last of consec. horizontals when matching with eMaxPair
			if (e == eMaxPair && isLastHorz)
			{
				if (horzEdge->outIdx >= 0)
					addLocalMaxPoly(horzEdge, eMaxPair, horzEdge->top());
				deleteFromAEL(horzEdge);
				deleteFromAEL(eMaxPair);
				return;
			}

			if (dir == RIGHT_TURN)
			{
				Point2 pt = Point2(e->cur.x(), horzEdge->cur.y());
				intersectEdges(horzEdge, e, pt);
			}
			else
			{
				Point2 pt = Point2(e->cur.x(), horzEdge->cur.y());
				intersectEdges(e, horzEdge, pt);
			}
			TEdge* eNext = getNextInAEL(e, dir);
			swapPositionsInAEL(horzEdge, e);
			e = eNext;
		} //end while(e)

		//Break out of loop if horzEdge.nextInLML is not also horizontal ...
		if (!horzEdge->nextInLML || !isHorizontal(*horzEdge->nextInLML))
			break;

		updateEdgeIntoAEL(horzEdge);
		if (horzEdge->outIdx >= 0) 
			addOutPt(horzEdge, horzEdge->bot());
		getHorzDirection(*horzEdge, dir, horzLeft, horzRight);

	} //end for (;;)

	if (horzEdge->outIdx >= 0 && !op1)
	{
		op1 = getLastOutPt(horzEdge);
		TEdge* eNextHorz = _sortedEdges;
		while (eNextHorz)
		{
			if (eNextHorz->outIdx >= 0 &&
				horzSegmentsOverlap(horzEdge->bot().x(), horzEdge->top().x(), eNextHorz->bot().x(), eNextHorz->top().x()))
			{
				OutPt* op2 = getLastOutPt(eNextHorz);
				addJoin(op2, op1, eNextHorz->top());
			}
			eNextHorz = eNextHorz->nextInSEL;
		}
		addGhostJoin(op1, horzEdge->top());
	}

	if (horzEdge->nextInLML)
	{
		if (horzEdge->outIdx >= 0)
		{
			op1 = addOutPt(horzEdge, horzEdge->top());
			updateEdgeIntoAEL(horzEdge);
			if (horzEdge->windDelta == 0)
				return;
			//nb: horzEdge is no longer horizontal here
			TEdge* ePrev = horzEdge->prevInAEL;
			TEdge* eNext = horzEdge->nextInAEL;
			if (ePrev && ePrev->cur.x() == horzEdge->bot().x() &&
				ePrev->cur.y() == horzEdge->bot().y() && ePrev->windDelta != 0 &&
				(ePrev->outIdx >= 0 && ePrev->cur.y() > ePrev->top().y() &&
					slopesEqual(*horzEdge, *ePrev)))
			{
				OutPt* op2 = addOutPt(ePrev, horzEdge->bot());
				addJoin(op1, op2, horzEdge->top());
			}
			else if (eNext && eNext->cur.x() == horzEdge->bot().x() &&
				eNext->cur.y() == horzEdge->bot().y() && eNext->windDelta != 0 &&
				eNext->outIdx >= 0 && eNext->cur.y() > eNext->top().y() &&
				slopesEqual(*horzEdge, *eNext))
			{
				OutPt* op2 = addOutPt(eNext, horzEdge->bot());
				addJoin(op1, op2, horzEdge->top());
			}
		}
		else
			updateEdgeIntoAEL(horzEdge);
	}
	else
	{
		if (horzEdge->outIdx >= 0)
			addOutPt(horzEdge, horzEdge->top());
		deleteFromAEL(horzEdge);
	}
}

template<typename K>
typename Clipper<K>::OutPt *
Clipper<K>::addOutPt(TEdge *e, const Point2 &pt)
{
    if (e->outIdx < 0)
    {
        OutRec *outRec = createOutRec();
        outRec->isOpen = (e->windDelta == 0);
        OutPt *newOp = new OutPt;
        outRec->pts = newOp;
        newOp->idx = outRec->idx;
        newOp->pt = pt;
        newOp->next = newOp;
        newOp->prev = newOp;
        if (!outRec->isOpen)
            setHoleState(e, outRec);
        e->outIdx = outRec->idx;
        return newOp;
    }
    else
    {
        OutRec *outRec = _polyOuts[e->outIdx];
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
        op->prev = newOp;
        if (tofront)
            outRec->pts = newOp;
        return newOp;
    }
    return nullptr;
}

template<typename K>
typename Clipper<K>::OutPt*
Clipper<K>::getLastOutPt(TEdge* e)
{
	OutRec* outRec = _polyOuts[e->outIdx];
	if (e->side == EdgeSide::esLeft)
		return outRec->pts;
	else
		return outRec->pts->prev;
}

template<typename K>
typename Clipper<K>::OutPt *
Clipper<K>::addLocalMinPoly(TEdge *e1, TEdge *e2, const Point2 &pt)
{
    OutPt *result = nullptr;
    TEdge *e, *prevE;
    if (isHorizontal(*e2) || (e1->dx > e2->dx))
    {
        result = addOutPt(e1, pt);
        e2->outIdx = e1->outIdx;
        e1->side = EdgeSide::esRight;
        e2->side = EdgeSide::esLeft;
        e = e1;
        if (e->prevInAEL == e2)
            prevE = e2->prevInAEL;
        else
            prevE = e->prevInAEL;
    }
    else
    {
        result = addOutPt(e2, pt);
        e1->outIdx = e2->outIdx;
        e1->side = EdgeSide::esLeft;
        e2->side = EdgeSide::esRight;
        e = e2;
        if (e->prevInAEL == e1)
            prevE = e1->prevInAEL;
        else
            prevE = e->prevInAEL;
    }

    if (prevE && prevE->outIdx >= 0 && prevE->top().y() < pt.y() && e->top().y() < pt.y())
    {
        CT xPrev = topX(*prevE, pt.y());
        CT xE = topX(*e, pt.y());
        if (xPrev == xE && (e->windDelta != 0) && (prevE->windDelta != 0) &&
            slopesEqual(Point2(xPrev, pt.y()), prevE->top(), Point2(xE, pt.y()), e->top()))
        {
            OutPt *outPt = addOutPt(prevE, pt);
            addJoin(result, outPt, e->top());
        }
    }
    return result;
}

template<typename K>
typename Clipper<K>::OutPt* 
Clipper<K>::dupOutPt(OutPt* outPt, bool InsertAfter)
{
	OutPt* result = new OutPt;
	result->pt = outPt->pt;
	result->idx = outPt->idx;
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
Clipper<K>::getBottomPt(Clipper<K>::OutPt* pp)
{
	OutPt* dups = 0;
	OutPt* p = pp->next;
	while (p != pp)
	{
		if (p->pt.y() > pp->pt.y())
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
				if (p->next != pp && p->prev != pp) dups = p;
			}
		}
		p = p->next;
	}
	if (dups)
	{
		//there appears to be at least 2 vertices at BottomPt so ...
		while (dups != p)
		{
			if (!firstIsBottomPt(p, dups))
				pp = dups;
			dups = dups->next;
			while (dups->pt != pp->pt) 
				dups = dups->next;
		}
	}
	return pp;
}

template<typename K>
void Clipper<K>::addLocalMaxPoly(TEdge *e1, TEdge *e2, const Point2 &pt)
{
    addOutPt(e1, pt);
    if (e2->windDelta == 0)
        addOutPt(e2, pt);
    if (e1->outIdx == e2->outIdx)
    {
        e1->outIdx = _unassigned;
        e2->outIdx = _unassigned;
    }
    else if (e1->outIdx < e2->outIdx)
        appendPolygon(e1, e2);
    else
        appendPolygon(e2, e1);
}

template<typename K>
typename Clipper<K>::OutRec *
Clipper<K>::getOutRec(int idx)
{
    OutRec* outrec = _polyOuts[idx];
    while (outrec != _polyOuts[outrec->idx])
    	outrec = _polyOuts[outrec->idx];
    return outrec;
}

template<typename K>
typename Clipper<K>::OutRec *
Clipper<K>::getLowermostRec(OutRec* outRec1, OutRec* outRec2)
{
	//work out which polygon fragment has the correct hole state ...
	if (!outRec1->botPt)
		outRec1->botPt = getBottomPt(outRec1->pts);
	if (!outRec2->botPt)
		outRec2->botPt = getBottomPt(outRec2->pts);
	OutPt* outPt1 = outRec1->botPt;
	OutPt* outPt2 = outRec2->botPt;
	if (outPt1->pt.y() > outPt2->pt.y())
		return outRec1;
	else if (outPt1->pt.y() < outPt2->pt.y()) 
		return outRec2;
	else if (outPt1->pt.x() < outPt2->pt.x())
		return outRec1;
	else if (outPt1->pt.x() > outPt2->pt.x())
		return outRec2;
	else if (outPt1->next == outPt1) 
		return outRec2;
	else if (outPt2->next == outPt2) 
		return outRec1;
	else if (firstIsBottomPt(outPt1, outPt2)) 
		return outRec1;
	else 
		return outRec2;
}

template<typename K>
void Clipper<K>::appendPolygon(TEdge *e1, TEdge *e2)
{
    //get the start and ends of both output polygons ...
    OutRec* outRec1 = _polyOuts[e1->outIdx];
    OutRec* outRec2 = _polyOuts[e2->outIdx];

    OutRec* holeStateRec;
    if (outRec1RightOfOutRec2(outRec1, outRec2))
    	holeStateRec = outRec2;
    else if (outRec1RightOfOutRec2(outRec2, outRec1))
    	holeStateRec = outRec1;
    else
    	holeStateRec = getLowermostRec(outRec1, outRec2);

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
    		reversePolyPtLinks(p2_lft);
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
    		reversePolyPtLinks(p2_lft);
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

    outRec1->botPt = 0;
    if (holeStateRec == outRec2)
    {
    	if (outRec2->firstLeft != outRec1)
    		outRec1->firstLeft = outRec2->firstLeft;
    	outRec1->isHole = outRec2->isHole;
    }
    outRec2->pts = 0;
    outRec2->botPt = 0;
    outRec2->firstLeft = outRec1;

    int OKIdx = e1->outIdx;
    int ObsoleteIdx = e2->outIdx;

    e1->outIdx = _unassigned; //nb: safe because we only get here via addLocalMaxPoly
    e2->outIdx = _unassigned;

    TEdge* e = _activeEdges;
    while (e)
    {
    	if (e->outIdx == ObsoleteIdx)
    	{
    		e->outIdx = OKIdx;
    		e->side = e1->side;
    		break;
    	}
    	e = e->nextInAEL;
    }

    outRec2->idx = outRec1->idx;
}

template<typename K>
void Clipper<K>::intersectEdges(TEdge *e1, TEdge *e2, Point2 &pt)
{
    bool e1Contributing = (e1->outIdx >= 0);
    bool e2Contributing = (e2->outIdx >= 0);

#ifdef use_xyz
    SetZ(pt, *e1, *e2);
#endif

#ifdef USE_LINES
    //if either edge is on an OPEN path ...
    if (e1->windDelta == 0 || e2->windDelta == 0)
    {
        //ignore subject-subject open path intersections UNLESS they
        //are both open paths, AND they are both 'contributing maximas' ...
        if (e1->windDelta == 0 && e2->windDelta == 0)
            return;

        //if intersecting a subj line with a subj poly ...
        else if (e1->polyType == e2->polyType &&
            e1->windDelta != e2->windDelta && _clipType == ClipType::ctUnion)
        {
            if (e1->windDelta == 0)
            {
                if (e2Contributing)
                {
                    addOutPt(e1, pt);
                    if (e1Contributing)
                        e1->outIdx = _unassigned;
                }
            }
            else
            {
                if (e1Contributing)
                {
                    addOutPt(e2, pt);
                    if (e2Contributing)
                        e2->outIdx = _unassigned;
                }
            }
        }
        else if (e1->polyType != e2->polyType)
        {
            //toggle subj open path outIdx on/off when std::abs(clip.WndCnt) == 1 ...
            if ((e1->windDelta == 0) && abs(e2->windCnt) == 1 &&
                (_clipType != ClipType::ctUnion || e2->windCnt2 == 0))
            {
                addOutPt(e1, pt);
                if (e1Contributing)
                    e1->outIdx = _unassigned;
            }
            else if ((e2->windDelta == 0) && (abs(e1->windCnt) == 1) &&
                (_clipType != ClipType::ctUnion || e1->windCnt2 == 0))
            {
                addOutPt(e2, pt);
                if (e2Contributing)
                    e2->outIdx = _unassigned;
            }
        }
        return;
    }
#endif

    //update winding counts...
    //assumes that e1 will be to the right of e2 ABOVE the intersection
    if (e1->polyType == e2->polyType) {
        if (isEvenOddFillType(*e1))
        {
            int oldE1WindCnt = e1->windCnt;
            e1->windCnt = e2->windCnt;
            e2->windCnt = oldE1WindCnt;
        }
        else
        {
            if (e1->windCnt + e2->windDelta == 0)
                e1->windCnt = -e1->windCnt;
            else
                e1->windCnt += e2->windDelta;

            if (e2->windCnt - e1->windDelta == 0)
                e2->windCnt = -e2->windCnt;
            else
                e2->windCnt -= e1->windDelta;
        }
    }
    else {
        if (isEvenOddFillType(*e2))
            e1->windCnt2 = (e1->windCnt2 == 0) ? 1 : 0;
        else
            e1->windCnt2 += e2->windDelta;

        if (isEvenOddFillType(*e1))
            e2->windCnt2 = (e2->windCnt2 == 0) ? 1 : 0;
        else
            e2->windCnt2 -= e1->windDelta;
    }

    PolyFillType e1FillType, e2FillType,
        e1FillType2, e2FillType2;
    if (e1->polyType == PolyType::ptSubject)
    {
        e1FillType = _subFillType;
        e1FillType2 = _clipFillType;
    }
    else
    {
        e1FillType = _clipFillType;
        e1FillType2 = _subFillType;
    }
    if (e2->polyType == PolyType::ptSubject)
    {
        e2FillType = _subFillType;
        e2FillType2 = _clipFillType;
    }
    else
    {
        e2FillType = _clipFillType;
        e2FillType2 = _subFillType;
    }

    int e1Wc, e2Wc;
    switch (e1FillType)
    {
    case PolyFillType::pftPositive:
        e1Wc = e1->windCnt;
        break;
    case PolyFillType::pftNegative:
        e1Wc = -e1->windCnt;
        break;
    default:
        e1Wc = std::abs(e1->windCnt);
    }
    switch (e2FillType)
    {
    case PolyFillType::pftPositive:
        e2Wc = e2->windCnt; 
		break;
    case PolyFillType::pftNegative:
        e2Wc = -e2->windCnt;
		break;
    default:
        e2Wc = std::abs(e2->windCnt);
    }

    if (e1Contributing && e2Contributing)
    {
        if ((e1Wc != 0 && e1Wc != 1) || (e2Wc != 0 && e2Wc != 1) ||
			(e1->polyType != e2->polyType && _clipType != ClipType::ctXor) ||
			(e1->polyType == e2->polyType && _clipType == ClipType::ctUnion))
        {
            addLocalMaxPoly(e1, e2, pt);
        }
        else
        {
            addOutPt(e1, pt);
            addOutPt(e2, pt);
            swapSides(*e1, *e2);
            swapPolyIndexes(*e1, *e2);
        }
    }
    else if (e1Contributing)
    {
        if (e2Wc == 0 || e2Wc == 1)
        {
            addOutPt(e1, pt);
            swapSides(*e1, *e2);
            swapPolyIndexes(*e1, *e2);
        }
    }
    else if (e2Contributing)
    {
        if (e1Wc == 0 || e1Wc == 1)
        {
            addOutPt(e2, pt);
            swapSides(*e1, *e2);
            swapPolyIndexes(*e1, *e2);
        }
    }
    else if ((e1Wc == 0 || e1Wc == 1) && (e2Wc == 0 || e2Wc == 1)) {
        if (e1->polyType != e2->polyType)
        {
            addLocalMinPoly(e1, e2, pt);
        }
		else if (e1Wc == 1 && e2Wc == 1) {
			CT e1Wc2, e2Wc2;
			switch (e1FillType2)
			{
			case PolyFillType::pftPositive:
				e1Wc2 = e1->windCnt2;
				break;
			case PolyFillType::pftNegative:
				e1Wc2 = -e1->windCnt2;
				break;
			default:
				e1Wc2 = std::abs(e1->windCnt2);
			}
			switch (e2FillType2)
			{
			case PolyFillType::pftPositive:
				e2Wc2 = e2->windCnt2;
				break;
			case PolyFillType::pftNegative:
				e2Wc2 = -e2->windCnt2;
				break;
			default:
				e2Wc2 = std::abs(e2->windCnt2);
			}

			switch (_clipType) {
			case ClipType::ctIntersection:
				if (e1Wc2 > 0 && e2Wc2 > 0)
					addLocalMinPoly(e1, e2, pt);
				break;
			case ClipType::ctUnion:
				if (e1Wc2 <= 0 && e2Wc2 <= 0)
					addLocalMinPoly(e1, e2, pt);
				break;
			case ClipType::ctDifference:
				if (((e1->polyType == PolyType::ptClip) && (e1Wc2 > 0) && (e2Wc2 > 0)) ||
					((e1->polyType == PolyType::ptSubject) && (e1Wc2 <= 0) && (e2Wc2 <= 0)))
					addLocalMinPoly(e1, e2, pt);
				break;
			case ClipType::ctXor:
				addLocalMinPoly(e1, e2, pt);
			}
		} else
            swapSides(*e1, *e2);
    }
}

template<typename K>
bool Clipper<K>::processIntersections(const CT topY)
{
    if (!_activeEdges)
        return true;
    try {
        buildIntersectList(topY);
        size_t ilSize = _intersectList.size();
        if (ilSize == 0)
            return true;
        if (ilSize == 1 || fixupIntersectionOrder())
            processIntersectList();
        else
            return false;
    }
    catch (...)
    {
        _sortedEdges = 0;
        disposeIntersectNodes();
        throw ClipperException("processIntersections error");
    }
    _sortedEdges = 0;
    return true;
}

template<typename K>
void Clipper<K>::buildIntersectList(const CT topY)
{
    if (!_activeEdges)
        return;

    //prepare for sorting ...
    TEdge *e = _activeEdges;
    _sortedEdges = e;
    while (e)
    {
        e->prevInSEL = e->prevInAEL;
        e->nextInSEL = e->nextInAEL;
        e->cur = Point2(topX(*e, topY), e->cur.y());
        e = e->nextInAEL;
    }

    //bubblesort ...
    bool isModified;
    do
    {
        isModified = false;
        e = _sortedEdges;
        while (e->nextInSEL)
        {
            TEdge *eNext = e->nextInSEL;
            if (e->cur.x() > eNext->cur.x())
            {
                Point2 pt = intersectPoint(*e, *eNext);
                if (pt.y() > topY)
                    pt = Point2(topX(*e, topY), topY);
                IntersectNode newNode;
                newNode.edge1 = e;
                newNode.edge2 = eNext;
                newNode.pt = pt;
                _intersectList.push_back(newNode);

                swapPositionInSEL(e, eNext);
                isModified = true;
            }
            else
                e = eNext;
        }
        if (e->prevInSEL)
            e->prevInSEL->nextInSEL = 0;
        else
            break;
    } while (isModified);
    _sortedEdges = 0; //important
}

template<typename K>
void Clipper<K>::processIntersectList()
{
    for (size_t i = 0; i < _intersectList.size(); ++i)
    {
        IntersectNode *iNode = &_intersectList[i];
        {
            intersectEdges(iNode->edge1, iNode->edge2, iNode->pt);
            swapPositionsInAEL(iNode->edge1, iNode->edge2);
        }
    }
    _intersectList.clear();
}

template<typename K>
void Clipper<K>::processEdgesAtTopOfScanbeam(const CT topY)
{
    TEdge *e = _activeEdges;
    while (e)
    {
        //1. process maxima, treating them as if they're 'bent' horizontal edges,
        //   but exclude maxima with horizontal edges. nb: e can't be a horizontal.
        bool isMaximaEdge = isMaxima(e, topY);

        if (isMaximaEdge)
        {
            TEdge *eMaxPair = getMaximaPairEx(e);
            isMaximaEdge = (!eMaxPair || !isHorizontal(*eMaxPair));
        }

        if (isMaximaEdge)
        {
            if (_strictSimple)
                _maxima.push_back(e->top().x());

            TEdge *ePrev = e->prevInAEL;
            doMaxima(e);
            if (!ePrev)
                e = _activeEdges;
            else
                e = ePrev->nextInAEL;
        }
        else
        {
            //2. promote horizontal edges, otherwise update cur.x() and cur.y() ...
            if (isIntermediate(e, topY) && isHorizontal(*e->nextInLML))
            {
                updateEdgeIntoAEL(e);
                if (e->outIdx >= 0)
                    addOutPt(e, e->bot());
                addEdgeToSEL(e);
            }
            else
            {
                e->cur = Point2(topX(*e, topY), topY);
#ifdef use_xyz
                e->cur.Z = topY == e->top().y() ? e->top().Z : (topY == e->bot().y() ? e->bot().Z : 0);
#endif
            }

            //When strictlySimple and 'e' is being touched by another edge, then
            //make sure both edges have a vertex here ...
            if (_strictSimple)
            {
                TEdge *ePrev = e->prevInAEL;
                if ((e->outIdx >= 0) && (e->windDelta != 0) && ePrev && (ePrev->outIdx >= 0) &&
                    (ePrev->cur.x() == e->cur.x()) && (ePrev->windDelta != 0))
                {
                    Point2 pt = e->cur;
#ifdef use_xyz
                    SetZ(pt, *ePrev, *e);
#endif
                    OutPt *op = addOutPt(ePrev, pt);
                    OutPt *op2 = addOutPt(e, pt);
                    addJoin(op, op2, pt); //strictlySimple (type-3) join
                }
            }

            e = e->nextInAEL;
        }
    }

    //3. Process horizontals at the top() of the scanbeam ...
    _maxima.sort();
    processHorizontals();
    _maxima.clear();

    //4. Promote intermediate vertices ...
    e = _activeEdges;
    while (e)
    {
        if (isIntermediate(e, topY))
        {
            OutPt *op = 0;
            if (e->outIdx >= 0)
                op = addOutPt(e, e->top());
            updateEdgeIntoAEL(e);

            //if output polygons share an edge, they'll need joining later ...
            TEdge *ePrev = e->prevInAEL;
            TEdge *eNext = e->nextInAEL;
            if (ePrev && ePrev->cur.x() == e->bot().x() &&
                ePrev->cur.y() == e->bot().y() && op &&
                ePrev->outIdx >= 0 && ePrev->cur.y() > ePrev->top().y() &&
                slopesEqual(e->cur, e->top(), ePrev->cur, ePrev->top()) &&
                (e->windDelta != 0) && (ePrev->windDelta != 0))
            {
                OutPt *op2 = addOutPt(ePrev, e->bot());
                addJoin(op, op2, e->top());
            }
            else if (eNext && eNext->cur.x() == e->bot().x() &&
                eNext->cur.y() == e->bot().y() && op &&
                eNext->outIdx >= 0 && eNext->cur.y() > eNext->top().y() &&
                slopesEqual(e->cur, e->top(), eNext->cur, eNext->top()) &&
                (e->windDelta != 0) && (eNext->windDelta != 0))
            {
                OutPt *op2 = addOutPt(eNext, e->bot());
                addJoin(op, op2, e->top());
            }
        }
        e = e->nextInAEL;
    }
}

template<typename K>
void Clipper<K>::buildResult(Paths &polys)
{
    polys.reserve(_polyOuts.size());
    for (PolyOutList::size_type i = 0; i < _polyOuts.size(); ++i)
    {
        if (!_polyOuts[i]->pts)
            continue;
        OutPt *p = _polyOuts[i]->pts->next;
        int cnt = pointCount(p);
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
void Clipper<K>::buildResult2(PolyTree<K>& polytree)
{
	polytree.clear();
	polytree._allNodes.reserve(_polyOuts.size());
	//add each output polygon/contour to polytree ...
	for (PolyOutList::size_type i = 0; i < _polyOuts.size(); i++)
	{
		OutRec* outRec = _polyOuts[i];
		int cnt = pointCount(outRec->pts);
		if ((outRec->isOpen && cnt < 2) || (!outRec->isOpen && cnt < 3)) continue;
		fixHoleLinkage(*outRec);
		PolyNode<K> * pn = new PolyNode<K>();
		//nb: polytree takes ownership of all the PolyNodes
		polytree._allNodes.push_back(pn);
		outRec->polyNode = pn;
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
	polytree._childs.reserve(_polyOuts.size());
	for (PolyOutList::size_type i = 0; i < _polyOuts.size(); i++)
	{
		OutRec* outRec = _polyOuts[i];
		if (!outRec->polyNode) continue;
		if (outRec->isOpen)
		{
			outRec->polyNode->_isOpen = true;
			polytree.addChild(outRec->polyNode);
		}
		else if (outRec->firstLeft && outRec->firstLeft->polyNode)
			outRec->firstLeft->polyNode->addChild(outRec->polyNode);
		else
			polytree.addChild(outRec->polyNode);
	}
}

template<typename K>
void Clipper<K>::setHoleState(TEdge *e, OutRec *outrec)
{
    TEdge *e2 = e->prevInAEL;
    TEdge *eTmp = 0;
    while (e2)
    {
        if (e2->outIdx >= 0 && e2->windDelta != 0)
        {
            if (!eTmp)
                eTmp = e2;
            else if (eTmp->outIdx == e2->outIdx)
                eTmp = 0;
        }
        e2 = e2->prevInAEL;
    }
    if (!eTmp)
    {
        outrec->firstLeft = 0;
        outrec->isHole = false;
    }
    else
    {
        outrec->firstLeft = _polyOuts[eTmp->outIdx];
        outrec->isHole = !outrec->firstLeft->isHole;
    }
}

template<typename K>
void Clipper<K>::disposeIntersectNodes()
{
    _intersectList.clear();
}

template<typename K>
bool Clipper<K>::fixupIntersectionOrder()
{
    //pre-condition: intersections are sorted Bottom-most first.
    //Now it's crucial that intersections are made only between adjacent edges,
    //so to ensure this the order of intersections may need adjusting ...
    copyAELToSEL();

    std::sort(_intersectList.begin(), _intersectList.end(),
        [](const IntersectNode &node1, const IntersectNode &node2) {
            return node2.pt.y() < node1.pt.y();
        });

    size_t cnt = _intersectList.size();
    for (size_t i = 0; i < cnt; ++i)
    {
        if (!edgesAdjacent(_intersectList[i]))
        {
            size_t j = i + 1;
            while (j < cnt && !edgesAdjacent(_intersectList[j]))
                j++;
            if (j == cnt)
                return false;
            std::swap(_intersectList[i], _intersectList[j]);
        }
        swapPositionInSEL(_intersectList[i].edge1, _intersectList[i].edge2);
    }
    return true;
}

template<typename K>
void Clipper<K>::fixupOutPolygon(OutRec &outrec)
{
    //fixupOutPolygon() - removes duplicate points and simplifies consecutive
    //parallel edges by removing the middle vertex.
    OutPt *lastOK = 0;
    outrec.botPt = 0;
    OutPt *pp = outrec.pts;
    bool preserveCol = _preserveCollinear || _strictSimple;

    for (;;)
    {
        if (pp->prev == pp || pp->prev == pp->next)
        {
            disposeOutPts(pp);
            outrec.pts = 0;
            return;
        }

        //test for duplicate points and collinear edges ...
        if ((pp->pt == pp->next->pt) || (pp->pt == pp->prev->pt) || (slopesEqual(pp->prev->pt, pp->pt, pp->next->pt) &&
                (!preserveCol || !pt2BetweenPt1AndPt3(pp->prev->pt, pp->pt, pp->next->pt))))
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
void Clipper<K>::fixupOutPolyline(OutRec &outrec)
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
        disposeOutPts(pp);
        outrec.pts = 0;
        return;
    }
}

template<typename K>
bool Clipper<K>::outRec1RightOfOutRec2(OutRec* outRec1, OutRec* outRec2)
{
	do
	{
		outRec1 = outRec1->firstLeft;
		if (outRec1 == outRec2) return true;
	} while (outRec1);
	return false;
}

template<typename K>
void Clipper<K>::fixHoleLinkage(OutRec &outrec)
{
    //skip OutRecs that (a) contain outermost polygons or
    //(b) already have the correct owner/child linkage ...
    if (!outrec.firstLeft || (outrec.isHole != outrec.firstLeft->isHole && outrec.firstLeft->pts)) 
		return;

    OutRec *orfl = outrec.firstLeft;
    while (orfl && ((orfl->isHole == outrec.isHole) || !orfl->pts))
        orfl = orfl->firstLeft;
    outrec.firstLeft = orfl;
}

template<typename K>
void Clipper<K>::addJoin(OutPt *op1, OutPt *op2, const Point2 offPt)
{
    Joint j;
    j.outPt1 = op1;
    j.outPt2 = op2;
    j.offPt = offPt;
    _joints.push_back(j);
}

template<typename K>
void Clipper<K>::clearJoints()
{
    _joints.clear();
}

template<typename K>
void Clipper<K>::clearGhostJoints()
{
    _ghostJoints.clear();
}

template<typename K>
void Clipper<K>::addGhostJoin(OutPt *op, const Point2 offPt)
{
    Joint j;
    j.outPt1 = op;
    j.outPt2 = 0;
    j.offPt = offPt;
    _ghostJoints.push_back(j);
}

template<typename K>
bool Clipper<K>::joinHorz(OutPt* op1, OutPt* op1b, OutPt* op2, OutPt* op2b, const Point2 pt, bool discardLeft)
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
		op1b = dupOutPt(op1, !discardLeft);
		if (op1b->pt != pt)
		{
			op1 = op1b;
			op1->pt = pt;
			op1b = dupOutPt(op1, !discardLeft);
		}
	}
	else
	{
		while (op1->next->pt.x() >= pt.x() && op1->next->pt.x() <= op1->pt.x() && op1->next->pt.y() == pt.y())
			op1 = op1->next;
		if (!discardLeft && (op1->pt.x() != pt.x()))
			op1 = op1->next;
		op1b = dupOutPt(op1, discardLeft);
		if (op1b->pt != pt)
		{
			op1 = op1b;
			op1->pt = pt;
			op1b = dupOutPt(op1, discardLeft);
		}
	}

	if (Dir2 == RIGHT_TURN)
	{
		while (op2->next->pt.x() <= pt.x() && op2->next->pt.x() >= op2->pt.x() && op2->next->pt.y() == pt.y())
			op2 = op2->next;
		if (discardLeft && (op2->pt.x() != pt.x())) 
			op2 = op2->next;
		op2b = dupOutPt(op2, !discardLeft);
		if (op2b->pt != pt)
		{
			op2 = op2b;
			op2->pt = pt;
			op2b = dupOutPt(op2, !discardLeft);
		};
	}
	else
	{
		while (op2->next->pt.x() >= pt.x() && op2->next->pt.x() <= op2->pt.x() && op2->next->pt.y() == pt.y())
			op2 = op2->next;
		if (!discardLeft && (op2->pt.x() != pt.x())) 
			op2 = op2->next;
		op2b = dupOutPt(op2, discardLeft);
		if (op2b->pt != pt)
		{
			op2 = op2b;
			op2->pt = pt;
			op2b = dupOutPt(op2, discardLeft);
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
bool Clipper<K>::joinPoints(Joint* j, OutRec* outRec1, OutRec* outRec2)
{
	OutPt* op1 = j->outPt1, * op1b;
	OutPt* op2 = j->outPt2, * op2b;

	//There are 3 kinds of joins for output polygons ...
	//1. Horizontal joins where Joint.outPt1 & Joint.outPt2 are vertices anywhere
	//along (horizontal) collinear edges (& Joint.offPt is on the same horizontal).
	//2. Non-horizontal joins where Joint.outPt1 & Joint.outPt2 are at the same
	//location at the Bottom of the overlapping segment (& Joint.offPt is above).
	//3. StrictSimple joins where edges touch but are not collinear and where
	//Joint.outPt1, Joint.outPt2 & Joint.offPt all share the same point.
	bool isHorizontal = (j->outPt1->pt.y() == j->offPt.y());

	if (isHorizontal && (j->offPt == j->outPt1->pt) && (j->offPt == j->outPt2->pt))
	{
		//Strictly Simple join ...
		if (outRec1 != outRec2)
			return false;
		op1b = j->outPt1->next;
		while (op1b != op1 && (op1b->pt == j->offPt))
			op1b = op1b->next;
		bool reverse1 = (op1b->pt.y() > j->offPt.y());
		op2b = j->outPt2->next;
		while (op2b != op2 && (op2b->pt == j->offPt))
			op2b = op2b->next;
		bool reverse2 = (op2b->pt.y() > j->offPt.y());
		if (reverse1 == reverse2)
			return false;
		if (reverse1)
		{
			op1b = dupOutPt(op1, false);
			op2b = dupOutPt(op2, true);
			op1->prev = op2;
			op2->next = op1;
			op1b->next = op2b;
			op2b->prev = op1b;
			j->outPt1 = op1;
			j->outPt2 = op1b;
			return true;
		}
		else
		{
			op1b = dupOutPt(op1, true);
			op2b = dupOutPt(op2, false);
			op1->next = op2;
			op2->prev = op1;
			op1b->prev = op2b;
			op2b->next = op1b;
			j->outPt1 = op1;
			j->outPt2 = op1b;
			return true;
		}
	}
	else if (isHorizontal)
	{
		//treat horizontal joins differently to non-horizontal joins since with
		//them we're not yet sure where the overlapping is. outPt1.pt & outPt2.pt
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
		if (!getOverlap(op1->pt.x(), op1b->pt.x(), op2->pt.x(), op2b->pt.x(), left, right))
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
		j->outPt1 = op1; j->outPt2 = op2;
		return joinHorz(op1, op1b, op2, op2b, pt, discardLeftSide);
	}
	else
	{
		//nb: For non-horizontal joins ...
		//    1. Jr.outPt1.pt.y() == Jr.outPt2.pt.y()
		//    2. Jr.outPt1.pt > Jr.offPt.y()

		//make sure the polygons are correctly oriented ...
		op1b = op1->next;
		while ((op1b->pt == op1->pt) && (op1b != op1)) 
			op1b = op1b->next;
		bool Reverse1 = ((op1b->pt.y() > op1->pt.y()) || !slopesEqual(op1->pt, op1b->pt, j->offPt));
		if (Reverse1)
		{
			op1b = op1->prev;
			while ((op1b->pt == op1->pt) && (op1b != op1)) 
				op1b = op1b->prev;
			if ((op1b->pt.y() > op1->pt.y()) || !slopesEqual(op1->pt, op1b->pt, j->offPt))
				return false;
		};
		op2b = op2->next;
		while ((op2b->pt == op2->pt) && (op2b != op2))
			op2b = op2b->next;
		bool Reverse2 = ((op2b->pt.y() > op2->pt.y()) || !slopesEqual(op2->pt, op2b->pt, j->offPt));
		if (Reverse2)
		{
			op2b = op2->prev;
			while ((op2b->pt == op2->pt) && (op2b != op2))
				op2b = op2b->prev;
			if ((op2b->pt.y() > op2->pt.y()) || !slopesEqual(op2->pt, op2b->pt, j->offPt)) 
				return false;
		}

		if ((op1b == op1) || (op2b == op2) || (op1b == op2b) || ((outRec1 == outRec2) && (Reverse1 == Reverse2)))
			return false;

		if (Reverse1)
		{
			op1b = dupOutPt(op1, false);
			op2b = dupOutPt(op2, true);
			op1->prev = op2;
			op2->next = op1;
			op1b->next = op2b;
			op2b->prev = op1b;
			j->outPt1 = op1;
			j->outPt2 = op1b;
			return true;
		}
		else
		{
			op1b = dupOutPt(op1, true);
			op2b = dupOutPt(op2, false);
			op1->next = op2;
			op2->prev = op1;
			op1b->prev = op2b;
			op2b->next = op1b;
			j->outPt1 = op1;
			j->outPt2 = op1b;
			return true;
		}
	}
	return false;
}

template<typename K>
void Clipper<K>::joinCommonEdges()
{
	for (JoinList::size_type i = 0; i < _joints.size(); i++)
	{
		Joint* join = &_joints[i];

		OutRec* outRec1 = getOutRec(join->outPt1->idx);
		OutRec* outRec2 = getOutRec(join->outPt2->idx);

		if (!outRec1->pts || !outRec2->pts) 
			continue;
		if (outRec1->isOpen || outRec2->isOpen)
			continue;

		//get the polygon fragment with the correct hole state (firstLeft)
		//before calling joinPoints() ...
		OutRec* holeStateRec;
		if (outRec1 == outRec2) 
			holeStateRec = outRec1;
		else if (outRec1RightOfOutRec2(outRec1, outRec2)) 
			holeStateRec = outRec2;
		else if (outRec1RightOfOutRec2(outRec2, outRec1))
			holeStateRec = outRec1;
		else 
			holeStateRec = getLowermostRec(outRec1, outRec2);

		if (!joinPoints(join, outRec1, outRec2))
			continue;

		if (outRec1 == outRec2)
		{
			//instead of joining two polygons, we've just created a new one by
			//splitting one polygon into two.
			outRec1->pts = join->outPt1;
			outRec1->botPt = 0;
			outRec2 = createOutRec();
			outRec2->pts = join->outPt2;

			//update all OutRec2.pts idx's ...
			updateOutPtIdxs(*outRec2);

			if (poly2ContainsPoly1(outRec2->pts, outRec1->pts))
			{
				//outRec1 contains outRec2 ...
				outRec2->isHole = !outRec1->isHole;
				outRec2->firstLeft = outRec1;

				if (_usingPolyTree)
					fixupFirstLefts2(outRec2, outRec1);

				if ((outRec2->isHole ^ _reverseOutput) == (area(outRec2->pts) > 0))
					reversePolyPtLinks(outRec2->pts);

			}
			else if (poly2ContainsPoly1(outRec1->pts, outRec2->pts))
			{
				//outRec2 contains outRec1 ...
				outRec2->isHole = outRec1->isHole;
				outRec1->isHole = !outRec2->isHole;
				outRec2->firstLeft = outRec1->firstLeft;
				outRec1->firstLeft = outRec2;

				if (_usingPolyTree) 
					fixupFirstLefts2(outRec1, outRec2);

				if ((outRec1->isHole ^ _reverseOutput) == (area(outRec1->pts) > 0))
					reversePolyPtLinks(outRec1->pts);
			}
			else
			{
				//the 2 polygons are completely separate ...
				outRec2->isHole = outRec1->isHole;
				outRec2->firstLeft = outRec1->firstLeft;

				//fixup firstLeft pointers that may need reassigning to OutRec2
				if (_usingPolyTree) 
					fixupFirstLefts1(outRec1, outRec2);
			}

		}
		else
		{
			//joined 2 polygons together ...

			outRec2->pts = 0;
			outRec2->botPt = 0;
			outRec2->idx = outRec1->idx;

			outRec1->isHole = holeStateRec->isHole;
			if (holeStateRec == outRec2)
				outRec1->firstLeft = outRec2->firstLeft;
			outRec2->firstLeft = outRec1;

			if (_usingPolyTree)
				fixupFirstLefts3(outRec2, outRec1);
		}
	}
}

template<typename K>
int Clipper<K>::pointCount(OutPt *pts)
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
void Clipper<K>::updateOutPtIdxs(OutRec& outrec)
{
	OutPt* op = outrec.pts;
	do
	{
		op->idx = outrec.idx;
		op = op->prev;
	} while (op != outrec.pts);
}

template<typename K>
inline void Clipper<K>::doSimplePolygons()
{
	PolyOutList::size_type i = 0;
	while (i < _polyOuts.size())
	{
		OutRec* outrec = _polyOuts[i++];
		OutPt* op = outrec->pts;
		if (!op || outrec->isOpen) 
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
					OutRec* outrec2 = createOutRec();
					outrec2->pts = op2;
					updateOutPtIdxs(*outrec2);
					if (poly2ContainsPoly1(outrec2->pts, outrec->pts))
					{
						//OutRec2 is contained by OutRec1 ...
						outrec2->isHole = !outrec->isHole;
						outrec2->firstLeft = outrec;
						if (_usingPolyTree) 
							fixupFirstLefts2(outrec2, outrec);
					}
					else
						if (poly2ContainsPoly1(outrec->pts, outrec2->pts))
						{
							//OutRec1 is contained by OutRec2 ...
							outrec2->isHole = outrec->isHole;
							outrec->isHole = !outrec2->isHole;
							outrec2->firstLeft = outrec->firstLeft;
							outrec->firstLeft = outrec2;
							if (_usingPolyTree) 
								fixupFirstLefts2(outrec, outrec2);
						}
						else
						{
							//the 2 polygons are separate ...
							outrec2->isHole = outrec->isHole;
							outrec2->firstLeft = outrec->firstLeft;
							if (_usingPolyTree) 
								fixupFirstLefts1(outrec, outrec2);
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
bool Clipper<K>::firstIsBottomPt(const OutPt* btmPt1, const OutPt* btmPt2)
{
	OutPt* p = btmPt1->prev;
	while ((p->pt == btmPt1->pt) && (p != btmPt1)) 
		p = p->prev;
	double dx1p = std::fabs(getDelta(btmPt1->pt, p->pt));
	p = btmPt1->next;
	while ((p->pt == btmPt1->pt) && (p != btmPt1)) 
		p = p->next;
	double dx1n = std::fabs(getDelta(btmPt1->pt, p->pt));

	p = btmPt2->prev;
	while ((p->pt == btmPt2->pt) && (p != btmPt2))
		p = p->prev;
	double dx2p = std::fabs(getDelta(btmPt2->pt, p->pt));
	p = btmPt2->next;
	while ((p->pt == btmPt2->pt) && (p != btmPt2))
		p = p->next;
	double dx2n = std::fabs(getDelta(btmPt2->pt, p->pt));

	if (std::max(dx1p, dx1n) == std::max(dx2p, dx2n) && std::min(dx1p, dx1n) == std::min(dx2p, dx2n))
		return area(btmPt1) > 0; //if otherwise identical use orientation
	else
		return (dx1p >= dx2p && dx1p >= dx2n) || (dx1n >= dx2p && dx1n >= dx2n);
}

template<typename K>
typename Clipper<K>:: OutRec* 
Clipper<K>::parseFirstLeft(OutRec* firstLeft)
{
	while (firstLeft && !firstLeft->pts)
		firstLeft = firstLeft->firstLeft;
	return firstLeft;
}

template<typename K>
bool Clipper<K>::poly2ContainsPoly1(OutPt* outPt1, OutPt* outPt2)
{
	OutPt* op = outPt1;
	do
	{
		//nb: PointInPolygon returns 0 if false, +1 if true, -1 if pt on polygon
		int res = true;// PointInPolygon(op->pt, outPt2);
		if (res >= 0) 
			return res > 0;
		op = op->next;
	} while (op != outPt1);
	return true;
}

template<typename K>
void Clipper<K>::fixupFirstLefts1(OutRec* oldOutRec, OutRec* newOutRec)
{
	//tests if newOutRec contains the polygon before reassigning firstLeft
	for (PolyOutList::size_type i = 0; i < _polyOuts.size(); ++i)
	{
		OutRec* outRec = _polyOuts[i];
		OutRec* firstLeft = parseFirstLeft(outRec->firstLeft);
		if (outRec->pts && firstLeft == oldOutRec)
		{
			if (poly2ContainsPoly1(outRec->pts, newOutRec->pts))
				outRec->firstLeft = newOutRec;
		}
	}
}

template<typename K>
void Clipper<K>::fixupFirstLefts2(OutRec* innerOutRec, OutRec* outerOutRec)
{
	//A polygon has split into two such that one is now the inner of the other.
	//It's possible that these polygons now wrap around other polygons, so check
	//every polygon that's also contained by outerOutRec's firstLeft container
	//(including 0) to see if they've become inner to the new inner polygon ...
	OutRec* orfl = outerOutRec->firstLeft;
	for (PolyOutList::size_type i = 0; i < _polyOuts.size(); ++i)
	{
		OutRec* outRec = _polyOuts[i];

		if (!outRec->pts || outRec == outerOutRec || outRec == innerOutRec)
			continue;
		OutRec* firstLeft = parseFirstLeft(outRec->firstLeft);
		if (firstLeft != orfl && firstLeft != innerOutRec && firstLeft != outerOutRec)
			continue;
		if (poly2ContainsPoly1(outRec->pts, innerOutRec->pts))
			outRec->firstLeft = innerOutRec;
		else if (poly2ContainsPoly1(outRec->pts, outerOutRec->pts))
			outRec->firstLeft = outerOutRec;
		else if (outRec->firstLeft == innerOutRec || outRec->firstLeft == outerOutRec)
			outRec->firstLeft = orfl;
	}
}

template<typename K>
void Clipper<K>::fixupFirstLefts3(OutRec* oldOutRec, OutRec* newOutRec)
{
	//reassigns firstLeft WITHOUT testing if newOutRec contains the polygon
	for (PolyOutList::size_type i = 0; i < _polyOuts.size(); ++i)
	{
		OutRec* outRec = _polyOuts[i];
		OutRec* firstLeft = parseFirstLeft(outRec->firstLeft);
		if (outRec->pts && firstLeft == oldOutRec)
			outRec->firstLeft = newOutRec;
	}
}

template<typename K>
typename Clipper<K>::Point2 
Clipper<K>::intersectPoint(TEdge &edge1, TEdge &edge2)
{
    auto result = intersection(edge1, edge2);
    Point2 ret= boost::get<Point2>(result.get());
    return ret;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

template<typename K>
void ClipperOffset<K>::reversePath(Path& p)
{
	std::reverse(p.begin(), p.end());
}

template<typename K>
void ClipperOffset<K>::reversePaths(Paths& p)
{
	for (Paths::size_type i = 0; i < p.size(); ++i)
		reversePath(p[i]);
}

template<typename K>
bool ClipperOffset<K>::orientation(const Path& p)
{
    return Clipper<K>::area(p) < 0;
}

template<typename K>
void ClipperOffset<K>::simplifyPolygon(const Path& in_poly, Paths& out_polys, PolyFillType fillType)
{
	Clipper<K> c;
	c.setStrictlySimple(true);
	c.addPath(in_poly, PolyType::ptSubject, true);
	c.execute(ClipType::ctUnion, out_polys, fillType, fillType);
}

template<typename K>
void ClipperOffset<K>::simplifyPolygons(const Paths& in_polys, Paths& out_polys, PolyFillType fillType)
{
	Clipper<K> c;
	c.setStrictlySimple(true);
	c.addPaths(in_polys, PolyType::ptSubject, true);
	c.execute(ClipType::ctUnion, out_polys, fillType, fillType);
}

template<typename K>
void ClipperOffset<K>::simplifyPolygons(Paths& polys, PolyFillType fillType)
{
	simplifyPolygons(polys, polys, fillType);
}

template<typename K>
bool ClipperOffset<K>::slopesNearCollinear(const Point2& pt1, const Point2& pt2, const Point2& pt3, double distSqrd)
{
	//this function is more accurate when the point that's geometrically
	//between the other 2 points is the one that's tested for distance.
	//ie makes it more likely to pick up 'spikes' ...
	if (std::abs(pt1.x() - pt2.x()) > std::abs(pt1.y() - pt2.y()))
	{
		if ((pt1.x() > pt2.x()) == (pt1.x() < pt3.x()))
			return distanceFromLineSqrd(pt1, pt2, pt3) < distSqrd;
		else if ((pt2.x() > pt1.x()) == (pt2.x() < pt3.x()))
			return distanceFromLineSqrd(pt2, pt1, pt3) < distSqrd;
		else
			return distanceFromLineSqrd(pt3, pt1, pt2) < distSqrd;
	}
	else
	{
		if ((pt1.y() > pt2.y()) == (pt1.y() < pt3.y()))
			return distanceFromLineSqrd(pt1, pt2, pt3) < distSqrd;
		else if ((pt2.y() > pt1.y()) == (pt2.y() < pt3.y()))
			return distanceFromLineSqrd(pt2, pt1, pt3) < distSqrd;
		else
			return distanceFromLineSqrd(pt3, pt1, pt2) < distSqrd;
	}
}

template<typename K>
bool ClipperOffset<K>::pointsAreClose(Point2 pt1, Point2 pt2, double distSqrd)
{
	double Dx = (double)pt1.x() - pt2.x();
	double dy = (double)pt1.y() - pt2.y();
	return ((Dx * Dx) + (dy * dy) <= distSqrd);
}

template<typename K>
void ClipperOffset<K>::cleanPolygon(const Path& in_poly, Path& out_poly, double distance)
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
		if (pointsAreClose(op->pt, op->prev->pt, distSqrd))
		{
			op = excludeOp(op);
			size--;
		}
		else if (pointsAreClose(op->prev->pt, op->next->pt, distSqrd))
		{
			excludeOp(op->next);
			op = excludeOp(op);
			size -= 2;
		}
		else if (slopesNearCollinear(op->prev->pt, op->pt, op->next->pt, distSqrd))
		{
			op = excludeOp(op);
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
void ClipperOffset<K>::cleanPolygon(Path& poly, double distance)
{
	cleanPolygon(poly, poly, distance);
}

template<typename K>
void ClipperOffset<K>::cleanPolygons(const Paths& in_polys, Paths& out_polys, double distance)
{
	out_polys.resize(in_polys.size());
	for (Paths::size_type i = 0; i < in_polys.size(); ++i)
		cleanPolygon(in_polys[i], out_polys[i], distance);
}

template<typename K>
void ClipperOffset<K>::cleanPolygons(Paths& polys, double distance)
{
	cleanPolygons(polys, polys, distance);
}

template<typename K>
void ClipperOffset<K>::minkowski(const Path& poly, const Path& path, Paths& solution, bool isSum, bool isClosed)
{
	int delta = (isClosed ? 1 : 0);
	size_t polyCnt = poly.size();
	size_t pathCnt = path.size();
	Paths pp;
	pp.reserve(pathCnt);
	if (isSum)
		for (size_t i = 0; i < pathCnt; ++i)
		{
			Path p;
			p.reserve(polyCnt);
			for (size_t j = 0; j < poly.size(); ++j)
				p.push_back(Point2(path[i].x() + poly[j].x(), path[i].y() + poly[j].y()));
			pp.push_back(p);
		}
	else
		for (size_t i = 0; i < pathCnt; ++i)
		{
			Path p;
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
			Path quad;
			quad.reserve(4);
			quad.push_back(pp[i % pathCnt][j % polyCnt]);
			quad.push_back(pp[(i + 1) % pathCnt][j % polyCnt]);
			quad.push_back(pp[(i + 1) % pathCnt][(j + 1) % polyCnt]);
			quad.push_back(pp[i % pathCnt][(j + 1) % polyCnt]);
			if (!orientation(quad)) 
				reversePath(quad);
			solution.push_back(quad);
		}
}

template<typename K>
void ClipperOffset<K>::minkowskiSum(const Path& pattern, const Path& path, Paths& solution, bool pathIsClosed)
{
	minkowski(pattern, path, solution, true, pathIsClosed);
	Clipper<K> c;
	c.addPaths(solution, PolyType::ptSubject, true);
	c.execute(ClipType::ctUnion, solution, PolyFillType::pftNonZero, PolyFillType::pftNonZero);
}

template<typename K>
void ClipperOffset<K>::translatePath(const Path& input, Path& output, const Point2 delta)
{
	//precondition: input != output
	output.resize(input.size());
	for (size_t i = 0; i < input.size(); ++i)
		output[i] = Point2(input[i].x() + delta.x(), input[i].y() + delta.y());
}

template<typename K>
void ClipperOffset<K>::minkowskiSum(const Path& pattern, const Paths& paths, Paths& solution, bool pathIsClosed)
{
	Clipper<K> c;
	for (size_t i = 0; i < paths.size(); ++i)
	{
		Paths tmp;
		minkowski(pattern, paths[i], tmp, true, pathIsClosed);
		c.addPaths(tmp, PolyType::ptSubject, true);
		if (pathIsClosed)
		{
			Path tmp2;
			translatePath(paths[i], tmp2, pattern[0]);
			c.addPath(tmp2, PolyType::ptClip, true);
		}
	}
	c.execute(ClipType::ctUnion, solution, PolyFillType::pftNonZero, PolyFillType::pftNonZero);
}

template<typename K>
void ClipperOffset<K>::minkowskiDiff(const Path& poly1, const Path& poly2, Paths& solution)
{
	minkowski(poly1, poly2, solution, false, true);
	Clipper<K> c;
	c.addPaths(solution, PolyType::ptSubject, true);
	c.execute(ClipType::ctUnion, solution, PolyFillType::pftNonZero, PolyFillType::pftNonZero);
}

template<typename K>
void ClipperOffset<K>::addPolyNodeToPaths(const PolyNode& polynode, NodeType nodetype, Paths& paths)
{
	bool match = true;
	if (nodetype == ntClosed) match = !polynode.isOpen();
	else if (nodetype == ntOpen) return;

	if (!polynode._contour.empty() && match)
		paths.push_back(polynode._contour);
	for (int i = 0; i < polynode.childCount(); ++i)
		addPolyNodeToPaths(*polynode._childs[i], nodetype, paths);
}

template<typename K>
void ClipperOffset<K>::polyTreeToPaths(const PolyTree& polytree, Paths& paths)
{
	paths.resize(0);
	paths.reserve(polytree.total());
	addPolyNodeToPaths(polytree, ntAny, paths);
}

template<typename K>
void ClipperOffset<K>::closedPathsFromPolyTree(const PolyTree& polytree, Paths& paths)
{
	paths.resize(0);
	paths.reserve(polytree.total());
	addPolyNodeToPaths(polytree, ntClosed, paths);
}

template<typename K>
void ClipperOffset<K>::openPathsFromPolyTree(PolyTree& polytree, Paths& paths)
{
	paths.resize(0);
	paths.reserve(polytree.total());
	//Open paths are top level only, so ...
	for (int i = 0; i < polytree.childCount(); ++i)
		if (polytree._childs[i]->isOpen())
			paths.push_back(polytree._childs[i]->_contour);
}

template<typename K>
double ClipperOffset<K>::distanceSqrt(const Point2& pt1, const Point2& pt2)
{
	double Dx = ((double)pt1.x() - pt2.x());
	double dy = ((double)pt1.y() - pt2.y());
	return (Dx * Dx + dy * dy);
}

template<typename K>
double ClipperOffset<K>::distanceFromLineSqrd( const Point2& pt, const Point2& ln1, const Point2& ln2)
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
ClipperOffset<K>::excludeOp(typename Clipper<K>::OutPt* op)
{
	Clipper<K>::OutPt* result = op->prev;
	result->next = op->next;
	op->next->prev = result;
	result->idx = 0;
	return result;
}

template<typename K>
typename ClipperOffset<K>:: DoublePoint 
ClipperOffset<K>::getUnitNormal(const Point2& pt1, const Point2& pt2)
{
	if (pt2.x() == pt1.x() && pt2.y() == pt1.y())
		return DoublePoint(0, 0);

	double Dx = (double)(pt2.x() - pt1.x());
	double dy = (double)(pt2.y() - pt1.y());
	double f = 1 * 1.0 / std::sqrt(Dx * Dx + dy * dy);
	Dx *= f;
	dy *= f;
	return DoublePoint(dy, -Dx);
}

//---------------------------------------------------------------------------

template<typename K>
ClipperOffset<K>::ClipperOffset(double miterLimit, double roundPrecision)
{
	_miterLimit = miterLimit;
	_arcTolerance = roundPrecision;
	_lowest = Point2(-1, _lowest.y());
}

template<typename K>
ClipperOffset<K>::~ClipperOffset()
{
	clear();
}

template<typename K>
void ClipperOffset<K>::addPath(const Path& path, JointType joinType, EndType endType) 
{
	int highI = (int)path.size() - 1;
	if (highI < 0)
		return;
	PolyNode* newNode = new PolyNode();
	newNode->_joinType = joinType;
	newNode->_endType = endType;

	//strip duplicate points from path and also get index to the lowest point ...
	if (endType == EndType::etClosedLine || endType == EndType::etClosedPolygon)
	{
		while (highI > 0 && path[0] == path[highI])
			highI--;
	}
	newNode->_contour.reserve(highI + 1);
	newNode->_contour.push_back(path[0]);
	int j = 0, k = 0;
	for (int i = 1; i <= highI; i++)
	{
		if (newNode->_contour[j] != path[i])
		{
			j++;
			newNode->_contour.push_back(path[i]);
			if (path[i].y() > newNode->_contour[k].y() ||
				(path[i].y() == newNode->_contour[k].y() &&
					path[i].x() < newNode->_contour[k].x())) k = j;
		}
	}
	if (endType == EndType::etClosedPolygon && j < 2)
	{
		delete newNode;
		return;
	}
	_polyNodes.addChild(*newNode);

	//if this path's lowest pt is lower than all the others then update _lowest
	if (endType != EndType::etClosedPolygon)
		return;
	if (_lowest.x() < 0)
		_lowest = Point2(_polyNodes.childCount() - 1, k);
	else
	{
		Point2 ip = _polyNodes._childs[(int)_lowest.x()]->_contour[(int)_lowest.y()];
		if (newNode->_contour[k].y() > ip.y() ||
			(newNode->_contour[k].y() == ip.y() &&
				newNode->_contour[k].x() < ip.x()))
			_lowest = Point2(_polyNodes.childCount() - 1, k);
	}

}

template<typename K>
void ClipperOffset<K>::addPaths(const Paths& paths, JointType joinType, EndType endType)
{
	for (Paths::size_type i = 0; i < paths.size(); ++i)
		addPath(paths[i], joinType, endType);
}

template<typename K>
void ClipperOffset<K>::execute(Paths& solution, double delta)
{
	solution.clear();
	fixOrientations();
	doOffset(delta);

	//now clean up 'corners' ...
	Clipper<K> clpr;
	clpr.addPaths(_destPolys, PolyType::ptSubject, true);
	if (delta > 0)
	{
		clpr.execute(ClipType::ctUnion, solution, PolyFillType::pftPositive, PolyFillType::pftPositive);
	}
	else
	{
		Rectangle2<K> r = clpr.getBounds();
		Path outer(4);
		outer[0] = Point2(r.left - 10, r.bottom + 10);
		outer[1] = Point2(r.right + 10, r.bottom + 10);
		outer[2] = Point2(r.right + 10, r.top - 10);
		outer[3] = Point2(r.left - 10, r.top - 10);

		clpr.addPath(outer, PolyType::ptSubject, true);
		clpr.setReverseSolution(true);
		clpr.execute(ClipType::ctUnion, solution, PolyFillType::pftNegative, PolyFillType::pftNegative);
		if (solution.size() > 0)
			solution.erase(solution.begin());
	}
}

template<typename K>
void ClipperOffset<K>::execute(PolyTree& solution, double delta)
{
	solution.clear();
	fixOrientations();
	doOffset(delta);

	//now clean up 'corners' ...
	Clipper<K> clpr;
	clpr.addPaths(_destPolys, PolyType::ptSubject, true);
	if (delta > 0)
	{
		clpr.execute(ClipType::ctUnion, solution, PolyFillType::pftPositive, PolyFillType::pftPositive);
	}
	else
	{
		Rectangle2 r = clpr.getBounds();
		Path outer(4);
		outer[0] = Point2(r.left - 10, r.bottom + 10);
		outer[1] = Point2(r.right + 10, r.bottom + 10);
		outer[2] = Point2(r.right + 10, r.top - 10);
		outer[3] = Point2(r.left - 10, r.top - 10);

		clpr.addPath(outer, PolyType::ptSubject, true);
		clpr.setReverseSolution(true);
		clpr.execute(ClipType::ctUnion, solution, PolyFillType::pftNegative, PolyFillType::pftNegative);
		//remove the outer PolyNode rectangle ...
		if (solution.childCount() == 1 && solution._childs[0]->childCount() > 0)
		{
			PolyNode* outerNode = solution._childs[0];
			solution._childs.reserve(outerNode->childCount());
			solution._childs[0] = outerNode->_childs[0];
			solution._childs[0]->_parent = outerNode->_parent;
			for (int i = 1; i < outerNode->childCount(); ++i)
				solution.addChild(*outerNode->_childs[i]);
		}
		else
			solution.clear();
	}
}

template<typename K>
void ClipperOffset<K>::clear()
{
	for (int i = 0; i < _polyNodes.childCount(); ++i)
		delete _polyNodes._childs[i];
	_polyNodes._childs.clear();
	_lowest = Point2(-1, _lowest.y());
}

template<typename K>
void ClipperOffset<K>::fixOrientations()
{
	//fixup orientations of all closed paths if the orientation of the
	//closed path with the lowermost vertex is wrong ...
	if (_lowest.x() >= 0 && !orientation(_polyNodes._childs[(int)_lowest.x()]->_contour))
	{
		for (int i = 0; i < _polyNodes.childCount(); ++i)
		{
			PolyNode& node = *_polyNodes._childs[i];
			if (node._endType == EndType::etClosedPolygon ||
				(node._endType == EndType::etClosedLine && orientation(node._contour)))
				reversePath(node._contour);
		}
	}
	else
	{
		for (int i = 0; i < _polyNodes.childCount(); ++i)
		{
			PolyNode& node = *_polyNodes._childs[i];
			if (node._endType == EndType::etClosedLine && !orientation(node._contour))
				reversePath(node._contour);
		}
	}
}

template<typename K>
void ClipperOffset<K>::doOffset(double delta)
{
	_destPolys.clear();
	_delta = delta;

	//if Zero offset, just copy any CLOSED polygons to m_p and return ...
	if (delta < 0 || delta > 0)
	{
		_destPolys.reserve(_polyNodes.childCount());
		for (int i = 0; i < _polyNodes.childCount(); i++)
		{
			PolyNode *node = _polyNodes._childs[i];
			if (node->_endType == EndType::etClosedPolygon)
				_destPolys.push_back(node->_contour);
		}
		return;
	}

	//see offset_triginometry3.svg in the documentation folder ...
	if (_miterLimit > 2) _miterLim = 2 / (_miterLimit * _miterLimit);
	else _miterLim = 0.5;

	double y;
	if (_arcTolerance <= 0.0) y = arcTolerance;
	else if (_arcTolerance > std::fabs(delta) * arcTolerance)
		y = std::fabs(delta) * arcTolerance;
	else y = _arcTolerance;
	//see offset_triginometry2.svg in the documentation folder ...
	double steps = pi / std::acos(1 - y / std::fabs(delta));
	if (steps > std::fabs(delta) * pi)
		steps = std::fabs(delta) * pi;  //ie excessive precision check
	_sin = std::sin(pi2 / steps);
	_cos = std::cos(pi2 / steps);
	_stepsPerRad = steps / pi2;
	if (delta < 0.0) _sin = -_sin;

	_destPolys.reserve(_polyNodes.childCount() * 2);
	for (int i = 0; i < _polyNodes.childCount(); i++)
	{
		PolyNode& node = *_polyNodes._childs[i];
		_srcPoly = node._contour;

		int len = (int)_srcPoly.size();
		if (len == 0 || (delta <= 0 && (len < 3 || node._endType != EndType::etClosedPolygon)))
			continue;

		_destPoly.clear();
		if (len == 1)
		{
			if (node._joinType == JointType::jtRound)
			{
				double x = 1.0, y = 0.0;
				for (int j = 1; j <= steps; j++)
				{
					_destPoly.push_back(Point2( _srcPoly[0].x() + x * delta, _srcPoly[0].y() + y * delta));
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
					_destPoly.push_back(Point2( _srcPoly[0].x() + x * delta, _srcPoly[0].y() + y * delta));
					if (x < 0)
						x = 1;
					else if (y < 0) 
						y = 1;
					else 
						x = -1;
				}
			}
			_destPolys.push_back(_destPoly);
			continue;
		}
		//build _normals ...
		_normals.clear();
		_normals.reserve(len);
		for (int j = 0; j < len - 1; ++j)
			_normals.push_back(getUnitNormal(_srcPoly[j], _srcPoly[j + 1]));
		if (node._endType == EndType::etClosedLine || node._endType == EndType::etClosedPolygon)
			_normals.push_back(getUnitNormal(_srcPoly[len - 1], _srcPoly[0]));
		else
			_normals.push_back(DoublePoint(_normals[len - 2]));

		if (node._endType == EndType::etClosedPolygon)
		{
			int k = len - 1;
			for (int j = 0; j < len; ++j)
				offsetPoint(j, k, node._joinType);
			_destPolys.push_back(_destPoly);
		}
		else if (node._endType == EndType::etClosedLine)
		{
			int k = len - 1;
			for (int j = 0; j < len; ++j)
				offsetPoint(j, k, node._joinType);
			_destPolys.push_back(_destPoly);
			_destPoly.clear();
			//re-build _normals ...
			DoublePoint n = _normals[len - 1];
			for (int j = len - 1; j > 0; j--)
				_normals[j] = DoublePoint(-_normals[j - 1].x(), -_normals[j - 1].y());
			_normals[0] = DoublePoint(-n.x(), -n.y());
			k = 0;
			for (int j = len - 1; j >= 0; j--)
				offsetPoint(j, k, node._joinType);
			_destPolys.push_back(_destPoly);
		}
		else
		{
			int k = 0;
			for (int j = 1; j < len - 1; ++j)
				offsetPoint(j, k, node._joinType);

			Point2 pt1;
			if (node._endType == EndType::etOpenButt)
			{
				int j = len - 1;
				pt1 = Point2((_srcPoly[j].x() + _normals[j].x() * delta), (_srcPoly[j].y() + _normals[j].y() * delta));
				_destPoly.push_back(pt1);
				pt1 = Point2((_srcPoly[j].x() - _normals[j].x() * delta), (_srcPoly[j].y() - _normals[j].y() * delta));
				_destPoly.push_back(pt1);
			}
			else
			{
				int j = len - 1;
				k = len - 2;
				_sinA = 0;
				_normals[j] = DoublePoint(-_normals[j].x(), -_normals[j].y());
				if (node._endType == EndType::etOpenSquare)
					doSquare(j, k);
				else
					doRound(j, k);
			}

			//re-build _normals ...
			for (int j = len - 1; j > 0; j--)
				_normals[j] = DoublePoint(-_normals[j - 1].x(), -_normals[j - 1].y());
			_normals[0] = DoublePoint(-_normals[1].x(), -_normals[1].y());

			k = len - 1;
			for (int j = k - 1; j > 0; --j) offsetPoint(j, k, node._joinType);

			if (node._endType == EndType::etOpenButt)
			{
				pt1 = Point2((_srcPoly[0].x() - _normals[0].x() * delta),
					(_srcPoly[0].y() - _normals[0].y() * delta));
				_destPoly.push_back(pt1);
				pt1 = Point2((_srcPoly[0].x() + _normals[0].x() * delta),
					(_srcPoly[0].y() + _normals[0].y() * delta));
				_destPoly.push_back(pt1);
			}
			else
			{
				k = 1;
				_sinA = 0;
				if (node._endType == EndType::etOpenSquare)
					doSquare(0, 1);
				else
					doRound(0, 1);
			}
			_destPolys.push_back(_destPoly);
		}
	}
}

template<typename K>
void ClipperOffset<K>::offsetPoint(int j, int& k, JointType jointype)
{
	//cross product ...
	_sinA = (_normals[k].x() * _normals[j].y() - _normals[j].x() * _normals[k].y());
	if (std::fabs(_sinA * _delta) < 1.0)
	{
		//dot product ...
		double cosA = (_normals[k].x() * _normals[j].x() + _normals[j].y() * _normals[k].y());
		if (cosA > 0) // angle => 0 degrees
		{
			_destPoly.push_back(Point2((_srcPoly[j].x() + _normals[k].x() * _delta),
				(_srcPoly[j].y() + _normals[k].y() * _delta)));
			return;
		}
		//else angle => 180 degrees   
	}
	else if (_sinA > 1.0) _sinA = 1.0;
	else if (_sinA < -1.0) _sinA = -1.0;

	if (_sinA * _delta < 0)
	{
		_destPoly.push_back(Point2((_srcPoly[j].x() + _normals[k].x() * _delta),
			(_srcPoly[j].y() + _normals[k].y() * _delta)));
		_destPoly.push_back(_srcPoly[j]);
		_destPoly.push_back(Point2((_srcPoly[j].x() + _normals[j].x() * _delta),
			(_srcPoly[j].y() + _normals[j].y() * _delta)));
	}
	else
		switch (jointype)
		{
		case JointType::jtMiter:
		{
			double r = 1 + (_normals[j].x() * _normals[k].x() + _normals[j].y() * _normals[k].y());
			if (r >= _miterLim)
				doMiter(j, k, r);
			else
				doSquare(j, k);
			break;
		}
		case JointType::jtSquare:
			doSquare(j, k);
			break;
		case JointType::jtRound:
			doRound(j, k);
			break;
		}
	k = j;
}

template<typename K>
void ClipperOffset<K>::doSquare(int j, int k)
{
	double dx = std::tan(std::atan2(_sinA,
		_normals[k].x() * _normals[j].x() + _normals[k].y() * _normals[j].y()) / 4);
	_destPoly.push_back(Point2(
		_srcPoly[j].x() + _delta * (_normals[k].x() - _normals[k].y() * dx),
		_srcPoly[j].y() + _delta * (_normals[k].y() + _normals[k].x() * dx)));
	_destPoly.push_back(Point2(
		_srcPoly[j].x() + _delta * (_normals[j].x() + _normals[j].y() * dx),
		_srcPoly[j].y() + _delta * (_normals[j].y() - _normals[j].x() * dx)));
}

template<typename K>
void ClipperOffset<K>::doMiter(int j, int k, double r)
{
	double q = _delta / r;
	_destPoly.push_back(Point2(_srcPoly[j].x() + (_normals[k].x() + _normals[j].x()) * q,
		_srcPoly[j].y() + (_normals[k].y() + _normals[j].y()) * q));
}

template<typename K>
void ClipperOffset<K>::doRound(int j, int k)
{
	double a = std::atan2(_sinA, _normals[k].x() * _normals[j].x() + _normals[k].y() * _normals[j].y());
	int steps = std::max<double>(_stepsPerRad * std::fabs(a), 1);

	double x = _normals[k].x(), y = _normals[k].y(), X2;
	for (int i = 0; i < steps; ++i)
	{
		_destPoly.push_back(Point2(_srcPoly[j].x() + x * _delta, _srcPoly[j].y() + y * _delta));
		X2 = x;
		x = x * _cos - _sin * y;
		y = X2 * _sin + y * _cos;
	}
	_destPoly.push_back(Point2(_srcPoly[j].x() + _normals[j].x() * _delta, _srcPoly[j].y() + _normals[j].y() * _delta));
}

template<typename T>
void test()
{
	CGAL::Polygon2Clipper::Clipper<T> cl;
	typedef Clipper<T>::Point2 Point2;

	int ex = 7;

	if (ex == 1) {
		Clipper<T>::Path poly, poly1;
		poly.push_back(Point2(0, 0));
		poly.push_back(Point2(10, 10));
		poly.push_back(Point2(0, 20));
		poly.push_back(Point2(0, 0));

		poly1.push_back(Point2(15, 0));
		poly1.push_back(Point2(15, 20));
		poly1.push_back(Point2(5, 10));
		poly1.push_back(Point2(15, 0));

		cl.addPath(poly, PolyType::ptSubject);
		cl.addPath(poly1, PolyType::ptClip);

		Clipper<K>::Paths out;
		cl.execute(ClipType::ctIntersection, out);
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

		cl.addPath(poly, PolyType::ptSubject);
		cl.addPath(poly1, PolyType::ptClip);

		Clipper<K>::Paths out;
		cl.execute(ClipType::ctIntersection, out);
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

		cl.addPath(poly, PolyType::ptSubject);
		cl.addPath(poly1, PolyType::ptClip);

		Clipper<K>::Paths out;
		cl.execute(ClipType::ctIntersection, out);
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

		cl.addPath(poly, PolyType::ptSubject);
		cl.addPath(poly1, PolyType::ptClip);

		Clipper<K>::Paths out;
		cl.execute(ClipType::ctIntersection, out);
	}
	else if (ex == 5)
	{
		Clipper<K>::Path poly;
		poly.push_back(Point2(0, 10));
		poly.push_back(Point2(0, -10));
		poly.push_back(Point2(10, 0));
		poly.push_back(Point2(-10, 0));
		poly.push_back(Point2(0, 10));
		cl.addPath(poly, PolyType::ptSubject, true);

		//cl.setStrictlySimple(true);
		PolyTree<K> out;
		cl.execute(ClipType::ctUnion, out, PolyFillType::pftPositive);
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
		cl.addPath(poly, PolyType::ptSubject, false);

		//cl.setStrictlySimple(true);
		PolyTree<K> out;
		cl.execute(ClipType::ctUnion, out, PolyFillType::pftPositive);
		printf("");
	} 
	else if (ex == 7)
	{
		Clipper<K>::Path poly;
		poly.push_back(Point2(0, 10));
		poly.push_back(Point2(0, -10));
		poly.push_back(Point2(10, 0));
		poly.push_back(Point2(-10, 5));
		cl.addPath(poly, PolyType::ptSubject, false);

		cl.setStrictlySimple(true);
		PolyTree<K> out;
		cl.execute(ClipType::ctUnion, out, PolyFillType::pftNegative);
		printf("");
	}
	else if (ex == 8)
	{
		Clipper<K>::Path poly;
		poly.push_back(Point2(0, 10));
		poly.push_back(Point2(0, -10));
		poly.push_back(Point2(10, 0));
		cl.addPath(poly, PolyType::ptSubject, false);

		//cl.setStrictlySimple(true);
		PolyTree<K> out;
		cl.execute(ClipType::ctUnion, out, PolyFillType::pftPositive);
		printf("");
	}
}

}}

#endif
