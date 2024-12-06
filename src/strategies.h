#ifndef STRATEGIES_H
#define STRATEGIES_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Polygon_2.h>

// Kernel and type definitions
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds, Itag> CDT;
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CDT::Face_handle Face_handle;

// Function declarations
int is_obtuse_angle(Point &A, Point &B, Point &C);
int count_Obtuse_Angles(CDT &cdt);
Point project_point(Point &A, Point &B, Point &P);
Polygon_2 find_convex_polygon_around_obtuse_triangle(CDT &cdt, Face_handle face);
Point insert_Steiner_point_in_convex_polygons(CDT &cdt, Polygon_2 &region_boundary);
Point select_steiner_point(Point &a, Point &b, Point &c, int strategy, CDT &cdt, Polygon_2 convex_hull);

#endif // STRATEGIES_H