#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/mark_domain_in_triangulation.h>
#include <CGAL/Delaunay_mesher_2.h>
#include "triangulation.h"
#include "strategies.h"
#include <CGAL/Polygon_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/centroid.h>
#include <CGAL/convex_hull_2.h>
#include <unordered_map>
#include <cmath>

//////////////////////////////////////////////////////////
// kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds, Itag> CDT;
typedef K::Line_2 Line;
typedef CDT::Point Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CDT::Face_handle Face_handle;
using namespace std;
/////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
int is_obtuse_angle(Point &A, Point &B, Point &C)
{
    return angle(A, B, C) == CGAL::OBTUSE;
}

int count_Obtuse_Angles(CDT &cdt)
{
    int count = 0;
    for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
    {
        Point a = fit->vertex(0)->point();
        Point b = fit->vertex(1)->point();
        Point c = fit->vertex(2)->point();
        // cout << "Checking triangle with points: (" << a.x() << ", " << a.y() << "), (" << b.x() << ", " << b.y() << "), (" << c.x() << ", " << c.y() << ")" << endl;

        if (is_obtuse_angle(a, b, c) || is_obtuse_angle(b, c, a) || is_obtuse_angle(c, a, b))
        {
            count++;
            // cout << "Obtuse angle found at vertex " << is_Obtuse << endl;
        }
    }
    return count;
}

Point project_point(Point &A, Point &B, Point &P)
{ // επιστρέφει την προβολή από ένα σημείο P στην πλευρά που σχηματίζουν τα Α-Β
    Line line(A, B);
    return line.projection(P);
}
///////////////////////////////////////////////////////////////

Polygon_2 find_convex_polygon_around_obtuse_triangle(CDT &cdt, Face_handle face)
{
    set<Face_handle> visited_faces;
    queue<Face_handle> face_queue;
    Polygon_2 convex_polygon;
    // Αρχικοποίηση με το αρχικό αμβλυγώνιο τρίγωνο
    face_queue.push(face);
    visited_faces.insert(face);

    while (!face_queue.empty())
    {
        Face_handle current_face = face_queue.front();
        face_queue.pop();

        // Προσθήκη σημείων του current_face στο polygon
        for (int i = 0; i < 3; i++)
        {
            convex_polygon.push_back(current_face->vertex(i)->point());
        }
        // Έλεγχος για τους γειτονικούς κόμβους
        for (int i = 0; i < 3; i++)
        {
            Face_handle neighbor = current_face->neighbor(i);
            if (!cdt.is_infinite(neighbor) && visited_faces.find(neighbor) == visited_faces.end())
            {
                Point a = neighbor->vertex(0)->point();
                Point b = neighbor->vertex(1)->point();
                Point c = neighbor->vertex(2)->point();

                if (is_obtuse_angle(a, b, c) || is_obtuse_angle(b, c, a) || is_obtuse_angle(c, a, b)) // Έλεγχος αν το γειτονικό τρίγωνο είναι αμβλυγώνιο
                {
                    face_queue.push(neighbor);
                    visited_faces.insert(neighbor);
                }
            }
        }
    }
    // Δημιουργία του κυρτού περιβλήματος των σημείων του πολυγώνου
    Polygon_2 convex_hull;
    convex_hull_2(convex_polygon.vertices_begin(), convex_polygon.vertices_end(), back_inserter(convex_hull));
    return convex_hull;
}

// Μέθοδος εισαγωγής Steiner points στο εσωτερικό κυρτών πολυγώνων που σχηματίζονται από αμβλυγώνια τρίγωνα
Point insert_Steiner_point_in_convex_polygons(CDT &cdt, Polygon_2 &region_boundary)
{
    int steiner_counter = 0;

    for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
    {
        Point a = fit->vertex(0)->point();
        Point b = fit->vertex(1)->point();
        Point c = fit->vertex(2)->point();
        // Ελέγχουμε αν το τρίγωνο έχει αμβλεία γωνία
        if (is_obtuse_angle(a, b, c) || is_obtuse_angle(b, c, a) || is_obtuse_angle(c, a, b))
        {
            // Δημιουργία κυρτού πολυγώνου γύρω από το αμβλυγώνιο τρίγωνο και τους γείτονές του
            Polygon_2 convex_polygon = find_convex_polygon_around_obtuse_triangle(cdt, fit);

            // Υπολογισμός του centroid του κυρτού πολυγώνου
            Point centroid = CGAL::centroid(convex_polygon.vertices_begin(), convex_polygon.vertices_end());
            return centroid;
        }
    }
}

///////////////////////////////////////////////////////////

// Συνάρτηση που επιστρέφει σημείο Steiner για μία από τις 5 στρατηγικές
Point select_steiner_point(Point &a, Point &b, Point &c, int strategy, CDT &cdt, Polygon_2 convex_hull)
{
    switch (strategy)
    {
    case 0:
    { // Περίκεντρο check
        return circumcenter(a, b, c);
    }
    case 1:
    { // Κέντρο βάρους ενός τριγώνου
        double cx = (a.x() + b.x() + c.x()) / 3.0;
        double cy = (a.y() + b.y() + c.y()) / 3.0;
        return Point(cx, cy);
    }
    case 2:
    { // Μέσο της μεγαλύτερης ακμής check
        double d_ab = squared_distance(a, b);
        double d_bc = squared_distance(b, c);
        double d_ca = squared_distance(c, a);
        if (d_ab >= d_bc && d_ab >= d_ca)
        {
            return midpoint(a, b);
        }
        else if (d_bc >= d_ab && d_bc >= d_ca)
        {
            return midpoint(b, c);
        }
        else
        {
            return midpoint(c, a);
        }
    }
    case 3:
    { // Προβολή της κορυφής της αμβλείας γωνίας στην απέναντι πλευρά
        if (is_obtuse_angle(b, a, c))
        {
            return project_point(b, c, a); // προβολή του A στην πλευρά B-C
        }
        else if (is_obtuse_angle(a, b, c))
        {
            return project_point(a, c, b); // του B στην πλευρά A-C
        }
        else if (is_obtuse_angle(a, c, b))
        {
            return project_point(a, b, c); // και του C στην πλευρά AB
        }
    }
    case 4:
    {
        return insert_Steiner_point_in_convex_polygons(cdt, convex_hull);
    }
    default:
        throw invalid_argument("Invalid strategy selected.");
    }
}
