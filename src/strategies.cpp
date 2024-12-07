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
    Point steiner_point;
    bool valid_point = false;
    int attempts = 0;
    const int max_attempts = 5;

    while (!valid_point && attempts < max_attempts)
    {
        try
        {
            switch (strategy)
            {
            case 0: // Circumcenter
                steiner_point = CGAL::circumcenter(a, b, c);
                break;
            case 1: // Centroid
            {
                double cx = (a.x() + b.x() + c.x()) / 3.0;
                double cy = (a.y() + b.y() + c.y()) / 3.0;
                steiner_point = Point(cx, cy);
                break;
            }
            case 2: // Midpoint of longest edge
            {
                double d_ab = CGAL::squared_distance(a, b);
                double d_bc = CGAL::squared_distance(b, c);
                double d_ca = CGAL::squared_distance(c, a);
                if (d_ab >= d_bc && d_ab >= d_ca)
                    steiner_point = CGAL::midpoint(a, b);
                else if (d_bc >= d_ab && d_bc >= d_ca)
                    steiner_point = CGAL::midpoint(b, c);
                else
                    steiner_point = CGAL::midpoint(c, a);
                break;
            }
            case 3: // Projection of obtuse angle vertex
            {
                if (is_obtuse_angle(b, a, c))
                    steiner_point = project_point(b, c, a);
                else if (is_obtuse_angle(a, b, c))
                    steiner_point = project_point(a, c, b);
                else if (is_obtuse_angle(a, c, b))
                    steiner_point = project_point(a, b, c);
                else
                    steiner_point = CGAL::centroid(a, b, c); // Fallback to centroid if no obtuse angle
                break;
            }
            case 4: // Convex polygon method
                steiner_point = insert_Steiner_point_in_convex_polygons(cdt, convex_hull);
                break;
            default:
                throw std::invalid_argument("Invalid strategy selected.");
            }

            // Check if the point is inside the triangle
            if (CGAL::side_of_bounded_circle(a, b, c, steiner_point) == CGAL::ON_BOUNDED_SIDE) 
            {
                valid_point = true;
            }
            else
            {
                // If not inside, fall back to centroid
                steiner_point = CGAL::centroid(a, b, c);
                valid_point = true;
            }
        }
        catch (const CGAL::Precondition_exception& e)
        {
            // If a CGAL precondition fails, try the next strategy
            strategy = (strategy + 1) % 5;
        }
        catch (const std::exception &e)
        {
            // For any other exception, fall back to centroid
            steiner_point = CGAL::centroid(a, b, c);
            valid_point = true;
        }

        attempts++;
    }

    if (!valid_point)
    {
        // If we've exhausted all attempts, use the centroid as a last resort
        steiner_point = CGAL::centroid(a, b, c);
    }

    return steiner_point;
}