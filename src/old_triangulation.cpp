#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/mark_domain_in_triangulation.h>
#include <CGAL/Delaunay_mesher_2.h>
#include "triangulation.h"
#include <CGAL/Polygon_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/centroid.h>
#include <CGAL/convex_hull_2.h>
#include <unordered_map>
#include <cmath>
#include "strategies.h"
#include "old_triangulation.h"

//////////////////////////////////////////////////////////

// kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds, Itag> CDT;
typedef CDT::Point Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CDT::Face_handle Face_handle;
using namespace std;

///////////////////////////////////////////////////////////

double calculate_angle_Old(const Point &A, const Point &B, const Point &C)
{
    double a2 = squared_distance(B, C); // απέναντι πλευρά από το A
    double b2 = squared_distance(A, C); // απέναντι πλευρά από το B
    double c2 = squared_distance(A, B); // απέναντι πλευρά από το C

    // νόμος των συνημιτόνων: cos(γωνία) = (b² + c² - a²) / (2*ρίζα(b)*ρίζα(c))
    double cos_A = (b2 + c2 - a2) / (2 * sqrt(b2) * sqrt(c2));

    // μεατροπή της γωνίας σε μοίρες
    return acos(cos_A) * 180.0 / M_PI;
}

// η collinear ελέγχει αν τα 3 σημεία που εισάγαμε είναι συγγραμικά, δηλαδή δεν σχηματίζουν τρίγωνο
// Η Has_Obuse_angle_old κάνει έλεγχο αν το τρίγωνο έχει έστω και μία αμβλεία γωνία
int Has_Obuse_angle_old(const Point &a, const Point &b, const Point &c)
{
    if (collinear(a, b, c))
    {
        return -1; // τα σημεία είναι συγγραμμικά, οπότε δεν σχηματίζουν τρίγωνο
    }
    // υπολογισμοσ  των γωνίων του τριγώνου
    double angle_A = calculate_angle_Old(a, b, c);
    double angle_B = calculate_angle_Old(b, a, c);
    double angle_C = calculate_angle_Old(c, a, b);

    // έλεγχος για αμβλεία γωνία
    if (angle_A > 90)
    {
        // cout << "Obtuse angle found at A! " << angle_A << endl;
        return 0; // Αμβλεία γωνία στο A
    }
    if (angle_B > 90)
    {
        // cout << "Obtuse angle found at B! " << angle_B << endl;
        return 1; // Αμβλεία γωνία στο B
    }
    if (angle_C > 90)
    {
        // cout << "Obtuse angle found at C! " << angle_C << endl;
        return 2; // Αμβλεία γωνία στο C
    }

    return -1; // Δεν υπάρχει αμβλεία γωνία
}

Point insert_Steiner(const Point &a, const Point &b) // επιστρέφει το μέσο της απέναντι πλευράς
{
    return Point((a.x() + b.x()) / 2, (a.y() + b.y()) / 2);
}

// επιστρέφει το πλήθος των αμβλείων γωνιών του γράφου
int count_Obtuse_Angles_old(CDT &cdt)
{
    int count = 0;
    for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
    {
        Point a = fit->vertex(0)->point();
        Point b = fit->vertex(1)->point();
        Point c = fit->vertex(2)->point();
        // cout << "Checking triangle with points: (" << a.x() << ", " << a.y() << "), (" << b.x() << ", " << b.y() << "), (" << c.x() << ", " << c.y() << ")" << endl;

        if (int is_Obtuse = Has_Obuse_angle_old(a, b, c) != -1)
        {
            count++;
            // cout << "Obtuse angle found at vertex " << is_Obtuse << endl;
        }
    }
    return count;
}

// Μέθοδος για εύρεση γειτονικών αμβλυγώνιων finite faces και δημιουργία του κυρτού περιβλήματος
Polygon_2 find_convex_polygon_around_obtuse_triangle_Old(CDT &cdt, Face_handle face)
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

                if (Has_Obuse_angle_old(a, b, c) != -1) // Έλεγχος αν το γειτονικό τρίγωνο είναι αμβλυγώνιο
                {
                    face_queue.push(neighbor);
                    visited_faces.insert(neighbor);
                }
            }
        }
    }

    // Δημιουργία του κυρτού περιβλήματος των σημείων του πολυγώνου
    Polygon_2 convex_hull;
    CGAL::convex_hull_2(convex_polygon.vertices_begin(), convex_polygon.vertices_end(), back_inserter(convex_hull));

    return convex_hull;
}

// Μέθοδος εισαγωγής Steiner points στο εσωτερικό κυρτών πολυγώνων που σχηματίζονται από αμβλυγώνια τρίγωνα
void insert_Steiner_points_in_convex_polygons_old(CDT &cdt, const Polygon_2 &region_boundary)
{
    int steiner_counter = 0;

    for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
    {
        Point a = fit->vertex(0)->point();
        Point b = fit->vertex(1)->point();
        Point c = fit->vertex(2)->point();

        // Ελέγχουμε αν το τρίγωνο έχει αμβλεία γωνία
        if (Has_Obuse_angle_old(a, b, c) != -1)
        {
            // Δημιουργία κυρτού πολυγώνου γύρω από το αμβλυγώνιο τρίγωνο και τους γείτονές του
            Polygon_2 convex_polygon = find_convex_polygon_around_obtuse_triangle_Old(cdt, fit);

            // Υπολογισμός του centroid του κυρτού πολυγώνου
            Point centroid = CGAL::centroid(convex_polygon.vertices_begin(), convex_polygon.vertices_end());

            // Έλεγχος αν το centroid βρίσκεται εντός του αρχικού κυρτού περιβλήματος
             if (region_boundary.bounded_side(centroid) == CGAL::ON_BOUNDED_SIDE||region_boundary.bounded_side(centroid) == CGAL::ON_BOUNDARY)
            {
                CDT temp_cdt = cdt;
                temp_cdt.insert(centroid);

                // Μέτρηση αμβλεών γωνιών πριν και μετά
                int original_obtuse_count = count_Obtuse_Angles_old(cdt);
                int new_obtuse_count = count_Obtuse_Angles_old(temp_cdt);
                if (new_obtuse_count < original_obtuse_count)
                {
                    cdt.insert(centroid); // Εισαγωγή του Steiner point στην τριγωνοποίηση
                    steiner_counter++;
                    //cout << "Inserted Steiner point in polygon centroid: (" << centroid.x() << ", " << centroid.y() << ")" << endl;
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// τριγωνοποίηση
void old_triangulate(const vector<int> &points_x, const vector<int> &points_y, const vector<int> &region_boundary, const vector<pair<int, int>> &additional_constraints)
{
    CDT cdt;

    // προσθήκη σημείων και ακμών του γράφου
    vector<Point> points;
    for (size_t i = 0; i < points_x.size(); ++i)
    {
        points.push_back(Point(points_x[i], points_y[i]));
    }
    // Δημιουργία του πολυγώνου του κυρτού περιβλήματος
    Polygon_2 convex_hull;
    for (size_t i : region_boundary)
    {
        convex_hull.push_back(points[i]);
    }

    for (int i = 0; i < points_x.size(); i++)
        cout << points_x[i] << " ";
    for (int i = 0; i < points_y.size(); i++)
        cout << points_y[i] << " ";

    // προσθήκη των ακμών του κυρτού περιβλήματος για όριο
    for (size_t i = 0; i < region_boundary.size(); ++i)
    {
        int next = (i + 1) % region_boundary.size();
        cdt.insert_constraint(points[region_boundary[i]], points[region_boundary[next]]);
    }

    // προσθήκη πρόσθετων constraints
    for (const auto &constraint : additional_constraints)
    {
        cdt.insert_constraint(points[constraint.first], points[constraint.second]);
    }

    int steiner_counter = 0;
    int original_graph_count = count_Obtuse_Angles_old(cdt); // μετρητής για το πόσες αμβλείες γωνίες υπάρχουν
    cout << "Number of obtuse angles in the triangulation (before using any methods):" << original_graph_count << endl;
    if (original_graph_count > 0)
    {
        bool found_steiner_point = true;
        while (found_steiner_point)
        {
            found_steiner_point = false; // aν δεν βρεθεί "καλό σημείο" για να βάλουμε steiner, το loop θα σταματήσει
            int original_graph_count = count_Obtuse_Angles_old(cdt);
            // Κάνουμε iterate στα finite faces της τριγωνοποίησης
            for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
            {
                Point a = fit->vertex(0)->point();
                Point b = fit->vertex(1)->point();
                Point c = fit->vertex(2)->point();

                int is_Obtuse = Has_Obuse_angle_old(a, b, c);
                if (is_Obtuse != -1)
                { // δηλαδή η γωνία ΔΕΝ είναι οξεία ή κάθετη -> αμβλεία
                    insert_Steiner_points_in_convex_polygons_old(cdt, convex_hull);
                    Point steiner;
                    if (is_Obtuse == 0)
                        steiner = insert_Steiner(b, c); // απέναντι πλευρά από το A
                    else if (is_Obtuse == 1)
                        steiner = insert_Steiner(a, c); // απέναντι πλευρά από το σημείο B
                    else                                // δλδ, is_Obtuse = 2
                        steiner = insert_Steiner(a, b); // απέναντι πλευρά από το C
                             if (convex_hull.bounded_side(steiner) == CGAL::ON_BOUNDED_SIDE|| convex_hull.bounded_side(steiner) ==CGAL::ON_BOUNDARY)
                    {
                        CDT temp_cdt = cdt;
                        temp_cdt.insert(steiner);
                        int new_obtuse_count = count_Obtuse_Angles_old(temp_cdt);
                        if (new_obtuse_count < original_graph_count)
                        {
                            cdt.insert(steiner);
                            found_steiner_point = true;
                            cout << "Inserted Steiner point at: (" << steiner.x() << ", " << steiner.y() << ")" << endl;
                            steiner_counter++;
                            break;
                        }
                    }
                    /////////////////////////////////////////////////////////////////// Μέθοδος 2
                    // Χρήση της έτοιμης συνάρτησης της CGAL circumcenter για υπολογισμό περικέντρου
                    steiner = circumcenter(a, b, c);
                    if (convex_hull.bounded_side(steiner) == CGAL::ON_BOUNDED_SIDE||convex_hull.bounded_side(steiner) == CGAL::ON_BOUNDARY )
                    {
                        CDT temp_cdt = cdt;
                        temp_cdt.insert(steiner);
                        int new_obtuse_count = count_Obtuse_Angles_old(temp_cdt);
                        if (new_obtuse_count < original_graph_count)
                        {
                            cdt.insert(steiner);
                            found_steiner_point = true;
                            cout << "Inserted Steiner point at circumcenter: (" << steiner.x() << ", " << steiner.y() << ")" << endl;
                            steiner_counter++;
                            break;
                        }
                        else
                        {
                            cout << "Skipped Steiner point at circumcenter (out of convex hull): (" << steiner.x() << ", " << steiner.y() << ")" << endl;
                        }
                    }
                }
            }
        }
        original_graph_count = count_Obtuse_Angles_old(cdt);
    }

        cout << "Number of obtuse angles in the triangulation: " << original_graph_count << endl;
    cout << "Number of steiner points added in the triangulation: " << steiner_counter << endl;
    draw(cdt);
}