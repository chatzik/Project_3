#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/mark_domain_in_triangulation.h>
#include <CGAL/Delaunay_mesher_2.h>
#include "triangulation.h"
#include "strategies.h"
#include <CGAL/Polygon_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/centroid.h>
#include <CGAL/convex_hull_2.h>
#include <unordered_map>
#include <random>
#include <cmath>
#include "old_triangulation.h"


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

// State definition
struct State
{
    CDT cdt;
    int obtuse_count;
    int steiner_points;
    vector<Point> steiner_locations;
    vector<int> strategies;

    bool operator==(const State &other) const
    {
        return obtuse_count == other.obtuse_count &&
               steiner_points == other.steiner_points &&
               steiner_locations == other.steiner_locations &&
               strategies == other.strategies;
    }
};
bool is_obtuse(Face_handle face) {
    Point a = face->vertex(0)->point();
    Point b = face->vertex(1)->point();
    Point c = face->vertex(2)->point();
    return (CGAL::angle(a, b, c) == CGAL::OBTUSE || 
            CGAL::angle(b, c, a) == CGAL::OBTUSE || 
            CGAL::angle(c, a, b) == CGAL::OBTUSE);
}

Face_handle choose_triangle(const vector<Face_handle>& obtuse_triangles, 
                            const map<Point, double>& pheromones,
                            double pheromone_importance, 
                            double heuristic_importance) {
    return obtuse_triangles[rand() % obtuse_triangles.size()];
}

vector<Point> generate_steiner_candidates(Face_handle triangle) {

    vector<Point> candidates;
    Point a = triangle->vertex(0)->point();
    Point b = triangle->vertex(1)->point();
    Point c = triangle->vertex(2)->point();
    
    candidates.push_back(CGAL::centroid(a, b, c));
    
    Point circumcenter = CGAL::circumcenter(a, b, c);
    if (CGAL::bounded_side_2(a, b, c, circumcenter) == CGAL::ON_BOUNDED_SIDE) {
        candidates.push_back(circumcenter);
    }
    
    candidates.push_back(CGAL::midpoint(a, b));
    candidates.push_back(CGAL::midpoint(b, c));
    candidates.push_back(CGAL::midpoint(c, a));
    
    return candidates;
}

Point choose_steiner_point(const vector<Point>& candidates, 
                           const map<Point, double>& pheromones,
                           double pheromone_importance, 
                           double heuristic_importance) {
    // Επιστρέφει τυχαίο candidate - πρέπει να αλλάξει με την ανάλογη πιθανότητα στις διαφάνειες ΔΕΣ ΤΟ!
    return candidates[rand() % candidates.size()];
}

void update_pheromones(map<Point, double>& pheromones, 
                       const vector<CDT>& ant_solutions, 
                       const vector<int>& ant_obtuse_counts,
                       double evaporation_rate) {
    // Evaporate pheromones
    for (auto& p : pheromones) {
        p.second *= (1 - evaporation_rate);
    }

    // βαζουμε pheromones με βαση τα μεμυγκια
    for (size_t i = 0; i < ant_solutions.size(); ++i) {
        double pheromone_deposit = 1.0 / (ant_obtuse_counts[i] + 1);  // Avoid division by zero
        for (auto vit = ant_solutions[i].finite_vertices_begin(); vit != ant_solutions[i].finite_vertices_end(); ++vit) {
            Point p = vit->point();
            pheromones[p] += pheromone_deposit;
        }
    }
}
State aco_triangulation(const vector<int>& points_x, const vector<int>& points_y, 
                        const vector<int>& region_boundary, 
                        const vector<pair<int, int>>& additional_constraints,
                        double alpha, double beta, int L, 
                        int ant_count, double pheromone_importance, 
                        double heuristic_importance, double evaporation_rate)
{
    CDT cdt;
    vector<Point> points;
    for (size_t i = 0; i < points_x.size(); ++i) {
        points.push_back(Point(points_x[i], points_y[i]));
    }

    // Add constraints
    for (size_t i = 0; i < region_boundary.size(); ++i) {
        int next = (i + 1) % region_boundary.size();
        cdt.insert_constraint(points[region_boundary[i]], points[region_boundary[next]]);
    }
    for (const auto& constraint : additional_constraints) {
        cdt.insert_constraint(points[constraint.first], points[constraint.second]);
    }

    // Initialize pheromone trails
    std::map<Point, double> pheromones;

    // αρχικοποίηση best
    State best_state;
    best_state.cdt = cdt;
    best_state.obtuse_count = count_Obtuse_Angles(cdt);
    

        // Main ACO loop
    for (int iteration = 0; iteration < L; ++iteration) {
        vector<CDT> ant_solutions(ant_count, cdt);
        vector<int> ant_obtuse_counts(ant_count);

        // καθε μεμυγκι φτιαχνει μια λυση
        for (int ant = 0; ant < ant_count; ++ant) {
            CDT& current_cdt = ant_solutions[ant];
            int current_obtuse_count = count_Obtuse_Angles(current_cdt);
            int steiner_count = 0;

            // το μερμυγκι προχωραει μεχρι να μην μπορει να βελτιωσει αλλο η δεν εχει αλλα βηματα
            for (int step = 0; step < 100; ++step) {
                // βρισκουμε τα αμβλεια τριγωνα
                vector<Face_handle> obtuse_triangles;
                for (auto fit = current_cdt.finite_faces_begin(); fit != current_cdt.finite_faces_end(); ++fit) {
                    if (is_obtuse(fit)) {
                        obtuse_triangles.push_back(fit);
                    }
                }

                if (obtuse_triangles.empty()) break;

                // Διαλέγει τυχαία αμβλεία - τρίγωνο
                Face_handle chosen_triangle = obtuse_triangles[rand() % obtuse_triangles.size()];

                // κανουμε genarate cadidates- thelei kiallo edw
                vector<Point> candidates = generate_steiner_candidates(chosen_triangle);

                // επιλεγουμε steiner point με βαση το pheromones 
                Point chosen_steiner = choose_steiner_point(candidates, pheromones, pheromone_importance, heuristic_importance);

                //β
                current_cdt.insert(chosen_steiner);
                steiner_count++;

                // 
                int new_obtuse_count = count_Obtuse_Angles(current_cdt);
                if (new_obtuse_count >= current_obtuse_count) break;
                current_obtuse_count = new_obtuse_count;
            }

            ant_obtuse_counts[ant] = current_obtuse_count;

            // ενημερωνουμε το  best solution αν χρειαζεται
            if (current_obtuse_count < best_state.obtuse_count) {
                best_state.cdt = current_cdt;
                best_state.obtuse_count = current_obtuse_count;
                
            }
        }

        // ενημερωνουμε pheromones
        update_pheromones(pheromones, ant_solutions, ant_obtuse_counts, evaporation_rate);
    }

    return best_state;}