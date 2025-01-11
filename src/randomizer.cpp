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

State local_triangulation(CDT &initial_cdt, Polygon_2 &convex_hull, int &best_obtuse, CDT &best_cdt, int max_iterations)
{
    queue<State> queue;
    unordered_set<State, StateHash> visited; // Custom hash for State

    // Initialize with the initial state
    State initial_state = {initial_cdt, count_Obtuse_Angles(initial_cdt), 0, {}, {}};
    State best_state = initial_state;
    queue.push(initial_state);
    visited.insert(initial_state);
    int iteration_count = 0;
    int stagnation_count = 0; // Counter for stagnation
    best_cdt = initial_cdt;

    // Exploration via local strategies
    while (!queue.empty() && iteration_count < max_iterations && best_state.obtuse_count > 0)
    {
        State current_state = queue.front();
        queue.pop();

        // If the current state is better, update the best solution
        if (current_state.obtuse_count < best_state.obtuse_count)
        {
            best_cdt = current_state.cdt;
            best_state = current_state;
            iteration_count -= 10; // Reset iteration count due to improvement
            stagnation_count = 0;  // Reset stagnation count
        }

        // If max iterations reached or no obtuse angles remain, stop
        if (iteration_count >= max_iterations || best_state.obtuse_count == 0)
            return best_state;

        // Explore all triangles with obtuse angles
        bool improvement_found = false;
        for (auto fit = current_state.cdt.finite_faces_begin(); fit != current_state.cdt.finite_faces_end(); ++fit)
        {
            Point a = fit->vertex(0)->point();
            Point b = fit->vertex(1)->point();
            Point c = fit->vertex(2)->point();

            // Test all strategies
            for (int strategy = 0; strategy < 5; ++strategy)
            {
                Point steiner = select_steiner_point(a, b, c, strategy, current_state.cdt, convex_hull);

                // Check if the point is within the convex hull
                if (convex_hull.bounded_side(steiner) == CGAL::ON_BOUNDED_SIDE || convex_hull.bounded_side(steiner) == CGAL::ON_BOUNDARY)
                {
                    CDT temp_cdt = current_state.cdt;
                    temp_cdt.insert(steiner);
                    int new_obtuse = count_Obtuse_Angles(temp_cdt);

                    // Create a new state
                    State new_state = {temp_cdt, new_obtuse, current_state.steiner_points + 1, {}, {}};

                    // If the new state is better and hasn't been visited
                    if (new_obtuse <= initial_state.obtuse_count && visited.find(new_state) == visited.end())
                    {
                        queue.push(new_state);
                        visited.insert(new_state);
                        improvement_found = true;
                    }
                }
            }
        }

        // Randomization: Trigger after 200 consecutive stagnation iterations
        if (!improvement_found)
        {
            stagnation_count++;
            if (stagnation_count >= 1)
            {
                Point random_point = generate_random_point_within_hull(convex_hull);
                cout << "Randomization triggered. Generated random point: (" << random_point.x() << ", " << random_point.y() << ")\n";

                if (convex_hull.bounded_side(random_point) == CGAL::ON_BOUNDED_SIDE || convex_hull.bounded_side(random_point) == CGAL::ON_BOUNDARY)
                {
                    // Insert the random point into the best CDT
                    best_cdt.insert(random_point);
                    int new_obtuse = count_Obtuse_Angles(best_cdt);

                    // Update the best state
                    best_state.cdt = best_cdt;
                    best_state.obtuse_count = new_obtuse;
                    best_state.steiner_points++;
                    best_state.steiner_locations.push_back(random_point);

                    stagnation_count = 0; // Reset stagnation count after randomization
                    cout << "Random point accepted and added to the best CDT. New obtuse count: " << new_obtuse << "\n";
                }
                else
                {
                    cout << "Random point is outside the convex hull.\n";
                }
            }
        }

        // Increment the iteration count if no improvement
        if (current_state.obtuse_count >= best_state.obtuse_count)
        {
            iteration_count++;
        }
    }

    return best_state;
}

State sa_triangulation(CDT &cdt, const Polygon_2 &convex_hull, int initial_obtuse, CDT &best_cdt, double alpha, double beta, int L)
{
    // Energy function for evaluating states
    auto energy = [alpha, beta](const CDT &triangulation, int initial_vertices)
    {
        int obtuse_count = count_Obtuse_Angles(const_cast<CDT &>(triangulation));
        int steiner_count = triangulation.number_of_vertices() - initial_vertices;
        return alpha * obtuse_count + beta * steiner_count;
    };

    // Initialize the current state
    State current_state;
    current_state.cdt = cdt;
    current_state.obtuse_count = initial_obtuse;
    current_state.steiner_points = 0;

    // Set the initial state as the best state
    State best_state = current_state;
    double current_energy = energy(cdt, cdt.number_of_vertices());
    double best_energy = current_energy;

    // Random number generation setup
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0.0, 100.0);

    // Temperature for simulated annealing
    double temperature = 1.0;

    // Counter for stagnation
    int stagnation_count = 0;

    // Main loop of simulated annealing
    while (temperature >= 0)
    {
        // Perform L iterations for each temperature
        for (int i = 0; i < L; ++i)
        {
            // Identify all obtuse triangles
            vector<CDT::Face_handle> obtuse_faces;
            for (auto fit = current_state.cdt.finite_faces_begin(); fit != current_state.cdt.finite_faces_end(); ++fit)
            {
                if (is_obtuse_triangle(fit))
                {
                    obtuse_faces.push_back(fit);
                }
            }

            // If no obtuse triangles remain, terminate
            if (obtuse_faces.empty())
                break;

            // Select a random obtuse triangle
            uniform_int_distribution<> face_dis(0, obtuse_faces.size() - 1);
            int selected_index = face_dis(gen);
            auto selected_face = obtuse_faces[selected_index];

            // Retrieve the vertices of the selected triangle
            Point a = selected_face->vertex(0)->point();
            Point b = selected_face->vertex(1)->point();
            Point c = selected_face->vertex(2)->point();

            // Choose a random strategy out of 5
            int strategy = gen() % 5;
            Point steiner = select_steiner_point(a, b, c, strategy, current_state.cdt, convex_hull);

            // Check if the Steiner point is inside or on the boundary
            if (convex_hull.bounded_side(steiner) == CGAL::ON_BOUNDED_SIDE ||
                convex_hull.bounded_side(steiner) == CGAL::ON_BOUNDARY)
            {
                // Create a new state with the Steiner point
                State new_state = current_state;
                CDT::Vertex_handle v = new_state.cdt.insert(steiner);

                // Proceed only if the insertion was successful
                if (v != CDT::Vertex_handle())
                {
                    new_state.obtuse_count = count_Obtuse_Angles(new_state.cdt);
                    new_state.steiner_points++;
                    new_state.steiner_locations.push_back(steiner);

                    // Compute the energy of the new state
                    double new_energy = energy(new_state.cdt, cdt.number_of_vertices());
                    double delta_energy = new_energy - best_energy;

                    // Metropolis criterion
                    if (delta_energy < 0 || exp(-delta_energy / temperature) > dis(gen))
                    {
                        current_state = new_state;
                        current_energy = new_energy;

                        // Update the best state if the new state is better
                        if (current_energy < best_energy)
                        {
                            best_state = current_state;
                            best_energy = current_energy;
                            best_cdt = current_state.cdt;
                            stagnation_count = 0; // Reset stagnation count
                        }
                    }
                }
            }

            // Remove the triangle from the obtuse_faces list
            obtuse_faces.erase(obtuse_faces.begin() + selected_index);
        }

        // Randomization: Insert a random Steiner point if stagnation persists for 100 iterations
        stagnation_count++;
        if (stagnation_count >= 8)
        {
            Point random_point = generate_random_point_within_hull(convex_hull);
            cout << "SA Randomization triggered. Generated random point: (" << random_point.x() << ", " << random_point.y() << ")\n";

            if (convex_hull.bounded_side(random_point) == CGAL::ON_BOUNDED_SIDE ||
                convex_hull.bounded_side(random_point) == CGAL::ON_BOUNDARY)
            {
                CDT::Vertex_handle v = current_state.cdt.insert(random_point);
                if (v != CDT::Vertex_handle())
                {
                    current_state.obtuse_count = count_Obtuse_Angles(current_state.cdt);
                    current_state.steiner_points++;
                    current_state.steiner_locations.push_back(random_point);

                    cout << "SA Random point accepted. New obtuse count: " << current_state.obtuse_count << "\n";

                    // Reset stagnation count after successful randomization
                    stagnation_count = 0;
                }
                else
                {
                    cout << "SA Random point insertion failed.\n";
                }
            }
            else
            {
                cout << "SA Random point is outside the convex hull.\n";
            }
        }

        // Decrease the temperature
        temperature -= 1.0 / 0.95;
    }

    // Return the best state found
    return best_state;
}

Point generate_random_point_within_hull(const Polygon_2 &convex_hull)
{
    // Generate a random point inside the bounding box of the convex hull
    auto bbox = CGAL::bbox_2(convex_hull.vertices_begin(), convex_hull.vertices_end());
    double x = CGAL::to_double(bbox.xmin()) + (rand() / (RAND_MAX + 1.0)) * (CGAL::to_double(bbox.xmax()) - CGAL::to_double(bbox.xmin()));
    double y = CGAL::to_double(bbox.ymin()) + (rand() / (RAND_MAX + 1.0)) * (CGAL::to_double(bbox.ymax()) - CGAL::to_double(bbox.ymin()));
    return Point(x, y);
}