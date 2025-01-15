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
#include <chrono>

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
bool randomization = true;
// State definition
struct State
{
    CDT cdt;
    int obtuse_count;
    int steiner_points;
    vector<Point> steiner_locations;
    vector<int> strategies;
    vector<int> obtuse_history;
    // Υπερφορτωση για το ==
    bool operator==(const State &other) const
    {
        return obtuse_count == other.obtuse_count &&
               steiner_points == other.steiner_points &&
               steiner_locations == other.steiner_locations &&
               strategies == other.strategies;
    }
};
double calculate_convergence_rate(int n, int obtuse_n, int obtuse_n_plus_1) {
    if (n == 0 || obtuse_n == obtuse_n_plus_1) {
        return 0.0;  // Avoid division by zero or log of 1
    }
    double ln_obtuse_ratio = log(static_cast<double>(obtuse_n_plus_1)) - log(static_cast<double>(obtuse_n));
    double ln_n_ratio = log(static_cast<double>(n + 1)) - log(static_cast<double>(n));
    return ln_obtuse_ratio / ln_n_ratio;
}
double calculate_energy(const State& state, double alpha, double beta) {
    return alpha * state.obtuse_count + beta * state.steiner_points;
}
Point generate_random_point_within_hull(const Polygon_2 &convex_hull)
{
    // Use the bounding box of the convex hull
    CGAL::Bbox_2 bbox = convex_hull.bbox();

    // Create a random number generator with a time-based seed
    static unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    static std::mt19937 gen(seed);
    std::uniform_real_distribution<> dis_x(bbox.xmin(), bbox.xmax());
    std::uniform_real_distribution<> dis_y(bbox.ymin(), bbox.ymax());

    while (true) {
        double x = dis_x(gen);
        double y = dis_y(gen);
        Point p(x, y);

        if (convex_hull.bounded_side(p) == CGAL::ON_BOUNDED_SIDE) {
            return p;
        }
    }
}

// εκτυπωνουμε τα στοιχεία της κατάστάσης
void printStateDetails(State &state)
{
    cout << "Obtuse Count: " << state.obtuse_count << endl;
    cout << "Steiner Points: " << state.steiner_points << endl;
    cout << endl;
}
// Function to compare two points with epsilon tolerance
bool points_equal(const Point& a, const Point& b, double epsilon = 1e-10) {
    return (std::abs(a.x() - b.x()) < epsilon) && (std::abs(a.y() - b.y()) < epsilon);
}
// Συνάρτηση για τη σύγκριση δύο καταστάσεων State
bool compareStates(const State &a, const State &b)
{
    return a.obtuse_count == b.obtuse_count &&
           a.steiner_points == b.steiner_points &&
           a.steiner_locations == b.steiner_locations &&
           a.strategies == b.strategies;
}

struct StateHash
{
    size_t operator()(const State &s) const
    {
        size_t hash_val = 0;
        // Hash για obtuse_count, steiner_points, και άλλα χαρακτηριστικά
        hash_val ^= hash<int>{}(s.obtuse_count) + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2);
        hash_val ^= hash<int>{}(s.steiner_points) + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2);
        for (const auto &loc : s.steiner_locations)
        {
            hash_val ^= hash<double>{}(loc.x()) + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2);
            hash_val ^= hash<double>{}(loc.y()) + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2);
        }
        return hash_val;
    }
};
// ελενχος αν τα τριγωνα ειναι αμβλειγωνια
bool is_obtuse_triangle(CDT::Face_handle face)
{
    Point a = face->vertex(0)->point();
    Point b = face->vertex(1)->point();
    Point c = face->vertex(2)->point();
    return is_obtuse_angle(a, b, c) || is_obtuse_angle(b, c, a) || is_obtuse_angle(c, a, b);
}

State local_triangulation(CDT &initial_cdt, Polygon_2 &convex_hull, int &best_obtuse, CDT &best_cdt, int max_iterations)
{
    queue<State> queue;
    unordered_set<State, StateHash> visited; // Custom hash for State
    auto start_time = std::chrono::high_resolution_clock::now();
    const int TIME_LIMIT = 500; // 8 minutes in seconds


    // Initialize with the initial state
    State initial_state = {initial_cdt, count_Obtuse_Angles(initial_cdt), 0, {}, {},{count_Obtuse_Angles(initial_cdt)}};
    State best_state = initial_state;
    queue.push(initial_state);
    visited.insert(initial_state);
    int iteration_count = 0;
    int stagnation_count = 0; // Counter for stagnation
    best_cdt = initial_cdt;
    

    // Exploration via local strategies
    while (!queue.empty() && iteration_count < max_iterations && best_state.obtuse_count > 0)
    {
        auto current_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time);
        if (duration.count() >= TIME_LIMIT)
        {
            std::cout << "Time limit reached. Stopping local search." << std::endl;
            break;
        }
        State current_state = queue.front();
        queue.pop();

        // If the current state is better, update the best solution
        if (current_state.obtuse_count < best_state.obtuse_count)
        {
            best_cdt = current_state.cdt;
            best_state = current_state;
            best_state.steiner_points ++;
            iteration_count -= 10; // Reset iteration count due to improvement
            stagnation_count = 0;  // Reset stagnation count
            best_state.obtuse_history.push_back(best_state.obtuse_count);
            
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
            {   stagnation_count++;
                Point steiner = select_steiner_point(a, b, c, strategy, current_state.cdt, convex_hull);
                
                // Check if the point is within the convex hull
                if (convex_hull.bounded_side(steiner) == CGAL::ON_BOUNDED_SIDE || convex_hull.bounded_side(steiner) == CGAL::ON_BOUNDARY)
                {
                    CDT temp_cdt = current_state.cdt;
                    temp_cdt.insert(steiner);
                    int new_obtuse = count_Obtuse_Angles(temp_cdt);

                    // Create a new state
                    State new_state = {temp_cdt, new_obtuse, current_state.steiner_points + 1, {}, {}, current_state.obtuse_history};
                    new_state.obtuse_history.push_back(new_obtuse);
                    // If the new state is better and hasn't been visited
                    if (new_obtuse <= initial_state.obtuse_count && visited.find(new_state) == visited.end())
                    {
                        queue.push(new_state);
                        visited.insert(new_state);
                        improvement_found = true;
                        stagnation_count = 0;
                    }
                }
            }
        }
        
        // Randomization: Trigger after 200 consecutive stagnation iterations
        if(!improvement_found && randomization)
        {
            stagnation_count++;
            if (stagnation_count >= 100)
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
                    cout << "New best solution found: Obtuse Count = " << best_state.obtuse_count << endl;
                    best_state.obtuse_history.push_back(new_obtuse);
                    cout<< "after the first update of best_state.obtuse_history"<<endl;



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
    auto start_time = std::chrono::high_resolution_clock::now();
    const int TIME_LIMIT = 500; // 8 minutes in seconds
    current_state.obtuse_count = initial_obtuse;
    current_state.steiner_points = 0;
    current_state.obtuse_history.push_back(initial_obtuse);

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
        auto current_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time);
        if (duration.count() >= TIME_LIMIT)
        {
            std::cout << "Time limit reached. Stopping local search." << std::endl;
            break;
        }
        // Perform L iterations for each temperature
        for (int i = 0; i < L; ++i)
        {
            //stagnation_count++;
    
            if (stagnation_count >= 500 && randomization == true){
                
                // Generate a random Steiner point
                Point random_point = generate_random_point_within_hull(convex_hull);
                cout << "Randomization triggered. Generated random point: (" << random_point.x() << ", " << random_point.y() << ")\n";
                // Check if the random Steiner point is inside or on the boundary
                if (convex_hull.bounded_side(random_point) == CGAL::ON_BOUNDED_SIDE ||
                    convex_hull.bounded_side(random_point) == CGAL::ON_BOUNDARY)
                {
                    // Create a new state with the random Steiner point
                    State new_state = current_state;
                    CDT::Vertex_handle v = new_state.cdt.insert(random_point);
                    
                    // Proceed only if the insertion was successful
                    if (v != CDT::Vertex_handle())
                    {
                        best_cdt.insert(random_point);
                    int new_obtuse = count_Obtuse_Angles(best_cdt);

                    // Update the best state
                    best_state.cdt = best_cdt;
                    best_state.obtuse_count = new_obtuse;
                    best_state.steiner_points++;
                    best_state.steiner_locations.push_back(random_point);
                    cout << "New best solution found: Obtuse Count = " << best_state.obtuse_count << endl;
                    best_state.obtuse_history.push_back(new_obtuse);
                    stagnation_count = 0; // Reset stagnation count
                    }
                    else
                    {
                        cout << "Random Steiner point insertion failed.\n";
                    }
                }
                else
                {
                    cout << "Random Steiner point is outside the convex hull.\n";
                }
                
                continue; 
            }
            else{
                stagnation_count++;
            }
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
                    new_state.steiner_locations.push_back(steiner);
                    new_state.obtuse_history.push_back(new_state.obtuse_count);
                    // Compute the energy of the new state
                    double new_energy = energy(new_state.cdt, cdt.number_of_vertices());
                    double delta_energy = new_energy - best_energy;
                    if (new_energy < best_energy)
                    {
                    best_state.steiner_points++;  // This will increment the number of Steiner points in the best state
                    best_energy = new_energy;
                    best_cdt = new_state.cdt;
                    best_state.steiner_locations.push_back(steiner);
                    best_state.obtuse_history.push_back(new_state.obtuse_count);
                    stagnation_count = 0; // Reset stagnation count
                    cout << "New best solution found with better energy: Obtuse Count = " << best_state.obtuse_count 
                    << ", Steiner Points = " << best_state.steiner_points << endl;
                    }
                    // Metropolis criterion
                    if (delta_energy < 0 || exp(-delta_energy / temperature) > dis(gen))
                    {
                        current_state = new_state;
                        current_energy = new_energy;
                       
                    }
                        // Update the best state if the new state is better
                        
                    
                }
            }

            // Remove the triangle from the obtuse_faces list
            obtuse_faces.erase(obtuse_faces.begin() + selected_index);
        }

        

        // Decrease the temperature
        temperature -= 1.0 / 0.95;
        
    }
    
    // Return the best state found
    return best_state;
}

double calculate_energy(const CDT& triangulation, const Point& steiner, double alpha, double beta) {
    CDT temp_cdt = triangulation;
    temp_cdt.insert(steiner);
    int obtuse_count = count_Obtuse_Angles(temp_cdt);
    int steiner_count = temp_cdt.number_of_vertices() - triangulation.number_of_vertices();
    return alpha * obtuse_count + beta * steiner_count;
}
State aco_triangulation(CDT &cdt, const Polygon_2 &convex_hull, int initial_obtuse,
                        double alpha, double beta, int L, int K,
                        double chi, double psi, double lambda)
{
    // Initialize random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    // Initialize pheromone trails
    std::map<Point, double> pheromones;
    for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit)
    {
        pheromones[vit->point()] = 1.0; // Initial pheromone level
    }

    // Initialize best state
    auto start_time = std::chrono::high_resolution_clock::now();
    const int TIME_LIMIT = 60; // 8 minutes in seconds
    State best_state;
    best_state.cdt = cdt;
    best_state.obtuse_count = initial_obtuse;
    best_state.steiner_points = 0;
    best_state.obtuse_history.push_back(initial_obtuse);  // Add initial obtuse count to history

    int stagnation_counter = 0;
    const int MAX_STAGNATION = 50; // Maximum number of iterations without improvement

    // Main ACO loop
    for (int cycle = 0; cycle < L; ++cycle)
    {
        auto current_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time);
        if (duration.count() >= TIME_LIMIT)
        {
            std::cout << "Time limit reached. Stopping aco search." << std::endl;
            break;
        }
        std::vector<State> ant_solutions(K);
        bool improved = false;

        // Each ant constructs a solution
        for (int ant = 0; ant < K; ++ant)
        {
            State current_state = best_state;
            int steps_without_improvement = 0;

            // Ant's solution construction
            for (int step = 0; step < 200 && steps_without_improvement < 50; ++step)
            {
                // Find obtuse triangles
                std::vector<CDT::Face_handle> obtuse_faces;
                for (auto fit = current_state.cdt.finite_faces_begin();
                     fit != current_state.cdt.finite_faces_end(); ++fit)
                {
                    if (is_obtuse_triangle(fit))
                    {
                        obtuse_faces.push_back(fit);
                    }
                }

                if (obtuse_faces.empty())
                    break;

                // Choose a random obtuse triangle
                CDT::Face_handle chosen_face = obtuse_faces[gen() % obtuse_faces.size()];

                // Generate Steiner point candidates
                std::vector<Point> candidates;
                for (int i = 0; i < 5; ++i)
                {
                    Point a = chosen_face->vertex(0)->point();
                    Point b = chosen_face->vertex(1)->point();
                    Point c = chosen_face->vertex(2)->point();
                    Point steiner = select_steiner_point(a, b, c, i, current_state.cdt, convex_hull);
                    if (convex_hull.bounded_side(steiner) != CGAL::ON_UNBOUNDED_SIDE)
                    {
                        candidates.push_back(steiner);
                    }
                }

                if (candidates.empty())
                    continue;

                // Choose Steiner point based on pheromone and heuristic
                std::vector<double> probabilities;
                double sum = 0.0;
                for (const Point &candidate : candidates)
                {
                    double pheromone = pheromones[candidate];
                    double heuristic = 1.0 / (1 + calculate_energy(current_state.cdt, candidate, alpha, beta));
                    double probability = std::pow(pheromone, chi) * std::pow(heuristic, psi);
                    probabilities.push_back(probability);
                    sum += probability;
                }

                double random_value = dis(gen) * sum;
                int chosen_index = 0;
                for (size_t i = 0; i < probabilities.size(); ++i)
                {
                    random_value -= probabilities[i];
                    if (random_value <= 0)
                    {
                        chosen_index = i;
                        break;
                    }
                }

                // Insert chosen Steiner point
                Point chosen_steiner = candidates[chosen_index];
                CDT::Vertex_handle v = current_state.cdt.insert(chosen_steiner);
                if (v != CDT::Vertex_handle())
                {
                    
                    int new_obtuse_count = count_Obtuse_Angles(current_state.cdt);
                    if (new_obtuse_count < current_state.obtuse_count)
                    {
                        current_state.steiner_points++;
                       current_state.steiner_locations.push_back(chosen_steiner);
                        current_state.obtuse_count = new_obtuse_count;
                        steps_without_improvement = 0;
                        best_state.obtuse_history.push_back(current_state.obtuse_count);
                        
                    }
                    else
                    {
                        steps_without_improvement++;
                    }
                }
    }

            

            ant_solutions[ant] = current_state;

            // Update best state if necessary
            if (current_state.obtuse_count < best_state.obtuse_count)
            {
                best_state = current_state;
                improved = true;
                stagnation_counter = 0;
                best_state.obtuse_history.push_back(current_state.obtuse_count);
            }
            if (steps_without_improvement >= 60 && randomization == true) {
                Point random_point = generate_random_point_within_hull(convex_hull);
                cout << "ACO Randomization triggered. Generated random point: (" << random_point.x() << ", " << random_point.y() << ")\n";

                if (convex_hull.bounded_side(random_point) == CGAL::ON_BOUNDED_SIDE ||
                convex_hull.bounded_side(random_point) == CGAL::ON_BOUNDARY)
                {
                 CDT::Vertex_handle v = current_state.cdt.insert(random_point);
                if (v != CDT::Vertex_handle())
                {
                int new_obtuse_count = count_Obtuse_Angles(current_state.cdt);
                best_state.steiner_points++;
                best_state.steiner_locations.push_back(random_point);
                best_state.obtuse_count = new_obtuse_count;
                steps_without_improvement = 0;
                best_state.obtuse_history.push_back(current_state.obtuse_count);
                
                cout << "steiner point counter: " << best_state.steiner_points << "\n";
                cout << "ACO Random point accepted. New obtuse count: " << current_state.obtuse_count << "\n";
                }
                else
                {
                cout << "ACO Random point insertion failed.\n";
                }
            }
            else
            {
            cout << "ACO Random point is outside the convex hull.\n";
            }
        }
        }

        // Update pheromones
        for (auto &p : pheromones)
        {
            p.second *= (1 - lambda); // Evaporation
        }

        for (const State &solution : ant_solutions)
        {
            double delta = 1.0 / (1 + calculate_energy(solution.cdt, Point(0, 0), alpha, beta));
            for (const Point &steiner : solution.steiner_locations)
            {
                pheromones[steiner] += delta;
            }
        }
        if (!best_state.obtuse_history.empty() && 
            best_state.obtuse_history.back() != best_state.obtuse_count)
        {
            
            best_state.obtuse_history.push_back(best_state.obtuse_count);
        }
        

        best_state.obtuse_history.push_back(best_state.obtuse_count);
        if (!improved)
        {
            stagnation_counter++;
            if (stagnation_counter >= MAX_STAGNATION)
            {
                break; // Terminate if no improvement for MAX_STAGNATION iterations
            }
        }
    }

    return best_state;
}

string recognize_input_category(const vector<int> &region_boundary, const vector<pair<int, int>> &additional_constraints, const vector<Point> &points)
{
    bool is_convex = true;
    bool has_open_constraints = false;
    bool has_closed_constraints = false;
    bool is_axis_aligned = true;

    CGAL::Orientation initial_orientation = CGAL::COLLINEAR;

    std::cout << "Debug: Checking convexity and axis alignment..." << std::endl;

    // Check convexity and axis alignment of the region boundary
    for (size_t i = 0; i < region_boundary.size(); ++i)
    {
        size_t j = (i + 1) % region_boundary.size();
        size_t k = (i + 2) % region_boundary.size();
        Point p = points[region_boundary[i]];
        Point q = points[region_boundary[j]];
        Point r = points[region_boundary[k]];

        CGAL::Orientation orientation = CGAL::orientation(p, q, r);

        if (orientation != CGAL::COLLINEAR)
        {
            if (initial_orientation == CGAL::COLLINEAR)
            {
                initial_orientation = orientation;
            }
            else if (orientation != initial_orientation)
            {
                is_convex = false;
                std::cout << "Debug: Non-convex shape detected at points " << i << ", " << j << ", " << k << std::endl;
                break;
            }
        }

        if (p.x() != q.x() && p.y() != q.y())
        {
            is_axis_aligned = false;
            std::cout << "Debug: Non-axis-aligned edge detected between points " << i << " and " << j << std::endl;
        }
    }

    std::cout << "Debug: Is convex: " << (is_convex ? "Yes" : "No") << std::endl;
    std::cout << "Debug: Is axis-aligned: " << (is_axis_aligned ? "Yes" : "No") << std::endl;

    std::cout << "Debug: Checking for closed and open constraints..." << std::endl;

    // Check for closed constraints
    if (!additional_constraints.empty())
    {
        std::cout << "Debug: Number of additional constraints: " << additional_constraints.size() << std::endl;
        
        std::unordered_map<int, std::vector<int>> constraint_graph;
        for (const auto &constraint : additional_constraints)
        {
            constraint_graph[constraint.first].push_back(constraint.second);
            constraint_graph[constraint.second].push_back(constraint.first);
        }

        std::function<bool(int, int, std::vector<bool>&)> dfs = [&](int start, int current, std::vector<bool>& visited) -> bool {
            if (current == start && visited[current]) return true;
            if (visited[current]) return false;
            
            visited[current] = true;
            for (int neighbor : constraint_graph[current]) {
                if (dfs(start, neighbor, visited)) return true;
            }
            return false;
        };

        for (const auto &entry : constraint_graph)
        {
            std::vector<bool> visited(points.size(), false);
            if (dfs(entry.first, entry.first, visited))
            {
                has_closed_constraints = true;
                std::cout << "Debug: Closed constraint detected starting from point " << entry.first << std::endl;
                if (entry.first == 0){
                    has_closed_constraints = false;
                }
                break;
            }
        }

        has_open_constraints = !has_closed_constraints || additional_constraints.size() > constraint_graph.size();
        
        std::cout << "Debug: Has closed constraints: " << (has_closed_constraints ? "Yes" : "No") << std::endl;
        std::cout << "Debug: Has open constraints: " << (has_open_constraints ? "Yes" : "No") << std::endl;
    }
    else
    {
        std::cout << "Debug: No additional constraints" << std::endl;
    }

    // Determine category
    std::string category;
    if (is_convex && additional_constraints.empty())
        category = "A";
    else if (is_convex && has_closed_constraints && !has_open_constraints)
        category = "C";
    else if (is_convex && has_open_constraints)
        category = "B";
    else if (!is_convex && is_axis_aligned && additional_constraints.empty())
        category = "D";
    else
        category = "E";

    std::cout << "Debug: Determined category: " << category << std::endl;
    return category;
}

// Κύρια συνάρτηση
TriangulationResult triangulate(const vector<int> &points_x, const vector<int> &points_y, const vector<int> &region_boundary, const vector<pair<int, int>> &additional_constraints, double alpha, double beta, int L, string &method, bool delaunay)
{
    auto start_time = std::chrono::high_resolution_clock::now();
    // Αρχεικοποιουμε το CDT
    CDT cdt;
    // φτιαχνουμε Vector για να αποθηκευσουμε ολα τα points
    vector<Point> points;
    

    for (size_t i = 0; i < points_x.size(); ++i)
    {
        points.push_back(Point(points_x[i], points_y[i]));
    }

    for (const Point& p : points) {
    cdt.insert(p);
    }
   // δημιουμε το convex hull
    Polygon_2 convex_hull;
    for (size_t i : region_boundary)
    {
        convex_hull.push_back(points[i]);
    }

    // Προσθήκη constraints
    for (size_t i = 0; i < region_boundary.size(); ++i)
    {
        int next = (i + 1) % region_boundary.size();

        cdt.insert_constraint(points[region_boundary[i]], points[region_boundary[next]]);
    }
    // Προσθηκη additional constraints
    for (const auto &constraint : additional_constraints)
    {
        cdt.insert_constraint(points[constraint.first], points[constraint.second]);
    }
    // Μετραμε τα αμβλειγωνια τριγωνα στο CDT
    if (method == "auto")
    {
        string category = recognize_input_category(region_boundary, additional_constraints, points);

        cout << "Input category: " << category << endl;
    }
    int best_obtuse = count_Obtuse_Angles(cdt);
    cout << "Initial obtuse angles: " << best_obtuse << endl;
    // αρχικοποιηση για το best_CDT και state
    CDT best_cdt;
    int max_depth = 6000;
    State best_overall_state;
    best_overall_state.obtuse_count = numeric_limits<int>::max();
    
    if (delaunay)
    {
        // κανουμε την μεθοδο sa
        if (method == "sa")
        {
            // για καλυτερα αποτελσματα την κανουμε 10 φορες και επιστεφουμε την καλυτερη
            for (int i = 0; i < 10; i++){
                cout<< "SA iteration " << i << endl;
                State current_best = sa_triangulation(cdt, convex_hull, best_obtuse, best_cdt, alpha, beta, L);

                if (current_best.obtuse_count < best_overall_state.obtuse_count)
                {
                    best_overall_state = current_best;
                    best_cdt = current_best.cdt;
                }
                if (current_best.obtuse_count == 0){
                    break;
                }
            }
            double total_convergence_rate = 0.0;
            int N = best_overall_state.obtuse_history.size();
            for (int n = 1; n < N; n++) {
                total_convergence_rate += calculate_convergence_rate(n, 
                                                             best_overall_state.obtuse_history[n-1], 
                                                             best_overall_state.obtuse_history[n]);
            }
            // Calculate average convergence rate
            double avg_convergence_rate = 0.0;
            if (N > 1) {
            avg_convergence_rate = total_convergence_rate / (N - 1);
            }
    if (std::isinf(avg_convergence_rate) || std::isnan(avg_convergence_rate)) {
    // Use energy calculation instead
    double initial_energy = alpha * best_overall_state.obtuse_history.front() + beta * 0; // Assuming no initial Steiner points
    double final_energy = alpha * best_overall_state.obtuse_history.back() + beta * best_overall_state.steiner_points;
    avg_convergence_rate = (initial_energy - final_energy) / initial_energy;
    cout << "Using energy-based convergence rate due to invalid average." << endl;
    cout << "final energy: " << final_energy << "\n";
    }else{
    cout << "Average convergence rate " << avg_convergence_rate << endl;
    }
            // καλουμε την μεθοδο local
        cout << "SA New best solution found: Obtuse Count = " << best_overall_state.obtuse_count << "\n";
        cout << "steiner pont counter: " << best_overall_state.steiner_points << "\n"; 
        }
        else if (method == "local")
        {
            
            State initial_state = {cdt, best_obtuse, 0, {}, {}, {count_Obtuse_Angles(cdt)}};
            best_overall_state = local_triangulation(cdt, convex_hull, best_obtuse, best_cdt, max_depth);
            best_cdt = best_overall_state.cdt;
            double total_convergence_rate = 0.0;
    int N = best_overall_state.obtuse_history.size();
    
    for (int n = 1; n < N; n++) {
        total_convergence_rate += calculate_convergence_rate(n, 
                                                             best_overall_state.obtuse_history[n-1], 
                                                             best_overall_state.obtuse_history[n]);
    }
    
    // Calculate average convergence rate
    double avg_convergence_rate = 0.0;
    if (N > 1) {
        avg_convergence_rate = total_convergence_rate / (N - 1);
    }
    if (std::isinf(avg_convergence_rate) || std::isnan(avg_convergence_rate)) {
    // Use energy calculation instead
    double initial_energy = alpha * best_overall_state.obtuse_history.front() + beta * 0; // Assuming no initial Steiner points
    double final_energy = alpha * best_overall_state.obtuse_history.back() + beta * best_overall_state.steiner_points;
    avg_convergence_rate = (initial_energy - final_energy) / initial_energy;
    cout << "Using energy-based convergence rate due to invalid average." << endl;
    }else{
    cout << "Average convergence rate " << avg_convergence_rate << endl;
    }
    cout<< "steiner pont counter: " << best_overall_state.steiner_points << "\n"; 
    cout << "Local New best solution found: Obtuse Count = " << best_overall_state.obtuse_count << "\n";
        
        } else if (method == "aco") {
            State initial_state = {cdt, best_obtuse, 0, {}, {}, {count_Obtuse_Angles(cdt)}};
            // ACO parameters
            int K = 10;  // Number of ants
            double chi = 1.0;  // Pheromone importance
            double psi = 3.0;  // Heuristic importance
            double lambda = 0.5;  // Evaporation rate

            best_overall_state = aco_triangulation(cdt, convex_hull, best_obtuse, 
                                                   alpha, beta, L, K, chi, psi, lambda);
            best_cdt = best_overall_state.cdt;
            double total_convergence_rate = 0.0;
    int N = best_overall_state.obtuse_history.size();
    
    for (int n = 1; n < N; n++) {
        total_convergence_rate += calculate_convergence_rate(n, 
                                                             best_overall_state.obtuse_history[n-1], 
                                                             best_overall_state.obtuse_history[n]);
    }
    cout<< "total convergence rate (p��): " << total_convergence_rate << endl;
    
    // Calculate average convergence rate
    double avg_convergence_rate = 0.0;
    if (N > 1) {
        avg_convergence_rate = total_convergence_rate / (N - 1);
    }
    if (std::isinf(avg_convergence_rate) || std::isnan(avg_convergence_rate)) {
    // Use energy calculation instead
    double initial_energy = alpha * best_overall_state.obtuse_history.front() + beta * 0; // Assuming no initial Steiner points
    double final_energy = alpha * best_overall_state.obtuse_history.back() + beta * best_overall_state.steiner_points;
    avg_convergence_rate = (initial_energy - final_energy) / initial_energy;
    cout << "Using energy-based convergence rate due to invalid average." << endl;
    }else{
    cout << "Average convergence rate " << avg_convergence_rate << endl;
    }
    cout << "steiner pont counter: " << best_overall_state.steiner_points << "\n"; 
    cout << "ACO New best solution found: Obtuse Count = " << best_overall_state.obtuse_count << "\n";
        } else {
        }
    }
    // αν το delaunay ειναι false τοτε τρεχουμε την πρωτη εργασια
    else
    {
        old_triangulate(points_x, points_y, region_boundary, additional_constraints);
    }
    // ετοιμαζουμε τα αποτελεσματα
    
    TriangulationResult results;
    results.obtuse_count = best_overall_state.obtuse_count;

    for (CDT::Finite_edges_iterator eit = best_overall_state.cdt.finite_edges_begin();
         eit != best_overall_state.cdt.finite_edges_end(); ++eit)
    {
        CDT::Face_handle face = eit->first;
        int index = eit->second;
        // παιρνουμε τις δυο κορυφες της ακμης
        CDT::Vertex_handle v1 = face->vertex((index + 1) % 3);
        CDT::Vertex_handle v2 = face->vertex((index + 2) % 3);

        // // βρειτε τις θεσεις αυτων των σημειων
        auto it1 = find(points.begin(), points.end(), v1->point());
        auto it2 = find(points.begin(), points.end(), v2->point());
        // εαν και τα δυο σημεια βρεθουν στην αρχικη εισοδο τα βαζουμε στο αποτελεσμα
        if (it1 != points.end() && it2 != points.end())
        {
            int index1 = distance(points.begin(), it1);
            int index2 = distance(points.begin(), it2);
            results.edges.push_back({index1, index2});
        }
    }
    // Παιρνουμε τα Steiner points
    set<Point> original_points(points.begin(), points.end());
    // ελενχουμε ολα τα τριγωνα και παιρνουμε το καλυτερο
    for (CDT::Finite_vertices_iterator vit = best_overall_state.cdt.finite_vertices_begin();
         vit != best_overall_state.cdt.finite_vertices_end(); ++vit)
    {
        Point p = vit->point();

        // Αν δεν βρισκεται στο original set τότε ειανι steiner point
        if (original_points.find(p) == original_points.end())
        {
            results.steiner_points_x.push_back(CGAL::to_double(p.x()));
            results.steiner_points_y.push_back(CGAL::to_double(p.y()));
        }
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Total execution time: " << duration.count() << " milliseconds" << std::endl;
   
    return results;
    
} 
//////////////////////////////////////////////////////////////////////////