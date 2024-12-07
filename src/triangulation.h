#ifndef TRIANGULATION_H
#define TRIANGULATION_H

#include <vector>
#include <utility>
using namespace std;
struct TriangulationResult {
    vector<double> steiner_points_x;
    vector<double> steiner_points_y;
    vector<pair<int, int>> edges;
    int obtuse_count;
};
// Δήλωση της συνάρτησης τριγωνοποίησης
TriangulationResult triangulate(const vector<int> &points_x, const vector<int> &points_y, const vector<int> &region_boundary, const vector<pair<int, int>> &additional_constraints);

#endif // TRIANGULATION_H