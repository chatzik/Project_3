#ifndef TRIANGULATION_H
#define TRIANGULATION_H

#include <vector>
#include <utility>
using namespace std;
// Δήλωση της συνάρτησης τριγωνοποίησης
void triangulate(const vector<int> &points_x, const vector<int> &points_y, const vector<int> &region_boundary, const vector<pair<int, int>> &additional_constraints);

#endif // TRIANGULATION_H