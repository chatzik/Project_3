#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "json.hpp"
#include "triangulation.h"

using json = nlohmann::json;
using namespace std;

// Συνάρτηση για φόρτωση των δεδομένων από το JSON αρχείο
void loadDataFromJSON(const string &filename, vector<int> &points_x, vector<int> &points_y, vector<int> &region_boundary, vector<pair<int, int>> &additional_constraints, string &instance_uid,string& method,double& alpha,double& beta ,int& L,bool& delaunay)
{
    // Άνοιγμα αρχείου JSON
    ifstream inputFile(filename);
    if (!inputFile.is_open())
    {
        cerr << "Error: Could not open the file " << filename << endl;
        return;
    }

    // Φόρτωση JSON δεδομένων
    json j;
    inputFile >> j;

    // Ανάγνωση των arrays από το JSON
    method = j["method"].get<string>();
    points_x = j["points_x"].get<vector<int>>();
    points_y = j["points_y"].get<vector<int>>();
    region_boundary = j["region_boundary"].get<vector<int>>();
    instance_uid = j["instance_uid"].get<string>();
    additional_constraints = j["additional_constraints"].get<vector<pair<int, int>>>();
    // Read parameters from items.json
    auto params = j["parameters"];
    alpha = params["alpha"].get<double>();
    beta = params["beta"].get<double>();
    L = params["L"].get<int>();
    delaunay = j["delaunay"].get<bool>();
    if (j.contains("delaunay")) {
        delaunay = j["delaunay"].get<bool>();
    }
}



int main()
{
    // Δεδομένα που θα φορτωθούν από το JSON αρχείο
    string instance_uid;
    string method;
    vector<int> points_x;
    vector<int> points_y;
    vector<int> region_boundary;
    vector<pair<int, int>> additional_constraints;
    double alpha, beta;
    int L;
    bool delaunay = true;

    // Κάλεσμα της συνάρτησης για φόρτωση δεδομένων
    loadDataFromJSON("data.json", points_x, points_y, region_boundary, additional_constraints, instance_uid, method,alpha ,beta ,L,delaunay);
    cout<< "alpha: " << alpha << ", beta: " << beta << ", L: " << L << ", method: " << method <<  endl;
    // Εκτέλεση τριγωνοποίησης
    cout<<delaunay<<endl;
    TriangulationResult result = triangulate(points_x, points_y, region_boundary, additional_constraints,alpha, beta, L, method,delaunay);
    json output;
    output["content_type"] = "CG_SHOP_2025_Solution";
    output["instance_uid"] = instance_uid;
    output["method"] = method;  
    output["obtuse_count"] = result.obtuse_count;
    output["steiner_points_x"] = result.steiner_points_x;
    output["steiner_points_y"] = result.steiner_points_y;
    json edges_json = json::array();
    for (const auto& edge : result.edges) {
        edges_json.push_back({edge.first, edge.second});
    }
    output["edges"] = edges_json;


    // Write to output.json file
    std::ofstream out_file("output.json");
    if (out_file.is_open()) {
        out_file << std::setw(4) << output << std::endl;
        cout << "Data written to output.json" << endl;
    } else {
        cerr << "Error: Unable to open output.json for writing" << endl;
    }

    // Debug: Print the JSON output
    //cout << "JSON output:" << endl;
   // cout << std::setw(4) << output << endl;

    
    return 0;
}