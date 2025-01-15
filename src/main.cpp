#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "json.hpp"
#include "triangulation.h"
#include <chrono>

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
    
    if (!j.contains("method")){
       method = "auto";
    }else{
        method = j["method"].get<string>();
    }
    points_x = j["points_x"].get<vector<int>>();
    points_y = j["points_y"].get<vector<int>>();
    region_boundary = j["region_boundary"].get<vector<int>>();
    instance_uid = j["instance_uid"].get<string>();
    additional_constraints = j["additional_constraints"].get<vector<pair<int, int>>>();
    // διαβαζουμε τα parameters απο το json
    // Check if parameters exist and read them
if (j.contains("parameters")) {
    auto params = j["parameters"];
    if (params.contains("alpha")) {
        alpha = params["alpha"].get<double>();
    }
    if (params.contains("beta")) {
        beta = params["beta"].get<double>();
    }
    if (params.contains("L")) {
        L = params["L"].get<int>();
    }
} else {
    // Set default values or handle the absence of parameters
    alpha = 5.0;  // Example default value
    beta = 0.5;   // Example default value
    L = 10000;     // Example default value
    cout << "Warning: Parameters not found in JSON. Using default values." << endl;
}
   // Check if delaunay exists
if (j.contains("delaunay")) {
    delaunay = j["delaunay"].get<bool>();
} else {
    delaunay = true;  // Default value if not specified
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
    std::chrono::seconds time_limit(600);
    // Κάλεσμα της συνάρτησης για φόρτωση δεδομένων
    loadDataFromJSON("data.json", points_x, points_y, region_boundary, additional_constraints, instance_uid, method,alpha ,beta ,L,delaunay);

    // Εκτέλεση τριγωνοποίησης
    // Εκτέλεση τριγωνοποίησης
    TriangulationResult result = triangulate(points_x, points_y, region_boundary, additional_constraints,alpha, beta, L, method,delaunay);
    json output;
    output["content_type"] = "CG_SHOP_2025_Solution";
    output["instance_uid"] = instance_uid;
    output["method"] = method;  
    output["obtuse_count"] = result.obtuse_count;
    vector<string> steiner_x, steiner_y;

    for (size_t i = 0; i < result.steiner_points_x.size(); ++i) {
        // Convert x-coordinate to string
        steiner_x.push_back(to_string(static_cast<int>(result.steiner_points_x[i])));
        
        // Convert y-coordinate to fraction string
        double y = result.steiner_points_y[i];
        int y_int = static_cast<int>(y);
        int y_frac = static_cast<int>((y - y_int) * 3 + 0.5);
        if (y_frac == 0) {
            steiner_y.push_back(to_string(y_int));
        } else if (y_frac == 3) {
            steiner_y.push_back(to_string(y_int + 1));
        } else {
            steiner_y.push_back(to_string(y_int) + " " + to_string(y_frac) + "/3");
        }
    }

    output["steiner_points_x"] = steiner_x;
    output["steiner_points_y"] = steiner_y;

    json edges_json = json::array();
    for (const auto& edge : result.edges) {
        edges_json.push_back({edge.first, edge.second});
    }
    output["edges"] = edges_json;

    
    // Γραφουμε τα αποτελεσμένα στο output.json
    ofstream out_file("output.json");
    if (out_file.is_open()) {
        out_file << setw(4) << output << endl;
        cout << "Data written to output.json" << endl;
    } else {
        cerr << "Error: Unable to open output.json for writing" << endl;
    }
    
    return 0;
}