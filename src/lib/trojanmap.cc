#include "trojanmap.h"

//-----------------------------------------------------
// TODO: Students should implement the following:
//-----------------------------------------------------
/**
 * GetLat: Get the latitude of a Node given its id. If id does not exist, return
 * -1.
 *
 * @param  {std::string} id : location id
 * @return {double}         : latitude
 */
double TrojanMap::GetLat(const std::string &id) { 
  if (data.find(id) != data.end()){
    return data[id].lat;
  }
  return -1;
}

/**
 * GetLon: Get the longitude of a Node given its id. If id does not exist,
 * return -1.
 *
 * @param  {std::string} id : location id
 * @return {double}         : longitude
 */
double TrojanMap::GetLon(const std::string &id) {
  if(data.find(id) != data.end()){
    return data[id].lon;
  }
  return -1;
}

/**
 * GetName: Get the name of a Node given its id. If id does not exist, return
 * "NULL".
 *
 * @param  {std::string} id : location id
 * @return {std::string}    : name
 */
std::string TrojanMap::GetName(const std::string &id) {
  if(data.find(id) != data.end()){
    return data[id].name;
  }
  return "NULL";
}

/**
 * GetNeighborIDs: Get the neighbor ids of a Node. If id does not exist, return
 * an empty vector.
 *
 * @param  {std::string} id            : location id
 * @return {std::vector<std::string>}  : neighbor ids
 */
std::vector<std::string> TrojanMap::GetNeighborIDs(const std::string &id) {
  if(data.find(id) != data.end()){
    return data[id].neighbors;
  }
  return {};
}

/**
 * GetID: Given a location name, return the id.
 * If the node does not exist, return an empty string.
 * The location name must be unique, which means there is only one node with the name.
 *
 * @param  {std::string} name          : location name
 * @return {std::string}               : id
 */
std::string TrojanMap::GetID(const std::string &name) {
  std::string res = "";
  for (const auto &pair : data){
    if (pair.second.name == name){
      return pair.first;
    }
  }
  return res;
}

/**
 * GetPosition: Given a location name, return the position. If id does not
 * exist, return (-1, -1).
 *
 * @param  {std::string} name          : location name
 * @return {std::pair<double,double>}  : (lat, lon)
 */
std::pair<double, double> TrojanMap::GetPosition(std::string name) {
  std::pair<double, double> results(-1, -1);
  for (const auto &pair : data){
    if (pair.second.name == name){
      return {pair.second.lat, pair.second.lon};
    }
  }
  return results;
}

/**
 * CalculateEditDistance: Calculate edit distance between two location names
 * @param  {std::string} a          : first string
 * @param  {std::string} b          : second string
 * @return {int}                    : edit distance between two strings
 */
int TrojanMap::CalculateEditDistance(std::string a, std::string b) {
  std::vector<std::vector<int>> dp(a.size() + 1, std::vector<int>(b.size() + 1));

  for (unsigned int i = 0; i <= a.size(); i++) {
    for (unsigned int j = 0; j <= b.size(); j++) {
      if (i == 0) {
        dp[i][j] = j;
      } else if (j == 0) {
        dp[i][j] = i;
      } else if (a[i - 1] == b[j - 1]) {
        dp[i][j] = dp[i - 1][j - 1];
      } else {
        dp[i][j] = 1 + std::min({dp[i - 1][j], dp[i][j - 1], dp[i - 1][j - 1]});
      }
    }
  }
  return dp[a.size()][b.size()];
}

/**
 * FindClosestName: Given a location name, return the name with the smallest edit
 * distance.
 *
 * @param  {std::string} name          : location name
 * @return {std::string} tmp           : the closest name
 */
std::string TrojanMap::FindClosestName(std::string name) {
  std::string closestName = "";
  int minDistance = INT_MAX;

  for (const auto &pair : data) {
    int distance = CalculateEditDistance(name, pair.second.name);
    if (distance < minDistance) {
      minDistance = distance;
      closestName = pair.second.name;
    }
  }
  return closestName;
}


/**
 * Autocomplete: Given a parital name return all the possible locations with
 * partial name as the prefix. The function should be case-insensitive.
 *
 * @param  {std::string} name          : partial name
 * @return {std::vector<std::string>}  : a vector of full names
 */
std::vector<std::string> TrojanMap::Autocomplete(std::string name) {
    std::vector<std::string> results;
    std::transform(name.begin(), name.end(), name.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    for (const auto& pair : data) {
        std::string nodeName = pair.second.name;
        std::transform(nodeName.begin(), nodeName.end(), nodeName.begin(),
                       [](unsigned char c){ return std::tolower(c); });
        if (nodeName.find(name) == 0) {
            results.push_back(pair.second.name);
        }
    }
    return results;
}

/**
 * GetAllCategories: Return all the possible unique location categories, i.e.
 * there should be no duplicates in the output.
 *
 * @return {std::vector<std::string>}  : all unique location categories
 */
std::vector<std::string> TrojanMap::GetAllCategories() {
  std::unordered_set<std::string> categorySet;
  for(const auto& pair : data){
    for (const auto& category : pair.second.attributes){
      categorySet.insert(category);
      }
  }
  std::vector<std::string> categories(categorySet.begin(), categorySet.end());
    return categories;
}

/**
 * GetAllLocationsFromCategory: Return all the locations of the input category (i.e.
 * 'attributes' in data.csv). If there is no location of that category, return
 * (-1, -1). The function should be case-insensitive.
 *
 * @param  {std::string} category         : category name (attribute)
 * @return {std::vector<std::string>}     : ids
 */
std::vector<std::string> TrojanMap::GetAllLocationsFromCategory(std::string category) {
    std::vector<std::string> res;
    bool categoryFound = false;

    // Convert the input category to lowercase for case-insensitive comparison
    std::transform(category.begin(), category.end(), category.begin(),
                   [](unsigned char c){ return std::tolower(c); });

    for (const auto& pair : data) {
        for (const auto& attr : pair.second.attributes) {
            std::string lowerAttr = attr;

            // Convert the attribute to lowercase
            std::transform(lowerAttr.begin(), lowerAttr.end(), lowerAttr.begin(),
                           [](unsigned char c){ return std::tolower(c); });

            if (lowerAttr == category) {
                res.push_back(pair.first); // Add the ID to the result
                categoryFound = true;
            }
        }
    }

    if (!categoryFound) {
        res.push_back("-1, -1"); // Add (-1, -1) if no location of that category
    }

    return res;
}


/**
 * GetLocationRegex: Given the regular expression of a location's name, your
 * program should first check whether the regular expression is valid, and if so
 * it returns all locations that match that regular expression.
 *
 * @param  {std::regex} location name      : the regular expression of location
 * names
 * @return {std::vector<std::string>}     : ids
 */
std::vector<std::string> TrojanMap::GetLocationRegex(std::regex location) {
  std::vector<std::string> matchingIDs;
  for (const auto& pair : data){
    if (std::regex_match(pair.second.name, location)){
      matchingIDs.push_back(pair.first);
    }
  }
  return matchingIDs;
}

/**
 * CalculateDistance: Get the distance between 2 nodes.
 * We have provided the code for you. Please do not need to change this function.
 * You can use this function to calculate the distance between 2 nodes.
 * The distance is in mile.
 * The distance is calculated using the Haversine formula.
 * https://en.wikipedia.org/wiki/Haversine_formula
 * 
 * @param  {std::string} a  : a_id
 * @param  {std::string} b  : b_id
 * @return {double}  : distance in mile
 */
double TrojanMap::CalculateDistance(const std::string &a_id,
                                    const std::string &b_id) {
  // Do not change this function
  Node a = data[a_id];
  Node b = data[b_id];
  double dlon = (b.lon - a.lon) * M_PI / 180.0;
  double dlat = (b.lat - a.lat) * M_PI / 180.0;
  double p = pow(sin(dlat / 2), 2.0) + cos(a.lat * M_PI / 180.0) *
                                           cos(b.lat * M_PI / 180.0) *
                                           pow(sin(dlon / 2), 2.0);
  double c = 2 * asin(std::min(1.0, sqrt(p)));
  return c * 3961;
}

/**
 * CalculatePathLength: Calculates the total path length for the locations
 * inside the vector.
 * We have provided the code for you. Please do not need to change this function.
 * 
 * @param  {std::vector<std::string>} path : path
 * @return {double}                        : path length
 */
double TrojanMap::CalculatePathLength(const std::vector<std::string> &path) {
  // Do not change this function
  double sum = 0;
  for (int i = 0; i < int(path.size()) - 1; i++) {
    sum += CalculateDistance(path[i], path[i + 1]);
  }
  return sum;
}

/**
 * CalculateShortestPath_Dijkstra: Given 2 locations, return the shortest path
 * which is a list of id. Hint: Use priority queue.
 *
 * @param  {std::string} location1_name     : start
 * @param  {std::string} location2_name     : goal
 * @return {std::vector<std::string>}       : path
 */
std::vector<std::string> TrojanMap::CalculateShortestPath_Dijkstra(
    std::string location1_name, std::string location2_name) {

    auto location1_id = GetID(location1_name);
    auto location2_id = GetID(location2_name);

    std::unordered_map<std::string, double> distances;
    std::unordered_map<std::string, std::string> previous;
    std::priority_queue<std::pair<double, std::string>,
                        std::vector<std::pair<double, std::string>>,
                        std::greater<std::pair<double, std::string>>> pq;

    // Initialize distances
    for (const auto& pair : data) {
        distances[pair.first] = std::numeric_limits<double>::infinity();
    }

    // Start with the first location
    distances[location1_id] = 0.0;
    pq.push({0.0, location1_id});

    while (!pq.empty()) {
        std::string current_id = pq.top().second;
        double current_distance = pq.top().first;
        pq.pop();

        // Early exit if we reached the destination
        if (current_id == location2_id) break;

        for (const auto& neighbor_id : data[current_id].neighbors) {
            double distance = current_distance + CalculateDistance(current_id, neighbor_id);
            if (distance < distances[neighbor_id]) {
                distances[neighbor_id] = distance;
                previous[neighbor_id] = current_id;
                pq.push({distance, neighbor_id});
            }
        }
    }
    // Reconstruct the shortest path
    std::vector<std::string> path;
    for (std::string at = location2_id; at != location1_id; at = previous[at]) {
        path.push_back(at);
        if (previous.find(at) == previous.end()) {
            return {}; // Path not found
        }
    }
    path.push_back(location1_id);
    std::reverse(path.begin(), path.end());

    return path;
}

/**
 * CalculateShortestPath_Bellman_Ford: Given 2 locations, return the shortest
 * path which is a list of id. Hint: Do the early termination when there is no
 * change on distance.
 *
 * @param  {std::string} location1_name     : start
 * @param  {std::string} location2_name     : goal
 * @return {std::vector<std::string>}       : path
 */

std::vector<std::string> TrojanMap::CalculateShortestPath_Bellman_Ford(
    std::string location1_name, std::string location2_name) {

    std::unordered_map<std::string, double> distances;
    std::unordered_map<std::string, std::string> predecessors;
    auto location1_id = GetID(location1_name);
    auto location2_id = GetID(location2_name);

    // Initialize distances
    for (const auto& pair : data) {
        distances[pair.first] = std::numeric_limits<double>::infinity();
    }
    distances[location1_id] = 0;

    // Relaxation step
    for (size_t i = 0; i < data.size() - 1; ++i) {
        bool updated = false;
        for (const auto& pair : data) {
            for (const auto& neighbor_id : pair.second.neighbors) {
                double distance = distances[pair.first] + CalculateDistance(pair.first, neighbor_id);
                if (distance < distances[neighbor_id]) {
                    distances[neighbor_id] = distance;
                    predecessors[neighbor_id] = pair.first;
                    updated = true;
                }
            }
        }
        if (!updated) {
            break; // Early termination if no updates
        }
    }

    // Reconstruct the shortest path
    std::vector<std::string> path;
    for (std::string at = location2_id; at != location1_id; at = predecessors[at]) {
        path.push_back(at);
        if (predecessors.find(at) == predecessors.end()) {
            return {}; // Path not found
        }
    }
    path.push_back(location1_id);
    std::reverse(path.begin(), path.end());

    return path;
}

/**
 * Traveling salesman problem: Given a list of locations, return the shortest
 * path which visit all the places and back to the start point.
 *
 * @param  {std::vector<std::string>} input : a list of locations needs to visit
 * @return {std::pair<double, std::vector<std::vector<std::string>>} : a pair of total distance and the all the progress to get final path, 
 *                                                                      for example: {10.3, {{0, 1, 2, 3, 4, 0}, {0, 1, 2, 3, 4, 0}, {0, 4, 3, 2, 1, 0}}},
 *                                                                      where 10.3 is the total distance, 
 *                                                                      and the first vector is the path from 0 and travse all the nodes and back to 0,
 *                                                                      and the second vector is the path shorter than the first one,
 *                                                                      and the last vector is the shortest path.
 */
// Please use brute force to implement this function, ie. find all the permutations.
std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravelingTrojan_Brute_force(
  std::vector<std::string> location_ids) {
    std::pair<double, std::vector<std::vector<std::string>>> records;
    double min_distance = std::numeric_limits<double>::max();

    // Start permutation from the second element to keep the first location fixed as starting point
    std::sort(location_ids.begin() + 1, location_ids.end());

    do {
        double current_distance = 0;
        std::vector<std::string> current_path = location_ids;
        current_path.push_back(location_ids[0]); // Return to start

        // Calculate the distance of the current path
        for (int i = 0; i < int(current_path.size()) - 1; i++) {
            current_distance += CalculateDistance(current_path[i], current_path[i + 1]);
        }

        // Update records if a new shorter path is found
        if (current_distance < min_distance) {
            min_distance = current_distance;
            records.first = min_distance;
            records.second.push_back(current_path); // Add the current path
        }

    } while (std::next_permutation(location_ids.begin() + 1, location_ids.end())); // Generate next permutation

    return records;
}


// Please use backtracking to implement this function
// We need a Helper function for backtracking
void TSP_Backtrack(int index, double current_distance, double &min_distance, 
                   std::vector<std::string> &current_path, std::vector<bool> &visited, 
                   const std::vector<std::string> &location_ids, TrojanMap &map, 
                   std::vector<std::vector<std::string>> &all_paths) {

    if (index == location_ids.size() && 
        current_distance + map.CalculateDistance(current_path.back(), current_path.front()) < min_distance) {
        // Completing the cycle by returning to the starting point
        min_distance = current_distance + map.CalculateDistance(current_path.back(), current_path.front());
        current_path.push_back(current_path.front()); // Add start to the end to complete the cycle
        all_paths.push_back(current_path);
        current_path.pop_back(); // Remove start from the end to allow further exploration
        return;
    }

    for (int i = 0; i < location_ids.size(); ++i) {
        if (!visited[i]) {
            visited[i] = true;
            current_path.push_back(location_ids[i]);
            
            double new_distance = current_distance;
            if (index > 0) {
                new_distance += map.CalculateDistance(current_path[index - 1], current_path[index]);
            }

            if (new_distance < min_distance) {
                TSP_Backtrack(index + 1, new_distance, min_distance, current_path, visited, location_ids, map, all_paths);
            }

            visited[i] = false;
            current_path.pop_back();
        }
    }
}

std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravelingTrojan_Backtracking(
                                    std::vector<std::string> location_ids) {
    std::pair<double, std::vector<std::vector<std::string>>> records;
    records.first = std::numeric_limits<double>::max();
    std::vector<std::string> current_path;
    std::vector<bool> visited(location_ids.size(), false);
    std::vector<std::vector<std::string>> all_paths;

    TSP_Backtrack(0, 0.0, records.first, current_path, visited, location_ids, *this, all_paths);

    records.second = all_paths;
    return records;
}

// Hint: https://en.wikipedia.org/wiki/2-opt
std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravelingTrojan_2opt(
      std::vector<std::string> location_ids) {
    std::pair<double, std::vector<std::vector<std::string>>> records;

    // Initialize with a random or specific path
    std::vector<std::string> current_path = location_ids;
    current_path.push_back(location_ids.front()); // Make it a cycle by returning to the start
    double min_distance = CalculatePathLength(current_path);
    records.second.push_back(current_path);

    bool improvement = true;
    while (improvement) {
        improvement = false;
        for (int i = 1; i < current_path.size() - 2; i++) {
            for (int j = i + 1; j < current_path.size() - 1; j++) {
                // Try to swap edges and check if it results in a shorter path
                std::vector<std::string> new_path = current_path;
                std::reverse(new_path.begin() + i, new_path.begin() + j + 1); // Reverse the segment between i and j

                double new_distance = CalculatePathLength(new_path);
                if (new_distance < min_distance) {
                    // Found a shorter path
                    min_distance = new_distance;
                    current_path = new_path;
                    records.second.push_back(new_path);
                    improvement = true;
                }
            }
        }
    }

    records.first = min_distance;
    return records;
}


// This is optional
std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravelingTrojan_3opt(
      std::vector<std::string> location_ids){
  std::pair<double, std::vector<std::vector<std::string>>> records;
  return records;
}

/**
 * Given CSV filename, it read and parse locations data from CSV file,
 * and return locations vector for topological sort problem.
 * We have provided the code for you. Please do not need to change this function.
 * Example: 
 *   Input: "topologicalsort_locations.csv"
 *   File content:
 *    Name
 *    Ralphs
 *    KFC
 *    Chick-fil-A
 *   Output: ['Ralphs', 'KFC', 'Chick-fil-A']
 * @param  {std::string} locations_filename     : locations_filename
 * @return {std::vector<std::string>}           : locations
 */
std::vector<std::string> TrojanMap::ReadLocationsFromCSVFile(
    std::string locations_filename) {
  std::vector<std::string> location_names_from_csv;
  std::fstream fin;
  fin.open(locations_filename, std::ios::in);
  std::string line, word;
  getline(fin, line);
  while (getline(fin, word)) {
    location_names_from_csv.push_back(word);
  }
  fin.close();
  return location_names_from_csv;
}

/**
 * Given CSV filenames, it read and parse dependencise data from CSV file,
 * and return dependencies vector for topological sort problem.
 * We have provided the code for you. Please do not need to change this function.
 * Example: 
 *   Input: "topologicalsort_dependencies.csv"
 *   File content:
 *     Source,Destination
 *     Ralphs,Chick-fil-A
 *     Ralphs,KFC
 *     Chick-fil-A,KFC
 *   Output: [['Ralphs', 'Chick-fil-A'], ['Ralphs', 'KFC'], ['Chick-fil-A', 'KFC']]
 * @param  {std::string} dependencies_filename     : dependencies_filename
 * @return {std::vector<std::vector<std::string>>} : dependencies
 */
std::vector<std::vector<std::string>> TrojanMap::ReadDependenciesFromCSVFile(
    std::string dependencies_filename) {
  std::vector<std::vector<std::string>> dependencies_from_csv;
  std::fstream fin;
  fin.open(dependencies_filename, std::ios::in);
  std::string line, word;
  getline(fin, line);
  while (getline(fin, line)) {
    std::stringstream s(line);
    std::vector<std::string> dependency;
    while (getline(s, word, ',')) {
      dependency.push_back(word);
    }
    dependencies_from_csv.push_back(dependency);
  }
  fin.close();
  return dependencies_from_csv;
}

/**
 * DeliveringTrojan: Given a vector of location names, it should return a
 * sorting of nodes that satisfies the given dependencies. If there is no way to
 * do it, return a empty vector.
 *
 * @param  {std::vector<std::string>} locations                     : locations
 * @param  {std::vector<std::vector<std::string>>} dependencies     : prerequisites
 * @return {std::vector<std::string>} results                       : results
 */
std::vector<std::string> TrojanMap::DeliveringTrojan(
    std::vector<std::string> &locations,
    std::vector<std::vector<std::string>> &dependencies) {
    std::vector<std::string> result;
    std::unordered_map<std::string, int> in_degree;
    std::unordered_map<std::string, std::vector<std::string>> graph;
    std::queue<std::string> process_queue;

    // Initialize in-degree of all locations
    for (const auto &loc : locations) {
        in_degree[loc] = 0;
    }

    // Create graph and calculate in-degrees
    for (const auto &dep : dependencies) {
        graph[dep[0]].push_back(dep[1]);
        in_degree[dep[1]]++;
    }

    // Add locations with no incoming edges to the queue
    for (const auto &pair : in_degree) {
        if (pair.second == 0) {
            process_queue.push(pair.first);
        }
    }

    // Perform topological sort
    while (!process_queue.empty()) {
        std::string u = process_queue.front();
        process_queue.pop();
        result.push_back(u);

        for (const auto &v : graph[u]) {
            in_degree[v]--;
            if (in_degree[v] == 0) {
                process_queue.push(v);
            }
        }
    }

    // Check if topological sort is possible
    if (result.size() != locations.size()) {
        return {}; // Return empty vector if a cycle is detected or dependencies are incomplete
    }

    return result;
}


/**
 * inSquare: Give a id retunr whether it is in square or not.
 *
 * @param  {std::string} id            : location id
 * @param  {std::vector<double>} square: four vertexes of the square area
 * @return {bool}                      : in square or not
 */
bool TrojanMap::inSquare(std::string id, std::vector<double> &square) {
    if (data.find(id) == data.end()) {
        return false; // Location id not found in the data
    }

    Node location = data[id];
    double lat = location.lat;
    double lon = location.lon;
    double min_lon = square[0]; // Left bound longitude
    double max_lon = square[1]; // Right bound longitude
    double max_lat = square[2]; // Upper bound latitude
    double min_lat = square[3]; // Lower bound latitude

    // Check if the location's coordinates are within the bounds
    return (lat >= min_lat && lat <= max_lat && lon >= min_lon && lon <= max_lon);
}


/**
 * GetSubgraph: Give four vertexes of the square area, return a list of location
 * ids in the squares
 *
 * @param  {std::vector<double>} square         : four vertexes of the square
 * area
 * @return {std::vector<std::string>} subgraph  : list of location ids in the
 * square
 */
//std::vector<std::string> TrojanMap::GetSubgraph(std::vector<double> &square) {
  // include all the nodes in subgraph
  //std::vector<std::string> subgraph;
  //return subgraph;
//}

std::vector<std::string> TrojanMap::GetSubgraph(std::vector<double> &square) {
    std::vector<std::string> subgraph;
    double min_lon = square[0];
    double max_lon = square[1];
    double max_lat = square[2];
    double min_lat = square[3];

    for (const auto &pair : data) {
        double lat = pair.second.lat;
        double lon = pair.second.lon;

        // Check if the location's coordinates are within the bounds
        if (lat >= min_lat && lat <= max_lat && lon >= min_lon && lon <= max_lon) {
            subgraph.push_back(pair.first);
        }
    }

    return subgraph;
}

/**
 * Cycle Detection: Given four points of the square-shape subgraph, return true
 * if there is a cycle path inside the square, false otherwise.
 *
 * @param {std::vector<std::string>} subgraph: list of location ids in the
 * square
 * @param {std::vector<double>} square: four vertexes of the square area
 * @return {bool}: whether there is a cycle or not
 */
// Helper function for DFS
bool DFS_Cycle(std::string current, std::string parent, std::unordered_set<std::string> &visited, 
               std::unordered_map<std::string, std::vector<std::string>> &graph) {
    visited.insert(current);

    for (const auto &neighbor : graph[current]) {
        if (visited.find(neighbor) == visited.end()) {
            if (DFS_Cycle(neighbor, current, visited, graph)) {
                return true;
            }
        } else if (neighbor != parent) {
            // Found a cycle
            return true;
        }
    }
    return false;
}

bool TrojanMap::CycleDetection(std::vector<std::string> &subgraph, std::vector<double> &square) {
    // Filter nodes that are actually inside the square
    std::unordered_set<std::string> nodes_in_square;
    for (const auto &id : subgraph) {
        if (inSquare(id, square)) {
            nodes_in_square.insert(id);
        }
    }

    // Create a graph from the nodes inside the square
    std::unordered_map<std::string, std::vector<std::string>> graph;
    for (const auto &id : nodes_in_square) {
        for (const auto &neighbor_id : data[id].neighbors) {
            if (nodes_in_square.find(neighbor_id) != nodes_in_square.end()) {
                graph[id].push_back(neighbor_id);
            }
        }
    }

    std::unordered_set<std::string> visited;
    for (const auto &node : nodes_in_square) {
        if (visited.find(node) == visited.end()) {
            if (DFS_Cycle(node, "", visited, graph)) {
                return true; // Cycle found
            }
        }
    }

    return false; // No cycle found
}

/**
 * FindNearby: Given a class name C, a location name L and a number r,
 * find all locations in class C on the map near L with the range of r and
 * return a vector of string ids
 *
 * @param {std::string} className: the name of the class
 * @param {std::string} locationName: the name of the location
 * @param {double} r: search radius
 * @param {int} k: search numbers
 * @return {std::vector<std::string>}: location name that meets the requirements
 */
std::vector<std::string> TrojanMap::FindNearby(std::string attributesName, std::string locationName, double r, int k) {
    std::vector<std::string> res;
    std::string reference_id = GetID(locationName);
    if (reference_id.empty() || data.find(reference_id) == data.end()) {
        return res; // Location name not found or ID is empty
    }

    Node reference_location = data[reference_id];
    std::vector<std::pair<double, std::string>> filtered_locations;

    for (const auto &pair : data) {
        if (pair.first != reference_id && pair.second.attributes.find(attributesName) != pair.second.attributes.end()) {
            double distance = CalculateDistance(reference_id, pair.first);
            if (distance <= r) {
                filtered_locations.push_back(std::make_pair(distance, pair.first));
            }
        }
    }

    // Sort by distance
    std::sort(filtered_locations.begin(), filtered_locations.end());

    // Collect up to k results
    for (int i = 0; i < std::min(k, static_cast<int>(filtered_locations.size())); i++) {
        res.push_back(filtered_locations[i].second);
    }

    return res;
}


/**
 * Shortest Path to Visit All Nodes: Given a list of locations, return the shortest
 * path which visit all the places and no need to go back to the start point.
 *
 * @param  {std::vector<std::string>} input : a list of locations needs to visit
 * @return {std::vector<std::string> }      : the shortest path
 */
std::vector<std::string> TrojanMap::TrojanPath(std::vector<std::string> &location_names) {
    std::vector<std::string> location_ids;
    for (const auto& name : location_names) {
        std::string id = GetID(name);
        if (!id.empty()) {
            location_ids.push_back(id);
        }
    }

    // Use the backtracking approach to find the shortest path
    auto result = TravelingTrojan_Backtracking(location_ids);

    // Extract the shortest path from the result
    std::vector<std::string> shortest_path;
    if (!result.second.empty()) {
        shortest_path = result.second.back(); // Assuming the last path in all_paths is the shortest
    }

    return shortest_path;
}


/**
 * Given a vector of queries, find whether there is a path between the two locations with the constraint of the gas tank.
 *
 * @param  {std::vector<std::pair<double, std::vector<std::string>>>} Q : a list of queries
 * @return {std::vector<bool> }      : existence of the path
 */
std::vector<bool> TrojanMap::Queries(const std::vector<std::pair<double, std::vector<std::string>>>& q) {
    std::vector<bool> ans(q.size());

    for (int i = 0; i < q.size(); ++i) {
        double tank = q[i].first;
        std::string start_name = q[i].second[0];
        std::string end_name = q[i].second[1];

        // Check if start and end locations exist
        if (data.find(GetID(start_name)) == data.end() || data.find(GetID(end_name)) == data.end()) {
            ans[i] = false;
            continue;
        }

        std::string start_id = GetID(start_name);
        std::string end_id = GetID(end_name);

        // BFS to find a path with the gas tank constraint
        std::queue<std::string> bfs_queue;
        std::unordered_set<std::string> visited;
        bfs_queue.push(start_id);
        visited.insert(start_id);
        bool path_found = false;

        while (!bfs_queue.empty() && !path_found) {
            std::string current = bfs_queue.front();
            bfs_queue.pop();

            if (current == end_id) {
                path_found = true;
                break;
            }

            for (const auto& neighbor_id : data[current].neighbors) {
                if (visited.find(neighbor_id) == visited.end()) {
                    double distance = CalculateDistance(current, neighbor_id);
                    if (distance <= tank) {
                        bfs_queue.push(neighbor_id);
                        visited.insert(neighbor_id);
                    }
                }
            }
        }

        ans[i] = path_found;
    }

    return ans;
}


/**
 * CreateGraphFromCSVFile: Read the map data from the csv file
 * We have provided the code for you. Please do not need to change this function.
 */
void TrojanMap::CreateGraphFromCSVFile() {
  // Do not change this function
  std::fstream fin;
  fin.open("src/lib/data.csv", std::ios::in);
  std::string line, word;

  getline(fin, line);
  while (getline(fin, line)) {
    std::stringstream s(line);

    Node n;
    int count = 0;
    while (getline(s, word, ',')) {
      word.erase(std::remove(word.begin(), word.end(), '\''), word.end());
      word.erase(std::remove(word.begin(), word.end(), '"'), word.end());
      word.erase(std::remove(word.begin(), word.end(), '{'), word.end());
      word.erase(std::remove(word.begin(), word.end(), '}'), word.end());
      if (count == 0)
        n.id = word;
      else if (count == 1)
        n.lat = stod(word);
      else if (count == 2)
        n.lon = stod(word);
      else if (count == 3)
        n.name = word;
      else {
        word.erase(std::remove(word.begin(), word.end(), ' '), word.end());
        if (isalpha(word[0])) n.attributes.insert(word);
        if (isdigit(word[0])) n.neighbors.push_back(word);
      }
      count++;
    }
    data[n.id] = n;
  }
  fin.close();
}
