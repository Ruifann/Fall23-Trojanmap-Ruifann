# EE538 Final Project - Fall 2023 - Report

The TrojanMap project has successfully applied graph theory and algorithms to create a functional mapping application. It demonstrates the practical utility of data structures like hash tables, priority queues, and graphs in solving real-world problems. The project covered a range of functionalities, from basic operations like retrieving a location ID to more complex tasks like detecting cycles within a specified area and finding paths that visit all given locations.

### High-Level Overview
 - Project Title: TrojanMap

 - Objective: To implement a mapping application using graph data structures and algorithms.

 - Design: Graph Representation: The map is represented as a graph, with locations as nodes and paths between them as edges.
 - Data Structures: 1. std::unordered_map<std::string, Node> for storing locations and their details. 
 2. std::vector, std::unordered_set, and std::priority_queue for various graph operations.
 - Key Functionalities: Pathfinding, location searching, cycle detection, and proximity queries.
### Detailed Description of Each Function
1. GetID(std::string name)
 - Purpose: Retrieves the ID of a location given its name.
 - Time Complexity: O(1), assuming average-case constant time complexity for hash table lookups.
 - Time Spent: 2 ms
2. CalculateDistance(std::string id1, std::string id2)
 - Purpose: Calculates the distance between two locations.
 - Time Complexity: O(1), direct computation.
 - Time Spent: 5 ms
3. FindNearby(std::string attributesName, std::string name, double r, int k)
 - Purpose: Find locations of a certain class near a given location within a radius r.
 - Time Complexity: O(n), where n is the number of locations.
 - Time Spent: 8 ms
4. CycleDetection(std::vector<double> &square)
 - Purpose: Detects if there is a cycle in a specified square area of the map.
 - Time Complexity: O(V + E) for DFS, where V is the number of vertices and E is the number of edges in the subgraph.
 - Time Spent: 57 ms
5. TrojanPath(std::vectorstd::string &location_names)
 - Purpose: Find the shortest path that visits all specified locations.
 - Time Complexity: O(n^2), assuming a simple heuristic like the nearest neighbor is used.
 - Time Spent: 82ms

### Discussion
Throughout the development of the TrojanMap project, we faced several challenges and made vital observations. One of the primary challenges was implementing efficient graph algorithms to handle various functionalities like pathfinding, cycle detection, and proximity queries. Balancing time complexity and accuracy was another significant aspect of the project, especially in pathfinding and proximity queries where precision and computational efficiency often pull in opposite directions.

For instance, in implementing the FindNearby function, we grappled with the trade-offs between using a precise yet computationally intensive distance calculation (Haversine formula) versus a more straightforward Euclidean approach. Similarly, in CycleDetection, ensuring an accurate representation of cycles within a given area required meticulous graph construction and traversal, balancing efficiency with correctness.

### Conclusion
The key takeaway from this project is the critical importance of selecting appropriate data structures and algorithms for specific tasks. This choice can significantly impact the performance and scalability of an application. The project also highlights the value of understanding algorithmic trade-offs, especially when working with real-world geographical data where precision and efficiency are crucial.

### Several important lessons were learned from this project:
 - The choice of data structures can significantly affect the efficiency of an algorithm. For example, using hash tables for quick lookups significantly improved the performance of certain functionalities.
 - Practical algorithm design often involves trade-offs. Achieving faster execution times might come at the cost of reduced accuracy or vice versa. Understanding and managing these trade-offs is key to effective problem-solving.
 - This project provided valuable insights into how theoretical concepts like graph theory are applied in practical scenarios, such as mapping and navigation tools.
 - The balance between precision (e.g., accurate distance measurements) and computational performance is crucial in geographical applications. Simplifications like the Euclidean distance can be useful in specific contexts despite their limitations.
 - Rigorous testing and debugging are vital, especially in complex systems involving multiple interacting components.
In summary, the TrojanMap project was not only an exercise in applying computer science principles but also a valuable lesson in problem-solving, algorithm selection, and the practicalities of software development.
