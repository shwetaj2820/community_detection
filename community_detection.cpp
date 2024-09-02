#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <omp.h>
#include <algorithm>

class Graph {
public:
    std::map<int, std::vector<std::pair<int, int>>> adjList; // adjacency list: node -> (neighbor, weight)
    std::map<int, int> degrees; // degree of each node
    int num_nodes;

    Graph() : num_nodes(0) {}

    void addNode(int node) {
        if (adjList.find(node) == adjList.end()) {
            adjList[node] = std::vector<std::pair<int, int>>(); //initalize adj list for new node with empty list of edges
            degrees[node] = 0; //initial node degree = 0
            num_nodes++;
        }
    }

// remove node (along with its connected edges)
    void removeNode(int node) {
        if (adjList.find(node) == adjList.end()) return; //check for node exist or not

        // Remove edges connected to this node by: the iteration to all the neighbors; updates the degrees 
        // this update done by minusing the weight of the edge
        for (const auto& [neighbor, weight] : adjList[node]) {
            degrees[neighbor] -= weight;
        }
        adjList.erase(node); //remove node
        degrees.erase(node); //erase degree of the node
        num_nodes--; 
    }

    void addEdge(int u, int v, int weight) {
        if (u == v) return; // Avoidment of self-loops

        addNode(u);
        addNode(v);

        // Add edge u -> v
        adjList[u].emplace_back(v, weight);
        degrees[u] += weight;

        // Add edge v -> u 
        adjList[v].emplace_back(u, weight);
        degrees[v] += weight;
    }

    // void removeEdge(int u, int v) {



    // }

    int degree(int u) const {
        auto it = degrees.find(u);
        if (it != degrees.end()) {
            return it->second;
        }
        return 0;
    }

    const std::vector<std::pair<int, int>>& neighbors(int u) const {
        auto it = adjList.find(u);
        if (it != adjList.end()) {
            return it->second;
        }
        return {};
    }

    double modularity(const std::map<int, int>& community, int m) const {
        double Q = 0.0;

        // TODO: can paralyze this loop:
        for (const auto& [u, neighbors] : adjList) {
            for (const auto& [v, weight] : neighbors) {
                if (community.at(u) == community.at(v)) {
                    Q += weight - (degree(u) * degree(v)) / (2.0 * m);
                }
            }
        }
        return Q / (2.0 * m);
    }

    double modularityGain(int u, int community_c, const std::map<int, int>& community, int m, const std::map<int, int>& community_totals) const {
        double gain = 0.0;
        int sum_in = 0;
        int sum_tot = community_totals.at(community_c);

        for (const auto& [v, weight] : neighbors(u)) {
            if (community.at(v) == community_c) {
                sum_in += weight;
            }
        }

        int k_i = degree(u);

        gain = (sum_in + k_i) / (2.0 * m) - std::pow(sum_tot + k_i, 2) / (2.0 * m * m)
               - sum_in / (2.0 * m) + std::pow(sum_tot, 2) / (2.0 * m * m) + std::pow(k_i, 2) / (2.0 * m * m);

        return gain;
    }

    void louvainAlgorithm(std::map<int, int>& community) {
        int m = 0; 

        // finding total weight of edges
        for (const auto& [node, neighbors] : adjList) {
            for (const auto& [neighbor, weight] : neighbors) {
                m += weight;
            }
        }
        m /= 2;

        std::map<int, int> community_totals;

        // iterate on nodes to add to community_totals
        for (const auto& [node, neighbors] : adjList) {
            if (degree(node) > 0) {
                community_totals[community[node]] += degree(node);
            }
        }

        bool improvement = true;

        // iteration until no improvement is found
        while (improvement) {
            improvement = false;
            for (const auto& [u, neighbors] : adjList) {
                int current_comm = community[u];
                double best_gain = 0.0;
                int best_community = current_comm;

                // gain by moving community u to community v:
                for (const auto& [v, weight] : neighbors) {
                    if (community[u] != community[v]) {
                        int neighbor_comm = community[v];
                        double gain = modularityGain(u, neighbor_comm, community, m, community_totals);
                        if (gain > best_gain) {
                            best_gain = gain;
                            best_community = neighbor_comm;
                        }
                    }
                }

                // updating community of node "u" if improvement happens:
                if (best_community != current_comm) {
                    community[u] = best_community;
                    improvement = true;
                    community_totals[current_comm] -= degree(u); //updating degrees
                    community_totals[best_community] += degree(u);
                }
            }
        }
    }
};

int main() {
    Graph g;

    // TODO: input large graph file
    g.addEdge(0, 1, 1);
    g.addEdge(0, 2, 2);
    g.addEdge(1, 2, 3);
    g.addEdge(2, 4, 8);
    g.addEdge(3, 4, 1);
    g.addEdge(4, 5, 3);
    g.addEdge(3, 5, 2);

    // Initialized communities
    std::map<int, int> community;
    for (const auto& [node, _] : g.adjList) {
        community[node] = node; // Each node starts in its own community
    }

    // Applay Louvain algorithm
    g.louvainAlgorithm(community);

    // Output results
    std::map<int, std::set<int>> communities;
    for (const auto& [node, comm] : community) {
        communities[comm].insert(node);
    }

    std::cout << "Communities:\n";
    for (const auto& [comm, nodes] : communities) {
        std::cout << "Community " << comm << ": ";
        for (int node : nodes) {
            std::cout << node << " ";
        }
        std::cout << "\n";
    }

// TODO: plan something to make add/remove in a loop directly so it becomes dynamic;
// current implementation: example for add/remove; same logic to be implemented if making dynamic (first need the graph)

    g.addEdge(6, 7, 5);  // Add new edge
    // g.removeEdge(0, 1); //Remove existing edge

    // Recompute communities after modifications
    g.louvainAlgorithm(community);

    // Output the updated results
    communities.clear();
    for (const auto& [node, comm] : community) {
        communities[comm].insert(node);
    }

    std::cout << "Updated Communities:\n";
    for (const auto& [comm, nodes] : communities) {
        std::cout << "Community " << comm << ": ";
        for (int node : nodes) {
            std::cout << node << " ";
        }
        std::cout << "\n";
    }

    return 0;
}

// TODO : calculation of modularity & communities(before --> after : add/remove edges and nodes)
