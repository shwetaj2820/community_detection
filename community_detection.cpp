#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <omp.h>

// Simple representation of a graph
class Graph {
public:
    std::vector<std::vector<int>> adjMatrix;
    std::vector<int> degrees;
    int num_nodes;

    Graph(int n) : num_nodes(n) {
        adjMatrix.resize(n, std::vector<int>(n, 0));
        degrees.resize(n, 0);
    }

    void addEdge(int u, int v, int weight) {
        adjMatrix[u][v] = weight;
        adjMatrix[v][u] = weight;
        degrees[u] += weight;
        degrees[v] += weight;
    }

    int degree(int u) const {
        return degrees[u];
    }

    const std::vector<int>& neighbors(int u) const {
        return adjMatrix[u];
    }
};

// Function to calculation modularity
double modularity(const Graph& g, const std::map<int, int>& community, int m) {
    double Q = 0.0;
    #pragma omp parallel for reduction(+:Q)
    for (int u = 0; u < g.num_nodes; ++u) {
        for (int v = 0; v < g.num_nodes; ++v) {
            if (g.adjMatrix[u][v] > 0 && community.at(u) == community.at(v)) {
                Q += g.adjMatrix[u][v] - (g.degree(u) * g.degree(v)) / (2.0 * m);
            }
        }
    }
    return Q / (2.0 * m);
}

// Function to calculation modularity gain
double modularityGain(const Graph& g, int u, int community_c, const std::map<int, int>& community, int m, const std::map<int, int>& community_totals) {
    double gain = 0.0;
    int sum_in = 0;
    int sum_tot = community_totals.at(community_c);

    for (int v = 0; v < g.num_nodes; ++v) {
        if (g.adjMatrix[u][v] > 0 && community.at(v) == community_c) {
            sum_in += g.adjMatrix[u][v];
        }
    }

    int k_i = g.degree(u);

    gain = (sum_in + k_i) / (2.0 * m) - std::pow(sum_tot + k_i, 2) / (2.0 * m * m)
           - sum_in / (2.0 * m) + std::pow(sum_tot, 2) / (2.0 * m * m) + std::pow(k_i, 2) / (2.0 * m * m);

    return gain;
}

// Function to detect communities using Louvain method
void louvainAlgorithm(Graph& g, std::map<int, int>& community) {
    int m = 0;
    for (int u = 0; u < g.num_nodes; ++u) {
        for (int v = 0; v < g.num_nodes; ++v) {
            m += g.adjMatrix[u][v];
        }
    }
    m /= 2;

    std::map<int, int> community_totals;
    for (int u = 0; u < g.num_nodes; ++u) {
        community_totals[community[u]] += g.degree(u);
    }

    bool improvement = true;
    while (improvement) {
        improvement = false;
        for (int u = 0; u < g.num_nodes; ++u) {
            int current_comm = community[u];
            double best_gain = 0.0;
            int best_community = current_comm;

            for (int v = 0; v < g.num_nodes; ++v) {
                if (g.adjMatrix[u][v] > 0 && community[u] != community[v]) {
                    int neighbor_comm = community[v];
                    double gain = modularityGain(g, u, neighbor_comm, community, m, community_totals);
                    if (gain > best_gain) {
                        best_gain = gain;
                        best_community = neighbor_comm;
                    }
                }
            }

            if (best_community != current_comm) {
                community[u] = best_community;
                improvement = true;
                community_totals[current_comm] -= g.degree(u);
                community_totals[best_community] += g.degree(u);
            }
        }
    }
}

int main() {
    // Initialize the graph
    Graph g(6);
    g.addEdge(0, 1, 1);  // ab
    g.addEdge(0, 2, 2);  // ac
    g.addEdge(1, 2, 3);  // bc
    g.addEdge(2, 4, 8);  // ce
    g.addEdge(3, 4, 1);  // de
    g.addEdge(4, 5, 3);  // ef
    g.addEdge(3, 5, 2);  // df

    // Initialize communities
    std::map<int, int> community;
    for (int i = 0; i < 6; ++i) {
        community[i] = i;  // Each node starts in its own community
    }

    // Apply Louvain algorithm
    louvainAlgorithm(g, community);

    // Output the results
    std::map<int, std::set<int>> communities;
    for (const auto& [node, comm] : community) {
        communities[comm].insert(node);
    }

    std::cout << "Communities:\n";
    for (const auto& [comm, nodes] : communities) {
        std::cout << "Community " << comm << ": ";
        for (int node : nodes) {
            std::cout << char('a' + node) << " ";
        }
        std::cout << "\n";
    }

    return 0;
}
