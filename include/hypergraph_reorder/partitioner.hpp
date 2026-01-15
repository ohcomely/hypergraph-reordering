// partitioner.hpp - KaHyPar integration for hypergraph partitioning
#ifndef HYPERGRAPH_REORDER_PARTITIONER_HPP
#define HYPERGRAPH_REORDER_PARTITIONER_HPP

#include "types.hpp"
#include "hypergraph.hpp"
#include "clique_cover.hpp"
#include <string>
#include <vector>

namespace hypergraph_reorder {

// Vertex partition (result of converting CNH partition to vertex separator)
struct VertexPartition {
    index_t n_parts;
    std::vector<std::vector<index_t>> parts;  // parts[i] = vertices in part i
    std::vector<index_t> separator;            // Separator vertices

    index_t separator_size() const { return separator.size(); }
    double separator_ratio(index_t n_total) const {
        return static_cast<double>(separator.size()) / n_total;
    }
};

// Hypergraph partitioner
class HypergraphPartitioner {
public:
    struct Options {
        index_t n_parts;
        double imbalance;
        std::string kahypar_config_path;
        int seed;
        bool suppress_output;

        Options()
            : n_parts(4),
              imbalance(0.03),
              kahypar_config_path(""),
              seed(-1),
              suppress_output(false) {}
    };

    explicit HypergraphPartitioner(const Options& opts = Options());
    ~HypergraphPartitioner();

    // Partition hypergraph using KaHyPar
    HypergraphPartition partition(const Hypergraph& hg);

    // Convert CNH partition to vertex separator partition
    VertexPartition create_vertex_partition(
        const HypergraphPartition& cnh_partition,
        const CliqueCover& cover,
        index_t n_vertices);

private:
    Options opts_;
    void* kahypar_context_;  // Opaque KaHyPar context

    void init_kahypar_context();
    void cleanup_kahypar_context();
};

} // namespace hypergraph_reorder

#endif // HYPERGRAPH_REORDER_PARTITIONER_HPP
