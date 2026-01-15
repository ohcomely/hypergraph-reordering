// partitioner.cpp - KaHyPar integration implementation
#include "hypergraph_reorder/partitioner.hpp"
#include <libkahypar.h>
#include <iostream>
#include <unordered_set>
#include <algorithm>

namespace hypergraph_reorder {

HypergraphPartitioner::HypergraphPartitioner(const Options& opts)
    : opts_(opts), kahypar_context_(nullptr) {
    init_kahypar_context();
}

HypergraphPartitioner::~HypergraphPartitioner() {
    cleanup_kahypar_context();
}

void HypergraphPartitioner::init_kahypar_context() {
    kahypar_context_ = kahypar_context_new();

    if (!opts_.kahypar_config_path.empty()) {
        kahypar_configure_context_from_file(
            static_cast<kahypar_context_t*>(kahypar_context_),
            opts_.kahypar_config_path.c_str()
        );
    }

    if (opts_.seed >= 0) {
        kahypar_set_seed(
            static_cast<kahypar_context_t*>(kahypar_context_),
            opts_.seed
        );
    }

    if (opts_.suppress_output) {
        kahypar_supress_output(
            static_cast<kahypar_context_t*>(kahypar_context_),
            true
        );
    }
}

void HypergraphPartitioner::cleanup_kahypar_context() {
    if (kahypar_context_) {
        kahypar_context_free(static_cast<kahypar_context_t*>(kahypar_context_));
        kahypar_context_ = nullptr;
    }
}

HypergraphPartition HypergraphPartitioner::partition(const Hypergraph& hg) {
    if (hg.n_nets() == 0) {
        std::cout << "Warning: Empty hypergraph, creating trivial partition" << std::endl;

        HypergraphPartition result;
        result.n_parts = opts_.n_parts;
        result.node_partition.resize(hg.n_nodes());
        for (index_t i = 0; i < hg.n_nodes(); ++i) {
            result.node_partition[i] = i % opts_.n_parts;
        }
        result.objective = 0;
        result.part_sizes.resize(opts_.n_parts, 0);
        for (auto p : result.node_partition) {
            result.part_sizes[p]++;
        }
        return result;
    }

    std::cout << "Partitioning CNH: " << hg.n_nodes() << " nodes, "
              << hg.n_nets() << " nets into " << opts_.n_parts << " parts" << std::endl;

    // Convert edges to unsigned int (KaHyPar expects uint, but we use int64_t internally)
    std::vector<unsigned int> edges_uint(hg.edges().size());
    for (size_t i = 0; i < hg.edges().size(); ++i) {
        edges_uint[i] = static_cast<unsigned int>(hg.edges()[i]);
    }

    // Convert edge_indices to unsigned int
    std::vector<size_t> edge_indices_copy = hg.edge_indices();

    // Create KaHyPar hypergraph
    kahypar_hypergraph_t* kahypar_hg = kahypar_create_hypergraph(
        opts_.n_parts,
        static_cast<unsigned int>(hg.n_nodes()),
        static_cast<unsigned int>(hg.n_nets()),
        edge_indices_copy.data(),
        edges_uint.data(),
        hg.edge_weights().data(),
        hg.node_weights().data()
    );

    if (!kahypar_hg) {
        throw HypergraphReorderError("Failed to create KaHyPar hypergraph");
    }

    // Partition
    std::vector<kahypar_partition_id_t> partition(hg.n_nodes());
    kahypar_hyperedge_weight_t objective = 0;

    kahypar_partition_hypergraph(
        kahypar_hg,
        opts_.n_parts,
        opts_.imbalance,
        &objective,
        static_cast<kahypar_context_t*>(kahypar_context_),
        partition.data()
    );

    // Convert to result
    HypergraphPartition result;
    result.n_parts = opts_.n_parts;
    result.node_partition.resize(hg.n_nodes());
    result.part_sizes.resize(opts_.n_parts, 0);

    for (index_t i = 0; i < hg.n_nodes(); ++i) {
        result.node_partition[i] = partition[i];
        result.part_sizes[partition[i]]++;
    }

    result.objective = objective;

    std::cout << "Partition objective (cutsize): " << objective << std::endl;
    std::cout << "Part sizes: ";
    for (auto sz : result.part_sizes) {
        std::cout << sz << " ";
    }
    std::cout << std::endl;

    // Cleanup
    kahypar_hypergraph_free(kahypar_hg);

    return result;
}

VertexPartition HypergraphPartitioner::create_vertex_partition(
    const HypergraphPartition& cnh_partition,
    const CliqueCover& cover,
    index_t n_vertices) {

    std::cout << "Creating vertex separator from CNH partition..." << std::endl;

    VertexPartition result;
    result.n_parts = cnh_partition.n_parts;
    result.parts.resize(result.n_parts);

    // For each vertex, determine if it belongs to a single part or separator
    for (index_t v = 0; v < n_vertices; ++v) {
        auto cliques = cover.get_vertex_cliques(v);

        if (cliques.empty()) {
            // Isolated vertex -> separator
            result.separator.push_back(v);
            continue;
        }

        // Get unique partitions this vertex belongs to
        std::unordered_set<index_t> vertex_parts;
        for (auto clique_id : cliques) {
            vertex_parts.insert(cnh_partition.node_partition[clique_id]);
        }

        if (vertex_parts.size() == 1) {
            // Vertex belongs to single partition
            index_t part_id = *vertex_parts.begin();
            result.parts[part_id].push_back(v);
        } else {
            // Vertex spans multiple partitions -> separator
            result.separator.push_back(v);
        }
    }

    // Sort for consistency
    for (auto& part : result.parts) {
        std::sort(part.begin(), part.end());
    }
    std::sort(result.separator.begin(), result.separator.end());

    std::cout << "Vertex partition completed:" << std::endl;
    std::cout << "  Part sizes: ";
    for (const auto& part : result.parts) {
        std::cout << part.size() << " ";
    }
    std::cout << std::endl;
    std::cout << "  Separator size: " << result.separator.size()
              << " (" << (100.0 * result.separator.size() / n_vertices) << "%)" << std::endl;

    return result;
}

} // namespace hypergraph_reorder
