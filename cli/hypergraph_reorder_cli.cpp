#include "hypergraph_reorder/reorderer.hpp"
#include <iostream>
#include <string>
#include <iomanip>
#include <filesystem>

using namespace hypergraph_reorder;

const char* DEFAULT_CONFIG = "config/fast_km1_kKaHyPar.ini";

std::string find_config_file(const std::string& argv0) {
    namespace fs = std::filesystem;

    std::vector<std::string> search_paths = {
        DEFAULT_CONFIG,
        std::string("../") + DEFAULT_CONFIG,
        fs::path(argv0).parent_path() / DEFAULT_CONFIG,
        fs::path(argv0).parent_path() / ".." / DEFAULT_CONFIG,
    };

    for (const auto& path : search_paths) {
        if (fs::exists(path)) {
            return path;
        }
    }

    return "";
}

void print_usage(const char* prog_name) {
    std::cout << "Usage: " << prog_name << " <input_matrix> <k>\n\n";
    std::cout << "Arguments:\n";
    std::cout << "  <input_matrix>     Input matrix file (format auto-detected)\n";
    std::cout << "  <k>                Number of diagonal blocks/parts\n\n";
    std::cout << "Supported formats:\n";
    std::cout << "  .mtx    - Matrix Market format\n";
    std::cout << "  .graph  - METIS graph format\n";
    std::cout << "  other   - Binary CSR format (auto-detected)\n\n";
    std::cout << "Output: <input>_reordered.mtx\n\n";
    std::cout << "Options:\n";
    std::cout << "  -h, --help         Show this help message\n";
}

void print_statistics(const Statistics& stats) {
    std::cout << "\n========================================\n";
    std::cout << "Statistics\n";
    std::cout << "========================================\n";
    std::cout << "Matrix: " << stats.n_rows << " x " << stats.n_cols
              << ", " << stats.nnz << " nonzeros\n";
    std::cout << "Graph: " << stats.n_vertices << " vertices, "
              << stats.n_edges << " edges\n";
    std::cout << "Clique cover: " << stats.n_cliques << " cliques "
              << "(max: " << stats.max_clique_size
              << ", avg: " << std::fixed << std::setprecision(1) << stats.avg_clique_size << ")\n";
    std::cout << "Hypergraph: " << stats.n_hypernodes << " nodes, "
              << stats.n_hyperedges << " nets, " << stats.total_pins << " pins\n";
    std::cout << "Partition: " << stats.n_parts << " parts, separator size: "
              << stats.separator_size << " (" << std::fixed << std::setprecision(2)
              << (stats.separator_ratio * 100) << "%)\n";
    std::cout << "Ordering: " << stats.blocks_ordered << " blocks ordered ("
              << stats.blocks_failed << " failed)\n";
    std::cout << "\nTotal time: " << std::fixed << std::setprecision(2)
              << stats.time_total_ms << " ms\n";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    if (argc == 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")) {
        print_usage(argv[0]);
        return 0;
    }

    if (argc != 3) {
        std::cerr << "Error: Expected 2 arguments, got " << (argc - 1) << "\n\n";
        print_usage(argv[0]);
        return 1;
    }

    std::string input_path = argv[1];
    int n_parts = std::atoi(argv[2]);

    if (n_parts < 2) {
        std::cerr << "Error: Number of parts must be >= 2\n";
        return 1;
    }

    size_t dot_pos = input_path.find_last_of('.');
    std::string output_path;
    if (dot_pos != std::string::npos) {
        output_path = input_path.substr(0, dot_pos) + "_reordered.mtx";
    } else {
        output_path = input_path + "_reordered.mtx";
    }

    std::string config_path = find_config_file(argv[0]);
    if (config_path.empty()) {
        std::cerr << "Warning: Could not find " << DEFAULT_CONFIG << "\n";
        std::cerr << "Proceeding with KaHyPar defaults (may be slower)\n\n";
    }

    SymmetricDBReorderer::Options opts;
    opts.n_parts = n_parts;
    opts.imbalance = 0.03;
    opts.kahypar_config_path = config_path;
    opts.seed = -1;
    opts.ordering_method = OrderingMethod::AMD;
    opts.use_openmp = true;
    opts.num_threads = -1;
    opts.suppress_kahypar_output = false;

    try {
        SymmetricDBReorderer reorderer(opts);
        auto result = reorderer.reorder_from_file(input_path, MatrixFormat::AUTO);

        std::cout << "\nSaving reordered matrix to " << output_path << std::endl;
        write_matrix(output_path, result.reordered_matrix);
        std::cout << "Matrix saved successfully" << std::endl;

        print_statistics(result.stats);

        std::cout << "\n========================================\n";
        std::cout << "Reordering completed successfully!\n";
        std::cout << "========================================\n";

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
        return 1;
    }
}
