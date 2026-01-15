#ifndef HYPERGRAPH_REORDER_IO_HPP
#define HYPERGRAPH_REORDER_IO_HPP

#include "sparse_matrix.hpp"
#include <string>

namespace hypergraph_reorder {

enum class MatrixFormat {
    MTX,
    METIS,
    CSR_BINARY,
    AUTO
};

MatrixFormat detect_format(const std::string& filepath);

CSRMatrix read_matrix(const std::string& filepath,
                      MatrixFormat format = MatrixFormat::AUTO);

void write_matrix(const std::string& filepath,
                  const CSRMatrix& matrix,
                  MatrixFormat format = MatrixFormat::MTX);

namespace mtx {
    CSRMatrix read(const std::string& filepath);
    void write(const std::string& filepath, const CSRMatrix& matrix);
}

namespace metis {
    CSRMatrix read(const std::string& filepath);
    void write(const std::string& filepath, const CSRMatrix& matrix);
}

namespace csr_binary {
    CSRMatrix read(const std::string& filepath);
    void write(const std::string& filepath, const CSRMatrix& matrix);
}

} // namespace hypergraph_reorder

#endif // HYPERGRAPH_REORDER_IO_HPP
