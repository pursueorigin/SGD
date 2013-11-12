#ifndef MATRIX_MARKET_READER_HPP
#define MATRIX_MARKET_READER_HPP

#include <utility>
#include "matrix_market/mmio.h"

struct read_exception: public std::exception {
  const char* reason;
  read_exception (const char* reason) : reason(reason) {}
  virtual const char* what() const throw() { return reason; }
};


/**
 * Read in from Market Matrix format. The value_type of the
 * RandomAccessContainer is a ForwardContainer. The implicit requirements are
 * also that the values that the ForwardContainer are of type std::pair<int,
 * double> >.
 *
 * \param[in] input_filename File name that acts as our input.
 * \param[in] has_to_be_symmetric Name says it all.
 * \param[inout] adjacency_list A container to hold the adjacencies.
 */
template <typename RandomAccessContainer>
static std::pair<int,int> 
read_matrix_market (const char* input_filename,
                    const bool& has_to_be_symmetric,
                    RandomAccessContainer& adjacency_list) {

  typedef typename RandomAccessContainer::value_type ForwardContainer;
  typedef typename ForwardContainer::value_type EdgeType;
  typedef typename EdgeType::first_type VertexIndexType;
  typedef typename EdgeType::second_type EdgeWeightType;

  int M;
  int N;
  int nnz;
  FILE* input_fp; 
  MM_typecode matcode;
 
  /* open the file */
  if (NULL==(input_fp=fopen (input_filename, "r"))) {
    fclose (input_fp);
    throw read_exception ("Could not open input file for reading");
  }
 
  /* read the banner in */
  if (0 != mm_read_banner (input_fp, &matcode)) {
    fclose (input_fp);
    throw read_exception("Could not process Matrix Market banner.\n");
  }
  
  /*  We only support sparse, real/integer-valued edge weights */
  if (!mm_is_matrix (matcode) || 
      mm_is_dense (matcode) ||
      mm_is_complex (matcode) ||
      mm_is_integer (matcode) ||
      mm_is_pattern (matcode) ||
      mm_is_hermitian (matcode)) {
    fclose (input_fp);
    throw read_exception("Reader does not support the MATRIX_MARKET file");
  }

  /* find out size of sparse matrix */
  if (0 != mm_read_mtx_crd_size(input_fp, &M, &N, &nnz)) {
    fclose (input_fp);
    throw read_exception("Could not determine size of the matrix");
  }
 
  /* make sure that M and N are the same */
  if (M!=N) {
    fclose (input_fp);
    throw read_exception("Input is not an adjacency matrix, M!=N");
  }
 
  /* 
   * Figure out the type of matrix we are dealing with here. If the matrix 
   * is skew-symmetric or symmetric, we need to store the forward and backward
   * edges separately, thereby increasing the number of non-zeros. The only 
   * distinction is that in the case of symmetric matrices, we may have 
   * non-zero diagonal entries and for these, we don't have to store forward
   * and backward directions; only one copy of the edge will suffice.
   */
  bool symmetric = (mm_is_symmetric(matcode))? true:
                     ((mm_is_skew(matcode))? true: false);

  if (has_to_be_symmetric && (!symmetric)) {
    fclose (input_fp);
    throw read_exception("Input is not a symmetric matrix");
  }
 
  /* Create space for the adjacency list */
  adjacency_list.resize (N);

  /* A variable that holds the actual number of non-zeros. */
  int actual_nnz = 0;
 
  /* 
   * Read in the matrix --- adjust for 1-based indexing of this format. 
   */
  for (int i=0; i<nnz; ++i) {
    int source;
    int target;
    double weight;
 
    /* scan the line to get source, target, weight */
    if (3 != fscanf (input_fp, "%d %d %lg\n", &source, &target, &weight)) {
      throw read_exception("Malformed edge");
    }
 
    /* adjust to make 0 indexed */
    --source; --target;

 
    /* Add the forward edge to the adjacency_list */
    adjacency_list [source].push_back (EdgeType(target, weight));
    ++actual_nnz;
 
    /* 
     * Add the reverse edge if needed. The 'source!=target' takes care
     * of symmetric matrices' special case. 
     */
    if (symmetric && (source!=target)) {
      adjacency_list [target].push_back (EdgeType(source, weight));
      ++actual_nnz;
    }
  }

  return std::pair<int,int> (N, actual_nnz);
}

#endif // MATRIX_MARKET_READER_HPP
