#ifndef WEIGHTED_CSR_GRAPH_HPP
#define WEIGHTED_CSR_GRAPH_HPP

#include <cassert>
#include <vector>
#include <algorithm>

template <typename RandomAccessContainer>
struct weighted_csr_graph {

  typedef typename RandomAccessContainer::value_type ForwardContainer;
  typedef typename ForwardContainer::const_iterator ForwardContainerIterator;
  typedef typename ForwardContainer::value_type EdgeType;
  typedef typename EdgeType::first_type VertexIndexType;
  typedef typename EdgeType::second_type EdgeWeightType;
  typedef typename std::pair<VertexIndexType, VertexIndexType> VertexPairType;

  private:
  int N;
  int nnz;
  bool sorted;
  std::vector<VertexIndexType> row_indices;
  std::vector<VertexIndexType> columns;
  std::vector<EdgeWeightType> weights;
  std::vector<int> integer_weights;

  public:
  /**
   * Construct a CSR format graph from an adjacency_list representation.
   * \param[in] adjacency_list An adjacency list representation.
   * \param[in] N The number of vertices.
   * \param[in] nnz The number of non-zero entries. 
   *
   * NOTE: We can deduce N and nnz from adjacency_list, but its a lot easier
   * to just get the whole thing passed as parameters. Avoids going through 
   * the adjacency_list twice.
   */
  weighted_csr_graph (const RandomAccessContainer& adjacency_list,
                      const int& N,
                      const int& nnz):
                      N (adjacency_list.size()), nnz(nnz), sorted (false) {
    row_indices.resize (N+1);
    columns.resize (nnz);
    weights.resize (nnz);

    /* Start copying over from the first vertex onwards */
    int offset = 0;
    for (int vertex_index=0; vertex_index<N; ++vertex_index) {
      row_indices[vertex_index] = offset;
      ForwardContainerIterator edge_iter=adjacency_list[vertex_index].begin();
      ForwardContainerIterator edge_end = adjacency_list[vertex_index].end();
      while (edge_iter!=edge_end) {
        columns[offset] = (*edge_iter).first;
        weights[offset] = (*edge_iter).second;
        ++offset;
        ++edge_iter;
      }
    }

    /* check that we have included every edge */
    assert (offset==nnz);

    /* set the offset of the sentinel */
    row_indices[N] = offset;
  }

  /**
   * Copy constructor
   * @param[in] original_graph The graph from which to copy.
   * @param[in] sample_counts An integer array containing information on edges.
   * @param[in] probabilities Normalized probabilities of an edge for weighting
   * @param[in] tau Used to reweight the edges.
   */
  weighted_csr_graph (const weighted_csr_graph& original_graph,
                      const std::vector<bool>& samples) :
           N (original_graph.get_num_vertices()), nnz(0), sorted(false) {
    /* Figure out the number of non zeros */
    for (int i=0; i<samples.size(); ++i) if (samples[i]) ++nnz;
  
    /* resize the containers */
    row_indices.resize (N+1);
    columns.resize (nnz);
    weights.resize (nnz);

    /* Copy over those edges that have been sampled and adjust weights */
    int new_target_index = 0;
    for (int source=0; source<N; ++source) {
      row_indices [source] = new_target_index;
      for (int target_index=original_graph.begin (source);
           target_index<original_graph.end (source);
           ++target_index) {
        if (samples [target_index]) {
          const int target = original_graph.get_target (target_index);
          columns[new_target_index] = target;
          const EdgeWeightType weight=original_graph.get_weight (target_index);
          weights[new_target_index] = weight;
          ++new_target_index;
        }
      }
    }

    assert (nnz == new_target_index);
    row_indices [N] = nnz;
  }

  /**
   * Copy constructor
   * @param[in] original_graph The graph from which to copy.
   * @param[in] sample_counts An integer array containing information on edges.
   * @param[in] probabilities Normalized probabilities of an edge for weighting
   * @param[in] tau Used to reweight the edges.
   */
  weighted_csr_graph (const weighted_csr_graph& original_graph,
                      const std::vector<int>& sample_counts,
                      const std::vector<double>& probabilities,
                      const double& tau) :
           N (original_graph.get_num_vertices()), nnz(0), sorted(false) {
    /* Figure out the number of non zeros */
    for (int i=0; i<sample_counts.size(); ++i) if (sample_counts[i]) ++nnz;
  
    /* resize the containers */
    row_indices.resize (N+1);
    columns.resize (nnz);
    weights.resize (nnz);

    /* Copy over those edges that have been sampled and adjust weights */
    int new_target_index = 0;
    for (int source=0; source<N; ++source) {
      row_indices [source] = new_target_index;
      for (int target_index=original_graph.begin (source);
           target_index<original_graph.end (source);
           ++target_index) {
        if (sample_counts [target_index]) {
          const int target = original_graph.get_target (target_index);
          columns[new_target_index] = target;
          const EdgeWeightType old_weight = 
                                  original_graph.get_weight (target_index);
          const double edge_probability = tau*probabilities [target_index];
          const int num_selections = sample_counts [target_index];
          weights[new_target_index] = 
                    num_selections * (old_weight/edge_probability);
          ++new_target_index;
        }
      }
    }

    assert (nnz == new_target_index);
    row_indices [N] = nnz;
  }

  /**
   * \return The number of vertices in this graph
   */
  int get_num_vertices () const { return N; }

  /**
   * \return The number of edges in this graph
   */
  int get_num_edges () const { return nnz; }

  /**
   * Part of the iterator interface. 
   * \param[in] vertex The vertex whose beginning position we seek.
   * \return The index that can be used to access columns and weights.
   */
  int begin (const VertexIndexType& vertex) const { 
    assert (N>vertex);
    return row_indices[vertex];
  }

  /**
   * Part of the iterator interface. 
   * \param[in] vertex The vertex whose ending position we seek.
   * \return The index that can be used to access columns and values.
   */
  int end (const VertexIndexType& vertex) const { 
    assert (N>vertex);
    return row_indices[vertex+1];
  }

  /**
   * Get the target of the edge stored in the index.
   * \param[in] The index in the CSR storage.
   * \return The destination of the edge.
   */
  VertexIndexType get_target (const int& index) const {
    assert (nnz>index);
    return columns[index];
  }

  /**
   * Get the weight of the edge stored in the index.
   * \param[in] The index in the CSR storage.
   * \return The destination of the edge.
   */
  EdgeWeightType get_weight (const int& index) const {
    assert (nnz>index);
    return weights[index];
  }

  /**
   * \param[in] vertex The vertex for which we want the out degree.
   * \return The outdegree of the vertex.
   */
  int get_out_degree (const VertexIndexType& vertex) const {
    assert (N>vertex);
    return (row_indices[vertex+1]-row_indices[vertex]); 
  }

  /**
   * @return the row_indices as an array.
   */
  VertexIndexType* get_row_indices () { return &(row_indices[0]); }

  /**
   * @return the row_indices as an array.
   */
  VertexIndexType* get_columns () { return &(columns[0]); }

  /**
   * @return the weights as an array.
   */
  EdgeWeightType* get_weights () { return &(weights[0]); }

  /**
   * @param[in] scaling_factor To scale the weights.
   * @return integral the weights as an array. Multiply by scaling factor.
   */
  int* get_integer_weights (const EdgeWeightType scaling_factor) { 
    integer_weights.resize (nnz);
    for (int i=0; i<nnz; ++i) 
      integer_weights[i] = static_cast<int> (weights[i]*scaling_factor);
    return &(integer_weights[0]); 
  }

  /**
   * Compute the aggregate weight of the edges cut based on the partitioning
   * specified. We don't really care about the number of partitions here.
   * @param[in] partition_map Specifies the partition for each vertex.
   */
  EdgeWeightType compute_edges_cut(const std::vector<int>& partition_map)const{
    assert (partition_map.size () == N);

    /* 
     * Go through the adjacency map and add up edges which are not in the 
     * same partition.
     */
    EdgeWeightType edges_cut = 0.0;
    for (int source=0; source<N; ++source) {
      for (int target_index=begin(source); 
           target_index<end(source);
           ++target_index) {
        const int target = get_target (target_index);
        if (partition_map[target] != partition_map[source]) {
          edges_cut += get_weight (target_index);
        }
      }
    }

    return edges_cut;
  }

  /**
   * Compute the maximum imbalance for the partitioning specified. The maximum
   * imbalance is the ratio of the heaviest partition to that of the average 
   * partition weight. For each vertex, the partition is a number between 0 and
   * (num_partitions-1).
   * @param[in] partition_map Specifies the partition for each vertex.
   * @param[in] num_partitions Each vertex belongs to [0,num_partitions).
   */
  double compute_max_imbalance (const std::vector<int>& partition_map,
                                const int& num_partitions) const {
    assert (partition_map.size () == N);
    std::vector<EdgeWeightType> partition_weights (num_partitions, 0.0);

    /* Go through the adjacency map and sum up the partition weights */
    for (int source=0; source<N; ++source) {
      /* If and when we have vertex weights, we can change this line */
      partition_weights [partition_map [source]] += 1.0;
    }

    /* Figure out the heaviest partition and the average weight */
    EdgeWeightType average_weight = 0.0;
    EdgeWeightType maximum_weight = 0.0;
    for (int partition=0; partition<num_partitions; ++partition) {
      const EdgeWeightType current_weight = partition_weights [partition];
      maximum_weight = (current_weight>maximum_weight) ?
                        current_weight:maximum_weight;
      average_weight += current_weight;
    }
    average_weight /= num_partitions;

    const double max_imbalance = 
      ((maximum_weight-average_weight)/average_weight) * 100.0;

    return static_cast<double>((maximum_weight/average_weight));
  }

  /**
   * Pretty print the graph
   */
  void pretty_print () const {
    for (int source = 0; source < N; ++source) {
      for (int target_index = begin (source);
           target_index < end (source);
           ++target_index) {
        std::cout << source << " ---> " << get_target (target_index) 
                  << " (weight=" << get_weight (target_index) << ")" 
                  << std::endl;
      }
    }
  }

  /**
   * Sort the edges in ascending order so that we can search for things 
   * far more easily later. This procedure repeatedly calls binary search.
   */
  void sort () {
    /*
     * For each vertex's adjacency, we want to sort the targets by their 
     * vertex index. For this, we can directly use std::sort(). The only 
     * catch is that we need to reorganize the weights based on the sorted
     * values. So, what we do instead is std::sort() on a copy of the 
     * adjacencies, and once then reorganize the weights to match the new 
     * values.
     */
    std::vector<VertexIndexType> columns_copy(columns.begin(), columns.end());
    std::vector<EdgeWeightType> weights_copy (nnz, 0.0);
    for (int source = 0; source < N; ++source) {
      const int target_begin = begin (source);
      const int target_end = end (source);
      std::sort ((columns_copy.begin()+target_begin), 
                 (columns_copy.begin()+target_end));

      for (int target_index = begin (source);
           target_index < end (source);
           ++target_index) {
        const VertexIndexType target = get_target (target_index);
        const EdgeWeightType weight = get_weight (target_index);

        const int target_new_position = std::lower_bound 
                                    ((columns_copy.begin()+target_begin),
                                     (columns_copy.begin()+target_end),
                                     target) - columns_copy.begin();

        if (target_new_position == 
            ((columns_copy.begin()+target_end)-columns_copy.begin())) {
          fprintf (stderr, "Could not locate (%d,%d) in graph\n", 
                                              source, target);
          exit (1);
        }

        weights_copy [target_new_position] = weight;
      }
    }
    columns.swap (columns_copy);
    weights.swap (weights_copy);
    sorted = true;
  }

  /**
   * Get the index at which this current edge is found. Only works on 
   * sorted graphs!
   * @param[in] source The source vertex.
   * @param[in] target The target vertex.
   * @return The index of the edge in columns array such that 
   *         columns[return_value]=(source,target).
   */
  int get_target_index (const VertexIndexType& source, 
                        const VertexIndexType& target) const {
    if (source>=N || target>=N) {
      fprintf (stderr, "(%d,%d) is an illegal edge\n", source, target);
      exit (1);
    } else if (sorted) {
      const int required_target_index = std::lower_bound 
                                ((columns.begin()+begin(source)),
                                 (columns.begin()+end(source)),
                                 target) - columns.begin();

      if (required_target_index == 
          ((columns.begin()+end(source))-columns.begin())) {
        fprintf (stderr, "Could not locate (%d,%d) in graph\n", 
                                            source, target);
        exit (1);
      }
      return required_target_index;
    } else {
      fprintf (stderr, "Graph is not sorted!\n");
      exit (1);
    }
  }

  /**
   * Get the source and the target given an index in the columns array.
   * @param[in] target_index Index of the edge we need to look up.
   * @return A pair which contains the source and the target.
   */
  VertexPairType get_edge (const int& target_index) const {
    if (0>target_index || target_index>=nnz) {
      fprintf (stderr, "Invalid target index: %d\n", target_index); 
      exit (1);
    } else if (!sorted) {
      fprintf (stderr, "This procedure only works on sorted graphs\n");
      exit (1);
    } else {
      const VertexIndexType source = std::upper_bound
                                (row_indices.begin(),
                                 row_indices.end(),
                                 target_index) - row_indices.begin() - 1; 
      const VertexIndexType target = columns [target_index];
      return std::pair<VertexIndexType,VertexIndexType> (source, target);
    }
  }

  /**
   * Given an index in the columns array, figure out the edge that is 
   * the exact reverse of this edge --- ONLY WORKS FOR SYMMETRIC GRAPHS!
   * @param[in] target_index Index of the target that we want to reverse.
   * @return The reverse edge index if it exits.
   */
  int get_reverse_target_index (const int& target_index) const {
    VertexPairType edge_pair = get_edge (target_index);
    return get_target_index (edge_pair.second, edge_pair.first);
  }

  void diag(std::vector<EdgeWeightType> &d) const {
    for (int i = 0; N > i; ++i) {
      int t;
      for (t = row_indices[i]; row_indices[i + 1] > t && columns[t] != i; ++t)
        ;
      d[i] = row_indices[i + 1] > t ? weights[t] : 0;
    } // for each row
  }

  /**
   * Computes y <- y + Ax
   */
  void matvec(std::vector<EdgeWeightType> &y,
              const std::vector<EdgeWeightType> &x) const {
    for (int i = 0; N > i; ++i) {
      for (int t = row_indices[i]; row_indices[i + 1] > t; ++t) {
        int j = columns[t];
        EdgeWeightType Aij = weights[t];
        y[i] = y[i] + Aij * x[j];
      } // for each nnz in row i
    } // for each row
  }

};

#endif /* WEIGHTED_CSR_GRAPH_HPP */
