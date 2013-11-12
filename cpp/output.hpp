#ifndef OUTPUT_HPP
#define OUTPUT_HPP

#include <cstdio>

/**
 * Shows the progress of the algorithm by printing messages to stdout.
 */
struct output_formatter_t {
public:
  void print_hdr() {
    std::printf("%-15s" "%5s" "%-15s\n", "epochs", "", "||r||_2 / ||b||_2");
  }

  /**
   * Prints the residual norm.
   *
   * @param[in] epochs num of completed epochs
   * @param[in] residual 2-norm of the residual
   */
  void print_progress(int epochs, double residual) {
    std::printf("%15d" "%5s" "%15.4e\n", epochs, "", residual);
  }
};

#endif // OUTPUT_HPP
