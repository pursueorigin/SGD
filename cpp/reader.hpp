#ifndef READER_HPP
#define READER_HPP

/**
 * @author pkambadu
 */

#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <iterator>
#include <fstream>
#include <vector>

template <typename OutputIterator>
void reader (const char* file_path,
             OutputIterator result) {
  typedef typename std::iterator_traits<OutputIterator>::value_type
                            integral_value_type;
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer_type;
  boost::char_separator<char> separators(" ");

  std::ifstream in(file_path);
  if (!in.is_open()) { printf ("Could not open %s\n", file_path); return; }

  std::string current_line;
  while (getline(in,current_line)) {
    tokenizer_type tokens(current_line, separators);
    for (tokenizer_type::iterator token_iter = tokens.begin();
         token_iter != tokens.end(); 
         ++token_iter) {
      *result++=boost::lexical_cast<integral_value_type>(token_iter->c_str());
    }
  }
}

#endif // READ_SNP_MATRIX_HPP
