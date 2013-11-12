#ifndef PARSER_HPP
#define PARSER_HPP

/**
 * @author pkambadu
 */

#include <string>
#include <cstdlib>
#include <ext/hash_map>

#include "utilities.hpp"

static const char* description = "Stochastic Gradient Descent";
typedef __gnu_cxx::hash_map<const char*, const char*> string_string_map_t;
typedef __gnu_cxx::hash_map<const char*, int> string_int_map_t;

static inline void print_help (const string_string_map_t& help_map,
                               const char* option=NULL,
                               const char* val=NULL) {
  if (NULL != option) printf ("%s: Invalid arguments \"%s\" is \"%s\"\n", 
                                                description, option, val);
  else printf ("%s\n", description);
  string_string_map_t::const_iterator iter = help_map.begin();

  while (iter != help_map.end()) {
    printf ("%s: %s\n", (*iter).first, (*iter).second);
    ++iter;
  }

  exit (3);
}

static inline void parse_parameters (int argc, char** argv) {
  /** Set up the command line arguments */
  string_string_map_t help_map;
  help_map["help"] = "produce help messages";
  help_map["max-epochs"] = "maximum number of epochs (default:10)";
  help_map["epoch-size"] = "number of iterations before sync (default:10)";
  help_map["M"] = "number of rows/cols in A, cols in X, and rows in Y";
  help_map["num-threads"] = "number of threads to start (default:2)";
  help_map["debug"] = "print helpful messages out (default:0)";
  help_map["seed"] = "seed for the random number generator (default:0)";
  help_map["A-file-path"] = "path to the input A (matrix market format)";
  help_map["Y-file-path"] = "path to the input Y (text format)";

  string_int_map_t int_options_map;
  int_options_map["max-epochs"]        = MAX_EPOCHS_INDEX;
  int_options_map["epoch-size"]        = EPOCH_SIZE_INDEX;;
  int_options_map["M"]                 = M_INDEX;
  int_options_map["debug"]             = DEBUG_INDEX;
  int_options_map["num-threads"]       = NUM_THREADS_INDEX;
  int_options_map["seed"]              = RAND_SEED_INDEX;

  string_int_map_t dbl_options_map;

  string_int_map_t chr_options_map;
  chr_options_map["A-file-path"]   = A_FILE_PATH_INDEX;
  chr_options_map["Y-file-path"]   = Y_FILE_PATH_INDEX;

  /* default initialize parameters */
  int_params[MAX_EPOCHS_INDEX]        = 10;
  int_params[EPOCH_SIZE_INDEX]        = 10;
  int_params[M_INDEX]                 = -1;
  int_params[DEBUG_INDEX]             = 0;
  int_params[NUM_THREADS_INDEX]       = 2;
  int_params[RAND_SEED_INDEX]         = 0;

  chr_params[A_FILE_PATH_INDEX] = "";
  chr_params[Y_FILE_PATH_INDEX] = "";

  /* parse the command line */
  if (!(argc&0x1)) print_help(help_map);

  for (int i=1; i<argc; i+=2) {
    char* option = argv[i];
    char* value = argv[i+1];

    if (0==strcmp(option,"help")) print_help (help_map);
    else if (int_options_map.end()!=int_options_map.find(option)) {
      int_params[int_options_map[option]] = atoi(value);
    } else if (dbl_options_map.end()!=dbl_options_map.find(option)) {
      dbl_params[dbl_options_map[option]] = atof(value);
    } else if (chr_options_map.end()!=chr_options_map.find(option)) {
      chr_params[chr_options_map[option]] = value;
    } else {
      print_help (help_map,option);
    }
  }

  /* Make sure that all the parameters are given and are correct */
  if (0==strcmp("",chr_params[A_FILE_PATH_INDEX])) {
    printf ("A matrix path has not been given\n");
    print_help (help_map, "A-file-path");
  } else if (0==strcmp("",chr_params[Y_FILE_PATH_INDEX])) {
    printf ("Y matrix path has not been given\n");
    print_help (help_map, "Y-file-path");
  } else if (0>int_params[M_INDEX]) {
    printf ("M is invalid\n");
    print_help (help_map, "M");
  } 
}

#endif // PARSER_HPP
