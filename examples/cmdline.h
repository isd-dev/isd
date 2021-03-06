/** @file cmdline.h
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.22.2
 *  http://www.gnu.org/software/gengetopt.
 *  DO NOT modify this file, since it can be overwritten
 *  @author GNU Gengetopt by Lorenzo Bettini */

#ifndef CMDLINE_H
#define CMDLINE_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef CMDLINE_PARSER_PACKAGE
/** @brief the program name (used for printing errors) */
#define CMDLINE_PARSER_PACKAGE "findlowweight"
#endif

#ifndef CMDLINE_PARSER_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#define CMDLINE_PARSER_PACKAGE_NAME "findlowweight"
#endif

#ifndef CMDLINE_PARSER_VERSION
/** @brief the program version */
#define CMDLINE_PARSER_VERSION "1.0"
#endif

/** @brief Where the command line options are stored */
struct gengetopt_args_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  int length_arg;	/**< @brief length of the code.  */
  char * length_orig;	/**< @brief length of the code original value given at command line.  */
  const char *length_help; /**< @brief length of the code help description.  */
  int dimension_arg;	/**< @brief dimension of the code.  */
  char * dimension_orig;	/**< @brief dimension of the code original value given at command line.  */
  const char *dimension_help; /**< @brief dimension of the code help description.  */
  int weight_arg;	/**< @brief weight threshold.  */
  char * weight_orig;	/**< @brief weight threshold original value given at command line.  */
  const char *weight_help; /**< @brief weight threshold help description.  */
  int seed_arg;	/**< @brief seed of the random number generator (default='24301').  */
  char * seed_orig;	/**< @brief seed of the random number generator original value given at command line.  */
  const char *seed_help; /**< @brief seed of the random number generator help description.  */
  int max_iter_arg;	/**< @brief maximum number of iteration (default='0').  */
  char * max_iter_orig;	/**< @brief maximum number of iteration original value given at command line.  */
  const char *max_iter_help; /**< @brief maximum number of iteration help description.  */
  int max_found_arg;	/**< @brief maximum number of different solutions (default='1').  */
  char * max_found_orig;	/**< @brief maximum number of different solutions original value given at command line.  */
  const char *max_found_help; /**< @brief maximum number of different solutions help description.  */
  int px_arg;	/**< @brief size of the tuples for X (default='1').  */
  char * px_orig;	/**< @brief size of the tuples for X original value given at command line.  */
  const char *px_help; /**< @brief size of the tuples for X help description.  */
  int py_arg;	/**< @brief size of the tuples for Y (default='1').  */
  char * py_orig;	/**< @brief size of the tuples for Y original value given at command line.  */
  const char *py_help; /**< @brief size of the tuples for Y help description.  */
  int kx_arg;	/**< @brief number of rows in X (default='0').  */
  char * kx_orig;	/**< @brief number of rows in X original value given at command line.  */
  const char *kx_help; /**< @brief number of rows in X help description.  */
  int ky_arg;	/**< @brief number of rows in Y (default='0').  */
  char * ky_orig;	/**< @brief number of rows in Y original value given at command line.  */
  const char *ky_help; /**< @brief number of rows in Y help description.  */
  int m_arg;	/**< @brief multiple widows before pivoting (default='1').  */
  char * m_orig;	/**< @brief multiple widows before pivoting original value given at command line.  */
  const char *m_help; /**< @brief multiple widows before pivoting help description.  */
  int l_arg;	/**< @brief length of the window (default='8').  */
  char * l_orig;	/**< @brief length of the window original value given at command line.  */
  const char *l_help; /**< @brief length of the window help description.  */
  int c_arg;	/**< @brief columns swapped for each pivot (default='1').  */
  char * c_orig;	/**< @brief columns swapped for each pivot original value given at command line.  */
  const char *c_help; /**< @brief columns swapped for each pivot help description.  */
  int r_arg;	/**< @brief size of the groups for pivoting (default='1').  */
  char * r_orig;	/**< @brief size of the groups for pivoting original value given at command line.  */
  const char *r_help; /**< @brief size of the groups for pivoting help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int length_given ;	/**< @brief Whether length was given.  */
  unsigned int dimension_given ;	/**< @brief Whether dimension was given.  */
  unsigned int weight_given ;	/**< @brief Whether weight was given.  */
  unsigned int seed_given ;	/**< @brief Whether seed was given.  */
  unsigned int max_iter_given ;	/**< @brief Whether max_iter was given.  */
  unsigned int max_found_given ;	/**< @brief Whether max_found was given.  */
  unsigned int px_given ;	/**< @brief Whether px was given.  */
  unsigned int py_given ;	/**< @brief Whether py was given.  */
  unsigned int kx_given ;	/**< @brief Whether kx was given.  */
  unsigned int ky_given ;	/**< @brief Whether ky was given.  */
  unsigned int m_given ;	/**< @brief Whether m was given.  */
  unsigned int l_given ;	/**< @brief Whether l was given.  */
  unsigned int c_given ;	/**< @brief Whether c was given.  */
  unsigned int r_given ;	/**< @brief Whether r was given.  */

} ;

/** @brief The additional parameters to pass to parser functions */
struct cmdline_parser_params
{
  int override; /**< @brief whether to override possibly already present options (default 0) */
  int initialize; /**< @brief whether to initialize the option structure gengetopt_args_info (default 1) */
  int check_required; /**< @brief whether to check that all required options were provided (default 1) */
  int check_ambiguity; /**< @brief whether to check for options already specified in the option structure gengetopt_args_info (default 0) */
  int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
} ;

/** @brief the purpose string of the program */
extern const char *gengetopt_args_info_purpose;
/** @brief the usage string of the program */
extern const char *gengetopt_args_info_usage;
/** @brief all the lines making the help output */
extern const char *gengetopt_args_info_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser (int argc, char * const *argv,
  struct gengetopt_args_info *args_info);

/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use cmdline_parser_ext() instead
 */
int cmdline_parser2 (int argc, char * const *argv,
  struct gengetopt_args_info *args_info,
  int override, int initialize, int check_required);

/**
 * The command line parser (version with additional parameters)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_ext (int argc, char * const *argv,
  struct gengetopt_args_info *args_info,
  struct cmdline_parser_params *params);

/**
 * Save the contents of the option struct into an already open FILE stream.
 * @param outfile the stream where to dump options
 * @param args_info the option struct to dump
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_dump(FILE *outfile,
  struct gengetopt_args_info *args_info);

/**
 * Save the contents of the option struct into a (text) file.
 * This file can be read by the config file parser (if generated by gengetopt)
 * @param filename the file where to save
 * @param args_info the option struct to save
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_file_save(const char *filename,
  struct gengetopt_args_info *args_info);

/**
 * Print the help
 */
void cmdline_parser_print_help(void);
/**
 * Print the version
 */
void cmdline_parser_print_version(void);

/**
 * Initializes all the fields a cmdline_parser_params structure 
 * to their default values
 * @param params the structure to initialize
 */
void cmdline_parser_params_init(struct cmdline_parser_params *params);

/**
 * Allocates dynamically a cmdline_parser_params structure and initializes
 * all its fields to their default values
 * @return the created and initialized cmdline_parser_params structure
 */
struct cmdline_parser_params *cmdline_parser_params_create(void);

/**
 * Initializes the passed gengetopt_args_info structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void cmdline_parser_init (struct gengetopt_args_info *args_info);
/**
 * Deallocates the string fields of the gengetopt_args_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void cmdline_parser_free (struct gengetopt_args_info *args_info);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int cmdline_parser_required (struct gengetopt_args_info *args_info,
  const char *prog_name);


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* CMDLINE_H */
