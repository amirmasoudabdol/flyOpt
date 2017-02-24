/** 
 * @file dSFMT_str_state.h 
 * @brief function to save and load the state of a dSFMT
 *
 * @author Andrea C G Mennucci (Scuola Normale Superiore)
 *
 * Copyright (C) 2010 Andrea C G Mennucci
 *
 * The new BSD License is applied to this software, see LICENSE.txt
 */

/**
 * This function returns a string that represents the state of the dSFMT.
 * The string is allocated and should be free-d after use.
 *
 * @param dsfmt dsfmt state vector.
 * @param prefix a prefix to start all lines; if it is NULL, it is set to "dsfmt_"
 */
char *dsfmt_state_to_str(dsfmt_t *dsfmt, char *prefix);

/**
 * This function reads a NULL terminated list of string that represents the state of the dSFMT, and fills the state.
 * It returns NULL if OK, or a string explaining the error, if any. (The error should be freed after use)
 *
 * @param dsfmt dsfmt state vector to be filled.
 * @param strlist the NULL terminated list of strings encoding the state
 * @param prefix the prefix that starts all lines; if it is NULL, it is set to "dsfmt_"
 */
char *dsfmt_strlist_to_state(dsfmt_t *dsfmt, char **strlist, char *prefix);

/**
 * This function reads a string that represents the state of the dSFMT, and fills the state.
 * It returns NULL if OK, or a string explaining the error, if any. (It should be freed after use)
 *
 * @param dsfmt dsfmt state vector.
 * @param origstr the string encoding the state
 * @param prefix the prefix that starts all lines; if it is NULL, it is set to "dsfmt_"
 */
char *dsfmt_str_to_state(dsfmt_t *dsfmt, char *origstr, char *prefix);

/**
 * This function reads a file descriptor that represents the state of the dSFMT, and fills the state.
 * It returns NULL if OK, or a string explaining the error, if any. (It should be freed after use)
 *
 * @param dsfmt dsfmt state vector.
 * @param origfile the file encoding the state
 * @param prefix the prefix that starts all lines; if it is NULL, it is set to "dsfmt_"
 */
//static inline int string_is_comment(char *s);

//static inline char *fgets_noncomment(char *s,int l,FILE *f);

char *dsfmt_file_to_state(dsfmt_t *dsfmt, FILE *origfile, char *prefix);
