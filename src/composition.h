// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#ifndef COMPOSITION_H
#define COMPOSITION_H

#include <groan.h>
#include <unistd.h>

/*! @brief Prints information about the supported command line arguments for this module.*/
void print_usage_composition(void);


/*! @brief Parses command line arguments for the composition module.
 * 
 * @return Zero, if parsing has been successful. Else returns non-zero.
 */
int get_arguments_composition(
        const int argc, 
        char **argv,
        char **gro_file,
        char **xtc_file,
        char **ndx_file,
        char **output_file,
        char **phosphates,
        float *dt);


/*! @brief Calculates number of lipids of different types either in a gro file or in an xtc trajectory (if provided).
 *
 * @paragraph xtc file is optional
 * If the xtc file is NULL, only gro file will be read and information about the membrane composition will be
 * printed into stdout. If the xtc file is NOT NULL, xvg file 'output_file' is written showing membrane composition
 * in time.
 * 
 * @paragraph Selecting frames for analysis
 * Only some frames will be analyzed based on the value of dt. For instance, if the dt is 10.0 (ns), only frames
 * every 10 ns will be analyzed.
 * 
 * @paragraph What lipids can scramblyzer recognize?
 * Be default scramblyzer is able to recognize all standard lipids of CG force-field Martini 2 (and probably also Martini 3).
 * That includes over a 200 lipid types. Scramblyzer also allows the user to add additional lipid types by writing them
 * into "lipids.txt" file and placing this file in a directory from which scramblyzer is being run.
 * "lipids.txt" must contain only one lipid type per line. Comments must start with '#'.
 * 
 * @param input_gro_file        gro file to read
 * @param input_xtc_file        xtc_file_to_read (not used if NULL)
 * @param output_file           output file (not used if input_xtc_file is NULL)
 * @param head_identifier       name of the atom identifying lipid phosphate/head
 * @param dt                    time interval between analyzed trajectory frames in ns
 * 
 * @return Zero, if the analysis was successful. Else non-zero.
 * 
 */
int calc_lipid_composition(
        const char *input_gro_file,
        const char *input_xtc_file,
        const char *ndx_file,
        const char *output_file,
        const char *head_identifier,
        const float dt);


#endif /* COMPOSITION_H */