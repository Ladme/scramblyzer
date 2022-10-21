// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#ifndef POSITIONS_H
#define POSITIONS_H

#include <groan.h>
#include <unistd.h>

/*! @brief Parses command line arguments for the positions module.
 * 
 * @paragraph Implementation details
 * Internally uses the get_arguments_rate function.
 * 
 * @return Zero, if parsing has been successful. Else returns non-zero.
 */
int get_arguments_positions(
        const int argc, 
        char **argv,
        char **gro_file,
        char **xtc_file,
        char **ndx_file,
        char **output_file,
        char **phosphates,
        float *dt);

/*! @brief Prints supported flags and arguments of this module */
void print_usage_positions(void);

/*! @brief Prints arguments that the program will use for the calculation. */
void print_arguments_positions(
        const char *gro_file,
        const char *xtc_file,
        const char *ndx_file,
        const char *output_file,
        const char *phosphates,
        const float timestep);

/*! @brief Analyzes and prints the positions of lipid heads during the simulation. */
int calc_lipid_positions(
        const char *input_gro_file,
        const char *input_xtc_file,
        const char *ndx_file,
        const char *output_file,
        const char *head_identifier,
        const float dt);

#endif /* POSITIONS_H */