// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#ifndef FLIPFLOPS_H
#define FLIPFLOPS_H

#include <groan.h>
#include <unistd.h>

/*! @brief Prints supported flags and arguments of this module */
void print_usage_flipflops(void);


/*! @brief Parses command line arguments for the flipflops module.
 * 
 * @return Zero, if parsing has been successful. Else returns non-zero.
 */
int get_arguments_flipflops(
        const int argc, 
        char **argv,
        char **gro_file,
        char **xtc_file,
        char **phosphates,
        float *spatial_limit,
        int *temporal_limit);

int calc_lipid_flipflops(
        const char *input_gro_file,
        const char *input_xtc_file,
        const char *head_identifier,
        const float spatial_limit,
        const int temporal_limit);

#endif /* FLIPFLOPS_H */