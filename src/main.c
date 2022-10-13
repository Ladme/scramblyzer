// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#include <stdio.h>
#include <string.h>
#include "composition.h"
#include "rate.h"
#include "flipflops.h"

const char VERSION[] = "v2022/07/17";


void print_usage(const char *program_name)
{
    printf("Usage: %s MODULE OPTIONS\n", program_name);
    printf("\nMODULES\n");
    printf("composition      calculates lipid composition of a membrane\n");
    //printf("positions        calculates position of each lipid phosphate in time [NOT IMPLEMENTED]\n");
    printf("rate             calculates percentage of scrambled lipids in time\n");
    printf("flipflops        calculates the number of flip-flop events\n");
    printf("\n");
}

int main(int argc, char **argv)
{
    printf("\n");
    if (argc <= 1) {
        fprintf(stderr, "Module not provided.\n");
        print_usage(argv[0]);
        return 1;
    }

    if (!strcmp(argv[1], "composition")) {
        char *gro_file = NULL;
        char *xtc_file = NULL;
        char *ndx_file = "index.ndx";
        char *output_file = "composition.xvg";
        char *phosphates = "name PO4";
        float dt = 1.0;

        if (get_arguments_composition(argc, argv, &gro_file, &xtc_file, &ndx_file, &output_file, &phosphates, &dt) != 0) {
            print_usage_composition();
            return 1;
        }

        //printf("\n>>> Lipid Composition Analysis by Scramblyzer %s <<<\n\n", VERSION);
        calc_lipid_composition(gro_file, xtc_file, ndx_file, output_file, phosphates, dt);
    
    } else if (!strcmp(argv[1], "rate")) {
        char *gro_file = NULL;
        char *xtc_file = NULL;
        char *ndx_file = "index.ndx";
        char *output_file = "rate.xvg";
        char *phosphates = "name PO4";
        float dt = 10.0;

        if (get_arguments_rate(argc, argv, &gro_file, &xtc_file, &ndx_file, &output_file, &phosphates, &dt) != 0) {
            print_usage_rate();
            return 1;
        }

        calc_scrambling_rate(gro_file, xtc_file, ndx_file, output_file, phosphates, dt);

    } else if (!strcmp(argv[1], "flipflops")) {
        char *gro_file = NULL;
        char *xtc_file = NULL;
        char *ndx_file = "index.ndx";
        char *phosphates = "name PO4";
        float spatial_limit = 1.5;
        int temporal_limit = 10;

        if (get_arguments_flipflops(argc, argv, &gro_file, &xtc_file, &ndx_file, &phosphates, &spatial_limit, &temporal_limit) != 0) {
            print_usage_flipflops();
            return 1;
        }

        calc_lipid_flipflops(gro_file, xtc_file, ndx_file, phosphates, spatial_limit, temporal_limit);
    
    } else if (!strcmp(argv[1], "-h")) {
        print_usage(argv[0]);
    } else {
        fprintf(stderr, "Unknown module %s\n", argv[1]);
        print_usage(argv[0]);
    }

    /*if (get_arguments(argc, argv, &gro_file, &xtc_file, &ndx_file, &output_file, &reference_atoms, &skip) != 0) {
        print_usage(argv[0]);
        return 1;
    }*/

    printf("\n");

    return 0;

    
}