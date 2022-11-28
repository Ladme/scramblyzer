// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#include "general.h"
#include "positions.h"
#include "rate.h"

// frequency of printing during the calculation
static const int PROGRESS_FREQ = 10000;

/*! @brief Prints supported flags and arguments of this module */
void print_usage_positions(void)
{
    printf("\nValid OPTIONS for the positions module:\n");
    printf("-h               print this message and exit\n");
    printf("-c STRING        gro file to read\n");
    printf("-f STRING        xtc file to read\n");
    printf("-n STRING        ndx file to read (optional, default: index.ndx)\n");
    printf("-o STRING        output file name (default: positions.xvg)\n");
    printf("-p STRING        selection of lipid head identifiers (default: name PO4)\n");
    printf("-t FLOAT         time interval between analyzed frames [in ns] (default: 1.0)");
    printf("\n");
}

int get_arguments_positions(
        const int argc, 
        char **argv,
        char **gro_file,
        char **xtc_file,
        char **ndx_file,
        char **output_file,
        char **phosphates,
        float *dt) 
{
    // we can reuse the get_arguments_rate function
    return get_arguments_rate(argc, argv, gro_file, xtc_file, ndx_file, output_file, phosphates, dt);
}

void print_arguments_positions(
        const char *gro_file,
        const char *xtc_file,
        const char *ndx_file,
        const char *output_file,
        const char *phosphates,
        const float timestep)
{
    printf("Parameters for Lipid Positions Analysis:\n");
    printf(">>> gro file:         %s\n", gro_file);
    printf(">>> xtc file:         %s\n", xtc_file);
    printf(">>> ndx file:         %s\n", ndx_file);
    printf(">>> output file:      %s\n", output_file);
    printf(">>> lipid heads:      %s\n", phosphates);
    printf(">>> time step:        %f ns\n", timestep);
    printf("\n");
}

int calc_lipid_positions(
        const char *input_gro_file,
        const char *input_xtc_file,
        const char *ndx_file,
        const char *output_file,
        const char *head_identifier,
        const float dt)
{
    print_arguments_positions(input_gro_file, input_xtc_file, ndx_file, output_file, head_identifier, dt);

    // read gro file
    system_t *system = load_gro(input_gro_file);
    if (system == NULL) return 1;

    atom_selection_t *all = select_system(system);

     // read ndx file
    dict_t *ndx_groups = read_ndx(ndx_file, system);

    // get lipid heads
    atom_selection_t *heads = smart_select(all, head_identifier, ndx_groups);
    if (heads == NULL || heads->n_atoms == 0) {
        fprintf(stderr, "No lipid headgroups ('%s') found.\n", head_identifier);
        dict_destroy(ndx_groups);
        free(all);
        free(heads);
        free(system);  
        return 1;    
    }

    dict_destroy(ndx_groups);
    free(all);

    // open output file
    FILE *output = fopen(output_file, "w");
    if (output == NULL) {
        fprintf(stderr, "Could not open output file %s\n", output_file);
        free(heads);
        free(system);
        return 1;
    }

    // write header for the output file
    fprintf(output, "# Generated with Scramblyzer Positions from file %s\n", input_xtc_file);
    fprintf(output, "@    title \"Positions of lipid heads in time\"\n");
    fprintf(output, "@    xaxis label \"time [ns]\"\n");
    fprintf(output, "@    yaxis label \"z-coordinate [nm]\"\n");
    for (size_t i = 0; i < heads->n_atoms; ++i) {
        fprintf(output, "@    s%zu legend \"index %d\"\n", i, heads->atoms[i]->atom_number);
    }

    // open xtc file for reading
    XDRFILE *xtc = xdrfile_open(input_xtc_file, "r");
    if (xtc == NULL) {
        fprintf(stderr, "File %s could not be read as an xtc file.\n", input_xtc_file);
        free(heads);
        free(system);
        fclose(output);
        return 1;
    }

    // check that the gro file and the xtc file match each other
    if (!validate_xtc(input_xtc_file, (int) system->n_atoms)) {
        fprintf(stderr, "Number of atoms in %s does not match %s.\n", input_xtc_file, input_gro_file);
        free(heads);
        free(system);
        xdrfile_close(xtc);
        fclose(output);
        return 1;
    }

    while (read_xtc_step(xtc, system) == 0) {
        // print info about the progress of reading and writing
        if ((int) system->time % PROGRESS_FREQ == 0) {
            printf("Step: %d. Time: %.0f ps\r", system->step, system->time);
            fflush(stdout);
        }

        if ((int) system->time % (int) roundf((dt * 1000)) != 0) {
            continue;
        }

        // loop through heads, get their positions and write them into output file
        fprintf(output, "%f ", system->time / 1000.0);
        for (size_t i = 0; i < heads->n_atoms; ++i) {
            fprintf(output, "%f ", heads->atoms[i]->position[2]);
        }
        fprintf(output, "\n");
    }

    free(heads);
    free(system);
    xdrfile_close(xtc);
    fclose(output);
    return 0;
}