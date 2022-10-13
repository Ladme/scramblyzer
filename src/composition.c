// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#include "general.h"
#include "composition.h"

// frequency of printing during the calculation
static const int PROGRESS_FREQ = 10000;

/*! @brief Identifier for all lipids in the classify_lipids() dictionaries
 * 
 * @paragraph Details
 * @@ is used to avoid any potential overlap with real lipid name.
 */
static const char ALL_LIPIDS_IDENTIFIER[50] = "@@TOTAL@@";

/*! @brief Assigns all lipids from the lipids dictionary into upper and lower leaflets. */
static void classify_lipids(
        lipid_composition_t *composition,
        vec_t membrane_center,
        box_t box,
        dict_t *upper_leaflet,
        dict_t *lower_leaflet)
{
    size_t total_upper = 0, total_lower = 0;
    // loop through all available lipid names
    for (size_t i = 0; i < composition->n_lipid_types; ++i) {
        atom_selection_t *selection = *((atom_selection_t **) dict_get(composition->lipids_dictionary, composition->lipid_types[i]));
        
        size_t upper = 0, lower = 0;
        // loop through the heads of the selection
        for (size_t j = 0; j < selection->n_atoms; ++j) {
            if (distance1D(selection->atoms[j]->position, membrane_center, z, box) > 0) {
                ++upper;
            } else {
                ++lower;
            }
        }

        total_upper += upper;
        total_lower += lower;
        // add the calculated numbers to the dictionaries
        dict_set(upper_leaflet, composition->lipid_types[i], &upper, sizeof(size_t));
        dict_set(lower_leaflet, composition->lipid_types[i], &lower, sizeof(size_t));
    }

    dict_set(upper_leaflet, ALL_LIPIDS_IDENTIFIER, &total_upper, sizeof(size_t));
    dict_set(lower_leaflet, ALL_LIPIDS_IDENTIFIER, &total_lower, sizeof(size_t));
}

/*! @brief Prints supported flags and arguments of this module */
void print_usage_composition(void)
{
    printf("\nValid OPTIONS for the composition module:\n");
    printf("-h               print this message and exit\n");
    printf("-c STRING        gro file to read\n");
    printf("-f STRING        xtc file to read (optional)\n");
    printf("-n STRING        ndx file to read (optional, default: index.ndx)\n");
    printf("-o STRING        output file name (default: composition.xvg)\n");
    printf("-p STRING        selection of lipid head identifiers (default: name PO4)\n");
    printf("-t FLOAT         time interval between analyzed trajectory frames in ns (default: 1.0)\n");
    printf("\n");
}

int get_arguments_composition(
        const int argc, 
        char **argv,
        char **gro_file,
        char **xtc_file,
        char **ndx_file,
        char **output_file,
        char **phosphates,
        float *dt) 
{
    int gro_specified = 0;

    int opt = 0;
    while((opt = getopt(argc - 1, argv + 1, "c:f:n:o:p:t:h")) != -1) {
        switch (opt) {
        // help
        case 'h':
            return 1;
        // gro file to read
        case 'c':
            *gro_file = optarg;
            gro_specified = 1;
            break;
        // xtc file to read
        case 'f':
            *xtc_file = optarg;
            break;
        // ndx file
        case 'n':
            *ndx_file = optarg;
            break;
        // output file name
        case 'o':
            *output_file = optarg;
            break;
        // phosphates identifier
        case 'p':
            *phosphates = optarg;
            break;
        // dt (time precision of the analysis)
        case 't':
            sscanf(optarg, "%f", dt);
            if (*dt <= 0) {
                fprintf(stderr, "dt must be positive.\n");
                return 1;
            }
            break;
        default:
            //fprintf(stderr, "Unknown command line option: %c.\n", opt);
            return 1;
        }
    }

    if (!gro_specified) {
        fprintf(stderr, "Gro file must always be supplied.\n");
        return 1;
    }
    return 0;
}

/* Prints arguments that the program will use for the calculation. */
void print_arguments_composition(
        const char *gro_file,
        const char *xtc_file,
        const char *ndx_file,
        const char *output_file,
        const char *phosphates,
        const float timestep)
{
    printf("Parameters for Composition Analysis:\n");
    printf(">>> gro file:         %s\n", gro_file);
    printf(">>> xtc file:         %s\n", xtc_file);
    printf(">>> ndx file:         %s\n", ndx_file);
    printf(">>> output file:      %s\n", output_file);
    printf(">>> lipid heads:      %s\n", phosphates);
    printf(">>> time step:        %f ns\n", timestep);
    printf("\n");
}

int calc_lipid_composition(
        const char *input_gro_file,
        const char *input_xtc_file,
        const char *ndx_file,
        const char *output_file,
        const char *head_identifier,
        const float dt)
{
    if (input_xtc_file != NULL) {
        print_arguments_composition(input_gro_file, input_xtc_file, ndx_file, output_file, head_identifier, dt);
    }

    // read gro file
    system_t *system = load_gro(input_gro_file);
    if (system == NULL) return 1;

    // read ndx file
    dict_t *ndx_groups = read_ndx(ndx_file, system);

    // get lipids present in the system
    lipid_composition_t *composition = get_lipid_composition(system, head_identifier, ndx_groups);
    if (composition == NULL) {
        free(system);
        dict_destroy(ndx_groups);
        return 1;
    }

    dict_destroy(ndx_groups);

    // if there are no lipids
    if (composition->n_lipid_types < 1) {
        fprintf(stderr, "No usable lipids detected.\n");
        lipid_composition_destroy(composition);
        free(system);
        return 1;
    }

    // if there is no xtc file, just analyze gro file and print to stdout
    if (input_xtc_file == NULL) {
        // get center of geometry of the membrane
        vec_t membrane_center = {0.0};
        center_of_geometry(composition->all_lipid_atoms, membrane_center, system->box);

        dict_t *upper_leaflet = dict_create();
        dict_t *lower_leaflet = dict_create();
        classify_lipids(composition, membrane_center, system->box, upper_leaflet, lower_leaflet);

        printf("Lipid | Upper | Lower | Full \n");
        for (size_t i = 0; i < composition->n_lipid_types; ++i) {
            size_t upper = *((size_t *) dict_get(upper_leaflet, composition->lipid_types[i]));
            size_t lower = *((size_t *) dict_get(lower_leaflet, composition->lipid_types[i]));
            printf("%-5s | %-5zu | %-5zu | %-5zu\n", composition->lipid_types[i], upper, lower, upper + lower);
        }
        // if there are 2 or more lipid types, also print TOTAL number of lipids
        if (composition->n_lipid_types > 1) {
            size_t total_upper = *((size_t *) dict_get(upper_leaflet, ALL_LIPIDS_IDENTIFIER));
            size_t total_lower = *((size_t *) dict_get(lower_leaflet, ALL_LIPIDS_IDENTIFIER));
            printf("-----------------------------\n");
            printf("%-5s | %-5zu | %-5zu | %-5zu\n", "TOTAL", total_upper, total_lower, total_upper + total_lower);
        }
        
        lipid_composition_destroy(composition);
        dict_destroy(upper_leaflet);
        dict_destroy(lower_leaflet);
        free(system);

        return 0;
    }

    // open output file
    FILE *output = fopen(output_file, "w");
    if (output == NULL) {
        fprintf(stderr, "Could not open output file %s\n", output_file);
        lipid_composition_destroy(composition);
        free(system);
        return 1;
    }

    // write header for the output file
    fprintf(output, "# Generated with Scramblyzer Composition from file %s\n", input_xtc_file);
    fprintf(output, "@    title \"Membrane composition in time\"\n");
    fprintf(output, "@    xaxis label \"time [ns]\"\n");
    fprintf(output, "@    yaxis label \"number of lipids\"\n");
    for (size_t i = 0; i < composition->n_lipid_types + 1; ++i) {
        // don't print TOTAL if there is only one lipid species
        if (composition->n_lipid_types < 2 && i == composition->n_lipid_types) break; 

        char *name = NULL;
        if (i < composition->n_lipid_types) {
            name = composition->lipid_types[i];
        } else {
            name = "TOTAL";
        }

        fprintf(output, "@    s%zu legend \"%s_upper\"\n", i * 3, name);
        fprintf(output, "@    s%zu legend \"%s_lower\"\n", i * 3 + 1, name);
        fprintf(output, "@    s%zu legend \"%s_full\"\n", i * 3 + 2, name);
    }

    fprintf(output, "@TYPE xy\n");

    // open xtc file for reading
    XDRFILE *xtc = xdrfile_open(input_xtc_file, "r");
    if (xtc == NULL) {
        fprintf(stderr, "File %s could not be read as an xtc file.\n", input_xtc_file);
        lipid_composition_destroy(composition);
        free(system);
        fclose(output);
        return 1;
    }

    // check that the gro file and the xtc file match each other
    if (!validate_xtc(input_xtc_file, (int) system->n_atoms)) {
        fprintf(stderr, "Number of atoms in %s does not match %s.\n", input_xtc_file, input_gro_file);
        lipid_composition_destroy(composition);
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

        // get center of geometry of the membrane
        vec_t membrane_center = {0.0};
        center_of_geometry(composition->all_lipid_atoms, membrane_center, system->box);

        dict_t *upper_leaflet = dict_create();
        dict_t *lower_leaflet = dict_create();
        classify_lipids(composition, membrane_center, system->box, upper_leaflet, lower_leaflet);

        fprintf(output, "%f     ", system->time / 1000.0);
        for (size_t i = 0; i < composition->n_lipid_types; ++i) {
            size_t upper = *((size_t *) dict_get(upper_leaflet, composition->lipid_types[i]));
            size_t lower = *((size_t *) dict_get(lower_leaflet, composition->lipid_types[i]));
            fprintf(output, "%zu      %zu      %zu      ", upper, lower, upper + lower);
        }

        if (composition->n_lipid_types > 1) {
            size_t total_upper = *((size_t *) dict_get(upper_leaflet, ALL_LIPIDS_IDENTIFIER));
            size_t total_lower = *((size_t *) dict_get(lower_leaflet, ALL_LIPIDS_IDENTIFIER));
            fprintf(output, "%zu      %zu      %zu      ", total_upper, total_lower, total_upper + total_lower);
        }
        fprintf(output, "\n");


        dict_destroy(upper_leaflet);
        dict_destroy(lower_leaflet);

    }

    printf("\nOutput file %s written.\n", output_file);

    lipid_composition_destroy(composition);
    free(system);
    fclose(output);
    xdrfile_close(xtc);

    return 0;
}
