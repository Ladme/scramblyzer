// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#include "general.h"
#include "composition.h"

// frequency of printing during the calculation
static const int PROGRESS_FREQ = 10000;

/*! @brief Assign lipids into individual leaflets and save this information into a dictionary. */
static dict_t *create_reference(
        const lipid_composition_t *composition,
        const vec_t membrane_center,
        const box_t box)
{
    dict_t *classified_lipids = dict_create();    
    // loop through all available lipid names
    for (size_t i = 0; i < composition->n_lipid_types; ++i) {
        atom_selection_t *selection = *((atom_selection_t **) dict_get(composition->lipids_dictionary, composition->lipid_types[i]));
        
        // an array in which each number corresponds to one lipid
        // 1 means that the lipid is in the upper leaflet, 0 means that the lipid is in the lower leaflet
        short *selection_ul = calloc(selection->n_atoms, sizeof(short));

        // loop through the heads of the selection
        for (size_t j = 0; j < selection->n_atoms; ++j) {
            if (distance1D(selection->atoms[j]->position, membrane_center, z, box) > 0) {
                selection_ul[j] = 1;
            // else 0, but that is already in the array via calloc
            }
        }

        // assign selection_ul to the dictionary
        dict_set(classified_lipids, composition->lipid_types[i], &selection_ul, sizeof(short *));
    }

    return classified_lipids;
}

/*! @brief Decide how many lipids have been scrambled by comparing their current positions with the reference. Print this information. */
static void classify_lipids(
        FILE *file,
        const lipid_composition_t *composition,
        const dict_t *reference,
        const vec_t membrane_center,
        const box_t box)
{
    // loop through lipid types
    size_t total_scrambled = 0;
    size_t total_lipids = 0;
    for (size_t i = 0; i < composition->n_lipid_types; ++i) {
        atom_selection_t *selection = *((atom_selection_t **) dict_get(composition->lipids_dictionary, composition->lipid_types[i]));
        short *reference_pos = *((short **) dict_get(reference, composition->lipid_types[i]));

        size_t scrambled = 0;
        // loop through all lipids, calculate their position and compare it with their reference position
        for (size_t j = 0; j < selection->n_atoms; ++j) {
            register float dist = distance1D(selection->atoms[j]->position, membrane_center, z, box);

            // lipid was in the lower leaflet, now is in the upper leaflet
            if (reference_pos[j] == 0 && dist > 0) ++scrambled;
            // lipid was in the upper leaflet, now is in the lower leaflet
            if (reference_pos[j] == 1 && dist < 0) ++scrambled;
        }

        fprintf(file, "%f     ", 100.0 * (float) scrambled / selection->n_atoms);

        total_scrambled += scrambled;
        total_lipids += selection->n_atoms;
    }

    if (composition->n_lipid_types > 1) {
        fprintf(file, "%f     ", 100.0 * (float) total_scrambled / total_lipids);
    }

    fprintf(file, "\n");
}


/*! @brief Prints supported flags and arguments of this module */
void print_usage_rate(void)
{
    printf("\nValid OPTIONS for the rate module:\n");
    printf("-h               print this message and exit\n");
    printf("-c STRING        gro file to read\n");
    printf("-f STRING        xtc file to read\n");
    printf("-n STRING        ndx file to read (optional, default: index.ndx)\n");
    printf("-o STRING        output file name (default: rate.xvg)\n");
    printf("-p STRING        selection of lipid head identifiers (default: name PO4)\n");
    printf("-t FLOAT         time interval between analyzed trajectory frames in ns (default: 10.0)\n");
    printf("\n");
}

int get_arguments_rate(
        const int argc, 
        char **argv,
        char **gro_file,
        char **xtc_file,
        char **ndx_file,
        char **output_file,
        char **phosphates,
        float *dt) 
{
    int gro_specified = 0, xtc_specified = 0;

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
            xtc_specified = 1;
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

    if (!gro_specified || !xtc_specified) {
        fprintf(stderr, "Gro and xtc file must always be supplied.\n");
        return 1;
    }
    return 0;
}

/* Prints arguments that the program will use for the calculation. */
void print_arguments_rate(
        const char *gro_file,
        const char *xtc_file,
        const char *ndx_file,
        const char *output_file,
        const char *phosphates,
        const float timestep)
{
    printf("Parameters for Scrambling Rate Analysis:\n");
    printf(">>> gro file:         %s\n", gro_file);
    printf(">>> xtc file:         %s\n", xtc_file);
    printf(">>> ndx file:         %s\n", ndx_file);
    printf(">>> output file:      %s\n", output_file);
    printf(">>> lipid heads:      %s\n", phosphates);
    printf(">>> time step:        %f ns\n", timestep);
    printf("\n");
}

int calc_scrambling_rate(
        const char *input_gro_file,
        const char *input_xtc_file,
        const char *ndx_file,
        const char *output_file,
        const char *head_identifier,
        const float dt)
{
    print_arguments_rate(input_gro_file, input_xtc_file, ndx_file, output_file, head_identifier, dt);

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

    // open output file
    FILE *output = fopen(output_file, "w");
    if (output == NULL) {
        fprintf(stderr, "Could not open output file %s\n", output_file);
        lipid_composition_destroy(composition);
        free(system);
        return 1;
    }

    // write header for the output file
    fprintf(output, "# Generated with Scramblyzer Rate from file %s\n", input_xtc_file);
    fprintf(output, "@    title \"Percentage of scrambled lipids in time\"\n");
    fprintf(output, "@    xaxis label \"time [ns]\"\n");
    fprintf(output, "@    yaxis label \"scrambled lipids [%%]\"\n");
    for (size_t i = 0; i < composition->n_lipid_types + 1; ++i) {
        // don't print TOTAL if there is only one lipid species
        if (composition->n_lipid_types < 2 && i == composition->n_lipid_types) break; 

        char *name = NULL;
        if (i < composition->n_lipid_types) {
            name = composition->lipid_types[i];
        } else {
            name = "TOTAL";
        }

        fprintf(output, "@    s%zu legend \"%s\"\n", i, name);
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

    int frame = 0;
    dict_t *reference = NULL;
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

        // if this is the first frame of the trajectory, create reference classification of lipids
        fprintf(output, "%f     ", system->time / 1000.0);
        if (frame == 0) {
            reference = create_reference(composition, membrane_center, system->box);
            for (size_t i = 0; i < composition->n_lipid_types; ++i) {
                fprintf(output, "0.0        ");
            }
            // total number of scrambled lipids
            if (composition->n_lipid_types > 1) fprintf(output, "0.0");
            fprintf(output, "\n");
            ++frame;
            continue;
        }

        // classify lipids in the current frame
        classify_lipids(output, composition, reference, membrane_center, system->box);
        ++frame;
    }

    printf("\nOutput file %s written.\n", output_file);

    // deallocate memory for the reference dictionary
    for (size_t i = 0; i < composition->n_lipid_types; ++i) {
        short *reference_pos = *((short **) dict_get(reference, composition->lipid_types[i]));
        free(reference_pos);
    }

    dict_destroy(reference);
    lipid_composition_destroy(composition);
    free(system);
    fclose(output);
    xdrfile_close(xtc);


    return 0;
}
