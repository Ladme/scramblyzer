// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#include "general.h"
#include "flipflops.h"

// frequency of printing during the calculation
static const int PROGRESS_FREQ = 10000;

/*! @brief Assigns all lipids into membrane leaflets and search for flipflops.*/
static void find_flipflops(
        const lipid_composition_t *composition,
        int **classified_lipids,
        size_t *flipflops_upper_lower,
        size_t *flipflops_lower_upper,
        const vec_t membrane_center,
        const box_t box,
        const float spatial_limit,
        const int time_frames)
{ 
    // loop through all available lipid names
    for (size_t i = 0; i < composition->n_lipid_types; ++i) {
        atom_selection_t *selection = *((atom_selection_t **) dict_get(composition->lipids_dictionary, composition->lipid_types[i]));

        int *assignment = classified_lipids[i];

        // loop through the heads of the selection
        for (size_t j = 0; j < selection->n_atoms; ++j) {
            float dist = distance1D(selection->atoms[j]->position, membrane_center, z, box);

            // UPPER LEAFLET
            if (dist > spatial_limit) {
                // this means that the lipid is stable in the upper leaflet; don't do anything
                if (assignment[j] > time_frames);
                // this means that the lipid flipped recently from the lower leaflet but has not yet stabilized in the upper leaflet
                else if (assignment[j] > 0) assignment[j]++; 
                // this means that the lipid just now flipped from the lower leaflet in which it was stable
                else if (assignment[j] <= -time_frames) assignment[j] = 1;
                // this means that the lipid moved here from the lower leaflet but it was not stable in it (no flip-flop)
                else if (assignment[j] < 0) assignment[j] = time_frames + 1;
                // at the start of the analysis
                else if (assignment[j] == 0) assignment[j] = time_frames + 1;

            // INTERMEDIATE UPPER LEAFLET
            } else if (dist > 0) {
                // this means that the lipid is stable in the upper leaflet; don't do anything
                if (assignment[j] > time_frames);
                // this means that the lipid flipped recently from the lower leaflet but has not yet stabilized in the upper leaflet
                else if (assignment[j] > 0) assignment[j]++;
                // this means that the lipid just now flipped from the lower leaflet in which it was stable
                // don't do anything because the lipid must flip to the true UPPER LEAFLET as defined by spatial limit
                else if (assignment[j] <= -time_frames);
                // this means that the lipid moved here from the lower leaflet but it was not stable in it (no flip-flop)
                else if (assignment[j] < 0) assignment[j] = time_frames + 1;
                // at the start of the analysis
                else if (assignment[j] == 0) assignment[j] = time_frames + 1;

            // LOWER LEAFLET
            } else if (dist < -spatial_limit) {
                // this means that the lipid just now flipped from the upper leaflet in which it was stable
                if (assignment[j] >= time_frames) assignment[j] = -1;
                // this means that the lipid moved here from the upper leaflet but it was not stable in it (no flip-flop)
                else if (assignment[j] > 0) assignment[j] = -time_frames - 1;
                // this means that the lipid is stable in the lower leaflet; don't do anything
                else if (assignment[j] < -time_frames);
                 // this means that the lipid flipped recently from the upper leaflet but has not yet stabilized in the lower leaflet
                else if (assignment[j] < 0) assignment[j]--;
                // at the start of the analysis
                else if (assignment[j] == 0) assignment[j] = -time_frames - 1;

            // INTERMEDIATE LOWER LEAFLET
            } else if (dist < 0) {
                // this means that the lipid just now flipped from the upper leaflet in which it was stable
                // don't do anything because the lipid must flip to the true LOWER LEAFLET as defined by spatial limit
                if (assignment[j] >= time_frames);
                // this means that the lipid moved here from the upper leaflet but it was not stable in it (no flip-flop)
                else if (assignment[j] > 0) assignment[j] = -time_frames - 1;
                // this means that the lipid is stable in the lower leaflet; don't do anything
                else if (assignment[j] < -time_frames);
                 // this means that the lipid flipped recently from the upper leaflet but has not yet stabilized in the lower leaflet
                else if (assignment[j] < 0) assignment[j]--;
                // at the start of the analysis
                else if (assignment[j] == 0) assignment[j] = -time_frames - 1;
            }

            // if the time_frames number is reached, increase the flip-flop counter
            if (assignment[j] == time_frames && dist > 0) {
                flipflops_lower_upper[i]++;
            } else if (assignment[j] == -time_frames && dist < 0) {
                flipflops_upper_lower[i]++;
            }
        }
    }
}

void print_usage_flipflops(void)
{
    printf("\nValid OPTIONS for the flipflops module:\n");
    printf("-h               print this message and exit\n");
    printf("-c STRING        gro file to read\n");
    printf("-f STRING        xtc file to read\n");
    printf("-n STRING        ndx file to read (optional, default: index.ndx)\n");
    printf("-o STRING        output file (default: positions.xvg)\n");
    printf("-p STRING        selection of lipid head identifiers (default: name PO4)\n");
    printf("-s FLOAT         how far into a leaflet must the head of the lipid move to count as flip-flop [in nm] (default: 1.5)\n");
    printf("-t INTEGER       how long must the lipid stay in a leaflet to count as flip-flop [in ns] (default: 10)\n");
    printf("\n");
}

int get_arguments_flipflops(
        const int argc, 
        char **argv,
        char **gro_file,
        char **xtc_file,
        char **ndx_file,
        char **phosphates,
        float *spatial_limit,
        int *temporal_limit) 
{
    int gro_specified = 0, xtc_specified = 0;

    int opt = 0;
    while((opt = getopt(argc - 1, argv + 1, "c:f:n:p:s:t:h")) != -1) {
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
        // phosphates identifier
        case 'p':
            *phosphates = optarg;
            break;
        // spatial limit
        case 's':
            if (sscanf(optarg, "%f", spatial_limit) != 1) {
                fprintf(stderr, "Could not read spatial limit.\n");
                return 1;
            }

            if (*spatial_limit < 0) {
                fprintf(stderr, "Spatial limit must be non-negative.\n");
                return 1;
            }
            break;
        // time limit
        case 't':
            if (sscanf(optarg, "%d", temporal_limit) != 1) {
                fprintf(stderr, "Could not read temporal limit.\n");
                return 1;
            }

            if (*temporal_limit < 1) {
                fprintf(stderr, "Temporal limit cannot be lower than 1 ns.\n");
                return 1;
            }
            break;
        default:
            //fprintf(stderr, "Unknown command line option: %c.\n", opt);
            return 1;
        }
    }

    if (!gro_specified || !xtc_specified) {
        fprintf(stderr, "Gro file and xtc file must always be supplied.\n");
        return 1;
    }
    return 0;
}

/* Prints arguments that the program will use for the calculation. */
void print_arguments_flipflops(
        const char *gro_file,
        const char *xtc_file,
        const char *ndx_file,
        const char *phosphates,
        const float spatial_limit,
        const int temporal_limit)
{
    printf("Parameters for FlipFlops Analysis:\n");
    printf(">>> gro file:         %s\n", gro_file);
    printf(">>> xtc file:         %s\n", xtc_file);
    printf(">>> ndx file:         %s\n", ndx_file);
    printf(">>> lipid heads:      %s\n", phosphates);
    printf(">>> spatial limit:    %f nm\n", spatial_limit);
    printf(">>> temporal limit:   %d ns\n", temporal_limit);
    printf("\n");
}

int calc_lipid_flipflops(
        const char *input_gro_file,
        const char *input_xtc_file,
        const char *ndx_file,
        const char *head_identifier,
        const float spatial_limit,
        const int temporal_limit)
{
    print_arguments_flipflops(input_gro_file, input_xtc_file, ndx_file, head_identifier, spatial_limit, temporal_limit);

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

    // open xtc file for reading
    XDRFILE *xtc = xdrfile_open(input_xtc_file, "r");
    if (xtc == NULL) {
        fprintf(stderr, "File %s could not be read as an xtc file.\n", input_xtc_file);
        lipid_composition_destroy(composition);
        free(system);
        return 1;
    }

    // check that the gro file and the xtc file match each other
    if (!validate_xtc(input_xtc_file, (int) system->n_atoms)) {
        fprintf(stderr, "Number of atoms in %s does not match %s.\n", input_xtc_file, input_gro_file);
        lipid_composition_destroy(composition);
        free(system);
        xdrfile_close(xtc);
        return 1;
    }

    // create an array for lipid classificiation
    int **classified = calloc(composition->n_lipid_types, sizeof(int *));
    for (size_t i = 0; i < composition->n_lipid_types; ++i) {
        atom_selection_t *selection = *((atom_selection_t **) dict_get(composition->lipids_dictionary, composition->lipid_types[i]));
        classified[i] = calloc(selection->n_atoms, sizeof(int));
    }
    // create arrays for flip-flop
    size_t *flipflops_upper_lower = calloc(composition->n_lipid_types, sizeof(size_t));
    size_t *flipflops_lower_upper = calloc(composition->n_lipid_types, sizeof(size_t));

    float prevtime = -1.0;
    while (read_xtc_step(xtc, system) == 0) {
        // print info about the progress of reading and writing
        if ((int) system->time % PROGRESS_FREQ == 0) {
            printf("Step: %d. Time: %.0f ps\r", system->step, system->time);
            fflush(stdout);
        }

        // only analyze every nanosecond
        if ((int) system->time % 1000 != 0) {
            continue;
        }

        // sanity check of the trajectory
        if (prevtime >= 0 && system->time - prevtime > 1000) {
            fprintf(stderr, "Scramblyzer flipflops expects trajectory time step not to be higher than 1 ns.\n");
            fprintf(stderr, "Times of concern: %f (current), %f (previous)\n", system->time, prevtime);
            for (size_t i = 0; i < composition->n_lipid_types; ++i) {
                free(classified[i]);
            }
            free(classified);

            lipid_composition_destroy(composition);
            free(system);
            free(flipflops_upper_lower);
            free(flipflops_lower_upper);
            xdrfile_close(xtc);
            return 1;

        } else {
            prevtime = system->time;
        }

        // get center of geometry of the membrane
        vec_t membrane_center = {0.0};
        center_of_geometry(composition->all_lipid_atoms, membrane_center, system->box);

        find_flipflops(composition, classified, flipflops_upper_lower, flipflops_lower_upper, membrane_center, system->box, spatial_limit, temporal_limit);
    }

    // printing output
    //printf("Detected flip-flops with spatial limit = %f nm and temporal limit = %d ns:\n", spatial_limit, temporal_limit);
    printf("\n\nLipid | U->L | L->U | All \n");
    size_t total_upper_lower = 0, total_lower_upper = 0;
    for (size_t i = 0; i < composition->n_lipid_types; ++i) {
        total_upper_lower += flipflops_upper_lower[i];
        total_lower_upper += flipflops_lower_upper[i];

        printf("%-5s | %-4zu | %-4zu | %-4zu\n", 
            composition->lipid_types[i], 
            flipflops_upper_lower[i], 
            flipflops_lower_upper[i],
            flipflops_upper_lower[i] + flipflops_lower_upper[i]);
    }
    
    // if there are 2 or more lipid types, also print TOTAL number of lipids
    if (composition->n_lipid_types > 1) {
        printf("-----------------------------\n");
        printf("TOTAL | %-4zu | %-4zu | %-4zu\n", total_upper_lower, total_lower_upper, total_upper_lower + total_lower_upper);
    }


    for (size_t i = 0; i < composition->n_lipid_types; ++i) {
        free(classified[i]);
    }
    free(classified);

    lipid_composition_destroy(composition);
    free(system);
    free(flipflops_upper_lower);
    free(flipflops_lower_upper);
    xdrfile_close(xtc);

    return 0;
}