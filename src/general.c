// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#include "general.h"

const char ALL_LIPIDS_IDENTIFIER[50] = "@@TOTAL@@";

/*! @brief Maximal length of a line in lipids.txt */
static const size_t MAX_LINE_LENGTH = 1024;
/*! @brief File to read user-defined lipid names/types from */
static const char LIPIDS_TXT[] = "lipids.txt";

void lipid_names_destroy(char **lipid_names, const size_t n_lipid_names)
{
    for (size_t i = 0; i < n_lipid_names; ++i) {
        free(lipid_names[i]);
    }

    free(lipid_names);
}

char **read_lipid_names(size_t *n_lipid_names)
{
    // all Martini lipids
    char default_lipid_names[] = "DAPC DBPC DFPC DGPC DIPC DLPC DNPC DOPC DPPC DRPC DTPC DVPC DXPC DYPC LPPC PAPC PEPC PGPC PIPC POPC PRPC PUPC DAPE DBPE DFPE DGPE DIPE DLPE DNPE DOPE DPPE DRPE DTPE DUPE DVPE DXPE DYPE LPPE PAPE PGPE PIPE POPE PQPE PRPE PUPE DAPS DBPS DFPS DGPS DIPS DLPS DNPS DOPS DPPS DRPS DTPS DUPS DVPS DXPS DYPS LPPS PAPS PGPS PIPS POPS PQPS PRPS PUPS DAPG DBPG DFPG DGPG DIPG DLPG DNPG DOPG DPPG DRPG DTPG DVPG DXPG DYPG JFPG JPPG LPPG OPPG PAPG PGPG PIPG POPG PRPG DAPA DBPA DFPA DGPA DIPA DLPA DNPA DOPA DPPA DRPA DTPA DVPA DXPA DYPA LPPA PAPA PGPA PIPA POPA PRPA PUPA DPP1 DPP2 DPPI PAPI PIPI POP1 POP2 POP3 POPI PUPI PVP1 PVP2 PVP3 PVPI PADG PIDG PODG PUDG PVDG TOG APC CPC IPC LPC OPC PPC TPC UPC VPC BNSM DBSM DPSM DXSM PGSM PNSM POSM PVSM XNSM DPCE DXCE PNCE XNCE DBG1 DPG1 DPG3 DPGS DXG1 DXG3 PNG1 PNG3 XNG1 XNG3 DFGG DFMG DPGG DPMG DPSG FPGG FPMG FPSG OPGG OPMG OPSG CHOA CHOL CHYO BOG DDM DPC EO5 SDS BOLA BOLB CDL0 CDL1 CDL2 CDL DBG3 ERGO HBHT HDPT HHOP HOPR ACA ACN BCA BCN LCA LCN PCA PCN UCA UCN XCA XCN RAMP REMP OANT";

    // get default lipid names
    char **default_names_split = NULL;
    size_t n_default = 0;
    if ((n_default = strsplit(default_lipid_names, &default_names_split, " ")) <= 0) {
        return NULL;
    };
    *n_lipid_names = n_default;

    // allocate new array for all the lipids
    size_t allocated_lipids = n_default;
    char **lipid_names = calloc(allocated_lipids, sizeof(char *));
    if (lipid_names == NULL) {
        *n_lipid_names = 0;
        free(default_names_split);
        return NULL;
    }

    // copy default lipids into the new array
    for (size_t i = 0; i < n_default; ++i) {
        lipid_names[i] = calloc(1, 10);
        strncpy(lipid_names[i], default_names_split[i], 10);
    }
    free(default_names_split);

    // try opening file with the user-defined lipids
    FILE *file = fopen(LIPIDS_TXT, "r");
    if (file == NULL) {
        return lipid_names;
    }

    // read the file with lipids    
    char line[MAX_LINE_LENGTH];
    while (fgets(line, MAX_LINE_LENGTH, file) != NULL) {
        // remove comments
        line[strcspn(line, "#")] = 0;
        // strip line
        strstrip(line);

        if (strlen(line) != 0) {
            // check that the lipid name does not already exist in the array
            int add = 1;
            for (size_t i = 0; i < *n_lipid_names; ++i) {
                if (!strcmp(lipid_names[i], line)) {
                    fprintf(stderr, "Warning. Lipid type %s from %s already exists in the default lipid set.\n\n", line, LIPIDS_TXT);
                    add = 0;
                    break;
                }
            }

            // if the lipid already exists in the array, skip to the next entry
            if (!add) continue;

            if (*n_lipid_names >= allocated_lipids) {
                allocated_lipids *= 2;
                lipid_names = realloc(lipid_names, allocated_lipids * sizeof(char *));
            }

            lipid_names[*n_lipid_names] = calloc(1, 10);
            strncpy(lipid_names[*n_lipid_names], line, 10);
            (*n_lipid_names)++;
        }

    }

    fclose(file);
    return lipid_names;
} 

void deallocate_lipid_types(dict_t *lipids_dictionary, char **lipid_names, size_t n_lipid_names)
{
    for (size_t i = 0; i < n_lipid_names; ++i) {
        atom_selection_t *selection = *((atom_selection_t **) dict_get(lipids_dictionary, lipid_names[i]));

        free(selection);
    }
}


lipid_composition_t *get_lipid_composition(
        system_t *system,
        const char *head_identifier) 
{
    // create lipid composition structure
    lipid_composition_t *composition = calloc(1, sizeof(lipid_composition_t));

    // select all atoms
    atom_selection_t *all = select_system(system);
    // select all head identifiers of lipids
    atom_selection_t *heads = smart_select(all, head_identifier, NULL); //select_atoms(all, head_identifier, &match_atom_name);


    // load lipid names from default and from lipids.txt
    size_t n_lipid_names = 0;
    char **lipid_names = read_lipid_names(&n_lipid_names);
    if (lipid_names == NULL) {
        fprintf(stderr, "Error obtaining lipid names.\n");
        free(all);
        free(heads);
        free(composition);
        return NULL;
    }

    // select lipid atoms corresponding to specific lipid types
    size_t all_lipids_allocated = 64;
    composition->all_lipid_atoms = selection_create(all_lipids_allocated);
    composition->lipids_dictionary = dict_create();
       
    for (size_t i = 0; i < n_lipid_names; ++i) {
        atom_selection_t *lipid_type = select_atoms(all, lipid_names[i], &match_residue_name);

        if (lipid_type == NULL || lipid_type->n_atoms == 0) {
            free(lipid_type);
            continue;
        }

        // add the selection to all lipid atoms
        selection_add(&composition->all_lipid_atoms, &all_lipids_allocated, lipid_type);

        // get only lipid heads of these lipids
        atom_selection_t *lipid_type_heads = selection_intersect(lipid_type, heads);

        // check that this selection is not empty
        if (lipid_type_heads->n_atoms == 0) {
            fprintf(stderr, "Warning. %zu atoms were found for %s lipids but none of these atoms was lipid head identifier %s.\n",
                lipid_type->n_atoms, lipid_names[i], head_identifier);
            fprintf(stderr, "Lipids of type %s will not be included in the analysis.\n\n", lipid_names[i]);
            free(lipid_type);
            free(lipid_type_heads);
            continue;
        }

        // add the selection to dictionary of lipid types
        dict_set(composition->lipids_dictionary, lipid_names[i], &lipid_type_heads, sizeof(atom_selection_t *));
        // deallocate lipid type
        free(lipid_type);
    }

    // deallocate unneeded selections
    free(all);
    free(heads);
    lipid_names_destroy(lipid_names, n_lipid_names);

    // get lipid types that are actually present in the system
    composition->n_lipid_types = dict_keys(composition->lipids_dictionary, &composition->lipid_types);

    return composition;
}

void lipid_composition_destroy(lipid_composition_t *composition)
{
    free(composition->all_lipid_atoms);
    deallocate_lipid_types(composition->lipids_dictionary, composition->lipid_types, composition->n_lipid_types);
    dict_destroy(composition->lipids_dictionary);
    free(composition->lipid_types);
    free(composition);
}