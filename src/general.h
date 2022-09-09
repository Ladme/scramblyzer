// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#ifndef GENERAL_H
#define GENERAL_H

#include <groan.h>

/*! @brief Identifier for all lipids in the classify_lipids() dictionaries
 * 
 * @paragraph Details
 * @@ is used to avoid any potential overlap with real lipid name.
 */
const char ALL_LIPIDS_IDENTIFIER[50];

/*! @brief Lipid composition of a membrane. See get_lipid_composition() for more details. */
typedef struct lipid_composition {
    atom_selection_t *all_lipid_atoms;
    dict_t *lipids_dictionary;
    char **lipid_types;
    size_t n_lipid_types;
} lipid_composition_t;


/*! @brief Deallocates memory allocated for the array of residue names of lipids in read_lipid_names().
 * 
 * @param lipid_names       pointer to the lipid names array to deallocate
 * @param n_lipid_names     number of lipid names saved in the array
 */
void lipid_names_destroy(char **lipid_names, const size_t n_lipid_names);


/*! @brief Gets residue names of lipids both from the default list and from the file lipids.txt (if present).
 * 
 * @paragraph Concerning deallocation
 * The returned pointer must later be deallocated using destroy_lipids().
 * 
 * @paragraph Concerning n_lipids argument
 * n_lipids can be pointer to zero. The number of lipid names loaded is saved into memory this pointer points at.
 * 
 * @param n_lipid_names          pointer to the number of lipid names read
 * 
 * @return Pointer to an array of lipid names. NULL in case of an error.
 */
char **read_lipid_names(size_t *n_lipid_names);


/*! @brief Loops through lipid selections from the lipids dictionary and deallocates them. */
void deallocate_lipid_types(dict_t *lipids_dictionary, char **lipid_names, size_t n_lipid_names);


/*! @brief Get lipid composition of a membrane. 
 *
 * @paragraph Lipid composition structure
 * This returns a pointer to lipid_composition structure containing the following information:
 * a) atom selection of all atoms that were identified as belonging to lipids (all_lipid_atoms),
 * b) dictionary of lipid name -> atom selection pairs; atom selection contains one headgroup atom for each selected lipid (lipids dictionary)
 * c) an array of lipid types present in the system (lipid_types)
 * d) number of lipid types present in the system (n_lipid_types)
 * 
 * @paragraph Note on deallocation
 * The memory pointed at by the returned pointer must be deallocated using lipid_composition_destroy().
 * 
 * @param system            pointer to system_t structure that should be read
 * @param head_identifier   string containing the names of atoms identifying lipid heads
 *
 * @return Pointer to lipid_composition structure. NULL in case of an error.
 */
lipid_composition_t *get_lipid_composition(
        system_t *system,
        const char *head_identifier);


/*! @brief Deallocates memory for lipid_composition_t strucutre */
void lipid_composition_destroy(lipid_composition_t *composition);

#endif /* GENERAL_H */