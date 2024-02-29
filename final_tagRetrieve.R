# Functions for process every cell and combine to a search term:
# Functions for extract terms from cell contents:

# Process function for Biomarker/Common/Scientific/Other names:
process_terms <- function(cell_content) {
  # Extract terms within parentheses
  terms_in_parentheses <- gregexpr("\\(([^)]+)\\)", cell_content, perl = TRUE)
  additional_terms <- regmatches(cell_content, terms_in_parentheses)[[1]]
  additional_terms <- gsub("\\(|\\)", "", additional_terms) # Remove parentheses
  
  # Remove parentheses and content within from the original string
  cell_content <- gsub("\\(.*?\\)", "", cell_content)
  
  # Split by semicolon 
  terms <- unlist(strsplit(cell_content, ";\\s*"))
  terms <- tolower(trimws(terms))
  additional_terms <- unlist(strsplit(additional_terms, ";\\s*"))
  additional_terms <- tolower(trimws(additional_terms))
  
  unique(c(terms, additional_terms))
}

# Process function for Mesh terms:
process_mesh <- function(cell_content) {
  cell_content <- as.character(cell_content)
  # Split by '//' 
  terms <- unlist(strsplit(cell_content, "//\\s*"))
  terms <- trimws(terms)
  
  unique(terms)
}

# Process function for entry terms:
process_entry_terms <- function(cell_content) {
  cell_content <- as.character(cell_content)
  # Split by new line
  terms <- unlist(strsplit(cell_content, "\n|//\\s*"))
  terms <- tolower(trimws(terms))
  unique(terms)
}

# Function to format terms in a cell
format_terms <- function(terms) {
  #unique_terms <- unique(terms) # Remove duplicated terms in every cell
  formatted_terms <- sapply(terms, function(term) paste0('"', term, '"'))
  paste(formatted_terms, collapse = "\n")
}

# Process and format terms for each row vector
# Since the chunks are smaller than before, we included Mesh/ Biomarker/ Scientific names/ Entry terms, some terms in common or other names are too broad
process_row_all_vector <- function(row_vector) {
  tibble(
    pathogen = row_vector["Biomarker"],
    #mesh_new = list(format_terms(process_mesh(unlist(strsplit(row_vector["MeSH terms"], "\n"))))),
    MeSH_new = list(format_terms(na.omit(unique(process_mesh(unlist(strsplit(row_vector["MeSH terms"], "\n"))))))),
    Entry_new = list(format_terms(na.omit(unique(c(
      process_terms(unlist(strsplit(row_vector["Biomarker"], "\n"))),
      #process_terms(unlist(strsplit(row_vector["Common Name of Disease"], "\n"))),
      process_terms(unlist(strsplit(row_vector["Scientific Name of Pathogen"], "\n"))),
      #process_terms(unlist(strsplit(row_vector["Other names"], "\n"))),
      process_entry_terms(unlist(strsplit(row_vector["Entry terms"], "\n")))
    )))))
  )
}

# Function to create search term from a row
create_search_term <- function(mesh_terms, entry_terms) {
  # Remove NA and empty strings
  mesh_terms <- mesh_terms[!is.na(mesh_terms) & mesh_terms != ""]
  entry_terms <- entry_terms[!is.na(entry_terms) & entry_terms != ""]
  
  # Create search terms with [Mesh] and [tw] if they are not empty
  mesh_search_with_Mesh <- if (length(mesh_terms) > 0) paste(paste0(mesh_terms, "[Mesh]"), collapse = " OR ") else NULL
  mesh_search_with_tw <- if (length(mesh_terms) > 0) paste(paste0(mesh_terms, "[tw]"), collapse = " OR ") else NULL
  entry_search_with_tw <- if (length(entry_terms) > 0) paste(paste0(entry_terms, "[tw]"), collapse = " OR ") else NULL
  
  # Combine mesh and entry searches with OR, if they exist
  full_search <- paste(na.omit(c(mesh_search_with_Mesh, mesh_search_with_tw, entry_search_with_tw)), collapse = " OR ")
  return(full_search)
}

################################################################################
#Function for retrieve specific tags
retrieve_tag <- function(tag1 = NULL, tag1_subgroup = NULL, tag2 = NULL, pathogen_name = NULL) { 
  # tag_cols for Tag1/Tag1 subgroup/Tag2/pathogens, tag_names for the specific tag names, add "" when using the function
  tag_name <- ""
  # Filter the targeted tags/pathogens
  
  # Initialize an empty vector to collect tag components
  tag_components <- c()
  
  # Filter and collect tag1 if not null
  if (!is.null(tag1)) {
    tag1 <- as.vector(tag1)
    tag_data <- tag_data %>% filter(`Tag1` %in% tag1)
    tag_components <- c(tag_components, gsub(" ", "_", tag1))
  }
  
  # Filter and collect tag1_subgroup if not null
  if (!is.null(tag1_subgroup)) {
    tag1_subgroup <- as.vector(tag1_subgroup) 
    tag_data <- tag_data %>% filter(`Tag1 subgroup` %in% tag1_subgroup)
    tag_components <- c(tag_components, gsub(" ", "_", tag1_subgroup))
  }
  
  # Filter and collect tag2 if not null
  if (!is.null(tag2)) {
    tag2 <- as.vector(tag2)
    tag_data <- tag_data %>% filter(`Tag2` %in% tag2)
    tag_components <- c(tag_components, gsub(" ", "_", tag2))
  }
  
  # Filter and collect pathogen_name if not null
  if (!is.null(pathogen_name)) {
    pathogen_name <- as.vector(pathogen_name)
    tag_data <- tag_data %>% filter(`Pathogen` %in% pathogen_name)
    tag_components <- c(tag_components, gsub(" ", "_", pathogen_name))
  }
  # Combine all components to form the tag_name
  tag_names <- paste(unique(tag_components), collapse = "_")
  
  tag_terms <- clean_term %>%
    semi_join(tag_data, by = c("Scientific Name of Pathogen" = "Pathogen"))
  
  # Get a new table for the specific tag with new Mesh and entry terms
  
  # Transfer the dataset to vector
  list_of_tag <- lapply(seq_len(nrow(tag_terms)), function(i) {
    unlist(tag_terms[i, ])
  })
  
  # Process all vectors
  processed_tag <- map_df(list_of_tag, process_row_all_vector)
  
  # Produce csv for search terms ("Mesh new" col with [Mesh] and "Entry_new" col with [tw])
  combined_tag <- left_join( processed_tag, clean_term, by = c( "pathogen" = "Biomarker"))
  combined_tag_formatted <- combined_tag %>%
    mutate(across(where(is.list), ~sapply(., function(x) paste(x, collapse = "\n"))))%>% 
    filter(!is.na(pathogen))
  write.csv(combined_tag_formatted, paste0("/Users/nynn/Library/CloudStorage/OneDrive-JohnsHopkins/Serologic scoping review/Search results/", tag_names, "_search_terms.csv"))
  
  
  # Get the search term
  
  # Apply the function to mesh and entry terms of each row
  tag_search_terms <- lapply(1:nrow(combined_tag_formatted), function(i) {
    mesh_terms <- strsplit(combined_tag_formatted$MeSH_new[i], "\n")[[1]]
    entry_terms <- strsplit(combined_tag_formatted$Entry_new[i], "\n")[[1]]
    create_search_term(mesh_terms, entry_terms)
  })
  
  # Flatten, de-duplicate, and concatenate all search terms
  flat_search_terms <- unique(na.omit(unlist(tag_search_terms)))
  tag_search_term <- paste(flat_search_terms, collapse = " OR ")
  
  # Print and save the final search term
  tag_final_search_term <- paste(Luminex_terms, " AND (", tag_search_term, ")")
  #cat(tag_final_search_term)
  writeLines(tag_final_search_term, paste0("/Users/nynn/Library/CloudStorage/OneDrive-JohnsHopkins/Serologic scoping review/Search results/", tag_names, "_search_terms.txt"))
  
  # Return a list containing both the formatted table and the search term
  return(list(combined_table = combined_tag_formatted, search_term = tag_final_search_term))
}


#To spurce the R script containing the R codes for the specific function
#source("final_tagRetrieve.R")