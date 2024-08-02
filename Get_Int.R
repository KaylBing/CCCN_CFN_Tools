# Load the STRINGdb package
library(STRINGdb)
#load("/home/mikhail/Code/PTMs_New/Temp/STRINGdb/R/rstring.R")

# Initialize the STRING database object
string_db <- STRINGdb$new( version="12.0", species=9606, score_threshold=0, link_data="detailed", input_directory="")

# Retrieve the proteins from the STRING database
string_proteins <- string_db$get_proteins()

# Check the dimensions of the retrieved proteins
print(dim(string_proteins))  # Expected output: 19566 4

# Read the dataset that you want to combine with the STRING database
filter_db <- read.table("/home/mikhail/Downloads/Stringdbfilter.txt", header = TRUE, sep = "")

# Check the column names of filter_db
print(colnames(filter_db))

# Assuming filter_db contains gene identifiers that match those in the STRING database,
# you may need to map your gene identifiers to STRING identifiers before proceeding.
# Replace 'gene_identifier_column' with the actual column name from your filter_db
mapped_genes <- string_db$map(filter_db, "experimental", removeUnmappedRows = TRUE)

# Print the mapped genes to check
print(head(mapped_genes))

# Retrieve the interactions for the mapped genes
interactions <- string_db$get_interactions(mapped_genes$STRING_id)

# Print the interactions to check
print((interactions))
#write.table(interactions[c(1, 2, 5), ], file = "interactions.txt", sep = "\t", row.names = TRUE, col.names = FALSE)
write.table(interactions, file = "interactions.txt", sep = " ")
# Assuming 'interactions' is your data frame or matrix


# Select and rename the columns to match the desired output
#output <- interactions[, c("protein1", "protein2", "experimental")]
#colnames(output) <- c("from", "to", "experimental_score")

# Print the final output
#print(output)

