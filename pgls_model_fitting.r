library(nlme)
library(ape)

#To construct residue-level pgls fitting of lifespan based on RES
run_pgls_analysis <- function(data_file, tree_file, out) {

  #data_file: csv file of residue-level RES
  #tree_file: phylogenetic tree of species in order with specified outgroup based on MEGA alignment
  #out: pgls model to be further processed to obtain results, see processing instructions below
  
  order <- read.csv(data_file, row.names = 1, check.names = FALSE)
  names(order) <- gsub("^X\\.|\\.+$", "", names(order)) 
  
  order_tree <- read.tree(tree_file)
  order_tree <- root(order_tree, outgroup = out, resolve.root = TRUE)
  order_tree <- multi2di(order_tree)
  
  matched_species <- intersect(rownames(order), order_tree$tip.label)
  if (length(matched_species) != nrow(order)) {
    warning("Some species in the data are not in the tree or vice versa. Subsetting the data and tree...")
  }

  order <- order[matched_species, ]
  order_tree <- drop.tip(order_tree, setdiff(order_tree$tip.label, matched_species))  # Tree
  
  order$spp <- rownames(order)

  resi_columns <- grep("resi_", names(order), value = TRUE)
  
  results <- data.frame(Residue = character(),
                        Coefficient = numeric(),
                        P_Value = numeric(),
                        T_Statistic = numeric(),
                        stringsAsFactors = FALSE)
  
  for (resi_col in resi_columns) {
    
    complete_data <- order[complete.cases(order[, c(resi_col, "lifespan")]), ]
    
    spp <- rownames(complete_data)
    corBM <- corBrownian(phy = order_tree, form = ~spp)
    
    if (nrow(complete_data) > 2) {
      model <- nlme::gls(as.formula(paste("lifespan ~", resi_col)), 
                         data = complete_data, 
                         correlation = corBM)
      
      # Extract the coefficient, p-value, and t-statistic from the summary
      model_summary <- summary(model)
      coefficient <- model_summary$tTable[2, "Value"]     
      p_value <- model_summary$tTable[2, "p-value"]     
      t_statistic <- model_summary$tTable[2, "t-value"]  
      
      results <- rbind(results, data.frame(Residue = resi_col,
                                           Coefficient = coefficient,
                                           P_Value = p_value,
                                           T_Statistic = t_statistic,
                                           stringsAsFactors = FALSE))
    } else {
      warning(paste("Not enough data to fit the model for", resi_col))
    }
  }
  return(list(results = results, corBM = corBM))
}

#**Processing pgls model**
# Obtain results using <variable_name_of_analysis_data>$results