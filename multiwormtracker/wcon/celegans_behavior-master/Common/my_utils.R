#
# My utils - Common general-purpose functions that are expected to be used by (nearly) all scripts
#
# Created on Oct 8th, 2013
# ...
# Modified on July 30th, 2014, to create function 'myRbindIterator'
#

# To clear the workspace:
#   rm(list = ls());

# Note that the order of the packages is important; if not properly set, dplyr seems not to work
library(parallel)
library(plyr)
library(reshape2)
library(ggplot2)
#library(gplots)
#library(grid)
#library(ClassDiscovery)
library(RColorBrewer)
library(dplyr)
#library(doSNOW)
#library(foreach)
#library(doParallel)
library(data.table)

getNamesLoadedPackages = function() {
    SS = search() %>% 
        strsplit(":", fixed = TRUE) %>% 
        sapply(function(S) { if ((length(S) == 2) && (S[[1]] == "package")) { return (S[[2]]) } else { return ("") }})
    return (SS[nchar(SS) > 0])
}

#getNamesAllCorePackages = function() 
#  c("plyr", "reshape2", "ggplot2", "gplots", "grid", "ClassDiscovery", "RColorBrewer", "dplyr", "doSNOW", "foreach", "doParallel")

# To clear the screen
clc = function() cat(rep("\n",50))

# To replace R's 'seq' function. Acutally, it wraps 'seq'
mySeq = function(first_val, last_val, incr = 1) {
  if (incr > 0) {
    if (first_val <= last_val) {
      return (seq(from = first_val, to = last_val, by = incr))
    } else {
      return (c())
    }
  } else if (incr < 0) {
    if (first_val >= last_val) {
      return (seq(from = first_val, to = last_val, by = incr))
    } else {
      return (c())
    }
  } else {
    return (c())
  }
}

#
# Check whether set A is contained in set B
#
isContainedIn = function(A, B) {
  return (isEmpty(setdiff(A, B)))
}

# To remove the names of a list or vector
noNames = function(SS) { names(SS) = NULL; return (SS) }

# To remove the column names of a matrix
noColNames = function(SS) { colnames(SS) = NULL; return (SS) }

# To remove the row names of a matrix
noRowNames = function(SS) { rownames(SS) = NULL; return (SS) }

# Helper function to check whether all elements of a vector/list are unique
areAllValuesUnique = function(v) { return (length(v) == length(unique(v))); }

# Helper MATLAB-like functions
myError = function(msg) { stop(msg) }
errorIf = function(b, msg) { if (b) myError(msg) }

myWarning = function(msg, flag_immediate = TRUE) { warning(msg, call. = FALSE, immediate. = flag_immediate) }
warningIf = function(b, msg, flag_immediate = TRUE) { if (b) myWarning(msg, flag_immediate = flag_immediate) }

isScalar = function(x) { return (length(x) == 1); }
isEmpty = function(x) { return (length(x) == 0); }
isString = function(s) { return (typeof(s) == "character"); }

addNames = function(X, names_str) {
  names(X) = names_str
  return (X)
}

setColNames = function(X, names_str, flag_check_uniq = FALSE) {
  errorIf(flag_check_uniq && !areAllValuesUnique(names_str), "The names for the columns must be unique")
  colnames(X) = names_str
  return (X)
}

setColClass = function(X, col_name, class_name) {
    errorIf(!is.data.frame(X) || is.null(colnames(X)), "")
    errorIf(!isSingleStr(col_name) || !(col_name %in% colnames(X)), "")
    errorIf(!isSingleStr(class_name), "")
    class(X[, col_name]) = class_name
    return (X)
}

applyFun2Cols = function(X, spec) {
    errorIf(!is.data.frame(X) || is.null(colnames(X)), "")
    errorIf(!isNamedList(spec) || !isContainedIn(names(spec), colnames(X)) || 
            !areAllValuesUnique(names(spec))|| any(!sapply(spec, is.function)), "")
    for (ss in names(spec)) {
        fun = spec[[ss]]
        X[, ss] = fun(X[, ss])
    }
    return (X)
}

colsAsChar = function(X, names_target_cols = NULL) {
    errorIf(!is.data.frame(X) || is.null(colnames(X)), "")
    if (is.null(names_target_cols)) names_target_cols = colnames(X)
    errorIf(!isStringVector(names_target_cols) || !areAllValuesUnique(names_target_cols), "")
    return (applyFun2Cols(X, lapplyNamesWrapper(names_target_cols, function(Z) as.character)))
}

charColsAsFactors = function(D) {
    errorIf(!is.data.frame(D), "")
    inds_target_cols = which(sapply(1:ncol(D), function(i_col) class(D[, i_col])) == "character")
    if (!isEmpty(inds_target_cols)) {
        for (i in inds_target_cols) {
            D[, i] = as.factor(D[, i])
        }
    }
    D
}

setRowNames = function(X, names_str, flag_check_uniq = FALSE) {
  errorIf(flag_check_uniq && !areAllValuesUnique(names_str), "The names for the rows must be unique")
  rownames(X) = names_str
  return (X)
}

# "Smart", general, implementation, supports 
# Created on April 21st, 2015
renameStrs = function(SS, naming_map, flag_check_if_unique = TRUE) {
  # Sanity checks
  errorIf(!isStringVector(SS), "")
  
  SS = sapply(SS, function(ss) {
    ind = which(names(naming_map) == ss)
    if (isScalar(ind)) {
      return (naming_map[[ss]])
    } else if (isEmpty(ind)) {
      return (ss)
    } else {
      myError(sprintf("Multiple matches in the renaming map for '%s'", ss))
    }
  })
  errorIf(flag_check_if_unique && !areAllValuesUnique(SS), "Output column names are not unique")
  return (SS)
}

renameCols = function(X, map_col_names) {
  col_names = lapply(colnames(X), function(ss) {  
    if (any(names(map_col_names) == ss)) {
      return (map_col_names[[ss]])
    } else {
      return (ss)
    } 
  })
  errorIf(!areAllValuesUnique(col_names), "Output column names are not unique")
  colnames(X) = col_names
  return (X)
}

renameRows = function(X, map_row_names) {
  row_names = lapply(rownames(X), function(ss) {  
    if (any(names(map_row_names) == ss)) {
      return (map_row_names[[ss]])
    } else {
      return (ss)
    }
  })
  rownames(X) = row_names
  return (X)
}

#
# To rename entries in a list
#
renameListEntries = function(lst, map_names) {
  # Sanity checks
  errorIf(!is.list(lst) || is.null(names(lst)), "The input must be a named list")
  errorIf(!is.vector(map_names) || is.null(names(map_names)), "The map must be a named vector")
  errorIf(!isContainedIn(names(map_names), names(lst)), "Excess names in the map")
  
  new_names = lapply(names(lst), function(ss) {  
    if (any(names(map_names) == ss)) {
      return (map_names[[ss]])
    } else {
      return (ss)
    }
  })
  names(lst) = new_names
  return (lst)
}

#
# To rename the entries in a vector
#
renameVector = function(X, naming_map) {
  # Sanity checks
  errorIf(!is.vector(X) || is.null(names(X)), "The input must be a named list")
  return (mySetNames(X, applyMap(names(X), naming_map, raise_err = FALSE, ret_input_on_err = TRUE), flag_check_uniq_names = TRUE))
}

#
# 
#
getStdString = function(ss) {
  if (nchar(ss)) {
    return (paste(toupper(substr(ss, 1, 1)), substr(ss, 2, nchar(ss)), sep = ""))
  } else {
    return ("")
  }
}

# "Convert" an environment to a list
# Use with caution
env2List = function(Env) {
  errorIf(!is.environment(Env), "The input must be an environment")
  f = ls(Env)
  n_f = length(f)
  lst = list()
  for (i in seq(1, n_f))
    lst[[i]] = Env[[f[[i]]]]
  names(lst) = f
  
  return (lst)
}
 
# Color-blind-friendly palette (see 'http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/')
myDefaultColorPalette = function() {
  color_palette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  return (color_palette)
}

createFolder = function(folder_path, recursive = TRUE, show_warnings = FALSE) { 
  dir.create(folder_path, recursive = recursive, showWarnings = show_warnings)
  return (folder_path)
}

# Helper function to apply a string-matching function ('pmatch') to a list of strings
lapply_pmatch = function(ss, C) { 
  pmatch_fun = function(ss_ref) { return (pmatch(ss, ss_ref, nomatch = -1)); }; 
  return (lapply(C, pmatch_fun)) 
} 
# Return the indices of string in a list ('C') that begin with the given string ('ss')
getIndsPMatchBeginningStr = function(ss, C) { return (which(unlist(lapply_pmatch(ss, C)) == 1)); }

# Determine if two vectors/lists are identical (have the same number of elements AND all elements are equal)
areIdenticalContainers = function (X, Y) {
  return ((length(X) == length(Y)) && (all(X == Y)))
}

#
# Load data from a file into a list (MATLAB-like format)
#
loadRData = function(file_path) {
  data_env = new.env()
  load("transcriptome_data.RData", data_env)
  field_names = ls(data_env)
  n_fields = length(field_names)
  res = list()
  for (i in seq(1, n_fields))
    res[[i]] = get(field_names[[i]], data_env)
  names(res) = field_names
  return (res)
}

# To fix the mess in R of a numeric matrix
# Used when extracting the results
toMatrix = function(A, col_names) { 
  NN = nrow(A); AA = matrix(nrow = NN, ncol = ncol(A))
  for (i in seq(1, NN)) AA[i, ] = unlist(A[i, ])
  rownames(AA) = rownames(A)
  colnames(AA) = col_names
  return (AA)
}

#
# Format date string 
# Example of output: "Oct_16_2013_14_48_50"
#
formatDateStr = function(val = Sys.time()) {
  return (gsub(" ", "_", gsub(":", "", format(val, "%b%d%Y_%X"))))
}

formatDateStr_msec = function(val = Sys.time()) 
    gsub(" ", "_", gsub(":", "", format(val, "%b%d%Y_%H%M%OS6")))

#
# Format filename with a string denoting the data
# For example
#
formatFilenameWithDateStr = function(filename_token_str, sep = "-", val = Sys.time()) {
  return (paste(filename_token_str, formatDateStr(val = val), sep = sep))
}

#
# Check if a given factor is "balanced", i.e., all its levels occur the same number of times
#
isBalancedFactor = function(x) {
  errorIf(typeof(x) != "integer", "The input must be a numeric vector of factors")
  return (isScalar(unique(sapply(unique(x), function(xx) { return (sum(x == xx)) }))))  
}

#
# Format a value of a floating point number as a string to be used in a file name or folder name
# For now, this basically involves replacing the decimal point by another character
#
formatNum4FileFolderName = function(x, fmt_str, new_char = "_") {
  return (gsub("\\.", new_char, sprintf(fmt_str, x)))
}

# map_vals is a list that is used as a map
replaceValuesInColumn = function(Data, col_name, map_vals) {
  errorIf(!any(colnames(Data) == col_name), "Invalid column name")
  errorIf(!areAllValuesUnique(names(map_vals)), "The names in the map must be unique")
  vals_not_in_map = setdiff(Data[, col_name], unlist(names(map_vals)))
  errorIf(!isEmpty(vals_not_in_map), "Incomplete map")
  vals = unlist(Data[, col_name])
  out_vals = unlist(map_vals[vals])
  Data[, col_name] = out_vals
  return (Data)
}

#
# Function for validating the data on predictors
# The input is a data.frame
# Raises an exception if at least one of the checks performed fails
#
validateBalancing = function(In_data) {
  errorIf(!is.data.frame(In_data), "The input data must be stored in a data frame")
  
  I_ancestral = In_data[, "exper_evol"] == "ancestral" # G0, by definition
  I_ngm = In_data[, "exper_evol"] == "ngm" # G50, by definition
  I_sudden_gradual = (In_data[, "exper_evol"] == "sudden") | (In_data[, "exper_evol"] == "gradual")
  I_G50 = (In_data[, "generation"] == 50)
  
  #
  # 'columns_groups' is a list of the names of the columns that define the various "groups"
  #
  coreValidateBalancing = function(IIn_data, columns_groups) {
    n_c = length(columns_groups)
    errorIf(!n_c, "Need to the names of the columns that define the 'groups'")
    if (n_c == 1) {
      Data_groups = as.matrix(IIn_data[, columns_groups])
    } else {
      Data_groups = IIn_data[, columns_groups]
    }
    
    un_group_vals = unique(Data_groups)
    n_un_rows = nrow(un_group_vals) # this is the number of (effective) groups
    
    # Returns a vector of boolean values
    matchRow2Matrix = function(row_data, Matrix_data) {
      return (sapply(1:nrow(Matrix_data), function(row_ii) { return (all(row_data == Matrix_data[row_ii, ])) }))
    }
    
    # Index of the group to which each row belows
    I_group = sapply(1:nrow(Data_groups), function(i_row) { return (which(matchRow2Matrix(Data_groups[i_row, ], un_group_vals))) })
    
    # Make sure number of samples in each (effective) group is the same
    n_elements_per_group = sapply(1:n_un_rows, function (i_group) { return (sum(i_group == I_group)) })
    
    return (isScalar(unique(n_elements_per_group)))
  }
  coreValidateBalancing_Ancestral_NGM_SuddenGradual = function(IIn_data, columns_groups) {
    return ((coreValidateBalancing(IIn_data[I_ancestral, ], columns_groups)) &&
              (coreValidateBalancing(IIn_data[I_ngm, ], columns_groups)) &&
              (coreValidateBalancing(IIn_data[I_sudden_gradual, ], columns_groups)) && 
              (coreValidateBalancing(IIn_data[I_G50, ], columns_groups)))
  }
  
  # {generation, exper_evol, NaCl, block, plate}
  # Brute-force approach for validating
  errorIf(!((coreValidateBalancing(In_data, c("NaCl"))) && 
              (coreValidateBalancing(In_data, c("block"))) && 
              (coreValidateBalancing(In_data, c("plate"))) && 
              (coreValidateBalancing_Ancestral_NGM_SuddenGradual(In_data, c("exper_evol"))) &&
              (coreValidateBalancing_Ancestral_NGM_SuddenGradual(In_data, c("generation"))) &&
              #
              (coreValidateBalancing(In_data, c("generation", "exper_evol"))) &&
              (coreValidateBalancing_Ancestral_NGM_SuddenGradual(In_data, c("generation", "NaCl"))) &&
              (coreValidateBalancing_Ancestral_NGM_SuddenGradual(In_data, c("generation", "block"))) &&
              (coreValidateBalancing_Ancestral_NGM_SuddenGradual(In_data, c("generation", "plate"))) &&
              (coreValidateBalancing_Ancestral_NGM_SuddenGradual(In_data, c("exper_evol", "NaCl"))) &&
              (coreValidateBalancing_Ancestral_NGM_SuddenGradual(In_data, c("exper_evol", "block"))) &&
              (coreValidateBalancing_Ancestral_NGM_SuddenGradual(In_data, c("exper_evol", "plate"))) &&
              #
              (coreValidateBalancing(In_data, c("generation", "exper_evol", "NaCl"))) &&
              (coreValidateBalancing(In_data, c("generation", "exper_evol", "block"))) &&
              (coreValidateBalancing(In_data, c("generation", "exper_evol", "plate"))) &&
              (coreValidateBalancing_Ancestral_NGM_SuddenGradual(In_data, c("generation", "NaCl", "block"))) &&
              (coreValidateBalancing_Ancestral_NGM_SuddenGradual(In_data, c("generation", "NaCl", "plate"))) &&
              (coreValidateBalancing_Ancestral_NGM_SuddenGradual(In_data, c("generation", "block", "plate"))) &&
              (coreValidateBalancing_Ancestral_NGM_SuddenGradual(In_data, c("exper_evol", "NaCl", "block"))) &&
              (coreValidateBalancing_Ancestral_NGM_SuddenGradual(In_data, c("exper_evol", "NaCl", "plate"))) &&
              (coreValidateBalancing_Ancestral_NGM_SuddenGradual(In_data, c("exper_evol", "block", "plate"))) &&
              (coreValidateBalancing(In_data, c("NaCl", "block", "plate"))) &&
              #
              (coreValidateBalancing(In_data, c("generation", "exper_evol", "NaCl", "block"))) &&
              (coreValidateBalancing(In_data, c("generation", "exper_evol", "NaCl", "plate"))) &&
              (coreValidateBalancing(In_data, c("generation", "exper_evol", "block", "plate"))) &&
              (coreValidateBalancing_Ancestral_NGM_SuddenGradual(In_data, c("generation", "NaCl", "block", "plate"))) &&
              (coreValidateBalancing_Ancestral_NGM_SuddenGradual(In_data, c("exper_evol", "NaCl", "block", "plate"))) &&
              #
              (coreValidateBalancing(In_data, c("generation", "exper_evol", "NaCl", "block", "plate")))),
          "Data appears to be unbalanced")
}

divideIntoSubLists = function(In_lst, max_n_elements) {
  N = length(In_lst)
  if (!max_n_elements) max_n_elements = N
  if (max_n_elements >= N) return (list(In_lst))
  n_parts = ceiling(N / max_n_elements)
  
  Out_lst = array(list(), dim = n_parts)
  for (i in seq(1, n_parts)) {
    if (i == n_parts) {
      rem = N %% max_n_elements
      if (rem > 0) {
        I_i = (i - 1) * max_n_elements + (1:(rem))
      } else {
        I_i = (i - 1) * max_n_elements + (1:max_n_elements)
      }
    } else {
      I_i = (i - 1) * max_n_elements + (1:max_n_elements)
    }
    errorIf(any((I_i < 1) | (I_i > N)), "Out of range")
    Out_lst[[i]] = In_lst[I_i]
  }
  return (Out_lst)
}

#
# "Initialize" field in a list that stores option values
#
setOptFieldIfNonExistent = function(Opts, field_name, field_val) {
  if (is.null(Opts[[field_name]])) Opts[[field_name]] = field_val
  return (Opts)
}

# Append data to a table, by matching the rows in the destination table with those in the source table
appendDataColumns = function(Dest_table, ref_field_name_dest_table, Src_table) {
  errorIf(!is.data.frame(Dest_table) || is.null(colnames(Dest_table)), "The destination table must be a data frame with named columns")
  errorIf(!is.character(ref_field_name_dest_table), "The reference field name in the destination table must be a string")
  errorIf(!is.data.frame(Src_table) || is.null(rownames(Src_table)) || is.null(colnames(Src_table)), 
          "The source table must be a data frame with named rows and columns")
  errorIf(all(colnames(Dest_table) != ref_field_name_dest_table), 
          sprintf("No column named '%s' in the destination table (columns are: {%s})", 
                  ref_field_name_dest_table, paste(colnames(Dest_table), collapse = ";")))
  
  ref_vals_in_src = rownames(Src_table)
  errorIf(!areAllValuesUnique(ref_vals_in_src), "The row names in the source table must be unique")
  
  errorIf(!isEmpty(intersect(colnames(Dest_table), colnames(Src_table))), "Columns with equal name(s)")
  inds_src = rep(NA, length = nrow(Dest_table))
  for (i in seq(1, length(ref_vals_in_src))) {
    I_i = which(Dest_table[, ref_field_name_dest_table] == ref_vals_in_src[[i]])
    inds_src[I_i] = i
  }
  errorIf(any(is.na(inds_src)), "Bug in the code: NA indices")
  errorIf(!areIdenticalContainers(Dest_table[, ref_field_name_dest_table], rownames(Src_table)[inds_src]), "Bug in the code")
  return (cbind(Dest_table, Src_table[inds_src, ]))
}

#
# General wrapper for the mapping: val = map[id]
#
# Modified on January 26th, 2014, to accept non-scalar values (in a recursive way)
#
applyMap = function(name, map_var, raise_err = TRUE, err_msg = "", ret_input_on_err = TRUE, default_val = NA) {
  if (length(name) > 1) {
    wrapper_fun = function (ss) { return (applyMap(ss, map_var, raise_err = raise_err, err_msg = err_msg, 
                                                   ret_input_on_err = ret_input_on_err, default_val = default_val)) }
    if (is.list(name)) { 
      return (lapply(name, wrapper_fun))
    } else if (is.vector(name)) {
      return (sapply(name, wrapper_fun))
    } else {
      stop("Unrecognized format for the input variable")
    }
  } else if (isEmpty(name)) {
    return (NULL)
  }
  
  if (any(name %in% names(map_var))) {
    return (map_var[[name]])
  } else {
    if (raise_err) { msg = paste("Unable to map input value: '", name, "'; ", err_msg, sep = ""); stop(msg) }
    else if (ret_input_on_err) { return (name) }
    else { return (default_val) }
  }
}

# New, optimized version, created on April 17th, 2015
applyMap_fast = function(SS, map_var) {
  # Sanity checks
  errorIf(!isStringVector(SS), "")
  errorIf(!isValidMap(map_var), "")
  errorIf(!isContainedIn(SS, names(map_var)), "")
  return (map_var[SS])
}

# To test:
#   SS = c("thiago", "guzella", "dos", "santos", "testing")
#   map_var = c("thiago" = "Thiago", "guzella" = "Guzella")
applyMap_fast_v2 = function(SS, map_var, raise_err = TRUE, err_msg = "", ret_input_on_err = TRUE, default_val = NA) {
  # Sanity checks
  errorIf(!isStringVector(SS), "")
  errorIf(!isValidMap(map_var), "")
  
  if (isEmpty(setdiff(SS, names(map_var)))) {
    return (map_var[SS])
  } else {
    errorIf(raise_err, err_msg)
    # Need to keep track of the ones that are not in the map
    bool_in_map = SS %in% names(map_var)
    out_SS = SS # will leave the input if it is not in the map
    if (any(bool_in_map)) out_SS[bool_in_map] = map_var[SS[bool_in_map]]
    if (!ret_input_on_err) out_SS[!bool_in_map] = default_val
    return (out_SS)
  }
}

mkMapFromTable = function(table_map, col_from, col_to, flag_col_to_as_character = FALSE) {
  errorIf(!is.data.frame(table_map), "The table map must be a data frame")
  errorIf(!isSingleStr(col_from) || !isSingleStr(col_to), "The names of the target columns must be single strings")
  errorIf(!isContainedIn(c(col_from, col_to), colnames(table_map)), sprintf("No columns named '%s' and/or '%s' in the input table", col_from, col_to))
  
  # Retrieve unique rows
  # Note that an error will be raised when generating the map if there are inconsistencies, so that it is not necessary to perform checks here
  table_map = unique(table_map[, c(col_from, col_to)])
  
  # Generate a map out of the data in the table
  if (flag_col_to_as_character) {
    data_map = mkMap(as.character(table_map[, col_to]), as.character(table_map[, col_from]))
  } else {
    data_map = mkMap(table_map[, col_to], as.character(table_map[, col_from]))
  }
  return (data_map)
}

# When the 'col_from' corresponds to the row names
mkMapFromTableRowNames = function(table_map, col_to, flag_col_to_as_character = FALSE) {
  errorIf(!is.data.frame(table_map), "The table map must be a data frame")
  errorIf(!isSingleStr(col_to), "The names of the target columns must be single strings")
  errorIf(!(col_to %in% colnames(table_map)), sprintf("No column named '%s' in the input table", col_to))
  
  # Generate a map out of the data in the table
  if (flag_col_to_as_character) {
    data_map = mkMap(as.character(table_map[, col_to]), rownames(table_map))
  } else {
    data_map = mkMap(table_map[, col_to], rownames(table_map))
  }
  return (data_map)
}

# Similar to applyMap, but when the map is a table
# applyMap = function(name, map_var, raise_err = TRUE, err_msg = "", ret_input_on_err = TRUE, default_val = NA) {
applyTableMap = function(name, table_map, col_from, col_to, flag_col_to_as_character = FALSE, 
                         raise_err = TRUE, err_msg = "", ret_input_on_err = TRUE, default_val = NA) {
  # Sanity checks
  errorIf(!isStringVector(name), "The name(s) must be stored in a string vector")
  errorIf(!is.data.frame(table_map), "The table map must be a data frame")
  errorIf(!isSingleStr(col_from) || !isSingleStr(col_to), "The names of the target columns must be single strings")
  errorIf(!isContainedIn(c(col_from, col_to), colnames(table_map)), sprintf("No columns named '%s' and/or '%s' in the input table", col_from, col_to))
  
  # Setup a map using the two columns supplied
  data_map = mkMapFromTable(table_map, col_from, col_to, flag_col_to_as_character = flag_col_to_as_character)
  
  return (applyMap(name, data_map, raise_err = raise_err, err_msg = err_msg, ret_input_on_err = ret_input_on_err, default_val = default_val))
}

applyTableMap_fast = function(name, table_map, col_from, col_to, flag_col_to_as_character = FALSE) {
  # Sanity checks
  errorIf(!isStringVector(name), "The name(s) must be stored in a string vector")
  errorIf(!is.data.frame(table_map), "The table map must be a data frame")
  errorIf(!isSingleStr(col_from) || !isSingleStr(col_to), "The names of the target columns must be single strings")
  errorIf(!isContainedIn(c(col_from, col_to), colnames(table_map)), sprintf("No columns named '%s' and/or '%s' in the input table", col_from, col_to))
  
  # Setup a map using the two columns supplied
  data_map = mkMapFromTable(table_map, col_from, col_to, flag_col_to_as_character = flag_col_to_as_character)
  
  return (applyMap_fast(name, data_map))
}

#
# Fixed a minor bug in this function on July 5th, 2014 (which would only lead to problems if 'default_val != NA')
# Modified on August 21st, 2014, to allow for empty vectors in the input list 'Lst'
# Modified on September 15th, 2014, to introduce the check for '!isEmpty(names_i) && !is.null(names_i)'
#
lst2ColsInMat = function(Lst, feature_names, default_val = NA) {
  # Sanity checks
  errorIf(!is.list(Lst) || any(sapply(Lst, function(X) { return (!is.vector(X) || (!isEmpty(X) && is.null(names(X)))) })), 
          "The input data must be stored in a list of named vectors")

  N_features = length(feature_names)
  NN = max(unique(sapply(Lst, length))) # use 'max' to allow for NA values
  errorIf(any(NN > N_features), "Incongruency in # features")
  N_cols = length(Lst)
  col_names = names(Lst)
  A = matrix(data = default_val, nrow = N_features, ncol = N_cols, dimnames = list(feature_names, col_names))
  for (i in seq(1, N_cols)) {
    a_i = Lst[[i]]
    names_i = names(a_i)[!is.na(a_i)]
    if (!isEmpty(names_i) && !is.null(names_i)) {
      errorIf(any(!nchar(names_i)), sprintf("Some of the names are empty: '%s'", paste(names_i, collapse = ", ")))
      errorIf(!isContainedIn(names_i, feature_names), sprintf("Unexpected feature(s) [names = {'%s'}]", paste(names_i, collapse = ", ")))
      A[names_i, i] = a_i
    }
  }
  return (A)
}

#
# Fixed a minor bug in the sanity checks on August 20th, 2014
#
matCols2Lst = function(Mat, flag_remove_NA = FALSE) {
  # Sanity checks
  errorIf(!is.matrix(Mat) || is.null(colnames(Mat)), "The input must be a matrix")
  
  lst = list()
  for (i in seq(1, ncol(Mat))) {
    if (flag_remove_NA) {
      a_i = Mat[, i]
      inds_rows = which(!is.na(a_i))
      lst[[i]] = setNames(a_i[inds_rows], names(inds_rows))
    } else {
      lst[[i]] = setNames(Mat[, i], rownames(Mat))
    }
  }
  return (setNames(lst, colnames(Mat)))
}
reOrderCols = function(A, new_order) {
  # Sanity checks
  errorIf((!is.data.frame(A) && !is.matrix(A)) || (is.null(colnames(A))), "")
  
  cols_1 = intersect(new_order, colnames(A))
  cols_2 = setdiff(colnames(A), new_order)
  return (A[, c(cols_1, cols_2)])
}

reOrderElements = function(X, new_order) {
  # Sanity checks
  errorIf((!is.vector(X) && !is.list(X)) || (is.null(names(X))), "")
  names_1 = intersect(new_order, names(X))
  names_2 = setdiff(names(X), new_order)
  return (X[c(names_1, names_2)])
}

canSaveFile = function(file_path, bool_allow_overwrite = NULL) {
  if (file.exists(file_path)) {
    cat(sprintf(sprintf("File '%s' already exists.\n", file_path)))
    if (is.null(bool_allow_overwrite)) {
      repeat {
        ss = tolower(readline(sprintf("File '%s' (in folder '%s') already exists. Overwrite (y/n)? ", basename(file_path), dirname(file_path))))
        if (ss == "y") {
          return (TRUE)
          break
        } else if (ss == "n") {
          return (FALSE)
          break
        } else {
          cat(sprintf("Unexpected reply."))
        }
      }
    } else {
      errorIf(!is.logical(bool_allow_overwrite), "The flag to allow overwriting a file must be a logical variable")
      return (bool_allow_overwrite)
    }
  } else {
    return (TRUE)
  }
}

#
# Format string according to the standard, in which the first letter is upper-case, and the rest are lower-case
#
getStdStr = function(ss) {
  nn = nchar(ss)
  return (paste(toupper(), tolower(), sep = ""))  
}

isFullyNamedDataFrame = function(D) {
  return (is.data.frame(D) && (!is.null(rownames(D))) && (!is.null(colnames(D))))
}

isDataFrameWithNamedCols = function(D) {
  return (is.data.frame(D) && (!is.null(colnames(D))))
}

isStringMatrix = function(D) {
  return (is.matrix(D) && is.character(D))
}

isFullyNamedMatrix = function(D) {
  return (is.matrix(D) && (!is.null(rownames(D))) && (!is.null(colnames(D))))
}

isFullyNamedStringMatrix = function(D) {
  return (isStringMatrix(D) && (!is.null(rownames(D))) && (!is.null(colnames(D))))
}

isFullyNamedNumericMatrix = function(D) {
  return (isFullyNamedMatrix(D) && is.numeric(D))
}

isFullyNamed3DArray = function(X) {
  return (is.array(X) && (length(dim(X)) == 3) && all(!sapply(dimnames(X), is.null)))
}

isFullyNamedNumeric3DArray = function(X) {
  return (isFullyNamed3DArray(X) && is.numeric(X))
}

isSingleStr = function(ss) {
  return (isScalar(ss) && is.character(ss))  
}

isStringVector = function(SS) { return (is.vector(SS) && is.character(SS)) }

isNamedStringVector = function(SS) { return (isStringVector(SS) && isNamedVector(SS)) }

# Sort based on the numeric value
sortGOIdStringsBasedOnNumValue = function(GO_ids) {
  errorIf(any(substr(GO_ids, 1, nchar("GO:")) != "GO:"), "Invalid GO id(s)")
  go_num_val_str = sapply(GO_ids, function (ss) { return (substr(ss, nchar("GO:") + 1, nchar(ss))) }) # helper, for sorting
  num_go_ids = as.numeric(go_num_val_str)
  errorIf(any(is.na(num_go_ids)), "Unable to convert the parsed GO ids to numeric values")
  return (GO_ids[sort(num_go_ids, index.return = TRUE)$ix])
}

invertMap = function(X) {
  errorIf(is.null(names(X)), "The input has no names")
  new_names = as.character(X)
  errorIf(!areAllValuesUnique(new_names), "The values of the input are not unique")
  return (setNames(names(X), new_names))
}

checkStringsStartingWith = function(Names, token) { return (substr(Names, 1, nchar(token)) == token) }
checkStringsEndWith = function(Names, token) { return (sapply(Names, function(ss) substr(ss, nchar(ss) - nchar(token) + 1, nchar(ss)) == token)) }

#
# A wrapper for the typically used line "X = X[sort(some_fun(X), index.return = TRUE)$ix]", but including a sanity check to make 
# sure that no elements are lost
#
# Raises an error if the input is empty
#
# Usage:
#   X = reorderBasedOnIndicesFromSorting(X, sort(some_fun(X), index.return = TRUE)$ix)
#
reorderBasedOnIndicesFromSorting = function(X, sorted_inds) {
  errorIf(!is.vector(X) && !is.list(X), "The input must be a vector or list")
  errorIf(isEmpty(X), "The input is empty")
  errorIf(any(is.na(sorted_inds)), "NA values in the sorted indices")
  errorIf(!areIdenticalContainers(sort(sorted_inds), 1:length(X)), "Indices are NOT complete")
  return (X[sorted_inds])
}

#
# Wrapper for lapply / mclapply
#
lapplyWrapper = function(lst, func, n_cores = 0) { 
  # Sanity checks
  errorIf(!is.numeric(n_cores) || !isScalar(n_cores), "The number of cores must be a numeric scalar")
  
  if (n_cores < 2) {
    return (lapply(lst, func))
  } else {
    errorIf(!require(parallel), "Unable to initialize 'parallel' package")
    return (mclapply(lst, func, mc.cores = n_cores))
  }
}

#
# Wrapper for lapply to set the names. Equivalent to: X = setNames(lapply(SS, fun), SS)
#
lapplyNamesWrapper = function(SS, fun) {
  errorIf(!is.character(SS), "The input is not a set of names")
  errorIf(!is.null(names(SS)), "The input is named already")
  errorIf(!areAllValuesUnique(SS), "The names are not unique")
  return (setNames(lapply(SS, fun), SS))
}

sapplyNamesWrapper = function(SS, fun) {
  errorIf(!is.character(SS), "The input is not a set of names")
  errorIf(!is.null(names(SS)), "The input is named already")
  errorIf(!areAllValuesUnique(SS), "The names are not unique")
  return (setNames(sapply(SS, fun), SS))
}

#
# Determine the overlap between pairs of lists / vectors
# Possible operating modes
#   "count": determine the number of overlapping items. Returns a matrix of integer values
#   "items": determine the actual overlapping items. Returns a matrix of lists
#   "both": determine both. Returns a list with entries 'count' and 'items'
#
determineOverlap = function(lst, op_mode = "count") {
  # Sanity checks
  errorIf(!is.list(lst) || is.null(names(lst)), "The input must be stored in a named list")
  errorIf(!areAllValuesUnique(names(lst)), "The names of the entries of the list must be unique")
  errorIf(!all(sapply(lst, function(X) { return (is.vector(X) || is.list(X)) })), "The input data must consist of vectors or lists")
  errorIf(!(op_mode %in% c("count", "items", "both")), "Invalid operating mode")
  
  labels = names(lst)
  N = length(lst) 
  errorIf(!N, "Data is empty")
  
  Items_matrix = matrix(data = list(), nrow = length(labels), ncol = length(labels), dimnames = list(labels, labels))
  Count_matrix = matrix(data = NA, nrow = length(labels), ncol = length(labels), dimnames = list(labels, labels))
  
  for (i in seq(1, N)) {
    X_1 = lst[[i]]
    for (j in seq(i, N)) {
      X_2 = lst[[j]]
      
      X_12 = intersect(X_1, X_2)
      n_12 = length(X_12)
      
      if (op_mode == "count") {
        Count_matrix[[i, j]] = n_12
        if (i != j)
          Count_matrix[[j, i]] = n_12
      } else if (op_mode == "items") {
        Items_matrix[[i, j]] = X_12
        if (i != j)
          Items_matrix[[j, i]] = X_12
      } else if (op_mode == "both") {
        Count_matrix[[i, j]] = n_12
        Items_matrix[[i, j]] = X_12
        if (i != j) {
          Count_matrix[[j, i]] = n_12
          Items_matrix[[j, i]] = X_12
        }
      } else {
        myError("Invalid operating mode")
      }
    }
  }
  
  if (op_mode == "count") {
    return (Count_matrix)
  } else if (op_mode == "items") {
    return (Items_matrix)
  } else if (op_mode == "both") {
    return (list(count = Count_matrix, items = Items_matrix))
  } else {
    myError("Invalid operating mode")
  }
}

mkOverlapMatrixStr = function(Overlap_matrix, results_test_overlap = NULL) {
  Overlap_matrix_str = Overlap_matrix
  class(Overlap_matrix_str) = "character"
  if (!is.null(results_test_overlap)) { # add significance information
    inds_signif = which(results_test_overlap$Table[, "signif_overlap"])
    for (ind in inds_signif) {
      label_1 = results_test_overlap$Table[[ind, "entry_1"]]
      label_2 = results_test_overlap$Table[[ind, "entry_2"]]
      Overlap_matrix_str[[label_1, label_2]] = paste(Overlap_matrix_str[[label_1, label_2]], "*", sep = "")
      Overlap_matrix_str[[label_2, label_1]] = paste(Overlap_matrix_str[[label_2, label_1]], "*", sep = "")
    }
    for (i in seq(1, nrow(Overlap_matrix_str))) {
      Overlap_matrix_str[[i, i]] = paste("[", Overlap_matrix_str[[i, i]], "]", sep = "")
    }
  }
  Overlap_matrix_str[lower.tri(Overlap_matrix_str)] = "" # erase the lower diagonal
  return (Overlap_matrix_str)
}

formatOverlapMatrixForOutput = function(Mat) {
  # OLD CODE, RENDERED DEPRECATED ON SEPT 9TH, 2014, IN LIEU OF FUNCTION 'mkOverlapMatrixStr'
  #   Mat_str = Mat
  #   class(Mat_str) = "character"
  #   Mat_str[lower.tri(Mat_str)] = "" # erase the lower diagonal
  #   return (Mat_str)
  return (mkOverlapMatrixStr(Mat))
}

#
# Take a data frame with multiple columns and "unwrap" it, given a set of "base" and "group" columns
# Example usage (where 'colnames(Overlap_count_matrix) = c("q_th", "count_gradual", "count_sudden", "count_union", "count_overlap")'): 
#
# Full_overlap_count_data = unwrapDataFrame(Overlap_count_matrix, c("q_th"), 
#                                           c("Gradual" = "count_gradual", "Sudden" = "count_sudden", 
#                                             "Union[Gradual, Sudden]" = "count_union", "Overlap[Gradual, Sudden]" = "count_overlap"), "count", "tp")
#
# will return a data frame with columns: "q_th", "count" and "tp", where "tp" will take the values: 
#   c("Gradual", "Sudden", "Union[Gradual, Sudden]", "Overlap[Gradual, Sudden]")
#
# This is particularly useful for "reshaping" data for use with ggplot2
#
# Function created on July 16th, 2014
#
unwrapDataFrame = function(D, base_col_names, group_col_names_map, data_col_name, group_name_col_name = "group") {
  # Sanity checks
  errorIf((!is.data.frame(D) && !is.matrix(D)) || is.null(colnames(D)), "The input data must a data frame with named columns")
  errorIf((!isEmpty(base_col_names) && !isStringVector(base_col_names)) || !isContainedIn(base_col_names, colnames(D)), "Invalid base column name(s)")
  errorIf(!isStringVector(group_col_names_map) || !isContainedIn(group_col_names_map, colnames(D)) || 
            !areAllValuesUnique(group_col_names_map) || !areAllValuesUnique(names(group_col_names_map)), "Invalid group column name(s)")
  errorIf(!isEmpty(intersect(base_col_names, group_col_names_map)), "Overlap between base and group column name(s)")
  errorIf(!isSingleStr(group_name_col_name) || (group_name_col_name %in% base_col_names), "Invalid group name column name")
  lst = list()
  n_groups = length(group_col_names_map)
  new_col_names = c(base_col_names, data_col_name, group_name_col_name)
  for (i in seq(1, n_groups)) {
    col_name = group_col_names_map[[i]]
    entry_name = names(group_col_names_map)[[i]]    
    
    if (!isEmpty(base_col_names)) {
      lst[[i]] = setColNames(cbind(D[, c(base_col_names, col_name)], entry_name), new_col_names)
    } else {
      lst[[i]] = setColNames(data.frame(D[, col_name], entry_name, stringsAsFactors = FALSE), c(data_col_name, group_name_col_name))
    }
  }
  return (do.call(rbind, lst))
}

#
# Manually tuned axis configs for the log count plots
# Works well for the two "typical" values of the number of features: 443 and 9505
#
getLogFeatureCountAxisNfo = function(N_features, flag_reduced = FALSE) {
  # Sanity checks
  errorIf(!isScalar(N_features) || !is.numeric(N_features) || (N_features <= 0), "Invalid number of features")
  
  y_min = 0.9
  y_max = round(N_features * 1.2 / 100) * 100
  if (!flag_reduced) {
    y_breaks = c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000)
  } else {
    y_breaks = c(1, 10, 100, 1000, 10000)
  }
  y_breaks = y_breaks[(y_min <= y_breaks) & (y_breaks <= y_max)]
  return (list(limits = c(y_min, y_max), breaks = y_breaks))
}

#
# Manually tuned configs for a frequency axis
#
getFreqAxisNfo = function() {
  return (list(limits = c(0, 1), breaks = c(0, 0.05, 0.15, 0.25, 0.5, 0.75, 1), labels = c("0", "0.05", "0.15", "0.25", "0.5", "0.75", "1")))
}

#
# Common function to simplify a data frame
# This can be used in the last step of converting a list of lists into a data frame
# Provides a much cleaner format for the typical command: 'do.call(rbind, lapply(..., function(...) { return (list) }))', 
#   in which "standard" columns are unwrapped, while the remaining are treated as lists
#
# Created on July 30th, 2014
# Modified on August 11th, 2014, to include 'integer' as a supported type of variable
# Modified on August 11th, 2014, to include 'double' as a supported type of variable
# Modified on September 2nd, 2014, to fix a bug in the check to see if all values are scalar (had to add the 'all' function)
# Modified on September 11th, 2014, to include 'logical' as a supported type of variable
# Modified on October 13th, 2014, to allow for processing NA values
#
simplifyDataFrame = function(D, supported_var_types = c("character", "numeric", "integer", "double", "logical")) {
  # Sanity checks
  errorIf(!is.data.frame(D), "The input must be stored in a data frame")
  
  if (is.null(D) || (!nrow(D)))
    return (NULL)
  
  # Process each column to see if it can be simplified
  for (i in seq(1, ncol(D))) {
    if ((class(D[[i]]) == "list") && (all(sapply(D[[i]], length) == 1))) { # Make sure the values in the list are all scalars
      uniq_var_types = unique(sapply(D[[i]], typeof))
      if (isScalar(uniq_var_types) && (uniq_var_types %in% supported_var_types)) { # "accepted" types of variable
        class(D[[i]]) = uniq_var_types
      } else if (length(uniq_var_types) == 2) {
        data_i = D[[i]]
        data_i = data_i[!sapply(data_i, is.na)]
        uniq_var_types = unique(sapply(data_i, typeof))
        if (all(sapply(data_i, length) == 1)) {
          if (isScalar(uniq_var_types) && (uniq_var_types %in% supported_var_types)) {
            class(D[[i]]) = uniq_var_types
          } else if (setequal(uniq_var_types, c("numeric", "integer")) || setequal(uniq_var_types, c("double", "integer"))) { # clause added on April 10th, 2015
            class(D[[i]]) = "numeric"
          }
        }
      }
    }
  }
  return (D)
}

#
# Wrapper for 'as.data.frame(do.call(rbind, lapply(...)))'
#
myDFRbindWrapper = function(lst, df_stringsAsFactors = FALSE) {
  # Sanity checks
  errorIf(!is.list(lst), "The input must be stored in a list")
  return (simplifyDataFrame(as.data.frame(do.call(rbind, lst), stringsAsFactors = df_stringsAsFactors)))
}

myHistFromDiscreteVals = function(X, mn = NULL, mx = NULL, increment_val = 1) {
  un_X = unique(X)
  if (is.null(mx)) 
    mx = max(un_X, na.rm = TRUE)
  if (is.null(mn))
    mn = min(un_X, na.rm = TRUE)
  errorIf(mn > mx, "Invalid minimum/maximum values")
  
  all_x = seq(mn, mx, by = increment_val)
  count = sapply(all_x, function(x) { return (sum(X == x)) })
  freq = count / sum(count)
  
  return (data.frame(x = all_x, count = count, freq = freq))
}

#
# Function that returns the colors used by ggplot by default
# Obtained from: http://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
#
ggplotColours = function(n = 6, h = c(0, 360) + 15) {
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

# Helper function
recycleVals = function(x, n_output, n_x_to_use = NULL) {
    # Sanity checks
    errorIf(!is.vector(x), "")
    errorIf(!isPositiveNumericScalar(n_output), "")
    
    NN = length(x)
    if (!is.null(n_x_to_use)) {
        errorIf(!isPositiveNumericScalar(n_x_to_use), "")
        errorIf(n_x_to_use > NN, "Invalid 'n_x_to_use'")
        NN = n_x_to_use
    }
    return (x[((0:(n_output - 1)) %% NN) + 1])
}

#
# Function that "shapes" to be used by ggplot-based code
# Obtained from: http://sape.inf.usi.ch/quick-reference/ggplot2/shape
#
getShapesForPlotting = function(n, base_vals = c(16, 17, 15, 3, 7, 8, 11, 18, 13), n_base_vals_to_use = NULL) {
    errorIf(!n, "The number of shapes must be positive")
    return (recycleVals(base_vals, n, n_x_to_use = n_base_vals_to_use))
}

getLinetypesForPlotting = function(n, base_vals = c("solid", "dashed", "dotted", "dotdash"), n_base_vals_to_use = NULL) {
    errorIf(!isPositiveNumericScalar(n), "The number of linetypes must be positive")
    return (recycleVals(base_vals, n, n_x_to_use = n_base_vals_to_use))
}

#
# Full comparison of two vectors, accounting for NA values
#
compareVectors = function(A, B, compare_names = TRUE) {
  errorIf(!is.vector(A) || !is.vector(B), "The inputs must be vectors")
  if (length(A) != length(B))
    return (FALSE)
  if (compare_names) {
    errorIf(is.null(names(A)) || is.null(names(B)), "No names in the vectors")
    if (!compareVectors(names(A), names(B), compare_names = FALSE)) # do a comparison of the names
      return (FALSE)
  }
  bool_NA_A = is.na(A)
  bool_NA_B = is.na(B)
  if (any(bool_NA_A != bool_NA_B))
    return (FALSE)
  results = areIdenticalContainers(A[!bool_NA_A], B[!bool_NA_B])
  errorIf(is.na(results), "Bug in the code for function 'compareVectors': NA value returned")
  return (results)
}

myWriteTableWrapper = function(Data, file_path, 
                               sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE) {
  write.table(Data, file = file_path, sep = sep, quote = quote, row.names = row.names, col.names = col.names)
}

myReadTableWrapper = function(file_path, 
                              header = TRUE, sep = "\t") {
  return (read.table(file_path, header = header, sep = sep))
}

myReadDelimWrapper = function(file_path, 
                              header = TRUE, sep = "\t") {
  return (read.delim(file_path, header = header, sep = sep))
}

# Helper funtion to write a set of strings to a file, and then a table
writeTableAfterText = function(file_path, text_SS, table_data, 
                               sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE) {
  # Sanity checks
  errorIf(!is.list(text_SS) || any(!sapply(text_SS, is.character)), "The text must be stored in a list of strings")
    
  file_conn = file(file_path, open = "wt")
  lapply(text_SS, function(ss) { return (writeLines(ss, file_conn)) })
  write.table(table_data, file = file_conn, sep = sep, quote = quote, row.names = row.names, col.names = col.names)
  close(file_conn)
}

# Helper function for appending a string to the names
appendTokenToNames = function(lst, token_str, flag_leading = TRUE) {
  errorIf(!is.vector(lst) || is.null(names(lst)), "The input strings must be stored in a named vector")
  orig_names = names(lst)
  if (flag_leading) {
    new_names = paste(token_str, orig_names, sep = "")
  } else {
    new_names = paste(orig_names, token_str, sep = "")
  }
  return (setNames(lst, new_names))
}

# Helper function for appending a string to the column names
# In principle, 'token_str' may be a vector of strings, with matching size
appendTokenToColNames = function(A, token_str, flag_leading = TRUE) {
  errorIf((!is.matrix(A) && (!is.data.frame(A))) || is.null(colnames(A)), "The input must be a matrix or data frame")
  errorIf(!isStringVector(token_str) || all(length(token_str) != c(1, ncol(A))), "Size mismatch between matrix/data.frame and token string")
  orig_names = colnames(A)
  if (flag_leading) {
    new_names = paste(token_str, orig_names, sep = "")
  } else {
    new_names = paste(orig_names, token_str, sep = "")
  }
  return (setColNames(A, new_names))
}

# Helper function for appending a string to the column names
# In principle, 'token_str' may be a vector of strings, with matching size
appendTokenToRowNames = function(A, token_str, flag_leading = TRUE) {
  errorIf((!is.matrix(A) && (!is.data.frame(A))) || is.null(colnames(A)), "The input must be a matrix or data frame")
  errorIf(!isStringVector(token_str) || all(length(token_str) != c(1, ncol(A))), "Size mismatch between matrix/data.frame and token string")
  orig_names = rownames(A)
  if (flag_leading) {
    new_names = paste(token_str, orig_names, sep = "")
  } else {
    new_names = paste(orig_names, token_str, sep = "")
  }
  return (setRowNames(A, new_names))
}

#
# Created a reference classification table, which can be used to classify features based on binary patterns (e.g., statistical significance 
# of evolution)
#
# Usage: mkRefClassifTable(list(W_tilde_L = c(FALSE, TRUE), W_tilde_H = c(FALSE, TRUE), Delta_t = c(FALSE, TRUE)))
#
mkRefClassifTable = function(spec_vars, new_order_vars = NULL) {
  # Sanity checks
  errorIf(!is.list(spec_vars) || is.null(names(spec_vars)) || any(sapply(spec_vars, function(A) { return (!is.vector(A)) })), 
          "The variables must be specified in a named list, with elements being vectors")
  
  Ref_classif_table = do.call(expand.grid, spec_vars)
  if (!is.null(new_order_vars)) # re-order the columns
    Ref_classif_table = reOrderCols(Ref_classif_table, new_order_vars)
  return (Ref_classif_table)
}

#
# Classify features based on a reference table
# 
# Example of usage (considering that the input corresponds to a model in which all generations were analyzed altogether):
#   All_var_groups = list(
#     "model_G15_G35_G50" = c("generation" = "generation.G15_G35_G50", "exper_evol" = "exper_evol.G15_G35_G50", "NaCl" = "NaCl.G15_G35_G50"), 
#      "model_G15" = c("exper_evol" = "exper_evol.G15", "NaCl" = "NaCl.G15"), 
#      "model_G35" = c("exper_evol" = "exper_evol.G35", "NaCl" = "NaCl.G35"), 
#      "model_G50" = c("exper_evol" = "exper_evol.G50", "NaCl" = "NaCl.G50"))
#   abbrev_var_names = c("generation" = "gen", "exper_evol" = "reg", "NaCl" = "NaCl") # OPTIONAL
#   results = classifyBinaryFeatureResultsUsingRefTable(Signif_matrix, All_var_groups, abbrev_var_names)
#
# Note that the code can be used to process the results from generations that were considered separately, by properly tuning the 
# mapping that is specified in the variable 'ref_classif_var_labels'
#
classifyBinaryFeatureResultsUsingRefTable = function(Signif_matrix, All_var_groups, abbrev_var_names = NULL) {
  # Sanity checks
  errorIf(!isFullyNamedMatrix(Signif_matrix), "The data on the features must be stored in a fully named matrix")
  errorIf(!is.list(All_var_groups) || is.null(names(All_var_groups)), "The specification of the groups and variables must be stored in a named list")
  errorIf(!isContainedIn(colnames(Signif_matrix), unlist(All_var_groups)), "Some of the variables in the data are not present in the specification of groups")
  
  # Create the (binary-based) reference table
  ref_classif_var_labels = applyMap(colnames(Signif_matrix), invertMap(unlist(noNames(All_var_groups), recursive = FALSE)))
  Ref_classif_table = mkRefClassifTable(lapplyNamesWrapper(noNames(ref_classif_var_labels[!duplicated(ref_classif_var_labels)]), 
                                                           function(SS) { return (c(FALSE, TRUE)) }))
  
  # Include an "id" column in the reference classification table, if the argument 'abbrev_var_names' has been supplied
  if (!is.null(abbrev_var_names))
    Ref_classif_table = cbind(Ref_classif_table, id = apply(Ref_classif_table, 1, function(BB) { 
      SS = colnames(Ref_classif_table)[unlist(BB)]; return (ifelse(isEmpty(SS), "none", paste(applyMap(SS, abbrev_var_names), collapse = "_")))
    }))
  
  return (classifyBinaryFeatureResultsFromRefTable(Signif_matrix, All_var_groups, Ref_classif_table, ref_classif_var_labels))
}

# See function 'classifyBinaryFeatureResultsUsingRefTable' for a detailed description
classifyBinaryFeatureResultsFromRefTable = function(Signif_matrix, All_var_groups, Ref_classif_table, ref_classif_var_labels) {
  # Helper function
  findRowInMat = function(row_data, Mat_data) {
    for (i in seq(1, nrow(Mat_data))) {
      if (all(row_data == Mat_data[i, ]))
        return (i)
    }
    stop("No match for row data")
  }
  
  # Select the var groups to be used
  Var_groups = All_var_groups[sapply(All_var_groups, function(X) { return (isContainedIn(X, colnames(Signif_matrix))) })]
  
  # Assign an id (the row in the reference table) to each feature that is present in the significance matrix
  Ids_of_features = lapply(Var_groups, function(var_data) {
    errorIf(!isContainedIn(var_data, colnames(Signif_matrix)), sprintf("Invalid column names '%s'", paste(col_names, collapse = "', '")))
    # Determine the id to which each feature belongs to
    return (apply(Signif_matrix[, var_data], 1, function(signif_data) { 
      errorIf(!isContainedIn(names(signif_data), names(ref_classif_var_labels)), 
              sprintf("Unable to match labels of the variables ('%s') to columns of the reference classif table ('%s')", 
                      paste(names(signif_data), collapse = "', '"), paste(names(ref_classif_var_labels), collapse = "', '")))
      ref_col_names = applyMap(names(signif_data), ref_classif_var_labels)
      return (findRowInMat(signif_data, Ref_classif_table[, ref_col_names]))
    }))
  })
  
  Count_table = cbind(Ref_classif_table, lapply(Ids_of_features, function(Ids) { 
    return (sapply(1:nrow(Ref_classif_table), function(i_row) { return (sum(Ids == i_row)) }))
  }))
  
  return (list(ref_classif_var_labels = ref_classif_var_labels, Ref_classif_table = Ref_classif_table, 
               Ids_of_features = Ids_of_features, Count_table = Count_table))
}

isValidMap = function(X) {
  return (!is.null(names(X)) && areAllValuesUnique(names(X)))
}

isValidNumMap = function(X) {
  return (isValidMap(X) && isNumericVector(X))
}

isNamedList = function(lst) {
  return (is.list(lst) && !is.null(names(lst)))
}

isNamedVector = function(X) {
  return (is.vector(X) && !is.null(names(X)))
}

isNumericVector = function(X) {
  return (is.vector(X) && is.numeric(X))
}

isNumericMatrix = function(X) {
  return (is.matrix(X) && is.numeric(X))
}

isNumericScalar = function(x) {
  return (isScalar(x) && is.numeric(x))
}

isPositiveNumericScalar = function(x) {
  return (isScalar(x) && is.numeric(x) && (x > 0))
}

isLogicalScalar = function(x) {
  return (isScalar(x) && is.logical(x))
}

isSingleCharacter = function(x) {
  return (isSingleStr(x) && (nchar(x) == 1))
}

isLogicalVector = function(x) { return (is.vector(x) && is.logical(x)) }

isNamedNumericVector = function(X) {
  return (isNumericVector(X) && !is.null(names(X)))
}

#
# My custom code for determining the ECDF out of a data vector
# Note that the 'discard_NA' flag changes the interpretation of the frequencies if there are NA values
#
myECDF = function(X, x_vals = NULL, n_x_vals = NULL, discard_NA = TRUE, x_min = NULL, x_max = NULL) {
  # Sanity checks
  errorIf(!isNumericVector(X), "The input must be a numeric vector")
  
  if (discard_NA) X = X[!is.na(X)]
  
  errorIf(isEmpty(X), "No data")
  
  if (is.null(x_vals)) {
    if (is.null(n_x_vals)) n_x_vals = 100
    if (is.null(x_min)) x_min = min(X, na.rm = TRUE)
    if (is.null(x_max)) x_max = max(X, na.rm = TRUE)
    errorIf(x_max < x_min, "")
    
    x_vals = seq(from = x_min, to = x_max, length.out = n_x_vals)
  }
  x_vals = unique(sort(x_vals))
  NN = length(X)
  
  X_count = sapply(x_vals, function(x_val) { return (length(which(X <= x_val))) }) # use 'which' to discard NA
  return (data.frame(val = x_vals, count = X_count, count_rev = NN - X_count, freq = X_count / length(X)))
}

myECDFFromUniqVals = function(X, discard_NA = TRUE, addit_vals = c()) {
  return (myECDF(X, x_vals = sort(c(addit_vals, unique(X))), discard_NA = TRUE))
}

# Convenience wrapper when using normalized inputs ([0, 1])
getECDFXValsForNormData = function() { return (c(1e-10, 2e-10, 4e-10, 8e-10, 
                                                 1e-9, 2e-9, 4e-9, 8e-9, 
                                                 1e-8, 2e-8, 4e-8, 8e-8, 
                                                 1e-7, 2e-7, 4e-7, 8e-7, 
                                                 1e-6, 2e-6, 4e-6, 8e-6, 
                                                 1e-5, 2e-5, 4e-5, 8e-5, 
                                                 1e-4, 2e-4, 4e-4, 8e-4, 
                                                 1e-3, 2e-3, 4e-3, 8e-3, 
                                                 seq(from = 0.01, to = 1.0, by = 0.01))) }

# Convenience wrapper to apply function 'myECDF' to each column of a matrix, and add a column corresponding to the column name
getECDFFromMatrixColData = function(M, x_vals = NULL, n_x_vals = NULL, discard_NA = TRUE) {
  # Sanity checks
  errorIf(!is.matrix(M) || is.null(colnames(M)), "The input must be a matrix with named columns")
  return (do.call(rbind, lapply(colnames(M), function(col_name) { 
    return (cbind(myECDF(X = M[, col_name], x_vals = x_vals, n_x_vals = n_x_vals, discard_NA = discard_NA), data_name = col_name))
  })))
}

getECDFFromMatrixColData_uniqVals = function(M, discard_NA = TRUE) {
  # Sanity checks
  errorIf(!is.matrix(M) || is.null(colnames(M)), "The input must be a matrix with named columns")
  return (do.call(rbind, lapply(colnames(M), function(col_name) { 
    return (cbind(myECDFFromUniqVals(M[, col_name], discard_NA = discard_NA), data_name = col_name))
  })))
}

# Convenience function
# set fun as 'scale_y_log10' or 'scale_y_continuous', for example
setAxisNfo = function(plt, axis_nfo, fun) {
  if (isContainedIn(c("limits", "breaks", "labels"), names(axis_nfo))) {
    plt = plt + fun(limits = axis_nfo$limits, breaks = axis_nfo$breaks, labels = axis_nfo$labels)
  } else if (isContainedIn(c("limits", "breaks"), names(axis_nfo))) {
    plt = plt + fun(limits = axis_nfo$limits, breaks = axis_nfo$breaks)
  } else if (!is.null(axis_nfo)) {
    myWarning("Ignoring axis info")
    plt = plt + fun()
  }
  return (plt)
}

#
# A more sophisticated version of function 'writeTableAfterText'
# Receives string(s) and table(s) as input, and writes them to files in the order provided
#
writeTextAndTableToFile = function(file_path, Data, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE) {
  # Sanity checks
  errorIf(!isSingleStr(file_path), "The file path must be stored in a single string")
  errorIf(!is.list(Data), "The data must be stored in a list")
  
  file_conn = file(file_path, open = "wt")
  lapply(Data, function(X) { 
    if (is.matrix(X) || is.data.frame(X)) {
      write.table(X, file = file_conn, sep = sep, quote = quote, row.names = row.names, col.names = col.names)
    } else if (is.character(X) && (isScalar(X))) {
      writeLines(X, file_conn) 
    } else {
      myError("Unrecognized data")
    }
  })  
  close(file_conn)
}

# 
# To force log labels for ggplot2 axis
# Function code from "http://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales"
#
ggplot2_scientific10 = function(x) { parse(text=gsub("e", " %*% 10^", scientific_format()(x))) }

# Helper script to set names in a list or vector
mySetNames = function(A, SS_names, flag_check_uniq_names = TRUE, err_msg_not_uniq = "The supplied names are not unique") {
  errorIf(!is.vector(A) && !is.list(A), "The input must be a list or vector")
  errorIf(!isStringVector(SS_names), "The names must be stored in a string vector")
  errorIf(length(A) != length(SS_names), "Size mismatch")
  errorIf(flag_check_uniq_names && (!areAllValuesUnique(SS_names)), err_msg_not_uniq)
  names(A) = SS_names
  return (A)
}

# Helper function for appending a string
appendToken = function(SS, token_str, flag_leading = TRUE) {
  errorIf(!is.vector(SS) || is.null(names(SS)), "The input strings must be stored in a named vector")
  orig_names = names(SS)
  if (flag_leading) {
    SS = paste(token_str, SS, sep = "")
  } else {
    SS = paste(SS, token_str, sep = "")
  }
  return (setNames(SS, orig_names))
}

mkMap = function(vals, SS_names) {
  X = mySetNames(vals, SS_names)
  errorIf(!isValidMap(X), "An invalid map has been created")
  return (X)
}

mkMapFromFun = function(SS_names, fun) mkMap(fun(length(SS_names)), SS_names)

# Create new map by applying a function over a map
mkNewMap = function(in_map, fun) {
  # Sanity checks
  errorIf(!isValidMap(in_map), "")
  errorIf(!is.function(fun), "")
  
  return (mkMap(sapply(in_map, fun), names(in_map)))
}

mkMapFromSingleVal = function(val, SS_names) {
  # Sanity checks
  errorIf(!isScalar(val), "")
  errorIf(!isStringVector(SS_names), "")
  
  X = mySetNames(rep(val, length.out = length(SS_names)), SS_names)
  errorIf(!isValidMap(X), "An invalid map has been created")
  return (X)
}

# returns the indices
sortStrIndsGivenOrderedSet = function(SS, ordered_set) {
  # Sanity checks
  errorIf(!isStringVector(SS) || !isStringVector(ordered_set), "The inputs must be string vectors")
  errorIf(!areAllValuesUnique(ordered_set), "The ordered set must be unique")
  errorIf(!isContainedIn(SS, ordered_set), "The ordered set is incomplete")
  val_map = setNames(1:length(ordered_set), ordered_set)
  return (sort(applyMap(SS, val_map), index.return = TRUE)$ix)
}

# Returns a list with all possible combinations of a set of elements
getAllPossibleCombs = function(X) {
  errorIf(!areAllValuesUnique(X), "The input values must be unique")
  return (unlist(lapply(1:length(X), function(n) { return (combn(X, n, simplify = FALSE)) }), recursive = FALSE))
}

mkNamedVector = function(SS_names, default_val = NA) {
  # Sanity checks
  errorIf(!isScalar(default_val), "The default value must be a scalar")
  errorIf(!isStringVector(SS_names), "The names must be stored in a string vector")
  return (mySetNames(rep(default_val, length.out = length(SS_names)), SS_names))
}

# General-purpose function for classifying features based on binary patterns, given a list of variables (typically a list of significant SNPs for a 
# set of variables)
classifyFeatsBasedOnBinaryPatterns = function(lst_signif_feats, names_all_feats, label_sep = "&") {
  # Sanity checks
  errorIf(!isNamedList(lst_signif_feats) || !all(sapply(lst_signif_feats, isStringVector)), "The data on the significant features must be stored in a named list of string vectors")
  errorIf(!isStringVector(names_all_feats), "The names of all features must be stored in a string vector")
  errorIf(any(sapply(lst_signif_feats, function(feat_names) { return (!isContainedIn(feat_names, names_all_feats)) })), 
          "Inconsistency between the input data and the names of all features")
    
  var_names = names(lst_signif_feats)
  errorIf(!areAllValuesUnique(var_names), "The names of the variables are expected to be unique")
  
  possible_patterns = do.call(expand.grid, c(lapplyNamesWrapper(var_names, function(var_name) { return (c(0, 1)) }), list(stringsAsFactors = FALSE)))
  labels_for_possible_patterns = apply(possible_patterns, 1, function(pattern) {
    col_names = names(which(pattern == 1))
    return (ifelse(isEmpty(col_names), "None", paste(col_names, collapse = label_sep)))
  })
  errorIf(!areAllValuesUnique(labels_for_possible_patterns), "Failed to generate unique labels for the patterns")
  rownames(possible_patterns) = labels_for_possible_patterns
  
  # Helper function
  findRowInMat = function(row_data, Mat_data) {
    errorIf(length(row_data) != ncol(Mat_data), sprintf("Size mismatch (%d versus %d)", length(row_data), ncol(Mat_data)))
    for (i in seq(1, nrow(Mat_data))) {
      if (all(row_data == Mat_data[i, ]))
        return (i)
    }
    stop("No match for row data")
  }
  
  # String vector storing, for each feature, the label of the category to which it belongs
  Categ_labels_for_features = mkNamedVector(names_all_feats, "")
  for (feat_name in names_all_feats) {
    feat_pattern = sapply(lst_signif_feats, function(feat_names) { return (feat_name %in% feat_names) }) + 0 # convert to 0 and 1
    id_feat = findRowInMat(feat_pattern, possible_patterns)
    Categ_labels_for_features[[feat_name]] = labels_for_possible_patterns[[id_feat]]
  }  
  lst_labels_for_feats = lapplyNamesWrapper(unique(Categ_labels_for_features), function(label) { return (names(which(Categ_labels_for_features == label))) })
  
  return (list(possible_patterns = possible_patterns, Categ_labels_for_features = Categ_labels_for_features, lst_labels_for_feats = lst_labels_for_feats))
}

#
# Wrapper to evaluate an expression and return the WARNINGS that were generated in the process of evaluating it
#
# Usage:
#   results = withWarnings(someFunction(X, Y))
#
# Obtained from 'http://stackoverflow.com/questions/3903157/how-can-i-check-whether-a-function-call-results-in-a-warning/4947528#4947528' on 
# September 30th, 2014
# 
withWarnings <- function(expr) {
  myWarnings <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, list(w))
    invokeRestart("muffleWarning")
  }
  val <- withCallingHandlers(expr, warning = wHandler)
  list(value = val, warnings = myWarnings)
} 

mkMatrix = function(row_names, col_names, data = NA, flag_sparse = FALSE) {
  # Sanity checks
  errorIf(!isStringVector(row_names) || !isStringVector(col_names), "The row and column names must be stored in string vectors")
  if (!flag_sparse) {
    return (matrix(data = data, nrow = length(row_names), ncol = length(col_names), dimnames = list(row_names, col_names)))
  } else {
    library(Matrix)
    return (Matrix(data = data, nrow = length(row_names), ncol = length(col_names), dimnames = list(row_names, col_names), sparse = TRUE))
  }
}

mkSparseMatrix = function(row_names, col_names, data = NA) {
  return (mkMatrix(row_names, col_names, data = data, flag_sparse = TRUE))
}

isListOfStringVectors = function(lst) {
  return (is.list(lst) && (all(sapply(lst, function(SS) { return (isStringVector(SS)) }))))
}

# Custom wrapper with sanity check
setAsFactorsGivenLevels = function(SS, lvs) {
  errorIf(!isStringVector(SS), "The input values must be stored in a string vector")
  errorIf(!isStringVector(lvs), "The levels must be stored in a string vector")
  errorIf(!areAllValuesUnique(lvs), "The levels must be unique strings")
  errorIf(!isContainedIn(SS, lvs), "Invalid levels/values")
  return (factor(SS, levels = lvs))
}

# Simple funciton for convenience
discardNullEntriesFromList = function(lst) {
  # Sanity checks
  errorIf(!is.list(lst), "The input must be a list")
  return (lst[!sapply(lst, is.null)])
}

isSorted = function(X) {
  # Sanity checks
  errorIf(!is.vector(X) && !is.list(X), "The input must be a vector or list")
  X = X[!is.na(X)] # exclude NA values
  N = length(X)
  if (N > 1) {
    return (all(1:N == sort(X, index.return = TRUE)$ix))
  } else {
    return (TRUE)
  }
}

#
# Calculate the standard error of the mean
#
sem = function(X) {
  # Sanity checks
  errorIf(!isNumericVector(X), "The input must be a numeric vector")
  X = X[!is.na(X)] # discard NA
  n = length(X)
  return (ifelse(n > 1, sd(X) / sqrt(n), NA))
}

# Adaptor for function 'apply' to a matrix of lists
applyToListMat = function(Mat, func) {
  return (apply(Mat, c(1, 2), function(X) {
    errorIf(!is.list(X) || !isScalar(X), "The elements of the table are expected to be singleton lists")
    X = X[[1]]
    return (func(X))
  }))
}

# Wrapper function to invoke function 'table, and return the results as a named vector
tableWrapperAsNamedVector = function(in_data) {
  res = table(in_data)
  # To get rid of attributes...
  return (setNames(as.vector(res), names(res)))
}

#
# Format of the "data selection function" ('sel_fun'): sel_fun(naming_data, val_data)
# Returns a named list
#
# *** DISABLED ON AUGUST 26th, 2015, BECAUSE THIS DOES NOT SEEM TO MAKE SENSE
if (0) {
  clusterValuesBasedOnNames = function(naming_data, val_data, sel_fun = NULL) {
    # Sanity checks
    errorIf(!isStringVector(naming_data), "The naming data must be stored in a string vector")
    errorIf(!is.vector(val_data) && !is.list(val_data), "The values must be stored in a vector or list")
    N = length(naming_data)
    errorIf(length(val_data) != N, "Size mismatch")
    
    if (!is.null(sel_fun)) {
      bool_use = sel_fun(naming_data, val_data)
      errorIf(!isLogicalVector(bool_use) || (length(bool_use) != N), "Invalid value(s) returned by the selection function")
      inds_use = which(bool_use)
      if (isEmpty(inds_use)) return (list())
      naming_data = naming_data[inds_use]
      val_data = val_data[inds_use]
    }
    
    uniq_names = unique(naming_data)
    return (lapplyNamesWrapper(uniq_names, function(name_str) { return (val_data[uniq_names == name_str]) }))
  }
}

myIfElse = function(b, val_true, val_false) {
  # Sanity check
  errorIf(!isLogicalScalar(b), "The input must be a logical scalar")
  if (b) {
    return (val_true)
  } else {
    return (val_false)
  }
}

#
# Create reverse map out of named list of values
#
# Example:
#   
#
mkRevMapFromNamedListOfVals = function(lst) {
  # Sanity checks
  errorIf(!isNamedList(lst) || !areAllValuesUnique(names(lst)), "The input must be a list with unique names")
  errorIf(!all(sapply(lst, isStringVector)), "The input must be a list of string vectors")
  # Pack the data in a table
  vals = unlist(lst)
  SS = unlist(mapply(function(ss, n) rep(ss, length.out = n), names(lst), sapply(lst, length), SIMPLIFY = FALSE))
  errorIf(!areAllValuesUnique(vals), "The values must be unique")
  return (mkMap(SS, vals))
}

#
# General function to sort the rows of a table, based on a given ordered set of columns to be considered, in which the string values are 
# assigned to numeric values
# The names of the columns to be considered are defined as the names of the argument
#
# For example: 
#   sorting_specs = list("treat" = c("AN" = 1, "FIX" = 2, "MOV" = 3, "NGM" = 4), 
#                        "gen" = c("G0" = 1, "G8" = 2, "G10" = 2, "G35" = 3, "G36" = 3, "G50" = 4), 
#                        "repl_pop" = c("M" = 1, "M1" = 2, "M2" = 3, "M3" = 4, "M5" = 5, "M6" = 6))
#
sortTableRows = function(Table, sorting_specs, flag_return_inds = FALSE) {
  # Sanity checks
  errorIf(!isDataFrameWithNamedCols(Table), "The input table must be a data frame with named columns")
  errorIf(!isNamedList(sorting_specs), "The sorting specs must be stored in a named list")
  errorIf(!isContainedIn(names(sorting_specs), colnames(Table)), "Invalid sorting specs")
  
  # For convenience, reverse the sorting specs
  N_specs = length(sorting_specs)
  errorIf(!N_specs, "No sorting specs")
  sorting_specs = sorting_specs[rev(seq(1, N_specs, by = 1))]
  
  # Build a numeric matrix with the index for each row and column
  N = nrow(Table)
  Ind_matrix = matrix(data = 0, nrow = N, ncol = N_specs)
  for (i in seq(1, N_specs)) {
    col_i = names(sorting_specs)[[i]]
    map_i = sorting_specs[[i]]
    Ind_matrix[, i] = applyMap(as.character(Table[, col_i]), map_i)
  }
  
  # Constant vector containing the set of bases
  K = cumprod(sapply(sorting_specs[1:(N_specs - 1)], length))
  K = as.matrix(c(1, K))
  
  sorted_inds = sort(Ind_matrix %*% K, index.return = TRUE)$ix
  errorIf(!setequal(sorted_inds, 1:N), "Problem in the sorting")
  
  if (flag_return_inds) {
    return (sorted_inds)
  } else {
    return (Table[sorted_inds, ])
  }
}

combineMaps = function(X, Y) {
  # Sanity checks
  errorIf(!isValidMap(X) || !isValidMap(Y), "")
  Z = c(X, Y)
  errorIf(!isValidMap(Z), "")
  return (Z)
}

#
# Sort strings based on their numeric value
# For each string, try to convert it to a numeric value, by removing all non-numeric characters
# This function returns the sorted INDICES; see function 'sortStrsByNumVal' for the actual strings
#
sortStrByNumVal_Inds = function(SS, flag_NA_leading = FALSE) {
  # Sanity checks
  errorIf(!isStringVector(SS), "The input must be a string vector")
  
  # Convert to numeric values
  num_vals = as.numeric(sapply(SS, function(ss) gsub("[^0123456789\\.]", "", ss)))
  bool_NA = is.na(num_vals)
  inds_not_NA = which(!bool_NA)
  inds_NA = which(bool_NA)
  
  inds_not_NA = inds_not_NA[sort(num_vals[inds_not_NA], index.return = TRUE)$ix]
  
  if (flag_NA_leading) {
    res = c(inds_not_NA, inds_NA)
  } else {
    res = c(inds_NA, inds_not_NA)
  }
  errorIf(!setequal(res, 1:length(SS)), "Invalid results obtained")
  
  return (res)
}

sortStrsByNumVal = function(SS, flag_NA_leading = FALSE) {
  return (SS[sortStrByNumVal_Inds(SS, flag_NA_leading = flag_NA_leading)])
}

# Returns a data frame
# This is useful for parsing the names of populations or individuals
parseStringsBySplitting = function(SS, token_names, sep = ".") {
  n_tokens = length(token_names)
  lst_tokens = strsplit(SS, split = sep, fixed = TRUE)
  errorIf(any(sapply(lst_tokens, length) != n_tokens), "Invalid tokens retrieved")
  return (setColNames(myDFRbindWrapper(lst_tokens), token_names))
}

# Returns a named list of string vectors
clusterNamesBasedOnVals = function(X, possible_vals = NULL) {
  # Sanity checks
  errorIf(!isNamedVector(X), "")
  errorIf(!is.null(possible_vals) && !is.vector(possible_vals), "")
  uniq_X = noNames(unique(X))
  
  if (is.null(possible_vals)) {
    possible_vals = uniq_X
  } else {
    errorIf(!isContainedIn(uniq_X, possible_vals), "Invalid possible values supplied")
  }
  
  lst = lapply(possible_vals, function(x) names(X)[X == x])
  return (setNames(lst, as.character(possible_vals)))
}

clusterNamesBasedOnValsFromTable = function(Table, col_from, col_to, flag_col_to_as_character = FALSE) {
  data_map = mkMapFromTable(Table, col_from, col_to, flag_col_to_as_character = flag_col_to_as_character)
  return (clusterNamesBasedOnVals(data_map))
}

# wrapper function for do.call(rbind, lst), with sanity check against duplicated row names
rowConcatLstRowNamedTables = function(lst) {
  # Sanity checks
  errorIf(!is.list(lst), "")
  errorIf(!areAllValuesUnique(unlist(lapply(lst, rownames))), "")
  return (do.call(rbind, noNames(lst)))
}

# Wrapper for retrieving colors
brewerPaletteWrapper = function(nn, palette_name = "Set1", thr = 9) {
  if (nn <= thr) {
    return (brewer.pal(nn, palette_name))
  } else {
    return (colorRampPalette(brewer.pal(thr, palette_name))(nn))
  }
}

brewerPaletteWrapper_fromNames = function(SS, palette_name = "Set1", ...) {
  errorIf(!isStringVector(SS) || !areAllValuesUnique(SS), "")
  return (mkMap(brewerPaletteWrapper(length(SS), palette_name = palette_name, ...), SS))
}

# A generalization of function 'tableWrapperAsNamedVector', but in which the "labels" are supplied by the used
# It is not required that the argument 'possible_labels' correspond to all labels present in the data
getCounts = function(X, possible_labels = NULL) {
  # Sanity checks
  errorIf(!isStringVector(X), "The input must be a string vector")
    
  if (is.null(possible_labels)) {
    possible_labels = unique(X)
  } else {
    errorIf(!isStringVector(possible_labels), "The possible labels must be stored in a string vector")
    errorIf(!areAllValuesUnique(possible_labels), "The possible labels supplied must be unique")
  }
  possible_labels = noNames(possible_labels) # to avoid problems
  
  return (sapplyNamesWrapper(possible_labels, function(label) sum(X == label, na.rm = TRUE)))
}

# Wrapper for function 'getCounts', returning the results in a data frame
getCountsAsDF = function(X, possible_labels = NULL, variable_col_name = "variable", count_col_name = "count") {
  # Sanity checks
  errorIf(!isSingleStr(variable_col_name), "")
  errorIf(!isSingleStr(count_col_name), "")
  errorIf(variable_col_name == count_col_name, "")
  count_data = getCounts(X, possible_labels = possible_labels)
  return (setColNames(data.frame(names(count_data), count_data), c(variable_col_name, count_col_name)))
}

# A generalization of 'getCountsAsDF' for named numeric vectors
getCountDataAsDf = function(X, value_label = "value", count_label = "count") {
  # Sanity checks
  errorIf(!isNamedNumericVector(X), "")
  errorIf(!isSingleStr(value_label), "")
  errorIf(!isSingleStr(count_label), "")
  errorIf(value_label == count_label, "")
  
  uniq_x = sort(noNames(unique(X)))
  return (setColNames(data.frame(uniq_x, sapply(uniq_x, function(x) sum(X == x)), row.names = NULL), c(value_label, count_label)))
}

compareValuesInTables = function(A, B, col_names, eps = 1e-10) {
  errorIf(!isFullyNamedDataFrame(A) || !isFullyNamedDataFrame(B), "")
  errorIf(!isContainedIn(col_names, colnames(A)) || !isContainedIn(col_names, colnames(B)), "")
  errorIf(!setequal(rownames(A), rownames(B)), "")
  B = B[rownames(A), , drop = FALSE]
  comp_vals = sapply(col_names, function(ss) max(abs(A[, ss] - B[, ss])))
  return (all(!is.na(comp_vals)) && all(comp_vals < eps))
}

# Can be used with either vectors or lists
# Format of the function: fun(name, value)
applyFunToNamedContainer = function(X, fun, ...) { # fun(name, value)
  errorIf(!isNamedList(X) && !isNamedVector(X), "")
  return (mapply(fun, names(X), X, ...))
}

# Wrapper for function 'applyFunToNamedContainer'
# Format of the function: fun(name, value)
applyFunToNamedContainer_lst = function(X, fun) { # function(name, value)
  errorIf(!isNamedList(X) && !isNamedVector(X), "")
  return (mapply(fun, names(X), X, SIMPLIFY = FALSE))
}

# For c("name_1" = 0, "name_2" = 1), produce: "name_1 = 0, name_2 = 1"
namedVec2Str = function(X, fmt_str = "%d") {
  errorIf(!isNamedVector(X), "")
  fmt_str = paste("%s = ", fmt_str, sep = "")
  return (applyFunToNamedContainer(X, function(name, val) sprintf(fmt_str, name, val)))
}

# For c("name_1" = 0, "name_2" = 1), produces: "name_1 = 0, name_2 = 1"
printNamedVectorVals = function(X, fmt_str = "%d", sep = ", ") {
  errorIf(!isNamedVector(X), "")
  return (paste(namedVec2Str(X, fmt_str = fmt_str), collapse = sep))
}

setupForFig = function(fig_file_path) {
  cat(sprintf("Figure will be stored in file '%s'\n", fig_file_path))
  createFolder(dirname(fig_file_path))
  return (fig_file_path)
}

setupForData = function(data_file_path, quiet = FALSE) {
  if (!quiet) cat(sprintf("Output data will be stored in file '%s'\n", data_file_path))
  createFolder(dirname(data_file_path))
  return (data_file_path)
}

setupFolderForFigs = function(fig_folder_path) {
  cat(sprintf("Figures will be stored in folder '%s'\n", fig_folder_path))
  createFolder(fig_folder_path)
  return (fig_folder_path)    
}

setupFolderForData = function(data_folder_path) {
  cat(sprintf("Data will be stored in folder '%s'\n", data_folder_path))
  createFolder(data_folder_path)
  return (data_folder_path)    
}

formatTokenStrForFilename = function(token, sep_str = "-", flag_leading_sep = TRUE) {
  if (nchar(token)) {
    if (flag_leading_sep) {
      return (paste(sep_str, token, sep = ""))
    } else {
      return (paste(token, sep_str, sep = ""))
    }
  } else {
    return ("")
  }
}

#
# Convert all factor columns in a data frame to strings
# Obtained from "http://stackoverflow.com/questions/2851015/convert-data-frame-columns-from-factors-to-characters" on April 17th, 2015
#
dataFrameFactors2Strs = function(D) {
  # Sanity checks
  errorIf(!is.data.frame(D), "")
  return (as.data.frame(rapply(D, as.character, classes = "factor", how = "replace"), stringsAsFactors = FALSE))
}

#
# Apply function to each row of a data frame, preserving the different formats of variables
# The function passed as a parameter will receive a single row, which can be converted to a list using 'as.list
# Example: func = function(row) do_stuff
# Returns a data frame
#
# Obtained from "http://stackoverflow.com/questions/1699046/for-each-row-in-an-r-dataframe" on April 17th, 2015
#
apply2DataFrameRows_df = function(D, func, flag_force_str = TRUE) {
  # Sanity checks
  errorIf(!is.data.frame(D), "")
  if (flag_force_str) D = dataFrameFactors2Strs(D)
  return (noRowNames(do.call(rbind, by(D, 1:nrow(D), func))))
  #return (myDFRbindWrapper(lapply(1:nrow(D), function(i_row) func(D[i_row, , drop = FALSE]))))
}

# Similar to 'apply2DataFrameRows_df', but the function should accept a list
apply2DataFrameRows_lst = function(D, func, flag_force_str = TRUE) {
  # Sanity checks
  errorIf(!is.data.frame(D), "")
  if (flag_force_str) D = dataFrameFactors2Strs(D)
  return (noRowNames(do.call(rbind, by(D, 1:nrow(D), function(row) func(as.list(row))))))
}

paletteAdaptorFun = function(fun, SS, n_before, n_after) {
  # Sanity checks
  errorIf(!isStringVector(SS), "")
  errorIf(!isNumericScalar(n_before) || !isNumericScalar(n_after), "")
  colors = fun(length(SS) + n_before + n_after)
  return (mkMap(colors[n_before + (1:length(SS))], SS))
}

seqFromRange = function(X, n = 100) {
  rng = range(X, na.rm = TRUE)
  return (seq(rng[[1]], rng[[2]], length.out = n))
}

# Given a matrix with counts, obtain a matrix with frequencies, by normalizing the values over the sum for each row
getFreqsFromCountMatrix_byRows = function(A) {
  # Sanity checks
  errorIf(!isNumericMatrix(A) && !(is.data.frame(A) && all(sapply(1:ncol(A), function(i_col) class(A[[i_col]]) %in% c("numeric", "integer", "double")))), "")
  return (A / as.matrix(rowSums(A))[, rep(1, length.out = ncol(A)), drop = FALSE])
}

# Test whether there is overlap between two sets
areDisjoint = function(A, B) {
  return (isEmpty(intersect(A, B)))
}

# Wrapper function for a combination of 'applyMap' and (typically applied to a column of factors in a data frame)
renameAndEnforceOrder = function(A, col_name, label_map) {
  # Sanity checks
  errorIf(!isDataFrameWithNamedCols(A), "")
  errorIf(!isSingleStr(col_name), "")
  errorIf(!(col_name %in% colnames(A)), "")
  
  A[, col_name] = setAsFactorsGivenLevels(applyMap(as.character(A[, col_name]), label_map), unique(label_map))
  return (A)
}

renameValsInDataCol = function(A, col_name, renaming_map) {
  # Sanity checks
  errorIf(!isDataFrameWithNamedCols(A), "")
  errorIf(!isSingleStr(col_name), "")
  errorIf(!(col_name %in% colnames(A)), "")
  errorIf(!isValidMap(renaming_map), "")
    
  A[, col_name] = applyMap_fast(as.character(A[, col_name]), renaming_map)
  return (A)
}

# Wrapper function, simplified version of function 'renameAndEnforceOrder' (no renaming is done)
enforceOrder = function(A, col_name, ordered_labels, flag_subset = FALSE) {
  # Sanity checks
  errorIf(!isDataFrameWithNamedCols(A), "")
  errorIf(!isSingleStr(col_name), "")
  errorIf(!(col_name %in% colnames(A)), "")
  errorIf(!isStringVector(ordered_labels), "")
  
  SS = as.character(A[, col_name])
  ordered_labels = unique(ordered_labels)
  if (flag_subset)
      ordered_labels = ordered_labels[ordered_labels %in% SS]
  
  A[, col_name] = setAsFactorsGivenLevels(SS, ordered_labels)
  return (A)
}

askForPressingEnter = function(msg = "Press <Enter> to continue ") {
  errorIf(!isSingleStr(msg), "")
  cat(msg)
  scan(n = 1, quiet = TRUE)
}

# Wrapper function for 'melt' when dealing with a fully named matrix
meltFullyNamedMatrix = function(A, row_id, variable_name = "variable", value_name = "value") {
  # Sanity checks
  errorIf(!isFullyNamedMatrix(A), "")
  errorIf(!isSingleStr(row_id), "")
  errorIf(row_id %in% colnames(A), "")
  
  library(reshape2)
  return (melt(setColNames(cbind(rownames(A), as.data.frame(A, stringsAsFactors = FALSE), stringsAsFactors = FALSE), c(row_id, colnames(A))), 
               id.vars = c(row_id), measure.vars = setdiff(colnames(A), row_id), variable.name = variable_name, value.name = value_name))
}

# Wrapper function for 'melt' when dealing with a fully named numeric matrix
# Kept here for backwards compatibility
meltFullyNamedNumericMatrix = meltFullyNamedMatrix

# do it by hand (this is faster)
meltFullyNamedNumericMatrix_fast = function(A, row_id, variable_name = "variable", value_name = "value", 
                                            ordered_rows = NULL, ordered_cols = NULL) {
    # Sanity checks
    errorIf(!isFullyNamedMatrix(A), "")
    errorIf(!isSingleStr(row_id), "")
    errorIf(row_id %in% colnames(A), "")
    
    if (!is.null(ordered_rows)) {
        errorIf(!isStringVector(ordered_rows) || !setequal(ordered_rows, rownames(A)), "")
        A = A[ordered_rows, , drop = FALSE]
    }
    if (!is.null(ordered_cols)) {
        errorIf(!isStringVector(ordered_rows) || !setequal(ordered_cols, colnames(A)), "")
        A = A[, ordered_cols, drop = FALSE]
    }
    
    return (setColNames(data.frame(
        rep(rownames(A), times = ncol(A)), rep(colnames(A), each = nrow(A)), as.vector(A)), 
        c(row_id, variable_name, value_name), flag_check_uniq = TRUE)
    )
}

ggplot2SideBySide = function(plot_1, plot_2, title_str, fig_file_path, 
                             rel_horiz_partition = 0.5, 
                             margin_top = 0.2, cex_title_str = 1, total_fig_width = 9, total_fig_height = 5, 
                             rel_lg_label_text_size = 0.6, rel_lg_title_text_size = 0.8) {
    
  # Pack the two plots
  gp_1 = ggplotGrob(plot_1)
  gp_2 = ggplotGrob(plot_2)
  
  # To make sure that the two figures align perfectly in terms of height (which could be a problem because of the difference in the
  # format of the x axis labels)
  # This solution was obtained from:
  #   http://stackoverflow.com/questions/13294952/left-align-two-graph-edges-ggplot/13295880#13295880 
  max_height = grid::unit.pmax(gp_1$heights[2:5], gp_2$heights[2:5])
  gp_1$heights[2:5] <- as.list(max_height)
  gp_2$heights[2:5] <- as.list(max_height)
  
  # Generate the joint figure
  pdf(setupForFig(fig_file_path), height = total_fig_height, width = total_fig_width)
  # Add a newline to the title string, to make space in the top
  grid.arrange(gp_1, gp_2, ncol = 2, 
               main = textGrob(sprintf("%s", title_str), gp = gpar(cex = cex_title_str), vjust = 1), 
               widths = c(rel_horiz_partition * total_fig_width, (1 - rel_horiz_partition) * total_fig_width))
  dev.off()
}

#
# Returns data frame sorted based on the input labels
#
getCumValsGivenOrder = function(vals, labels, ordered_labels, label_col_name = NULL, output_col_name_token = "cum_freq") {
  # Sanity checks
  errorIf(!isNumericVector(vals), "Values must be numeric")
  errorIf(!isStringVector(labels) || !areAllValuesUnique(labels), "Invalid labels")
  errorIf(!isStringVector(ordered_labels), "Invalid ordered labels")
  n = length(vals)
  errorIf(n != length(labels), "Size mismatch")
  errorIf(!isContainedIn(labels, ordered_labels), "Mismatch in the labels")
  
  if (0) {
    cat(sprintf("Vals are: %s\n", paste(vals, collapse = ", ")))
    cat(sprintf("Labels are: %s\n", paste(labels, collapse = ", ")))
    cat(sprintf("Ordered labels are: %s\n", paste(ordered_labels, collapse = ", ")))
    cat(sprintf("\n"))
  }
  
  # Build a map with the input values and the labels
  in_vals = mkMap(vals, labels)
  
  # Sort based on the ordered labels, and determine the cumulative frequencies
  in_vals = in_vals[intersect(ordered_labels, labels)]
  cum_vals = cumsum(in_vals)
  if (n > 1) {
    cum_vals_start = c(0, cum_vals[1:(n - 1)])
  } else {
    cum_vals_start = 0
  }
  cum_vals_start = setNames(cum_vals_start, names(cum_vals))
  cum_vals_end = cum_vals
  
  # Sort back
  cum_vals = cum_vals[labels]
  cum_vals_start = cum_vals_start[labels]
  cum_vals_end = cum_vals_end[labels]
  
  D_out = setColNames(data.frame(X = cum_vals, Y = cum_vals_start, Z = cum_vals_end, row.names = names(cum_vals)), 
                      paste(output_col_name_token, c("", ".start", ".end"), sep = ""), flag_check_uniq = FALSE)
  
  if (!is.null(label_col_name)) {
    return (setColNames(cbind(name = rownames(D_out), D_out), c(label_col_name, colnames(D_out)), flag_check_uniq = FALSE))
  } else {
    return (D_out)
  }
}

# Wrapper function for using 'ddply' when each output data frame grows
ddply_growRowsWrapper = function(D_in, id_col_names, func, target_col_names = NULL) {
  # Sanity checks
  errorIf(!isDataFrameWithNamedCols(D_in), "")
  errorIf(!isStringVector(id_col_names), "")
  errorIf(!is.function(func), "")
  errorIf(!is.null(target_col_names) && !isStringVector(target_col_names), "")
  errorIf(!isContainedIn(id_col_names, colnames(D_in)), "")
  errorIf(!isContainedIn(target_col_names, colnames(D_in)), "")
  if (is.null(target_col_names)) target_col_names = colnames(D_in)
  
  return (ddply(D_in, id_col_names, function(D) {
    D_new = func(D)
    D_subset = D[, target_col_names, drop = FALSE]
    
    bool_bad_cols = sapply(1:ncol(D_subset), function(i_col) length(unique(D_subset[, i_col]))) != 1
    names_bad_cols = colnames(D_subset)[bool_bad_cols]
    vals_bad_cols = sapply(names_bad_cols, function(ss) sprintf("(%s)", paste(unique(D_subset[, ss]), collapse = ", ")))
    
    errorIf(any(bool_bad_cols), sprintf("Bad columns: '%s' => {%s}", paste(names_bad_cols, collapse = ", "), paste(vals_bad_cols, collapse = "; ")))
    return (cbind(D_subset[rep(1, length.out = nrow(D_new)), , drop = FALSE], D_new))
  }))
}

# A wrapper for convenience
mustBeScalar = function(X, msg = "The value is not a scalar") {
  errorIf(length(X) != 1, msg)
  return (X)
}

mustBeNumericScalar = function(X, msg = "The value is not a numeric scalar") {
  errorIf(!isNumericScalar(X), msg)
  return (X)
}

mustBeNonEmpty = function(X, msg = "The value is empty") {
  errorIf(isEmpty(X), msg)
  return (X)
}

mustBeAllUnique = function(X, msg = "The values are not unique") {
  errorIf(!areAllValuesUnique(X), msg)
  return (X)
}

mustBeAllUniqCols = function(D, msg = "The column names are not unique (or have not been set)") {
    errorIf(is.null(colnames(D)) || !areAllValuesUnique(colnames(D)), msg)
    return (D)
}

mustBePositiveNumScalar = function(X, msg = "The value is not a positive numeric scalar") {
  errorIf(!isScalar(X) || !is.finite(X) || (X <= 0), msg)
  return (X)
}

# Patched on Oct 4th, 2016
mustBe = function(X, chk_fun, msg_fun = function(XX) "") {
    chk_val = chk_fun(X)
    if (is.na(chk_val)) chk_val = FALSE
    
    if (is.function(msg_fun)) {
        errorIf(!chk_val, msg_fun(X))
    } else if (isSingleStr(msg_fun)) {
        errorIf(!chk_val, msg_fun)
    } else {
        myError("Unrecognized case")
    }
    return (X)
}

mustBeSingleStr = function(X, msg = "The value is not a single string") mustBe(X, isSingleStr)

mustBeFiniteNumericVector = function(X, msg = "The values must be all finite numeric")
    mustBe(X, function(x) isNumericVector(x) && all(is.finite(x)))

mustBeValidIntegerVector = function(X, msg = "The values must be all valid integers")
    mustBe(X, function(x) is.vector(x) && is.integer(x) && all(!is.na(x)))

mustBeFinitePositiveNumericScalar = function(X, msg = "The value must be finite numeric positive")
    mustBe(X, function(x) isNumericScalar(x) && is.finite(x) && (x > 0))


# Format of the function: fun(name, value)
# applyFunToNamedContainer_lst = function(X, fun) {

# Wrapper function for using 'dlply', keeping the values
# Example: lst_col_specs = list("col_name_1" = as.character, "col_name_2" = as.numeric, ...)
# Set the value as NULL to avoid type cast
dlply_valsWrapper = function(Data, lst_col_specs, id_col_names = NULL, data_entry_name = "Data", ...) {
  # Sanity checks
  errorIf(!isDataFrameWithNamedCols(Data), "")
  errorIf(!isNamedList(lst_col_specs), "")
  errorIf(!isSingleStr(data_entry_name), "")
  errorIf(!is.null(id_col_names) && !isStringVector(id_col_names), "")
  errorIf(!isContainedIn(names(lst_col_specs), colnames(Data)), "")
  errorIf(data_entry_name %in% names(lst_col_specs), "")
  
  if (is.null(id_col_names)) {
    id_col_names = names(lst_col_specs)
  } else {
    errorIf(!isContainedIn(id_col_names, names(lst_col_specs)), "")
  }
  
  dlply(Data, id_col_names, function(D) {
    # Obtain the "header"
    hdr = applyFunToNamedContainer_lst(lst_col_specs, function(col_name, conv_fun) {
      val = unique(D[, col_name])
      if (!is.null(conv_fun)) val = conv_fun(val)
      return (mustBeScalar(val))
    })
    
    return (c(hdr, setNames(list(D[, setdiff(colnames(D), names(lst_col_specs)), drop = FALSE]), data_entry_name)))
  }, ...)
}

# Similar to 'dlply_valsWrapper', but framed in terms of a function:
#   func(D_subset, hdr)
# where hdr is a list containing the header data
dlply_valsWrapper_fun = function(Data, lst_col_specs, func, id_col_names = NULL, flag_include_header = TRUE, ...) { # func(D_subset, hdr)
  # Sanity checks
  errorIf(!isDataFrameWithNamedCols(Data), "")
  errorIf(!isNamedList(lst_col_specs), "")
  errorIf(!is.function(func), "")
  errorIf(!is.null(id_col_names) && !isStringVector(id_col_names), "")
  errorIf(!isContainedIn(names(lst_col_specs), colnames(Data)), "")
  
  if (is.null(id_col_names)) {
    id_col_names = names(lst_col_specs)
  } else {
    errorIf(!isContainedIn(id_col_names, names(lst_col_specs)), "")
  }
  
  dlply(Data, id_col_names, function(D) {
    # Obtain the "header"
    hdr = applyFunToNamedContainer_lst(lst_col_specs, function(col_name, conv_fun) {
      val = unique(D[, col_name])
      if (!is.null(conv_fun)) val = conv_fun(val)
      return (mustBeScalar(val))
    })
    
    lst_D = func(D[, setdiff(colnames(D), names(lst_col_specs)), drop = FALSE], hdr)
    errorIf(!isNamedList(lst_D), "The output of the function must be a named list")
        
    if (flag_include_header) {
      errorIf(!areDisjoint(names(hdr), names(lst_D)), "Naming conflict")
      return (c(hdr, lst_D))
    } else {
      return (lst_D)
    }
  }, ...)
}

# Analogous to 'dlply_valsWrapper', but the core function returns a data frame
#   func(D_subset, hdr)
ddply_valsWrapper_fun = function(Data, lst_col_specs, func, id_col_names = NULL, flag_include_header = TRUE, ...) { # func(D_subset, hdr)
  # Sanity checks
  errorIf(!isDataFrameWithNamedCols(Data), "")
  errorIf(!isNamedList(lst_col_specs), "")
  errorIf(!is.function(func), "")
  errorIf(!is.null(id_col_names) && !isStringVector(id_col_names), "")
  errorIf(!isContainedIn(names(lst_col_specs), colnames(Data)), "")
  
  if (is.null(id_col_names)) {
    id_col_names = names(lst_col_specs)
  } else {
    errorIf(!isContainedIn(id_col_names, names(lst_col_specs)), "")
  }
  
  ddply(Data, id_col_names, function(D) {
    # Obtain the "header"
    hdr = buildValsWrapperHeader(D, lst_col_specs)
    D_out = func(D[, setdiff(colnames(D), names(lst_col_specs)), drop = FALSE], hdr)
    if (flag_include_header) {
      return (do.call(cbind, c(hdr, list(D_out))))
    } else {
      return (D_out)
    }
  }, ...)
}

# Given a named list, return a single data frame, adding a column corresponding to the names of the elements of the list 
# myLdply = function(lst, id_col_names) {
#   ldply(lst, function(D) {
#     # Obtain the "header"
#     hdr = buildValsWrapperHeader(D, lst_col_specs)
#     D_out = func(D[, setdiff(colnames(D), names(lst_col_specs)), drop = FALSE], hdr)
#     if (flag_include_header) {
#       return (do.call(cbind, c(hdr, list(D_out))))
#     } else {
#       return (D_out)
#     }
#   }, ...)
# }

# Rbind, with padding of columns of the second argument
rbind_padCols = function(A, B, B_padding_col_names, default_val = NA) {
  # Sanity checks
  errorIf(!isNumericMatrix(A) && !isDataFrameWithNamedCols(A), "")
  errorIf(!isNumericMatrix(B) && !isDataFrameWithNamedCols(B), "")
  errorIf(is.null(colnames(A)) || is.null(colnames(B)), "")
  errorIf(!isContainedIn(colnames(B), colnames(A)), "")
  
  if (isStringVector(B_padding_col_names)) {
    errorIf(!areDisjoint(B_padding_col_names, colnames(B)), "")
    errorIf(!setequal(c(colnames(B), B_padding_col_names), colnames(A)), "")
  } else if (is.null(B_padding_col_names)) {
    B_padding_col_names = setdiff(colnames(A), colnames(B))
  } else {
    myError("Invalid specification of the columns for padding")
  }
  
  return (rbind(A, do.call(cbind, c(as.list(mkNamedVector(B_padding_col_names, default_val = default_val)), list(B)))))
}

# Helper function used by functions 'dlply_valsWrapper', 'dlply_valsWrapper_fun', ...
buildValsWrapperHeader = function(D, lst_col_specs) {
  # Sanity checks
  errorIf(!isDataFrameWithNamedCols(D), "")
  errorIf(!isNamedList(lst_col_specs), "")
  hdr = applyFunToNamedContainer_lst(lst_col_specs, function(col_name, conv_fun) {
    val = unique(D[, col_name])
    if (!is.null(conv_fun)) val = conv_fun(val)
    return (mustBeScalar(val))
  })
  return (hdr)
}

# Analogous to 'dlply_valsWrapper', but the core function returns nothing
#   func(D_subset, hdr)
d_ply_valsWrapper_fun = function(Data, lst_col_specs, func, id_col_names = NULL, ...) { # func(D_subset, hdr)
  # Sanity checks
  errorIf(!isDataFrameWithNamedCols(Data), "")
  errorIf(!isNamedList(lst_col_specs), "")
  errorIf(!is.function(func), "")
  errorIf(!is.null(id_col_names) && !isStringVector(id_col_names), "")
  errorIf(!isContainedIn(names(lst_col_specs), colnames(Data)), "")
  
  if (is.null(id_col_names)) {
    id_col_names = names(lst_col_specs)
  } else {
    errorIf(!isContainedIn(id_col_names, names(lst_col_specs)), "")
  }
  
  d_ply(Data, id_col_names, function(D) {
    # Obtain the "header"
    hdr = buildValsWrapperHeader(D, lst_col_specs)
    func(D[, setdiff(colnames(D), names(lst_col_specs)), drop = FALSE], hdr)
  }, ...)
}

# A simple wrapper
discardCols_names = function(A, col_names) {
  # Sanity checks
  errorIf(!is.matrix(A) && !is.data.frame(A), "")
  errorIf(is.null(colnames(A)), "")
  errorIf(!isStringVector(col_names), "")
  errorIf(!isContainedIn(col_names, colnames(A)), "")
  return (A[, setdiff(colnames(A), col_names), drop = FALSE])
}

discardRows_names = function(A, row_names) {
  # Sanity checks
  errorIf(!is.matrix(A) && !is.data.frame(A), "")
  errorIf(is.null(rownames(A)), "")
  errorIf(!isStringVector(row_names), "")
  errorIf(!isContainedIn(row_names, rownames(A)), "")
  return (A[setdiff(rownames(A), row_names), , drop = FALSE])
}

discardElements_names = function(X, SS_names) {
  # Sanity checks
  errorIf(!is.list(X) && !is.vector(X), "")
  errorIf(is.null(names(X)), "")
  errorIf(!isStringVector(SS_names), "")
  errorIf(!isContainedIn(SS_names, names(X)), "")
  
  return (X[setdiff(names(X), SS_names)])
}

# Wrapper function, for convenience
useColForNamingRows = function(A, col_name) {
  errorIf(!isDataFrameWithNamedCols(A) && !is.matrix(A), "")
  errorIf(is.null(colnames(A)), "")
  errorIf(!isSingleStr(col_name), "")
  errorIf(!(col_name %in% colnames(A)), "")
  return (setRowNames(discardCols_names(A, col_name), as.character(A[, col_name]), flag_check_uniq = TRUE))
}

# 
# A function providing "strict mode" (as in Perl)
# See http://stackoverflow.com/questions/6216968/r-force-local-scope
# 
# Usage:
#   f <- useStrict(function(myvar) { return(myVar); } )
#
useStrict = function(f, pos = 2) eval(substitute(f), as.environment(pos))

#
# Helper functions for 'doSNOW': 'clusterExport_doSNOW' and 'createCluster_doSNOW'
# See: http://www.r-bloggers.com/parallelization-using-plyr-loading-objects-and-packages-into-worker-nodes/
#
clusterExport_doSNOW <- local({
  gets <- function(n, v) { assign(n, v, envir = .GlobalEnv); NULL }
  function(cl, list, envir = .GlobalEnv) {
    ## do this with only one clusterCall--loop on slaves?
    for (name in list) {
      clusterCall(cl, gets, name, get(name, envir = envir))
    }
  }
})
createCluster_doSNOW = function(noCores, logfile = "/dev/null", export = NULL, lib = NULL) {
  require(doSNOW)
  cl <- makeCluster(noCores, type = "SOCK", outfile = logfile)
  if(!is.null(export)) clusterExport_doSNOW(cl, export)
  if(!is.null(lib)) {
    l_ply(lib, function(dum) { 
      clusterExport_doSNOW(cl, "dum", envir = environment())
      clusterEvalQ(cl, library(dum, character.only = TRUE))
    })
  }
  registerDoSNOW(cl)
  return(cl)
}

# For convenience
createCluster_doSNOW_exportAll = function(N_cores, ...)
    createCluster_doSNOW(N_cores, export = ls(), lib = getNamesLoadedPackages(), ...)

# Function copied from file 'my_utils.R' from folder 'Genotyping_data-NaCl_experiment' on July 15th, 2015
getIndsSepGivenGroupAssignment = function(group_assignment) {
  # Sanity checks
  errorIf(!isStringVector(group_assignment), "")
  errorIf(isEmpty(group_assignment), "")
  
  NN = length(group_assignment)
  return (1 + noNames(which(group_assignment[1:(NN - 1)] != group_assignment[2:NN])))
}

# Wrapper to apply binary function
seqPairwiseComparisons = function(X, fun) {
  errorIf(!is.list(X) || !is.vector(X), "")
  errorIf(isScalar(X), "")
  return (sapply(X[2:length(X)], function(x) fun(x, X[[1]])))
}

# Helper function to make sure that a set of objects are all identical
areAllIdenticalSets = function(lst) {
    # Sanity checls
    errorIf(!is.list(lst), "The input must be a list")
    return (all(seqPairwiseComparisons(lst, identical)))
}

allAllEqualSets = function(lst) {
    # Sanity checls
    errorIf(!is.list(lst), "The input must be a list")
    return (all(seqPairwiseComparisons(lst, setequal)))
}

# Wrapper functions for convenience
quantileLower = function(X, alpha = 0.05, ...) quantile(X, probs = 0.5 * alpha, ...)
quantileUpper = function(X, alpha = 0.05, ...) quantile(X, probs = 1 - 0.5 * alpha, ...)

# Wrapper for the chi-square bootstrap-based test
chisqBootpTest_wrapper = function(count_data, p_ref, n_bootp = 10000) {
  # Sanity checks
  errorIf(!isNamedNumericVector(count_data), "")
  errorIf(!isNamedNumericVector(p_ref), "")
  errorIf(!setequal(names(count_data), names(p_ref)), "")
  p_ref = p_ref[names(count_data)] # enforce the order
  return (chisq.test(count_data, p = p_ref, simulate.p.value = TRUE, B = n_bootp))
}

chisqTestFromTable_wrapper = function(Table, count_col_name, ref_freq_col_name, ...) {
  # Sanity checks
  errorIf((!is.data.frame(Table) && !is.matrix(Table)) || is.null(colnames(Table)), "")
  errorIf(!isSingleStr(count_col_name), "")
  errorIf(!isSingleStr(ref_freq_col_name), "")
  errorIf(!isContainedIn(c(count_col_name, ref_freq_col_name), colnames(Table)))
  count_data = Table[, count_col_name]
  p_ref = Table[, ref_freq_col_name]
  errorIf(any(is.na(count_data)) || any(is.na(p_ref)), "")
  return (chisq.test(count_data, p = p_ref, ...))
}

# relies on function 'chisqTestFromTable_wrapper'
chisqBootpTestFromTable_wrapper = function(Table, count_col_name, ref_freq_col_name, n_bootp = 10000) 
  chisqTestFromTable_wrapper(Table, count_col_name, ref_freq_col_name, simulate.p.value = TRUE, B = n_bootp)

# The names are set as: [row_name][sep][col_name]
# Note that the matrix returned has no row or column names
getNamingMatrix = function(X, sep = ".", flag_by_row = TRUE) {
  # Sanity checks
  errorIf(!isFullyNamedMatrix(X), "")
  row_names = rownames(X)
  col_names = colnames(X)
  
  if (flag_by_row) {
    return (do.call(rbind, lapply(row_names, function(ss) paste(ss, col_names, sep = sep))))
  } else {
    return (do.call(cbind, lapply(col_names, function(ss) paste(row_names, ss, sep = sep))))
  }
}

# Convert a matrix to a (named) vector, by combining the names of the rows and columns of the matrix
mat2Vec_combineNames = function(X, ...) {
  # Sanity checks
  errorIf(!isFullyNamedMatrix(X), "")
  
  XX = as.vector(X) # linearize
  SS = as.vector(getNamingMatrix(X, ...)) # linearize a matching matrix with the names of the elements
  
  return (mySetNames(XX, SS, flag_check_uniq_names = TRUE))
}

# Wrapper for function 'mat2Vec_combineNames', if the names need to be combined
mat2Vec_wrapper = function(X, flag_combine_names, ...) {
  # Sanity checks
  errorIf(!isFullyNamedMatrix(X), "")
  errorIf(!isLogicalScalar(flag_combine_names), "")
  
  if (flag_combine_names) {
    return (mat2Vec_combineNames(X, ...))
  } else {
    return (as.vector(X))
  }
}

# Wrapper function
checkIfProcessingShouldBeDone = function(output_data_file_path, flag_allow_overwrite) {
  # Sanity checks
  errorIf(!isSingleStr(output_data_file_path), "")
  errorIf(!isLogicalScalar(flag_allow_overwrite), "")
  
  # Check if the output data file already exists  
  if (file.exists(output_data_file_path)) {
    cat(sprintf("*** The output file '%s' already exists\n", output_data_file_path))
    errorIf(!flag_allow_overwrite, "Not allowed to overwrite the output file")
  }
}

# returns a numeric vector
myInterp = function(x, y, x_out, y_left = 0, y_right = NULL) {
  # Sanity checks
  errorIf(!isNumericVector(x), "Argument must be a numeric vector (x)")
  errorIf(!isNumericVector(y), "Argument must be a numeric vector (y)")
  errorIf(!isNumericVector(x_out), "Argument must be a numeric vector (x_out)")
  errorIf(length(x) != length(y), "Size mismatch")
  errorIf(isEmpty(x), "No data")
  
  if (is.null(y_left) || is.null(y_right)) {
    inds = order(x, decreasing = FALSE)
    if (is.null(y_left)) y_left = y[[inds[[1]]]] # left-most value
    if (is.null(y_right)) y_right = y[[inds[[length(inds)]]]] # right-most value
  }
  
  return (approx(x, y, xout = x_out, method = "linear", rule = 2, yleft = y_left, yright = y_right)$y)
}

# Wrapper function for convenience
# Allows to wait if the file is not available yet
loadDataFromFile = function(input_data_file_path, input_file_mode = "default", sleep_time_secs = 300, max_N_attempts = -1) {
  # Sanity checks
  errorIf(!isSingleStr(input_data_file_path), "")
  errorIf(!isSingleStr(input_file_mode), "")
  errorIf(!isNumericScalar(max_N_attempts), "")
  
  wait_msg = sprintf("Failed to load the input file '%s'. Sleeping for %d seconds...\n", input_data_file_path, sleep_time_secs)
  
  i = 0
  while (TRUE) {
    if (!file.exists(input_data_file_path)) {
      # File does not exist
      if (input_file_mode == "wait_if_no_input_file") {
        cat(wait_msg)
      } else {
        myError(sprintf("The input file '%s' does not exist", input_data_file_path))
      }
    } else {
      # Try loading the file with the bootstrap results
      loading_res = tryCatch({
        cat(sprintf("Loading input data from file '%s'\n", input_data_file_path))
        in_data = new.env(); load(input_data_file_path, envir = in_data)
        cat(sprintf("DONE\n"))        
        list(results = in_data, status = "OK")
      }, warning = function(w) {
        msg_str = sprintf("[WARNING] %s", w)
        myWarning(msg_str)
        return (list(results = NULL, status = msg_str))
      }, error = function(e) {
        msg_str = sprintf("[ERROR] %s", e)
        myWarning(msg_str)
        return (list(results = NULL, status = msg_str))
      })
      
      if (!is.null(loading_res$results)) {
        return (loading_res$results)
      } else {
        cat(wait_msg)
      }
    }
    i = i + 1
    if ((max_N_attempts > 0) && (i >= max_N_attempts))
      break
    Sys.sleep(sleep_time_secs) # sleep for 5 minutes
  }
  myError(sprintf("Failed to load file '%s'", input_data_file_path))
}

# Given a data frame having a column with values ('in_col_name'), add a new character column ('new_col_name'), setting NA's as empty strings
addValueLabelColToDataFrame = function(D, in_col_name, new_col_name, format_str = NULL) {
  # Sanity checks
  errorIf(!isDataFrameWithNamedCols(D), "")
  errorIf(!isSingleStr(in_col_name), "")
  errorIf(!isSingleStr(new_col_name), "")
  errorIf(!(in_col_name %in% colnames(D)), "")
  errorIf(new_col_name %in% colnames(D), "")
  
  X = D[, in_col_name]
  if (is.null(format_str)) {
    label_data = as.character(X)
  } else if (isSingleStr(format_str)) {
    label_data = sprintf(format_str, X)
  } else {
    myError("")
  }
  
  D = do.call(cbind, c(list(D), setNames(list(label_data), new_col_name), list(stringsAsFactors = FALSE)))
  D[, new_col_name] = as.character(D[, new_col_name]) # just to be sure
  D[is.na(D[, new_col_name]), new_col_name] = ""
  return (D)
}

stepwiseStrConcat = function(SS, n_per_block, within_block_sep = ", ", among_block_sep = "\n") {
  # Sanity checks
  errorIf(!isStringVector(SS), "")
  errorIf(!isPositiveNumericScalar(n_per_block), "")
  errorIf(!isSingleStr(within_block_sep), "")
  errorIf(!isSingleStr(among_block_sep), "")
  
  return (paste(sapply(divideIntoSubLists(SS, n_per_block), function(SS) paste(SS, collapse = within_block_sep)), collapse = among_block_sep))
}

# do nothing if the input is a standard data frame
convertDPlyrDF = function(D) {
  # Sanity checks
  errorIf(!is.data.frame(D), "")
  if ("tbl_df" %in% class(D)) D = as.data.frame(D)
  return (D)
}

castDFCols2Strings = function(D, col_names) {
  # Sanity checks
  errorIf(!isDataFrameWithNamedCols(D), "")
  errorIf(!isStringVector(col_names), "")
  errorIf(!isContainedIn(col_names, colnames(D)), "")
  errorIf(!isContainedIn(sapply(col_names, function(ss) class(D[, ss])), c("factor", "character")))
  
  for (ss in col_names) D[, ss] = as.character(D[, ss])
  return (D)
}

#
# General function to setup specs for plotting populations
# Specs = color, linetype, shape
# The default: color per generation, shape and linetype per replicate population
#
mkSpecsForPlotting = function(Sample_info_table, source_col_name = "name_sample", 
                              color_col_name = "generation", color_fun = ggplotColours, 
                              shape_col_name = "repl_pop", shape_fun = getShapesForPlotting, 
                              linetp_col_name = "repl_pop", linetp_fun = getLinetypesForPlotting) {
  # Sanity checks
  errorIf(!isDataFrameWithNamedCols(Sample_info_table), "")
  errorIf(!isContainedIn(c(source_col_name, color_col_name, shape_col_name, linetp_col_name), colnames(Sample_info_table)), "")
  
  mkSpecMap = function(fun, labels) mkMap(fun(length(labels)), labels) 
  
  label_2_color_map = mkSpecMap(color_fun, as.character(unique(Sample_info_table[, color_col_name])))
  label_2_shape_map = mkSpecMap(shape_fun, as.character(unique(Sample_info_table[, shape_col_name])))
  label_2_linetp_map = mkSpecMap(linetp_fun, as.character(unique(Sample_info_table[, linetp_col_name])))
  
  All_sample_names = as.character(Sample_info_table[, source_col_name])
  Plot_specs_table = data.frame(
    name_sample = All_sample_names, 
    sample_color = label_2_color_map[applyTableMap(All_sample_names, Sample_info_table, source_col_name, color_col_name, flag_col_to_as_character = TRUE)], 
    sample_shape = label_2_shape_map[applyTableMap(All_sample_names, Sample_info_table, source_col_name, shape_col_name, flag_col_to_as_character = TRUE)], 
    sample_linetp = label_2_linetp_map[applyTableMap(All_sample_names, Sample_info_table, source_col_name, linetp_col_name, flag_col_to_as_character = TRUE)], 
    row.names = NULL
  )
  
  return (list(Plot_specs_table = Plot_specs_table, 
               color_map = mkMapFromTable(Plot_specs_table, "name_sample", "sample_color", flag_col_to_as_character = TRUE), 
               shape_map = mkMapFromTable(Plot_specs_table, "name_sample", "sample_shape", flag_col_to_as_character = FALSE), 
               linetp_map = mkMapFromTable(Plot_specs_table, "name_sample", "sample_linetp", flag_col_to_as_character = TRUE)))
}

#
# Helper function to "expand" subsets of a table
#
# For example, consider a data frame 'Results', with columns: 
#   "Model_stats" (data frame), "test_sample_name" (string)
# then return a data frame obtained by concatenating: [test_sample_name, Model_stats]
#
expandTable = function(D, id_col_name, data_col_name, flag_id_as_char = FALSE) {
    # Sanity checks
    errorIf(!isDataFrameWithNamedCols(D), "")
    errorIf(!isSingleStr(id_col_name), "")
    errorIf(!isSingleStr(data_col_name), "")
    errorIf(!isContainedIn(c(id_col_name, data_col_name), colnames(D)), "")
    errorIf(id_col_name == data_col_name, "")
    
    Ids = D[, id_col_name]
    if (flag_id_as_char)
        Ids = as.character(Ids)
    
    Res = do.call(rbind, mapply(function(.id, .data) cbind(.id, .data), Ids, D[, data_col_name], SIMPLIFY = FALSE))
    
    colnames(Res)[[1]] = id_col_name
    return (mustBeAllUniqCols(Res))
}

# Convenience function for generating a figure with log10(abs(x))
setupDataForLog10AbsPlot = function(Data, col_names, min_abs_val = 1e-8) {
    # Sanity checks
    errorIf(!isDataFrameWithNamedCols(Data), "")
    errorIf(!isStringVector(col_names), "")
    errorIf(!isPositiveNumericScalar(min_abs_val), "")
    errorIf(!areAllValuesUnique(col_names), "")
    
    for (name in col_names) {
        x = Data[, name]
        errorIf(!isNumericVector(x), "")
        x = abs(x)
        x[x < min_abs_val] = min_abs_val
        Data[, name] = log10(x)
    }
    
    return (Data)    
}

getHeatmapBaseEdgesFromStats_log10Abs = function(X, p = c(0.01, 0.25, 0.5, 0.75, 0.99)) {
    # Sanity checks
    errorIf(!isNumericMatrix(X), "")
    return (getHeatmapBaseEdgesFromStats_log10Abs(as.vector(X), p = p))
}

getHeatmapBaseEdgesFromStats_log10Abs_vec = function(x, p = c(0.01, 0.25, 0.5, 0.75, 0.99)) {
    # Sanity checks
    errorIf(!isNumericVector(x), "")
    x = abs(x)
    x = x[!is.na(x)]
    x = x[x > 0]
    xx = log10(x)
    return (10 ^ (c(min(xx), quantile(xx, probs = p, names = FALSE), max(xx))))
}

# Helper function to return the value of R2
getR2Val = function(Y_true, Y_pred) {
    # Sanity checks
    errorIf(!isNumericVector(Y_true), "")
    errorIf(!isNumericVector(Y_pred), "")
    return (summary(lm(y_true ~ -1 + y_pred, data = data.frame(y_true = Y_true, y_pred = Y_pred)))$r.squared)
}

# Fast version to melt a (large) matrix, e.g., a correlation matrix
# Will return only the upper diagonal, by construction
# Code obtained from file 'Postdoc/Andrei/Hapl_corr_matrix/helper_process_plot.R' 
#   (function 'meltCorMat_upperDiag_fast') on November 17th, 2015
meltMat_upperDiag_fast = function(Rho_mat) {
    # Sanity checks
    errorIf(!isFullyNamedNumericMatrix(Rho_mat) || (!identical(rownames(Rho_mat), colnames(Rho_mat))), "")
    names_features = rownames(Rho_mat)
    N = length(names_features)
    
    library(data.table)
    DD = ldply(2:N, function(i_col) { # start from 2
        inds_rows = seq(1, i_col - 1)
        return (data.frame(i = inds_rows, 
                           j = i_col, 
                           rho_val = Rho_mat[inds_rows, i_col], 
                           row.names = NULL))
    }, .progress = "text")
    DD = mutate(DD, feature_i = names_features[i], feature_j = names_features[j])
    return (DD)
}

# MATLAB-like function
logspace = function(x, n = 200, x_thr = 1) 10 ^ seq(log10(max(x_thr, min(x))), log10(max(x)), length.out = n)

#
# Usage:
# For:
#   > D
#        name x y
#   1  thiago 0 1
#   2  santos 1 1
#   3 guzella 2 2
#
#   > DD = setValues(D %>% colsAsChar("name"), D$y == 1, list("name" = "???", x = NA))
#
# The second argument may be also be a function
#
setValues = function(D, bool_rows, lst_val_data) {
    # Sanity checks
    errorIf(!isDataFrameWithNamedCols(D), "Invalid data frame")
    if (is.function(bool_rows)) bool_rows = bool_rows(D)
    errorIf(!isLogicalVector(bool_rows) || (length(bool_rows) != nrow(D)), "Invalid boolean vector")
    errorIf(!isNamedList(lst_val_data) || !isContainedIn(names(lst_val_data), colnames(D)), "Invalid specification of the new values")
    
    for (i in seq(1, length(lst_val_data))) {
        col_name = names(lst_val_data)[[i]]
        val = lst_val_data[[i]]
        D[bool_rows, col_name] = val
    }
    return (D)
}

# Wrapper for do.call(rbind, lst), so that it can be used with the %>% operator
# Force with as.data.frame just to deal with the case of 'NULL'
doCallRbind_wrapper = function(lst) as.data.frame(do.call(rbind, lst))

# Sanity check for rbindlist (which does not check that the columns have the same order)
# Returns the argument itself if ok, otherwise raises an error
sanitizeForRbindList = function(lst) {
    # Sanity checks
    errorIf(!is.list(lst) || any(!sapply(lst, function(X) (isDataFrameWithNamedCols(X) || is.null(X)))), 
            "The input must be a list of data frames with valid column names")
    # Get rid of empty entries
    inds_discard = which(sapply(lst, function(X) is.null(X) || !nrow(X)))
    lst = lst[-inds_discard]
    if (isEmpty(lst))
        return (lst)
    errorIf(any(!sapply(lst, function(X) areAllValuesUnique(colnames(X)))), "Column names must be unique")
    target_col_names = colnames(lst[[1]])
    errorIf(is.null(target_col_names), "No column names")
    
    lapply(lst, function(D) {
        errorIf(!setequal(target_col_names, colnames(D)), 
                sprintf("Mismatch in column names: got '%s', expected '%s'", 
                        paste(colnames(D)), paste(target_col_names, collapse = ", ")))
        D[, target_col_names, drop = FALSE]
    })
}

# A shortcut for dcast() + melt()
padTable = function(D, header_cols, padding_col, value_col, fill_val = 0, ...) {
    # Sanity checks
    errorIf(!isDataFrameWithNamedCols(D), "")
    errorIf(!isStringVector(header_cols), "")
    errorIf(!isSingleStr(padding_col), "")
    errorIf(!isSingleStr(value_col), "")
    errorIf(!isScalar(fill_val), "")
    all_col_names = c(header_cols, padding_col, value_col)
    errorIf(!areAllValuesUnique(all_col_names) || !isContainedIn(all_col_names, colnames(D)), "")
    
    D %>%
        dcast(as.formula(sprintf("%s ~ %s", paste(header_cols, collapse = " + "), padding_col)), 
              value.var = value_col, fill = fill_val, ...) %>% # Padding
        melt(id.vars = header_cols, variable.name = padding_col, value.name = value_col)
}

# Write compressed data to file
writeTable_bz2 = function(Data, file_path, 
                          n_threads = NULL, 
                          compression = 6, 
                          sep = "\t", use_quote = FALSE, row.names = FALSE, col.names = TRUE) {
    ff = bzfile(file_path, open = "wb", compression = compression)
    write.table(Data, file = ff, sep = sep, quote = use_quote, row.names = row.names, col.names = col.names)
    close(ff)
}

# Read compressed data from file
readTable_bz2 = function(file_path, 
                         header = TRUE, sep = "\t") {
    read.table(file_path, header = header, sep = sep)
}

reOrderCols = function(A, new_order) {
    cols_1 = intersect(new_order, colnames(A))
    cols_2 = setdiff(colnames(A), new_order)
    return (A[, c(cols_1, cols_2)])
}

getECDF <- function(X) {
    errorIf(!isNumericVector(X), "")
    X = X[!is.na(X)]
    data.frame(val = sort(unique(X))) %>% 
        mutate(n_points = sapply(val, function(thr) sum(X <= thr)), 
               freq_points = n_points / length(X))
}

getRevECDF <- function(X) {
    errorIf(!isNumericVector(X), "")
    X = X[!is.na(X)]
    data.frame(val = sort(unique(X))) %>% 
        mutate(n_points = sapply(val, function(thr) sum(X > thr)), 
               freq_points = n_points / length(X))
}

# Convenience function, equivalent to 'ldply(lst, .id = col_name)', but 
# inputs are data tables, and returns a data table
# No format sanity checks are done...
namedDTList2DataTable = function(lst, col_name) {
    # Basic sanity checks
    errorIf(!isNamedList(lst), "")
    errorIf(!isSingleStr(col_name), "")
    mapply(function(.name, DT) DT[, c(col_name) := list(.name)], 
           names(lst), lst, SIMPLIFY = FALSE) %>%
        rbindlist(use.names = TRUE)
}

# Turn a map into a data farme
map2DF = function(x, name_col, val_col) {
    stopifnot(isValidMap(x))
    stopifnot(isSingleStr(name_col))
    stopifnot(isSingleStr(val_col))
    stopifnot(name_col != val_col)
    data.frame(names(x), noNames(x)) %>%
        setColNames(c(name_col, val_col))
}

getCol = function(D, col_name, ...) {
    stopifnot(is.data.frame(D) || is.data.table(D))
    stopifnot(isSingleStr(col_name))
    stopifnot(col_name %in% colnames(D))
    if (is.data.table(D)) D[[col_name]] else D[, col_name, ...]
}

mustBeNoNA = function(D, msg = "Missing data is present") {
    stopifnot(is.data.frame(D) || is.data.table(D))
    D %>%
        mustBe(function(DD)
            all(apply(DD, 2, function(x) !any(is.na(x)))), msg)
}
