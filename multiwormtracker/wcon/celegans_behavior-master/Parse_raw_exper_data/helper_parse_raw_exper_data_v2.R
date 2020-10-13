#
# New version of the code for parsing the raw data, out of the folders for the new experiments, 
# in which the Bash script 'run_choreography.sh' was used to run choreography (i.e., 
# in which the 'phenotype_names.txt' file is generated)
#
# This new version assumes a minimalistic format: no "covariates" included in the actual data
#
# Created on June 15th, 2016, out of the devel code from script 'scrap_parseRawData_jun15th2016.R'
#

source(file.path("Common", "my_utils.R"))

#
# Main processing function for a single experiment folder
# Assumes that the input path is in the format "[...]/[plate_folder_name]", i.e., without ending in "Unified raw data"
#
# Arguments:
#	'Data_location_str': string used to describe where the data is located
#
parseDataFromFilesInExperimentFolder = function(base_folder_path, verbose = TRUE) {
    # Sanity checks
    errorIf(!isSingleStr(base_folder_path), "The base folder path must be a single string")
    
    # Look for the file with the names of the phenotypes
    phen_names_file_path = file.path(base_folder_path, "phenotype_names.txt")
    errorIf(!file.exists(phen_names_file_path), "No file with the names of the phenotypes")
    phen_names = readLines(phen_names_file_path)
    input_col_names = c("time", strsplit(phen_names, ",")[[1]])
    
    # For backwards compatibility
    col_renaming_map = c(
        'time' = 'Time', 'id' = 'IndWorm', 'persistence' = 'Persistence', 'area' = 'Area', 
        'speed' = 'Velocity', 'angular' = 'TurningRate', 'length' = 'RegressedLength', 
        'width' = 'RegressedWidth', 'aspect' = 'Aspect', 'midline' = 'Length', 
        'morphwidth' = 'Width', 'bias' = 'Bias', 'pathlen' = 'NetDistanceTravelled', 
        'curve' = 'Curvature', 'loc_x' = 'XPosition', 'loc_y' = 'Yposition', 
        'orient' = 'BodyOrientation', 'area:jitter' = 'AreaJitter', 
        'midline:jitter' = 'LengthJitter', 'morphwidth:jitter' = 'WidthJitter', 
        'loc_x:jitter' = 'XPositionJitter', 'loc_y:jitter' = 'YPositionJitter', 'kink' = 'kink'
    )
    errorIf(!isValidMap(col_renaming_map) || !areAllValuesUnique(col_renaming_map), "Invalid column renaming map")
    
    data_col_names = noNames(applyMap(input_col_names, col_renaming_map))
    
    file_list = list.files(path = base_folder_path, pattern = glob2rx("*.dat"), full.names = TRUE)
    file_list = file_list[!sapply(file_list, function (file_path) { return (file.info(file_path)$isdir) } )]
    n_files = length(file_list)
    
    if (verbose) cat(sprintf("Will process %d files in folder '%s'\n", n_files, base_folder_path))
    
    if (n_files) {
        # First, parse all the filenames, to retrieve information on the line name
        # and track number (supposedly)
        # This is a N x 2 data frame, with the values as strings
        base_file_names = basename(file_list)
        errorIf(!areAllValuesUnique(base_file_names), "The base filenames should necessarily be unique")
        file_names_token_info = setRowNames(do.call(rbind, lapply(base_file_names, function (file_name_i) {
            # kludge to allow '.' in file prefixes
            i = gsub('.dat', '', file_name_i)
            tokens = strsplit(i, '.', fixed=T)[[1]]
            file_id = tokens[[length(tokens)]]
            data_id = gsub(paste0('.', file_id), '', i)
            #tokens = strsplit(file_name_i, ".", fixed = TRUE)[[1]]
            #errorIf(length(tokens) != 3, sprintf("Expected to read 3 tokens from filename '%s', but got %d", file_name_i, length(tokens)))
            #data_id = tokens[[1]] # line id and environment
            #file_id = tokens[[2]]
            errorIf(is.null(data_id) || is.null(file_id), sprintf("Unable to parse patterns from filename '%s'", file_name_i))
            return (c(data_id = data_id, file_id = file_id)) 
        } )), base_file_names)
        errorIf(nrow(file_names_token_info) != n_files, "Size mismatch in the data storing info on the file names")
        
        # Make sure all the data ids are the same
        uniq_data_ids = unique(file_names_token_info[, "data_id"])
        errorIf(!isScalar(uniq_data_ids), sprintf("Expected a single file id, but got: {'%s'}", paste(uniq_data_ids, collapse = "', '")))
        data_id = uniq_data_ids
        
        file_ids = file_names_token_info[, "file_id"] %>%
            as.character() %>%
            as.numeric() %>%
            mustBe(function(x) all(is.finite(x)))
        
        library(data.table)
        Data = mapply(function(file_path, file_id) {
            # Load the file data
            File_data = read.table(file_path, header = FALSE)
            errorIf(ncol(File_data) != length(data_col_names), "Invalid data")
            File_data = setColNames(File_data, data_col_names, flag_check_uniq = TRUE)
            id = mustBeScalar(unique(File_data$IndWorm))
            errorIf(id != file_id, "Mismatch in the id")
            File_data
        }, file_list, file_ids, SIMPLIFY = FALSE) %>%
            rbindlist() %>%
            as.data.frame()
        
        # Replace non-finite values with NA
        for (i in seq(1, ncol(Data))) {
            if (class(Data[, i]) %in% c("numeric", "integer", "logical"))
                Data[!is.finite(Data[, i]), i] = NA
        }
        
        return (list(file_list = file_list, 
                     data_id = data_id, 
                     Data = Data))
    } else {
        myWarning(sprintf("No files found in folder '%s'...", base_folder_path))
        return (NULL)
    }
}

#
# Perform all operations related to a folder:
# 	1) parse the path
# 	2) process the files and the data
# 	3) save the results to a file
#
parseDataFromExperimentFolder = function(folder_path, base_output_folder_path, verbose = TRUE, print_single_summary_msg = FALSE, 
                                         save_results = TRUE, return_results = FALSE) {
    # Get information on the lower-level folder
    upper_dir_name = basename(folder_path)
    desc_folder_name = basename(dirname(folder_path))
    
    data_name_token = sprintf("%s-%s", gsub("-", "_", desc_folder_name), gsub("-", "_", upper_dir_name))
    if (verbose) cat(sprintf("*** Processing folder '%s' ('%s')\n", folder_path, data_name_token))
    
    # Generate the output data file name
    output_data_file_path = file.path(base_output_folder_path, sprintf("Parsed_data-%s.RData", data_name_token))
    if (verbose) cat(sprintf("Processing folder '%s' ('%s')\n", folder_path, data_name_token))
    
    # Parse the date and time
    folder_name_tokens = strsplit(basename(folder_path), "_")[[1]]
    errorIf(length(folder_name_tokens) != 2, 
            sprintf("Expected to read 2 tokens from folder name '%s', but got %d", 
                    basename(folder_path), length(folder_name_tokens)))
    data_date = folder_name_tokens[[1]]
    data_time = folder_name_tokens[[2]]
    
    if (!file.exists(output_data_file_path) || (!save_results)) {
        processing_start = Sys.time()
        if (verbose) cat(sprintf("Processing started on '%s'\n", processing_start))
        results = parseDataFromFilesInExperimentFolder(folder_path, verbose)
        processing_end = Sys.time()
        if (verbose) cat(sprintf("Processing finished on '%s'\n", processing_end))
        
        if (save_results) {
            if (verbose) cat(sprintf("Saving the data to file '%s'...\n", output_data_file_path))
            save(folder_path, base_output_folder_path, upper_dir_name, desc_folder_name, data_name_token, 
                 data_date, data_time, 
                 results, 
                 processing_start, processing_end, 
             file = setupForData(output_data_file_path, quiet = TRUE), compression_level = 9, precheck = TRUE)
            if (verbose) cat(sprintf("Data saved to file '%s'\n", output_data_file_path))
        }
        if (print_single_summary_msg) cat(sprintf("Folder '%s' -> %s -> %s [%s, %s]\n", folder_path, data_name_token, output_data_file_path, processing_start, processing_end))
        
        if (return_results) {
            return (list(
                upper_dir_name = upper_dir_name, 
                desc_folder_name = desc_folder_name, 
                data_name_token = data_name_token, 
                data_date = data_date, 
                data_time = data_time, 
                results = results
            ))
        }
    } else {
        if (print_single_summary_msg) cat(sprintf("Skipping processing: output file '%s' already exists\n", folder_path, output_data_file_path))
        if (return_results) {
            return (list(
                upper_dir_name = upper_dir_name, 
                desc_folder_name = desc_folder_name, 
                data_name_token = data_name_token, 
                data_date = data_date, 
                data_time = data_time, 
                results = NULL
            ))
        }
    }
}
