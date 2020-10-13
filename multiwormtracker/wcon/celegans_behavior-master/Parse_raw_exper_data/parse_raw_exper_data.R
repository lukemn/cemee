#
# Parse raw experimental data on behavior data
#
# The experimental data must have been "generated" via the new 
# 'run_choreography.sh' script (i.e., the file 'phenotype_names.txt' must be 
# available; see helper script 'helper_parse_raw_exper_data_v2.R')
#
# The code has a slightly different "flow" than normal, because of the need to set 
# the working directory via the command line argument '--project-folder-path'
#
# Created on July 5th, 2016, based on the code from 
#   script 'Test_scripts/parse_test_data_acquired_jun10th2016.R'
# Modified on September 18th, 2016, to accept command line arguments
#

rm(list = ls())

prj_folder_path = NULL
input_folder_path = NULL
root_output_folder = NULL

library(optparse)
arg_opts = parse_args(OptionParser(option_list = list(
    make_option(c("-p", "--project-folder-path")), 
    make_option(c("-i", "--input-folder-path")),
    make_option(c("-o", "--output-root-folder-path"))
)))

if ("project-folder-path" %in% names(arg_opts)) prj_folder_path = arg_opts[["project-folder-path"]]
if ("input-folder-path" %in% names(arg_opts)) input_folder_path = arg_opts[["input-folder-path"]]
if ("output-root-folder-path" %in% names(arg_opts)) root_output_folder = arg_opts[["output-root-folder-path"]]

# ---

# Set the working directory, and move on to parse the data
if (!is.null(prj_folder_path)) setwd(prj_folder_path)

sys_info = Sys.info()
source(file.path("Parse_raw_exper_data", "helper_parse_raw_exper_data_v2.R"))
session_info = sessionInfo()

errorIf(is.null(input_folder_path), "No input folder")
errorIf(is.null(root_output_folder), "No output folder")

# Parse the data, and return the results
parsed_data = parseDataFromExperimentFolder(
    input_folder_path,
    "/tmp",
    save_results = FALSE,
    return_results = TRUE)

#
# Save the data
#
upper_dir_name = parsed_data$upper_dir_name
desc_folder_name = parsed_data$desc_folder_name
data_name_token = parsed_data$data_name_token
data_date = parsed_data$data_date
data_time = parsed_data$data_time
results = parsed_data$results

data_token = parsed_data$data_name_token
data_id = parsed_data$results$data_id
errorIf(!isSingleStr(data_token) || !isSingleStr(data_id), "Invalid format")

output_data_file_path = file.path(
    root_output_folder,
    sprintf("Parsed_data-%s-%s.RData", data_token, data_id))

if (!file.exists(output_data_file_path)){
    save(output_data_file_path, sys_info, session_info,
     input_folder_path,
     upper_dir_name, desc_folder_name, data_name_token,
     data_date, data_time,
     results,
     file = setupForData(output_data_file_path),
     compression_level = 9, precheck = TRUE)
}
