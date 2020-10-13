#
# Helper code with definitions related to the assay locations (Paris or Lisbon)
#
# Created on May 6th, 2017
#

source(file.path("Common", "my_utils.R"))

getLocationDefs = function()
    data.frame(location_label = c("Paris", "Lisbon")) %>%
        mutate(location_val = seq(0, n() - 1))

# For the data from Bruno's dataset
getDataFromBrunoLocationYearMap = function()
    c("2012" = "Lisbon", 
      "2013" = "Lisbon", 
      "2014" = "Lisbon", 
      "2015" = "Paris", 
      "2016" = "Paris")

mkExperimentLocationSpecsTable = function(Header_data, exper_index_col_name = "experiment_ind") {
    stopifnot(isDataFrameWithNamedCols(Header_data))
    stopifnot(isSingleStr(exper_index_col_name))
    stopifnot(all(c(exper_index_col_name, "dataset_name", "date_str") %in% colnames(Header_data)))
    
    Header_data[, c(exper_index_col_name, "dataset_name", "date_str")] %>% 
        mustBe(function(D) areAllValuesUnique(D[, exper_index_col_name])) %>%
        # Retrieve the year
        mutate(year_str = substr(as.character(date_str), 1, 4)) %>%
        distinct() %>%
        mutate(location_label = ifelse(
            dataset_name == "Data_from_Bruno", 
            applyMap(as.character(year_str), getDataFromBrunoLocationYearMap(), 
                     raise_err = FALSE, ret_input_on_err = FALSE, default_val = as.character(NA)), 
            "Paris"
        )) %>%
        mustBe(function(D) all(!is.na(unlist(D))))
}
