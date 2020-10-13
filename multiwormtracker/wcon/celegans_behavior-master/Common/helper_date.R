#
# Helper code for dealing with dates (such as the acquisition date)
#
# Created on April 12th, 2016
#

source(file.path("Common", "my_utils.R"))

parseAcqDate = function(SS) {
    # Parse the date strings
    sapply(strsplit(as.character(SS), ""), function(SS) {
        errorIf(length(SS) != 8, sprintf("Invalid data: {%s}", paste(SS, collapse = ";")))
        year_str = paste(SS[1:4], collapse = "")
        month_str = paste(SS[5:6], collapse = "")
        day_str = paste(SS[7:8], collapse = "")
        as.Date(sprintf("%s-%s-%s", year_str, month_str, day_str))
    })
}

# New version, created because function 'parseAcqDate' converts the Date values to numeric values, 
# complicating subsequent processing steps
# This version returns as a vector of strings, which can be cast as dates via as.Date
parseAcqDate_v2 = function(SS) {
    # Parse the date strings
    sapply(strsplit(as.character(SS), ""), function(SS) {
        errorIf(length(SS) != 8, sprintf("Invalid data: {%s}", paste(SS, collapse = ";")))
        year_str = paste(SS[1:4], collapse = "")
        month_str = paste(SS[5:6], collapse = "")
        day_str = paste(SS[7:8], collapse = "")
        format(as.Date(sprintf("%s-%s-%s", year_str, month_str, day_str)))
    })
}
