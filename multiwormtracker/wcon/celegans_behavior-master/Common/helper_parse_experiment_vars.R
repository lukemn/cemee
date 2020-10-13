#
# Helper code for parsing experiment "variables" (time, date, etc)
#
# Created on October 31st, 2016
# Modified JAn 2018, added convertDateStrToMWTFormat_FM

source(file.path("Common", "my_utils.R"))

# Convert date string formatted as '21/11/16' to '20161121'
convertDateStrToMWTFormat = function(SS) {
    errorIf(!isStringVector(SS), "")
    SS %>%
        strsplit("/", fixed = TRUE) %>%
        mustBe(function(lst) all(sapply(lst, length) == 3)) %>%
        sapply(function(X) paste0("20", paste(rev(X), collapse = "")))
}


# Convert date string formatted as '21/11/16' to '20161121'
convertDateStrToMWTFormat_FM = function(SS) {
    errorIf(!isStringVector(SS), "")
    SS %>%
        strsplit("/", fixed = TRUE) %>%
        mustBe(function(lst) all(sapply(lst, length) == 3)) %>%
        sapply(function(X) paste0(paste(rev(X), collapse = "")))
}


# Returns a numeric value (in seconds)
parseTimeValsFromNumTuple = function(X) {
    if (isNumericVector(X) && (length(X) == 3) && (all(is.finite(X)))) {
        t_h = X[[1]]
        t_m = X[[2]]
        t_s = X[[3]]
        errorIf(!(t_h %in% 0:23) || !(t_m %in% 0:59) || !(t_s %in% 0:60), "")
        t_h * 3600 + t_m * 60 + t_s
    } else {
        NA
    }
}

parseExperSummDataTimeStr = function(SS) {
    if (length(SS) >= 6) {
        list(SS[1:2], SS[3:4], SS[5:6]) %>%
            sapply(function(SSS) as.integer(paste(SSS, collapse = ""))) %>%
            parseTimeValsFromNumTuple()
    } else {
        NA
    }
}
