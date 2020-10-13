#
# Helper code for formatting the data
#
# Created on April 12th, 2016
#

source(file.path("Common", "my_utils.R"))

getEnvRenamingMap = function() c("NGM" = "NGM", "NGM+NaCl" = "NaCl")

formatEnvironment = function(D, ren_map = getEnvRenamingMap()) 
    D %>% 
        renameAndEnforceOrder("Environment", ren_map)

standardizeLineNames = function(x, specs = c("G140A6" = "A6140", "G50MOVA1" = "GA150", "G50MOVA2" = "GA250", "G50MOVA4" = "GA450")) {
    errorIf(!isStringVector(x), "")
    errorIf(!isValidMap(specs), "")
    for (i in seq_along(specs)) {
        x = gsub(paste0("^", names(specs)[[i]]), specs[[i]], x)
    }
    x
}
