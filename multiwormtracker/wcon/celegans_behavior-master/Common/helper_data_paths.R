#
# Path specifications, so that data can be "shared" across computers
#

source(file.path("Common", "my_utils.R"))

getRootDataPathSpecs = function() 
    c("gevpc05" = file.path("/", "users", "gev", "mallard", "behavioral_analysis","celegans_behavior-master"), 
      "gevpc10" = file.path("/", "users", "gev", "mallard", "behavioral_analysis","celegans_behavior-master"))


