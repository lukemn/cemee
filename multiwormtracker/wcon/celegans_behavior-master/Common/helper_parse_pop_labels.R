#
# Helper code to parse the population labels (standard format, e.g., GA150L56)
#
# Created on October 31st, 2016
# Modified 20/11/2018 to add all the new populations

source(file.path("Common", "my_utils.R"))

getPopReprSysData = function(){

vect_A <- c("A0","A00",paste0(c("A4","A5","A6"),rep(c(60,140),each=3)))
vect_CA <- paste0(c("CA1","CA2","CA3"),rep(c(50,100),each=3))
vect_GA <- paste0(c("GA1","GA2","GA4"),rep(c(50),each=3))
vect_GT <- paste0(c("GT1","GT2"),rep(c(50),each=3))
vect_GM <- paste0(c("GM1","GM3"),rep(c(50),each=3))
vect_WI <- c("N2","JU400","N2anc","JU345","AB1","CB4858","MY1","PX174","JU319","CB4507","CB4852","CB4856","PB306","MY16","RC301","PX179","CB4855")

all_pops <- c(vect_A, vect_CA, vect_GA, vect_GT, vect_GM)
all_pops_noM <- paste0(all_pops,"noM")

all_pops_rep <- rep(c("A","A","A","T","M"), c(
length(vect_A), length(vect_CA), length(vect_GA), length(vect_GT), length(vect_GM)))

    data.frame(population = c(all_pops, all_pops_noM, vect_WI) ,
               repr_system = c(all_pops_rep,all_pops_rep,rep("wild_isolate",length(vect_WI))))
}
getPopTreatData = function(){

vect_A <- c("A0","A00",paste0(c("A4","A5","A6"),rep(c(60,140),each=3)))
vect_CA <- paste0(c("CA1","CA2","CA3"),rep(c(50,100),each=3))
vect_GA <- paste0(c("GA1","GA2","GA4"),rep(c(50),each=3))
vect_GT <- paste0(c("GT1","GT2"),rep(c(50),each=3))
vect_GM <- paste0(c("GM1","GM3"),rep(c(50),each=3))
vect_WI <- c("N2","JU400","N2anc","JU345","AB1","CB4858","MY1","PX174","JU319","CB4507","CB4852","CB4856","PB306","MY16","RC301","PX179","CB4855")

all_pops <- c(vect_A, vect_CA, vect_GA, vect_GT, vect_GM)
all_pops_noM <- paste0(all_pops,"noM")

all_pops_treatment <- rep(c("ancestral_0","control","gradual","gradual","gradual"), c(
length(vect_A), length(vect_CA), length(vect_GA), length(vect_GT), length(vect_GM)))


    data.frame(population = c(all_pops, all_pops_noM, vect_WI), 
               treatment = c(all_pops_treatment,all_pops_treatment,rep("wild_isolate", length(vect_WI))))
}
#getReplPopData = function(){	
#    data.frame(population = c("A6140", "CA150", "CA250", "CA350", "GA150", "GA250", "GA450", "GT150", "GT250", "GM150", "GM350", "CA1100", "CA2100", "CA3100"), 
#              repl_pop_std = c("A6", "A1", "A2", "A3", "A1", "A2", "A4", "T1", "T2", "M1", "M3", "A1", "A2", "A3"), 
#               repl_pop_v2 = c("A6", "CA1", "CA2", "CA3", "GA1", "GA2", "GA4", "GT1", "GT2", "GM1", "GM3", "CA1", "CA2", "CA3"))
#}

# Parse strings of the format: 'GA150L56'
parsePopLabel = function(Pop_labels, pop_treatment_data = getPopTreatData(), flag_relaxed = FALSE) {
    # Sanity checks
    errorIf(!isStringVector(Pop_labels), "")
    
    # Parse the 'pop_label' from variable 'final_Header_data'
    lst_pop_label_tokens = strsplit(Pop_labels, "L", fixed = TRUE)
    errorIf(!isContainedIn(sapply(lst_pop_label_tokens, length), c(1, 2)), "Invalid pop_label tokens")
    
    Pop_header_table = data.frame(
        pop_label = Pop_labels, 
        population = sapply(lst_pop_label_tokens, function(S) S[[1]]), 
        line_label_id = sapply(lst_pop_label_tokens, function(S)
            ifelse(length(S) == 2, S[[2]], as.character(NA)))
    ) %>%
        mutate(line_num_id = as.integer(as.character(line_label_id)))
    if (!flag_relaxed)
        Pop_header_table = Pop_header_table %>%
            mustBe(function(D) all(is.na(D$line_num_id) == is.na(D$line_label_id)))
    Pop_header_table %>%
        mutate(line_label_id = ifelse(is.na(line_label_id), as.character(NA), paste0("L", line_label_id))) %>%
        left_join(pop_treatment_data) %>%
        mustBe(function(D) all(!is.na(D$treatment)))
}

