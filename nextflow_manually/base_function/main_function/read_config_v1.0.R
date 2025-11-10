# this script is for R to read config file
# author: laojp
# time: 2023.08.08
# position: SYSUCC bioinformatic platform
# workflow: diffbind---chipseeker
# usage
#   process_read_config(config_atacseq,[string: the script name to read parameter])


process_read_config <- function(config_atacseq, script_name){
    # read config
    print("read configure file and create parameter")
    if(length(commandArgs(trailingOnly = TRUE))<1){
        stop("please input config file")
    } else if(length(commandArgs(trailingOnly = TRUE))<2){
        stop("please input the script name to read")
    }
    readLines(commandArgs(trailingOnly = TRUE)[1]) %>%
        {
            ifelse(all(grep(x = .,pattern = "^#+")<=grep(x=.,pattern = script_name)),
                .[grep(x=.,pattern = script_name):length(.)],
                .[grep(x=.,pattern = script_name):(grep(x = .,pattern = "^#+")[grep(x = .,pattern = "^#+")>grep(x=.,pattern = script_name)][1])])
        } %>%
        grep("^[^#[:space:]]", ., value = TRUE) %>%
        for (line in .) {
            if (length(strsplit(line, " = ")[[1]]) == 2) {
                assign(trimws(strsplit(line, " = ")[[1]][1]), 
                    trimws(strsplit(line, " = ")[[1]][2]), 
                    envir = globalenv())
            cat(paste(trimws(strsplit(line, " = ")[[1]][1]), " reading as ", 
                trimws(strsplit(line, " = ")[[1]][2]), "\n", sep = ""))
        }
        else{stop(paste(trimws(strsplit(line, " = ")[[1]][1]), " not defination \n", sep = ""))}
    }
}


