#' @title An R Package for Assigning Results from MALDI-TOF Real-Time Classification (RTC) of Bacterial Isolates
#'
#' @description This is a function for RTC result assignment. MALDI-TOF results for microbial isolates are generated in a PDF format, with embedded table of results.
#' Four lines of results for each isolate is provided.
#' This package takes the result table in a text file format, applies predetermined logic for processing those results and assignes classes, and tabulates the assignment.
#'
#' What this program does:
#'
#' Step1: Accepts the raw data and cleans it up.
#'
#' Step2: Creates a flat file (each id on one line). This is the intermed file object.
#'
#' Step3: Uses part of the intermediate file as final input, to generate the result file, based on logic provided.
#'
#' Step4: Writes a text file as output with the results
#'
#' Input: A text file (tax separated) with 3 columns, id, organism, and score.
#'
#' Note: Orgnasim id should be in this format: Genus species_strain
#'
#' To run this function from command line, ensure you have the source code file in your working directory & run the following to load:
#' source("rtcassignment.R")
#'
#'
#' Computational Requirements
#' This program is written for easy implementation, not speed.
#' It uses base R code, no libraries required, so can be run on any machine with R installed.
#'
#' @param rtc_input_name Path & name your input file.
#' @param result_name Path & name for your desired output file
#' @return Writes a text file as output with the assignemnt results
#'
#' @examples
#' rtc_to_results("inst/extdata/example_data.txt", "example_assignment.txt")
#' @export
#'
rtc_to_results <- function(rtc_input_name, result_name){

  # Load data externally
  rtc_data <- read.table(file = rtc_input_name, header = TRUE, sep="\t", stringsAsFactors = FALSE)

  # Clean column names
  colnames(rtc_data) <- c("name", "id", "score")

  # Clean names
  rtc_data$id <- gsub("Non", "No", rtc_data$id)
  rtc_data$id <- gsub("no", "No", rtc_data$id)

  # Flatten rtc_data file

  # define matrix
  iterations <- length(unique(rtc_data$name))
  variables = 9 # (1 name + 4 pairs)
  intermed <- matrix(ncol=variables, nrow=iterations)

  # Loop to flatten file
  for (i in 1:iterations){
    intermed[i,1] <- rtc_data$name[4*i-3]
    intermed[i,2] <- rtc_data$id[4*i-3]
    intermed[i,3] <- rtc_data$score[4*i-3]
    intermed[i,4] <- rtc_data$id[4*i-3+1]
    intermed[i,5] <- rtc_data$score[4*i-3+1]
    intermed[i,6] <- rtc_data$id[4*i-3+2]
    intermed[i,7] <- rtc_data$score[4*i-3+2]
    intermed[i,8] <- rtc_data$id[4*i-3+3]
    intermed[i,9] <- rtc_data$score[4*i-3+3]
  }

  # Save as df
  intermed <- data.frame(intermed, stringsAsFactors=FALSE)
  colnames(intermed) <- c("name", "id1", "score1", "id2", "score2", "id3", "score3", "id4", "score4")
  #View(intermed)

  # Extract genus (split string by space and select first part)
  intermed$gen1 <- sapply(strsplit(intermed$id1, " "), "[", 1)
  intermed$gen2 <- sapply(strsplit(intermed$id2, " "), "[", 1)
  intermed$gen3 <- sapply(strsplit(intermed$id3, " "), "[", 1)
  intermed$gen4 <- sapply(strsplit(intermed$id4, " "), "[", 1)

  # Extract species (split string by space and select second part)
  intermed$sp1 <- sapply(strsplit(intermed$id1, " "), "[", 2)
  intermed$sp2 <- sapply(strsplit(intermed$id2, " "), "[", 2)
  intermed$sp3 <- sapply(strsplit(intermed$id3, " "), "[", 2)
  intermed$sp4 <- sapply(strsplit(intermed$id4, " "), "[", 2)

  # Clean species name (remove strain info) (remove part after underscore)
  intermed$sp1 <- sapply(strsplit(intermed$sp1, "_"), "[", 1)
  intermed$sp2 <- sapply(strsplit(intermed$sp2, "_"), "[", 1)
  intermed$sp3 <- sapply(strsplit(intermed$sp3, "_"), "[", 1)
  intermed$sp4 <- sapply(strsplit(intermed$sp4, "_"), "[", 1)

  # Subset rtc_data file (select essential fields)
  rtc_data2 <- intermed[, c("name", "gen1", "sp1", "score1", "gen2", "sp2", "score2", "gen3", "sp3", "score3", "gen4", "sp4", "score4")]
  #View(rtc_data2)


  # Define final output matrix
  output <- matrix(ncol=3, nrow=nrow(rtc_data2))

  # loop for each row
  for (i in 1:nrow(rtc_data2)){

    # Set name
    output[i,1] <- as.character(rtc_data2$name[i])

    # Extract species info
    sp <- vector()
    sp <- c(rtc_data2$sp1[i], rtc_data2$sp2[i], rtc_data2$sp3[i], rtc_data2$sp4[i])

    # Extract genus info
    gen <- vector()
    gen <- c(rtc_data2$gen1[i], rtc_data2$gen2[i], rtc_data2$gen3[i], rtc_data2$gen4[i])

    # Extract score info
    score <- vector()
    score <- c(rtc_data2$score1[i], rtc_data2$score2[i], rtc_data2$score3[i], rtc_data2$score4[i])

    # Create otu_ids
    otu_id <- vector()
    otu_id <- c(paste(gen[1], sp[1], sep="_"), paste(gen[2], sp[2], sep="_"), paste(gen[3], sp[3], sep="_"), paste(gen[4], sp[4], sep="_"))

    # Check for no peak (sp field will be peak)
    if (length(which(sp == "peaks")) > 2){
      output[i,2] <- "No Peaks"
      output[i,3] <- "NA"
    }

    # Check for no reliable (sp field will be Reliable)
    else if (length(which(sp == "Reliable")) > 1){
      output[i,2] <- "Non Reliable ID"
      output[i,3] <- "NA"
    }

    # Check others condtions
    # 2/more scores =/> 2.2, gen field not "No"
    else if (length(which(score > 2.19)) > 2 & length(unique(gen[gen != "No"])) == 1){

      if (!("sp." %in% sp)){

        if (length(unique(otu_id)) == 1){
          output[i,2] <- "Perfect ID"
          output[i,3] <- unique(otu_id)
        }

        else {
          output[i,2] <- "Taxonomic Group"
          output[i,3] <- paste(unique(gen), paste(unique(sp)[1], unique(sp)[2], sep="_"), sep=" ")
        }
      }

      else {
        output[i,2] <- "Novel Species"
        output[i,3] <- paste(unique(unique(gen[gen != "No"])), "sp.", sep=" ")
      }

    }

    # Check for genus match
    else if (length(which(score > 1.99)) > 2){   # 2/more scores =/> 2.0
      if (length(unique(gen)) == 1){ # unique gen
        output[i,2] <- "Genus Match"
        output[i,3] <- unique(gen)
      }
      else {
        output[i,2] <- "Check" # For anything that doesnt fit
        output[i,3] <- "NA"
      }

    }

  }

  # Set column names
  colnames(output) <- c("name", "result" ,"taxon")

  # Output final result file as text file
  write.table(output, file = result_name, sep="\t", row.names = FALSE, quote=FALSE)

}
