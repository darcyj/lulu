#' Post Clustering Curation of Amplicon Data.
#' @description This algorithm \code{lulu} consumes an OTU table and a matchlist, and
#'   evaluates cooccurence of 'daughters' (potential analytical artefacts) and
#'   their 'parents' (≈ real biological species/OTUs). The algorithm requires an
#'   OTU table (species/site matrix), and a match list. The OTU table can be
#'   made with various r-packages (e.g. \code{DADA2}) or external pipelines (\code{VSEARCH, USEARCH, QIIME}, etc.), and the
#'   match-list can be made with external bioinformatic tools like \code{VSEARCH, USEARCH, BLASTN} or another algorithm
#'   for pair-wise sequence matching.
#' @param otutable a data.frame with with an OTU table that has sites/samples as
#'   columns and OTUs (unique OTU id's) as rows, and observations as read
#'   counts.
#' @param matchlist a data.frame containing three columns: (1) OTU id of
#'   potential child, (2) OTU id of potential parent, (3) match - \% identiti
#'   between the sequences of the potential parent and potential child OTUs.
#'   \strong{NB: The matchlist is the product of a mapping of OTU sequences against each other. This is
#'   currently carried out by an external script in e.g. Blastn or VSEARCH, prior to running lulu!}
#' @param minimum_ratio_type sets whether a potential error must have lower
#'   abundance than the parent in all samples \code{min} (default), or if an error
#'   just needs to have lower abundance on average \code{avg}. Choosing lower
#'   abundance on average over globally lower abundance will greatly increase
#'   the number of designated errors. This option was introduced to make it
#'   possible to account for non-sufficiently clustered intraspecific variation,
#'   but is not generally recommended, as it will also increase the potential of
#'   cluster well-separated, but co-occuring, sequence similar species.
#' @param minimum_ratio sets the minimim abundance ratio between a potential error
#'   and a potential parent to be identified as an error. If the \code{minimum_ratio_type} is
#'   set to \code{min} (default), the \code{minimum_ratio} applies to the lowest observed
#'   ration across the samples.  If the \code{minimum_ratio_type} is
#'   set to \code{avg} (default), the \code{minimum_ratio} applies to the mean of observed
#'   ration across the samples.\code{avg}. (default is 1).
#' @param minimum_match minimum threshold of sequence similarity
#'   for considering any OTU as an error of another can be set (default 84\%).
#' @param minimum_relative_cooccurence minimum co-occurrence rate – i.e. the
#'   lower rate of occurrence of the potential error explained by co-occurrence
#'   with the potential parent for considering error state.
#' @param log_name name for log file (default is "auto"). If "auto", a name will
#'   be generated starting with "lulu.log_" and ending in the system date and time.
#'   If "none", no log will be written. Any other value will cause a log to be written
#'   under that name. 
#' @param output_type the type of output desired. If "otutable", default list output
#'  will be overridden, and only an otutable (class data.frame) will be returned.
#' @return Function \code{lulu} returns a list of results based on the input OTU
#'   table and match list.
#'   \enumerate{
#'   \item \code{curated_table} - a curated
#'   OTU table with daughters merged with their matching parents.
#'   \item \code{curated_count} - number of curated (parent) OTUs.
#'   \item \code{curated_otus} - ids of the OTUs that were accepted as valid OTUs.
#'   \item \code{discarded_count} - number of discarded (merged with parent) OTUs.
#'   \item \code{discarded_otus} - ids of the OTUs that were identified as
#'   errors (daughters) and merged with respective parents.
#'   \item \code{runtime} - time used by the script.
#'   \item \code{minimum_match} - the id threshold
#'   (minimum match \% between parent and daughter) for evaluating co-occurence (set
#'   by user).
#'   \item \code{minimum_relative_cooccurence} - minimum ratio of
#'   daughter-occurences explained by co-occurence with parent (set by user).
#'   \item \code{otu_map} - information of which daughters were mapped to which
#'   parents.
#'   \item \code{original_table} - original OTU table. }
#' @examples
#' lulu(my_table, my_matchlist)
#' lulu(my_table, my_matchlist, minimum_ratio_type = "avg", minimum_match = 0.8, minimum_ratio = 10, minimum_relative_cooccurence = 0.9)
#' @details The matchlist is the product of a mapping of OTU sequences against each other. This is
#'   currently carried out by an external script in e.g. BLASTN or VSEARCH, prior to running \code{lulu}!
#'   Producing the match list requires a file with all the OTU sequences (centroids) - e.g. \code{OTUcentroids.fasta}. The matchlist can be produced by mapping all OTUs against each other with an external algorithm like VSEARCH or BLASTN. In \code{VSEARCH} a matchlist can be produced e.g. with the following command: \code{vsearch --usearch_global OTUcentroids.fasta --db OTUcentroids.fasta --strand plus --self --id .80 --iddef 1 --userout matchlist.txt --userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10}. In \code{BLASTN} a matchlist can be produces e.g. with the following commands. First we produce a blast-database from the fasta file: \code{makeblastdb -in OTUcentroids.fasta -parse_seqids -dbtype nucl}, then we match the centroids against that database: \code{blastn -db OTUcentoids.fasta -num_threads 10 -outfmt'6 qseqid sseqid pident' -out matchlist.txt -qcov_hsp_perc .90 -perc_identity .84 -query OTUcentroids.fasta}
#' @author Tobias Guldberg Frøslev
#' @export
lulu <- function (otutable, matchlist, minimum_ratio_type = "min", minimum_ratio = 1, minimum_match = 84, minimum_relative_cooccurence = 0.95, log_name="auto", output_type="full"){
  require(dplyr)
  start.time <- Sys.time()
  colnames(matchlist) <- c("OTUid", "hit", "match")
  # remove no hits (vsearch)
  matchlist = matchlist[which(matchlist$hit != "*"), ]
  # remove self-hits
  matchlist = matchlist[which(matchlist$hit != matchlist$OTUid), ]

  # Making a separate table with stats (total readcount and spread).
  statistics_table <- otutable[, 0]
  statistics_table$total <- rowSums(otutable)

  # calculating spread (number of presences (samples with 1+ read) pr OTU)
  statistics_table$spread <- rowSums(otutable > 0)
  statistics_table <- statistics_table[with(statistics_table, order(spread, total, decreasing = TRUE)), ]
  otutable <- otutable[match(row.names(statistics_table), row.names(otutable)), ]
  # add parent_id column to statistics_table, it will store result for simplification
  
  # open log for writing, use new argument I added for log name. Default ("auto") preserves original functionality.
  logTF <- TRUE
  if(log_name == "auto"){
    log_fp <- paste0("lulu.log_", format(start.time, "%Y%m%d_%H%M%S"))
  }else if(log_name=="none"){
    logTF <- FALSE
  }else{
    log_fp <- log_name
  }
  if(logTF == TRUE){ log_con <- file(log_fp, open="a") }
  
  # define function for writing to log. This tidies up the code a bit, and enables user to turn off log.
  write2log <- function(text, writeTF=logTF, con=log_con){
    if(writeTF == TRUE){
      cat(text, file=con)
    }
  }
  
  # start up progress bar
  pb <- txtProgressBar(min=0, max=nrow(statistics_table), style=3)

  # for each otu (in statistics_table / otu_table), 
  # there was an extra "seq" in this line, it's not necessary, I removed it.
  for (line in 1:nrow(statistics_table)) {
      # progress bar will be updated at the end of the loop.
      potential_parent_id <- row.names(otutable)[line]
      write2log(paste0("\n", "#### Processing: ", potential_parent_id, " #####"))
      daughter_samples <- otutable[line, ]
      hits <- matchlist[which(matchlist$OTUid == potential_parent_id &  matchlist$match > minimum_match), "hit"]
      write2log(paste0("\n", "---hits: ", hits))
      last_relevant_entry <- sum(statistics_table$spread >= statistics_table$spread[line])
      potential_parents <- which(row.names(otutable)[1:last_relevant_entry] %in% hits)
      write2log(paste0("\n", "---potential parent: ", row.names(statistics_table)[potential_parents]))
      
      success <- FALSE
      if (length(potential_parents) > 0) {
          for (line2 in potential_parents) {
              write2log(paste0("\n", "------checking: ", row.names(statistics_table)[line2]))
              if (!success) {
                relative_cooccurence <- sum((daughter_samples[otutable[line2, ] > 0]) > 0)/sum(daughter_samples > 0)
                write2log(paste0("\n", "------relative cooccurence: ", relative_cooccurence))
                if (relative_cooccurence >= minimum_relative_cooccurence) {
                  write2log(paste0(" which is sufficient!"))
                  if (minimum_ratio_type == "avg") {
                    relative_abundance <- mean(otutable[line2, ][daughter_samples > 0]/daughter_samples[daughter_samples > 0])
                    write2log(paste0("\n", "------mean avg abundance: ", relative_abundance))
                  } else {
                    relative_abundance <- min(otutable[line2, ][daughter_samples > 0]/daughter_samples[daughter_samples > 0])
                    write2log(paste0("\n", "------min avg abundance: ", relative_abundance))
                  }
                  if (relative_abundance > minimum_ratio) {
                    write2log(paste0(" which is OK!"))
                    if (line2 < line) {
                      statistics_table$parent_id[line] <- statistics_table[row.names(otutable)[line2], "parent_id"]
                      write2log(paste0("\n", "SETTING ", potential_parent_id, " to be an ERROR of ", (statistics_table[row.names(otutable)[line2], "parent_id"]), "\n"))
                    } else {
                      statistics_table$parent_id[line] <- row.names(otutable)[line2]
                      write2log(paste0("\n", "SETTING ", potential_parent_id, " to be an ERROR of ", (row.names(otutable)[line2]), "\n"))
                    }
                    success <- TRUE
                  }
                }
              }
          }
      }
      if (!success) {
          statistics_table$parent_id[line] <- row.names(statistics_table)[line]
          write2log(paste0("\n", "No parent found!", "\n"))
      }
      
      # update progress bar
      setTxtProgressBar(pb, line)
  }


  total_abundances <- rowSums(otutable)
  curation_table <- cbind(nOTUid = statistics_table$parent_id, otutable)
  statistics_table$curated <- "merged"
  curate_index <- row.names(statistics_table) == statistics_table$parent_id
  statistics_table$curated[curate_index] <- "parent"
  statistics_table <- transform(statistics_table, rank = ave(total, FUN = function(x) rank(-x, ties.method = "first")))

  # changed the following line to avoid deprecation error message
  # curation_table2 <- as.data.frame(curation_table %>% group_by(nOTUid) %>% summarise_each(funs(sum)))
  curation_table <- as.data.frame(curation_table %>% group_by(nOTUid) %>% summarise_all(funs(sum)))
  # quick check to see if they actually do the same thing, made the summarize_each one "curation_table_old"
  # all(curation_table_old == curation_table)
  # it's TRUE, we're good here
                                                             
  row.names(curation_table) <- as.character(curation_table$nOTUid)
  
  # remove redundant name column from curation_table
  curation_table <- curation_table[, -1]
  
  curated_otus <- names(table(statistics_table$parent_id))
  curated_count <- length(curated_otus)
  discarded_otus <- setdiff(row.names(statistics_table), curated_otus)
  discarded_count <- length(discarded_otus)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  
  # I made it so the result can be the bare-minimum output, the curated OTU table. 
  # This is useful if your output is HUGE - it speeds things up a bit.
  if(output_type == "otutable"){
      result <- curation_table
  }else{
      result <- list(curated_table = curation_table, curated_count = curated_count, 
          curated_otus = curated_otus, discarded_count = discarded_count, 
          discarded_otus = discarded_otus, runtime = time.taken, 
          minimum_match = minimum_match, minimum_relative_cooccurence = minimum_relative_cooccurence, 
          otu_map = statistics_table, original_table = otutable)
  }
  
  # write runtime to log, then close it
  time_in_mins <- round(as.numeric(time.taken, units="mins"), 4)
  write2log(paste0("\n", "#### FINISHED; took ", time_in_mins, " minutes."))
  close(log_con)

  return(result)
}


