#' @import GO.db
#' @import GSEABase
#' @import org.Dm.eg.db
#' @import hash  
#' @export
#' @title Function to Count the Number of Genes with an Annotation to  
#' an Ontology Term for each Term within an Ontology.
#' @description This function counts the number of genes that are associated 
#' with each term 
#' of a specified ontology within a given geneset. Furthermore it counts the 
#' it counts the number of genes in the geneset that are characterized by the 
#' ontology, that is genes that have an annotation in the ontology.
#' @param  geneset geneset of interest e.g. up-regulated genes. Must be supplied
#'  as character vector of FbgnIDs.
#' @param  modus "GO" or "PANTHER". Defaults to "GO".
#' @param  cagetory GO category. Can be set to "BP", "MF", or "CC". Defaults to 
#' "BP".
#' @param  GO_slim  logical. If TRUE, the function only uses PANTHER GO slim 
#' terms. Defaults to FALSE.
#' 
#' 
#' @return \code{term_2_gene_count} returns a list of four elements.\cr 
#' The first element identified by "countTable" is a dataframe containing a 
#' column with 
#' the ontology identifers and the number of genes in the geneset that are 
#' annotated to it.\cr The second
#' element identified by "num_char" is an integer giving the number of 
#' characterized genes in the geneset.\cr
#' The third element identified by "num_all"
#' is an integer giving the number of all genes that are present in 
#' the geneset.\cr
#' The fourth element identified by "gen_annot_2_ID" is a dataframe with the 
#' Entrez IDs that are annotated to each Ontology ID.
#' 
#' @details This function could potentially be extended to use other ontologies,
#'  if a mapping between genes and terms of the desired ontology is 
#'  retrievable from some database. Please use up-to-date Fbgn IDs when using
#'  this function.
term_2_gene_count <- function(geneset, 
                            modus ="GO",
                            category ="BP",
                            GO_slim = FALSE
                            ){ 
PANTHER.db::pthOrganisms(PANTHER.db) <- "FLY" #  sets panther.db to fly
# gets panther go slim ids from obo file format

  
# functions to remove NULL objects from list of lists ==========================
is.NullOb <- function(x) is.null(x) | all(sapply(x, is.null))
      
rmNullObs <- function(x) {
  x <- Filter(Negate(is.NullOb), x)
  lapply(x, function(x) if (is.list(x)) rmNullObs(x) else x)
  }                                                         
#===============================================================================

# retrieves entrez ids that are required for mapping to goids    
EIDsgen <- unlist(AnnotationDbi::mget(geneset, org.Dm.egFLYBASE2EG, ifnotfound = NA)) 
    
if (modus == "GO"){
# character vecors of all GOIDs in respective categories; used to find 
#characterized genes
## "GO:0008150" is term "BP"
  BP_terms <- c(AnnotationDbi::as.list(GO.db::GOBPOFFSPRING)[["GO:0008150"]], "GO:0008150")  
  MF_terms <- c(AnnotationDbi::as.list(GO.db::GOMFOFFSPRING)[["GO:0003674"]], "GO:0003674")  
  CC_terms <- c(AnnotationDbi::as.list(GO.db::GOCCOFFSPRING)[["GO:0005575"]], "GO:0005575")  
      
#' map ezids to goids; no NAs in input allowed;  output is list of lists
#' with evidence source information and other information
  GOIDsgen_ann <- AnnotationDbi::mget(EIDsgen[!is.na(EIDsgen)], org.Dm.egGO, ifnotfound = NA)
   
      
#' make list of GO_IDs for every Entrez_ids -> list of lists
#' unique to remove replicates from different evidence sources
  GOIDsgen <- lapply(lapply(GOIDsgen_ann, names), unique) 
  
  # rmNUllObs prerequisite for use of relist()    
  GOIDsgen <- as.relistable(rmNullObs(GOIDsgen))
  
  # flatten list for intersection    
  GO_term_list <- unlist(GOIDsgen)  
      
  if (category == "BP"){
    
    #check if genes are characterized by ontology category
    only_BP_terms <- GO_term_list %in% BP_terms 
    
    
    # relist  --> entrez to go term, back to old structure    
    GO_term_list_BP <- relist(flesh = ifelse(only_BP_terms == TRUE, GO_term_list, NA),
                              skeleton = attr(GO_term_list, "skeleton")) 
    
    # remove nas
    GO_term_list_BP <- lapply(GO_term_list_BP, function(x) x[!is.na(x)])
    
    
    # removes empty lists, genes that are not characterized by "BP" are deleted
    GO_term_list_BP <- GO_term_list_BP[lapply(GO_term_list_BP,length)>0]  
    GOIDsgen <- GO_term_list_BP
    
    
    ancestor_hash <- hash(AnnotationDbi::as.list(GO.db::GOBPANCESTOR)) #  produce lookup table
    my_keys = hash::keys(ancestor_hash)
    helper.list <- lapply(GO_term_list_BP, function(x) ancestor_hash[x])
    helper.list = lapply(helper.list, AnnotationDbi::as.list)
    
    # puts list identifier in the list and removes replicate GOIDs
    GOIDsgen <- lapply(helper.list, function(x) unique(c(names(x), unlist(x)))) 
   
    
    }
      
  else if (category == "MF"){
    only_MF_terms <- GO_term_list %in% MF_terms
        
    GO_term_list_MF <-  relist(flesh = ifelse(only_MF_terms == TRUE, GO_term_list, NA),
                               skeleton= attr(GO_term_list, "skeleton")) # relist  --> entrez to go term
    GO_term_list_MF <- lapply(GO_term_list_MF, function(x) x[!is.na(x)])  # remove nas
    GO_term_list_MF <- GO_term_list_MF[lapply(GO_term_list_MF,length)>0]  # removes empty lists
    GOIDsgen <- GO_term_list_MF
    
    ancestor_hash <- hash(AnnotationDbi::as.list(GO.db::GOMFANCESTOR)) #  produce lookup table
    my_keys = hash::keys(ancestor_hash)
    helper.list <- lapply(GO_term_list_MF, function(x) ancestor_hash[x])
    helper.list = lapply(helper.list, AnnotationDbi::as.list)
    
    # puts list identifier in the list and removes replicate goids
    GOIDsgen <- lapply(helper.list, function (x) unique(c(names(x), unlist(x))))
    
    }
      
  else if (category == "CC"){
    only_CC_terms <- GO_term_list %in% CC_terms
        
    GO_term_list_CC <-  relist(flesh = ifelse(only_CC_terms==TRUE, GO_term_list, NA),
                               skeleton= attr(GO_term_list, "skeleton"))  # relist  --> entrez to go term
    GO_term_list_CC <- lapply(GO_term_list_CC, function(x) x[!is.na(x)])  # remove nas
    GO_term_list_CC <- GO_term_list_CC[lapply(GO_term_list_CC,length)>0]  # removes empty lists
    GOIDsgen <- GO_term_list_CC
    
    ancestor_hash <- hash(AnnotationDbi::as.list(GO.db::GOCCANCESTOR)) #  produce lookup table
    my_keys = hash::keys(ancestor_hash)
    helper.list <- lapply(GO_term_list_CC, function(x) ancestor_hash[x])
    helper.list = lapply(helper.list, AnnotationDbi::as.list)
    
    # puts list identifier in the list and removes replicate goids
    GOIDsgen <- lapply(helper.list, function (x) unique(c(names(x), unlist(x))))
  }
  
  chargenes <- length(GOIDsgen) # number of characterized genes in geneset
  GOIDsgen <- lapply(GOIDsgen, function (x) x[which(x != "all")] )  # removes annoying 
                                                           # "all" entries that
                                                           # cannot be passed to
                                                           # AnnotationDbi::mget, "all" 
                                                           # probably marks
                # genes that can be assigned to a gene ontology category but not to a specific term.
                                                                    
  
  # that function I found late it reverses the mapping from gene_to_ids to
  # to id_to_genes
  genes_4_ID <- Biobase::reverseSplit(GOIDsgen)
  
  # table of goids with respective number of genes that are annotated to it
  # here the counting happens; lapply(genes_4_ID, length) would give the same
  # counting
  countgen <- as.data.frame(table(factor(unlist(GOIDsgen)))) 
  names(countgen)<- c("ID","number_of_genes_with_annotation") 
  
    
  
#______________________Panther GOslim________________________________________  
  if (GO_slim == TRUE){
  intsec_index <- countgen$ID %in% Panther_GO_Slim     
  countgen <- countgen[intsec_index,]
  
  
  names(countgen)<- c("ID", "number_of_genes_with_annotation")
  
  #  number of genes characterized by GO_slim terms
  chargenes <- sum(unlist(lapply(GOIDsgen, function(x) any(x %in% Panther_GO_Slim))))
          
  }
  }
#_______________________________________________________________________________

  else if (modus=="PANTHER"){
    #  I m not sure about the structure of panther pathways; do they have a hierarchical structure like GO
    #  or are they independent? 
    
    # maps entrez ids with removed NAs to Panther Pathways and removes EZs with no annotation  
    Path_IDs <- na.omit(PANTHER.db::select(PANTHER.db, EIDsgen[!is.na(EIDsgen)],
                                           "PATHWAY_ID", "ENTREZ")) 
    # makes a list of  with Entrez_ID as identifier and content the assigned Pathway ids
    Path_IDs <- lapply(split(Path_IDs, Path_IDs$ENTREZ), function(x) x$PATHWAY_ID) 
    
    genes_4_ID <- Biobase::reverseSplit(Path_IDs)
    
    # table with pathway ids with respective number of gene annotations to it
    countgen <- as.data.frame(table(factor(unlist(Path_IDs))))  
    
    names(countgen) <- c("ID", "number_of_genes_with_annotation")
    
    
    chargenes <-length(Path_IDs) # this would only count the number
            #genes that are characterized by panther pathway which is very few
            #and therefore maybe to strict?, too pessimistic?
    
  
  }

  else if (modus == "KEGG"){
    keggPaths <-  AnnotationDbi::mget(EIDsgen[!is.na(EIDsgen)], org.Dm.egPATH, ifnotfound = NA)
    keggPaths <- keggPaths[!is.na(keggPaths)] #remove entrez ids that do not have an annotation for pathway
    genes_4_ID <- Biobase::reverseSplit(keggPaths) # list pathway IDs and the genes that are annotated to it
    
    countgen <- as.data.frame(table(factor(unlist(keggPaths))))  # counts numbe rof annotations for every ID
    names(countgen) <- c("ID", "number_of_genes_with_annotation") 
    
    chargenes <- length(keggPaths)
    
  }
    
# number of all genes charecterized and non characterized; this value is taken as number for 
# the gene ratio in ggplot
numbAll <- length(geneset)
  
# Return outputList: contains table with GO_ID and number 
parameter_counts <- list(countgen, chargenes, numbAll, genes_4_ID)
names(parameter_counts) <-  c("countTable", "num_char", "num_all", "gen_annot_2_ID")


return(parameter_counts)

} 




#' @export                
#' @title Testing for Enrichment and Depletion of Ontology Terms in a Geneset
#' @description \code{term_2_gene_count}  tests for enrichment or depletion of
#' ontology term as compared to randomly sample a set of genes in the reference
#' set of \code{length(your_favourite_experimental_set)}.     
#' @param parameters_list_exp the list that is created by 
#' \code{term_2_gene_count()} for your \emph{gene set of interest} e.g. 
#' down-regulated genes
#' @param parameters_list_ref the list that is created by 
#' \code{term_2_gene_count()} for your \emph{reference gene set}.
#' @param stage specify name of experiment.
#' @param modus "GO", "KEGG" or "PANTHER" analysis. Defaults to "GO".
#' @param twosidetest "midpval" or "simpledouble" according to PMID: 17182697.
#' Default is "midpval".
#' @param multCorrect "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#' "fdr", "none". Default is "bonferroni".

#' @return \code{term_enrichment} creates a dataframe with information of 
#' the statistical analysis.\cr
#' The column named "k" represents the number of genes
#' in the experimental set with an annotation to the row ID.\cr 
#' The column named "m" represents the number 
#' of genes in the reference set with an annotation to the row ID.\cr
#' The columns "kcorrected", "whichmin", "pvl_R" and "pval_L" are essential
#' for the twosided test. Please see PMID: 17182697 and source code.\cr
#' The column "k_to_n" ratio gives the proportion of genes in the experimental
#' set that have an annotation to the ID specified in the row.\cr
#' The column "twosided.logpval" is the negative log10 of the twosided pvalue.
#' \cr
#' "E" stands for enriched and "D" for depleted. 
#' @details Please be careful when you set the background or reference set.
#' The background has to be set in a reasonable fashion to get an appropriate
#' answer to a predefined biological question.  
#' 
term_enrichment <- function(parameters_list_exp, # contains parameters for experimental set
                            parameters_list_ref, # contains parameers for reference set
                            stage,    # number of experimental stage
                            modus = "GO", # "GO" "KEPANTHER" analysis
                            twosidetest = "midpval",
                            multCorrect = "bonferroni"){ # "twoside "midpval" or "simpledouble"
                                                    # according to PMID: 17182697
  
  PANTHER.db::pthOrganisms(PANTHER.db) <- "FLY" #  sets panther.db to fly  
  N <- parameters_list_ref$num_char # num of characterized genes in refset   
  n <- parameters_list_exp$num_char # num of chararacterized (with annotation) genes in expset
  
  num_exp <- stage # number of experimental stage, is added to the table as column
  
  
  tabks <- merge(parameters_list_exp$countTable, parameters_list_ref$countTable, by ="ID", all = TRUE)  
  tabks$ID <-  as.character(tabks$ID)  # to make the GO_ID column a character  vector from factor type for sorting
  tabks[is.na(tabks)] <- 0            ## replaces nas in dataframe (GO ids from experimental and reference data set that have no intersection)
  tabks <- tabks[order(tabks[1]), ] # sorts dataframe by GO ID
  
  
  
  
  
  ls_pval_enriched <- vector("numeric", length = dim(tabks)[1]) #  right tail test                    
  names(ls_pval_enriched) <- vector("character", length = dim(tabks)[1]) 
    
  ls_pval_depleted <- vector("numeric", length = dim(tabks)[1])# left tail test                    
  names(ls_pval_depleted) <- vector("character", length = dim(tabks)[1])
  
  
  #### rename colums
  names(tabks)[c(2,3)] <- c("k","m")
  
  ### assign gene ontology category("BP", "CC", "MF") to each GO_ID in tabk and add col #########
  if(modus=="GO"){
  tabks$GO_categories <- sapply(AnnotationDbi::mget(tabks$ID, GO.db::GOTERM, ifnotfound = NA),
                                function(x) x@Ontology)
  
  }
  
  
  # add extra col for correction when plug in phyper (issue with margin value included or excluded  P[X ≤ x]  P[X > x]) 
  tabks$kcorrected <- ifelse(tabks$k== 0 , 0 , tabks$k-1) 
  
  numtests <- dim(tabks)[1] # for Bonferroni correction
  index <- 0 #  initializing value
  
  for (i in tabks$ID){ #=============================================================================================================
    # statistical heart of the function; for every GOID the twosided pvalue is
    # calculated from geometric distribution using folowing parameters
    # k = #genes in expset with annotation to the specified GOID
    # m = #genes in refset with annotation to the specified GOID
    # N = #allgenes in reference set that are characterized meaning they have some annotation to a GOID
    # n = #allgenes in expset that are characterized
    index <- index + 1
    names(ls_pval_enriched)[index] <- i 
    names(ls_pval_depleted)[index] <- i    
    
    if(twosidetest=="simpledouble"){
   #gives back named list (names are corresponding GO_IDs) of p values for
   # 
   # 
      if (tabks[tabks$ID == i, "kcorrected"] == 0 & (tabks[tabks$ID == i, "k"] == 0)) {
        ls_pval_enriched[index] <-  phyper(tabks[tabks$ID == i, "kcorrected"], #  k
                                       tabks[tabks$ID == i,"m"], #  m
                                       (N - tabks[tabks$ID == i, "m"]), # N-m
                                       n, lower.tail = FALSE)  # logical; if TRUE (default), 
                                                               # probabilities are 
                                                               # P[X ≤ x], otherwise, P[X > x]
                                                               
        + dhyper(tabks[tabks$ID == i, "k"],
                 tabks[tabks$ID == i,"m"],
                (N - tabks[tabks$ID == i, "m"]),
                 n) 
      }
      else ls_pval_enriched[index] <-  phyper(tabks[tabks$ID == i, "kcorrected"],
                                              tabks[tabks$ID == i,"m"],
                                              (N - tabks[tabks$ID == i, "m"]),
                                              n, lower.tail = FALSE) 
    
        ls_pval_depleted[index] <- phyper(tabks[tabks$ID == i, "k"],
                                      tabks[tabks$ID == i,"m"],
                                      (N - tabks[tabks$ID == i, "m"]),
                                       n)          
   }
    
    else if (twosidetest=="midpval"){
    ls_pval_enriched[index] <-  phyper(tabks[tabks$ID == i, "k"],
                                       tabks[tabks$ID == i,"m"],
                                       (N - tabks[tabks$ID == i, "m"]),
                                       n, lower.tail = FALSE) +
     (0.5 * dhyper(tabks[tabks$ID == i, "k"],
                    tabks[tabks$ID == i,"m"],
                    (N - tabks[tabks$ID == i, "m"]),
                    n))
    
      if ((tabks[tabks$ID == i, "kcorrected"] == 0) & (tabks[tabks$ID == i, "k"] == 0)){
        ls_pval_depleted[index] <- (0.5 * dhyper(tabks[tabks$ID == i, "k"],
                                                 tabks[tabks$ID == i,"m"],
                                                 (N - tabks[tabks$ID == i, "m"]),
                                                  n))
      }
    
      else ls_pval_depleted[index] <- phyper(tabks[tabks$ID == i, "kcorrected"],
                                           tabks[tabks$ID == i,"m"],
                                           (N - tabks[tabks$ID == i, "m"]),
                                            n)  +
      (0.5 * dhyper(tabks[tabks$ID == i, "k"],
                    tabks[tabks$ID == i,"m"],
                    (N - tabks[tabks$ID == i, "m"]),
                     n))
    }
    
  }
  #========================================================================================================
  
  ##### add cols for p vals enriched and depleted
  tabks$pvl_R <- ls_pval_enriched ## right tail P(X>=k)
  tabks$pvl_L <- ls_pval_depleted ## left tail  P(X<=k)
   
  ##### add col for two-sided test
  tabks$whichmin <- apply(tabks[,c("pvl_R","pvl_L")], 1, which.min)
  
  #  p_doublingtwo (k) = 2 * min[P(X>=k),P(X<=k)] PMID: 17182697
  tabks$twosidepval <-ifelse(tabks$whichmin ==1, 2*tabks$pvl_R, 2* tabks$pvl_L)
  
  # add col for k to n(all genes characterized and not characterized) ratio
  tabks$k_n_ratio <- unlist(lapply(tabks$k, function(x) x/parameters_list_exp$num_all))  # unlist to make it numeric vector
  
  if (modus == "GO"){
    
       # adds column with the adjusted pvalues
       tabks$padj <- p.adjust(tabks$twosidepval, method = multCorrect, n = numtests)
       
       
       tabks$TERM <- as.character(sapply(AnnotationDbi::mget(tabks$ID, GO.db::GOTERM,
                                              ifnotfound = NA),
                                              function(x) x@Term))
       
       # marks enriched or depleted GOIDs
       tabks$E_or_D <- factor(ifelse(tabks[ ,"whichmin"] == 1, "E","D")) 
       
       tabks$twosided.logpval <-  -log10(tabks[ ,"twosidepval"])
       
       # needed for ggplot scale_gradient
       tabks$twosided.logpval <- ifelse(tabks$E_or_D=="E",
                                    tabks$twosided.logpval,
                                    -1* tabks$twosided.logpval)
       
       tabks$log.padj <- -log10(tabks$padj)
       tabks$log.padj <- ifelse(tabks$E_or_D == "E",
                                tabks$log.padj,
                                -1* tabks$log.padj)
       
       tabks$state <-  replicate(dim(tabks)[1], num_exp)
       tabks <- rename(tabks, c(k_n_ratio="k_to_n_ratio")) 
       tabks$kdivn <-  base::paste(tabks[ ,"k"], "/", parameters_list_exp$num_all, sep="")
       
   }
  
  else if (modus =="PANTHER"){
    
  temp.tabks <- select(PANTHER.db, tabks[ ,"ID"],"PATHWAY_TERM", "PATHWAY_ID")

  temp.tabks <-  temp.tabks[order(temp.tabks["PATHWAY_ID"]),] # order by pathway id
  
  temp.tabks$padj <- p.adjust(tabks$twosidepval, method = multCorrect, n = numtests)
  
 
  temp.tabks$E_or_D <- factor(ifelse(tabks[, "whichmin"] == 1, "E","D"))
  temp.tabks$twosided.logpval <- -log10(tabks[, "twosidepval"])
  temp.tabks$twosided.logpval <- ifelse(temp.tabks$E_or_D == "E",
                                    temp.tabks$twosided.logpval,
                                    -1* temp.tabks$twosided.logpval)
  
  temp.tabks$log.padj <- -log10(temp.tabks$padj) 
  temp.tabks$log.padj <- ifelse(temp.tabks$E_or_D == "E",
                                temp.tabks$log.padj,
                                -1* temp.tabks$log.padj)
  
   
  temp.tabks$k_to_n_ratio <- tabks[, "k_n_ratio"]
  temp.tabks$kdivn        <- base::paste(tabks[, "k"],"/", parameters_list_exp$num_all, sep ="")
  names(temp.tabks)[names(temp.tabks)=="PATHWAY_TERM"]   <-  "TERM"
  names(temp.tabks)[names(temp.tabks)=="PATHWAY_ID"] <- "ID"
  temp.tabks$state        <- replicate(dim(temp.tabks)[1], num_exp)
  tabks                   <- temp.tabks
    
  }
  
  else if (modus == "KEGG"){
  tabks$padj <- p.adjust(tabks$twosidepval, method = multCorrect, n = numtests)
    
  temp.ids  <-  sapply(tabks$ID, function(x) paste0('map', x))
  
  # split list because keggGet can take a vector of length 10
  temp.ids  <-  base::split(temp.ids, ceiling(seq_along(unlist(tabks$ID))/10)) 

  
  
  res <-  lapply(temp.ids, KEGGREST::keggGet) # gets a list of the kegg input
  res2 <-  unlist(res, recursive = FALSE) # l
  names(res2) <-  sapply(res2, function(x) x$ENTRY)
  terms <-  lapply(res2, function(x) x$NAME)
  
  try(if(identical(names(res2), unname(unlist(temp.ids))) == FALSE) stop("there are gaps in the mapping between ID and term name, IDs originated from org.Dm.eg.db"))
  tabks$TERM <- base::unname(base::unlist(terms))                                        
                                      
  # marks enriched or depleted GOIDs
  tabks$E_or_D <- factor(ifelse(tabks[ ,"whichmin"] == 1, "E","D")) 
  
  tabks$log.padj <- -log10(tabks$padj)
  
  tabks$log.padj <- ifelse(tabks$E_or_D == "E",
                           tabks$log.padj,
                           -1* tabks$log.padj) 
    
 
    
  tabks$twosided.logpval <-  -log10(tabks[ ,"twosidepval"])
    
  # needed for ggplot scale_gradient
  tabks$twosided.logpval <- ifelse(tabks$E_or_D=="E",
                                     tabks$twosided.logpval,
                                     -1* tabks$twosided.logpval)
    
  tabks$state <-  replicate(dim(tabks)[1], num_exp)
  tabks <- rename(tabks, c(k_n_ratio="k_to_n_ratio")) 
  tabks$kdivn <-  base::paste(tabks[ ,"k"], "/", parameters_list_exp$num_all, sep="")
  }
  
  
  return(tabks)
  
}




#'@export
#'@title Prune the GO Graph Beginning from the Roots
#'@description The \code{cut_GO()} function takes a character vector of GOIDS as
#'input and removes the terms to a specified level in the DAG.
#'@param GOids character vector of GOIDs
#'@param level up to which level the terms should be deleted ? Defaults to 2.
#'@param startingnodes set the starting nodes. The specified level is relative 
#'to the starting nodes. The default is set to the ID representing the terms:
#''BP', 'MF', 'CC'.
#'@return The output is a truncated character vector of the entered GOIDs. 
#'The more general terms are removed.

cut_GO <- function(GOids, level = 2, startingnodes=c("GO:0008150","GO:0003674","GO:0005575")){         
  
  # empty list for init.
  whole <- vector(mode = "list")
  
  for (i in seq_len(level-1)){
    whole[[i]] <- startingnodes
    
    startingnodes <- c(unlist(AnnotationDbi::mget(startingnodes, GO.db::GOBPCHILDREN, ifnotfound = NA)), unlist(AnnotationDbi::mget(startingnodes, GO.db::GOMFCHILDREN, ifnotfound = NA)),  unlist(AnnotationDbi::mget(startingnodes, GO.db::GOCCCHILDREN, ifnotfound = NA)))[!is.na(c(unlist(AnnotationDbi::mget(startingnodes, GO.db::GOBPCHILDREN, ifnotfound = NA)), unlist(AnnotationDbi::mget(startingnodes, GO.db::GOMFCHILDREN, ifnotfound = NA)),  unlist(AnnotationDbi::mget(startingnodes, GO.db::GOCCCHILDREN, ifnotfound = NA))))]
    
    
  }
  whole <- unlist(whole)  
  cutgo_index <- GOids %in% whole
  pruned_IDs <- GOids[!cutgo_index]
  return(pruned_IDs)
}


                              
#' @export
#' @title  Merge the dataframes of different stages and cut out general GO 
#' terms
#' @description \code{merge_n_clean()} takes the the data output of 
#' \code{term_enrichment()} for different stages, merges them, cuts out more
#' general GO terms and subsets the merged tables by the GO terms that are 
#' significantly enriched or depleted in any of the entered stages.
#' @param diff_stage_tables put in dataframes obtained by \code{term_enrichment} for several
#' stages as list
#' @param cut_level  up to which level starting from the root the terms should 
#' be deleted ? Defaults to 2.
#' @param sort_by how should the terms be sorted for the ggplot script ? The
#' default is "pvalSum" and it sorts the terms according to the sum of p-values
#' over the different stages. The second option is "variance" which sorts the 
#' terms according to the biggest spread in p-values over the stages.
#' @param mode "GO", "PANTHER" or "KEGG"
#' @return The output is a list of two items.\cr
#' The first item identified by "table_all_stages_sign_IDs" is a dataframe which
#' contains all Ontology IDs that are significantly enriched or depleted in any
#' of the experimental stages. This dataframe is the input for the ggplot script
#' \cr
#' The second item identified by "order_by_pval" specifies the order of terms 
#' for the visualization in the ggplot script.
merge_n_clean <- function(diff_stage_tables,
                          cut_level = 2, # cutting general GO terms top down
                          pval  =  0.01, # desired pval threshold
                          sort_by = "pvalSum",
                          mode = "GO"){
  
  PANTHER.db::pthOrganisms(PANTHER.db) <- "FLY"  # sets panther database to fly
  
  mergedTabs <- diff_stage_tables[[1]]
  for (i in 2:length(diff_stage_tables)) mergedTabs <- rbind(mergedTabs, diff_stage_tables[[i]])
  
  if(!any(mergedTabs$padj < pval)) return(NULL) #FIXME extremely ugly solution to prevent going on if there are no significant terms.

  sig_ID_all_states <- unique(mergedTabs[mergedTabs$padj < pval, "ID"])  # significant go terms from all different states ; --> cut out general terms with the cut_go function
  
  if(mode == "GO"){
  sig_ID_all_states <- cut_GO(sig_ID_all_states, level = cut_level) # trims away super general terms like "biological process"
  }
  
  tabs <- mergedTabs[mergedTabs$ID %in% sig_ID_all_states, ]
  
  if(sort_by == "pvalSum"){
  
  sumpval <- vector("numeric")                     
  names(sumpval) <- vector("character") 
  
  initz = 0
  for (i in sig_ID_all_states) {
    
    initz = initz + 1
    
    sumpval[initz] <- sum(tabs[tabs$ID == i, "twosided.logpval"])
    names(sumpval)[initz] <- i # assign ID to pval
    
  }
  
  sumpval <- base::sort(sumpval, decreasing =TRUE)
  
  if(mode == "GO"){
  term_order <- unname(unlist(lapply(AnnotationDbi::mget(names(sumpval), GO.db::GOTERM), function(x) x@Term)))   # order of terms for ggplot; defined in scale_x_discrete
  }
  
  else if (mode == "PANTHER"){
  term_order <- PANTHER.db::select(PANTHER.db::PANTHER.db, keys = names(sumpval), columns = "PATHWAY_TERM", keytype = "PATHWAY_ID")[,2]
  }
  
  else if (mode == "KEGG"){
    
    ordered_ids <- base::split(names(sumpval), ceiling(seq_along(names(sumpval))/10)) # split ids for KEGGREST::keggGet
    
    ordered_ids <- sapply(ordered_ids, function(x) paste0('map', x))
    
    res <-  lapply(ordered_ids, KEGGREST::keggGet) # gets a list of the kegg input
    res2 <-  unlist(res, recursive = FALSE) # l
    names(res2) <-  sapply(res2, function(x) x$ENTRY)
    terms <-  lapply(res2, function(x) x$NAME)

    term_order <- base::unname(base::unlist(terms)) 
    
  }
  
  }
  
  else if(sort_by == "variance"){
    
  variance <- vector("numeric")                     
  names(variance) <- vector("character") 
    
    initz = 0
    for (i in sig_ID_all_states) {
      
      initz = initz + 1
      
      variance[initz] <- var(tabs[tabs$ID == i, "twosided.logpval"])
      names(variance)[initz] <- i
      
    }
    
    
    variance <- base::sort(variance, decreasing =TRUE)
    
    if(mode == "GO"){
    term_order <- unname(unlist(lapply(AnnotationDbi::mget(names(variance), GO.db::GOTERM), function(x) x@Term)))   # order of terms for ggplot; defined in scale_x_discrete
    }
    
    else if (mode == "PANTHER"){
    term_order <- PANTHER.db::select(PANTHER.db::PANTHER.db, keys = names(variance), columns = "PATHWAY_TERM", keytype = "PATHWAY_ID")[,2]
    
    }
    
    else if (mode == "KEGG"){
      ordered_ids <- base::split(names(variance), ceiling(seq_along(names(variance))/10)) # split ids for keggget
      ordered_ids <- sapply(ordered_ids, function(x) paste0('map', x))
      
      res <-  lapply(ordered_ids, KEGGREST::keggGet) # gets a list of the kegg input
      res2 <-  unlist(res, recursive = FALSE) # l
      names(res2) <-  sapply(res2, function(x) x$ENTRY)
      terms <-  lapply(res2, function(x) x$NAME)
      term_order <- base::unname(base::unlist(terms)) 
    }
    
  }
  
  output_list <- list(tabs, rev(term_order))
  names(output_list) <-  c("table_all_stages_sign_IDs", "order_by_pval")
  return(output_list)
  
}
 
  




  








