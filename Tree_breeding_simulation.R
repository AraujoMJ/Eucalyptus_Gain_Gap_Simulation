#' @title Simulate Breeding Population with Pedigree, Genomic, and Haplotypic Analyses
#' 
#' @description 
#' This function simulates a breeding population with options for open-pollination (OP) or controlled-pollination (CP).
#' It generates a pedigree, simulates genotypic and phenotypic data, estimates variance components, heritabilities, and breeding values (BLUPs).
#' The function also calculates haplotype similarities to classify self-half-sibs (SHS) and self-full-sibs (SFS).
#' 
#' @param generations Integer. Number of generations to simulate (default = 2).
#' @param ids Integer or vector. Number of individuals per generation (default = 1000).
#' @param animals Logical. Whether the population is animal-based (default = FALSE).
#' @param nP Numeric. Selfing rate (proportion of inbreeding in F1 generation, default = 0.1).
#' @param nFam_F1 Integer. Number of families in F1 generation (default = 100).
#' @param num_cores Integer. Number of cores for parallel computation (default = NULL).
#' @param use_cores Logical. Whether to enable parallel computation (default = NULL).
#' @param n_Fam Integer. Number of families in the base generation.
#' @param n_Block Integer. Number of replication blocks for experiments.
#' @param n_Tree Integer. Number of trees per plot.
#' @param VarPlot Numeric. Standard deviation for plot effect (default = NULL).
#' @param nSNPPerChr Integer. Number of SNPs per chromosome.
#' @param nQTLPerChr Integer. Number of quantitative trait loci (QTL).
#' @param mu Numeric. Mean additive effect.
#' @param Vp Numeric. Phenotypic variance.
#' @param Heritability Numeric. Narrow-sense heritability.
#' @param Va Numeric. Additive genetic variance.
#' @param AddMat Character. Method for additive genomic relationship matrix (e.g., "VanRaden").
#' @param MeanDom Numeric. Mean dominance variance.
#' @param DomMat Character. Method for dominance genomic relationship matrix.
#' @param Vd Numeric. Dominance variance.
#' @param nProg Integer. Number of progenies per cross.
#' @param TypePop Character. Population type: "OP" (open-pollination) or "CP" (controlled-pollination).
#' @param Samp Numeric. Proportion of genotyped individuals for genomic analyses (default = NULL).
#' @param Text_track Character. Text description for tracking simulations.
#' @param pullSNPs Logical. Whether to pull SNPs (default = NULL).
#' @param pullHaplotype Logical. Whether to extract haplotypes (default = TRUE).
#' @param Get_only_data Logical. Whether to stop after generating data (default = TRUE).
#' 
#' @return 
#' A list containing:
#' \item{TypePop}{Character. Population type ("OP" or "CP").}
#' \item{SelfRate}{Character. Observed selfing rate in the F1 generation.}
#' \item{RealParameters}{Data frame. Simulated genetic parameters.}
#' \item{VarComp}{List. Variance components for different models.}
#' \item{Herit.All}{Data frame. Heritability estimates (narrow-sense, dominance, broad-sense).}
#' \item{BV}{List. Breeding values (BLUPs) from pedigree and genomic models.}
#' \item{BV_ALL}{Data frame. Breeding and genotypic values from all models.}
#' \item{VarComp.ssGBLUP}{List. Variance components from ssGBLUP models.}
#' \item{BV.ssGBLUP}{List. Breeding values from ssGBLUP models.}
#' \item{BV_ALL.ssGBLUP}{Data frame. Combined BLUPs and dominance effects from ssGBLUP.}
#' \item{Mean.Inbreeding}{Numeric. Mean inbreeding coefficient.}
#' \item{Data}{Data frame. Experimental design and phenotypic values.}
#' \item{DiagG}{Data frame. Diagonal elements (1 + Fi) from the G matrix.}
#' \item{snps}{Matrix. SNP genotypes (if \code{pullSNPs = TRUE}).}
#' 
#' @details 
#' The function performs the following steps:
#' 1. Generates a pedigree for specified generations.
#' 2. Simulates phenotypes, genotypes, and genetic effects for base and progeny populations.
#' 3. Estimates variance components, heritabilities, and breeding values (BLUPs) using ASReml-R models.
#' 4. Performs genomic-based (G and GD matrices) and single-step GBLUP (ssGBLUP) analyses.
#' 5. Calculates haplotype similarities for selfing classification (SHS vs. SFS).
#' 
#' @examples
#' \dontrun{
#' # Simulate a CP population with 1000 individuals, 10% selfing rate
#' result <- TB_SimulR(
#'   generations = 2,
#'   ids = 1000,
#'   nP = 0.1,
#'   nFam_F1 = 100,
#'   n_Block = 3,
#'   n_Tree = 4,
#'   VarPlot = 0.1,
#'   Seg_Sites = 500,
#'   nQTL = 50,
#'   mu = 0,
#'   Vp = 1,
#'   Heritability = 0.3,
#'   TypePop = "CP"
#' )
#' 
#' # Access results
#' print(result$VarComp)
#' print(result$Herit.All)
#' }
#' 
#' @import dplyr
#' @import asreml
#' @import AlphaSimR
#' @import pedigreemm
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom parallel makeCluster stopCluster
#' @importFrom scales percent
#' @importFrom stats var rnorm
#' @importFrom tibble rownames_to_column
#' @export

##################################################################
#                                                               #
# Function for simulations                                      #
##################################################################

## Sub function `pin` to calculate the heritabilities and the SE's
# Load pin function
pin <- function(object, transform) {
  pframe <- as.list(summary(object)$varcomp[, 1])
  names(pframe) <- paste("V", seq(1, length(pframe)), sep = "")
  tvalue <-
    eval(deriv(transform[[length(transform)]], names(pframe)),
         pframe)
  X <- as.vector(attr(tvalue, "gradient"))
  tname <- if (length(transform) == 3)
    transform[[2]]
  else
    ""
  n <- length(pframe)
  i <- rep(1:n, 1:n)
  j <- sequence(1:n)
  k <- 1 + (i > j)
  Vmat <- object$ai
  toext <- upper.tri(Vmat)
  diag(toext) <- TRUE
  Vmat <- Vmat[which(toext, arr.ind = TRUE)]
  se <- sqrt(abs(sum(Vmat * X[i] * X[j] * k)))
  toreturn2 <- data.frame(Estimate = tvalue, SE = se)
  rownames(toreturn2) <- tname
  class(toreturn2) <- c("vpredict.asr", "data.frame")
  return(toreturn2)
}


##### Function to calculate Self-Half_Sibs (SHS) and Self-Full-Sibs (SFS) ######
# Function to calculate haplotype correspondence per individual
calculate_haplotype_similarity <- function(haplo_data, pedigree_data) {
  Sele_id <- pedigree_data |> 
    filter(
      gene == 1 & (father == mother)
    ) |> 
    pull(id)
  # Filter individuals selfing
  haplo_data <- haplo_data[gsub("_.*", "", rownames(haplo_data)) %in% Sele_id,]
  
  # Extract individual IDs
  individuals <- unique(gsub("_.*", "", rownames(haplo_data)))
  results <- data.frame(Individual = individuals, Similarity = NA)
  
  for (ind in individuals) {
    # Select haplotypes for the current individual
    haplo1 <- haplo_data[paste0(ind, "_1"), ]
    haplo2 <- haplo_data[paste0(ind, "_2"), ]
    
    # Calculate the proportion of matching loci
    matching_loci <- sum(haplo1 == haplo2)
    total_loci <- length(haplo1)
    similarity <- matching_loci / total_loci
    
    # Store the result
    results$Similarity[results$Individual == ind] <- round(similarity, 2)
  }
  
  # Classify individuals based on similarity threshold
  results$Classification <- ifelse(results$Similarity >= 0.95, "SFS", "SHS")
  return(results)
}

############# Function to Calculate homozygosity and heterozygosity ######################
calculate_homo_and_hetero <- function(geno_matrix) {
  apply(geno_matrix, 1, function(genotype) {
    n_0_loci <- sum(genotype == 0)
    n_1_loci <- sum(genotype == 1)
    n_2_loci <- sum(genotype == 2)
    homozygous_loci <- sum(genotype == 0 | genotype == 2)
    heterozygous_loci <- sum(genotype == 1)
    total_loci <- length(genotype)
    total_loci2 <- sum(genotype == 1 | genotype == 2)
    homozygosity <- homozygous_loci / total_loci
    heterozygosity <- heterozygous_loci / total_loci
    homozygosity2 <- n_2_loci / total_loci2
    heterozygosity2 = n_1_loci / total_loci2
    return(c(
      n_0_loci = n_0_loci,
      n_1_loci = n_1_loci,
      n_2_loci = n_2_loci,
      Homozygosity = homozygosity,
      Heterozygosity = heterozygosity,
      Homozygosity2 = homozygosity2,
      Heterozygosity2 = heterozygosity2))
  })
}
############### Set up Parameters ################################

TB_SimulR <- function(
    # Set up parameters for simulate pedigree
    generations = 2,
    ids = 1000,
    animals = F,
    #familySize = n_Block * n_Tree,
    #nPlott = ids * n_Block,
    ## Desirable selfing rate
    nP = 0.1,
    ## Set up number of families at F1
    nFam_F1 = 100,
    # Set up parameters for simulation
    # Parameters to set up the machine
    ## Set up the core numbers
    num_cores = NULL,
    ## set up if use N cores
    use_cores = NULL,
    # Parameters to get the pedigree
    ## Number of families: Founders in the base population
    n_Fam = NULL,
    ## Number of replication/blocks
    n_Block = NULL,
    ## Number of trees on the plot
    n_Tree = NULL,
    ## Plot sd
    VarPlot = NULL,
    # Parameter for simulations
    ## Number of the QTL affecting the trait
    nQTLPerChr = NULL,
    ## Number of SNPs in the population
    nSNPPerChr = NULL,
    ## Set up the Ne
    Ne = NULL,
    ## mean of additive effect
    mu = NULL,
    ## Phenotypic variance
    Vp = NULL,
    ## Heritability
    Heritability = NULL,
    ## Additive variance: Vp * Heritability
    Va = NULL,
    AddMat = NULL,
    ## Mean of the dominance variance
    MeanDom = NULL,
    DomMat = NULL,
    ## Dominance effec
    #Domeff = 0.2,
    ## Variance of the dominance: Vp * Domeff
    Vd = NULL,
    ## Set up number of progeny per crosses to be simulated
    nProg = NULL,
    # Choose between open-pollination or control-polination
    ## OP = Open-pollination
    ## CP = Control-pollination
    TypePop = NULL,
    # Set up parameter for extract parameters simulated
    # Get the sample of genotyped individuals from Gmatrix in the ssGBLUP models
    Samp = NULL,
    Text_track = "",
    # Want ssGBLUP model?
    ssGBLUP = TRUE,
    # Option to pull qtl, snps and haplotypes
    pullQTL = NULL,
    pullSNPs = NULL,
    pullHaplotype = TRUE,
    Get_only_data = TRUE
    
    ) {
##--------- Part 1: Function for simulate the pedigree #################
  # Set up number of progenies by Families and number of progenies by crosses if OP or CP
  cat("Population type:\n", TypePop, ": ", Text_track, "\n")
  if (!TypePop %in% c("OP", "CP")) {
    stop("TypePop have to be OP or CP")
  }
  
  if (TypePop == "OP") {
    n_progenie_by_family <- nProg
    n_Progeny_by_crosses <- 1
  }
  
  if (TypePop == "CP") {
    n_progenie_by_family <- 1
    n_Progeny_by_crosses <- nProg
  }
  
  # The parameters used to build the pedigree
  
  familySize <- n_progenie_by_family
  
  cat("\n\nstart part 1: build the pedigree\n\n")
  # if only one value is given, set same number for all generations
  if (length(ids) == 1) {
    ids <- rep(ids, times = generations)
    
    # initialisation
    gener <- rep(1:generations, times = ids)
    ID <- 1:sum(ids)
    father <- mother <- rep(0, length(ID))
  }
  # random mating for plants (inbreeds are likely)
  if (!animals) {
    for (i in 2:generations) {
      father[gener == i] <-
        ID[sample(ID[gener == i - 1], size = ids[i], replace = F)]
      mother[gener == i] <-
        ID[sample(ID[gener == i - 1], size = ids[i], replace = F)]
    }
    pedi <-
      data.frame(
        ID = ID,
        father = father,
        mother = mother,
        gener = gener - 1,
        sex = NA
      )
    # define sire and dams for animals (no inbreeds)
  } else {
    # define sex for 1st generation
    # 0 = female
    # 1 = male
    sex <- rep(0, length(ID))
    sex[gener == 1] <-
      sample(rep(c(0, 1), length = sum(gener == 1)), sum(gener == 1), replace =
               FALSE)
    
    for (i in 2:generations) {
      sex[gener == i] <-
        sample(rep(c(0, 1), length = sum(gener == i)), sum(gener == i), replace =
                 FALSE)
      
      father[gener == i] <-
        ID[sample(ID[(gener == i - 1)], size = ids[i], replace = F)]
      mother[gener == i] <-
        ID[sample(ID[(gener == i - 1)], size = ids[i], replace = F)]
    }
    pedi <-
      data.frame(
        ID = ID,
        father = father,
        mother = mother,
        gener = gener - 1,
        sex = sex
      )
  }
  
  # create open-polinated families in the last generation
  pedTemp <- pedi[pedi$gener == generations - 1, ]
  Parents <- pedi[pedi$gener < generations - 1, ]
  
  # 1st randomization of fathers
  pedFamily <- pedTemp |>
    dplyr::mutate(father = sample(father, 
                           size = length(father),
                           replace = F))
  rownames(pedFamily) <- NULL
  
  # Select number of families declared in`n_fam` argument for F1 generation
  pedFamily2 <- pedFamily[sample(
    1:nrow(pedFamily),
    size = nFam_F1,
    replace = F
  ),]
  rownames(pedFamily2) <- NULL
  
  # lenght of the F1 data
  LenF1 <- length(1:(nrow(pedFamily2) * familySize) + pedFamily$ID[1] - 1)
  # Get the F1 generation
  F1 <-
    data.frame(
      ID = rep(1:(nrow(pedFamily2) * familySize) + pedFamily$ID[1] - 1),
      mother = rep(pedFamily2$mother, each = familySize),
      # 2nd randomization of fathers
      father = sample(pedFamily$father,
                      size = LenF1,
                      replace = T),
      gener = rep(pedFamily2$gener, each = familySize),
      sex = rep(pedFamily2$sex, each = familySize)
    )
  
  # randomization of inbreeding
  Rand_inbr <- sample(1:nrow(F1),
                      size = round(nrow(F1) * nP),
                      replace = F)
  # Simulate the inbreeding
  F1[Rand_inbr, "father"] <- 
    F1[Rand_inbr, "mother"]
  
  # Filter parents present in the pedigree
  Parents2 <- Parents[Parents$ID %in% unique(c(F1$mother, F1$father)), ] 
  rownames(Parents2) <- NULL
  # Set up Generation 0 (Parents) and Generation 2 (F1)
  pedi <- rbind(Parents2, F1)
  # Include class
  class(pedi) <- c("pedigr", "data.frame")
  #rownames(pedi)
  
  # Ordering by mother
  ped2 <- pedi[order(pedi$mother),]
  rownames(ped2) <- NULL
  # Correct ID
  ped2$ID <- seq(from = 1,
                 to = length(row.names(ped2)),
                 by = 1)
  # Remove Gener and sex columns
  Pedigree <- ped2[, -c(4, 5)]
  
##------------ Parte 2: Function for simulations ####################
  cat("\n\nstart part 2: Simulation\n\n")
  # Set up number of progenies by Families and number of progenies by crosses if OP or CP
  if (!TypePop %in% c("OP", "CP")) {
    stop("TypePop have to be OP or CP")
  }
  
  if (TypePop == "OP") {
    n_progenie_by_family <- nProg
    n_Progeny_by_crosses <- 1
  }
  
  if (TypePop == "CP") {
    n_progenie_by_family <- 1
    n_Progeny_by_crosses <- nProg
  }
  
  # check inbreeding
  Inbr <- Pedigree |>
    filter(father != 0 & mother != 0 & father == mother) |>
    nrow()
  nRowPed <- Pedigree |>
    filter(father != 0 & mother != 0) |>
    nrow()
  Inbreeding <- scales::percent(Inbr / nRowPed, accuracy = 0.1)
  
  # Getting crosses from pedigree
  CROSS.PED <- Pedigree |>
    filter(father != 0 & mother != 0) |>
    dplyr::select(father, mother) |>
    as.matrix()
  
  # Get the founders --> Allow to use runMacs2
  system.time(
    founderGenomes <- runMacs2(
      nInd = n_Fam,
      nChr = 11,
      segSites =  (nSNPPerChr + nQTLPerChr),
      Ne = Ne,
      # Myburg et al., 2014
      bp = 5.82e+07,
      # average from 0.91 - 1.09
      genLen = 1,
      # Myburg et al., 2014
      mutRate = 1.2e-9,
      # Set up the historical effective population size (Ne)
      histNe = c(1000),
      # Set up the number of generations ago for effective population sizes
      histGen = c(1),
      inbred = FALSE,
      split = NULL,
      ploidy = 2L,
      returnCommand = FALSE,
      nThreads = NULL
    )
  )
  # ~ 137s
  
  # start the new population
  SP <<- SimParam$new(founderGenomes)
  # Define the trait
    # Set up additive and dominance effects
  SP$addTraitAD(
    nQtlPerChr = nQTLPerChr,
    mean = mu,
    var = Va,
    meanDD = MeanDom,
    varDD = Vd,
    useVarA = T
  )
  ## Remove trac sex
  SP$setSexes("no")
  ## Track recombinations (pedigree is also tracked)
  SP$setTrackRec(isTrackRec = TRUE)
  # Create base population
  basePop = newPop(founderGenomes)
  
  # Phenotype base population individuals
  basePop = setPheno(pop = basePop,
                     h2 = Heritability)
  # Get the crosses
  Crosses <- makeCross(pop = basePop,
                       crossPlan = CROSS.PED[, c(2, 1)],
                       nProgeny = n_Progeny_by_crosses)
  
  # Set up the phenotype for the progenies
  Crosses = setPheno(pop = Crosses,
                     h2 = Heritability,
                     H2 = (Va + Vd) / Vp)
  
  # Additive variance
  Va_basePop <- varA(basePop)
  Va_Crosses <- varA(Crosses)
  # Dominance variance
  Vd_basePop <- varD(basePop)
  Vd_Crosses <- varD(Crosses)
  # Genotypic variance
  Vg_basePop <- varG(basePop)
  Vg_Crosses <- varG(Crosses)
  # Phenotypic variance
  Vp_basePop <- varP(basePop)
  Vp_Crosses <- varP(Crosses)
  # Residual variances
  Ve_basePop <- Vp_basePop - Vg_basePop
  Ve_Crosses <- Vp_Crosses - Vg_Crosses
  # Narrow sense heritabilities
  h2a_basePop <- Va_basePop / Vp_basePop
  h2a_Crosses <- Va_Crosses / Vp_Crosses
  # Coefficient of Dominance
  h2d_basePop <- Vd_basePop / Vp_basePop
  h2d_Crosses <- Vd_Crosses / Vp_Crosses
  # Broad sense heritability
  h2g_basePop <- Vg_basePop / Vp_basePop
  h2g_Crosses <- Vg_Crosses / Vp_Crosses
  
  # Get the table with estimates
  GenVars <- data.frame(
    Pop = rep(c("basePop", "F1"), each = 8),
    VarType = rep(c(
      "Va", "Vd", "Ve", "Vg", "Vp", "h2a", "h2d", "h2g"
    ), 2),
    Estimates = c(
      Va_basePop,
      Vd_basePop,
      Ve_basePop,
      Vg_basePop,
      Vp_basePop,
      h2a_basePop,
      h2d_basePop,
      h2g_basePop,
      Va_Crosses,
      Vd_Crosses,
      Ve_Crosses,
      Vg_Crosses,
      Vp_Crosses,
      h2a_Crosses,
      h2d_Crosses,
      h2g_Crosses
    )
  )
  
  # Get the pedigree of progenies - F1
  Pedigree.F1 <- data.frame(id = Crosses@id,
                            father = Crosses@father,
                            mother = Crosses@mother)
  
  Pedigree.F1 <- with(
    Pedigree.F1,
    pedigreemm::editPed(
      sire = as.numeric(father),
      dam = as.numeric(mother),
      label = as.numeric(id)
    )
  ) |>
    arrange(as.numeric(label)) |>
    rename(id = label,
           father = sire,
           mother = dam) |>
    dplyr::mutate(
      id = ifelse(is.na(id), "0", id),
      father = ifelse(is.na(father), "0", father),
      mother = ifelse(is.na(mother), "0", mother)
    )
  
  # Pull QTL?
  if (pullQTL) {
    # Get the QTL
    QTL <- pullQtlGeno(Crosses)
  } else {
    QTL <- pullQTL
  }
  
  # Add snp chip
  SP$addSnpChip(nSnpPerChr = nSNPPerChr) 
  
  SNPs <- pullSegSiteGeno(Crosses)
  
  if (pullSNPs) {
    # Get the SNPS
    SNPs_export <- SNPs
  } else {
    SNPs_export <- pullSNPs
  }
  
  if (pullHaplotype) {
    
    # Get the haplotypes
    haplotypes <- pullSegSiteHaplo(Crosses)
  } else {
    haplotypes <- pullHaplotype
  }
  
  # Get data for:
  # Get the Genetic values
  GV_basePop <- gv(basePop)
  GV <- gv(Crosses)
  # Breeding values
  BV_basePop <- bv(basePop)
  BV <- bv(Crosses)
  # Phenotype values
  PhenoValues_basePop <- pheno(basePop)
  PhenoValues <- pheno(Crosses)
  # Additive effect
  A_eff_basePop <- aa(basePop)
  A_eff <- aa(Crosses)
  # Dominance effect
  D_eff_basePop <- dd(basePop)
  D_eff <- dd(Crosses)
  # Pull QTL map
  QTL_Map <- getQtlMap(trait = 1,
                       sex = "A",
                       simParam = NULL)
  # Simulating data set
  ## Block effect
  if (TypePop == "CP") {
    Data <- Pedigree.F1 |>
      filter(gene == 1) |>
      dplyr::mutate(
        Fam = paste(mother, father, sep = "_"),
        Block = rep(1:n_Block, each = n_Tree, times = nFam_F1),
        Tree = rep(1:n_Tree, times = nFam_F1 * n_Block)
      )
  } else {
    Data <- Pedigree.F1 |>
      filter(gene == 1) |>
      dplyr::mutate(
        Fam = paste0(mother, "_"),
        Block = rep(1:n_Block, each = n_Tree, times = nFam_F1),
        Tree = rep(1:n_Tree, times = nFam_F1 * n_Block)
      )
  }
  
  # Include the Plot effect
  if (TypePop == "CP") {
    Data[["Plot"]] <-
      paste0("fam", Data$Fam, "block", Data$Block)
    
  } else {
    Data[["Plot"]] <-
      paste0("fam", Data$mother,"block", Data$Block)
  }
  
  ## Insert values for Plot effect
  PlotID <- data.frame(Plot = unique(Data$Plot),
                       Plot_eff = rnorm(nFam_F1 * n_Block, mean = 0, sd = VarPlot))
  Data <- left_join(Data,
                    PlotID,
                    by = "Plot") |>
    # Include all effects
    dplyr::mutate(
      Pheno.initial = as.vector(PhenoValues),
      GV = as.vector(GV),
      BV = as.vector(BV),
      E_eff = Pheno.initial - GV,
      Dom_eff = as.vector(D_eff),
      Pheno = Pheno.initial + Plot_eff,
      # Check the Pheno as the difference between (mean + BV + Dom + E_eff) - Pheno
      Pheno_check = mean(Pheno.initial, na.rm = T) + BV + Dom_eff + E_eff + Plot_eff - Pheno
    )
  
  # Stop simulations if wanna only data
  
  if (Get_only_data) {
    # Run the function on IBD haplotypes
    similarity_results <- calculate_haplotype_similarity(
      haplo_data = haplotypes,
      pedigree_data = Pedigree.F1
    )
    
    Results <- list(
      TypePop = TypePop,
      `% of genotyped individuals` = scales::percent(as.numeric(Samp)),
      # Results Pedigree and SNP based models
      SelfRate = Inbreeding,
      RealParameters = GenVars,
      Data = Data,
      QTL = QTL,
      snps = SNPs_export,
      haplotypes = haplotypes,
      Self_relationships = similarity_results,
      Pedigree = Pedigree.F1
    )
    
  } else {
    
  # Calculating the genomic matrices: Gmatrix and GD matrix
  if (use_cores == TRUE) {
    require(parallel)
    require(foreach)
    require(doParallel)
    # Set up the machine to run in parallel
    cl <- makeCluster(num_cores) 
    registerDoParallel(cl)
    # Export the Gmatrix function to the parallel workers
    clusterExport(cl, c("Gmatrix"))
    # Get the Gmatrix
    system.time({
      G <- foreach(i = 1:num_cores, .combine = cbind) %dopar% {
        Gmatrix(
          SNPs[, ((i - 1) * (ncol(SNPs) / num_cores) + 1):(i * (ncol(SNPs) / num_cores))],
          method = AddMat,
          maf = 0.01,
          thresh.missing = 0.01
        )
      }
    })
    
    # ~ 24s and the file size in global environment is 256.4 MB with 8 core
    # Get the Genomic-Dominance matrix: GDmatrix
    system.time({
      GD <- foreach(i = 1:num_cores, .combine = cbind) %dopar% {
        Gmatrix(
          SNPs[, ((i - 1) * (ncol(SNPs) / num_cores) + 1):(i * (ncol(SNPs) / num_cores))],
          method = DomMat,
          maf = 0.01,
          thresh.missing = 0.01
        )
      }
    })
    # ~ 24s and the file size in global environment is 256.4 MB with 8 core
    
    # Stop the parallel backend
    stopCluster(cl)
  } else {
    # Get the Gmatrix
    system.time(G <- Gmatrix(
      SNPs,
      method = AddMat,
      maf = 0.01,
      thresh.missing = 0.01
    ))
    
    # ~ 48s and the file size in global environment is 32.3 MB
    # Get the Genomic-Dominance matrix: GDmatrix
    system.time(GD <- Gmatrix(
      SNPs,
      method = DomMat,
      maf = 0.01,
      thresh.missing = 0.01
    ))
    # ~ ~ 49s and the file size in global environment is 32.3 MB
  }
  
    # A matrix: pedigree-based
    ## For Half-sibs: excluding information of fathers
    A.HS <- Amatrix(Pedigree.F1[, 1:3] |>
                      dplyr::mutate(father = "0"))
    ## For Full-sibs: with information of fathers
    if (TypePop == "OP") {
      A.FS <- Amatrix(Pedigree.F1[, 1:3] |> 
                        dplyr::mutate(father = "0"))
    } else {
      A.FS <- Amatrix(Pedigree.F1[, 1:3])
    }
  
  # Calculatin the Dominance matrix (pedigree-based)
  if (use_cores == TRUE) {
    # Register parallel backend
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    
    # Export the Amatrix function to the parallel workers
    clusterExport(cl, c("Amatrix"))
    
    # Parallelize the Amatrix function for full-sibs
    if (TypePop == "CP") {
      system.time({
        D.FS <- foreach(i = 1:num_cores, .combine = cbind) %dopar% {
          Amatrix(Pedigree.F1[i, 1:3],
                  dominance = TRUE)
        }
      })
    } else {
      D.FS <- diag(
        x = 1,
        nrow = dim(A.HS)[1],
        ncol = dim(A.HS)[2]
      )
      # Include dim
      dimnames(D.FS) <- dimnames(A.HS)
      
    }
    
    # Stop the parallel backend
    stopCluster(cl)
  } else {
    if (TypePop == "CP") {
      ## For Full-sibs: with information of fathers
      system.time(D.FS <- Amatrix(Pedigree.F1[, 1:3],
                                  dominance = T))
    } else {
      D.FS <- diag(
        x = 1,
        nrow = dim(A.HS)[1],
        ncol = dim(A.HS)[2]
      )
      # Include dim
      dimnames(D.FS) <- dimnames(A.HS)
    }
    
    # ~94s
  }
  
  ##--------------- Parte 3: Extract the paremeter simulated #################
  cat("\n\nstart part 3: Extract the parameter simulated\n\n")
  # Half-sibs
  ## Amatrix
  if (TypePop == "OP") {
    A.inv.HS <<- ainverse(Pedigree.F1[, 1:3] |>
                            dplyr::mutate(father = "0"))
  }
  
  # Full-sibs
  ## Amatrix
  if (TypePop == "OP") {
    ### Remember, if OP population, the FS model will be the same as HS model because the father it's supposed unknown
    A.inv.FS <<- ainverse(Pedigree.F1[, 1:3] |>
                        dplyr::mutate(father = "0"))
  } else {
    A.inv.FS <<- ainverse(Pedigree.F1[, 1:3])
  }
  
  
  #unique(c(Simul$D.FS))
  ## Dmatrix
  D.inv <<- G.inverse(G = D.FS, sparseform = TRUE)$Ginv
  ##
  #GD <<- GD
  # Gmatrix
  ## Alignment
  G2A <-
    match.G2A(
      A = A.FS,
      G = G,
      clean = TRUE,
      ord = TRUE,
      mism = TRUE,
      RMdiff = TRUE
    )
  
  G2D <-
    match.G2A(
      A = D.FS,
      G = GD,
      clean = TRUE,
      ord = TRUE,
      mism = TRUE,
      RMdiff = TRUE
    )
  ## Get tuned Gmatrix
  G_tuned <<- G.tuneup(G = G2A$Gclean,
                       A = G2A$Aclean,
                       align = T,
                       bend = F
                       )$Gb
  
  GD_tuned <<- G.tuneup(
    G = G2D$Gclean,
    A = G2D$Aclean,
    align = T,
    bend = F,
    blend = F
  )$Gb
  
  ## Gmatrix
  G.inv <- G.inverse(G = G_tuned,
                      sparseform = FALSE,
                      rcn.thr = 1e-30)$Ginv
  # GDmatrix
  GD.inv <- G.inverse(G = G2D$Gclean,
                      sparseform = FALSE,
                      rcn.thr = 1e-30)$Ginv

  ## H
  system.time(H_full <<-
                H.inverse(
                  A = A.HS,
                  G = G.inv,
                  #tau = 2,
                  #omega = 0.1,
                  lambda = 0.9,  
                  sparseform = TRUE
                ))
  
  
  ## HD
  system.time(HD_full <<-
                H.inverse(
                  A = D.FS,
                  G = GD.inv,
                  #tau = 2,
                  #omega = 0.1,
                  lambda = 0.9,  
                  sparseform = FALSE
                ))
  
  # Set up model if single-tree plot
  if (n_Tree == 1) {
    # Pedigree-based model: Half-sibs
    ## Family model
    HS.Fam.model <- asreml(
      fixed = Pheno ~ 1 + Block,
      random = ~ Fam,
      data = Data |>
        dplyr::mutate(
          Block = as.factor(Block),
          Plot = as.factor(Plot),
          Fam = as.factor(Fam)
        )
    )
    HS.Fam.model <- update.asreml(HS.Fam.model)
    if (TypePop == "OP") {
      ## Additive model
      HS.Add.model <- try({asreml(
        fixed = Pheno ~ 1 + Block,
        random = ~ vm(id, A.inv.HS),
        data = Data |>
          dplyr::mutate(
            Block = as.factor(Block),
            Plot = as.factor(Plot),
            id = as.factor(id)
          )
      )})
      HS.Add.model <- try({update.asreml(HS.Add.model)})
    } else {
      HS.Add.model <- NA
    }
    
    
    # Pedigree-based model: Full-sibs
    FS.Add.model <- try({asreml(
      fixed = Pheno ~ 1 + Block,
      random = ~ vm(id, A.inv.FS),
      data = Data |>
        dplyr::mutate(
          Block = as.factor(Block),
          Plot = as.factor(Plot),
          id = as.factor(id)
        )
    )})
    FS.Add.model <- try({update.asreml(FS.Add.model)})
    # AD
    if (TypePop == "OP") {
      FS.AddDom.model <- try({asreml(
        fixed = Pheno ~ 1 + Block,
        random = ~ vm(id, A.inv.FS),
        data = Data |>
          dplyr::mutate(
            Block = as.factor(Block),
            Plot = as.factor(Plot),
            id = as.factor(id),
            id2 = as.factor(id)
          )
      )})
    } else {
      FS.AddDom.model <- try({asreml(
        fixed = Pheno ~ 1 + Block,
        random = ~ vm(id, A.inv.FS) + vm(id2, D.inv),
        data = Data |>
          dplyr::mutate(
            Block = as.factor(Block),
            Plot = as.factor(Plot),
            id = as.factor(id),
            id2 = as.factor(id)
          )
      )})
    }
    
    FS.AddDom.model <- try({update.asreml(FS.AddDom.model)})
    
    # Genomic-based model
    ## GBlup model
    GBlup.model <- asreml(
      fixed = Pheno ~ 1 + Block,
      random = ~ vm(id, G_tuned) ,
      data = Data |>
        dplyr::mutate(
          Block = as.factor(Block),
          Plot = as.factor(Plot),
          id = as.factor(id),
          id2 = as.factor(id)
        ),
      workspace = 200e05
    )
    #GBlup.model <- update.asreml(GBlup.model)
    ## GDBlup model
    GDBlup.model <- asreml(
      fixed = Pheno ~ 1 + Block,
      random = ~ vm(id, G_tuned) + vm(id2, GD),
      data = Data |>
        dplyr::mutate(
          Block = as.factor(Block),
          Plot = as.factor(Plot),
          id = as.factor(id),
          id2 = as.factor(id)
        ),
      workspace = 300e06
    )
    
    ## HDBlup model
    HDBlup.model <- asreml(
      fixed = Pheno ~ 1 + Block,
      random = ~ vm(id, H_full) + vm(id2, HD_full),
      data = Data |>
        dplyr::mutate(
          Block = as.factor(Block),
          Plot = as.factor(Plot),
          id = as.factor(id),
          id2 = as.factor(id)
        ),
      workspace = 300e06
    )
    
    # If not a single-tree-plot model
  } else {
    # Pedigree-based model: Half-sibs
    ## Family model
    HS.Fam.model <- asreml(
      fixed = Pheno ~ 1 + Block,
      random = ~ Plot + Fam,
      data = Data |>
        dplyr::mutate(
          Block = as.factor(Block),
          Plot = as.factor(Plot),
          Fam = as.factor(Fam)
        )
    )
    HS.Fam.model <- update.asreml(HS.Fam.model)
    if (TypePop == "OP") {
      ## Additive model
      HS.Add.model <- try({asreml(
        fixed = Pheno ~ 1 + Block,
        random = ~ Plot + vm(id, A.inv.HS),
        data = Data |>
          dplyr::mutate(
            Block = as.factor(Block),
            Plot = as.factor(Plot),
            id = as.factor(id)
          )
      )})
      HS.Add.model <- try({update.asreml(HS.Add.model)})
    } else {
      HS.Add.model <- NA
    }
    
    
    # Pedigree-based model: Full-sibs
    FS.Add.model <- try({asreml(
      fixed = Pheno ~ 1 + Block,
      random = ~ Plot + vm(id, A.inv.FS),
      data = Data |>
        dplyr::mutate(
          Block = as.factor(Block),
          Plot = as.factor(Plot),
          id = as.factor(id)
        )
    )})
    FS.Add.model <- try({update.asreml(FS.Add.model)})
    # AD
    if (TypePop == "OP") {
      FS.AddDom.model <- try({asreml(
        fixed = Pheno ~ 1 + Block,
        random = ~ Plot + vm(id, A.inv.FS),
        data = Data |>
          dplyr::mutate(
            Block = as.factor(Block),
            Plot = as.factor(Plot),
            id = as.factor(id),
            id2 = as.factor(id)
          )
      )})
    } else {
      FS.AddDom.model <- try({asreml(
        fixed = Pheno ~ 1 + Block,
        random = ~ Plot + vm(id, A.inv.FS) + vm(id2, D.inv),
        data = Data |>
          dplyr::mutate(
            Block = as.factor(Block),
            Plot = as.factor(Plot),
            id = as.factor(id),
            id2 = as.factor(id)
          )
      )})
    }
    
    FS.AddDom.model <- try({update.asreml(FS.AddDom.model)})
    
    # Genomic-based model
    ## GBlup model
    GBlup.model <- asreml(
      fixed = Pheno ~ 1 + Block,
      random = ~ Plot + vm(id, G_tuned) ,
      data = Data |>
        dplyr::mutate(
          Block = as.factor(Block),
          Plot = as.factor(Plot),
          id = as.factor(id),
          id2 = as.factor(id)
        ),
      workspace = 200e05
    )
    #GBlup.model <- update.asreml(GBlup.model)
    ## GDBlup model
    GDBlup.model <- asreml(
      fixed = Pheno ~ 1 + Block,
      random = ~ Plot + vm(id, G_tuned) + vm(id2, GD_tuned),
      data = Data |>
        dplyr::mutate(
          Block = as.factor(Block),
          Plot = as.factor(Plot),
          id = as.factor(id),
          id2 = as.factor(id)
        ),
      workspace = 300e06
    )
    ## HDBlup model
    HDBlup.model <- asreml(
      fixed = Pheno ~ 1 + Block,
      random = ~ Plot + vm(id, H_full) + vm(id2, HD_full),
      data = Data |>
        dplyr::mutate(
          Block = as.factor(Block),
          Plot = as.factor(Plot),
          id = as.factor(id),
          id2 = as.factor(id)
        ),
      workspace = 300e06
    )
  }
  
  
  if (TypePop == "OP" & class(HS.Add.model) != "try-error") {
    HS.Fam <- summary(HS.Fam.model)$varcomp
    HS.Add <- summary(HS.Add.model)$varcomp
  } else {
    HS.Fam <- NA
    HS.Add <- NA
  }
  
  if (class(FS.Add.model) != "try-error") {
    FS.Add <- summary(FS.Add.model)$varcomp
  } else {
    FS.Add <- NA
  }
  
  if (class(FS.AddDom.model) != "try-error") {
    FS.AddDom <- summary(FS.AddDom.model)$varcomp
  } else {
    FS.AddDom <- NA
  } 
  
  GBlup <- summary(GBlup.model)$varcomp
  GDBlup <- summary(GDBlup.model)$varcomp
  HDBlup <- summary(HDBlup.model)$varcomp
  
  VarComp <- list(
    HS.Fam = HS.Fam,
    HS.Add = HS.Add,
    FS.Add = FS.Add,
    FS.AddDom = FS.AddDom,
    GBlup = GBlup,
    GDBlup = GDBlup,
    HDBlup = HDBlup
  )
  
  ## Extract heritabilities if single-tree-plot
  if (n_Tree == 1) {
    if (TypePop == "OP") {
      Herit <- rbind(
        # h2a
        data.frame(
          h2.type = "h2a",
          Estimate = round(GenVars[which(GenVars$Pop == "F1" &
                                     GenVars$VarType == "h2a"),
                             "Estimates"], 2
          ),
          SE = NA,
          TypePop = "True"
        ),
        pin(HS.Fam.model, h2a ~ 4 * V1 / (V1 + V2)) |>
          dplyr::mutate(TypePop = "Fam") |>
          rownames_to_column(var = "h2.type"),
        pin(HS.Add.model, h2a ~ V1 / (V1 + V2)) |>
          dplyr::mutate(TypePop = "HS") |>
          rownames_to_column(var = "h2.type"),
        pin(FS.Add.model, h2a ~ V1 / (V1 + V2)) |>
          dplyr::mutate(TypePop = "FS.Add") |>
          rownames_to_column(var = "h2.type"),
        pin(FS.AddDom.model, h2a ~ V1 / (V1 + V2)) |>
          dplyr::mutate(TypePop = "FS.AddDom") |>
          rownames_to_column(var = "h2.type"),
        pin(GBlup.model, h2a ~ V1 / (V1 + V2)) |>
          dplyr::mutate(TypePop = "GBlup") |>
          rownames_to_column(var = "h2.type"),
        pin(GDBlup.model, h2a ~ V1 / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "GDBlup") |>
          rownames_to_column(var = "h2.type"),
        pin(HDBlup.model, h2a ~ V1 / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "ssGBlup100") |>
          rownames_to_column(var = "h2.type"),
        # h2d
        data.frame(
          h2.type = "h2d",
          Estimate = round(GenVars[which(GenVars$Pop == "F1" &
                                           GenVars$VarType == "h2d"),
                                   "Estimates"], 2
          ),
          SE = NA,
          TypePop = "True"
        ),
        pin(FS.AddDom.model, h2d ~ V1 / (V1 + V2)) |>
          dplyr::mutate(TypePop = "FS.AddDom") |>
          rownames_to_column(var = "h2.type"),
        pin(GDBlup.model, h2d ~ V2 / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "GDBlup") |>
          rownames_to_column(var = "h2.type"),
        pin(HDBlup.model, h2d ~ V2 / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "ssGBlup100") |>
          rownames_to_column(var = "h2.type"),
        # h2g
        data.frame(
          h2.type = "h2g",
          Estimate = round(GenVars[which(GenVars$Pop == "F1" &
                                           GenVars$VarType == "h2g"),
                                   "Estimates"], 2
          ),
          SE = NA,
          TypePop = "True"
        ),
        pin(HS.Fam.model, h2g ~ V1 / (V1 + V2)) |>
          dplyr::mutate(TypePop = "HS") |>
          rownames_to_column(var = "h2.type"),
        pin(FS.AddDom.model, h2g ~ (V1) / (V1 + V2)) |>
          dplyr::mutate(TypePop = "FS.AddDom") |>
          rownames_to_column(var = "h2.type"),
        pin(GDBlup.model, h2g ~ (V1 + V2) / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "GDBlup") |>
          rownames_to_column(var = "h2.type"),
        pin(HDBlup.model, h2g ~ (V1 + V2) / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "ssGBlup100") |>
          rownames_to_column(var = "h2.type")
      ) |>
        dplyr::mutate(Estimate = round(Estimate, 2),
                      SE = round(SE, 3)) |> 
        # Include Selfing and number of genotyped individuals
        dplyr::mutate(
          SelfRate = nP,
          `% genotyped` = Samp
        )
      # If CP
    } else {
      Herit <- rbind(
        # h2a
        data.frame(
          h2.type = "h2a",
          Estimate = round(GenVars[which(GenVars$Pop == "F1" &
                                           GenVars$VarType == "h2a"),
                                   "Estimates"], 2
          ),
          SE = NA,
          TypePop = "True"
        ),
        pin(HS.Fam.model, h2a ~ 4 * V1 / (V1 + V2)) |>
          dplyr::mutate(TypePop = "Fam") |>
          rownames_to_column(var = "h2.type"),
        pin(FS.Add.model, h2a ~ V1 / (V1 + V2)) |>
          dplyr::mutate(TypePop = "FS.Add") |>
          rownames_to_column(var = "h2.type"),
        pin(FS.AddDom.model, h2a ~ V1 / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "FS.AddDom") |>
          rownames_to_column(var = "h2.type"),
        pin(GBlup.model, h2a ~ V1 / (V1 + V2)) |>
          dplyr::mutate(TypePop = "GBlup") |>
          rownames_to_column(var = "h2.type"),
        pin(GDBlup.model, h2a ~ V1 / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "GDBlup") |>
          rownames_to_column(var = "h2.type"),
        pin(HDBlup.model, h2a ~ V1 / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "ssGBlup100") |>
          rownames_to_column(var = "h2.type"),
        # h2d
        data.frame(
          h2.type = "h2d",
          Estimate = round(GenVars[which(GenVars$Pop == "F1" &
                                           GenVars$VarType == "h2d"),
                                   "Estimates"], 2
          ),
          SE = NA,
          TypePop = "True"
        ),
        pin(FS.AddDom.model, h2d ~ V2 / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "FS.AddDom") |>
          rownames_to_column(var = "h2.type"),
        pin(GDBlup.model, h2d ~ V2 / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "GDBlup") |>
          rownames_to_column(var = "h2.type"),
        pin(HDBlup.model, h2d ~ V2 / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "ssGBlup100") |>
          rownames_to_column(var = "h2.type"),
        # h2g
        data.frame(
          h2.type = "h2g",
          Estimate = round(GenVars[which(GenVars$Pop == "F1" &
                                           GenVars$VarType == "h2g"),
                                   "Estimates"], 2
          ),
          SE = NA,
          TypePop = "True"
        ),
        pin(HS.Fam.model, h2g ~ V1 / (V1 + V2)) |>
          dplyr::mutate(TypePop = "HS") |>
          rownames_to_column(var = "h2.type"),
        pin(FS.AddDom.model, h2g ~ (V1 + V2) / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "FS.AddDom") |>
          rownames_to_column(var = "h2.type"),
        pin(GDBlup.model, h2g ~ (V1 + V2) / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "GDBlup") |>
          rownames_to_column(var = "h2.type"),
        pin(HDBlup.model, h2g ~ (V1 + V2) / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "ssGBlup100") |>
          rownames_to_column(var = "h2.type")
      ) |>
        dplyr::mutate(Estimate = round(Estimate, 2),
                      SE = round(SE, 3)) |> 
        # Include Selfing and number of genotyped individuals
        dplyr::mutate(
          SelfRate = nP,
          `% genotyped` = Samp
        )
    }
    ## Extract heritabilities if not a single-tree-plot
  } else {
    if (TypePop == "OP") {
      Herit <- rbind(
        # h2a
        data.frame(
          h2.type = "h2a",
          Estimate = round(GenVars[which(GenVars$Pop == "F1" &
                                           GenVars$VarType == "h2a"),
                                   "Estimates"], 2
          ),
          SE = NA,
          TypePop = "True"
        ),
        pin(HS.Fam.model, h2a ~ 4 * V2 / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "Fam") |>
          rownames_to_column(var = "h2.type"),
        pin(HS.Add.model, h2a ~ V2 / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "HS") |>
          rownames_to_column(var = "h2.type"),
        pin(FS.Add.model, h2a ~ V2 / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "FS.Add") |>
          rownames_to_column(var = "h2.type"),
        pin(FS.AddDom.model, h2a ~ V2 / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "FS.AddDom") |>
          rownames_to_column(var = "h2.type"),
        pin(GBlup.model, h2a ~ V2 / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "GBlup") |>
          rownames_to_column(var = "h2.type"),
        pin(GDBlup.model, h2a ~ V2 / (V1 + V2 + V3 + V4)) |>
          dplyr::mutate(TypePop = "GDBlup") |>
          rownames_to_column(var = "h2.type"),
        pin(HDBlup.model, h2a ~ V2 / (V1 + V2 + V3 + V4)) |>
          dplyr::mutate(TypePop = "ssGBlup") |>
          rownames_to_column(var = "h2.type"),
        # h2d
        data.frame(
          h2.type = "h2d",
          Estimate = round(GenVars[which(GenVars$Pop == "F1" &
                                           GenVars$VarType == "h2d"),
                                   "Estimates"], 2
          ),
          SE = NA,
          TypePop = "True"
        ),
        pin(FS.AddDom.model, h2d ~ V2 / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "FS.AddDom") |>
          rownames_to_column(var = "h2.type"),
        pin(GDBlup.model, h2d ~ V3 / (V1 + V2 + V3 + V4)) |>
          dplyr::mutate(TypePop = "GDBlup") |>
          rownames_to_column(var = "h2.type"),
        pin(HDBlup.model, h2d ~ V3 / (V1 + V2 + V3 + V4)) |>
          dplyr::mutate(TypePop = "ssGBlup") |>
          rownames_to_column(var = "h2.type"),
        # h2g
        data.frame(
          h2.type = "h2g",
          Estimate = round(GenVars[which(GenVars$Pop == "F1" &
                                           GenVars$VarType == "h2g"),
                                   "Estimates"], 2
          ),
          SE = NA,
          TypePop = "True"
        ),
        pin(HS.Fam.model, h2g ~ V2 / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "HS") |>
          rownames_to_column(var = "h2.type"),
        pin(FS.AddDom.model, h2g ~ (V2) / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "FS.AddDom") |>
          rownames_to_column(var = "h2.type"),
        pin(GDBlup.model, h2g ~ (V2 + V3) / (V1 + V2 + V3 + V4)) |>
          dplyr::mutate(TypePop = "GDBlup") |>
          rownames_to_column(var = "h2.type"),
        pin(HDBlup.model, h2g ~ (V2 + V3) / (V1 + V2 + V3 + V4)) |>
          dplyr::mutate(TypePop = "ssGBlup") |>
          rownames_to_column(var = "h2.type")
      ) |>
        dplyr::mutate(Estimate = round(Estimate, 2),
                      SE = round(SE, 3)) |> 
        # Include Selfing and number of genotyped individuals
        dplyr::mutate(
          SelfRate = nP,
          `% genotyped` = Samp
        )
    } else {
      Herit <- rbind(
        # h2a
        data.frame(
          h2.type = "h2a",
          Estimate = round(GenVars[which(GenVars$Pop == "F1" &
                                           GenVars$VarType == "h2a"),
                                   "Estimates"], 2
          ),
          SE = NA,
          TypePop = "True"
        ),
        pin(HS.Fam.model, h2a ~ 4 * V2 / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "Fam") |>
          rownames_to_column(var = "h2.type"),
        pin(FS.Add.model, h2a ~ V2 / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "FS.Add") |>
          rownames_to_column(var = "h2.type"),
        pin(FS.AddDom.model, h2a ~ V2 / (V1 + V2 + V3 + V4)) |>
          dplyr::mutate(TypePop = "FS.AddDom") |>
          rownames_to_column(var = "h2.type"),
        pin(GBlup.model, h2a ~ V2 / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "GBlup") |>
          rownames_to_column(var = "h2.type"),
        pin(GDBlup.model, h2a ~ V2 / (V1 + V2 + V3 + V4)) |>
          dplyr::mutate(TypePop = "GDBlup") |>
          rownames_to_column(var = "h2.type"),
        pin(HDBlup.model, h2a ~ V2 / (V1 + V2 + V3 + V4)) |>
          dplyr::mutate(TypePop = "ssGBlup100") |>
          rownames_to_column(var = "h2.type"),
        # h2d
        data.frame(
          h2.type = "h2d",
          Estimate = round(GenVars[which(GenVars$Pop == "F1" &
                                           GenVars$VarType == "h2d"),
                                   "Estimates"], 2
          ),
          SE = NA,
          TypePop = "True"
        ),
        pin(FS.AddDom.model, h2d ~ V3 / (V1 + V2 + V3 + V4)) |>
          dplyr::mutate(TypePop = "FS.AddDom") |>
          rownames_to_column(var = "h2.type"),
        pin(GDBlup.model, h2d ~ V3 / (V1 + V2 + V3 + V4)) |>
          dplyr::mutate(TypePop = "GDBlup") |>
          rownames_to_column(var = "h2.type"),
        pin(HDBlup.model, h2d ~ V3 / (V1 + V2 + V3 + V4)) |>
          dplyr::mutate(TypePop = "ssGBlup100") |>
          rownames_to_column(var = "h2.type"),
        # h2g
        data.frame(
          h2.type = "h2g",
          Estimate = round(GenVars[which(GenVars$Pop == "F1" &
                                           GenVars$VarType == "h2g"),
                                   "Estimates"], 2
          ),
          SE = NA,
          TypePop = "True"
        ),
        pin(HS.Fam.model, h2g ~ V2 / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "HS") |>
          rownames_to_column(var = "h2.type"),
        pin(FS.AddDom.model, h2g ~ (V2 + V3) / (V1 + V2 + V3 + V4)) |>
          dplyr::mutate(TypePop = "FS.AddDom") |>
          rownames_to_column(var = "h2.type"),
        pin(GDBlup.model, h2g ~ (V2 + V3) / (V1 + V2 + V3 + V4)) |>
          dplyr::mutate(TypePop = "GDBlup") |>
          rownames_to_column(var = "h2.type"),
        pin(HDBlup.model, h2g ~ (V2 + V3) / (V1 + V2 + V3 + V4)) |>
          dplyr::mutate(TypePop = "ssGBlup100") |>
          rownames_to_column(var = "h2.type")
      ) |>
        dplyr::mutate(Estimate = round(Estimate, 2),
                      SE = round(SE, 3)) |> 
        # Include Selfing and number of genotyped individuals
        dplyr::mutate(
          SelfRate = nP,
          `% genotyped` = Samp
        )
    } 
  }
  
  
  if (TypePop == "OP") {
    # BLUPs from HS
    Ranef.HS.Add <-
      summary(HS.Add.model, coef = T)$coef.random[grep("vm", rownames(summary(HS.Add.model, coef = T)$coef.random)),] %>%
      `rownames<-`(gsub("*.*_", "", rownames(.))) |>
      data.frame() |>
      rownames_to_column(var = "id") |>
      filter(id %in% Data$id) |> 
      # include PEV
      dplyr::mutate(
        PEV = std.error ^ 2,
        # Calculate accuracy
        raa = sqrt(1 - PEV / VarComp$HS.Add[str_detect(rownames(VarComp$HS.Add),
                                                       "id,"), "component"])
      )
    
  } else {
    Ranef.HS.Add <- NA
  }
  
  
  # BLUPs from FS
  Ranef.FS.Add <-
    summary(FS.Add.model, coef = T)$coef.random[grep("vm", rownames(summary(FS.Add.model, coef = T)$coef.random)),] %>%
    `rownames<-`(gsub("*.*_", "", rownames(.))) |>
    data.frame() |>
    rownames_to_column(var = "id") |>
    filter(id %in% Data$id) |> 
    # include PEV
    dplyr::mutate(
      PEV = std.error ^ 2,
      # Calculate accuracy
      raa = sqrt(1 - PEV / VarComp$FS.Add[str_detect(rownames(VarComp$FS.Add),
                                                     "id,"), "component"])
    )
  
  # from FS AddDom
  Ranef.FS.AddDom <-
    summary(FS.AddDom.model, coef = T)$coef.random[grep("A.inv.FS", rownames(summary(FS.AddDom.model, coef = T)$coef.random)),] %>%
    `rownames<-`(gsub("*.*_", "", rownames(.))) |>
    data.frame() |>
    rownames_to_column(var = "id") |>
    filter(id %in% Data$id) |> 
    # include PEV
    dplyr::mutate(
      PEV = std.error ^ 2,
      # Calculate accuracy
      raa = sqrt(1 - PEV / VarComp$FS.AddDom[str_detect(rownames(VarComp$FS.AddDom),
                                                     "id,"), "component"])
    )
  
  # BLUPs from GBlup
  Ranef.GBlup <-
    summary(GBlup.model, coef = T)$coef.random[grep("vm", rownames(summary(GBlup.model, coef = T)$coef.random)),] %>%
    `rownames<-`(gsub("*.*_", "", rownames(.))) |>
    data.frame() |>
    rownames_to_column(var = "id") |>
    filter(id %in% Data$id) |> 
    # include PEV
    dplyr::mutate(
      PEV = std.error ^ 2,
      # Calculate accuracy
      raa = sqrt(1 - PEV / VarComp$GBlup[str_detect(rownames(VarComp$GBlup),
                                                     "id,"), "component"])
    )
  
  # BLUPs from GDBlup
  Ranef.GDBlup <-
    summary(GDBlup.model, coef = T)$coef.random[grep("G_tuned", rownames(summary(GDBlup.model, coef = T)$coef.random)),] %>%
    `rownames<-`(gsub("*.*_", "", rownames(.))) |>
    data.frame() |>
    rownames_to_column(var = "id") |>
    filter(id %in% Data$id) |> 
    # include PEV
    dplyr::mutate(
      PEV = std.error ^ 2,
      # Calculate accuracy
      raa = sqrt(1 - PEV / VarComp$GDBlup[str_detect(rownames(VarComp$GDBlup),
                                                     "id,"), "component"])
    )
  
  # Blups from HDBlup
  Ranef.HDBlup <-
    summary(HDBlup.model, coef = T)$coef.random[grep("H_full", rownames(summary(HDBlup.model, coef = T)$coef.random)),] %>%
    `rownames<-`(gsub("*.*_", "", rownames(.))) |>
    data.frame() |>
    rownames_to_column(var = "id") |>
    filter(id %in% Data$id) |> 
    # include PEV
    dplyr::mutate(
      PEV = std.error ^ 2,
      # Calculate accuracy
      raa = sqrt(1 - PEV / VarComp$HDBlup[str_detect(rownames(VarComp$HDBlup),
                                                     "id,"), "component"])
    )
  
  # Dominance effects
  ## From FS
  if (TypePop == "CP") {
    Dom.eff.FS <-
      summary(FS.AddDom.model, coef = T)$coef.random[grep("D.inv", rownames(summary(FS.AddDom.model, coef = T)$coef.random)),] %>%
      `rownames<-`(gsub("*.*_", "", rownames(.))) |>
      data.frame() |>
      rownames_to_column(var = "id") |>
      filter(id %in% Data$id)
  } else {
    Dom.eff.FS <- 0
  }
  
  
  ## From GDBlup
  Dom.eff.GD <-
    summary(GDBlup.model, coef = T)$coef.random[grep("GD", rownames(summary(GDBlup.model, coef = T)$coef.random)),] %>%
    `rownames<-`(gsub("*.*_", "", rownames(.))) |>
    data.frame() |>
    rownames_to_column(var = "id") |>
    filter(id %in% Data$id)
  
  ## From HDBlup
  Dom.eff.HD <-
    summary(HDBlup.model, coef = T)$coef.random[grep("HD", rownames(summary(HDBlup.model, coef = T)$coef.random)),] %>%
    `rownames<-`(gsub("*.*_", "", rownames(.))) |>
    data.frame() |>
    rownames_to_column(var = "id") |>
    filter(id %in% Data$id)
  
  if (TypePop == "OP") {
  BV <- list(
    BV.HS = Ranef.HS.Add |> 
      # Include Selfing and number of genotyped individuals
      dplyr::mutate(
        SelfRate = nP,
        `% genotyped` = Samp
      ),
    BV.FS = Ranef.FS.Add |> 
      # Include Selfing and number of genotyped individuals
      dplyr::mutate(
        SelfRate = nP,
        `% genotyped` = Samp
      ),
    BV.FS.Dom = Ranef.FS.AddDom |> 
      # Include Selfing and number of genotyped individuals
      dplyr::mutate(
        SelfRate = nP,
        `% genotyped` = Samp
      ),
    BV.GBlup = Ranef.GBlup |> 
      # Include Selfing and number of genotyped individuals
      dplyr::mutate(
        SelfRate = nP,
        `% genotyped` = Samp
      ),
    BV.GDBlup = Ranef.GDBlup |> 
      # Include Selfing and number of genotyped individuals
      dplyr::mutate(
        SelfRate = nP,
        `% genotyped` = Samp
      ),
    BV.HDBlup = Ranef.HDBlup |> 
      # Include Selfing and number of genotyped individuals
      dplyr::mutate(
        SelfRate = nP,
        `% genotyped` = Samp
      )
  )
  } else {
    BV <- list(
      BV.FS = Ranef.FS.Add |> 
        # Include Selfing and number of genotyped individuals
        dplyr::mutate(
          SelfRate = nP,
          `% genotyped` = Samp
        ),
      BV.FS.Dom = Ranef.FS.AddDom |> 
        # Include Selfing and number of genotyped individuals
        dplyr::mutate(
          SelfRate = nP,
          `% genotyped` = Samp
        ),
      BV.GBlup = Ranef.GBlup |> 
        # Include Selfing and number of genotyped individuals
        dplyr::mutate(
          SelfRate = nP,
          `% genotyped` = Samp
        ),
      BV.GDBlup = Ranef.GDBlup |> 
        # Include Selfing and number of genotyped individuals
        dplyr::mutate(
          SelfRate = nP,
          `% genotyped` = Samp
        ),
      BV.HDBlup = Ranef.HDBlup |> 
        # Include Selfing and number of genotyped individuals
        dplyr::mutate(
          SelfRate = nP,
          `% genotyped` = Samp
        )
    )
  }
  
  if (TypePop == "OP") {
    # Merge BV's
    BV_ALL <- Data |>
      dplyr::select(id, Pheno, BV, Dom_eff) |>
      dplyr::mutate(GV = BV + Dom_eff) |>
      left_join(BV$BV.HS |>
                  dplyr::select(id, solution) |>
                  rename(BV.HS = solution),
                by = c("id")) |>
      left_join(BV$BV.FS |>
                  dplyr::select(id, solution) |>
                  rename(BV.FS = solution),
                by = c("id")) |>
      left_join(BV$BV.FS.Dom |>
                  dplyr::select(id, solution) |>
                  rename(BV.FS.Dom = solution),
                by = c("id")) |>
      left_join(BV$BV.GBlup |>
                  dplyr::select(id, solution) |>
                  rename(BV.GBlup = solution),
                by = c("id")) |>
      left_join(BV$BV.GDBlup |>
                  dplyr::select(id, solution) |>
                  rename(BV.GDBlup = solution),
                by = c("id")) |>
      left_join(BV$BV.HDBlup |>
                  dplyr::select(id, solution) |>
                  rename(BV.HDBlup = solution),
                by = c("id")) |> 
      left_join(Dom.eff.GD |>
                  # Set up solution to NA
                  mutate(
                    solution = NA
                  ) |> 
                dplyr::select(id, solution) |>
                  rename(DomEff.FS = solution),
                by = c("id")) |>
      left_join(Dom.eff.GD |> 
                dplyr::select(id, solution) |>
                  rename(DomEff.GD = solution),
                by = c("id")) |>
      left_join(Dom.eff.HD |> 
                    dplyr::select(id, solution) |>
                    rename(DomEff.HD = solution),
                  by = c("id")) |> 
      dplyr::mutate(GV.FS = BV.FS + DomEff.FS,
                    GV.GD = BV.GDBlup + DomEff.GD,
                    GV.HD = BV.HDBlup + DomEff.HD) |> 
      # Include Selfing and number of genotyped individuals
      dplyr::mutate(
        SelfRate = nP,
        `% genotyped` = Samp
      )
  } else {
    BV_ALL <- Data |>
      dplyr::select(id, Pheno, BV, Dom_eff) |>
      dplyr::mutate(GV = BV + Dom_eff) |>
      left_join(BV$BV.FS |>
                  dplyr::select(id, solution) |>
                  rename(BV.FS = solution),
                by = "id") |>
      left_join(BV$BV.FS.Dom |>
                  dplyr::select(id, solution) |>
                  rename(BV.FS.Dom = solution),
                by = c("id")) |>
      left_join(BV$BV.GBlup |>
                  dplyr::select(id, solution) |>
                  rename(BV.GBlup = solution),
                by = c("id")) |>
      left_join(BV$BV.GDBlup |>
                  dplyr::select(id, solution) |>
                  rename(BV.GDBlup = solution),
                by = c("id")) |>
      left_join(BV$BV.HDBlup |>
                  dplyr::select(id, solution) |>
                  rename(BV.HDBlup = solution),
                by = c("id")) |>
      left_join(Dom.eff.FS |>
                  dplyr::select(id, solution) |>
                  rename(DomEff.FS = solution),
                by = c("id")) |>
      left_join(Dom.eff.GD |>
                  dplyr::select(id, solution) |>
                  rename(DomEff.GD = solution),
                by = c("id")) |>
      left_join(Dom.eff.HD |>
                  dplyr::select(id, solution) |>
                  rename(DomEff.HD = solution),
                by = c("id")) |>
      dplyr::mutate(GV.FS = BV.FS + DomEff.FS,
                    GV.GD = BV.GDBlup + DomEff.GD,
                    GV.HD = BV.HDBlup + DomEff.HD) |> 
      
      # Include Selfing and number of genotyped individuals
      dplyr::mutate(
        SelfRate = nP,
        `% genotyped` = Samp
      )
  }
  
  # Pull snps
  if (pullSNPs) {
    snps <- SNPs
  } else {
    snps <- NULL
  }
  
  # Extract 1+Fi from diagonal of Gmatrix
  DiagG <- diag(G) |>
    data.frame() |>
    `colnames<-`("1+Fi") |>
    rownames_to_column(var = "id")
  
  ## Calculate the mean inbreeding of population
  Mean.Inbreeding <- round(mean(DiagG$`1+Fi`, na.rm = T), 2) - 1
  
  # ------------------------------ Extract `ssGBLUP` ---------------------------------- #
  
  if (ssGBLUP) {
    # Get the data sets and matrix to be used in the ssGBLUP analysis
    if (TypePop == "OP") {
      A <- A.HS
      # D <- Simul[["D.FS"]]
    } else {
      A <- A.FS
    }
    Samp <<- Samp
    #Samp2 <- "0.10"
    SampleG <- c(sample(
      row.names(G),
      size = nrow(G) * as.numeric(Samp),
      replace = F
    )) |> 
      as.numeric() |> 
      sort() |> 
      as.character()
    G_Sample <- G[SampleG, SampleG]
    
    GD_Sample <- GD[SampleG, SampleG]
    # unique(c(D))
    # unique(c(A))
    
    #D_test <- D[SampleG, SampleG]
    # Tunning a Gmatrix
    G2A <-
      match.G2A(
        A = A,
        G = G_Sample,
        clean = TRUE,
        ord = TRUE,
        mism = TRUE,
        RMdiff = TRUE
      )
    # Tunning a GDmatrix
    G2D <-
      match.G2A(
        A = D.FS,
        G = GD_Sample,
        clean = TRUE,
        ord = TRUE,
        mism = TRUE,
        RMdiff = TRUE
      )
    ## Get tuned Gmatrix
    G_tuned <<- G.tuneup(G = G2A$Gclean,
                         A = G2A$Aclean,
                         #bend = T,
                         align = T)$Gb
    # G_tuned2 <<- G.tuneup(G = G2A$Gclean,
    #                      A = G2A$Aclean,
    #                      bend = T,
    #                      align = F)$Gb
    ## Get tuned GDmatrix
    GD_tuned <<- G.tuneup(G = G2D$Gclean,
                          A = G2D$Aclean,
                          #bend = T,
                          align = T)$Gb
    
    # Get G and GD inverse matrices
    Ginv <<- G.inverse(G = G_tuned, sparseform = FALSE)$Ginv
    #Ginv2 <<- G.inverse(G = G2A$Gclean, sparseform = FALSE)$Ginv
    GDinv <<- G.inverse(G = G2D$Gclean, sparseform = FALSE)$Ginv
    
    # Get H and HD inverse matrices
    ## H
    system.time(H <<-
                  H.inverse(
                    A = A,
                    G = Ginv,
                    
                    #tau = 2,
                    #omega = 0,
                    lambda = 0.9,
                    sparseform = TRUE
                  ))
    # ~22s
    
    # system.time(H2 <<-
    #               H.inverse(
    #                 A = A,
    #                 G = Ginv2,
    #                 
    #                 #tau = 2,
    #                 #omega = 0,
    #                 lambda = 0.9,
    #                 sparseform = TRUE
    #               ))
    ## HD
    system.time(HD <<-
                  H.inverse(
                    A = D.FS,
                    G = GDinv,
                    #tau = 2,
                    #omega = 0.1,
                    lambda = 0.9,  
                    sparseform = TRUE
                  ))
    # ~22s
    
    
    # Set up single tree plot models
    if (n_Tree == 1) {
      # Fit ssGBlup
      system.time(
        ssGBlup <- try({asreml(
          Pheno ~ 1 + Block,
          random = ~ vm(id, H),
          data = Data |>
            dplyr::mutate(
              Block = as.factor(Block),
              Plot = as.factor(Plot),
              id = as.factor(id),
              id2 = as.factor(id)
            ),
          workspace = 100e05
        )})
      )
      # ~ 7s
      
      # Fit ssGDBlup
      system.time(
        ssGDBlup <- try({asreml(
          Pheno ~ 1 + Block,
          random = ~ vm(id, H) + vm(id2, HD),
          data = Data |>
            dplyr::mutate(
              Block = as.factor(Block),
              Plot = as.factor(Plot),
              id = as.factor(id),
              id2 = as.factor(id)
            ),
          workspace = 200e05
        )})
      )
      # ~ 50s
      
      # Set up non single tree plot models
    } else {
      # Fit ssGBlup
      system.time(
        ssGBlup <- try({asreml(
          Pheno ~ 1 + Block,
          random = ~ Plot + vm(id, H),
          data = Data |>
            dplyr::mutate(
              Block = as.factor(Block),
              Plot = as.factor(Plot),
              id = as.factor(id),
              id2 = as.factor(id)
            ),
          workspace = 100e05
        )})
      )
      # ~ 4s
      
      # Fit ssGDBlup
      system.time(
        ssGDBlup <- try({asreml(
          Pheno ~ 1 + Block,
          random = ~ Plot + vm(id, H2) + vm(id2, HD),
          data = Data |>
            dplyr::mutate(
              Block = as.factor(Block),
              Plot = as.factor(Plot),
              id = as.factor(id),
              id2 = as.factor(id)
            ),
          workspace = 200e05
        )})
      )
      # ~ 118.52s
    }
    
    # Extract the variance 
    if (class(ssGBlup) == "try-error") {
      ssG.varcomp <- NA
    } else {
      #ssGBlup <- update.asreml(ssGBlup)
      ssG.varcomp <- summary(ssGBlup)$varcomp 
    }
    
    if (class(ssGDBlup) == "try-error") {
      ssGD.varcomp <- NA
    } else {
      #ssGDBlup <- update.asreml(ssGDBlup)
      ssGD.varcomp <- summary(ssGDBlup)$varcomp
    }
    
    VarComp.ssGBLUP <- list(ssGBlup = ssG.varcomp,
                            ssGDBlup = ssGD.varcomp)
    
    # IF single tree plot
    if (n_Tree == 1) {
      # Heritabilities
      ## h2a
      if (class(ssGBlup) == "try-error") {
        h2a.ssGBlup <- data.frame(
          h2.type = "h2a",
          Estimate = NA,
          SE = NA,
          TypePop = "ssGBlup"
        )
      } else {
        h2a.ssGBlup <- pin(ssGBlup, h2a ~ V1 / (V1 + V2)) |>
          dplyr::mutate(TypePop = "ssGBlup") |>
          rownames_to_column(var = "h2.type") |> 
          # Include Selfing and number of genotyped individuals
          dplyr::mutate(
            SelfRate = nP,
            `% genotyped` = Samp
          )
      }
      
      if (class(ssGDBlup) == "try-error") {
        h2a.ssGDBlup <- data.frame(
          h2.type = "h2a",
          Estimate = NA,
          SE = NA,
          TypePop = "ssGDBlup"
        )
      } else {
        h2a.ssGDBlup <- pin(ssGDBlup, h2a ~ V1 / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "ssGDBlup") |>
          rownames_to_column(var = "h2.type") |> 
          # Include Selfing and number of genotyped individuals
          dplyr::mutate(
            SelfRate = nP,
            `% genotyped` = Samp
          )
      }
      
      ## h2d
      if (class(ssGDBlup) == "try-error") {
        h2d.ssGDBlup <- data.frame(
          h2.type = "h2d",
          Estimate = NA,
          SE = NA,
          TypePop = "ssGDBlup"
        )
      } else {
        h2d.ssGDBlup <- pin(ssGDBlup, h2d ~ V2 / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "ssGDBlup") |>
          rownames_to_column(var = "h2.type") |> 
          # Include Selfing and number of genotyped individuals
          dplyr::mutate(
            SelfRate = nP,
            `% genotyped` = Samp
          )
      }
      ## h2g
      if (class(ssGDBlup) == "try-error") {
        h2g.ssGDBlup <- data.frame(
          h2.type = "h2g",
          Estimate = NA,
          SE = NA,
          TypePop = "ssGDBlup"
        )
      } else {
        h2g.ssGDBlup <- pin(ssGDBlup, h2g ~ (V1 + V2) / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "ssGDBlup") |>
          rownames_to_column(var = "h2.type") |> 
          # Include Selfing and number of genotyped individuals
          dplyr::mutate(
            SelfRate = nP,
            `% genotyped` = Samp
          )
      }
      # If not a single tree plot
    } else {
      # Heritabilities
      ## h2a
      if (class(ssGBlup) == "try-error") {
        h2a.ssGBlup <- data.frame(
          h2.type = "h2a",
          Estimate = NA,
          SE = NA,
          TypePop = "ssGBlup"
        )
      } else {
        h2a.ssGBlup <- pin(ssGBlup, h2a ~ V2 / (V1 + V2 + V3)) |>
          dplyr::mutate(TypePop = "ssGBlup") |>
          rownames_to_column(var = "h2.type")
      }
      
      if (class(ssGDBlup) == "try-error") {
        h2a.ssGDBlup <- data.frame(
          h2.type = "h2a",
          Estimate = NA,
          SE = NA,
          TypePop = "ssGDBlup"
        )
      } else {
        h2a.ssGDBlup <- pin(ssGDBlup, h2a ~ V2 / (V1 + V2 + V3 + V4)) |>
          dplyr::mutate(TypePop = "ssGDBlup") |>
          rownames_to_column(var = "h2.type")
      }
      #names(ssGBLUP.CP[["S1"]][["0.01"]][["0.25"]][["Herit"]])
      ## h2d
      if (class(ssGDBlup) == "try-error") {
        h2d.ssGDBlup <- data.frame(
          h2.type = "h2d",
          Estimate = NA,
          SE = NA,
          TypePop = "ssGDBlup"
        )
      } else {
        h2d.ssGDBlup <- pin(ssGDBlup, h2d ~ V3 / (V1 + V2 + V3 + V4)) |>
          dplyr::mutate(TypePop = "ssGDBlup") |>
          rownames_to_column(var = "h2.type")
      }
      ## h2g
      if (class(ssGDBlup) == "try-error") {
        h2g.ssGDBlup <- data.frame(
          h2.type = "h2g",
          Estimate = NA,
          SE = NA,
          TypePop = "ssGDBlup"
        )
      } else {
        h2g.ssGDBlup <- pin(ssGDBlup, h2g ~ (V2 + V3) / (V1 + V2 + V3 + V4)) |>
          dplyr::mutate(TypePop = "ssGDBlup") |>
          rownames_to_column(var = "h2.type")
      }
    }
    
    
    Herit.ssGBlup <- rbind(
      # h2a
      h2a.ssGBlup,
      h2a.ssGDBlup,
      # h2d
      h2d.ssGDBlup,
      # h2g
      h2g.ssGDBlup
    ) |>
      dplyr::mutate(Estimate = round(Estimate, 2),
                    SE = round(SE, 3)) |> 
      # Include Selfing and number of genotyped individuals
      dplyr::mutate(
        SelfRate = nP,
        `% genotyped` = Samp
      )
    
    
    # Extract Blups
    # BLUPs from ssGBlup
    if (class(ssGBlup) == "try-error") {
      Ranef.ssGBlup <- data.frame(
        id = Data$id,
        solution = NA,
        std.error = NA,
        z.ratio = NA,
        PEV = NA,
        raa = NA
      )
    } else {
      Ranef.ssGBlup <-
        summary(ssGBlup, coef = T)$coef.random[grep("vm", rownames(summary(ssGBlup, coef = T)$coef.random)), ] %>%
        `rownames<-`(gsub("*.*_", "", rownames(.))) |>
        data.frame() |>
        rownames_to_column(var = "id") |>
        filter(id %in% Data$id) |> 
        # include PEV
        dplyr::mutate(
          PEV = std.error ^ 2,
          # Calculate accuracy
          raa = sqrt(1 - PEV / VarComp.ssGBLUP$ssGBlup[str_detect(rownames(VarComp.ssGBLUP$ssGBlup),
                                                                  "id,"), "component"])
        )
    }
    
    # BLUPs from ssGDBlup
    if (class(ssGDBlup) == "try-error") {
      Ranef.ssGDBlup <- data.frame(
        id = Data$id,
        solution = NA,
        std.error = NA,
        z.ratio = NA,
        PEV = NA,
        raa = NA
      )
    } else {
      Ranef.ssGDBlup <-
        summary(ssGDBlup, coef = T)$coef.random[grep("id, H", rownames(summary(ssGDBlup, coef = T)$coef.random)), ] %>%
        `rownames<-`(gsub("*.*_", "", rownames(.))) |>
        data.frame() |>
        rownames_to_column(var = "id") |>
        filter(id %in% Data$id) |> 
        # include PEV
        dplyr::mutate(
          PEV = std.error ^ 2,
          # Calculate accuracy
          raa = sqrt(1 - PEV / VarComp.ssGBLUP$ssGDBlup[str_detect(rownames(VarComp.ssGBLUP$ssGDBlup),
                                                                   "id,"), "component"])
        )
    }
    
    # Dominance effects From ssGDBlup
    if (class(ssGDBlup) == "try-error") {
      Dom.eff.GD.ssGDBLUP <- data.frame(
        id = Data$id,
        solution = NA,
        std.error = NA,
        z.ratio = NA,
        PEV = NA,
        raa = NA
      ) 
    } else {
      Dom.eff.GD.ssGDBLUP <-
        summary(ssGDBlup, coef = T)$coef.random[grep("id2, HD", rownames(summary(ssGDBlup, coef = T)$coef.random)), ] %>%
        `rownames<-`(gsub("*.*_", "", rownames(.))) |>
        data.frame() |>
        rownames_to_column(var = "id") |>
        filter(id %in% Data$id) |> 
        # include PEV
        dplyr::mutate(
          PEV = NA,
          # Calculate accuracy
          raa = NA
        )
    }
    
    BV.ssGBLUP <- list(
      BV.ssGBlup = Ranef.ssGBlup |> 
        # Include Selfing and number of genotyped individuals
        dplyr::mutate(
          SelfRate = nP,
          `% genotyped` = Samp
        ),
      BV.ssGDBlup = Ranef.ssGDBlup |> 
        # Include Selfing and number of genotyped individuals
        dplyr::mutate(
          SelfRate = nP,
          `% genotyped` = Samp
        ),
      Dom.ssGDBlup = Dom.eff.GD.ssGDBLUP |> 
        # Include Selfing and number of genotyped individuals
        dplyr::mutate(
          SelfRate = nP,
          `% genotyped` = Samp
        )
    )
    # Merge BV's
    BV_ALL.ssGBLUP <- Data |>
      dplyr::select(id, Pheno, BV, Dom_eff, GV) |>
      dplyr::mutate(GV_true = GV,
                    GV = BV + Dom_eff) |>
      left_join(DiagG,
                by = "id") |>
      left_join(BV.ssGBLUP$BV.ssGBlup |>
                  dplyr::select(id, solution) |>
                  rename(BV.ssGBlup = solution),
                by = "id") |>
      left_join(BV.ssGBLUP$BV.ssGDBlup |>
                  dplyr::select(id, solution) |>
                  rename(BV.ssGDBlup = solution),
                by = "id") |>
      left_join(
        BV.ssGBLUP$Dom.ssGDBlup |>
          dplyr::select(id, solution) |>
          rename(Dom.ssGDBlup = solution),
        by = "id"
      ) |>
      dplyr::mutate(GV.ssGDBlup = BV.ssGDBlup + Dom.ssGDBlup) |> 
      # Include Selfing and number of genotyped individuals
      dplyr::mutate(
        SelfRate = nP,
        `% genotyped` = Samp
      )
    
    # Merge Herit
    Herit.All <- full_join(
      Herit,
      Herit.ssGBlup,
      by = c(
        "h2.type", 
        "Estimate", 
        "SE", 
        "TypePop", 
        "SelfRate", 
        "% genotyped"
      )
    ) |> 
      # Arrange
      arrange(
        h2.type
      )
    Results <- list(
      TypePop = TypePop,
      `% of genotyped individuals` = scales::percent(as.numeric(Samp)),
      # Results Pedigree and SNP based models
      SelfRate = Inbreeding,
      RealParameters = GenVars,
      VarComp = VarComp,
      Herit.All = Herit.All,
      BV = BV,
      BV_ALL = BV_ALL,
      # Results ssGBLUP models
      Mean.Inbreeding = Mean.Inbreeding,
      VarComp.ssGBLUP = VarComp.ssGBLUP,
      BV.ssGBLUP = BV.ssGBLUP,
      BV_ALL.ssGBLUP = BV_ALL.ssGBLUP,
      Data = Data,
      # 1 + Fi
      DiagG = DiagG,
      snps = SNPs,
      haplotype = haplotypes
    )
  } else {
    Results <- list(
      TypePop = TypePop,
      `% of genotyped individuals` = scales::percent(as.numeric(Samp)),
      # Results Pedigree and SNP based models
      SelfRate = Inbreeding,
      RealParameters = GenVars,
      VarComp = VarComp,
      Herit.All = Herit,
      BV = BV,
      BV_ALL = BV_ALL,
      # Results ssGBLUP models
      Mean.Inbreeding = Mean.Inbreeding,
      Data = Data,
      # 1 + Fi
      DiagG = DiagG,
      snps = SNPs,
      haplotype = haplotypes
    )
  }
  
  
  
  }
  
  return(Results)
}

# test_dom0.2_GD_vitezica <- Results
# test_dom0.2_GD_tuned_vitezica <- Results
# test_dom0.2_GD_tuned_Su <- Results
# test_dom0.2_GD_Su <- Results
# In these tests, the "GD" matrix from the "Su" Method outperforms the "Vitezica" method and should be used in the analysis.
# test_dom0.2_G.Yang_GD.Su <- Results
# test_dom0.2_G.VanRaden_GD.Su <- Results
# test_dom0.2_G.VanRaden_GD.Su_nP0.5 <- Results
# test_dom0.2_G.Yang_GD.Su_nP0.5 <- Results
# 
# test_dom0.3_G.Vanraden_GD.Su <- Results
# test_dom0.3_G.Vanraden_GD.Vitezica <- Results
