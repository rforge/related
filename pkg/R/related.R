#--------------------------#
# Function - coancestry
#--------------------------#
# This function implements Jinlinag Wang's Fortran code for the
# Coancestry program.  It provides seven estimators of relatedness
# based on codominant molecular data, and can accomodate
# inbreeding and genotyping errors.
# *****  IT REQUIRES *****
# A data set on which to calculate relatedness, and specification
# of options.
#--------------------------#

coancestry <- function(genotype.data, error.rates = 0, allele.freqs = NULL, 
                    trioml = 0L, wang = 0L, lynchli = 0L, lynchrd = 0L,
                    ritland = 0L, quellergt = 0L, dyadml = 0L,
                    ci95.num.bootstrap = 100L, trioml.num.reference = 100L,
                    allow.inbreeding = FALSE, rng.seed = NULL,
                    working.directory = tempdir(),
                    output.file = FALSE) {

  # Constants consistent with the requirements of Related's Fortran source.
  kMaxBoots         <- 100000L
  kMaxRefInds       <- 10000L

  # Filenames that the Related binary reads data from depend upon
  # the prefix of long names given for its output file.
  temp.filename <- 'output.txt'
  file.prefix <- substr(temp.filename, 1, nchar(output.file)-10)
  genotype.filename <- paste0(file.prefix, 'GenotypeData.Txt')
  freqs.filename    <- paste0(file.prefix, 'FrequencyData.Txt')
  config.filename   <- paste0(file.prefix, 'TrioR.DAT')

  # These parameters are simpler to handle as a vector.
  # Order in this vector matters.
  estimates <- c(trioml, wang, lynchli, lynchrd, ritland, quellergt, dyadml)

  # An absolute path makes it easier to work with output.file.
  if (is.character(output.file)) {
    output.file <- file_path_as_absolute(output.file)
  }

  # Do some basic input validation.
  storage.mode(error.rates)          <- 'double'
  storage.mode(estimates)            <- 'integer'
  storage.mode(ci95.num.bootstrap)   <- 'integer'
  storage.mode(trioml.num.reference) <- 'integer'
  if (!is.null(rng.seed)) storage.mode(rng.seed) <- 'integer'
  else rng.seed <- sample.int(.Machine$integer.max, 1)

  if (any(estimates < 0) || any(estimates > 2))
    stop(paste('Set an estimator to 0, 1 or 2 to skip, perform,',
               'or perform with CI95 via bootstrapping, respectively.'))

  if (!any(estimates > 0)) stop('No estimators selected.')

  if (any(2 == estimates)
  && (ci95.num.bootstrap < 1 || ci95.num.bootstrap > kMaxBoots))
    stop(sprintf('The number of bootstrapping samples is outside [1,%d].',
                 kMaxBoots))
  else if (!any(2 == estimates))  # just leave a safe default value
    ci95.num.bootstrap <- 100L

  if (trioml > 0
  && (trioml.num.reference < 1 || trioml.num.reference > kMaxRefInds))
    stop(paste(sprintf('The number of reference individuals for the',
                       'Triadic ML estimate is outside [1,%d].',
                       kMaxRefInds)))
  else if (0 == trioml)  # just leave a safe default value
    trioml.num.reference <- 100L

  if (!is.logical(allow.inbreeding))
    stop('allow.inbreeding should be of type "logical".')

  # genotype.data may be a filename or a data.frame.
  # column names are currently ignored in a data.frame, and should
  # be arranged
  # ID, locus_1, ..., locus_n
  # where each locus contains a pair of integers.
  # If a filename is provided, the file should formatted according to
  # the coancestry User Guide, Section 3.2.2.1. The file is co-located
  # with the working directory of the Related binary when it is
  # invoked.
  ret <- readgenotypedata(genotype.data)
  gdata    <- ret$gdata
  nloci    <- ret$nloci
  nalleles <- ret$nalleles
  ninds    <- ret$ninds
  freqs    <- ret$freqs

  # Validate allele.freqs now that we know nloci and nalleles.
  # This user-provided allele.freqs may override the calculates freqs.
  # nalleles may also be altered by this code block.
  if (!is.null(allele.freqs)) {
    if (!is.list(allele.freqs)
    ||  length(allele.freqs) != nloci)
      stop('allele.freqs should be a list of nloci data.frames')

    # Each data.frame should have at least as many rows as nalleles[i].
    for (i in 1:nloci) {
      if (!is.data.frame(allele.freqs[[i]])
      ||  length(allele.freqs[[i]]) != 2
      ||  !is.integer(allele.freqs[[i]][,1])
      ||  !is.numeric(allele.freqs[[i]][,2]))
        stop(paste('each data.frame in allele.freqs should have',
                   'integer alleles and numeric frequencies'))

      gdata.i <- c(gdata[[2*i]], gdata[[2*i+1]])
      gdata.i <- unique(gdata.i[gdata.i > 0])

      if (length(allele.freqs[[i]][[1]]) < nalleles[i]
      ||  length(setdiff(gdata.i, allele.freqs[[i]][,1])) != 0)  # subset test
        stop('allele.freqs must have a frequency for each allele in data')
      if (length(allele.freqs[[i]][[1]])
      !=  length(unique(allele.freqs[[i]][[1]])))
        stop('cannot have two frequencies for one allele')
      if (any(allele.freqs[[i]][,2] < 0) || any(1 < allele.freqs[[i]][,2]))
        stop('frequencies in allele.freqs should be in [0,1]')

      # It is possible to have frequencies for more alleles than appear
      # in the sample.
      nalleles[i] <- max(nalleles[i], length(allele.freqs[[i]][[1]]))

      # Normalize the frequencies.
      allele.freqs[[i]][,2] <- allele.freqs[[i]][,2] / sum(allele.freqs[[i]][,2])
    }

    # User-provided allele.freqs should be valid; give them to Related.
    freqs <- allele.freqs
  }

  # Write out the data.frame the Related binary's working directory.
  # This is done even when the user provides genotype.data as a filename.
  gdatafile <- paste0(working.directory, '/', genotype.filename)
  write.table(gdata, gdatafile, quote = F, na = '0', row.names = F,
              col.names = F)

  # Also write out the per-locus allele frequencies to the Related
  # binary's working directory.
  freqs.lines <- vector(mode = 'character', length = 2*nloci)
  for (i in 1:nloci) {
    freqs.lines[2*i-1] <- paste(freqs[[i]][,1], collapse = ' ')
    freqs.lines[ 2*i ] <- paste(freqs[[i]][,2], collapse = ' ')
  }
  freqsfile <- paste0(working.directory, '/', freqs.filename)
  writeLines(freqs.lines, freqsfile)

  # error.rates are the mistyping rates of alleles at each locus.
  # These form line 10 of the configuration file.
  if (1 == length(error.rates))
    error.rates <- rep(error.rates, nloci)
  else if (nloci != length(error.rates))
    stop('error.rates must be either a scalar or a vector of length nloci.')
  if (any(error.rates < 0) || any(1 < error.rates))
    stop('error.rates must be between 0 and 1.')

  # Build the coancestry configuration file, which contains twelve lines.
  # Please see the coancestry User Guide for more information.
  config.lines <- vector(mode = 'character', length = 12)
  # 01: name of output file
  config.lines[1] <- temp.filename
  # 02: number of loci considered
  config.lines[2] <- nloci
  # 03: number of alleles observed, per locus
  config.lines[3] <- paste(nalleles, collapse = ' ')
  # 04: 1 if allele frequencies are available in FrequencyData.Txt, else 0
  config.lines[4] <- 1L
  # 05: 1 if inbreeding is NOT allowed, 0 if it is allowed
  config.lines[5] <- ifelse(allow.inbreeding, 1, 0)
  # 06: seed for Related's pseudo random number generator
  config.lines[6] <- rng.seed
  # 07: number of bootstrapping samples for CI95 estimates
  config.lines[7] <- ci95.num.bootstrap
  # 08: number of reference individuals for the Triadic ML estimator
  config.lines[8] <- trioml.num.reference
  # 09: number of individuals in the dataset
  config.lines[9] <- ninds
  # 10: allele mistyping rate, per locus
  config.lines[10] <- paste(error.rates, collapse = ' ')
  # 11: 1 to estimate r for all dyads - 0 is buggy 
  config.lines[11] <- 1L
  # 12: 0,1,2 for each of the seven estimator selections
  config.lines[12] <- paste(estimates, collapse = ' ')

  # Write the configuration file to where the Fortran code will expect it.
  cfgfile <- paste0(working.directory, '/', config.filename)
  writeLines(config.lines, cfgfile)

  # Invoke the Related code using .Fortran.
  old.wd <- getwd()
  setwd(working.directory)
  tryCatch({
    print(system.time(.Fortran('related', PACKAGE = 'related')))
  }, interrupt = function(ex) {
    # If the user interrupts the Related code, memory will be leaked.
    # Clean it up here.
    .Fortran('clean_mem', PACKAGE = 'related')
    setwd(old.wd)
    stop('coancestry execution interrupted')
  })

  # Read files output by the Related code into a list of data.frames.
  cat('\nReading output files into data.frames... ')
  output <- readrelatedoutput(working.directory, file.prefix,
                              any(2 == estimates), allow.inbreeding)
  cat('Done!\n')

  # Insert the frequency data into the output list.
  output[['freqs']] <- freqs

  # If the output.file parameter is set, copy the Related code's output
  # file to that location.
  if (is.character(output.file))
    file.copy(temp.filename, output.file)

  # Reset R's working directory.
  setwd(old.wd)

  # Return the list of output data.frame objects to the user.
  return(output)
}

Interleave <- function(v1, v2) {
  m <- matrix(c(v1, v2), ncol = 2)
  return(as.vector(t(m)))
}

#--------------------------#
# Function - readgenotypedata
#--------------------------#
# This function will read a genotype data file into R
# and format as a dataframe as needed for the coancestry
# function.  It will also format the allele frequency data
# in a way that can be used by the simulation functions.
# *****  IT REQUIRES *****
# A file containing genotype data that is to be read
# into R for analyses
#--------------------------#
readgenotypedata <- function(genotype.data) {

  # Constants consistent with the requirements of Related's Fortran source.
  kMaxAlleles <- 127L

  # If genotype.data is given as a filename, read it to
  # a) check it's formatting
  # b) count nloci, nalleles and ninds
  if (is.character(genotype.data))
    genotype.data <- read.table(genotype.data, header = F, strip.white = T,
                                comment.char = '', check.names = F,
                                colClasses = c(V1='character'))
  else if (!is.data.frame(genotype.data))
    stop('genotype.data must be either a data.frame or a string.')

  # The first column of the data.frame holds IDs for the sampled individuals.
  storage.mode(genotype.data[[1]]) <- 'character'
  # Remaining columns contain integers for alleles, and must come in pairs
  # (one pair per locus).
  if (ncol(genotype.data)%%2L != 1L)
    stop('genotype.data must have two columns per locus for alleles.')
  nloci <- as.integer((ncol(genotype.data)-1)/2)
  if (nloci < 1) stop('There data at at least one locus.')

  # Count the number of unique alleles at each locus and determine
  # allele frequencies.
  freqs <- list()
  nalleles <- vector(mode = 'integer', length = nloci)
  for (i in 1:nloci) {
    storage.mode(genotype.data[[ 2*i ]]) <- 'integer'
    storage.mode(genotype.data[[2*i+1]]) <- 'integer'

    v <- c(genotype.data[[2*i]], genotype.data[[2*i+1]])
    v <- v[v > 0]  # 0 means unknown allele. Discount it.
    nalleles[i] <- length(unique(v))

    if (nalleles[i] > kMaxAlleles)
      stop(sprintf(
        'Number of distinct alleles at locus %d is too large: %d > %d',
        i, nalleles[i], kMaxAlleles))

    # Use table() to calculate allele frequencies and then load the
    # results into a more convenient data.frame.
    freq_table <- table(v, dnn = 'allele_freqs') / length(v)
    freq_df <- data.frame(allele = integer(nalleles[i]),
                          frequency = numeric(nalleles[i]))
    freq_df$allele <- as.integer(paste(attr(freq_table,'dimnames')[[1]],
                                        sep = ' '))
    freq_df$frequency <- as.numeric(freq_table)
    freqs[[paste0('locus',i)]] <- freq_df
  }

  ninds <- nrow(genotype.data)
  if (ninds < 1) stop('The sample size should be positive.')

  return(list(gdata = genotype.data,
              nloci = nloci,
              nalleles = nalleles,
              ninds = ninds,
              freqs = freqs))
}

#--------------------------#
# Function - readrelatedoutput
#--------------------------#
# This function reads the output from the coancestry
# code into appropriate R data frames.  This 
# function is not seen by the user, but is used within
# the functions provided.
# *****  IT REQUIRES *****
# Output from the coancestry code.
#--------------------------#
readrelatedoutput <- function(working.directory = tempdir(),
                                 file.prefix = '', any.ci95 = F,
                                 allow.inbreeding = F) {

  # Filenames that the Related binary writes data to depend upon
  # the prefix of long names given for its main output file.
  out.prefix <- paste0(working.directory, '/', file.prefix)
  re.filename         <- paste0(out.prefix, 'RelatednessEstimates.Txt')
  re.ci95.filename    <- paste0(out.prefix, 'RelatednessCI95.Txt')
  d7.filename         <- paste0(out.prefix, 'Delta7Estimates.Txt')
  d7.ci95.filename    <- paste0(out.prefix, 'Delta7CI95.Txt')
  d8.filename         <- paste0(out.prefix, 'Delta8Estimates.Txt')
  d8.ci95.filename    <- paste0(out.prefix, 'Delta8CI95.Txt')
  inbrd.filename      <- paste0(out.prefix, 'InbreedingEstimates.Txt')
  inbrd.ci95.filename <- paste0(out.prefix, 'InbreedingCI95.Txt')

  # Read files output by the Related binary into several data.frame objects.

  # Point estimate of relatedness, r
  enames      <- c('trioml', 'wang', 'lynchli', 'lynchrd',
                   'ritland', 'quellergt', 'dyadml')
  cnames      <- c('pair.no', 'ind1.id', 'ind2.id', 'group', enames)
  cclasses    <- c('integer', rep('character', 3),
                   rep('double', length(enames)))
  relatedness <- read.csv(re.filename, header = F, strip.white = T,
                          col.names = cnames, colClasses = cclasses)

  if (any.ci95) {
    # CI-95 estimate of relatedness, r
    cnames           <- c('pair.no', 'ind1.id', 'ind2.id', 'group',
                          Interleave(paste0(enames, '-low'),
                                     paste0(enames, '-high')))
    cclasses         <- c('integer', rep('character',3),
                          rep('double',length(enames)*2))
    relatedness.ci95 <- read.csv(re.ci95.filename, header = F, strip.white = T,
                                 col.names = cnames, colClasses = cclasses)
  }

  # Point estimate of delta7
  enames   <- c('trioml', 'wang', 'lynchrd', 'dyadml')
  cnames   <- c('pair.no', 'ind1.id', 'ind2.id', 'group', enames)
  cclasses <- c('integer', rep('character', 3), rep('double', length(enames)))
  delta7   <- read.csv(d7.filename, header = F, strip.white = T,
                       col.names = cnames, colClasses = cclasses)

  # Point estimate of delta8
  delta8 <- read.csv(d8.filename, header = F, strip.white = T,
                     col.names = cnames, colClasses = cclasses)

  if (any.ci95) {
    # CI-95 estimate of delta7
    cnames      <- c('pair.no', 'ind1.id', 'ind2.id', 'group',
                     Interleave(paste0(enames, '-low'),
                                paste0(enames, '-high')))
    cclasses    <- c('integer', rep('character', 3),
                     rep('double', length(enames)*2))
    delta7.ci95 <- read.csv(d7.ci95.filename, header = F, strip.white = T,
                            col.names = cnames, colClasses = cclasses)

    # CI-95 estimate of delta8
    delta8.ci95 <- read.csv(d8.ci95.filename, header = F, strip.white = T,
                            col.names = cnames, colClasses = cclasses)
  }

  # Point estimate of inbreeding coefficient, F
  if (allow.inbreeding)
    enames   <- c('LH', 'LR')
  else
    enames   <- c('L3', 'LH', 'LR', 'L2')
  cnames     <- c('ind.id', enames)
  cclasses   <- c('character', rep('double', length(enames)))
  inbreeding <- read.table(inbrd.filename, header = F, strip.white = T,
                           comment.char = '',
                           col.names = cnames, colClasses = cclasses)

  if (any.ci95) {
  # CI-95 estimate of inbreeding coefficient, F
    cnames          <- c('ind.id', Interleave(paste0(enames, '-low'),
                                              paste0(enames, '-high')))
    cclasses        <- c('integer', rep('double', length(enames)*2))
    inbreeding.ci95 <- read.table(inbrd.ci95.filename, header = F,
                                  strip.white = T, comment.char = '',
                                  col.names = cnames, colClasses = cclasses)
  }

  # Return a list of data.frame objects to the user.
  if (any.ci95)
    return(list(relatedness      = relatedness,
                relatedness.ci95 = relatedness.ci95,
                delta7           = delta7,
                delta7.ci95      = delta7.ci95,
                delta8           = delta8,
                delta8.ci95      = delta8.ci95,
                inbreeding       = inbreeding,
                inbreeding.ci95  = inbreeding.ci95))
  else
    return(list(relatedness      = relatedness,
                delta7           = delta7,
                delta8           = delta8,
                inbreeding       = inbreeding))
}


#--------------------------#
# Function - familysim
#--------------------------#
# This function will generate a user-defined number of pairs of 
# individuals of known relatedness (parent-offspring, full-sib,
# half-sib, and unrelated), that can be used to test the performance
# of different estimators, or to assess the expected resolution for
# a given data set.
# *****  IT REQUIRES *****
# Allele frequency data, formatted as by the readgenotypedata function, and
# an integer for the number of pairs to simulate at each relatedness value.
#--------------------------#
familysim <- function(freqs, ninds = 100L) {

  # Extract information from freqs, which is a named list of named tables,
  # to prepare for the .External call.

  storage.mode(ninds) <- 'integer'
  nloci <- length(freqs)

  nalleles <- vector(length = nloci, mode = 'integer')
  for (i in 1:nloci) nalleles[i] <- length(freqs[[i]]$allele)

  alleles     <- vector(length = sum(nalleles), mode = 'integer')
  frequencies <- vector(length = sum(nalleles), mode = 'numeric')
  start <- 0L
  for (i in 1:nloci) {
    alleles[(start+1):(start+nalleles[i])]     <- freqs[[i]]$allele
    frequencies[(start+1):(start+nalleles[i])] <- freqs[[i]]$frequency
    start <- start + nalleles[i]
  }

  # Use family_sim to generate an integer matrix of simulated genotype data.
  sim_matrix <- .External("family_sim", ninds, nloci, nalleles,
                          alleles, frequencies)

  # make labels for individuals (required by Related).
  # these labels reflect relationship categories:
  #   1st ninds pairs: simulated parent-offspring
  #   2nd ninds pairs: simulated full siblings
  #   3rd ninds pairs: simulated half siblings
  #   4th ninds pairs: simulated unrelated individuals
  IDs <- vector(mode = 'character', length = 8*ninds)
  for (i in 1:ninds) {
    IDs[2*i-1] <- sprintf('PO%.7dp',i)
    IDs[ 2*i ] <- sprintf('PO%.7do',i)
  }
  for (i in (ninds+1):(2*ninds)) {
    IDs[2*i-1] <- sprintf('SB%.7da',i-ninds)
    IDs[ 2*i ] <- sprintf('SB%.7db',i-ninds)
  }
  for (i in (2*ninds+1):(3*ninds)) {
    IDs[2*i-1] <- sprintf('HS%.7da',i-2*ninds)
    IDs[ 2*i ] <- sprintf('HS%.7db',i-2*ninds)
  }
  for (i in (3*ninds+1):(4*ninds)) {
    IDs[2*i-1] <- sprintf('UR%.7da',i-3*ninds)
    IDs[ 2*i ] <- sprintf('UR%.7db',i-3*ninds)
  }

  # Convert the integer matrix into a more convenient data.frame.
  sim <- data.frame(IDs, sim_matrix, stringsAsFactors = F)
  colnames(sim)[1] <- 'ID'
  for (i in 1:nloci) {
    colnames(sim)[ 2*i ] <- paste0('locus',i,'a')
    colnames(sim)[2*i+1] <- paste0('locus',i,'b')
  }

  return(sim)
}

#--------------------------#
# Function - cleanuprvals
#--------------------------#
# This function removes the unnecessary pairwise
# relatedness values from a data frame generated by
# the Related function on a simulated data set.
# *****  IT REQUIRES *****
# 1. The name of the data frame to be cleaned up
# 2. The number of individuals that were simulated for each type of 
#    relatedness (same number for all relatedness values).
#------------------------------------------------------------#

cleanuprvals <- function(simdata, ninds) {
    skip <- ((ninds * 4) * 4) - 3
    desired <- 1
    pair <- 1
    simcleaned <- simdata[1, ]
    while (skip > 0) {
        simcleaned[pair, ] <- simdata[desired, ]    
        desired <- desired + skip
        skip <- skip - 4
        pair <- pair + 1
    }
    return(simcleaned)
}
#--------------------------#
# Function - compareestimators
#--------------------------#
# This function analyses simulated data with all moment-based relatedness
# estimators and then plots the data in way where estimates can be easily 
# compared.
# It takes 2 arguments:
# 1. A data frame of genotype data
# 2. The number of individuals to simulate
#------------------------------------------------------------#

compareestimators <- function(filename, ninds) {
	#------------------------------------------------#
	# Generate simulated data
	#------------------------------------------------#
	simdata <- familysim(filename$freqs, ninds)
	output <- coancestry(simdata, lynchli=1, lynchrd=1, quellergt=1, wang=1)
	simrel <- cleanuprvals(output$relatedness, ninds)
	
	#------------------------------------------------#
	# Parse out data based on relatedness type and
	# estimator used
	#------------------------------------------------#	
	wangpo <- simrel[1:ninds, 6]
	wangfs <- simrel[(ninds+1):(2*ninds), 6]
	wanghs <- simrel[((2*ninds) + 1):(3*ninds), 6]
	wangur <- simrel[((3*ninds)+1):(4*ninds), 6]
	lynchlipo <- simrel[1:ninds, 7]
	lynchlifs <- simrel[(ninds+1):(2*ninds), 7]
	lynchlihs <- simrel[((2*ninds) + 1):(3*ninds), 7]
	lynchliur <- simrel[((3*ninds)+1):(4*ninds), 7]
	lynchrdpo <- simrel[1:ninds, 8]	
	lynchrdfs <- simrel[(ninds+1):(2*ninds), 8]
	lynchrdhs <- simrel[((2*ninds) + 1):(3*ninds), 8]
	lynchrdur <- simrel[((3*ninds)+1):(4*ninds), 8]
	quellergtpo <- simrel[1:ninds, 10]
	quellergtfs <- simrel[(ninds+1):(2*ninds), 10]
	quellergths <- simrel[((2*ninds) + 1):(3*ninds), 10]
	quellergtur <- simrel[((3*ninds)+1):(4*ninds), 10]
	
	#-------------------------------------------------------------#
	# Create a list of labels for the different estimators, with
	# each repeated the appropriate number of times
	#-------------------------------------------------------------#
	wang <- rep("W", ninds)
	lynchli <- rep("L & L", ninds)
	lynchrd <- rep("L & R", ninds)
	quellergt <- rep("Q & G", ninds)
	estimator2 <- c(wang, lynchli, lynchrd, quellergt)
	Estimator <- rep(estimator2, 4)
	
	#-------------------------------------------------------------#
	# Create a list of labels for the different relatedness types
	#-------------------------------------------------------------#
	po <- rep("Parent-Offspring", (4*ninds))
	fs <- rep("Full-Sibs", (4*ninds))
	hs <- rep("Half-Sibs", (4*ninds))
	ur <- rep("Unrelated", (4*ninds))
	relationship <- c(po, fs, hs, ur)
	
	#-------------------------------------------------------------#
	# Combine the different values for each estimator based
	# on relatedness type, as lists
	#-------------------------------------------------------------#	
	relatednesspo <- c(wangpo, lynchlipo, lynchrdpo, quellergtpo)
	relatednessfs <- c(wangfs, lynchlifs, lynchrdfs, quellergtfs)
	relatednesshs <- c(wanghs, lynchlihs, lynchrdhs, quellergths)
	relatednessur <- c(wangur, lynchliur, lynchrdur, quellergtur)
	Relatedness_Value <- c(relatednesspo, relatednessfs, relatednesshs, relatednessur)
	
	#-------------------------------------------------------------#
	# Combine the data
	#-------------------------------------------------------------#		
	combineddata <- as.data.frame(cbind(Estimator, relationship, Relatedness_Value))
	combineddata$Relatedness_Value <- as.numeric(as.character(combineddata$Relatedness_Value))
	
	#-------------------------------------------------------------#
	# Plot the data
	#-------------------------------------------------------------#	
	ggplot(combineddata, aes(x = Estimator, y = Relatedness_Value), ylim = c(-0.5, 1.0)) +
	geom_boxplot() +
	facet_wrap(~ relationship)	
}
