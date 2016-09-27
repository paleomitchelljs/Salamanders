## Run BioGeoBEARS on salamander data
# Load the package (after installation, see above).
library(optimx)         # You need to have some version of optimx available
                        # as it is a BioGeoBEARS dependency; however, if you
                        # don't want to use optimx, and use optim() (from R core) 
                        # you can set:
                        # BioGeoBEARS_run_object$use_optimx = FALSE
                        # ...everything should work either way -- NJM 2014-01-08
library(FD)       # for FD::maxent() (make sure this is up-to-date)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)
library(BioGeoBEARS)
library(phytools)


##########################################################################################
#### Source all of the annoying BGB stuff
source("http://phylo.wdfiles.com/local--files/biogeobears/cladoRcpp.R") # (needed now that traits model added; source FIRST!)
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_add_fossils_randomly_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_basics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_calc_transition_matrices_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_classes_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_detection_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_DNA_cladogenesis_sim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_extract_Qmat_COOmat_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_generics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_models_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_on_multiple_trees_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_plots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_readwrite_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_simulate_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_makePlots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stochastic_mapping_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stratified_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_univ_model_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_uppass_probs_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_loglike_sp_v01.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/get_stratified_subbranch_top_downpass_likelihoods_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/runBSM_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/stochastic_map_given_inputs.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/summarize_BSM_tables_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_traits_v1.R") # added traits model
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)
    # slight speedup hopefully
# However, you can find the extdata directory like this:
extdata_dir <- np(system.file("extdata", package="BioGeoBEARS"))
scriptdir <- np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
list.files(extdata_dir)
##########################################################################################
Path_base <- "~/Documents/Salamanders/"
setwd(paste(Path_base, "Rscripts/", sep=""))
source("plot_BGB.R")

setwd(paste(Path_base, "datafiles", sep=""))

tree <- read.tree("fossilTree.tre")
inputGeog <- read.table("sal_geog.data", stringsAsFactors=F, header=T, row.names=1)
inputGeog <- cbind(tree$tip.label, inputGeog[tree$tip.label,])
Name <- "softFos"

# Toggle this for hardFos (on) or softFos (off)
inputGeog[,2] <- gsub("\\?", "0", inputGeog[,2])
Name <- "hardFos"

inputGeog <- cbind(inputGeog, rep("", nrow(inputGeog)))
colnames(inputGeog) <- c(Ntip(tree), 6, "(A W S Er C E)")

write.table(inputGeog, file="input_geog.data", quote=F, row.names=F, col.names=T, sep="\t")


# d = dispersal, e= extinction (these two are LaGrange)
z_dec <- define_BioGeoBEARS_run()

z_dec$geogfn <- "input_geog.data"
z_dec$trfn <- "fossilTree.tre"
z_dec$include_null_range <- FALSE	# This makes it DEC*xx, which seems to perform better based on doi: http://dx.doi.org/10.1101/026914
z_dec$use_optimx <- T
z_dec$num_cores_to_use <- 1
z_dec$max_range_size <- 6
z_dec$force_sparse <- FALSE
z_dec$useAmbiguities <- TRUE
z_dec <- readfiles_BioGeoBEARS_run(z_dec)
z_dec$return_condlikes_table <- TRUE
z_dec$calc_TTL_loglike_from_condlikes_table <- TRUE
z_dec$calc_ancprobs <- TRUE
tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn=z_dec$geogfn)
check_BioGeoBEARS_run(z_dec)

y_dec <- bears_optim_run(z_dec)

# add in j = founder-event speciation (bears_3param), and v = vicariance proportion. Possible to manipulate max_ent_constraint_05 which impacts the probability of different range sizes (low=small sizes, hi=big sizes)
z_decjv <- z_dec
z_decjv$BioGeoBEARS_model_object@params_table["d", "est"] <- y_dec$outputs@params_table["d", "est"]
z_decjv$BioGeoBEARS_model_object@params_table["e", "est"] <- y_dec$outputs@params_table["e", "est"]
z_decjv$BioGeoBEARS_model_object@params_table["j", "type"] <- "free"
z_decjv$BioGeoBEARS_model_object@params_table["j", "init"] <- 0.001
z_decjv$BioGeoBEARS_model_object@params_table["j", "est"] <- 0.001
z_decjv$BioGeoBEARS_model_object@params_table["v", "type"] <- "free"
z_decjv$BioGeoBEARS_model_object@params_table["v", "init"] <- 0.001
z_decjv$BioGeoBEARS_model_object@params_table["v", "est"] <- 0.001
#z_decjv <- readfiles_BioGeoBEARS_run(z_decjv)
check_BioGeoBEARS_run(z_decjv)

y_decjv <- bears_optim_run(z_decjv)

# Save output to load in analyzeBGB.R
save(y_dec, file=paste(Path_base, "output/sal_dec*_", Name, ".Rdata", sep=""))
save(y_decjv, file=paste(Path_base, "output/sal_dec*jv_", Name, ".Rdata", sep=""))