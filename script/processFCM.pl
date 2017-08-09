#!/usr/bin/perl

use strict;
use warnings;
use feature 'say';
use String::Util qw(trim);

=pod
 
=head1 INTRODUCTION

This is the pipeline for the FlowCytometry data processing scripts. 

The code generates five R scripts and one SLURM batch script. 

Depending upon the command line options, the code will submit the job on compute nodes or will
run the R scripts on the command line of the login node. 

Running on compute nodes is recommended.

To run the script, simply invoke it like so:

./ProcessFCM.pl [cluster|clean]

Note that the cluster parameter is optional and if not given, the generated R
scripts will run on login node. 

If the clean parameter is used the generated files will get cleaned.

=head1 CONFIGURATION

Configuration parameters for the runs are specified in a file called conf.conf.

The file simply has the list of names and their values separated by '='. 

There should be no space in the either side of '='.

Comments may be added to the configuration file by putting '#' in the beginning of the line. 

Partial line comments are not supported at this time.

=head1 COMPUTATION

The current script accomplishes the following tasks:

1. Provides a built-in documentation that you are reading now.

2. Parses the configuration file and plugs the values from configuration file into generated R scripts.

3. Generates the R and SLURM scripts.

4. Runs the R scripts on login node or submits the SLURM scripts to compute nodes depending on the command line option provided.

So, the actual computation is done by running the R scripts.

=head1 DATA

Make sure that the input data is present in the path that is specified in the configuration file as mentioned above. 

Just the path of directory is required for the FCS file.

To get this path cd into the directory containing the FCS file and run the command pwd.

Copy the output of this command and paste it into the configuration file at right place.

For other files, such as the gates file and compensation matrix file, full path of the file is needed.

To get the full path append the filename to the path of the directory in which the file is located.

=head1 TECHNICAL DETAILS

This one and the next section are mostly for the future maintainers of this code, including myself.

The details are for code maintenance. The script is written in Perl. 

It use the heredoc feature to embed R and SLURM scripts in the single file.

There are several parameters that are in the scripts that could be modified.

If you identify such parameters, it is straightforward to pick them up from the R script within the code added to the conf.conf file and used just as other such parameters are being used.

=head1 TROUBLESHOOT

The code may fail at a few places depending upon configuration of the system, data and R packages.

Often times one or more individual R scripts may fail. In such cases, it is wise to look into the generated failed R script outside the current script and remedy it.

Once the failed R script is fixed, the fixes may be incorporated back into the current script.

Some of the items to verify are shown in the following checklist:

1. Make sure the required R packages are installed and available. They are FlowSOM, FlowCore, etc.

2. Make sure the data is present in the place that you are telling the script to find. In other words, the values in conf.conf must be in agreement to those that form the parameters for processing.

3. Make sure you have access and allocation to CRC systems.

=cut

# If the command line parameter is clean--deleted generated files and quit
if (@ARGV and $ARGV[0] eq "clean"){
    my @filestoclean = ('plotsettings.R', 'preprocessing.R', 
		        'Rtsne.R',        'Flowsom.R',      
		        'Spade.R',        'Spade_10000.R', 
		        'jobdef.sbatch');
    
    my $status = system('rm', @filestoclean);
    exit 0;
    say "This should never be said";
}

# Parse the conf.conf file to extract the parameters
open FILE, "conf.conf" or die "Can't open file";
my %hash;
my @val;
while (<FILE>){
    next if (/^\s+$/); 
    next if (/^#/);
    @val=split(/=/,$_);
    $hash{$val[0]}=trim($val[1]);
}
close FILE;        


# Fetch parameters into local variables
# The parameters will be plugged into the R script that are generated
my $fcspath          = $hash{"fcspath"};
my $gatespath        = $hash{"gatespath"};
my $compmatpath      = $hash{"compmatpath"};
my $resultspath      = $hash{"resultspath"};
my $plotsettingspath = $hash{"plotsettingspath"};
my $fcsrdatapath     = $hash{"fcsrdatapath"};

#Start building the plotsettings.R
my $plotsettings = <<"END_PLOT_SETTINGS";

cellTypeColors <- c("total MBC" = "#e74c3c",
"activated MBC subset" = "#ff3399",
"classical MBC subset" = "#3498db",
"total naive B" = "#a65628", 
"atypical naive B subset" = "#ff7f00", 
"classical naive B subset" = "#fecc08",
"plasmablast" = "#9b59b6")

markerNames <- c("Comp-FITC-A" = "CD45RB",
"Comp-BV421-A" = "CD27",
"Comp-PE-Cy5-A" = "IgM",
"Comp-PerCP-Cy5-5-A" = "CD21",
"Comp-PE-Cy7-A" = "IgD",
"Comp-BV510-A" = "CD69",
"Comp-APC-A"="CD11a",
"Comp-PE-A" = "CD35",
"Comp-BUV 395-A" = "CD38",
"Comp-APC-Cy7-A" = "CD19")

markerColors <- c("Comp-FITC-A" = "#fdae6b",
"Comp-BV421-A" = "#e41a1c",
"Comp-PE-Cy5-A" = "#ff7f00",
"Comp-PerCP-Cy5-5-A" = "#fecc08",
"Comp-PE-Cy7-A" = "#ae017e",
"Comp-BV510-A" = "#9b59b6",
"Comp-APC-A" = "#3498db",
"Comp-PE-A" = "#1f78b4",
"Comp-BUV 395-A" = "#a4cc2e",
"Comp-APC-Cy7-A" = "#fb9a99")

circular_markerOrder = c(7,8,9,10,11,12,13,14,15,16)
grid_markerOrder = c(7,8,9,10,11,12,13,14,15,16)

save(cellTypeColors,markerColors,markerNames,circular_markerOrder,grid_markerOrder, file="$plotsettingspath")
END_PLOT_SETTINGS

open my $plotsettings_fh, '>', 'plotsettings.R';
say $plotsettings_fh $plotsettings;
close $plotsettings_fh;

#========================================================

#Start building the preprocessing.R
my $preprocessing = <<"END_PREPROCESSING";
library(flowCore)
library(methods)

#dataID <- "../data"
dataID <- "$fcspath"

######################
# Parameter settings #
######################

# Files to use
files = list.files(dataID, pattern=".fcs", full.names = TRUE)

# Compensation matrix
compensationFile = file.path("$compmatpath")
colsToCompensate = c(7:17)

# GatingML file
gatingFile <- file.path("$gatespath")

gate_ids <- c("Live single CD19+ cells" = 4,
"activated MBC subset" = 7,
"classical MBC subset" = 8,
"atypical naive B subset" = 10,
"classical naive B subset" = 11,
"plasmablast" = 5)

cellTypes <- c("activated MBC subset",
               "classical MBC subset",
               "atypical naive B subset",
               "classical naive B subset",
               "plasmablast"
	       )

# Columns to use for analysis
colsToCluster = c(7:16)

###################
# Helper function #
###################

ProcessGatingML <- function(flowFrame,gatingFile,gateIDs, cellTypes,silent=FALSE){
    gating_xml <- XML::xmlToList(XML::xmlParse(gatingFile))
    flowEnv <- new.env()
    flowUtils::read.gatingML(gatingFile, flowEnv)
    #  A. Read the gates from xml
    filterList <- list()
    for(cellType in names(gateIDs)){
        filterList[[cellType]] <-  flowEnv[[ as.character(gating_xml[[gateIDs[cellType]]]\$.attrs["id"]) ]]
    }
    #  B. Process the fcs file for all specified gates
    results <- matrix(NA,nrow=nrow(flowFrame),ncol=length(gateIDs), dimnames = list(NULL,names(gateIDs)))
	    
    for(cellType in names(gateIDs)){
        if(!silent){message(paste0("Processing celltype ",cellType))}
        results[,cellType] <- flowCore::filter(flowFrame, filterList[[cellType]])\@subSet
    }
    #  C. Assign one celltype to each cell
    manual <- rep("Unknown",nrow(flowFrame))
    for(celltype in cellTypes){
        manual[results[,celltype]] <- celltype
    }
    manual <- factor(manual,levels = c("Unknown",cellTypes))
    
    list("matrix"=results,"manual"=manual)
}

########################
# Preprocess the data  #
########################

# Read the compensation matrix
comp <- read.csv(compensationFile, row.names=1, check.names = FALSE)
colnames(comp) <- rownames(comp) <- gsub(" ::.*","",colnames(comp))

for(file in files){
    message(paste0("Processing file ",file))
    # Load the raw data
    ff <- read.FCS(file)
    # and compensate
    ff <- compensate(ff,comp)
    ff\@description\$SPILL <- comp
    colnames(ff)[colsToCompensate] <- paste("Comp-", colnames(ff)[colsToCompensate],sep="")
    
    # Extract manual gating
    gatingRes <- ProcessGatingML(ff, gatingFile, gate_ids, cellTypes)
    gatingMatrix <- gatingRes\$matrix
    manual <- gatingRes\$manual
    
    # Cells to use as input for the algorithms
    selected <- gatingMatrix[,"Live single CD19+ cells"]
    # Save selected cells to fcs file for algorithms which can only take
    # a file as input. File is already compensated, so identity matrix as
    # compensation matrix
    new_comp <- diag(length(colnames(ff)[7:16]))
    colnames(new_comp) <- colnames(ff)[7:16]
    ff\@description\$SPILL <- new_comp
    write.FCS(ff[selected,], file=gsub(".fcs","_selected.fcs",file))
    # Save also a subset of only the first 10.000 cells for faster processing
    write.FCS(ff[selected,][1:10000,], file=gsub(".fcs","_selected_10000.fcs",file))
    
    # Transform the data
    ff_t <- transform(ff,transformList(colnames(ff)[7:16],logicleTransform()))
    
    # Save results so this step can be skipped next time
    save(ff,ff_t,selected,manual,gatingMatrix,colsToCluster,
         file=gsub(".fcs",".Rdata",file))
}
END_PREPROCESSING

open my $preprocessing_fh, '>', 'preprocessing.R';
say $preprocessing_fh $preprocessing;
close $preprocessing_fh;

#========================================================

#Start building the Rtsne.R
my $Rtsne = <<"END_RTSNE";
# Load the preprocessed data:
# See script_preprocessing.R
#    ff:            Compensated flowFrame
#    ff_t:          Compensated and logicle transformed flowFrame
#    manual:        Array with label for each cell
#    selected:      Array with TRUE/FALSE whether cell falls in single live 
#                   cells
#    gatingMatrix:  Matrix with rows corresponding to cells and a column for 
#                   each manual gate. Each column contains TRUE/FALSE values
#                   indicating whether the cells fall in the specific gate
#    colsToCluster: Columns to use for clustering
load("$fcsrdatapath") #../data/D197_B_ST1_Tube_001_015.Rdata

# Load the plot settings
# See script_plotSettings.R
#    cellTypeColors,markerColors,markerNames,
#    circular_markerOrder,grid_markerOrder
load("plotSettings.Rdata")

# Load the FlowSOM library
library(Rtsne)
library(flowCore)

# Set some parameters
tSNE_subsample <- 10000

# Set seed for reproducable results
set.seed(42)
# Record start time
start <- Sys.time()
    # Run the rtsne algorithm
    rtsne_res <- Rtsne(exprs(ff_t[selected,colsToCluster])[seq_len(tSNE_subsample),])
    # Plot the results, marker plots in separate png because pdf is too big
    
    # Plot manual
    pdf("Rtsne_manual.pdf",useDingbats = FALSE)
    plot(rtsne_res\$Y,col=c("#888888",cellTypeColors)[manual[selected][seq_len(tSNE_subsample)]],pch=19,
         bty="n",axes=F,xlab="",ylab="")
    dev.off()
    
    # Plot individual markers
    png("Rtsne_markers.png",width = 1200,height=1200)
    par(mfrow=c(4,4))
    for(m in grid_markerOrder){
        channel <- colnames(ff)[m]
        plot(rtsne_res\$Y,
             col=colorRampPalette(c("#dddddd",markerColors[channel]))(100)[
                 as.numeric(cut(exprs(ff_t[selected,channel])[seq_len(tSNE_subsample),],
                                breaks = 100))],
             main=markerNames[channel],bty="n",axes=F,xlab="",ylab="",pch=19)
    }
    par(mfrow=c(1,1))
    dev.off()
    
# Record end time
t_Rtsne_10000<- Sys.time() - start
# Save results
save(t_Rtsne_10000, rtsne_res, file="rtsne_10000.Rdata")
END_RTSNE

open my $Rtsne_fh, '>', 'Rtsne.R';
say $Rtsne_fh $Rtsne;
close $Rtsne_fh;

#===========================================================

# Start building the Flowsom.R
my $flowsom = << "END_FLOWSOM";
# Load the preprocessed data:
# See script_preprocessing.R
#    ff:            Compensated flowFrame
#    ff_t:          Compensated and logicle transformed flowFrame
#    manual:        Array with label for each cell
#    selected:      Array with TRUE/FALSE whether cell falls in single live 
#                   cells
#    gatingMatrix:  Matrix with rows corresponding to cells and a column for 
#                   each manual gate. Each column contains TRUE/FALSE values
#                   indicating whether the cells fall in the specific gate
#    colsToCluster: Columns to use for clustering
#load("FR-FCM-ZZQY/21-10-15_Tube_028.Rdata")
load("$fcsrdatapath") #../data/D197_B_ST1_Tube_001_015.Rdata

# Load the plot settings
# See script_plotSettings.R
#    cellTypeColors,markerColors,markerNames,
#    circular_markerOrder,grid_markerOrder
load("$plotsettingspath")

# Load the FlowSOM library
library(FlowSOM)
library(flowCore)
library(methods)

# Set some parameters
flowSOM_metaClusters = 10
flowSOM_xdim = 7
flowSOM_ydim = 7

# Set seed for reproducable results
set.seed(42)
# Record start time
start <- Sys.time()
# Run the flowSOM algorithm
fsom <- FlowSOM(ff_t[selected,],
        compensate = FALSE, transform = FALSE, scale = TRUE,
        colsToUse = colsToCluster, xdim=flowSOM_xdim, ydim=flowSOM_xdim,
        nClus = flowSOM_metaClusters)
# Plot the results
pdf("FlowSOM.pdf", useDingbats = FALSE)
PlotPies(UpdateNodeSize(fsom[[1]],reset=T), manual[selected], colorPalette = colorRampPalette(c('#FFFFFF',cellTypeColors)))
PlotStars(UpdateNodeSize(fsom[[1]],reset=T),
#           view = "MST",
#           backgroundValues = as.factor(fsom[[2]]),
           markers = circular_markerOrder,
           colorPalette = colorRampPalette(markerColors[
           colnames(ff)[circular_markerOrder]]),
#             backgroundColor = c("#ff7f0055","#FECC0844",
#                                 "#bdc9e133","#9b59b644",
#                                 '#00000022',"#de2d2655",
#                                 "#e74c3c22","#3498db44",
#                                 "#addd8e44","#FF000066") 
                  # Colors adapted to manual annotation
        )
PlotStars(fsom[[1]],
#          view = "MST",
#          backgroundValues = as.factor(fsom[[2]]),
          markers = circular_markerOrder,
          colorPalette = colorRampPalette(markerColors[
          colnames(ff)[circular_markerOrder]]),
#              backgroundColor = c("#ff7f0055","#FECC0844",
#                                  "#bdc9e133","#9b59b644",
#                                  '#00000022',"#de2d2655",
#                                  "#e74c3c22","#3498db44",
#                                  "#addd8e44","#FF000066") 
	)
dev.off()
	
# Record end time
t_flowSOM<- Sys.time() - start
# Save results
save(t_flowSOM, fsom, file="flowSOM.Rdata")
END_FLOWSOM


open my $flowsom_fh, '>', 'Flowsom.R';
say $flowsom_fh $flowsom;
close $flowsom_fh;

#============================================================

# Start building the Spade.R
my $spade = << "END_SPADE";
# Load the preprocessed data:
# See script_preprocessing.R
#    ff:            Compensated flowFrame
#    ff_t:          Compensated and logicle transformed flowFrame
#    manual:        Array with label for each cell
#    selected:      Array with TRUE/FALSE whether cell falls in single live 
#                   cells
#    gatingMatrix:  Matrix with rows corresponding to cells and a column for 
#                   each manual gate. Each column contains TRUE/FALSE values
#                   indicating whether the cells fall in the specific gate
#    colsToCluster: Columns to use for clustering
load("$fcsrdatapath") #../data/D197_B_ST1_Tube_001_015.Rdata

# Load the plot settings
# See script_plotSettings.R
#    cellTypeColors,markerColors,markerNames,
#    circular_markerOrder,grid_markerOrder
load("$plotsettingspath")

# Load the SPADE library
library(Biobase)
library(spade)
library(flowCore)

# Set some parameters
SPADE_nClusters = 50

# Set seed for reproducable results
set.seed(42)
# Record start time
start <- Sys.time()
# Run the SPADE algorithm
out_dir <- "SPADE"
SPADE.driver("../data/D197_B_ST1_Tube_001_015_selected.fcs", out_dir=out_dir, cluster_cols=colnames(ff)[colsToCluster], comp=FALSE, transforms=logicleTransform(), k=SPADE_nClusters)

# Plot the results
mst_graph <- igraph:::read.graph(paste(out_dir,"mst.gml", sep=.Platform\$file.sep), format="gml")
layout <- read.table(paste(out_dir,"layout.table",sep=.Platform\$file.sep))
SPADE.plot.trees(mst_graph,out_dir,layout=as.matrix(layout),
                 out_dir=paste(out_dir,"pdf",sep=.Platform\$file.sep))
# Record end time
t_SPADE<- Sys.time() - start


# Fit the SPADE results in a FlowSOM object for the visualization 
# used in the paper.
library(FlowSOM)
spade_fcs <- read.FCS(paste0(out_dir,"/D197_B_ST1_Tube_001_015_selected.fcs.density.fcs.cluster.fcs"))
spade_res <- list("MST"=list(),"map"=list(),"data"=exprs(ff_t[selected,]))
spade_res\$MST\$graph <- mst_graph
spade_res\$MST\$l <- as.matrix(layout)
spade_res\$map\$mapping <- exprs(spade_fcs[,"cluster"])
vsize <- table(spade_res\$map\$mapping)/nrow(spade_res\$map\$mapping)
spade_res\$MST\$size <- (vsize/(max(vsize, na.rm = TRUE)) * 3 + 2) * 4 # same resizing as in SPADE library
#spade_res\$map\$codes <- spade_res\$map\$medianValues <- t(sapply(seq_len(SPADE_nClusters), function(i) {
spade_res\$map\$codes <- spade_res\$map\$meanValues <- t(sapply(seq_len(SPADE_nClusters), function(i) {
    apply(subset(spade_res\$data, spade_res\$map\$mapping[, 1] == i), 
          2, median)
}))
spade_res\$map\$medianValues[is.nan(spade_res\$map\$medianValues)] <- 0
spade_res\$map\$grid <- matrix(1:SPADE_nClusters,nrow=SPADE_nClusters,ncol=1)

# Save results
save(t_SPADE, spade_res, file="SPADE.Rdata")

# Plot results
pdf("SPADE.pdf",useDingbats = FALSE)
PlotPies(UpdateNodeSize(spade_res,reset=T),manual[selected], colorPalette = colorRampPalette(c("#FFFFFF",cellTypeColors)))

# Plot individual markers
par(mfrow=c(3,3))
for(m in grid_markerOrder){
    channel <- colnames(ff)[m]
    PlotMarker(spade_res, channel, colorPalette = colorRampPalette(c("#dddddd", markerColors[channel])), main=markerNames[channel])
}
par(mfrow=c(1,1))
dev.off()
END_SPADE

open my $spade_fh, '>', 'Spade.R';
say $spade_fh $spade;
close $spade_fh;

#========================================================

# Start building the Spade_10000.R
my $spade_10000 = << "END_SPADE_10000";
# Load the preprocessed data:
# See script_preprocessing.R
#    ff:            Compensated flowFrame
#    ff_t:          Compensated and logicle transformed flowFrame
#    manual:        Array with label for each cell
#    selected:      Array with TRUE/FALSE whether cell falls in single live 
#                   cells
#    gatingMatrix:  Matrix with rows corresponding to cells and a column for 
#                   each manual gate. Each column contains TRUE/FALSE values
#                   indicating whether the cells fall in the specific gate
#    colsToCluster: Columns to use for clustering
#load("FR-FCM-ZZQY/21-10-15_Tube_028.Rdata")
load("$fcsrdatapath") #../data/D197_B_ST1_Tube_001_015.Rdata

# Load the plot settings
# See script_plotSettings.R
#    cellTypeColors,markerColors,markerNames,
#    circular_markerOrder,grid_markerOrder
load("$plotsettingspath")

library(Biobase)
library(flowCore)
# Load the SPADE library
library(spade)

# Set some parameters
SPADE_nClusters = 50

# Set seed for reproducable results
set.seed(42)
# Record start time
start <- Sys.time()
    # Run the SPADE algorithm
    out_dir <- "SPADE_10000"
    SPADE.driver("../data/D197_B_ST1_Tube_001_015_selected_10000.fcs",
                 out_dir=out_dir,
                 cluster_cols=colnames(ff)[colsToCluster], 
                 comp=FALSE, transforms=logicleTransform(),
                 k=SPADE_nClusters)
    
    # Plot the results
    mst_graph <- igraph:::read.graph(paste(out_dir,"mst.gml",
                                           sep=.Platform\$file.sep),
                                           format="gml")
    layout <- read.table(paste(out_dir,"layout.table",sep=.Platform\$file.sep))
    SPADE.plot.trees(mst_graph,out_dir,layout=as.matrix(layout),
                     out_dir=paste(out_dir,"pdf",sep=.Platform\$file.sep))
# Record end time
t_SPADE_10000<- Sys.time() - start

# Fit the SPADE results in a FlowSOM object for the visualization 
# used in the paper.
library(FlowSOM)
spade_fcs <- read.FCS(paste0(out_dir,"/D197_B_ST1_Tube_001_015_selected_10000.fcs.density.fcs.cluster.fcs"))
spade_res <- list("MST"=list(),"map"=list(),"data"=exprs(ff_t[selected,][1:10000,]))
spade_res\$MST\$graph <- mst_graph
spade_res\$MST\$l <- as.matrix(layout)
spade_res\$map\$mapping <- exprs(spade_fcs[,"cluster"])
vsize <- table(spade_res\$map\$mapping)/nrow(spade_res\$map\$mapping)
spade_res\$MST\$size <- (vsize/(max(vsize, na.rm = TRUE)) * 3 + 2) * 4 # same resizing as in SPADE library
#spade_res\$map\$codes <- spade_res\$map\$medianValues <- t(sapply(seq_len(SPADE_nClusters), function(i) {
spade_res\$map\$codes <- spade_res\$map\$meanValues <- t(sapply(seq_len(SPADE_nClusters), function(i) {
    apply(subset(spade_res\$data, spade_res\$map\$mapping[, 1] == i), 
          2, median)
}))
spade_res\$map\$medianValues[is.nan(spade_res\$map\$medianValues)] <- 0
spade_res\$map\$grid <- matrix(1:SPADE_nClusters,nrow=SPADE_nClusters,ncol=1)


# Save results
save(t_SPADE_10000, spade_res, file=file.path(out_dir, "SPADE_10000.Rdata"))

# Plot results
pdf(file.path(out_dir, "SPADE_10000.pdf"),useDingbats = FALSE)
    PlotPies(UpdateNodeSize(spade_res,reset=T),manual[selected][0:10000],
         colorPalette = colorRampPalette(c("#FFFFFF",cellTypeColors)))
    
    # Plot individual markers
    par(mfrow=c(3,3))
    for(m in grid_markerOrder){
        channel <- colnames(ff)[m]
        PlotMarker(spade_res,channel,
                   colorPalette = colorRampPalette(c("#dddddd",markerColors[channel])),
                   main=markerNames[channel])
    }
    par(mfrow=c(1,1))
dev.off()
END_SPADE_10000

open my $spade_10000_fh, '>', 'Spade_10000.R';
say $spade_10000_fh $spade_10000;
close $spade_10000_fh;

#========================================================

# Start building the job script.
my $jobscript = <<"END_JOBSCRIPT";
\#!/bin/bash
\#SBATCH -N 1 
\#SBATCH --job-name=set1
\#SBATCH --output=../results/set1/set1.out 
\#SBATCH --error=../results/set1/set1.err 
\#SBATCH --time=00:50:00 
\#SBATCH --cpus-per-task=8 
\#SBATCH --mem=32g 

module purge
module load R
module load gsl 

Rscript plotsettings.R
echo "Finished plotsettings.R"

Rscript preprocessing.R
echo "Finished preprocessing.R"
echo "Set 1 finished"

Rscript Rtsne.R
echo "Finished Rtsne.R"

Rscript Flowsom.R
echo "Finished Flowsom.R"

Rscript Spade.R
echo "Finished Spade.R"

Rscript Spade_10000.R
echo "Finished Spade_10000.R"
echo "Set 2 finished"
END_JOBSCRIPT

open my $jobscript_fh, '>', 'jobdef.sbatch';
say $jobscript_fh $jobscript;
close $jobscript_fh;

if (@ARGV and $ARGV[0] eq "cluster"){
	say "Run over COMPUTE node.";
	my $status=system('sbatch', 'jobdef.sbatch');
}else{
	say "Run over LOGIN node.";
	my $status=system('Rscript', 'plotsettings.R');
           $status=system('Rscript', 'preprocessing.R');
           $status=system('Rscript', 'Rtsne.R');
           $status=system('Rscript', 'Flowsom.R');
           $status=system('Rscript', 'Spade.R');
           $status=system('Rscript', 'Spade_10000.R');
	   say $status;
}

# XXXXXXXXXXXXXXXXXXXXXX THE END XXXXXXXXXXXXXXXXXXXXXXXXXXXXX

