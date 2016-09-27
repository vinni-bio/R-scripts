#! /usr/bin/env Rscript

step_bar <- 50


### PACKAGES

list.of.packages <- c("ggplot2", "seqinr", "viridis", "ggthemes")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

suppressMessages(require(seqinr))
suppressMessages(require(ggplot2))
suppressMessages(require(viridis))
suppressMessages(require(ggthemes))

### READING USER INPUT
cat("\n")
cat("Please provide the path to your parsed RELIC-match file: \n")
RELIC.file <- readLines(con=file("stdin"),n=1)
RELIC.file <- gsub("\\s", "", RELIC.file)
cat("\n")
cat("Please provide the path to your protein sequence FASTA file: \n")
PROTEIN.file <- readLines(con=file("stdin"),n=1)
PROTEIN.file <- gsub("\\s", "", PROTEIN.file)
cat("\n")
cat("Please provide the path for the folder to store your output heatmap graph\n(If nothing is provided, the file will be saved into the current folder): \n")
output.path <- readLines(con=file("stdin"),n=1)
output.path <- gsub("\\s", "", output.path)
closeAllConnections()

if(nchar(output.path) == 0){
	output.path <- paste(getwd(), "/",sep="")
} else{
	if(file.info(output.path)$isdir==T){
		if(substr(output.path,nchar(output.path),nchar(output.path)) != "/"){
			output.path <- paste(output.path,"/",sep="")
		}
	}
	else{
		stop("Sorry, but I can't find your output folder. Please try again.")		
	}
} 

cat("\n")
cat("RELIC_heatmap.png will be saved to the following folder:\n")
cat(output.path)
cat("\n")

#RELIC.file <- "~/Documents/2016/EBOLA/NORM/RELIC/RELIC_parsed.txt"
#PROTEIN.file <- "~/Documents/2016/EBOLA/EBOV_GP.faa"

ebolaseq <- read.fasta(PROTEIN.file, seqtype ="AA", as.string = F, set.attributes = FALSE)

dat <- read.table(RELIC.file, header=F, sep="\t", stringsAsFactors = F)
names(dat) <- c("AB", "Start","End")

counts <- data.frame(Pos=1:length(ebolaseq[[1]]),AA=ebolaseq[[1]], counts=rep(0,length(ebolaseq[[1]])))

for(i in 1:nrow(dat)){
	seq_num <- dat$Start[i]: dat$End[i]
	counts[seq_num,3] <- counts[seq_num,3]+1
}

counts$hundreds <- (counts$Pos-1)%/%100*100
counts$numbs <- counts$Pos - counts$hundreds
counts$hundreds <- as.factor(counts$hundreds)
counts$numbs <- as.factor(counts$numbs)

x <- max(counts$counts)

gg <- ggplot(counts, aes(x=hundreds,y= numbs,fill=counts))
gg <- gg + geom_tile(color="white", size=0.01)
gg <- gg+ scale_fill_viridis(name="# of HITS", limits=c(0,x),breaks=seq(0,x, step_bar), guide=guide_colourbar(ticks=F, barheight=unit(0.8,"npc"),label.theme = element_text(colour = "blue", angle = 0, size=12), title.theme=element_text(colour = "blue", angle = 0, size=14, face="bold")))+ labs(x=NULL, y=NULL) + scale_y_discrete(breaks=seq(2,100,2))
gg <- gg + theme_tufte(base_family="Arial") 
gg <- gg + theme(axis.text.x=element_text(size=14), axis.ticks=element_blank(), legend.title=element_text(size=14),plot.title = element_text(size = 14, face = "bold"))
gg <- gg+ggtitle(paste(c("Sites with RELIC matches (n = ",nrow(dat),") on protein sequence"),collapse=""))


png(file=paste(output.path,"RELIC_heatmap.png", sep=""), res=300, width=10,height=10,units="in")
gg
invisible(dev.off())
