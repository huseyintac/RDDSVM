#installs required packages
packages <- c("BSgenome.Hsapiens.UCSC.hg19", "e1071", "optparse")

if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

#calls required packages
suppressMessages(library("BSgenome.Hsapiens.UCSC.hg19"))
suppressMessages(library ( "e1071" ))
suppressMessages(library("optparse"))


#Argument options
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="input file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)



if (is.null(opt$file)){
  message_help(opt_parser)
  stop("An input file name should be provided.", call.=FALSE)
}

#function for retrieving sequences for query

get_query_sequence <- function(x ,length ){
  
  if(query_data[x,3]=="+") {
    
    as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19,
                        paste0("chr",query_data[x,1]), start = query_data[x,2]-length ,
                        #query_data[x,1], start = query_data[x,2]-length ,
                        end = query_data[x,2] + length, strand = "+"))
  } else {
    
    as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19,
                        paste0("chr",query_data[x,1]), start = query_data[x,2]-length ,
                        #query_data[x,1], start = query_data[x,2]-length ,
                        end = query_data[x,2] + length, strand = "-"))
  }
}

#reading query file



 message("Reading input file..")
 message (opt$file)


 query_data <- read.table(opt$file)




 
message(paste("Completed reading",nrow(query_data),"editing points"))
message("Obtaining sequences for input file..")

for (i in 1:nrow(query_data)) {
    query_data[i,4] <- get_query_sequence (i,50)
}


message(paste("Sequence obtained for ",nrow(query_data),"editing points"))

SVM_query <- matrix(ncol = 65, nrow = 0)
SVM_query <- as.data.frame(SVM_query)
SVM_query[1:nrow(query_data),] <- 0


triplet_columns <- c("group","UUU", "UUC", "UUA", "UUG", "UCU", "UCC", "UCA", "UCG","UAU", "UAC", "UAA","UAG","UGU","UGC", "UGA","UGG","CUU","CUC","CUA","CUG","CCU", "CCC", "CCA", "CCG", "CAU", "CAC", "CAA",
                     
                     "CAG", "CGU", "CGC", "CGA", "CGG", "AUU", "AUC", "AUA","AUG","ACU", "ACC","ACA", "ACG", "AAU", "AAC", "AAA", "AAG", "AGU","AGC","AGA","AGG","GUU","GUC","GUA", "GUG","GCU", "GCC",
                     
                     "GCA", "GCG", "GAU", "GAC", "GAA", "GAG","GGU","GGC","GGA","GGG") 

colnames(SVM_query) <- triplet_columns


#Calculation of triplet information

message("Calculating triplet information in sequences")

for (i in 1:nrow(query_data)) {
  
  query_data[i,4] <- gsub("T","U",query_data[i,4])
  
  triplet_arr <- unlist(strsplit(query_data[i,4],""))
  
  for (k in 1: (length(triplet_arr)-2)) {
    
    
    triplet <- paste (triplet_arr[k],triplet_arr[k+1],triplet_arr[k+2],sep = "")
    
    SVM_query[i,triplet] <- SVM_query[i,triplet]+1
    
  }
  
}

# Normalization of triplet information
message("Normalizing triplet information")

for (i in 1:nrow(SVM_query)) {
  rowSum <- rowSums(SVM_query[i,2:65])
  
  for (k in 2:65) {
    SVM_query[i,k] <- round(SVM_query[i,k]/rowSum,2)
    
  }
  
}

# Reading of training data
SVM_query$group <- 1

message("Reading training data")
structure_data <- read.table("SVM_data.txt",header = TRUE,check.names = FALSE, stringsAsFactors = FALSE)
SVM_data <- structure_data
SVM_data <- structure_data[,-(1:3)]
SVM_data <- SVM_data [,-2]
SVM_data <- SVM_data [,-(2:33)]


SVM_data[SVM_data$group=="editing",]$group <- 1  
SVM_data[SVM_data$group=="control",]$group <- -1  
SVM_data$group <- as.factor(SVM_data$group)


train <- sample (10000,5000)

training_data <- SVM_data[train,]

message("Training SVM")
tune.out3 = tune ( svm , group~. , data =training_data, kernel ="polynomial" ,
                   ranges = list ( cost = c (1) ,
                                   degree = c (3) ) )
message("SVM training results:")
summary(tune.out3)

bestmod3 = tune.out3$best.model

ypred3 = predict (bestmod3,SVM_query)
test <- table ( predict = ypred3, truth = SVM_query$group)
message("Query results:")
message(paste("Number of sites predicted as true editing sites ->",sum(ypred3==1)))
message(paste("Number of sites predicted as false editing sites ->",sum(ypred3==-1)))

query_data$V5 <- ypred3

message("Creating output file:")
message(opt$out)
write.table(query_data, opt$out, sep = "\t" , col.names= FALSE ,row.names = FALSE, quote=FALSE)
