# STEP 1: DOWNLOAD
options(stringsAsFactors = F)
project="TCGA-THCA"
library(stringr)
project='TCGA-THCA'
if(!dir.exists("clinical"))dir.create("clinical")
if(!dir.exists("expdata"))dir.create("expdata")
dir()

command1<-"./gdc-client download -m gdc_manifest_clinical.txt -d clinical"
command2<-"./gdc-client download -m gdc_manifest_expdata.txt -d expdata"

system(command = command1) #download clinical data
system(command = command2) #download expression data

length(dir("./clinical/"))
length(dir("./expdata/"))

# STEP 2: CLINICAL DATA
library(XML)
xmls = dir("clinical/",pattern = "*.xml$",recursive = T)
cl = list()
for (i in 1:length(xmls)) {
  result = xmlParse(paste0("clinical/",xmls[[i]]))
  rootnode = xmlRoot(result)
  cl[[i]] = xmlToDataFrame(rootnode[2])
}
clinical = do.call(rbind,cl)
clinical[1:3,1:3]

# STEP 3: EXPRESSION MATRIX
count_files = dir("expdata/",pattern = ".tsv$",recursive = T)

exp = list()
for (i in 1:length(count_files)) {
  exp[[i]] =
    read.table(paste0("expdata/",count_files[[i]]),header = T,sep = "\t")
  exp[[i]] = exp[[i]][-(1:4),] # the first 4 rows contain unneeded information
  exp[[i]] = exp[[i]]$unstranded # the fourth column (unstranded) is the count value
}
exp = as.data.frame(do.call(cbind,exp))

dim(exp)

# TCGA ID and file name match
meta = jsonlite::fromJSON("metadata.cart.json")
ID = sapply(meta$associated_entities,
            function(x){x$entity_submitter_id})
file2id = data.frame(file_name = meta$file_name, 
                     ID = ID)

count_files2 = stringr::str_split(count_files,"/",simplify = T)[,2]
table(count_files2 %in% file2id$file_name)

file2id = file2id[match(count_files2,file2id$file_name),]
identical(file2id$file_name,count_files2)

colnames(exp) = file2id$ID

gene_name = data.table::fread(paste0("expdata/",count_files[1]))$gene_name
gene_name = gene_name[-seq(1,4)] # the first 4 rows contain unneeded information

exp = cbind(gene_name=gene_name,exp)
dim(exp)

exp = exp[!duplicated(exp$gene_name),]
rownames(exp) = exp$gene_name
exp = exp [,-1]
dim(exp)
# [1] 60660   571

# STEP 4: GENE FILTRATION
# keep the genes expressed in more than 50% samples
exp = exp[apply(exp, 1, function(x) sum(x>0) > 0.5*ncol(exp)),]
dim(exp)
# [1] 30937   570

# STEP 5: GROUP INFORMATION
# identify the normal and tumor samples
library(stringr)
table(str_sub(colnames(exp),14,15))
Group = ifelse(as.numeric(str_sub(colnames(exp),14,15)) < 10,'tumor','normal')
Group = factor(Group,levels = c("normal","tumor"))
table(Group)
# normal  tumor 
#  59    511 

# STEP 6: SAVING
if(!dir.exists("data"))dir.create("data")
save(exp,clinical,Group,project,file = paste0("data/",project,"_gdc.Rdata"))
