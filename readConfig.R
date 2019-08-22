#### title:  readConfig.R
####
#### purpose:  code to parse and eval config file
####

#### read config file
daConfig <- read.table(configFile,as.is=TRUE,header=TRUE)

#### generate two sets of R-commands with different formatting
#### why do we need two?  i don't know, but some rows will not eval correctly without them
commands1 <- paste(daConfig$Variable," <- ",daConfig$Value,sep="")
commands2 <- paste(daConfig$Variable," <- ","\"",daConfig$Value,"\"",sep="")

#### execute commands
for (ii in 1:length(commands1)) {
  if (inherits(try(eval(parse(text=commands1[ii])),silent=TRUE),"try-error")) {
    eval(parse(text=commands2[ii]))
  }
}

#### load libraries
if (exists("libRequests")) {
  libraries <- unlist(strsplit(libRequests,","))
  if (length(libraries)>0) {
    for (ii in 1:length(libraries)) {
      eval(parse(text=paste("require(",libraries[ii],")",sep="")))
    }
  }
}

#### source functions

if (exists("dirNameFuncs")) {
  if (length(unlist(strsplit(fileNameFuncs,",")))>0) {
    fileFuncs <- paste(dirNameFuncs,unlist(strsplit(fileNameFuncs,",")),".R",sep="")
    for (ii in 1:length(fileFuncs)) {
      source(paste(fileFuncs[ii]))
    }
  }
}
