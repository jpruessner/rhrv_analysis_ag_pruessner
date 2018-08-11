###############################################################
### R script for                                          #####
### 1 - automatic split in subfiles of main HRV datafile  #####
### 2 - automatic analysis of subfiles for HRV            #####
### 3 - automatic saving in datafile for subsequent ANOVA #####
###                                                       #####
### EMPRA SS 2018                                         #####
###############################################################

# libraries required:
# RHRV
library(RHRV)
# parsedate
library(parsedate)

# define function for loading data, splitting into fragments, and save each fragment into a new file 

splitRawData <- function(code,datum) {
  
  eventsFileName <- paste(datum, '_Events', '_', code, '.csv',  sep = '');
  RRFileName <- paste(datum, '_RR', '_', code, '.csv',  sep = '');
  KubiosFileName <- paste(datum, '_RR_KUBIOS', '_', code, '.txt',  sep = '');
  
  events <- read.csv2(eventsFileName, sep = ",", dec = ".", header = T)     
  RR <- read.csv2(RRFileName, sep = ",", dec = ".", header = T)
  kubios <- read.csv2(KubiosFileName, sep = ",", dec = ".", header = F)
  
  colnames(kubios) <- c('kubios'); 
  RRkubios <- cbind(RR, kubios);
  colnames(RRkubios)[1] <- 'dates';
  colnames(events)[1] <- 'dates';
  colnames(RR)[1] <- 'dates';
  colnames(RRkubios)[3] <- 'rrcum';
  
  RRkubios <- within(RRkubios, {
    dates <- parse_date(as.POSIXct(dates)); 
  });
  
  events <- within(events, {
    dates <- parse_date(as.POSIXct(dates)); 
  });
  
  # tos_event = time of start (date and time column 1) of presentation in events file
  # we have identical entry date and time column 1) in RRkubios, so this can be matched to select starting time 
  
  tos_event <- (events[1,1])
  
  # we have 45 pictures for which we know the valence and the order
  # can be defined through three vectors for saving the rr data into respective files indexed by valence and number
  
  pos_picts <- c(1,5,12,18,20,21,23,28,30,33,34,35,38,39,42)
  neut_picts <- c(2,7,9,10,11,14,17,22,27,29,31,36,37,40,44)
  neg_picts <- c(3,4,6,8,13,15,16,19,24,25,26,32,41,43,45)
  
  # now we define a loop for all 45 pictures that selects the respective sections from RRkubios
  # and writes them into separate files in the current directory all indexed by number
  # and valence of the picture, to prepare for subsequent statistical analysis
  
  # event.nr assigns a new variable to RRkubios based on current picture number
  # non-relevant rr sections are assigned 0
  # the rest is assigned 1 to 15 based on start time of presentation (tos_event + 20)
  # and duration of black screen period (30 secs)
  # code could be changed here to adjust analysis for other aspects of dataset, e.g. HRV during pictures themselves
  
  RRkubiosFinal <- within(RRkubios, {
    event.nr <- rep(0,nrow(RRkubios));
})
  
  for (indvdl in 1:45) {
    start_time <- tos_event + 20
    low_tshld <- start_time + ((indvdl-1) * 40)
    upp_tshld <- start_time + ((indvdl * 40)-10)
    RRkubiosFinal$event.nr[RRkubiosFinal$dates > low_tshld & RRkubiosFinal$dates < upp_tshld ] <- indvdl
  }
  
  # now we have assigned picture numbers to each row in the datafile, indicating
  # which rr interval is associated with which picture

  # select each segment and write to file, indicating picture number and valence
  # content of file is RR interval only (all that is needed for HRV analysis)
  
  for(index in 1:45) {
    
    pospic <- index %in% pos_picts
    neutpic <- index %in% neut_picts
    negpic <- index %in% neg_picts
    if (pospic == "TRUE") {val <-"pos"}
    if (neutpic == "TRUE") {val <-"neu"}
    if (negpic == "TRUE") {val <-"neg"}
    
    filename <- paste(code,'_',index,'_',val,'.csv', sep = '')
    write.table(
      RRkubiosFinal$kubios[RRkubiosFinal$event.nr == index], 
      file = filename, 
      row.names = F,
      col.names = F
    )
    
  }
  
  }
  
# define function for hrv analyis of each of the created files 

processHRVData <- function(code) {
  
  
  pos_picts <- c(1,5,12,18,20,21,23,28,30,33,34,35,38,39,42)
  neut_picts <- c(2,7,9,10,11,14,17,22,27,29,31,36,37,40,44)
  neg_picts <- c(3,4,6,8,13,15,16,19,24,25,26,32,41,43,45)
  
  result = list(); 
  
  for(index in 1:45) {
    
      pospic <- index %in% pos_picts
      neutpic <- index %in% neut_picts
      negpic <- index %in% neg_picts
      if (pospic == "TRUE") {val <-"pos"}
      if (neutpic == "TRUE") {val <-"neu"}
      if (negpic == "TRUE") {val <-"neg"}
      filename <- paste(code,'_',index,'_',val,'.csv', sep = '')
      
      hrv.data <- CreateHRVData()
      hrv.data <- LoadBeatRR(hrv.data, RecordName = filename, scale =0.001, verbose = NULL)
      hrv.data <- BuildNIHR(hrv.data)
      hrv.data <- FilterNIHR(hrv.data)
      hrv.data <- InterpolateNIHR(hrv.data, freqhr = 4)
      # show the plot for the HR data to get an idea of the results
      PlotHR(hrv.data, main = c("HR data of event",index))
      
      # size is interval length within epochs for calculation of hrv, should be 30 <= size <= 300
      # we are exceptionally using 10 second intervals to allow calculation of SDs
      hrv.data = CreateTimeAnalysis(hrv.data, size = 10, interval = 7.8125)
      # we are not doing frequency analysis in this analysis
      #hrv.data = CreateFreqAnalysis(hrv.data)
      
      result[[index]] <- hrv.data
      
    }
  
  return(result);     
  
}

# define function for saving the obtained results to file 

writeFile <- function(code, datum) {
  
  
  pos_picts <- c(1,5,12,18,20,21,23,28,30,33,34,35,38,39,42)
  neut_picts <- c(2,7,9,10,11,14,17,22,27,29,31,36,37,40,44)
  neg_picts <- c(3,4,6,8,13,15,16,19,24,25,26,32,41,43,45)
  
  # once data is analysed, write it to file so it can be used subsequently for statistical analysis
  # 
  # remove any leftover files or content in case this is executed repeatedly for several subjects
  
  endResultR= c();
  endResultP= c();
  if (exists("endResult.df")) {rm(endResult.df)}
  
  # load behavioral data file for allowing to include condition in endResult.df
  tempdata.df <- read.csv("~/Dropbox/ExPra_SS2018/ExPra_daten.csv", sep = ",")
  tempdata.df <- tempdata.df[-1,]
  
  for(index in 1:length(processedData) ) {
    endResultR <- c(endResultR, processedData[[index]]$TimeAnalysis[[1]]$rMSSD) 
  }
  
  for(index in 1:length(processedData) ) {
    endResultP <- c(endResultP, processedData[[index]]$TimeAnalysis[[1]]$pNN50) 
  }
  
  endResult.df <- data.frame(cbind(endResultP,endResultR))
  endResult.df$pict_nr <- c(1:length(endResultP))
  endResult.df$code[1:length(endResultP)] <- code
  endResult.df$beding <- NA
  #endResult.df$beding <- as.factor(endResult.df$beding)
  # retrieve the group from the behavioral data file and add it to the current file as variable for use in subsequent ANOVA
  endResult.df$beding <-  tempdata.df$group[which(tempdata.df$Code == endResult.df$code[1])]
  
  for(index in 1:45) {
    
    pospic <- index %in% pos_picts
    neutpic <- index %in% neut_picts
    negpic <- index %in% neg_picts
    if (pospic == "TRUE") {val <-"pos"}
    if (neutpic == "TRUE") {val <-"neu"}
    if (negpic == "TRUE") {val <-"neg"}
    
    endResult.df$valence[index] <- val
  }
  
  # sort the data by valence to create time series index for ANOVA
  endResult.df <- endResult.df[with(endResult.df, order(endResult.df$valence, endResult.df$pict_nr)), ]
  
  # create a new variable for repeat picture appearance
  endResult.df$time <- (rep(1:15, 3))
  
  # give proper names to your variables
  colnames(endResult.df) <- c("pnn50","rmssd","picture", "code", "beding", "valence", "time")
  
  # get rid of picture variable since it has served its purpose
  endResult.df$picture <- NULL
  
  # prepare file name
  resultsFileName <- paste(datum, '_results', '_', code, '.csv',  sep = '');
  
  # and write to file
  write.csv(endResult.df, file = resultsFileName)
  
}


#############################################################
#########  define all necessary variables for execution #####
#############################################################

datum <- '2018-6-11';
code <- 'GB26';

# set working directory to where your data is
# sewd("/path/to/where/your/data/is")

splitRawData(code,datum); 
processedData <- processHRVData(code); 
writeFile(code,datum)

####################################################
####  end of RHRV analysis                     #####
####################################################
