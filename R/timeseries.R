#' @title Produce a time series that counts up or otherwise summarises the intensity of observed archaeological/historical sites or site phases per year  
#'
#' @description Given a set of archaeological/historical sites or site phases with absolute start and end dates, count up how many sites are definitely or possible present for each year in a specified time range. 
#' 
#' @param startdates A vector of numeric values specifying the approximate start dates of each site or site phase.
#' @param enddates A vector of numeric values specifying the approximate end dates of each site or site phase.
#' @param siteprobs A vector of numeric values between 0 and 1 specifying the probability that the site in question was in use between the specified start and end date (often a quantification of how confidently an archaeological survey or excavation project were able to assign a given site to the specified chronological period based on recovered finds or other evidence).
#' @param siteareas A vector of numeric values >=0 that provide a size estimate of the sites or site phases (e.g. in hectares).
#' @param consensus A two column matrix or data.frame with a first column of years in the specified calendar, and a second column with an alternative set of settlement or population intensities (e.g. a summed radiocarbon distribution) for which a consensus model might be specified. This dataset must cover the full range of start and end dates, as well as of timeRange if this is given.
#' @param conprop Controls the importance of the consensus datset relative to the main sitecount dataset (default is 0.5 which imposes an equal weighting, but higher values would prioritise the consensus dataset more). Either a single number between 0 and 1 or a numeric vector of the same numebr of rows as consensus (each value between 0 and 1)
#' @param timeRange Numeric vector of length 2, specifying the earliest and latest date to include in the results (in the year units of the calendar argument). The default option will calculate this automatically from the input data.
#' @param calendar Specifies whether the input calendar is in years BC (-ve)/AD (+ve) or BP (default is "BCAD"). Note there should be no year 0 in "BCAD".
#' @param aoristic Logical value specifying whether an \"aoristic\" sum should be calculated.
#' @param verbose A logical variable indicating whether extra information on progress should be reported. Default is TRUE.
#' 
#' @return A data.frame with the maximum number of sites or site phases counted for each year of the time range. Optionally, this dataframe may also have columns for a minimum count and middle-range, probability-weighted count, calculated based on the probability of site existence (siteprobs argument). Equivalent measures of maximum, minimum and probability-weighted total site area can also be calculated if the siteareas argument is specified. If aoristicsum is set to TRUE, then equivalent measures for an aoristic sum of sites can also be added. If a consensus dataset is provided, such as a radiocarbon summeed rpobability distribution, then an attempt is made (for site count, and if present, site area) to achieve a consensus between the settlement data and the consensus dataset. For example, for each site or site phases, it would assign the shape of the consensus dataset intensity for that timespan.
#'
#' @references #' Ratcliffe, J.H. 2000. Aoristic analysis: the spatial interpretation of unspecifed temporal events, International Journal of Geographical Information Science 14: 669–679.x#'
#' Johnson, I. 2004. Aoristic analysis: seeds of a new approach to mapping archaeological distributions through time, In K.F. Ausserer, W. Börner, M. Goriany and L. Karlhuber-Vöckl (eds.), [Enter the past]. The E-way into the Four Dimensions of Cultural Heritage (CAA 2003): 448–452. Oxford: Archaeopress.x#'
#' Palmisano, A., Bevan, A. and S. Shennan 2017. Comparing archaeological proxies for long-term population patterns: An example from central Italy, Journal of Archaeological Science 87: 59-72 448–452.x#'
#' @export
siteCount <- function(startdates, enddates, siteprobs=NA, siteareas=NA, consensus=NA, conprop=0.5, timeRange=NA, calendar="BCAD", aoristic=FALSE, verbose=TRUE){

    ## Initial checks    
    if (!is.na(timeRange[1])){
        if (length(timeRange)!=2){
            stop("timeRange must be either NA (default) or a numeric vector of length 2 specifying the oldest and youngest date in the desired range.")
        }
    }
    if (calendar=="BCAD"){
        years <- seq(min(startdates)+1,max(enddates),1)
        years <- years[years !=0]
        if (!is.na(timeRange[1])){
            years <- years[years >= timeRange[1] & years <= timeRange[2]]
        }
        sitecounts <- data.frame(YearBCAD=years, MaxCount=0)
    } else if (calendar=="BP"){
        years <- seq(max(startdates)-1,min(enddates),1)
        if (!is.na(timeRange[1])){
            years <- years[years <= timeRange[1] & years >= timeRange[2]]
        }
        sitecounts <- data.frame(YearBP=years, MaxCount=0)
    } else {
        stop("calendar argument must be either \"BCAD\" or \"BP\".")
    }
    if (length(startdates) != length(enddates)){
        stop("The arguments startdates and enddates must be of the same length.")
    }
    if (!is.na(siteprobs[1])){
        if (length(startdates) != length(siteprobs)){
            stop("The argument siteprobs must the same length as startdates.")
        }
        if (max(siteprobs) > 1 | min(siteprobs) < 0){
            stop("siteprobs argument must fall within the range [0,1].")
        }
        sitecounts$MinCount <- 0
        sitecounts$MidCount <- 0
    }
    if (!is.na(siteareas[1])){
        if (length(startdates) != length(siteareas)){
            stop("The argument siteprobs must the same length as startdates.")
        }
        if (min(siteareas) < 0){
            stop("siteareas argument must be greater than 0.")
        }
        sitecounts$MaxArea <- 0
        if (!is.na(siteprobs[1])){
            sitecounts$MinArea <- 0
            sitecounts$MidArea <- 0
        }
    }
    if (aoristic){
        sitecounts$MaxAorSum <- 0
        if (!is.na(siteprobs[1])){
            sitecounts$MinAorSum <- 0
            sitecounts$MidAorSum <- 0
        }
    }
    if (any(!is.na(consensus[1]))){
        if (length(conprop)==1){
            conprop <- rep(conprop,nrow(consensus))
        } else {
            if (length(conprop) != nrow(consensus)){
                stop("conprop must be either a single number or a numeric vector of the same length as start dates.")
            }
        }
        if (all(conprop >1 | conprop < 0)){
            stop("conprop must be a fraction between 0 and 1.")
        }
        cvals <- as.data.frame(consensus)
        names(cvals) <- c("Year","Intensity")
        cvals$ConProp <- conprop
        cvals <- merge(sitecounts, cvals, by.x=names(sitecounts)[1], by.y="Year", sort=FALSE)
        cvals <- cvals[,c(names(sitecounts)[1],"Intensity","ConProp")]
        if (nrow(cvals) != nrow(sitecounts)){
            stop("The consensus argument must entirely cover the range of the start and end dates.")
        } else {
            sitecounts$ConsMaxCount <- 0
            if (!is.na(siteprobs[1])){
                sitecounts$ConsMinCount <- 0
                sitecounts$ConsMidCount <- 0
            }
            if (!is.na(siteareas[1])){
                sitecounts$ConsMaxArea <- 0
            }
            if (!is.na(siteprobs[1]) & !is.na(siteareas[1])){
                sitecounts$ConsMinArea <- 0
                sitecounts$ConsMidArea <- 0
            }
        }
    }
    if (length(startdates)>1 & verbose){
        print("Working on main site calculations...")
        flush.console()
        pb <- txtProgressBar(min=1, max=length(startdates), style=3)
    }
    ## Main calculations per site or site phase
    for (a in 1:length(startdates)){
        if (length(startdates)>1 & verbose){ setTxtProgressBar(pb, a) }
        siteyears <- years <= enddates[a] & years > startdates[a]
        sitecounts[siteyears, "MaxCount"] <- sitecounts[siteyears, "MaxCount"] + 1
        if (!is.na(siteprobs[1])){
            if (siteprobs[a]==1){
                sitecounts[siteyears, "MinCount"] <- sitecounts[siteyears, "MinCount"] + 1
            }
            sitecounts[siteyears, "MidCount"] <- sitecounts[siteyears, "MidCount"] + siteprobs[a]
        }
        if (!is.na(siteareas[1])){
            sitecounts[siteyears, "MaxArea"] <- sitecounts[siteyears, "MaxArea"] + siteareas[a]
            if (!is.na(siteprobs[1])){
                if (siteprobs[a]==1){
                    sitecounts[siteyears, "MinArea"] <- sitecounts[siteyears, "MinArea"] + siteareas[a]
                }
                sitecounts[siteyears, "MidArea"] <- sitecounts[siteyears, "MidArea"] + (siteareas[a]*siteprobs[a])
            }
        }
        if (aoristic){
            yearrange <- abs(enddates[a] - startdates[a])
            sitecounts[siteyears, "MaxAorSum"] <- sitecounts[siteyears, "MaxAorSum"] +(1/yearrange)
            if (!is.na(siteprobs[1])){
                if (siteprobs[a]==1){
                    sitecounts[siteyears, "MinAorSum"] <- sitecounts[siteyears, "MinAorSum"] +(1/yearrange)
                }
                sitecounts[siteyears, "MidAorSum"] <- sitecounts[siteyears, "MidAorSum"] + ((1/yearrange) * siteprobs[a])
            }
        }
    }
    if (length(startdates)>1 & verbose){
        close(pb)
    }
    ## Consensus model if required
    if (any(!is.na(consensus[1]))){
        if (length(startdates)>1 & verbose){
            print("Calculating consensus values...")
            flush.console()
            pb <- txtProgressBar(min=1, max=length(startdates), style=3)
        }
        for (b in 1:length(startdates)){
            if (length(startdates)>1 & verbose){ setTxtProgressBar(pb, b) }
            siteyears <- years <= enddates[b] & years > startdates[b]
            sitevals <- cvals[siteyears,]
            consweight <- sum(sitevals$Intensity)/sum(cvals$Intensity)
            siteweight <- sum(sitecounts[siteyears, "MaxCount"])/sum(sitecounts[, "MaxCount"])
            finalweight <- (consweight*sitevals$ConProp)+(siteweight*(1-sitevals$ConProp))
            sitevals1 <- (sitevals$Intensity/sum(sitevals$Intensity))*finalweight
            sitecounts[siteyears, "ConsMaxCount"] <- sitecounts[siteyears, "ConsMaxCount"] + sitevals1
            if (!is.na(siteprobs[1])){
                if (siteprobs[b]==1){
                    siteweight <- sum(sitecounts[siteyears, "MinCount"])/sum(sitecounts[, "MinCount"])
                    finalweight <- (consweight*sitevals$ConProp)+(siteweight*(1-sitevals$ConProp))
                    sitevals1 <- (sitevals$Intensity/sum(sitevals$Intensity))*finalweight
                    sitecounts[siteyears, "ConsMinCount"] <- sitecounts[siteyears, "ConsMinCount"] + sitevals1
                }
                siteweight <- sum(sitecounts[siteyears, "MidCount"])/sum(sitecounts[, "MidCount"])
                finalweight <- (consweight*sitevals$ConProp)+(siteweight*(1-sitevals$ConProp))
                sitevals1 <- (sitevals$Intensity/sum(sitevals$Intensity))*finalweight
                sitecounts[siteyears, "ConsMidCount"] <- sitecounts[siteyears, "ConsMidCount"] + (siteprobs[b]*sitevals1)
            }
            if (!is.na(siteareas[1])){
                siteweight <- sum(sitecounts[siteyears, "MaxArea"])/sum(sitecounts[, "MaxArea"])
                finalweight <- (consweight*sitevals$ConProp)+(siteweight*(1-sitevals$ConProp))
                if (any(sitevals$ConProp==0)){
                    sitevals1 <- (sitecounts[siteyears, "MaxArea"]/sum(sitecounts[siteyears, "MaxArea"]))*finalweight
                } else {
                    sitevals1 <- (sitevals$Intensity/sum(sitevals$Intensity))*finalweight
                }
                sitecounts[siteyears, "ConsMaxArea"] <- sitecounts[siteyears, "ConsMaxArea"] + siteareas[b] *sitevals1
                if (!is.na(siteprobs[1])){
                    if (siteprobs[b]==1){
                        siteweight <- sum(sitecounts[siteyears, "MinArea"])/sum(sitecounts[, "MinArea"])
                        finalweight <- (consweight*sitevals$ConProp)+(siteweight*(1-sitevals$ConProp))
                        if (any(sitevals$ConProp==0)){
                            sitevals1 <- (sitecounts[siteyears, "MinArea"]/sum(sitecounts[siteyears, "MinArea"]))*finalweight
                        } else {
                            sitevals1 <- (sitevals$Intensity/sum(sitevals$Intensity))*finalweight
                        }
                        sitecounts[siteyears, "ConsMinArea"] <- sitecounts[siteyears, "ConsMinArea"] + (siteareas[b]*sitevals1)
                    }
                    siteweight <- sum(sitecounts[siteyears, "MidArea"])/sum(sitecounts[, "MidArea"])
                    finalweight <- (consweight*sitevals$ConProp)+(siteweight*(1-sitevals$ConProp))
                    if (any(sitevals$ConProp==0)){
                        sitevals1 <- (sitecounts[siteyears, "MidArea"]/sum(sitecounts[siteyears, "MidArea"]))*finalweight
                    } else {
                        sitevals1 <- (sitevals$Intensity/sum(sitevals$Intensity))*finalweight
                    }
                    sitecounts[siteyears, "ConsMidArea"] <- sitecounts[siteyears, "ConsMidArea"] + (siteareas[b]*siteprobs[b]*sitevals1)
                }
            }
        }
        if (length(startdates)>1 & verbose){
            close(pb)
        }
    }
    class(sitecounts) <- append(class(sitecounts),"siteCounts")
    if (length(startdates)>1 & verbose){
        print("Done.")
    }
    return(sitecounts)
}

