

#=============================================================
#=============================================================


meand.ipr <- function(phi = NULL, u = NULL, v = NULL, obs.start, obs.end = F, obs.length = F,
                      wd.start = F, wd.end = F, wd.length = F, dt.steps, percentile = .05,
                      plot = F, path = "./", break.step = 15, output = F){
  
  # This function calculates percentiles of wind direction changes within
  # defiend time windows to allow an estimation of meandering. To also
  # consider slower direction shifts, extending over several observations,
  # the maximum shift for all time scales from input resolution up to
  # input resolution * 'dt.steps' is considered and their maximum used.
  # For each time window, defined by 'wd.start/-end/-length' the distribution
  # of these direction shifts is used to fit a Students t distribution, of
  # which the desired percentiles are calculated. If the fit should fail,
  # the normal distribution will be used instead, which will be prompted in
  # the console.
  # Output (which can also be written into a file, see 'output') contains
  # the following values for each window:
  # start and end of the window ('wd.start/-end'), the lower and upper
  # percentile ('min' and 'max'), the mean of the absolute percentiles
  # ('abs_mean') and the absolute maximum within the window ('abs_max').
  #
  # The following packages must be installed:
  # - MASS (install.packages("MASS"))
  # - plyr (install.packages("plyr"))
  # - ncdf4 (install.packages("ncdf4")) (If output as NetCDF desired)
  
  ##################
  ### Parameters ###
  ##################
  #
  # phi:        - vector of wind direction in degree with same
  #               length as 'time.start'. Instead, 'u' and 'v'
  #               can be defined
  # u:          - wind from east, ignored, when 'phi' defined
  # v:          - wind from north, ignored, when 'phi' defined
  # obs.start:  - time vector in POSIXct format with same length
  #               as 'phi', defining the start of each observation
  #               interval of 'phi'
  # obs.end:    - same as 'time.start', defining the end of each
  #               observation interval. If FALSE, next 'time.start'
  #               is used.
  # obs.length: - used, if no 'obs.end' defined. Useful, if end of
  #               an observation varies from start of the next
  # wd.start:   - time vector in POSIXct format, defining the start
  #               of each window to estimate meandering. If not
  #               defined, 'wd.length' must be given. In that case,
  #               first window will start from beginning of 'obs.start',
  #               second from there + 'wd.length' and so on
  # wd.end:     - as 'wd.start', defining the end and hence length of
  #               each window. If not defined, calculated by 'wd.length'
  #               if available. Else, next 'wd.start' will define the
  #               end of a window.
  # wd.length:  - as described in 'wd.start' and 'wd.end',
  #               time in second
  # dt.steps:   - number of Delta times to include in dphi/dt,
  #               starting at input resolution.
  #               Example: If input data has 2 min resoltion and
  #               'dt.steps' is 3, a Delta phi for every time step
  #               will be computed for Delta time = 2 min (1* input
  #               resolution), 4 min (2* nput resolution) and 6 min
  #               ('dt.steps' = 3* input resolution). Of these three
  #               Delta phi the largest will be selected for each time
  #               step for further use
  # percentile: - fraction to be cut from both ends of the distribution
  #               to obtain the interpercentile range. 0.05 corresponds
  #               to 5%, hence, the 5th and 95th percentiles would be
  #               used
  # plot:       - if TRUE, the distribution of Delta phi including the
  #               fitted curves, fit parameters and ipr values is
  #               plotted for each window and stored in a new
  #               directory in 'path' (if defined, else in working
  #               directory)
  # path:       - defines, where a new directory for plots and/or
  #               result files will be created, if any of these are
  #               selected. Has to end with a '/'. If not defined,
  #               directory will be created in the current working
  #               directory.
  # break.step: - width of bars (in degree) in histograms, if
  #               'plot' is TRUE
  # output:     - if, FALSE, no output written t file. Options:
  #               'csv' for .csv, 'nc' for NetCDF. Destination:
  #               'path'. If not defined, working directory.
  #
  
  
  ### manage wind direction
  if(is.null(phi)){
    if(any(is.null(u), is.null(v))){
      stop("either 'phi' or 'u' and 'v' must be defined")
    }else{
      phi <- WindVector(wind_east = u, wind_north = v)[,2]
    }
  }
  
  ### manage observation time
  {
    if(obs.end[1] == F){
      if(obs.length == F){
        # if no observation end or length defineed, use time step of
        # 'obs.start' as observation length
        if(length(table(diff(obs.start))) > 1)
          stop("'obs.start' has no regular time step -> 'obs.end' or 'obs.length' must be defined")
        obs.end <- obs.start + diff(obs.start)[1]
      }else{
        obs.end <- obs.start + obs.length
      }
      
    }
  }
  
  ### manage window time
  {
    # ERROR if neither window start nor length are defined
    if(wd.start[1] == F & wd.length == F)
      stop("at least one of 'wd.start' and 'wd.length' must be defined")
    
    if(wd.start[1] == F & is.null(running.wd.center)){
      # if no window start defined, start at beginning of 'obs.start'
      wd.start <- seq(obs.start[1], obs.start[length(obs.start)], wd.length)
      wd.end <- wd.start + wd.length
    }else{
      if(wd.end[1] == F){
        if(wd.length == F){
          # if no window end or length defineed, use time step of
          # 'wd.start' as window length
          if(length(table(diff(wd.start))) > 1)
            stop("'wd.start' has no regular time step -> 'wd.end' or 'wd.length' must be defined")
          wd.end <- wd.start + diff(wd.start)[1]
        }else{
          wd.end <- wd.start + wd.length
        }
      }
    }
  }
  
  t_ipr <- data.frame(wd.start, wd.end)
  t_ipr$abs_max <- t_ipr$abs_mean <- t_ipr$max <- t_ipr$min <- NA
  
  # load required packages
  library(MASS)
  library(plyr)
  
  ##################################
  ### calculate and correct dphi ###
  ##################################
  
  # create matrix to store dphi/dt for all dt values
  diff <- matrix(nrow = length(phi), ncol = dt.steps)
  # iterate over dt steps from raw resolution until
  # raw resolution * 'dt.steps'
  for(i in 1:dt.steps){
    st <- floor(i/2)
    en <- ceiling(i/2)
    # iterate over time raw steps
    for(j in (1+st):(length(phi)-en)){
      # if |dphi| > 180°
      if(abs(phi[j-st] - phi[j+en]) > 180){
        # check if any of the respective steps on raw resolution
        # has larger absolute value than 180° (only applies for
        # dt > raw resolution)
        if(max(abs(diff(phi[c((j-st):(j+en))]))) > 180){
          # if yes, check if its < -180° or > 180° and
          # compensate
          d <- phi[j+en] - phi[j-st]
          if(d < -180){
            diff[j,i] <- (phi[j+en]+360) - phi[j-st]
          }else if(d > 180){
            diff[j,i] <- phi[j+en] - (phi[j-st]+360)
          }
          # else, there is no single |jump| > 180° -> keep value 
        }else diff[j,i] <- phi[j+en] - phi[j-st]
        # else, no need to adjust
      }else diff[j,i] <- phi[j+en] - phi[j-st]
    }
  }; rm(i,j,st,en,d)
  
  
  ############################
  ### iterate over windows ###
  ############################
  
  ### create new directory to contain plots
  if(plot | output != F)
    dir.create(paste0(path,"distr/",dt.steps*2,"min_",percentile), showWarnings = F, recursive = T)
  # dummy vector to store failed fits
  excl <- c(); f <- T
  
  # start progress bar
  cat("compute IPR\n")
  pb <- txtProgressBar(min = 0, max = length(wd.start), style = 3, width = 75, char = "=")
  
  for(i in seq_along(wd.start)){
    ### subset for current window
    dphi <- data.frame(diff[which(obs.start >= wd.start[i] & obs.end <= wd.end[i]),])
    # update progress bar
    setTxtProgressBar(pb, i)
    ### remove all values that overlap with the previous/next window
    for(j in 1:dt.steps){
      st <- floor(j/2)
      en <- ceiling(j/2)
      dphi[c((0:0+st),((nrow(dphi)+1-en):(nrow(dphi)))),j] <- NA
    }
    ### select positive or negative maximum of all dts
    x <- c()
    for(j in 1:nrow(dphi)){
      tryCatch({
        x <- c(x, dphi[j,which(abs(dphi[j,]) == max(abs(dphi[j,]), na.rm = T))[1]])
      }, warning = function(w){}, error = function(e){})
    }
    dphi <- as.vector(na.omit(x))
    rm(j,x)
    
    t_ipr$abs_max[i] <- max(abs(dphi), na.rm = T)
    
    if(max(abs(dphi)) <= 180){ xlim <- c(-180,180)
    }else{ xlim <- c(-round_any(max(abs(dphi), na.rm = T), break.step, ceiling), round_any(max(abs(dphi), na.rm = T), break.step, ceiling))}
    
    if(plot){
      png(paste0(path,"distr/",dt.steps*2,"min_",percentile,"/hist_fit_",i,".png"), width = 16, height = 12, res = 400, units = "cm")
      par(mar = c(3.5,4,1.2,4), oma = c(0,0,0,0), mfrow = c(1,1))
      hist(dphi, breaks = seq(round_any(min(dphi, na.rm = T), break.step, floor), round_any(max(dphi, na.rm = T), break.step, ceiling), break.step), xlim = xlim, freq = F, main = "", xaxt = "n", xlab = "", ylab = "", las = 1)
      axis(side = 1, at = seq(-360,360,45))
      mtext(side = 1, line = 2.3, expression(paste(Delta, phi, " (°)")))
      mtext(side = 2, line = 3, "density")
      lines(density(dphi, na.rm = T), lwd = 2, col = "orange")
      curve(dnorm(x, mean=mean(dphi), sd=sd(dphi)), add=TRUE, lwd = 2, col = "darkgreen")
    }
    
    ###==================###
    ### fit distribution ###
    ###==================###
    
    # in case optimization should fail, the next section, where
    # normal distribution is used to calculate ipr, won't be skipped
    skip <- F
    
    tryCatch({
      ### fit distribution
      fail <- T
      k <- 1
      dphi <- na.omit(dphi)
      # try to fit t distribution with start parameter for df (degrees of
      # freedom) = k = 1. If successfull, 'fail' is overwritten. Else, try again
      # with df = current k + 2. Continue until successfull or k > 500. If latter,
      # use normal distribution instead.
      while(fail){
        tryCatch({
          # fit t distribution to data
          par <- fitdistr(dphi, "t", start = list(m = mean(dphi), s = sd(dphi), df = k))
          fail <- F
        }, error = function(e){})
        k <- k+2
        # if optimization fails continously, stop here and jump to
        # where parameters are calculated based on normal distribution
        if(k > 500){
          if(f){cat("fitting t distribution failed for window: "); f <- F}
          stop("fit not successfull")
        }
      }
      
      ### following bracket only executed when fit successfull
      {
        # if fit successfull -> skips the section, where normal distribution
        # is be used to calulating the ipr
        skip <- T
        
        # fit-parameters of t distribution
        m <- par$estimate[1] # location
        s <- par$estimate[2] # scale
        df <- par$estimate[3] # degrees of freedom
        
        ### calculate cumulative density of the t function
        t_cum <- pt((seq(xlim[1],xlim[2],1) - m)/s, df)
        
        ### calculate interpercentile range
        # all values within range
        ipr <- which(t_cum >= percentile & t_cum <= (1-percentile))
        # range od values
        ra <- range(seq(xlim[1],xlim[2],1)[ipr])
        
        if(plot){
          # plot t function
          curve(dt((x - m)/s, df)/s, add=TRUE, col = "darkblue", lwd = 2)
          
          par(new = T)
          ### plot cumulative density
          plot(seq(xlim[1],xlim[2],1), t_cum, type = "l", lty = 2, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "darkblue", lwd = 2, xlim = xlim)
          axis(side = 4, las = 1)
          mtext(side = 4, line = 2.8, "cumulative density")
          # plot line at ipr range
          abline(v = ra, lwd = 2, col = "grey30")
          
          legend("right", inset = c(.03,0), c("t dist params.:", paste0("location = ",round(m,2)), paste0("scale   =  ", round(s,2)), paste0("freed.  =  ",round(df,2)), paste0("ipr  = ",ra[1],"° to ",ra[2],"°")), bty = "n", text.font = c(2,1,1,1,1))
          legend("topleft", inset = c(.05,0), c("density", "normal", "t", "t cumul."), lwd = 2, col = c("orange", "darkgreen", "darkblue", "darkblue"), lty = c(1,1,1,2), bty = "n")
        }
      }
      # print each 'i' where fit failed in the console
    }, error = function(e){cat(paste0(i,","))})
    
    ### if optimization failed:
    ###========================
    if(!skip){
      ### calculate ipr based on normal distribution
      ipr <- which(pnorm(seq(-180,180,1), mean=mean(dphi, na.rm = T), sd=sd(dphi, na.rm = T)) >= percentile & pnorm(seq(-180,180,1), mean=mean(dphi), sd=sd(dphi)) <= (1-percentile))
      ra <- range(seq(-180,180,1)[ipr])
      
      if(plot){
        par(new = T)
        ### plot cumulative density
        plot(seq(xlim[1],xlim[2],1), pnorm(seq(xlim[1],xlim[2],1), mean=mean(dphi, na.rm = T), sd=sd(dphi, na.rm = T)), type = "l", lty = 2, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "darkgreen", lwd = 2, xlim = xlim)
        axis(side = 4, las = 1)
        mtext(side = 4, line = 2.8, "cumulative density")
        # plot line at ipr range
        abline(v = ra, lwd = 2, col = "grey30")
        
        legend("right", inset = c(.03,0), c("norm dist params.:", paste0("mean = ",round(mean(dphi),2)), paste0("sd   =  ", round(sd(dphi),2)), paste0("ipr  = ",ra[1],"° to ",ra[2],"°")), bty = "n", text.font = c(2,1,1,1))
        legend("topleft", inset = c(.05,0), c("density", "normal", "norm cumul."), lwd = 2, col = c("orange", "darkgreen", "darkgreen"), lty = c(1,1,2), bty = "n")
      }
      # write to vector with failed fits
      # (currently not used/not in output)
      excl <- c(excl,i)
    }
    
    ## write to data frame
    t_ipr$min[i] <- ra[1]
    t_ipr$max[i] <- ra[2]
    t_ipr$abs_mean[i] <- mean(abs(ra))
    
    if(plot){
      legend("left", inset = c(-0.02,0), c(paste0("dt = ",dt.steps*2," min"), paste0("thresh. = ",percentile*100,"%")), bty = "n")
      mtext(side = 3, paste0(wd.start[i], " to ", wd.end[i]))
      
      box()
      dev.off()
    }
  }
  
  ####################
  ### write output ###
  ####################
  
  if(output == "csv"){
    write.csv(t_ipr, file = paste0(path,"distr/",dt.steps*2,"min_",percentile,"/results.csv"), row.names = F)
  }else if(output == "nc"){
    library(ncdf4)
    # create dimensions
    origin <- wd.start[1]
    
    time.start <- ncdim_def(name = "wd.start", vals = as.numeric(difftime(wd.start, wd.start[1], units = "secs")), units = paste0("seconds since ", origin), calendar = "proleptic_gregorian", create_dimvar = T)
    
    # create variables
    time.end <- ncvar_def(name = "wd.end", units = paste0("seconds since ", origin), dim = time.start, prec = "double")
    min <- ncvar_def(name = "min", units = "degree", dim = time.start, prec = "double")
    max <- ncvar_def(name = "max", units = "degree", dim = time.start, prec = "double")
    abs_mean <- ncvar_def(name = "abs_mean", units = "degree", dim = time.start, prec = "double")
    abs_max <- ncvar_def(name = "abs_max", units = "degree", dim = time.start, prec = "double")
    
    # create netCDF
    nc <- nc_create(filename = paste0(path,"distr/",dt.steps*2,"min_",percentile,"/results.nc"),
                    vars = list(time.end, min, max, abs_mean, abs_max),
                    force_v4 = T)
    
    # enter values
    ncvar_put(nc = nc, varid = time.end, vals = as.numeric(difftime(wd.end, wd.start[1], units = "secs")))
    ncvar_put(nc = nc, varid = min, vals = t_ipr$min)
    ncvar_put(nc = nc, varid = max, vals = t_ipr$max)
    ncvar_put(nc = nc, varid = abs_mean, vals = t_ipr$abs_mean)
    ncvar_put(nc = nc, varid = abs_max, vals = t_ipr$abs_max)
    
    # add attributes
    ncatt_put(nc = nc, varid = "wd.start", attname = "window length", attval = paste0(as.numeric(difftime(wd.end[1],wd.start[1], units = "mins"))," min"))
    
    nc_close(nc)
  }else if(output != F) warning("no output written, must be 'csv' or 'nc'")
  
  return(t_ipr)
  
  # close progress bar
  close(pb)
  
}# end function meand.ipr()


#=============================================================
#=============================================================


ipr.categ <- function(ipr, categ, start = NULL, end = NULL, tz = "UTC", plot = T, path = "./", phi = NULL, u = NULL, v = NULL, obs.time = NULL, output = F){
  
  ##################
  ### Parameters ###
  ##################
  #
  # ipr:        - data.frame, being output of 'meand.ipr()'
  # categ:      - numeric or vector of numerics, indicating the
  #               thresholds between categories in degree (must
  #               be increasing)
  # start:      - start time to analyze/plot in format
  #               YYYY-mm-dd HH:MM:SS
  # end:        - as 'start' for the end of period to analyze
  # tz:         - time zone of 'start' and 'end' input
  # plot:       - if TRUE, the distribution of Delta phi including the
  #               fitted curves, fit parameters and ipr values is
  #               plotted for each window and stored in a new
  #               directory in 'path' (if defined, else in working
  #               directory)
  # path:       - defines, where a new directory for plots and/or
  #               result files will be created, if any of these are
  #               selected. Has to end with a '/'. If not defined,
  #               directory will be created in the current working
  #               directory.
  # phi:        - vector of wind direction in degree with same
  #               length as 'time.start'. Instead, 'u' and 'v'
  #               can be defined
  # u:          - wind from east, ignored, when 'phi' defined
  # v:          - wind from north, ignored, when 'phi' defined
  # obs.time:   - time vector in POSIXct format with same length
  #               as 'phi', defining the start of each observation
  #               interval of 'phi'
  # output:     - if, FALSE, no output written t file. Options:
  #               'csv' for .csv, 'nc' for NetCDF. Destination:
  #               'path'. If not defined, working directory.
  #
  
  if(is.null(start)) start <- ipr$wd.start[1]
  if(is.null(end)) end <- ipr$wd.end[nrow(ipr)]
  start <- as.POSIXct(start, tz = tz)
  end <- as.POSIXct(end, tz = tz)
  
  ipr <- ipr[which(ipr$wd.start >= start & ipr$wd.end <= end),]
  
  ### assign categories according to input
  ### (and create legend for plotting)
  ipr$cat <- NA
  n_cats <- length(categ)
  legend <- rep(NA, n_cats+1)
  if(n_cats == 1){
    ipr$cat[which(ipr$abs_mean <= categ)] <- 1
    ipr$cat[which(ipr$abs_mean > categ)] <- 2
    legend <- c(paste0("ipr \u2264 ",categ),paste0("ipr > ",categ))
  }else{
    ipr$cat[which(ipr$abs_mean <= categ[1])] <- 1
    ipr$cat[which(ipr$abs_mean > categ[n_cats])] <- n_cats + 1
    legend[1] <- paste0("ipr \u2264 ",categ[1]); legend[n_cats+1] <- paste0("ipr > ",categ[length(categ)])
    for(c in 2:n_cats){
      ipr$cat[which(ipr$abs_mean > categ[c-1] & ipr$abs_mean <= categ[c])] <- c
      legend[c] <- paste0(categ[c-1]," \u2264 ipr < ",categ[c])
    }
  }
  
  ### create new directory to contain plots and/or output
  if(plot | output != F){
    nam <- "categ"
    for(c in seq_along(categ)){nam <- paste0(nam,"_", categ[c])}
    dir.create(paste0(path,"distr/",nam), showWarnings = F, recursive = T)
  }
  
  ### plot result as time series
  if(plot){
    
    ### manage wind direction
    if(is.null(phi)){
      if(any(is.null(u), is.null(v))){
        stop("no plot was made. Either 'phi' or 'u' and 'v' must be defined")
      }else{
        phi <- WindVector(wind_east = u, wind_north = v)[,2]
      }
    }
    if(is.null(obs.time)) stop("no plot was made since 'obs.time' is not defined")
    
    library(viridis)
    col <- rev(inferno(n_cats+3))
    
    {
      png(paste0(path,"distr/",nam,"/meand_ipr_",substr(start,1,10),"-",substr(end,1,10),".png"), width = 22, height = 14, res = 400, units = "cm")
      
      par(mar = c(.5,4,.2,.5), mfrow = c(2,1), oma = c(1.5,0,1.2,0))
      
      ## wind dir ##
      plot(NULL, NULL, xlab = "", ylab = "wind direction (°)", las = 1, xlim = c(start,end), ylim = c(0,360), xaxt = "n", yaxt = "n", xaxs = "i")
      axis.POSIXct(side = 1, at = c(ipr$wd.start, ipr$wd.end[nrow(ipr)]), format = "%H:%M" , labels = F)
      axis(side = 2, seq(0,360,45), labels = F)
      axis(side = 2, seq(0,360,90), las = 1)
      # plot categories as background color
      for(c in 1:(n_cats+1)){
        cat <- c
        for(i in which(ipr$cat == c)){
          polygon(x = c(ipr$wd.start[i],ipr$wd.end[i],ipr$wd.end[i],ipr$wd.start[i]), y = c(400,400,-50,-50), col = col[c], border = F)}
      }
      abline(h = seq(0,360,45), lwd = .8, col = "grey60")
      abline(v = c(ipr$wd.start, ipr$wd.end), lwd = .8, col = "grey60")
      lines(obs.time, phi, col = "black", lwd = 2)
      #mtext(side = 3, line = .2, paste0("ipr <= ",thr_1,"°: n = ",length(con)," | ",thr_1," > ipr >= ",thr_2,": n = ",length(med_1)," | ",thr_2," > ipr > ",thr_3,": n = ",length(med_2)," |  ipr >= ",thr_3,"°: n = ",length(con)," | dt = ",dt*2," min | pr = ",pr), adj = 0, cex = .8)
      mtext(side = 3, line = .2, paste0(format(start, "%Y-%m-%d %H:%M")," to ",format(end, "%Y-%m-%d %H:%M")), adj = 1, cex = .8)
      box()
      
      ## ipr ##
      plot(NULL, NULL, xlab = "", ylab = expression(paste("ipr ",Delta,phi," (°)")), las = 1, ylim = c(0,max(ipr$abs_mean)), xlim = c(start,end), xaxt = "n", yaxt = "n", xaxs = "i")
      axis.POSIXct(side = 1, at = c(ipr$wd.start, ipr$wd.end[nrow(ipr)]), format = "%H:%M")
      axis(side = 2, seq(0,360,15), las = 1)
      # plot categories as background color
      for(c in 1:(n_cats+1)){
        cat <- c
        for(i in which(ipr$cat == c)){
          polygon(x = c(ipr$wd.start[i],ipr$wd.end[i],ipr$wd.end[i],ipr$wd.start[i]), y = c(400,400,-50,-50), col = col[c], border = F)}
      }
      abline(h = seq(0,360,15), lwd = .8, col = "grey60")
      abline(v = c(ipr$wd.start, ipr$wd.end), lwd = .8, col = "grey60")
      abline(h = categ, lwd = 2, lty = 2)
      lines(ipr$wd.start+(ipr$wd.end-ipr$wd.start)/2, ipr$abs_mean, lwd = 2, type = "b")
      legend("topright", fill = col[1:(n_cats+1)], legend = legend, box.lty = 0, bg = adjustcolor("white", alpha.f = .4))
      box()
      
      dev.off()
    }# end plot
  }# end if(plot)
  
  ### output
  if(output == "csv"){
    write.csv(ipr, file = paste0(path,"distr/",nam,"/meand_ipr_",start,"-",end,".csv"), row.names = F)
  }else if(output == "nc"){
    library(ncdf4)
    # create dimensions
    origin <- ipr$wd.start[1]
    
    time.start <- ncdim_def(name = "wd.start", vals = as.numeric(difftime(ipr$wd.start, ipr$wd.start[1], units = "secs")), units = paste0("seconds since ", origin), calendar = "proleptic_gregorian", create_dimvar = T)
    
    # create variables
    time.end <- ncvar_def(name = "wd.end", units = paste0("seconds since ", origin), dim = time.start, prec = "double")
    min <- ncvar_def(name = "min", units = "degree", dim = time.start, prec = "double")
    max <- ncvar_def(name = "max", units = "degree", dim = time.start, prec = "double")
    abs_mean <- ncvar_def(name = "abs_mean", units = "degree", dim = time.start, prec = "double")
    abs_max <- ncvar_def(name = "abs_max", units = "degree", dim = time.start, prec = "double")
    cat <- ncvar_def(name = "cat", units = paste0("category 1-",n_cats+1), dim = time.start, prec = "double")
    
    # create netCDF
    nc <- nc_create(filename = paste0(path,"distr/",nam,"/meand_ipr_",start,"-",end,".nc"),
                    vars = list(time.end, min, max, abs_mean, abs_max, cat),
                    force_v4 = T)
    
    # enter values
    ncvar_put(nc = nc, varid = time.end, vals = as.numeric(difftime(ipr$wd.end, ipr$wd.start[1], units = "secs")))
    ncvar_put(nc = nc, varid = min, vals = ipr$min)
    ncvar_put(nc = nc, varid = max, vals = ipr$max)
    ncvar_put(nc = nc, varid = abs_mean, vals = ipr$abs_mean)
    ncvar_put(nc = nc, varid = abs_max, vals = ipr$abs_max)
    ncvar_put(nc = nc, varid = cat, vals = ipr$cat)
    
    # add attributes
    ncatt_put(nc = nc, varid = "wd.start", attname = "window length", attval = paste0(as.numeric(difftime(ipr$wd.end[1],ipr$wd.start[1], units = "mins"))," min"))
    
    nc_close(nc)
  }else if(output != F) warning("no output written, must be 'csv' or 'nc'")
  
  return(ipr)
  
}# end function ipr.categ()


#=============================================================
#=============================================================


# by Wolfgang Babel
# compute the mean wind velocity and vector wind direction from the pair:
# wind from east and wind from north
WindVector = function(wind_east,wind_north = NULL){
  # wind_east, x: wind FROM east
  # wind_north, y: wind FROM north
  if (length(wind_north) == 0){
    wind_north = wind_east[,2]
    wind_east  = wind_east[,1]
  }
  dir = rep(NA, length(wind_east))
  for (i in 1:length(wind_east)){
    x = wind_east[i]
    y = wind_north[i]
    if (!is.na(x) & !is.na(y)){
      if (y == 0) {
        if (x > 0) dir[i] = 90 else 
          if (x < 0) dir[i] = 270 
      } else {
        dir[i] = (atan2(x,y) * 180 / pi) %% 360
      }
    }
  }
  umean = sqrt(wind_east^2 + wind_north^2)
  out = data.frame("WindSpd" = umean,
                   "WindDir" = dir)
  return(out)
}


#=============================================================
#=============================================================