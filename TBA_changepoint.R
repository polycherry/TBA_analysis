library(changepoint)
library(tidyverse)
library(dplyr)
library(ggplot2)

#Mention path to folder where CSV files are being exported
path_to_folder = "/Users/polisev/Desktop/TBA/TBAOpt15/"

#Read the CSV files
csv_files = list.files(path = path_to_folder,
                       pattern = "*.csv",
                       full.names = TRUE)

#Output Dataframe
df_summary = data.frame()

#Loop for CSV files one by one
for (files in csv_files) {
  file_name <- substring(tools::file_path_sans_ext(basename(files)),13)
  
  d = read_csv(files)
  
  
  x = d$X
  dapi = d$DAPI
  green = d$Green
  test  = d$Test
  
  
  #y_smooth = filter(butter(3,0.1),y)
  #Butter filter didnt work as well so switched to running median
  
  #Running median filter with window of 81, can be changed
  med_window = 161
  dapi_smooth = runmed(dapi, k = med_window)
  
  #Can set number of layers - works well with 5
  numLayers <- 5
  cpt <- cpt.mean(dapi_smooth, method = "BinSeg", Q = 4)
  changePoints <- cpts(cpt)
  
  test_smooth = runmed(test, k = 81)
  
  #segmenting the test and test_smooth arrays by the 5 layers
  DG_layers = c(
    "Molecular_ventral",
    "Granular_ventral",
    "Hilus",
    "Granular_dorsal",
    "Molecular_dorsal"
  )
  
  #These 2 lists will carry the segmented 647 and 647_smooth values
  segments_test <- list()
  segments_test_smooth <- list()
  
  #These 2 vectors will carry just the mean of each segment
  mean_test = c()
  mean_test_smooth = c()
  median_test = c()
  
  #Discarding the first 150 and last 150 microns of whole and leaving 40 for buffer
  startIdx <- 150
  
  for (i in changePoints) {
    endIdx <- i - 40
    
    segments_test[[length(segments_test) + 1]] <- test[startIdx:endIdx]
    mean_test[length(segments_test)] = mean(segments_test[[length(segments_test)]])
    median_test[length(segments_test)] = median(segments_test[[length(segments_test)]])
    
    segments_test_smooth[[length(segments_test_smooth) + 1]] <- test_smooth[startIdx:endIdx]
    mean_test_smooth[length(segments_test_smooth)] = mean(segments_test_smooth[[length(segments_test_smooth)]])
    
    startIdx <- endIdx + 80
  }
  segments_test[[length(segments_test) + 1]] <- test[startIdx:(length(test) - 150)]
  mean_test[length(segments_test)] = mean(segments_test[[length(segments_test)]])
  median_test[length(segments_test)] = median(segments_test[[length(segments_test)]])
  
  segments_test_smooth[[length(segments_test_smooth) + 1]] <- test_smooth[startIdx:(length(test) - 150)]
  mean_test_smooth[length(segments_test_smooth)] = mean(segments_test_smooth[[length(segments_test_smooth)]])
  
  
  #Calculating max value from all 3 layers
  max_value = max(c(max(green), max(dapi), max(test)))
  
  x1 = seq(1,length(dapi))
  
  # Plot the original and smoothed data
  plot(
    x1,
    dapi,
    ylim = c(350, 40000),
    #can use max_value insted of 40k for plotting
    type = 'l',
    col = 'blue',
    main = file_name,
    xlab = 'Distance across Hippocampus',
    ylab = 'Intensity'
  )
  lines(x1, dapi_smooth, col = 'cyan')
  lines(x1, green, col = "green")
  lines(x1, test, col = "red")
  lines(x1, test_smooth, col = "magenta")
  
  startIdx <- 0
  count = 1
  # Add vertical lines at change points
  for (i in changePoints) {
    abline(v = x1[i],
           col = 'black',
           lty = 2)
    #drawing the buffer zones
    abline(v = x1[i + 40],
           col = 'grey',
           lty = 2)
    abline(v = x1[i - 40],
           col = 'grey',
           lty = 2)
    
    segments(startIdx,
             mean_test_smooth[count],
             x1[i],
             mean_test_smooth[count],
             col = 'brown')
    
    
    count = count + 1
    
    startIdx = x1[i] + 1
    
  }
  
  segments(startIdx,
           mean_test_smooth[count],
           max(x1),
           mean_test_smooth[count],
           col = 'brown')
  
  # Add a legend
  legend(
    'bottomright',
    legend = c(
      'DAPI',
      '488',
      '647 test',
      'smoothed DAPI',
      'smoothed 647',
      'Mean 647',
      'Mean 647-smoothed'
    ),
    col = c('blue', 'green', 'red', 'cyan', 'magenta', 'black', 'brown'),
    lty = 1,
    xpd = TRUE,
    # Allow drawing outside plot region
    cex = 0.7,
    # Smaller legend text
    
  )
  
  #populating dataframe with the mean values
  new_row = c(file_name, mean_test, median_test, mean_test_smooth)
  
  df_summary = rbind(df_summary, as.data.frame(t(new_row), stringsAsFactors = FALSE))
  
}

#Converting chr to numeric
df_summary[, 2:16] <- lapply(df_summary[, 2:16], as.numeric)

#Renaming the summary dataframe
names(df_summary) <- c(
  'slide',
  'mol1mean',
  'gran1mean',
  'hilusmean',
  'gran2mean',
  'mol2mean',
  'mol1median',
  'gran1median',
  'hilusmedian',
  'gran2median',
  'mol2median',
  'mol1momean',
  'gran1momean',
  'hilusmomean',
  'gran2momean',
  'mol2momean'
)

df_condensed = df_summary[, 1, drop = FALSE]
colnames(df_condensed) = "slide"
df_condensed$layer1meanratio = df_summary$mol1mean / df_summary$gran1mean
df_condensed$layer2meanratio = df_summary$mol2mean / df_summary$gran2mean

df_condensed$layer1medianratio = df_summary$mol1median / df_summary$gran1median
df_condensed$layer2medianratio = df_summary$mol2median / df_summary$gran2median

df_condensed$layer1momratio = df_summary$mol1momean / df_summary$gran1momean
df_condensed$layer2momratio = df_summary$mol2momean / df_summary$gran2momean


#Which row is the negative control in
normalize_row = 7

#We normalize based on subtracting the layers background first
normalization_values = unlist(df_summary[normalize_row, -1])
norm_mat = as.matrix(df_summary[, -1])
norm_mat = sweep(norm_mat,
                 MARGIN = 2,
                 STATS = normalization_values,
                 FUN = "-")

#And then dividing by the mean intensity across the negative control
NCAvgmean = mean(unlist(df_summary[normalize_row, 2:6]))
NCAvgmedian =  mean(unlist(df_summary[normalize_row, 7:11]))
NCAvgmomean =  mean(unlist(df_summary[normalize_row, 12:16]))

norm_mat[, 1:5] = norm_mat[, 1:5] / NCAvgmean
norm_mat[, 6:10] = norm_mat[, 6:10] / NCAvgmedian
norm_mat[, 11:15] = norm_mat[, 11:15] / NCAvgmomean

#stitching up this matrix with the slide name
df_normalized = cbind(df_summary$slide, as.data.frame(norm_mat))

#Calculating the ratios layerwise for mean, median, MoMF
df_condensed1 = as.data.frame(df_summary$slide)
colnames(df_condensed1) = "slide"

df_condensed1$layer1meanratio = df_normalized$mol1mean / df_normalized$gran1mean
df_condensed1$layer2meanratio = df_normalized$mol2mean / df_normalized$gran2mean

df_condensed1$layer1medianratio = df_normalized$mol1median / df_normalized$gran1median
df_condensed1$layer2medianratio = df_normalized$mol2median / df_normalized$gran2median

df_condensed1$layer1momratio = df_normalized$mol1momean / df_normalized$gran1momean
df_condensed1$layer2momratio = df_normalized$mol2momean / df_normalized$gran2momean

df_condensed1[is.na(df_condensed1)] <- 0


df_long1 <- df_condensed1 %>%
  pivot_longer(cols = -slide,
               # Keep 'slide' column as is
               names_to = "variable",
               # Create a new column 'variable' for column names
               values_to = "value"          # Create a new column 'value' for values
               ) %>%
               
               mutate(group = case_when(
                 # Create a 'group' variable for pairs of columns
                 variable %in% c("layer1meanratio", "layer2meanratio") ~ "Mean",
                 variable %in% c("layer1medianratio", "layer2medianratio") ~ "Median",
                 variable %in% c("layer1momratio", "layer2momratio") ~ "MoMF"
               ))
               
               # Plot the dot plot using ggplot2
               ggplot(df_long1, aes(x = group, y = value, color = slide)) +
                 geom_point(size = 3, alpha = 0.7) +  # Create dot plot
                 #facet_wrap( ~ group, scales = "free_y") +  # Separate plots for each group
                 theme_minimal() +  # Minimal theme for better visualization
                 labs(
                   title = "Ratio of Mol/Gran ratio in DG for diff avergaing measures",
                   x = "Measure of Central tendency",
                   y = "Ratio Mol/Gran",
                   color = "Slide"
                 ) + scale_y_continuous(limits = c(0, 2.4)) +
                 theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels
          
               #calculating momean (which seems best) for both layers together
               df_momean = as.data.frame(df_summary$slide)
               colnames(df_momean) = "slide"
               df_momean$avg_ratio = (df_normalized$mol1mean+df_normalized$mol2mean)/(df_normalized$gran1mean+df_normalized$gran2mean)
               
               #Calculating absolute difference with subtraction and normalizing it with gran of neg control
               df_subtract = as.data.frame(df_summary$slide)
               colnames(df_subtract) = "slide"
               df_subtract$diff = (df_summary$mol1momean + df_summary$mol2momean - df_summary$gran1momean - df_summary$gran2momean)/2
               df_subtract$norm_diff = df_subtract$diff/df_summary$gran1momean[normalize_row]
               df_subtract$norm_diff2 = (df_subtract$diff*(df_summary$mol1momean + df_summary$mol2momean))/2*df_summary$mol1momean[normalize_row]
               plot(df_subtract$norm_diff)
               plot(df_subtract$norm_diff2)               