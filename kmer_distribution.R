# slide window---------------
## require packages---------
library(dplyr)
library(ggplot2)
dir <- '/Users/user/Desktop/2018_TPC'
groups <- list.files(dir,
                        pattern="*",
                        full.names=FALSE)
strands <- c('anti', 'both', 'sen')

# Create an Empty DataFrame with 0 rows and n columns 
columns = c("bin", "tp_tn", "zscore", "group", "strand", "cut.off") 
com.res = data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(com.res) = columns

for (x in 1:length(groups)) {
  for (y in 1:3) {
    file_name <- list.files(paste0(dir, '/', groups[x], '/', strands[y]),
                            pattern="*_sites.genemap",
                            full.names=FALSE)
    for (z in 1:2) {
      df <- read.delim(paste0(dir, '/', groups[x], '/', strands[y], '/', file_name[z]))
      # clarify genes to positive or negative ones based on Class
      pos <- df %>% filter(Class == 1)
      neg <- df %>% filter(Class == 0)
      
      # positive--------
      pos <- pos[,-2] # remove Class column
      pos.id <- data.frame(pos$ID) # extract positive gene ID
      kmer.name <- names(pos[,-1]) # extract all K-mer name
      
      for (i in 1:length(pos[,-1])) {
        if(i ==1){
          number <- i
          kmer <- kmer.name[number]
          pos.1 <- data.frame(str_split(pos[ ,kmer], ',', n =Inf, simplify = TRUE))
          pos.1 <- cbind(pos.id, pos.1)
          pos.1 <- gather(pos.1, key = 'column', value = 'kmers', -pos.ID, na.rm = TRUE)
          pos.1 <- pos.1[,-2]
          pos.1 <- pos.1 %>% filter(kmers != '')
          names(pos.1) <- c('ID', kmer)
          
          pos.1 <- gather(pos.1, key = 'kmers', value = 'location', -'ID', na.rm = TRUE)
          
        }else{
          number <- i
          kmer <- kmer.name[number]
          pos.2 <- data.frame(str_split(pos[ ,kmer], ',', n =Inf, simplify = TRUE))
          pos.2 <- cbind(pos.id, pos.2)
          pos.2 <- gather(pos.2, key = 'column', value = 'kmers', -pos.ID, na.rm = TRUE)
          pos.2 <- pos.2[,-2]
          pos.2 <- pos.2 %>% filter(kmers != '')
          names(pos.2) <- c('ID', kmer)
          
          pos.2 <- gather(pos.2, key = 'kmers', value = 'location', -'ID', na.rm = TRUE)
          
          pos.1 <- rbind(pos.1, pos.2) # gene-kmer-location long data.frame
        }
        
        dis <- pos.1 %>%
          select(3)
        dis$location <- as.numeric(dis$location) # extract location info
        
        # Basic density
        # give a screen of original K-mer distribution on gene promoter
        p <- ggplot(dis, aes(x=location)) + 
          geom_density() +
          theme_classic()
        p
        # Add mean line
        p +
          geom_vline(aes(xintercept=mean(location)),
                     color="blue", linetype="dashed", linewidth=1) +
          labs(title = expression(bolditalic('positive K-mer distribution'))) +
          scale_x_continuous('Location',
                             limits = c(0,1500),
                             breaks = seq(0, 1500, 500),
                             labels = c('-1000','-500', 'TSS', '500'))
        
        
        # count appearence one each site
        dis <- pos.1 %>%
          select(3) %>% 
          group_by(location) %>% 
          summarise(count = n())
        
        dis$location <- as.numeric(dis$location)
        
        # prepare window
        window.pos <- data.frame(location = 1:1500)
        window.pos <- full_join(window.pos, dis, by = 'location')
        window.pos[is.na(window.pos$count),]$count <- 0
        
        
        
        
        
        
        # negative-----------
        neg <- neg[,-2] # remove Class column
        neg.id <- data.frame(neg$ID) # extract positive gene ID
        kmer.name <- names(neg[,-1]) # extract all K-mer name
        
        
        for (i in 1:length(neg[,-1])) {
          if(i ==1){
            number <- i
            kmer <- kmer.name[number]
            neg.1 <- data.frame(str_split(neg[ ,kmer], ',', n =Inf, simplify = TRUE))
            neg.1 <- cbind(neg.id, neg.1)
            neg.1 <- gather(neg.1, key = 'column', value = 'kmers', -neg.ID, na.rm = TRUE)
            neg.1 <- neg.1[,-2]
            neg.1 <- neg.1 %>% filter(kmers != '')
            names(neg.1) <- c('ID', kmer)
            
            neg.1 <- gather(neg.1, key = 'kmers', value = 'location', -'ID', na.rm = TRUE)
            
          }else{
            number <- i
            kmer <- kmer.name[number]
            neg.2 <- data.frame(str_split(neg[ ,kmer], ',', n =Inf, simplify = TRUE))
            neg.2 <- cbind(neg.id, neg.2)
            neg.2 <- gather(neg.2, key = 'column', value = 'kmers', -neg.ID, na.rm = TRUE)
            neg.2 <- neg.2[,-2]
            neg.2 <- neg.2 %>% filter(kmers != '')
            names(neg.2) <- c('ID', kmer)
            
            neg.2 <- gather(neg.2, key = 'kmers', value = 'location', -'ID', na.rm = TRUE)
            
            neg.1 <- rbind(neg.1, neg.2) # gene-kmer-location long data.frame
          }
          
        }
        
        
        dis2 <- neg.1 %>%
          select(3)
        dis2$location <- as.numeric(dis2$location) # extract location info
        
        
        
        # Basic density
        # give a screen of original K-mer distribution on gene promoter
        p <- ggplot(dis2, aes(x=location)) + 
          geom_density() +
          theme_classic()
        p
        # Add mean line
        p +
          geom_vline(aes(xintercept=mean(location)),
                     color="blue", linetype="dashed", size=1) +
          labs(title = expression(bolditalic('positive K-mer distribution'))) +
          scale_x_continuous('Location',
                             limits = c(0,1500),
                             breaks = seq(0, 1500, 500),
                             labels = c('-1000','-500', 'TSS', '500'))
        
        # count appearence one each site
        dis2 <- neg.1 %>%
          select(3) %>% 
          group_by(location) %>% 
          summarise(count = n())
        
        dis2$location <- as.numeric(dis2$location)
        
        # prepare window
        window.neg <- data.frame(location = 1:1500)
        window.neg <- full_join(window.neg, dis2, by = 'location')
        window.neg[is.na(window.neg$count),]$count <- 0
        
        
        
        
        
        # combine window---------
        colnames(window.pos) <- c('location', 'pos')
        colnames(window.neg) <- c('location', 'neg')
        window <- merge(window.pos, window.neg, by = 'location')
        
        
        
        
        
        # calculate-----------
        # bin:100; sliding window:25; median
        library(zoo)
        pos.win <- rollapply(window$pos, width = 100, median, by = 25, partial = FALSE)
        neg.win <- rollapply(window$neg, width = 100, median, by = 25, partial = FALSE)
        
        
        tp.tn <- pos.win/neg.win
        res <- data.frame(bin = seq(from = 100, to = 1500, by = 25))
        res$tp_tn <- tp.tn
        
        res <- res %>% 
          mutate(zscore = (tp_tn - mean(tp_tn))/sd(tp_tn))
        
        cut_off <- str_sub(file_name[z], start = 1, end = 5)
        res$group <- groups[x]
        res$strand <- strands[y]
        res$cut.off <- cut_off
        
        com.res <- rbind(com.res, res)
      }
    }
  }
}

# first kind
ggplot(com.res, aes(x = bin,
                    y = zscore)) +
  geom_line(aes(colour = group), linewidth = 2)+
  # geom_point(size=2,shape=21) +
  labs(title = expression(italic(bold("K-mers distribution"))),
       subtitle = expression(italic("bin = 100; sliding window = 25; median")),
       x = 'Distance to TSS (bp)',
       y = 'TP/TN (z-score)') +
  scale_x_continuous('Distance to TSS (kb)',
                     limits = c(0,1500),
                     breaks = seq(0, 1500, 500),
                     labels = c('-1.0','-0.5', 'TSS', '0.5')) +
  facet_grid(cut.off ~ group) +
  theme_classic()