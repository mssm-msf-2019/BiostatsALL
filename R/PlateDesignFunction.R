#' Function to randomize and plot samples on a 96-well plate 
#' @description Given number of unique subjects, visits, replicates, buffers, iterplate controllers (IPC) this function randomizes and plots samples on a 96-well plate 
#' @param n.subjects: total number of unique subjects
#' @param n.visits: planned number of visits (min=1)
#' @param n.replicates: number of repliactes (min=1)
#' @param n.buffers: number of buffers (min=0)
#' @param n.ipc: number of interplate controllers (min=0)
#' @param ht: height of the pdf plot
#' @param wd: width of the pdf plot
#' @param wd: name and directory to save the plot 
#' @examples
#' Image.Randomized.Plates(n.subjects=25, n.visits=2, n.replicates=1, n.buffers=2, n.ipc=1)

Image.Randomized.Plates <- function(n.subjects=25, n.visits=2, n.replicates=2, 
                                    n.buffers=2, n.ipc=1, ht=6.5, wd=8, 
                                    fname="Image.Plate.Randomization.") {
  
  require(reshape2)
  require(magrittr)
  require(plyr)
  require(ggplot2)

  n.wells <- 96
  
  # 1. replicates for the same subject are on the same plate
  # 2. replicates are randomly distributed within a plate
  
  ###### 1. defining parameters ----
  subject.ids <- seq(1, n.subjects)
  replicate.ids <- seq(1, n.replicates)
  visit.ids <- seq(1, n.visits)
  if (n.buffers > 0) {buffer.ids <- paste0("Buffer.", seq(1, n.buffers))} else {buffer.ids <- NULL}
  if (n.ipc > 0) {ipc.ids <- paste0("IPC.", seq(1, n.ipc))} else {ipc.ids <- NULL}
  ur <- sort(as.vector(outer(subject.ids, visit.ids, paste, sep=".")))  # unit of replication (Subject*Visit)
  ur.replicate <- sort(as.vector(outer(ur, replicate.ids, paste, sep=".")))  # unit of replication technical replicates
  
  # number of plates needed
  n.plates<-ceiling((length(subject.ids)*n.visits*n.replicates)/(n.wells-sum(n.buffers,n.ipc)))
  
  # total number of emtry wells
  n.empty <- n.wells*n.plates-(n.subjects*n.visits*n.replicates + n.buffers*n.plates + n.ipc*n.plates)
  empty.per.plate <- floor(n.empty/n.plates)
  
  # max number of subjects on 1 plate given number of replicates, visits, and empty wells 
  max.subjects.plate <- (n.wells-sum(n.buffers,n.ipc,empty.per.plate))%/%(n.replicates*n.visits) + 1
  
  ##### 2. randomly selecting subject ids to be distributed among plates -----
  plates.w.subjects<-list()
  
  for(i in 1:n.plates) {
    
    if (i == 1) {
      used.samples <- NA
      left.samples <- subject.ids
    } else {
      used.samples <- melt(plates.w.subjects)[,1]
      left.samples <- setdiff(subject.ids, used.samples)
    }
    
    if(length(left.samples) > max.subjects.plate) {
      c <- sample(left.samples, max.subjects.plate, replace=FALSE, prob=NULL)
    } else {
      c <- sample(left.samples, length(left.samples), replace=FALSE, prob=NULL)
    }
    
    name <- paste0('Plate.', i)
    tmp <- list(subj.ids = c)
    plates.w.subjects[[name]] <- tmp
    
    rm(tmp,name,c,used.samples,left.samples)
  }
  
  
  ##### 3. expanding subjectid to include visit and replicate ------
  plates.w.ur <- list()
  
  for (i in 1:length(plates.w.subjects)) {
    
    v <- ldply(plates.w.subjects[[i]], data.frame)[,2]
    
    for (j in 1:length(v)) {
      name <- paste0("ur.ids.",j)
      z <- ur.replicate[grep(paste0("^",v[j],"\\."),ur.replicate)]
      tmp <- list(ur = z)
      plates.w.ur[[paste0('Plate.', i)]][name] <- tmp
    }
    
    rm(tmp,i,j,name,z,v)
  }
  
  
  ##### 4. randomizing plates ------
  db <- ldply(plates.w.ur, data.frame)
  plates.random <- list()
  
  for (i in names(plates.w.subjects)) {
    
    # selecting one plate at a time
    v <- c(as.vector(unlist(subset(db, .id==i,select=-c(.id)))), buffer.ids, ipc.ids)
    v <- v[!is.na(v)]
    
    # counting and assigning ids to emplty wells 
    n.empty.plate <- (n.wells-length(v))
    empty.ids <- paste0("None.", seq(1, n.empty.plate))  # wells that have nothig in them
    
    # randomizing samples on a plate 
    v.s <- sample(c(v,empty.ids))
    m <- matrix(v.s, nrow=8, ncol=12)
    
    plates.random[[i]] <- m
    
  }
  
  
  all.colors<-c('aliceblue', 'antiquewhite2', 'aquamarine', 'lightblue4', 'thistle2', 
                'slategray2','yellowgreen', 'darkorange2', 'lightsalmon1', 'violetred3',
                'orchid4','royalblue3','palevioletred','orangered3','olivedrab4',
                'lightpink4','lightslateblue','yellow2','thistle3','magenta3',
                'orange3','firebrick2','darkslategray','darkseagreen1','firebrick4',
                'darkgreen')
  
  if (max.subjects.plate > 25) {
    all.colors <- setdiff(colors()[grep("1|4",colors())],
                          colors()[grep("gray|grey|peach|beige|goldenrod|plum|seashell|steelblue|wheat
                                        |navajowhite|thistle|lavenderblush|chartreuse|sienna|antiquewhite
                                        |darkolivegreen|wheat|antiquewhite|chocolate|turquoise|lightsalmon
                                        |maroon|maroon|coral|darkorchid|rosybrown|palegreen|darkseagreen
                                        |lightyellow|lightyellow|^red|azure|lightskyblue|ivory|red",colors())])
  }
  
  
  all.colors.s <- sample(all.colors, size=max.subjects.plate+1)
  
  
  pdf(paste0(fname,Sys.Date(),".pdf"), onefile = T, width = wd, height = ht)
  for (i in 1:length(plates.random)) {
    
    extra.colors<-c("gray25","gray50","#000000")
    names(extra.colors)<-c("Buffer","IPC","..")
    
    db <- as.data.frame(plates.random[[i]])
    db <- data.frame(lapply(db, as.character), stringsAsFactors=FALSE)
    db$LETTER <- factor(LETTERS[1:8], levels=LETTERS[8:1])
    image.db <- reshape2::melt(db, id="LETTER") %>% 
      mutate(Number = factor(gsub("[[:alpha:]]", "", variable), levels=1:12),
             Patient = sapply(strsplit(as.character(value), ".", fixed=T), "[", 1),
             Visit = sapply(strsplit(value, ".", fixed=T), "[", 2),
             Replicate = sapply(strsplit(value, ".", fixed=T), "[", 3),
             Label2 = ifelse(!Patient%in%c("None","IPC","Buffer"),
                             paste0("vs", Visit, ".rp", Replicate), Visit))
    
    pat.ord <- as.character(sort(as.numeric(setdiff(unique(image.db$Patient), c("Buffer","IPC","None")))))
    extra.ord <- intersect(c("Buffer","IPC","None"),unique(image.db$Patient) )
    
    image.db$Label2 <- ifelse(image.db$Patient=="None", " ", image.db$Label2)
    image.db$Patient <- factor(mapvalues(image.db$Patient, "None", ".."), levels=c(pat.ord,mapvalues(extra.ord,"None","..")))
    
    pat.colors <- sample(all.colors.s, size=length(pat.ord))
    pat.colors <- as.vector(c(pat.colors, extra.colors[mapvalues(extra.ord,"None","..")]))
    
    p <- ggplot(image.db,aes(x = Number,y = LETTER)) +
      geom_tile(aes(fill = Patient),colour = "black") + 
      geom_text(aes(label = Patient, fontface = "bold", vjust = 0)) +
      geom_text(aes(label = Label2, vjust = 2)) +
      theme_bw() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "bottom") + 
      scale_fill_manual(values = alpha(pat.colors, 0.9), breaks=c(pat.ord,extra.ord),
                        guide = guide_legend(nrow=5)) +
      labs(title = paste0('Plate # ',i, ". N subjects = ", length(pat.ord)),
           x = "", y = "") 
    
    plot(p)
    
  }
  dev.off()
  
}

