# Analysis of Luis's simulation results

# Read in results_general_simulations.csv

ldat <- read.csv(file.choose())
       qlogis(het_mean_OBS),
          pch = '.')
     abline(0,1)
     }
     )

Fdiff <- with(ldat,
              qlogis(Fstp_EXP) - qlogis(Fstp_OBS)
              )

Hdiff <- with(ldat,
              qlogis(het_mean_OBS) - qlogis(het_mean_EXP))

ldat <- data.frame(ldat,Fdiff, Hdiff)
rm(Fdiff,Hdiff)

rb <- rainbow(5)[as.numeric(as.factor(ldat$h))]
symb <- as.numeric(as.factor(ldat$Ne)) 
sseq <- symb
sseq[symb == 3] <- 0
sseq[symb == 4] <- 5

ldat <- data.frame(ldat,rb,sseq)
rm(rb,symb,sseq)

with(ldat,{ 
     plot(qlogis(Fstp_EXP),
          qlogis(Fstp_OBS),
          pch = '.')
    abline(0,1)
    }
    )

ldat <- data.frame(ldat,
                   Hdiff = ldat$het_mean_OBS - ldat$het_mean_EXP,
                   Fdiff = ldat$Fst_OBS - ldat$Fst_EXP)


xvals <- 1/(log(ldat$Ne * ldat$cM))
hfac <- as.factor(ldat$h)
sfac <- as.factor(ldat$s)
ldat <- data.frame(ldat,xvals)

mod1a <- lm(Hdiff ~ xvals*sfac*hfac, data = ldat)

cvals <- rep(1/seq(4,12,length.out = 1000),16)
ch <- rep(rep(c(0,0.1,0.3,0.5), each = 4), each = 1000)
cs <- rep(rep(c(.005,.001,.0005,.0001), 4), each = 1000)
cframe <- data.frame(
  xvals = cvals,
  sfac = as.factor(cs),
  hfac = as.factor(ch),
  cs, 
  ch,
  cr = rainbow(5)[as.numeric(as.factor(ch))]
)


fvalsa <- predict(mod1a, 
                 newdata = cframe)



cframe <- data.frame(cframe,fvalsa)
rm(fvalsa)

par(mfrow = c(2,2))
hr <- range(ldat$Hdiff)
for (i in c(0.005, 0.001, 0.0005, 0.0001)){
  with(subset(ldat,  s == i),
       {plot(log(Ne * cM), Hdiff,
             col = rb,
             pch = sseq,
             main = paste('Selection, ',sprintf(i, fmt = '%#.4f')),
             ylab = 'Deviation in logit(H)',
             ylim = hr
       )
         abline(0,0, lty = 3)
         text(8,1,'h=')
         text(8.5,1,'0.0,',col='blue')
         text(9.2,1,'0.1,',col='green')
         text(9.9,1,'0.3,',col='red')
         text(10.6,1,'0.5')
         for (j in c(0,0.1,0.3,0.5)) with(
           subset(cframe, cs == i & ch == j),{
             lines(1/xvals,fvalsa, col = cr)
           }
         )
       }
  )
}


mod1 <- lm(Fdiff ~ xvals*sfac*hfac, data = ldat)



fvals <- predict(mod1, 
                 newdata = cframe)
cframe <- data.frame(cframe,fvals)

ldat <- data.frame(ldat,fvals)
rm(fvals)

par(mfrow = c(2,2))
for (i in c(0.005, 0.001, 0.0005, 0.0001)){
  with(subset(ldat,  s == i),
       {plot(log(Ne * cM), Fdiff,
             col = rb,
             pch = sseq,
             main = paste('Selection, ',sprintf(i, fmt = '%#.4f')),
             ylim = hr,
             ylab = 'Deviation in logit(Fst)'
       )
         text(8,1,'h=')
         text(8.5,1,'0.0,',col='blue')
         text(9.2,1,'0.1,',col='green')
         text(9.9,1,'0.3,',col='red')
         text(10.6,1,'0.5')
         abline(0,0, lty = 3)
         
       }
  )
  for (j in c(0,0.1,0.3,0.5)) with(
    subset(cframe, cs == i & ch == j),{
    lines(1/xvals,fvals, col = cr)
      }
    )
}

fitmat <- 1:1000
for (i in c(0.005, 0.001, 0.0005, 0.0001)){
  for (j in c(0,0.1,0.3,0.5)){
    fitmat <- cbind(fitmat,
                    subset(cframe, cs == i & ch == j,
                           fvals)
                    )
    }
}

range(fitmat[,-1])

