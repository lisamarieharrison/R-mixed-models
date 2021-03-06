Fluoro GLMs 
========================================================

This is an R Markdown document. Markdown is a simple formatting syntax for authoring web pages (click the **Help** toolbar button for more details on using R Markdown).

When you click the **Knit HTML** button a web page will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

```
setwd(dir = "C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/Data")
dat <- read.csv(file = "thinCTD.csv", header = T)
attach(dat)
library(nlme)


glm.dat <- dat[, c("fluoro", "sal", "temp", "par", "oxy", "stn", "profile.depth")]
glm.dat$stn <- as.factor(glm.dat$stn)
l.fluoro <- log(glm.dat$fluoro)
l.fluoro[is.nan(l.fluoro)] <- NA
glm.dat <- cbind(glm.dat, l.fluoro)
glm.dat$sal[glm.dat$sal == -9] <- NA
glm.dat$temp[glm.dat$temp == -9] <- NA

glm.dat$fluoro[glm.dat$fluoro <= 0] <- NA
glm.dat$l.fluoro[glm.dat$l.fluoro == -Inf] <- NA

fluoro.glm <- glm(fluoro ~ sal + par + temp + oxy + profile.depth, data = na.omit(glm.dat), family = Gamma(link = "inverse"), start = c(1, 1, 1, 1, 1, 1))
summary(fluoro.glm)
plot(fluoro.glm)
```

```r
summary(fluoro.glm)
```

```
## Error: object 'fluoro.glm' not found
```

```r
plot(fluoro.glm)
```

```
## Error: object 'fluoro.glm' not found
```


```
fluoro.lme <- lme(fluoro ~ sal + par + temp + oxy + profile.depth, random =~ stn + profile.depth|stn, data = na.omit(glm.dat))
```

```r
summary(fluoro.lme)
```

```
## Error: object 'fluoro.lme' not found
```

```r
plot(fluoro.lme)
```

```
## Error: object 'fluoro.lme' not found
```


```
glm.dat$l.fluoro <- log(glm.dat$fluoro)
fluoro.lme <- lme(l.fluoro ~ sal + par + temp + oxy + profile.depth, random =~ stn, data = na.omit(glm.dat))
```

```r
summary(fluoro.lme)
```

```
## Error: object 'fluoro.lme' not found
```







