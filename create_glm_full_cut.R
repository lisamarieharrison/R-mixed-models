#create a cut down matrix that only includes every 4th station
glm.full.cut <- glm.full.sort[1, ]
for(i in 1:nrow(glm.full.sort)){
  if(any(seq(2, 118, , by = 4) == as.numeric(glm.full.sort$stn[i])))
    glm.full.cut <- rbind(glm.full.cut, glm.full.sort[i, ])
}
glm.full.cut <- glm.full.cut[-1, ]
glm.full.cut$profile.depth <- as.factor(glm.full.cut$profile.depth)