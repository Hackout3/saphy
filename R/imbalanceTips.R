#return suspicious tips
imbalanceTips<-function(tree,reps=100){
  treemetrics<-imbalanceMetrics(tree)
  if(length(which(as.vector(unlist(treemetrics[1:12]))<threthold(tree,reps=100)[2,]
                  |as.vector(unlist(treemetrics[1:12])>threthold(tree,reps=100)[1,])))>1)
  {
    ImbalanceTips<-NULL
    for(i in 1:tree$Nnode){
      if(!any(unlist(imbalanceMetrics(timeprune(tree)$trees[[i]]))[1:12]<threthold(timeprune(tree)$trees[[i]],reps=100)[2,]|
              unlist(imbalanceMetrics(timeprune(tree)$trees[[i]]))[1:12]>threthold(timeprune(tree)$trees[[i]],reps=100)[1,]))
      {
        ImbalanceTips<-c(ImbalanceTips,timeprune(tree)$trees[[i]]$tip.label[1])
      }
    }
    return(ImbalanceTips)}
  else (return("This tree is balanced"))}
