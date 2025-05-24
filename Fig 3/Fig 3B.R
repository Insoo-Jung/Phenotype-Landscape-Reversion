##### library #####

library(BoolNet)
library(ggplot2)
library(tidyverse)

##### Figure 3B #####

# Change the file path #
network_folder = "C:/Users/Insoo Jung/Desktop/Revision/PLR_project/Phenotype-Landscape-Reversion/"
RBN_6_1_2_folder = "C:/Users/Insoo Jung/Desktop/RBN_Nets/RBN(6,1,2)/"
RBN_6_1_4_folder = "C:/Users/Insoo Jung/Desktop/RBN_Nets/RBN(6,1,4)/"
RBN_6_5_6_folder = "C:/Users/Insoo Jung/Desktop/RBN_Nets/RBN(6,5,6)/"
RBN_8_1_2_folder = "C:/Users/Insoo Jung/Desktop/RBN_Nets/RBN(8,1,2)/"
RBN_8_1_6_folder = "C:/Users/Insoo Jung/Desktop/RBN_Nets/RBN(8,1,6)/"
RBN_8_7_8_folder = "C:/Users/Insoo Jung/Desktop/RBN_Nets/RBN(8,7,8)/"
RBN_10_1_2_folder = "C:/Users/Insoo Jung/Desktop/RBN_Nets/RBN(10,1,2)/"
RBN_10_1_6_folder = "C:/Users/Insoo Jung/Desktop/RBN_Nets/RBN(10,1,6)/"
RBN_10_9_10_folder = "C:/Users/Insoo Jung/Desktop/RBN_Nets/RBN(10,9,10)/"

set.seed(123)

load(paste0(network_folder, "Fig 3/Exact_answer_RBN.RData"))

Compare_Ranking = function(rank1, data1){
  reordered_data = data1[order(data1$Score, decreasing = FALSE),]
  dist = 0
  for(i in 1:length(rank1)){
    score_rank = which(reordered_data$Score == reordered_data$Score[i])
    ctrl_target = paste0(reordered_data$Node[i], "-", reordered_data$Value[i])
    exact_rank = which(rank1 == ctrl_target)
    dist_vector = c()
    for(j in 1:length(score_rank)){
      for(k in 1:length(exact_rank)){
        dist_vector = c(dist_vector, abs(score_rank[j] - exact_rank[k]))
      }
    }
    dist = dist + min(dist_vector)
  }
  return(dist * 2 / (length(rank1) * (length(rank1) - 1)))
}

RBN_folder = c(RBN_6_1_2_folder,RBN_6_1_4_folder,RBN_6_5_6_folder,RBN_8_1_2_folder,RBN_8_1_6_folder,RBN_8_7_8_folder,RBN_10_1_2_folder,RBN_10_1_6_folder,RBN_10_9_10_folder)
RBN_name = c("RBN(6,1,2)","RBN(6,1,4)","RBN(6,5,6)","RBN(8,1,2)","RBN(8,1,6)","RBN(8,7,8)","RBN(10,1,2)","RBN(10,1,6)","RBN(10,9,10)")

whole_value = c()
for(i in 1:9){
  print(i)
  RBN_path = RBN_folder[i]
  RBN_Networks = list.files(RBN_path, pattern = ".txt")
  data = summary_per_rbn[[RBN_name[i]]]
  for(j in 1:length(RBN_Networks)){
    net = loadNetwork(paste0(RBN_path, RBN_name[i], j, ".txt"))
    attr = Large_Attr(net, 2^length(net$genes))
    scores = Large_ARC_Search(net, attr$Attr_Mat, attr$Attr_Size, as.integer(data$mut_node[which(data$network == paste0(RBN_name[i], j, ".txt"))]), as.integer(data$mut_value[which(data$network == paste0(RBN_name[i], j, ".txt"))]), length(net$genes), 15)
    rank = unlist(data[which(data$network == paste0(RBN_name[i], j, ".txt")), 4:dim(data)[2]])
    names(rank) = c()
    real_rank = c()
    for(k in 1:length(rank)){
      if(str_detect(rank[k], as.character(length(net$genes))) == FALSE){
        real_rank = c(real_rank, rank[k])
      }
    }
    similarity = Compare_Ranking(real_rank, scores[2:dim(scores)[1],])
    whole_value = c(whole_value, similarity)
  }
}
whole_value = 1 - whole_value


accuracy = c(whole_value[1:50], whole_value[101:150], whole_value[151:200], whole_value[251:300], whole_value[301:350], whole_value[401:450])
accuracy[which(accuracy < 0)] = 0
property = rep(c(rep("sparse",50),rep("dense",50)),3)
Fig_3B_data = data.frame(accuracy, property)
save(Fig_3B_data, file = "Fig_3B_data.RData")

ggplot(Fig_3B_data, aes(x = property, y = accuracy, fill = property)) + geom_violin() + geom_boxplot(width = 0.1) + theme_bw() + labs(x = "Network Structural Characteristics") + labs(y = "Target Identification Accuracy") + ylim(c(0,1)) + theme(legend.position = "none")

wilcox.test(subset(Fig_3B_data, property == "sparse")$accuracy, subset(Fig_3B_data, property == "dense")$accuracy)

