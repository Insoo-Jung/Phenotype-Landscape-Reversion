##### library #####

library(BoolNet)
library(ggplot2)
library(readxl)
library(dplyr)

##### Figure 3D #####

# Change the file path #
network_folder = "C:/Users/Insoo Jung/Desktop/network_insoo/"

set.seed(123)

load("Fig 3/network_controllers.RData")

Score_Alteration = c()
Score_LDOI = c()
Score_GS = c()
Score_PLR = c()
for(i in 1:dim(network_controllers)[1]){
  print(i)
  net = loadNetwork(paste0(network_folder, network_controllers$O_net[1]))
  N_Num = length(net$genes)
  Adj = matrix(0, length(net$genes), length(net$genes))
  for(j in 1:length(net$genes)){
    inputs = net$interactions[[j]]$input
    Adj[j,inputs] = 1
  }
  P_Nodes = which(colSums(Adj) == 0)
  
  approx_attr = getAttractors(net, method = "random", startStates = 10000)
  initial_states_matrix = approx_attr$stateInfo$initialStates
  
  score_A = 0
  score_L = 0
  score_G = 0
  score_P = 0
  count = 0
  for(j in 1:dim(initial_states_matrix)[2]){
    state = c()
    check = TRUE
    for(k in 1:dim(initial_states_matrix)[1]){
      if(is.na(initial_states_matrix[k,j])){
        check = FALSE
      }
    }
    if(check){
      count = count + 1
      for(k in 1:dim(initial_states_matrix)[1]){
        node_num = min(32, N_Num - 32*(k-1))
        state = c(state, rev(Binary_Vector(2^32 + initial_states_matrix[k,j], node_num)))
      }
      path = getPathToAttractor(net, state)
      reached = unlist(path[dim(path)[1],])
      names(reached) = c()
      
      a_state = state
      a_state[network_controllers$Mut_node[i]] = network_controllers$Mut_Value[i]
      a_path = getPathToAttractor(fixGenes(net, network_controllers$Mut_node[i], network_controllers$Mut_Value[i]), a_state)
      a_reached = unlist(a_path[dim(a_path)[1],])
      names(a_reached) = c()
      
      LDOI_state = state
      LDOI_state[c(network_controllers$Mut_node[i], network_controllers$Target_LDOI[i])] = c(network_controllers$Mut_Value[i], network_controllers$Target_LDOI_Value[i])
      LDOI_path = getPathToAttractor(fixGenes(net, c(network_controllers$Mut_node[i], network_controllers$Target_LDOI[i]), c(network_controllers$Mut_Value[i], network_controllers$Target_LDOI_Value[i])), LDOI_state)
      LDOI_reached = unlist(LDOI_path[dim(LDOI_path)[1],])
      names(LDOI_reached) = c()
      
      PLR_state = state
      PLR_state[c(network_controllers$Mut_node[i], as.integer(network_controllers$Target_PLR[i]))] = c(network_controllers$Mut_Value[i], network_controllers$Target_PLR_Value[i])
      PLR_path = getPathToAttractor(fixGenes(net, c(network_controllers$Mut_node[i], as.integer(network_controllers$Target_PLR[i])), c(network_controllers$Mut_Value[i], network_controllers$Target_PLR_Value[i])), PLR_state)
      PLR_reached = unlist(PLR_path[dim(PLR_path)[1],])
      names(PLR_reached) = c()
      
      GS_state = state
      if(Network_Controllers$Target_GS[i] == "None"){
        GS_control = c()
        GS_value = c()
      } else {
        GS_control = as.integer(network_controllers$Target_GS[i])
        GS_value = network_controllers$Target_GS_Value[i]
        GS_state[c(network_controllers$Mut_node[i], GS_control)] = c(network_controllers$Mut_Value[i], GS_value)
        GS_path = getPathToAttractor(fixGenes(net, c(network_controllers$Mut_node[i], GS_control), c(network_controllers$Mut_Value[i], GS_value)), GS_state)
        GS_reached = unlist(GS_path[dim(GS_path)[1],])
        names(GS_reached) = c()
        score_G = score_G + sum(abs(reached[P_Nodes] - GS_reached[P_Nodes]))/length(P_Nodes)
      }
      
      score_A = score_A + sum(abs(reached[P_Nodes] - a_reached[P_Nodes]))/length(P_Nodes)
      score_L = score_L + sum(abs(reached[P_Nodes] - LDOI_reached[P_Nodes]))/length(P_Nodes)
      score_P = score_P + sum(abs(reached[P_Nodes] - PLR_reached[P_Nodes]))/length(P_Nodes)
    }
  }
  Score_Alteration = c(Score_Alteration, score_A/count)
  Score_LDOI = c(Score_LDOI, score_L/count)
  Score_GS = c(Score_GS, score_G/count)
  Score_PLR = c(Score_PLR, score_P/count)
}
Score_Alteration = 1 - Score_Alteration
Score_LDOI = 1 - Score_LDOI
Score_GS = 1 - Score_GS
Score_PLR = 1 - Score_PLR
total_data = data.frame(Score_Alteration, Score_GS, Score_LDOI, Score_PLR)

save(total_data, file = "Fig_3D_data.RData")

Fig3_plot_data = total_data[order(total_data$Score_PLR, decreasing = TRUE),]
rownames(Fig3_plot_data) = c()

scores = c(Fig3_plot_data$Score_PLR, Fig3_plot_data$Score_GS, Fig3_plot_data$Score_LDOI)
methods = c(rep("ARC",dim(Fig3_plot_data)[1]),rep("FVS",dim(Fig3_plot_data)[1]),rep("LDOI",dim(Fig3_plot_data)[1]))
index = rep(seq(1,dim(Fig3_plot_data)[1]),3)
Fig_3B_data = data.frame(scores, methods, index)
ggplot(Fig_3B_data, aes(x = index, y = scores, color = methods)) + geom_line() + theme_bw() + geom_point(shape = 15, size = 3) + labs(x = "Network Index") + labs(y = "Phenotype Landscape Reversion Effectiveness") + scale_x_continuous(breaks = seq(1,33,1)) + theme(legend.background = element_rect(color = "black"))


##### Figure 3C #####

PLR_data = total_data[,c(1,4)]
Distortion_Degree = 1 - PLR_data$Score_Alteration
PLR_data[["Distortion_Degree"]] = Distortion_Degree

save(PLR_data, file = "Fig_3C_data.RData")

ggplot(PLR_data, aes(x = Distortion_Degree, y = Score_PLR)) + geom_point() + theme_bw() + geom_smooth(method = "lm") + labs(x = "Distortion Degree after Alteration") + labs(y = "Phenotype Landscape Reversion Effectiveness")
model = lm(Distortion_Degree ~ Score_PLR, data = PLR_data)
summary(model)