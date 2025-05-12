##### library #####

library(progress) # Version 1.2.2
library(BoolNet) # Version 2.1.5
library(ggplot2) # Version 3.4.2
library(ggsignif) # Version 0.6.4
library(igraph) # Version 1.2.6
library(visNetwork) # Version 2.1.2

##### 1. Basic Functions #####

Binary_Vector = function(Number, Length){
  # Input : decimal number, length of the binary vector
  # Output : binary vector
  Binary = c()
  if((Number+1) > 2^Length){
    Number = Number %% (2^Length)
  }
  for(i in 1:Length){
    Quotient = Number %/% (2^(Length-i))
    Binary = c(Binary, Quotient)
    if(Quotient == 1){
      Number = Number - (2^(Length-i))
    }
  }
  return(Binary)
}

Decimal_Vector = function(Vector){
  # Input : binary vector
  # Output : decimal number
  Decimal = 0
  for(i in 1:length(Vector)){
    Value = (2 ^ (length(Vector)-i)) * Vector[i]
    Decimal = Decimal + Value
  }
  names(Decimal) = c()
  return(Decimal)
}

MM = function(Matrix, Length){
  # Input : matrix, number of the multiplication
  # Output : computed matrix
  M = diag(1, dim(Matrix)[1])
  if(Length == 0){
    return(M)
  } else {
    for(i in 1:Length){
      M =  M %*% Matrix
    }
    return(M)
  }
}

##### 2. STP-based functions #####

E_mat = function(Mutation_Indices, Mutation_Values, N_Num){
  # Input : Fixed nodes indices & values, number of network nodes
  # Output : effect matrix of given nodes fixation
  I_mat = diag(1, (2^N_Num))
  for(i in 1:(2^N_Num)){
    state_vector = Binary_Vector((2^N_Num - i), N_Num)
    mutated_state_vector = state_vector
    mutated_state_vector[Mutation_Indices] = Mutation_Values
    new_column_vector = rep(0, (2^N_Num))
    new_column_vector[2^N_Num - Decimal_Vector(mutated_state_vector)] = 1
    I_mat[,i] = new_column_vector
  }
  return(I_mat)
}

Fast_E = function(Mutation_Indices, Mutation_Values, N_Num){
  # Input : Fixed nodes indices & values, number of network nodes
  # Output : effect matrix of given nodes fixation / can be only used in fast calculation process
  I_mat = c()
  for(i in 1:(2^N_Num)){
    state_vector = Binary_Vector((2^N_Num - i), N_Num)
    state_vector[Mutation_Indices] = Mutation_Values
    I_mat = c(I_mat, (2^N_Num - Decimal_Vector(state_vector)))
  }
  return(I_mat)
}

L_mat = function(Network, Mutation_Indices = c(), Mutation_Values = c()){
  # Input : Fixed nodes indices & values, network
  # Output : state transition matrix with the fixation of given nodes
  N_Num = length(Network$interactions)
  pb = progress_bar$new(total = 2^N_Num)
  for(i in 1:(2^N_Num)){
    pb$tick()
    state = Binary_Vector(((2^N_Num) - i), N_Num)
    next_state_dec = (2^N_Num) - Decimal_Vector(stateTransition(Network, state))
    vector = rep(0, 2^N_Num)
    vector[next_state_dec] = 1
    if(i == 1){
      l_mat = matrix(vector, ncol = 1)
    } else {
      l_mat = cbind(l_mat, vector)
    }
  }
  if(length(Mutation_Indices) > 0){
    E = E_mat(Mutation_Indices, Mutation_Values, N_Num)
    return(E %*% l_mat %*% E)
  }
  return(l_mat)
}

Fast_L = function(Network){
  # Input : Fixed nodes indices & values, network
  # Output : state transition matrix in the nominal network / can be only used in fast calculation process
  N_Num = length(Network$interactions)
  pb = progress_bar$new(total = 2^N_Num)
  l_mat = c()
  for(i in 1:(2^N_Num)){
    pb$tick()
    state = Binary_Vector(((2^N_Num) - i), N_Num)
    next_state_dec = (2^N_Num) - Decimal_Vector(stateTransition(Network, state))
    l_mat = c(l_mat, next_state_dec)
  }
  return(l_mat)
}

A_mat = function(Network, Mutation_Indices = c(), Mutation_Values = c()){
  # Input : Fixed nodes indices & values, network
  # Output : attractor matrix with the fixation of given nodes
  L = L_mat(Network, Mutation_Indices, Mutation_Values)
  eig = eigen(L)
  e_values = eig$values
  e_vectors = eig$vectors
  A = sign(abs(e_vectors[,which(((round(Re(e_values),5) == 1) & ((round(Im(e_values),5)) == 0)))]))
  return(matrix(A, nrow = dim(L)[1]))
}

B_mat = function(Network, Mutation_Indices = c(), Mutation_Values = c()){
  # Input : Fixed nodes indices & values, network
  # Output : basin matrix with the fixation of given nodes
  L = t(L_mat(Network, Mutation_Indices, Mutation_Values))
  eig = eigen(L)
  e_values = eig$values
  e_vectors = eig$vectors
  B = sign(abs(e_vectors[,which(((round(Re(e_values),5) == 1) & ((round(Im(e_values),5)) == 0)))]))
  return(matrix(B, nrow = dim(L)[1]))
}

P_mat = function(N_Num, Phenotype_Nodes){
  # Input : number of network nodes, output nodes indices
  # Output : attractor to phenotype converter matrix
  P = matrix(0, nrow = 2^length(Phenotype_Nodes), ncol = 2^N_Num)
  for(i in 1:(2^N_Num)){
    original_vector = Binary_Vector((2^N_Num - i), N_Num)
    pheno_vector = original_vector[Phenotype_Nodes]
    pheno_index = 2^length(Phenotype_Nodes) - Decimal_Vector(pheno_vector)
    P[pheno_index,i] = 1
  }
  return(P)
}

Fast_P = function(N_Num, P_Nodes){
  # Input : number of network nodes, output nodes indices
  # Output : attractor to phenotype converter matrix / can be only used in fast calculation process
  P = c()
  for(i in 1:(2^N_Num)){
    P_vector = Binary_Vector((2^N_Num) - i, N_Num)[P_Nodes]
    P = c(P, (1+Decimal_Vector(P_vector)))
  }
  return(P)
}

##### 3. Approximated STP-based functions #####

L_E_mat = function(Mutation_Indices, Mutation_Values, N_Num){
  # Input : Fixed nodes indices & values, number of network nodes
  # Output : effect matrix of given nodes fixation in large network
  I_mat = diag(1, (1 + N_Num))
  for(i in 1:length(Mutation_Indices)){
    if(Mutation_Values[i] == 0){
      I_mat[1+Mutation_Indices[i], 1+Mutation_Indices[i]] = 0
    } else {
      I_mat[1+Mutation_Indices[i], 1+Mutation_Indices[i]] = 0
      I_mat[1+Mutation_Indices[i], 1] = 1
    }
  }
  return(I_mat)
}

Large_L_F = function(Network){
  # Input : network
  # Output : approximated nominal forward state transition matrix in large network
  N_Num = length(Network$genes)
  Adj = matrix(0, N_Num, N_Num)
  
  ## Adj Matrix ##
  for(i in 1:N_Num){
    inputs = Network$interactions[[i]]$input
    Adj[i,inputs] = 1
  }
  Inputs = c()
  for(i in 1:N_Num){
    input_vector = rep(0,N_Num)
    input_vector[i] = 1
    if(identical(Adj[i,], input_vector)){
      Inputs = c(Inputs, i)
    }
  }
  Outputs = which(colSums(Adj) == 0)
  
  ## L Matrix ##
  L = matrix(0, 1 + N_Num, 1 + N_Num)
  L[1,1] = 1
  for(i in 1:N_Num){
    inputs = Network$interactions[[i]]$input
    if(identical(inputs,0)){
      L[1+i, 1] = Network$interactions[[i]]$func[1]
    } else {
      one_index = which(Network$interactions[[i]]$func == 1)
      if(length(one_index) > 0){
        for(j in 1:length(one_index)){
          input_vector = Binary_Vector(one_index[j]-1, length(inputs))
          for(k in 1:length(inputs)){
            L[1+i, 1+inputs[k]] = L[1+i, 1+inputs[k]] + 2 * sign(input_vector[k] - 0.5)
            L[1+i, 1] = L[1+i, 1] - sign(input_vector[k] - 0.5)
          }
        }
        L[1+i,] = L[1+i,] / (2^length(inputs))
      }
    }
  }
  
  ## Bias Correction ##
  Correct_L = matrix(0, 1+N_Num, 1+N_Num)
  Correct_L[1,1] = 1
  for(i in 1:N_Num){
    if(identical(L[1+i,c(2:(1+N_Num))],rep(0,N_Num)) == FALSE){
      minn = 0
      maxx = 0
      for(j in 1:N_Num){
        minn = minn + min(0, L[1+i,1+j])
        maxx = maxx + max(0, L[1+i,1+j])
      }
      Correct_L[1+i,2:(1 + N_Num)] = L[1+i,2:(1 + N_Num)] * (1 / (maxx - minn))
      Correct_L[1+i,1] = -minn / (maxx - minn)
    } else {
      Correct_L[1+i,1] = L[1+i,1]
    }
  }
  return(Correct_L)
}

Large_L_B = function(Network){
  # Input : network
  # Output : approximated nominal backward state transition matrix in large network
  N_Num = length(Network$genes)
  Adj = matrix(0, N_Num, N_Num)
  
  ## Adj Matrix ##
  for(i in 1:N_Num){
    inputs = Network$interactions[[i]]$input
    Adj[i,inputs] = 1
  }
  Inputs = c()
  for(i in 1:N_Num){
    input_vector = rep(0,N_Num)
    input_vector[i] = 1
    if(identical(Adj[i,], input_vector)){
      Inputs = c(Inputs, i)
    }
  }
  Outputs = which(colSums(Adj) == 0)
  
  ## B Matrix ##
  B = matrix(0, 1 + N_Num, 1 + N_Num)
  B[1,1] = 1
  for(i in 1:N_Num){
    #print(i)
    if((i %in% Outputs) == FALSE){
      targets = which(Adj[,i] > 0)
      one_prob = matrix(0, nrow = 2, ncol = length(targets))
      for(j in 1:length(targets)){
        inputs = Network$interactions[[targets[j]]]$input
        for(k in 1:(2^length(inputs))){
          i_value = Binary_Vector((k-1), length(inputs))[which(inputs == i)]
          one_prob[1+Network$interactions[[targets[j]]]$func[k],j] = one_prob[1+Network$interactions[[targets[j]]]$func[k],j] + i_value
        }
        if(length(unique(Network$interactions[[targets[j]]]$func)) == 2){
          one_prob[1,j] = one_prob[1,j] / length(which(Network$interactions[[targets[j]]]$func == 0))
          one_prob[2,j] = one_prob[2,j] / length(which(Network$interactions[[targets[j]]]$func == 1))
        } else {
          if(unique(Network$interactions[[targets[j]]]$func)[1] == 0){
            one_prob[1,j] = 1
            one_prob[2,j] = 0
          } else {
            one_prob[1,j] = 0
            one_prob[2,j] = 1
          }
        }
      }
      for(j in 1:(2^length(targets))){
        target_vector = Binary_Vector((j-1), length(targets))
        prob = 1
        for(k in 1:length(targets)){
          prob = prob * one_prob[1+target_vector[k],k]
        }
        for(k in 1:length(targets)){
          B[1+i, 1+targets[k]] = B[1+i, 1+targets[k]] + 2 * sign(target_vector[k] - 0.5) * prob
          B[1+i, 1] = B[1+i, 1] - sign(target_vector[k] - 0.5) * prob
        }
      }
      if(identical(B[1+i,], rep(0, 1+N_Num))){
        B[1+i, 1] = 0.5
      }
    } else {
      B[1+i,1+i] = 1
    }
  }
  
  ## Bias Correction ##
  Correct_B = matrix(0, 1+N_Num, 1+N_Num)
  Correct_B[1,1] = 1
  for(i in 1:N_Num){
    if((identical(B[1+i,2:(1+N_Num)], rep(0, N_Num))) == FALSE){
      minn = 0
      maxx = 0
      for(j in 1:N_Num){
        minn = minn + min(0, B[1+i,1+j])
        maxx = maxx + max(0, B[1+i,1+j])
      }
      Correct_B[1+i,2:(1 + N_Num)] = B[1+i,2:(1 + N_Num)] * (1 / (maxx - minn))
      Correct_B[1+i,1] = -minn / (maxx - minn)
    } else {
      Correct_B[1+i,] = B[1+i,]
    }
  }
  return(Correct_B)
}

##### 4. Control Target Identification in Small Network #####

Score = function(Network, Mutation_Indices = c(), Mutation_Values = c(), P_Nodes){
  # Input : Fixed nodes indices & values, network, output nodes indices
  # Output : distortion degree of phenotype landscape for given nodes fixations
  L = L_mat(Network)
  E = E_mat(Mutation_Indices, Mutation_Values, length(Network$genes))
  P = P_mat(length(Network$genes), P_Nodes)
  
  eig = eigen(t(L))
  e_values = eig$values
  e_vectors = eig$vectors
  B = sign(abs(e_vectors[,which(((round(Re(e_values),5) == 1) & ((round(Im(e_values),5)) == 0)))]))
  B = matrix(B, nrow = dim(L)[1])
  k = 0
  check_B = B
  while(TRUE){
    prev_B = check_B
    check_B = sign(L %*% check_B)
    k = k + 1
    if(identical(prev_B, check_B)){
      k = k - 1
      break
    }
  }
  
  nom_mat = P %*% MM(L, k) %*% B
  mat = P %*% MM(E %*% L, k) %*% E %*% B
  print(mat)
  return(sum(abs(nom_mat - mat)) / (2^(length(Network$genes)+1)))
}

BF_Search = function(Network, Mutation_Indices, Mutation_Values, P_Nodes){
  # Input : Altered nodes indices & values, network, output nodes indices
  # Output : distortion degrees of controlled phenotype landscapes for every single node control
  L = L_mat(Network)
  E = E_mat(Mutation_Indices, Mutation_Values, length(Network$genes))
  P = P_mat(length(Network$genes), P_Nodes)
  
  eig = eigen(t(L))
  e_values = eig$values
  e_vectors = eig$vectors
  B = sign(abs(e_vectors[,which(((round(Re(e_values),5) == 1) & ((round(Im(e_values),5)) == 0)))]))
  B = matrix(B, nrow = dim(L)[1])
  k = 0
  check_B = B
  while(TRUE){
    prev_B = check_B
    check_B = sign(L %*% check_B)
    k = k + 1
    if(identical(prev_B, check_B)){
      k = k - 1
      break
    }
  }
  
  nom_mat = P %*% MM(L, k) %*% B
  mat = P %*% MM(E %*% L, k) %*% E %*% B
  score = sum(abs(nom_mat - mat)) / (2^(length(Network$genes)+1))
  
  Node = c("None")
  Value = c("M")
  Score = c(score)
  
  for(i in 1:length(Network$genes)){
    for(j in 1:2){
      print(c(i,j))
      if((i %in% Mutation_Indices) == FALSE){
        C = E_mat(c(Mutation_Indices, i), c(Mutation_Values, (j-1)), length(Network$genes))
        mat = P %*% MM(C %*% L, k) %*% C %*% B
        score = sum(abs(nom_mat - mat)) / (2^(length(Network$genes)+1))
        Node = c(Node, i)
        Value = c(Value, (j-1))
        Score = c(Score, score)
      }
    }
  }
  return(data.frame(Node, Value, Score))
}

BF_Simul = function(Network, Mutation_Indices, Mutation_Values, P_Nodes){
  # Input : Altered nodes indices & values, network, output nodes indices
  # Output : distortion degrees of controlled phenotype landscapes for every single node control calculated through Brute-Force simulation
  N_Num = length(Network$genes)
  Attr = getAttractors(Network)
  length = max(Attr$stateInfo$stepsToAttractor)
  
  diff = 0
  for(i in 1:2^N_Num){
    vector = Binary_Vector((i-1), N_Num)
    m_vector = vector
    m_vector[Mutation_Indices] = Mutation_Values
    for(j in 1:length){
      vector = stateTransition(Network, vector)
      m_vector = stateTransition(fixGenes(Network, Mutation_Indices, Mutation_Values), m_vector)
    }
    if(identical(vector[P_Nodes], m_vector[P_Nodes]) == FALSE){
      diff = diff + 1/(2^N_Num)
    }
  }
  
  Node = c("None")
  Value = c("M")
  Score = c(diff)
  for(i in 1:N_Num){
    if(i != Mutation_Indices){
      for(j in 1:2){
        diff = 0
        for(k in 1:2^N_Num){
          vector = Binary_Vector((k-1), N_Num)
          c_vector = vector
          c_vector[c(i, Mutation_Indices)] = c((j-1), Mutation_Values)
          for(l in 1:length){
            vector = stateTransition(Network, vector)
            c_vector = stateTransition(fixGenes(Network, c(i,Mutation_Indices), c((j-1),Mutation_Values)), c_vector)
          }
          if(identical(vector[P_Nodes], c_vector[P_Nodes]) == FALSE){
            diff = diff + 1/(2^N_Num)
          }
        }
        Node = c(Node, i)
        Value = c(Value, (j-1))
        Score = c(Score, diff)
      }
    }
  }
  return(data.frame(Node, Value, Score))
}

##### 5. Control Target Identification in Large Network #####

Large_Attr = function(Network, startStates = 100000){
  # Input : Network, number of initial states
  # Output : major attractor states & basin size ratio
  attr = getAttractors(Network, method = "random", startStates = startStates)
  N_Num = length(Network$genes)
  matrix_elements = c()
  attr_size = c()
  for(i in 1:length(attr$attractors)){
    attr_vector = rep(0, N_Num)
    for(j in 1:dim(attr$attractors[[i]]$involvedStates)[2]){
      attr_partial_vector = c()
      for(k in 1:dim(attr$attractors[[i]]$involvedStates)[1]){
        node_num = min(32, N_Num-32*(k-1))
        attr_partial_vector = c(attr_partial_vector, rev(Binary_Vector(2^32+attr$attractors[[i]]$involvedStates[k,j],node_num)))
      }
      attr_vector = attr_vector + attr_partial_vector
    }
    attr_vector = attr_vector / dim(attr$attractors[[i]]$involvedStates)[2]
    matrix_elements = c(matrix_elements, 1)
    for(j in 1:length(attr_vector)){
      matrix_elements = c(matrix_elements, attr_vector[j])
    }
    attr_size = c(attr_size, attr$attractors[[i]]$basinSize)
  }
  
  result = list()
  result[["Attr_Mat"]] = matrix(matrix_elements, ncol = length(attr$attractors))
  result[["Attr_Size"]] = attr_size
  return(result)
}

Large_Score = function(Network, Attr_mat, Attr_size, Mutation_Indices = c(), Mutation_Values = c(), P_Nodes = c(), Length = 15){
  # Input : Network, major attractor states & basin size ratio, fixed nodes indices & values, output nodes indices, state transition length
  # Output : distortion degree of the phenotype landscape for given nodes fixations
  N_Num = length(Network$genes)
  #N = L_STP_Taylor(Network)
  L_F_M = Large_L_F(Network)
  L_B_M = Large_L_B(Network)
  N = list()
  N[["Forward"]] = L_F_M
  N[["Backward"]] = L_B_M
  E = L_E_mat(Mutation_Indices, Mutation_Values, N_Num)
  #N_P = MM(N$Forward, Length) %*% MM(N$Backward, Length) %*% Attr_mat
  N_P = Attr_mat
  M_P = MM(E %*% N$Forward, Length) %*% E %*% MM(N$Backward, Length) %*% Attr_mat
  for(i in 1:length(Attr_size)){
    N_P[,i] = N_P[,i] * Attr_size[i]
    M_P[,i] = M_P[,i] * Attr_size[i]
  }
  
  diff = abs(N_P[P_Nodes+1,] - M_P[P_Nodes+1,])
  return(sum(diff) / (sum(Attr_size) * length(P_Nodes)))
}

Large_BF_Search = function(Network, Attr_mat, Attr_size, Mutation_Indices, Mutation_Values, P_Nodes = c(), Length = 15){
  # Input : Network, major attractor states & basin size ratio, fixed nodes indices & values, output nodes indices, state transition length
  # Output : distortion degrees of the phenotype landscape for every single node control
  N_Num = length(Network$genes)
  #N = L_STP_Taylor(Network)
  L_F_M = Large_L_F(Network)
  L_B_M = Large_L_B(Network)
  N = list()
  N[["Forward"]] = L_F_M
  N[["Backward"]] = L_B_M
  E = L_E_mat(Mutation_Indices, Mutation_Values, N_Num)
  #N_P = MM(N$Forward, Length) %*% MM(N$Backward, Length) %*% Attr_mat
  N_P = Attr_mat
  M_P = MM(E %*% N$Forward, Length) %*% E %*% MM(N$Backward, Length) %*% Attr_mat
  for(i in 1:length(Attr_size)){
    N_P[,i] = N_P[,i] * Attr_size[i]
    M_P[,i] = M_P[,i] * Attr_size[i]
  }
  
  diff = abs(N_P[P_Nodes+1,] - M_P[P_Nodes+1,])
  s = (sum(diff) / (sum(Attr_size) * length(P_Nodes)))
  
  Node = c("M")
  Value = c(-1)
  Score = c(s)
  for(i in 1:N_Num){
    if((i %in% c(Mutation_Indices,P_Nodes)) == FALSE){
      for(j in 1:2){
        #print(c(i,j))
        C = L_E_mat(c(i,Mutation_Indices), c((j-1),Mutation_Values), N_Num)
        C_P = MM(C %*% N$Forward, Length) %*% C %*% MM(N$Backward, Length) %*% Attr_mat
        for(k in 1:length(Attr_size)){
          C_P[,k] = C_P[,k] * Attr_size[k]
        }
        diff = abs(N_P[P_Nodes+1,] - C_P[P_Nodes+1,])
        s = (sum(diff) / (sum(Attr_size) * length(P_Nodes)))
        Node = c(Node, i)
        Value = c(Value, (j-1))
        Score = c(Score, s)
      }
    }
  }
  result = data.frame(Node,Value,Score)
  return(result)
}

Large_BF_Score = function(Network, Mutation_Indices, Mutation_Values, Control_Indices, Control_Values, P_Nodes, startStates = 100000){
  # Input : Network, fixed nodes indices & values, Controlled node indices & values, output nodes indices, number of random initial states
  # Output : distortion degree of the phenotype landscape for given nodes fixations computed by Brute-Force simulation
  N_Num = length(Network$genes)
  startStates = min(startStates, 2^N_Num)
  Attr = getAttractors(Network, method = "random", startStates = startStates)
  length = max(Attr$stateInfo$stepsToAttractor)
  
  diff = 0
  count = 0
  for(i in 1:dim(Attr$stateInfo$initialStates)[2]){
    check = TRUE
    for(j in 1:dim(Attr$stateInfo$initialStates)[1]){
      if(is.na(Attr$stateInfo$initialStates[j,i])){
        check = FALSE
      }
    }
    if(check){
      count = count + 1
      vector = c()
      for(j in 1:dim(Attr$stateInfo$initialStates)[1]){
        node_num = min(32, N_Num - 32*(j-1))
        if((is.na(Attr$stateInfo$initialStates[j,i])) == FALSE){
          
        }
        vector = c(vector, rev(Binary_Vector(2^32+Attr$stateInfo$initialStates[j,i], node_num)))
      }
      m_vector = vector
      m_vector[Mutation_Indices] = Mutation_Values
      m_vector[Control_Indices] = Control_Values
      for(j in 1:length){
        vector = stateTransition(Network, vector)
        m_vector = stateTransition(fixGenes(Network, c(Mutation_Indices, Control_Indices), c(Mutation_Values, Control_Values)), m_vector)
      }
      diff = diff + sum(abs(vector[P_Nodes] - m_vector[P_Nodes]))/length(P_Nodes)
    }
  }
  diff = diff / count
  return(diff)
}

Large_BF_Simul = function(Network, Mutation_Indices, Mutation_Values, P_Nodes, startStates = 100000){
  # Input : Network, fixed nodes indices & values, output nodes indices, number of random initial states
  # Output : distortion degrees of the phenotype landscape for every single node control computed by Brute-Force simulation
  N_Num = length(Network$genes)
  startStates = min(startStates, 2^N_Num)
  Attr = getAttractors(Network, method = "random", startStates = startStates)
  length = max(Attr$stateInfo$stepsToAttractor)
  
  diff = 0
  count = 0
  for(i in 1:dim(Attr$stateInfo$initialStates)[2]){
    check = TRUE
    for(j in 1:dim(Attr$stateInfo$initialStates)[1]){
      if(is.na(Attr$stateInfo$initialStates[j,i])){
        check = FALSE
      }
    }
    if(check){
      count = count + 1
      vector = c()
      for(j in 1:dim(Attr$stateInfo$initialStates)[1]){
        node_num = min(32, N_Num - 32*(j-1))
        if((is.na(Attr$stateInfo$initialStates[j,i])) == FALSE){
          
        }
        vector = c(vector, rev(Binary_Vector(2^32+Attr$stateInfo$initialStates[j,i], node_num)))
      }
      m_vector = vector
      m_vector[Mutation_Indices] = Mutation_Values
      for(j in 1:length){
        vector = stateTransition(Network, vector)
        m_vector = stateTransition(fixGenes(Network, Mutation_Indices, Mutation_Values), m_vector)
      }
      if(identical(vector[P_Nodes], m_vector[P_Nodes]) == FALSE){
        diff = diff + 1
      }
    }
  }
  diff = diff / count
  
  Node = c("None")
  Value = c("M")
  Score = c(diff)
  for(i in 1:N_Num){
    if(((i %in% Mutation_Indices) == FALSE) & ((i %in% P_Nodes) == FALSE)){
      for(j in 1:2){
        # print(c(i,j))
        diff = 0
        count = 0
        for(k in 1:dim(Attr$stateInfo$initialStates)[2]){
          check = TRUE
          for(l in 1:dim(Attr$stateInfo$initialStates)[1]){
            if(is.na(Attr$stateInfo$initialStates[l,k])){
              check = FALSE
            }
          }
          if(check){
            count = count + 1
            vector = c()
            for(l in 1:dim(Attr$stateInfo$initialStates)[1]){
              node_num = min(32, N_Num - 32*(l-1))
              vector = c(vector, rev(Binary_Vector(2^32+Attr$stateInfo$initialStates[l,k], node_num)))
            }
            c_vector = vector
            c_vector[c(i, Mutation_Indices)] = c((j-1), Mutation_Values)
            for(l in 1:length){
              vector = stateTransition(Network, vector)
              c_vector = stateTransition(fixGenes(Network, c(i,Mutation_Indices), c((j-1),Mutation_Values)), c_vector)
            }
            if(identical(vector[P_Nodes], c_vector[P_Nodes]) == FALSE){
              diff = diff + 1
            }
          }
        }
        diff = diff / count
        Node = c(Node, i)
        Value = c(Value, (j-1))
        Score = c(Score, diff)
      }
    }
  }
  return(data.frame(Node, Value, Score))
}