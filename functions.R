
Block_Evaluation = function(X, LMT){
  
  X_Seg_Block = X
  
  K = X_Seg_Block %>% select(IDENTIFIER, AH, RLS, Iab, Ian, FL, Pil2, EM, Pil, LL)
  
  H = melt(K, id.vars = 'IDENTIFIER')
  
  H$variable = as.character(H$variable)
  
  H2 = H %>% left_join(LMT, by = 'variable')
  
  H3 = H2 %>% mutate(SINGLE_NOK_FLAG = (value < LCL) | (value > UCL))
  
  S = dcast(H3, IDENTIFIER~variable, value.var = 'SINGLE_NOK_FLAG', fun.aggregate = function(x) x[1])
  
  S2 = S %>% mutate(NOK = AH | RLS | Iab | Ian | FL | Pil2 | EM | Pil | LL) %>%
    mutate(NOK_ASS = AH | RLS | Iab | Ian) %>% 
    mutate(NOK_EMI =  FL | Pil2 | EM | Pil | LL)
  
  Act_FPY = 1 - sum(S2$NOK)/length(S2$NOK)
  Act_FPY_ASS = 1 - sum(S2$NOK_ASS)/length(S2$NOK_ASS)
  Act_FPY_EMI = 1 - sum(S2$NOK_EMI)/length(S2$NOK_EMI)
  
  
  S_EST = H3 %>% group_by(variable) %>% dplyr::summarise(
    MEAN = mean(value, na.rm = TRUE),
    SD = sd(value, na.rm = TRUE),
    LCL = unique(LCL),
    UCL = unique(UCL)
  ) %>% mutate(
    LEFT_TAIL = pnorm(LCL, MEAN, SD),
    RIGHT_TAIL = 1 - pnorm(UCL, MEAN, SD)
  ) %>% mutate(
    TAILS = LEFT_TAIL + RIGHT_TAIL
  ) %>% mutate(
    EST_FPY = 1 - TAILS
  )
  
  Est_FPY = prod(S_EST$EST_FPY)
  Est_FPY_EMI = prod(S_EST$EST_FPY[S_EST$variable %in% c('FL', 'Pil2', 'EM', 'Pil', 'LL')])
  Est_FPY_ASS = prod(S_EST$EST_FPY[! S_EST$variable %in% c('FL', 'Pil2', 'EM', 'Pil', 'LL')])
  
  AH_MEAN = mean(K$AH, na.rm = TRUE)
  Iab_MEAN = mean(K$Iab, na.rm = TRUE)
  AH_MEDIAN = median(K$AH, na.rm = TRUE)
  Iab_MEDIAN = median(K$Iab, na.rm = TRUE)
  
  R = data.frame(Act_FPY, Act_FPY_ASS, Act_FPY_EMI, Est_FPY, Est_FPY_ASS, Est_FPY_EMI, AH_MEAN, AH_MEDIAN, Iab_MEAN, Iab_MEDIAN)
  
  return(R)
}


interval_segmentator = function(x, MAX_INTERVAL_TIME = 60*60, RETURN = 'SEG', TYPE = 'time',VERBOSE = FALSE){
  N = length(x)
  if(TYPE == 'time'){
    x = lubridate::ymd_hms(x)
    DIFF = as.numeric(as.duration(x[2:N] - x[1:(N-1)]))
  } else {
    DIFF = diff(x)
  }
  FIRSTS = which(DIFF > MAX_INTERVAL_TIME) + 1
  if(length(FIRSTS) == 0) {
    L = list(1:N)
    L2 = rep(1,N)
  } else {
    m = length(FIRSTS)
    L = list()
    for(i in 0:m){
      if(i == 0){
        INDEX = 1:(FIRSTS[1]-1)
      } else if( i == m) {
        INDEX = FIRSTS[m]:N
      } else {
        INDEX = FIRSTS[i]:(FIRSTS[i+1]-1)
      }
      L = rlist::list.append(L, INDEX)
    }
    LEN = length(L)
    L2 = unlist(lapply(1:LEN, function(i) {rep(i, length(L[[i]]))}))
  }
  if(VERBOSE){
    print('Return: SEG: segmentation vector, others: INDEX in a list')
    print(paste('You have chosen Returning', RETURN))
  }
  if(RETURN == 'SEG') R = L2 else R = L
  return(R)
}
