discard
adopath + "C:\Users\a-e235920k\Documents\Stage Victor Rechard\ado"

rcall clear

import delimited "C:\Users\a-e235920k\Documents\Stage Victor Rechard\code\data_sim_dif_suite\sim_N(500)_M(4)_J(5)_grEFF(-.2)_var(1)_h",clear
keep if replication==10

mat test=J(2,2,1)
matlist test



rcall data=st.data()
rcall test=st.matrix(test)
rcall : source("C:/Users/a-e235920k/Documents/Stage Victor Rechard/code/toStataWLE.R")
rcall st.load(data)

