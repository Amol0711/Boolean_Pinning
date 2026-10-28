import math
import numpy as np
import pandas as pd
import random
import time
import matplotlib.pyplot as plt
# 记录开始时间
start_time = time.time()

return_list = []
return_list1 = []
return_list2 = []
return_list3 = []
return_list4 = []
return_list5 = []
return_list6 = []
return_list7 = []
return_list8 = []
return_list9 = []
return_list10 = []
return_list11 = []
return_list12 = []
return_list13 = []
return_list14 = []
return_list15 = []
return_list16 = []
return_list17 = []
return_list18 = []
return_list19 = []
return_list20 = []
return_list21 = []
return_list22 = []
return_list23 = []
return_list24 = []
return_list25 = []
return_list26 = []
return_list27 = []
return_list28 = []
return_list29 = []
A1 = []
A2 = []
A3 = []
A4 = []
A5 = []
A6 = []
A7 = []
A8 = []
A9 = []
A10 = []
A11 = []
A12 = []
A13 = []
A14 = []
A15 = []
A16 = []
A17 = []
A18 = []
A19 = []
A20 = []
A21 = []
A22 = []
A23 = []
A24 = []
A25 = []
A26 = []
A27 = []
A28 = []
A29 = []
pd.set_option('display.max_columns',None)
pd.set_option('display.max_rows',None)
np.random.seed(1)
ACTIONS = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29]
Q_list = ['22222222222222222222222222222']
epsilon = 0.9
alpha = 0.05
gamma = 0.9
num_episodes = 100000
S = [0.5] * 29
S_counter = 0
attractor = ['1','1','1','1','1','1','1','1','0','0','0','0','0','0','1','1','1','0','0','0','0','1','1','1','0','0','0','0','0']
def Entropy(S):
    h = []
    for i in range(29):
        if S[i] == 0 or S[i] == 1:
            h.append(0)
        else:
            h.append(-S[i]*np.log2(S[i])-(1-S[i])*np.log2(1-S[i]))
    return np.sum(h)/len(S)
def action_func(action_num,list):
    if action_num == 1:
        list[0] = 1
    elif action_num == 2:
        list[1] = 1
    elif action_num == 3:
        list[2] = 1
    elif action_num == 4:
        list[3] = 1
    elif action_num == 5:
        list[4] = 1
    elif action_num == 6:
        list[5] = 1
    elif action_num == 7:
        list[6] = 1
    elif action_num == 8:
        list[7] = 1
    elif action_num == 9:
        list[8] = 0
    elif action_num == 10:
        list[9] = 0
    elif action_num == 11:
        list[10] = 0
    elif action_num == 12:
        list[11] = 0
    elif action_num == 13:
        list[12] = 0
    elif action_num == 14:
        list[13] = 0
    elif action_num == 15:
        list[14] = 1
    elif action_num == 16:
        list[15] = 1
    elif action_num == 17:
        list[16] = 1
    elif action_num == 18:
        list[17] = 0
    elif action_num == 19:
        list[18] = 0
    elif action_num == 20:
        list[19] = 0
    elif action_num == 21:
        list[20] = 0
    elif action_num == 22:
        list[21] = 1
    elif action_num == 23:
        list[22] = 1
    elif action_num == 24:
        list[23] = 1
    elif action_num == 25:
        list[24] = 0
    elif action_num == 26:
        list[25] = 0
    elif action_num == 27:
        list[26] = 0
    elif action_num == 28:
        list[27] = 0
    elif action_num == 29:
        list[28] = 0
    return list
#计算H的值
def calculate_H(S1,S2,S3,S4,S5,S6, S7, S8, S9, S10, S11, S12, S13, S14, S15, S16, S17, S18, S19, S20, S21, S22, S23, S24, S25, S26, S27, S28, S29):
    global Hh
    S = list(map(lambda x: x, [S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15, S16, S17, S18, S19, S20, S21, S22, S23, S24, S25, S26, S27, S28, S29]))
    list1 = []
    for i in range(29):
        if S[i] == 1 or S[i] == 0:
            list1.append(i)

    for itn in range(10):
        new_S = np.zeros(29)
        new_S[0] = S[0]
        new_S[1] = S[0]
        new_S[2] = S[1]
        new_S[3] = S[0]
        new_S[4] = S[0]
        new_S[5] = S[3]
        new_S[6] = S[4] * (1 - S[5]) + S[4] * S[5] + (1 - S[4]) * (1 - S[5])
        new_S[7] = S[13] * (S[5] * S[2] * S[4] + S[5] * S[2] * (1 - S[4]) + S[5] * (1 - S[2]) * S[4] + S[5] * (1 - S[2]) * (1 - S[4]) + (1 - S[5]) * S[2] * S[4] + (1 - S[5]) * (1 - S[2]) * (1 - S[4]) + (1 - S[5]) * S[2] * (1 - S[4]) + (1 - S[5]) * (1 - S[2]) * S[4])+(1 - S[13]) * S[5] * (S[2] * S[4] + (1 - S[2]) * S[4] + S[2] * (1 - S[4]))
        new_S[8] = S[8]
        new_S[9] = S[8]
        new_S[10] = S[9]
        new_S[11] = (1 - S[3]) * (1 - S[10])
        new_S[12] = (1 - S[3]) * (1 - S[10])
        new_S[13] = S[10]
        new_S[14] = S[10] * S[15] + S[10] * (1 - S[15]) + (1 - S[10]) * S[15]
        new_S[15] = S[14]
        new_S[16] = S[14]
        new_S[17] = (1 - S[16]) * (S[0] * S[10] + S[0] * (1 - S[10]) + (1 - S[0]) * S[10] + (1 - S[0]) * (1 - S[10]))
        new_S[18] = S[17]
        new_S[19] = (1 - S[0]) * S[18]
        new_S[20] = S[19]
        new_S[21] = S[0]
        new_S[22] = S[1]
        new_S[23] = S[3]
        new_S[24] = S[11]
        new_S[25] = S[13]
        new_S[26] = S[13]
        new_S[27] = S[13]
        new_S[28] = S[10]
        for j in range(len(list1)):
            new_S[list1[j]] = S[list1[j]]
        S = new_S
    if Entropy(S) == 0:
        Hh = 0
    else:
        Hh = 1
    return Hh

def calculate_len(list):
    L = []
    for i in range(29):
        if list[i] == 1 or list[i] == 0:
            L.append(i)
    return len(L)

def build_qtable(n_states, actions):
    table = pd.DataFrame(
        np.zeros((n_states, len(actions))),
        columns=actions,
    )
    return table
def choose_action(state, q_table, action,epsilon):
    state_actions = q_table.iloc[state,:]
    if(np.random.uniform()<epsilon):
        action_name = np.random.choice(action)
    else:
        action_name = state_actions.idxmax()
    return action_name

def reward(H,C):
    if H == 0:
        return 10*(29-C)
    else:
        return -1

def next_S(state,action,number, H):
    S_q = list(Q_list[state])
    A = int(action)-1
    if number > 5 or H == 0:
        Sn = 'terminal'
    else:
        S_q[A] = attractor[A]
        S_l = ''.join(S_q)
        if S_l in Q_list:
            Sn = Q_list.index(S_l)
        else:
            Sn = S_counter+1
            Q_list.append(S_l)
    return Sn

def ql(epsilon):
    q_table = build_qtable(100000, ACTIONS)
    Reward_list = []
    rewards = []
    for episode in range(num_episodes):
        Sq = 0
        S = [0.5] * 29
        if episode % 1000 == 0:
            epsilon = max(0, epsilon * 0.95)
        is_terminated = False
        while not is_terminated:
            A = choose_action(Sq,q_table,ACTIONS,epsilon)
            S = action_func(A,S)
            S1 = S[0]
            S2 = S[1]
            S3 = S[2]
            S4 = S[3]
            S5 = S[4]
            S6 = S[5]
            S7 = S[6]
            S8 = S[7]
            S9 = S[8]
            S10 = S[9]
            S11 = S[10]
            S12 = S[11]
            S13 = S[12]
            S14 = S[13]
            S15 = S[14]
            S16 = S[15]
            S17 = S[16]
            S18 = S[17]
            S19 = S[18]
            S20 = S[19]
            S21 = S[20]
            S22 = S[21]
            S23 = S[22]
            S24 = S[23]
            S25 = S[24]
            S26 = S[25]
            S27 = S[26]
            S28 = S[27]
            S29 = S[28]
            H = calculate_H(S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15, S16, S17, S18, S19, S20, S21, S22, S23, S24, S25, S26, S27, S28, S29)
            C = calculate_len(S)
            R = reward(H,C)
            Reward_list.append(R)
            Sn = next_S(Sq, A, C, H)
            q_predict = q_table.loc[Sq, A]
            A1 = q_table.iloc[645, 0]
            A2 = q_table.iloc[0, 1]
            A3 = q_table.iloc[0, 2]
            A4 = q_table.iloc[0, 3]
            A5 = q_table.iloc[0, 4]
            A6 = q_table.iloc[0, 5]
            A7 = q_table.iloc[0, 6]
            A8 = q_table.iloc[0, 7]
            A9 = q_table.iloc[224, 8]
            A10 = q_table.iloc[0, 9]
            A11 = q_table.iloc[0, 10]
            A12 = q_table.iloc[0, 11]
            A13 = q_table.iloc[0, 12]
            A14 = q_table.iloc[0, 13]
            A15 = q_table.iloc[49, 14]
            A16 = q_table.iloc[49, 15]
            A17 = q_table.iloc[0, 16]
            A18 = q_table.iloc[0, 17]
            A19 = q_table.iloc[0, 18]
            A20 = q_table.iloc[0, 19]
            A21 = q_table.iloc[0, 20]
            A22 = q_table.iloc[0, 21]
            A23 = q_table.iloc[0, 22]
            A24 = q_table.iloc[0, 23]
            A25 = q_table.iloc[0, 24]
            A26 = q_table.iloc[0, 25]
            A27 = q_table.iloc[0, 26]
            A28 = q_table.iloc[0, 27]
            A29 = q_table.iloc[0, 28]
            if Sn != 'terminal':
                q_target = R+gamma*q_table.iloc[Sn,:].max()
            else:
                q_target = R
                is_terminated = True
            q_table.loc[Sq,A] += alpha*(q_target-q_predict)
            Sq = Sn
            return_list1.append(A1)
            return_list2.append(A2)
            return_list3.append(A3)
            return_list4.append(A4)
            return_list5.append(A5)
            return_list6.append(A6)
            return_list7.append(A7)
            return_list8.append(A8)
            return_list9.append(A9)
            return_list10.append(A10)
            return_list11.append(A11)
            return_list12.append(A12)
            return_list13.append(A13)
            return_list14.append(A14)
            return_list15.append(A15)
            return_list16.append(A16)
            return_list17.append(A17)
            return_list18.append(A18)
            return_list19.append(A19)
            return_list20.append(A20)
            return_list21.append(A21)
            return_list22.append(A22)
            return_list23.append(A23)
            return_list24.append(A24)
            return_list25.append(A25)
            return_list26.append(A26)
            return_list27.append(A27)
            return_list28.append(A28)
            return_list29.append(A29)
        if episode % 100 == 0:
            m = np.mean(Reward_list)
            rewards.append(m)
            Reward_list = []
            print(m)
    plt.plot(range(len(rewards)), rewards)
    plt.xlabel('Episode')
    plt.ylabel('Average Reward')
    plt.title('Convergence of Q-learning')
    plt.show()
    return q_table

q_table = ql(epsilon)
print(q_table)

# 记录结束时间
end_time = time.time()
# 计算运行时间
elapsed_time = end_time - start_time
print(f"\n代码运行时间：{elapsed_time}秒")


episodes_list=list(range(len(return_list)))
plt.figure(1)
plt.plot(episodes_list,return_list)
'''plt.legend(ncol=1, loc='upper left')'''
plt.title('')
plt.xlabel('Episodes')
plt.ylabel('Returns')
plt.show()

plt.figure(2)
episodes_list=list(range(len(return_list1)))
plt.plot(episodes_list,return_list1, label='X1')
episodes_list=list(range(len(return_list2)))
plt.plot(episodes_list,return_list2, label='X2')
episodes_list=list(range(len(return_list3)))
plt.plot(episodes_list,return_list3, label='X3')
episodes_list=list(range(len(return_list4)))
plt.plot(episodes_list,return_list4, label='X4')
episodes_list=list(range(len(return_list5)))
plt.plot(episodes_list,return_list5, label='X5')
episodes_list=list(range(len(return_list6)))
plt.plot(episodes_list,return_list6, label='X6')
episodes_list=list(range(len(return_list7)))
plt.plot(episodes_list,return_list7, label='X7')
episodes_list=list(range(len(return_list8)))
plt.plot(episodes_list,return_list8, label='X8')
episodes_list=list(range(len(return_list9)))
plt.plot(episodes_list,return_list9, label='X9')
episodes_list=list(range(len(return_list10)))
plt.plot(episodes_list,return_list10, label='X10')
episodes_list=list(range(len(return_list11)))
plt.plot(episodes_list,return_list11, label='X11')
episodes_list=list(range(len(return_list12)))
plt.plot(episodes_list,return_list12, label='X12')
episodes_list=list(range(len(return_list13)))
plt.plot(episodes_list,return_list13, label='X13')
episodes_list=list(range(len(return_list14)))
plt.plot(episodes_list,return_list14, label='X14')
episodes_list=list(range(len(return_list15)))
plt.plot(episodes_list,return_list15, label='X15')
episodes_list=list(range(len(return_list16)))
plt.plot(episodes_list,return_list16, label='X16')
episodes_list=list(range(len(return_list17)))
plt.plot(episodes_list,return_list17, label='X17')
episodes_list=list(range(len(return_list18)))
plt.plot(episodes_list,return_list18, label='X18')
episodes_list=list(range(len(return_list19)))
plt.plot(episodes_list,return_list19, label='X19')
episodes_list=list(range(len(return_list20)))
plt.plot(episodes_list,return_list20, label='X20')
episodes_list=list(range(len(return_list21)))
plt.plot(episodes_list,return_list21, label='X21')
episodes_list=list(range(len(return_list22)))
plt.plot(episodes_list,return_list22, label='X22')
episodes_list=list(range(len(return_list23)))
plt.plot(episodes_list,return_list23, label='X23')
episodes_list=list(range(len(return_list24)))
plt.plot(episodes_list,return_list24, label='X24')
episodes_list=list(range(len(return_list25)))
plt.plot(episodes_list,return_list25, label='X25')
episodes_list=list(range(len(return_list26)))
plt.plot(episodes_list,return_list26, label='X26')
episodes_list=list(range(len(return_list27)))
plt.plot(episodes_list,return_list27, label='X27')
episodes_list=list(range(len(return_list28)))
plt.plot(episodes_list,return_list28, label='X28')
episodes_list=list(range(len(return_list29)))
plt.plot(episodes_list,return_list29, label='X29')
episodes_list=list(range(len(return_list30)))
plt.legend(ncol=3, loc='upper left', fontsize=8)
plt.title('')
plt.xlabel('Episodes')
plt.ylabel('Returns')
plt.show()



