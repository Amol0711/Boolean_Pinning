import math
import numpy as np
import pandas as pd
import random
import time
import matplotlib.pyplot as plt
import csv
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
return_list30 = []
return_list31 = []
return_list32 = []
return_list33 = []
return_list34 = []
return_list35 = []
return_list36 = []
return_list37 = []
return_list38 = []
return_list39 = []
return_list40 = []
return_list41 = []
return_list42 = []
return_list43 = []
return_list44 = []
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
A30 = []
A31 = []
A32 = []
A33 = []
A34 = []
A35 = []
A36 = []
A37 = []
A38 = []
A39 = []
A40 = []
A41 = []
A42 = []
A43 = []
A44 = []
pd.set_option('display.max_columns',None)
pd.set_option('display.max_rows',None)
np.random.seed(2)
ACTIONS = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44]
Q_list = ['22222222222222222222222222222222222222222222']
epsilon = 0.9
alpha = 0.05
gamma = 0.9
num_episodes = 100000
S = [0.5] * 44
S_counter = 0
attractor = ['1','1','1','1','1','1','0','0','0','0','0','0','0','0','0','0','0','1','0','0','0','0','1','1','1','0','1','0','1','0','0','1','0','0','1','0','0','0','1','1','0','0','1','1']
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
        list[6] = 0
    elif action_num == 8:
        list[7] = 0
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
        list[14] = 0
    elif action_num == 16:
        list[15] = 0
    elif action_num == 17:
        list[16] = 0
    elif action_num == 18:
        list[17] = 1
    elif action_num == 19:
        list[18] = 0
    elif action_num == 20:
        list[19] = 0
    elif action_num == 21:
        list[20] = 0
    elif action_num == 22:
        list[21] = 0
    elif action_num == 23:
        list[22] = 1
    elif action_num == 24:
        list[23] = 1
    elif action_num == 25:
        list[24] = 1
    elif action_num == 26:
        list[25] = 0
    elif action_num == 27:
        list[26] = 1
    elif action_num == 28:
        list[27] = 0
    elif action_num == 29:
        list[28] = 1
    elif action_num == 30:
        list[29] = 0
    elif action_num == 31:
        list[30] = 0
    elif action_num == 32:
        list[31] = 1
    elif action_num == 33:
        list[32] = 0
    elif action_num == 34:
        list[33] = 0
    elif action_num == 35:
        list[34] = 1
    elif action_num == 36:
        list[35] = 0
    elif action_num == 37:
        list[36] = 0
    elif action_num == 38:
        list[37] = 0
    elif action_num == 39:
        list[38] = 1
    elif action_num == 40:
        list[39] = 1
    elif action_num == 41:
        list[40] = 0
    elif action_num == 42:
        list[41] = 0
    elif action_num == 43:
        list[42] = 1
    elif action_num == 44:
        list[43] = 1
    return list
#计算H的值
def Entropy(S):
    h = []
    for i in range(44):
        if S[i] == 0 or S[i] == 1:
            h.append(0)
        else:
            h.append(-S[i]*np.log2(S[i])-(1-S[i])*np.log2(1-S[i]))
    return np.sum(h)/len(S)

def calculate_len(list):
    L = []
    for i in range(44):
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
    if(np.random.uniform()<epsilon) or ((state_actions==0).all()):
        action_name = np.random.choice(action)
    else:
        action_name = state_actions.idxmax()
    return action_name

def reward(H,C):
    if H == 0:
        return 10*(44-C)
    else:
        return -1
def next_S(state,action,number, H):
    S_q = list(Q_list[state])
    A = int(action)-1
    if H == 0:
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

def calculate_H(S1,S2,S3,S4,S5,S6, S7, S8, S9, S10, S11, S12, S13, S14, S15, S16, S17, S18, S19, S20, S21, S22, S23, S24, S25, S26, S27, S28, S29, S30, S31, S32, S33, S34, S35, S36, S37, S38, S39, S40, S41, S42, S43, S44):
    global Hh
    S = list(map(lambda x: x, [S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15, S16, S17, S18, S19, S20, S21, S22, S23, S24, S25, S26, S27, S28, S29, S30, S31, S32, S33, S34, S35, S36, S37, S38, S39, S40, S41, S42, S43, S44]))
    list1 = []
    for i in range(44):
        if S[i] == 1 or S[i] == 0:
            list1.append(i)

    for itn in range(10):
        new_S = np.zeros(44)
        new_S[0] = S[0]
        new_S[1] = S[1]
        new_S[2] = S[2]
        new_S[3] = S[3]
        new_S[4] = (1 - S[32]) * S[34] * (1 - S[41])
        new_S[5] = S[11] * (1 - S[37] + S[37]) + (1 - S[11]) * (1 - S[37])
        new_S[6] = S[29]
        new_S[7] = S[19]
        new_S[8] = S[4] * S[11] * S[34] + (1 - S[4]) * (S[11] * (1 - S[34]) + S[11] * S[34] + (1 - S[11]) * S[34])
        new_S[9] = (1 - S[4]) * S[31] * S[34] * S[40]
        new_S[10] = S[11]
        new_S[11] = (1 - S[10]) * (S[13] * S[15] + S[13] * (1 - S[15]) + (1 - S[13]) * S[15])
        new_S[12] = S[6]
        new_S[13] = (1 - S[17]) * (
                S[41] * (S[0] * S[1] + (1 - S[0]) * (1 - S[1]) + S[0] * (1 - S[1]) + (1 - S[0]) * S[1]) + (
                1 - S[41]) * (1 - S[0]) * (S[1] + (1 - S[1])) + (1 - S[41]) * S[0] * (1 - S[1]))
        new_S[14] = S[18]
        new_S[15] = S[12] * S[14] * (
                S[21] * S[22] + (1 - S[21]) * (1 - S[22]) + (1 - S[21]) * S[22] + S[21] * (1 - S[22])) + S[21] * S[
                        22] * (S[12] * S[14] + (1 - S[12]) * (1 - S[14]) + (1 - S[12]) * S[14] + S[12] * (1 - S[14]))

        new_S[16] = S[8] * S[5] * (1 - S[27]) * (S[26] * (1 - S[24]) + S[24] * (1 - S[26]) + S[24] * S[26])

        new_S[17] = S[21] * (
                (1 - S[4]) * (1 - S[16]) * (1 - S[22]) * (1 - S[7]) + (1 - S[4]) * (1 - S[16]) * (1 - S[22]) * S[
            7] + (1 - S[4]) * (1 - S[16]) * S[22] * (1 - S[7]) + (1 - S[4]) * (1 - S[16]) * S[22] * S[7] + (
                        1 - S[4]) * S[16] * (1 - S[22]) * (1 - S[7]) + (1 - S[4]) * S[16] * (1 - S[22]) * S[
                    7] + (1 - S[4]) * S[16] * S[22] * (1 - S[7]) + (1 - S[4]) * S[16] * S[22] * S[7] + S[4] * (
                        1 - S[16]) * (1 - S[22]) * (1 - S[7]) + S[4] * (1 - S[16]) * (1 - S[22]) * S[7] + S[
                    4] * (1 - S[16]) * S[22] * (1 - S[7]) + S[4] * (1 - S[16]) * S[22] * S[7] + S[4] * S[16] * (
                        1 - S[22]) * (1 - S[7]) + S[4] * S[16] * (1 - S[22]) * S[7] + S[4] * S[16] * S[22] * (
                        1 - S[7]) + S[4] * S[16] * S[22] * S[7]) + (1 - S[21]) * S[4] * (
                            (1 - S[16]) * (1 - S[22]) * (1 - S[7]) + (1 - S[16]) * (1 - S[22]) * S[7] + (
                            1 - S[16]) * S[22] * (1 - S[7]) + (1 - S[16]) * S[22] * S[7] + S[16] * (
                                    1 - S[22]) * (1 - S[7]) + S[16] * (1 - S[22]) * S[7] + S[16] * S[22] * (
                                    1 - S[7]) + S[16] * S[22] * S[7]) + (1 - S[21]) * (1 - S[4]) * (
                            1 - S[16]) * (
                            (1 - S[22]) * (1 - S[7]) + (1 - S[22]) * S[7] + S[22] * (1 - S[7]) + S[22] * S[7]) + (
                            1 - S[21]) * (1 - S[4]) * S[16] * (1 - S[22]) * (1 - S[7] + S[7]) + (1 - S[21]) * (
                            1 - S[4]) * S[16] * S[22] * S[7]

        new_S[18] = S[29]
        new_S[19] = S[7] * (S[42] * (1 - S[2]) + S[42] * S[2] + (1 - S[42]) * (1 - S[2]))
        new_S[20] = (1 - S[11]) * (1 - S[34]) * (1 - S[41])
        new_S[21] = S[35]
        new_S[22] = S[23]
        new_S[23] = S[3]
        new_S[24] = S[17] * ((1 - S[34]) * S[11] + (1 - S[34]) * (1 - S[11]) + (1 - S[11]) * S[34])
        new_S[25] = S[11]
        new_S[26] = S[17] * (
                S[34] * (S[41] * S[29] + (1 - S[41]) * S[29] + (1 - S[41]) * (1 - S[29]) + S[41] * (1 - S[29])) + (
                1 - S[34]) * (1 - S[41]) * (1 - S[29] + S[29]) + (1 - S[34]) * S[41] * (1 - S[29]))
        new_S[27] = S[33] * (1 - S[3]) * (1 - S[8])
        new_S[28] = S[38]
        new_S[29] = S[30] * S[28]
        new_S[30] = S[11]
        new_S[31] = S[3]
        new_S[32] = S[36]
        new_S[33] = (1 - S[3])
        new_S[34] = S[3]
        new_S[35] = S[3] * S[11]
        new_S[36] = S[19]
        new_S[37] = (1 - S[3]) * (1 - S[4])
        new_S[38] = S[3]
        new_S[39] = S[0]
        new_S[40] = S[32]
        new_S[41] = S[34] * S[3] * S[32]
        new_S[42] = S[43]
        new_S[43] = S[3]
        for j in range(len(list1)):
            new_S[list1[j]] = S[list1[j]]
        S = new_S
    if S[0] == 1 and S[1] == 1 and S[2] == 1 and S[3] == 1 and S[4] == 1 and S[5] == 1 and S[6] == 0 and S[7] == 0 and \
            S[8] == 0 and S[9] == 0 and S[10] == 0 \
            and S[11] == 0 and S[12] == 0 and S[13] == 0 and S[14] == 0 and S[15] == 0 and S[16] == 0 and S[17] == 1 and \
            S[18] == 0 and S[19] == 0 and S[20] == 0 \
            and S[21] == 0 and S[22] == 1 and S[23] == 1 and S[24] == 1 and S[25] == 0 and S[26] == 1 and S[27] == 0 and \
            S[28] == 1 and S[29] == 0 and S[30] == 0 \
            and S[31] == 1 and S[32] == 0 and S[33] == 0 and S[34] == 1 and S[35] == 0 and S[36] == 0 and S[37] == 0 and \
            S[38] == 1 and S[39] == 1 and S[40] == 0 and S[41] == 0 and S[42] == 1 and S[43] == 1:
        Hh = 0
    else:
        Hh = 1
    return Hh




def ql(epsilon):
    q_table = build_qtable(200000, ACTIONS)
    Reward_list = []
    rewards = []
    for episode in range(num_episodes):
        Sq = 0
        S = [0.5] * 44
        if episode % 500 == 0:
            epsilon = max(0.05, epsilon * 0.95)
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
            S30 = S[29]
            S31 = S[30]
            S32 = S[31]
            S33 = S[32]
            S34 = S[33]
            S35 = S[34]
            S36 = S[35]
            S37 = S[36]
            S38 = S[37]
            S39 = S[38]
            S40 = S[39]
            S41 = S[40]
            S42 = S[41]
            S43 = S[42]
            S44 = S[43]
            H = calculate_H(S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15, S16, S17, S18, S19, S20, S21, S22, S23, S24, S25, S26, S27, S28, S29, S30, S31, S32, S33, S34, S35, S36, S37, S38, S39, S40, S41, S42, S43, S44)
            C = calculate_len(S)
            R = reward(H,C)
            Reward_list.append(R)
            Sn = next_S(Sq, A, C, H)
            q_predict = q_table.loc[Sq, A]
            A1 = q_table.iloc[0, 0]
            A2 = q_table.iloc[0, 1]
            A3 = q_table.iloc[0, 2]
            A4 = q_table.iloc[0, 3]
            A5 = q_table.iloc[0, 4]
            A6 = q_table.iloc[0, 5]
            A7 = q_table.iloc[0, 6]
            A8 = q_table.iloc[58, 7]
            A9 = q_table.iloc[0, 8]
            A10 = q_table.iloc[0, 9]
            A11 = q_table.iloc[0, 10]
            A12 = q_table.iloc[100, 11]
            A13 = q_table.iloc[0, 0]
            A14 = q_table.iloc[0, 13]
            A15 = q_table.iloc[0, 14]
            A16 = q_table.iloc[74, 15]
            A17 = q_table.iloc[0, 16]
            A18 = q_table.iloc[0, 17]
            A19 = q_table.iloc[0, 18]
            A20 = q_table.iloc[19, 19]
            A21 = q_table.iloc[0, 20]
            A22 = q_table.iloc[0, 21]
            A23 = q_table.iloc[0, 22]
            A24 = q_table.iloc[0, 23]
            A25 = q_table.iloc[0, 24]
            A26 = q_table.iloc[0, 25]
            A27 = q_table.iloc[0, 26]
            A28 = q_table.iloc[0, 27]
            A29 = q_table.iloc[0, 28]
            A30 = q_table.iloc[0, 29]
            A31 = q_table.iloc[0, 30]
            A32 = q_table.iloc[0, 31]
            A33 = q_table.iloc[0, 32]
            A34 = q_table.iloc[0, 33]
            A35 = q_table.iloc[0, 34]
            A36 = q_table.iloc[0, 35]
            A37 = q_table.iloc[0, 36]
            A38 = q_table.iloc[0, 37]
            A39 = q_table.iloc[0, 38]
            A40 = q_table.iloc[0, 39]
            A41 = q_table.iloc[0, 40]
            A42 = q_table.iloc[0, 41]
            A43 = q_table.iloc[0, 42]
            A44 = q_table.iloc[0, 43]
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
            return_list30.append(A30)
            return_list31.append(A31)
            return_list32.append(A32)
            return_list33.append(A33)
            return_list34.append(A34)
            return_list35.append(A35)
            return_list36.append(A36)
            return_list37.append(A37)
            return_list38.append(A38)
            return_list39.append(A39)
            return_list40.append(A40)
            return_list41.append(A41)
            return_list42.append(A42)
            return_list43.append(A43)
            return_list44.append(A44)
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

