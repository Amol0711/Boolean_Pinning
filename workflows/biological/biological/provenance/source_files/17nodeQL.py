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
pd.set_option('display.max_columns',None)
pd.set_option('display.max_rows',None)
np.random.seed()
ACTIONS = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
Q_list = ['22222222222222222']
epsilon = 0.9
alpha = 0.05
gamma = 0.9
num_episodes = 80000
S = [0.5] * 17
S_counter = 0
attractor = ['1','1','0','0','0','0','0','0','0','0','1','0','0','1','1','0','1']
def calculate_len(list):
    L = []
    for i in range(10):
        if list[i] == 1 or list[i] == 0:
            L.append(i)
    return len(L)
def action_func(action_num,list):
    if action_num == 1:
        list[0] = 1
    elif action_num == 2:
        list[1] = 1
    elif action_num == 3:
        list[2] = 0
    elif action_num == 4:
        list[3] = 0
    elif action_num == 5:
        list[4] = 0
    elif action_num == 6:
        list[5] = 0
    elif action_num == 7:
        list[6] = 0
    elif action_num == 8:
        list[7] = 0
    elif action_num == 9:
        list[8] = 0
    elif action_num == 10:
        list[9] = 0
    elif action_num == 11:
        list[10] = 1
    elif action_num == 12:
        list[11] = 0
    elif action_num == 13:
        list[12] = 0
    elif action_num == 14:
        list[13] = 1
    elif action_num == 15:
        list[14] = 1
    elif action_num == 16:
        list[15] = 0
    elif action_num == 17:
        list[16] = 1
    return list
#计算H的值
def calculate_H(S1,S2,S3,S4,S5,S6, S7, S8, S9, S10, S11, S12, S13, S14, S15, S16, S17):
    global Hh
    S = list(map(lambda x: x, [S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15, S16, S17]))
    list1 = []
    for i in range(17):
        if S[i] == 1 or S[i] == 0:
            list1.append(i)
    for itn in range(10):
        new_S = np.zeros(17)
        new_S[0] = S[0]
        new_S[1] = S[1]
        new_S[2] = S[2]
        new_S[3] = (1 - S[16]) * (
                    S[0] * (1 - S[3]) * S[15] + S[0] * S[3] * S[15] + S[0] * S[3] * (1 - S[15]) + (1 - S[0]) * S[3] * S[
                15])
        new_S[4] = S[3]
        new_S[5] = S[1] * (1 - S[0])
        new_S[6] = S[5]
        new_S[7] = S[6] * (1 - S[16])
        new_S[8] = S[7]
        new_S[9] = S[15] * (1 - S[6]) * (1 - S[16])
        new_S[10] = S[10] * (1 - S[9]) * (1 - S[2]) + S[9] * (
                    (1 - S[10]) * (1 - S[2]) + (1 - S[10]) * S[2] + S[10] * (1 - S[2]) + S[10] * S[2])
        new_S[11] = S[10] * S[2]
        new_S[12] = S[2] * (1 - S[10]) + S[2] * S[10] + (1 - S[10]) * (1 - S[2])
        new_S[13] = 1 - S[6]
        new_S[14] = S[13]
        new_S[15] = S[14] * (S[12] * S[2] + S[12] * (1 - S[2]) + (1 - S[12]) * S[2])
        new_S[16] = S[14] * (1 - S[12]) * (1 - S[2])
        for j in range(len(list1)):
            new_S[list1[j]] = S[list1[j]]
        S = new_S

    if S[0] == 1 and S[1] == 1 and S[2] == 0 and S[3] == 0 and S[4] == 0 and S[5] == 0 and S[6] == 0 and S[7] == 0 and \
            S[8] == 0 and S[9] == 0 and S[10] == 1 \
            and S[11] == 0 and S[12] == 0 and S[13] == 1 and S[14] == 1 and S[15] == 0 and S[16] == 1:
    #if Entropy(S) == 0:
        H = 0
    else:
        H = 1
    return H

def Entropy(S):
    h = []
    for i in range(17):
        if S[i] == 0 or S[i] == 1:
            h.append(0)
        else:
            h.append(-S[i]*np.log2(S[i])-(1-S[i])*np.log2(1-S[i]))
    return np.sum(h)/len(S)
def calculate_len(list):
    L = []
    for i in range(17):
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
        return 10*(17-C)
    else:
        return -1

def next_S(state,action,number, H):
    S_q = list(Q_list[state])
    A = int(action)-1
    if H == 0 or number > 5:
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
    q_table = build_qtable(50000,ACTIONS)
    Reward_list = []
    rewards = []
    for episode in range(num_episodes):
        Sq = 0
        S = [0.5] * 17
        if episode % 500 == 0:
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
            H = calculate_H(S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15, S16, S17)
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
            A8 = q_table.iloc[0, 7]
            A9 = q_table.iloc[0, 8]
            A10 = q_table.iloc[0, 9]
            A11 = q_table.iloc[0, 10]
            A12 = q_table.iloc[0, 11]
            A13 = q_table.iloc[0, 12]
            A14 = q_table.iloc[0, 13]
            A15 = q_table.iloc[0, 14]
            A16 = q_table.iloc[0, 15]
            A17 = q_table.iloc[0, 16]
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
plt.legend(ncol=1, loc='upper right')
plt.title('')
plt.xlabel('Episodes')
plt.ylabel('Rewards')
plt.show()


