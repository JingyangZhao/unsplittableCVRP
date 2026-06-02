import math
from gurobipy import *
import itertools
import sys
import time
################################################################################################################################
def lp(N): # the factor-revealing LP in Theorem 3
    ############################################################################################################################
    m = Model('calculate ratio')
    y = m.addVar(vtype=GRB.CONTINUOUS, name="y")
    m.setParam("Threads", 4)
    m.setParam("Method", 2)
    m.setParam("BarOrder", 1)
    m.setParam("BarHomogeneous", 1)
    m.setParam("NumericFocus", 2)
    m.setParam("Crossover", 0)
    ############################################################################################################################
    M = 1000
    N = int(N)
    nums = [N//4, N//3, N//2]
    ############################################################################################################################
    xi0 = m.addVars(range(1, N+1), vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)
    xi1 = m.addVars(range(1, N+1), vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)  
    ############################################################################################################################
    ### Bounds in (3)
    m.addConstr(xi1[1] <= 1/N * xi0[1])
    for i in range(2, N+1):
        m.addConstr(xi0[i] - xi0[i-1] >= 0)
        m.addConstr(xi1[i] - xi1[i-1] >= (i-1)/N * (xi0[i] - xi0[i-1]))
        m.addConstr(xi1[i] - xi1[i-1] <= i/N * (xi0[i] - xi0[i-1]))
    m.addConstr(xi1[N] <= 1)
    ############################################################################################################################
    ### Vertex sets of canonical tours
    V = [
        [],
        [0, 1],             # Q1: <1>
        [0, 2],             # Q2: <2>
        [0, 3],             # Q3: <3>
        [0, 4, 5],          # Q4: <1,1>
        [0, 6, 7],          # Q5: <1,2>
        [0, 8, 9],          # Q6: <1,3>
        [0, 10, 11],        # Q7: <2,1>
        [0, 12, 13],        # Q8: <2,2>
        [0, 14, 15],        # Q9: <2,3>
        [0, 16, 17],        # Q10: <3,1>
        [0, 18, 19],        # Q11: <3,2>
        [0, 20, 21, 22],    # Q12: <1,1,1>
        [0, 23, 24, 25],    # Q13: <1,1,2>
        [0, 26, 27, 28],    # Q14: <1,2,1>
        [0, 29, 30, 31],    # Q15: <1,2,2>
        [0, 32, 33, 34],    # Q16: <2,1,1>
        [0, 35, 36, 37],    # Q17: <2,1,2>
        [0, 38, 39, 40]     # Q18: <2,2,1>
    ]
    ############################################################################################################################
    ### Edge sets of canonical tours
    ET = [[]]
    for Vi in V[1:]:
        ETi = []
        for i in range(len(Vi) - 1):
            ETi.append((Vi[i], Vi[i + 1]))
        ETi.append((Vi[-1], Vi[0]))  # return to the depot
        ET.append(ETi)
    Npat = len(V) - 1
    ############################################################################################################################
    ### Edge set of the graph \mathbb{G}
    edge = [(0, 0)]
    for i in range(1, 41):
        edge.append((0, i))
        edge.append((i, 0))
    for Vi in V[1:]:
        for i in Vi[1:]:
            for j in Vi[1:]:
                edge.append((i, j))
                
    chi = m.addVars(edge, vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)
    ############################################################################################################################
    ### Bounds in Lemma 13
    for Vi in V[1:]:
        for i in Vi: 
            m.addConstr(chi[i, i] == 0)
            for j in Vi:
                for k in Vi:
                    m.addConstr(chi[i, j] + chi[j, k] >= chi[i, k])
             
    m.addConstr(quicksum(quicksum(chi[i, j] for i, j in ET[k]) for k in range(1, Npat+1)) <= 1)
    ############################################################################################################################
    ### Bounds in Lemma 16
    J1 = {1, 4, 5, 6, 8, 11, 17, 20, 21, 22, 23, 24, 26, 28, 29, 33, 34, 36, 40}  # Type-1 vertices 
    J2 = {2, 7, 10, 12, 13, 14, 19, 25, 27, 30, 31, 32, 35, 37, 38, 39}          # Type-2 vertices
    J3 = {3, 9, 15, 16, 18}                                                      # Type-3 vertices

    m.addConstr(xi0[N] - xi0[N//2] >= quicksum(chi[0, i] + chi[i, 0] for i in J3))
    m.addConstr(xi0[N//2] - xi0[N//3] >= quicksum(chi[0, i] + chi[i, 0] for i in J2))
    m.addConstr(xi0[N//3] - xi0[N//4] >= quicksum(chi[0, i] + chi[i, 0] for i in J1))
    ############################################################################################################################
    ### Bound under \mu=1/2 in (7): Note that phi[N//2] represents \phi_{\gamma}(N/2)/\gamma in (7)
    phi = m.addVars(nums, vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)
    m.addConstr(phi[N//2] <= xi0[N] - xi0[N//2])
    ############################################################################################################################
    ### Lemma 14: \mu = 1/3
    m.addConstr(phi[N//3] <= 
        (chi[0, 2] + chi[2, 0]) + 
        (chi[0, 3] + chi[3, 0]) + 
        (chi[0, 7] + chi[7, 0]) + 
        (chi[0, 9] + chi[9, 0]) + 
        (chi[0, 10] + chi[10, 0]) + 
        (chi[0, 12] + chi[12, 13] + chi[13, 0]) + 
        (chi[0, 14] + chi[14, 15] + chi[15, 0]) + 
        (chi[0, 16] + chi[16, 0]) + 
        (chi[0, 18] + chi[18, 19] + chi[19, 0]) + 
        (chi[0, 25] + chi[25, 0]) + 
        (chi[0, 27] + chi[27, 0]) + 
        (chi[0, 30] + chi[30, 31] + chi[31, 0]) + 
        (chi[0, 32] + chi[32, 0]) + 
        (chi[0, 35] + chi[35, 37] + chi[37, 0]) + 
        (chi[0, 38] + chi[38, 39] + chi[39, 0])
    )
    ############################################################################################################################
    ### Lemma 15: constraints for \mu = 1/4 (linearization of the min operator via auxiliary variables and additional constraints)
    Gamma4 = m.addVars(range(1, 19), vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)
    
    # for the first 11 patterns，Gamma_i = the chi-cost of the corresponding canonical tour
    for i in range(1, 12):
        m.addConstr(Gamma4[i] == quicksum(chi[u, v] for u, v in ET[i]))
        
    # for the remaining patterns，Gamma_i is upper bounded by the chi-cost of the 3 permutations of the corresponding canonical tour
    for i in range(12, 19):
        p1, p2, p3 = V[i][1], V[i][2], V[i][3]
        # (p1, p2, p3)
        m.addConstr(Gamma4[i] <= chi[0, p1] + chi[p1, 0] + chi[0, p2] + chi[p2, p3] + chi[p3, 0])
        # (p2, p1, p3)
        m.addConstr(Gamma4[i] <= chi[0, p2] + chi[p2, 0] + chi[0, p1] + chi[p1, p3] + chi[p3, 0])
        # (p3, p1, p2)
        m.addConstr(Gamma4[i] <= chi[0, p3] + chi[p3, 0] + chi[0, p1] + chi[p1, p2] + chi[p2, 0])
        
    m.addConstr(phi[N//4] <= quicksum(Gamma4[i] for i in range(1, 19)))
    ############################################################################################################################
    xi0epsilonN = m.addVar(vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)
    xi1epsilonN = m.addVar(vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)
    m.addConstr(xi0[1] - xi0epsilonN >= 0)
    m.addConstr(xi1[1] - xi1epsilonN >= 0)
    m.addConstr(xi1[1] - xi1epsilonN <= 1/N * (xi0[1] - xi0epsilonN))
    
    for iM in range(M + 1): ### iN denotes I_M
        gamma = iM / M
        p = math.e**(-gamma) ### p denotes p_\gamma
        for iN in set(nums) | {-1}: ### iN denotes I_\delta and "iN = -1" means that I_\delta = \varepsilon\cdot N
            delta = iN / N
            if delta <= 1/2:
                for iiN in (x for x in set(nums) if x >= iN): ### iiN denotes I_\mu
                    if iN > 0:
                        if iiN == N/2:
                            m.addConstr(y <= gamma + (1/(1-delta))*xi1epsilonN + p*((1/(1-delta))*(xi1[iN]-xi1epsilonN)
                                            + 2/(1-delta)*(xi1[iiN]-xi1[iN]) - delta/(1-delta)*(xi0[iiN]-xi0[iN]))
                                            + p*phi[iiN])
                        else:
                            m.addConstr(y <= gamma + (1/(1-delta))*xi1epsilonN + p*((1/(1-delta))*(xi1[iN]-xi1epsilonN)
                                            + 2/(1-delta)*(xi1[iiN]-xi1[iN]) - delta/(1-delta)*(xi0[iiN]-xi0[iN]))
                                            + phi[iiN])
                    else:
                        if iiN == N/2:
                            m.addConstr(y <= gamma + xi1epsilonN + p*(2*(xi1[iiN]-xi1epsilonN)) + p*phi[iiN])
                        else:
                            m.addConstr(y <= gamma + xi1epsilonN + p*(2*(xi1[iiN]-xi1epsilonN)) + phi[iiN])
    ############################################################################################################################                        
    m.setObjective(y, GRB.MAXIMIZE)
    m.optimize()
    
    if m.status == GRB.OPTIMAL:
        return m.objVal
    else: return -1
    ############################################################################################################################

if __name__ == "__main__":
    N = int(sys.argv[1])

    start = time.time()
    value = lp(N)
    end = time.time()
    runtime = end - start

    print(f"N: {N}, value: {value}, runtime: {runtime} seconds")