import matplotlib.pyplot as plt
import numpy as np
import math
from fractions import Fraction
from gurobipy import *
import itertools
import os
import sys
import time
import fcntl
################################################################################################################################
def lp(generalflag, N, randomflag, metricflag):
    ############################################################################################################################
    m = Model('calculate ratio')
    y = m.addVar(vtype=GRB.CONTINUOUS, name="y"); 
    # m.setParam('OutputFlag', 0)
    m.setParam("Threads", 4)
    m.setParam("Method", 2)
    m.setParam("BarOrder", 1)
    m.setParam("BarHomogeneous", 1)
    m.setParam("NumericFocus", 3)
    m.setParam("Crossover", 0)
    # m.setParam("BarConvTol", 1e-12) # default 1e-8
    # m.setParam("FeasibilityTol", 1e-9) # default 1e-6
    # m.setParam("OptimalityTol", 1e-9) # default 1e-6
    ############################################################################################################################
    M = 1000
    Nmax = 4 # N should be divisible by (Nmax!)
    ############################################################################################################################
    generalflag = int(generalflag)
    N = int(N); nums = range(N//Nmax, N//2 + 1)
    randomflag = int(randomflag)
    metricflag = int(metricflag)
    ############################################################################################################################   
    xi0 = m.addVars(range(1, N+1), vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)
    xi1 = m.addVars(range(1, N+1), vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)  
    ############################################################################################################################   
    m.addConstr( xi1[1] <= 1/N * xi0[1] )
    for i in range(2, N+1):
        m.addConstr( xi0[i]-xi0[i-1] >= 0 )
        m.addConstr( xi1[i]-xi1[i-1] >= (i-1)/N * (xi0[i]-xi0[i-1]) )
        m.addConstr( xi1[i]-xi1[i-1] <= i/N * (xi0[i]-xi0[i-1]) )
    m.addConstr( xi1[N] <= 1 )
    ############################################################################################################################
    ### Enumerate all feasible patterns
    patterns = []
    for i in range(1, Nmax):
        if i == 1:
            for num in nums:
                patterns.append((num,))
        else:
            for val in itertools.combinations_with_replacement(nums, i):
                if sum(val) < N:
                    for p in itertools.permutations(val):
                        if metricflag == 1:
                            patterns.append(min(p, p[::-1])) # metric case
                        else: 
                            patterns.append(p) # asymmetric case
    patterns = list(dict.fromkeys(patterns))
    ############################################################################################################################
    vertex = [0]; nextvertex = 1
    mapvertex = [0]
    V=[[]]
    for perm in patterns:
        Vi=[0]
        for x in perm:
            Vi.append(nextvertex)
            vertex.append(nextvertex); mapvertex.append(x)
            nextvertex += 1
        V.append(Vi)
    ### Edge sets of canonical tours
    ET = [[]]
    for Vi in V[1:]:
        ETi = []
        for i in range(len(Vi)):
            ETi.append((Vi[i], Vi[(i+1)%(len(Vi))]))
        ET.append(ETi)
    Npat = len(V) - 1
    ############################################################################################################################
    edge = [(0,0)]
    for i in vertex[1:]:
        edge.append((0,i)); edge.append((i,0))
    for Vi in V[1:]:
        for i in Vi[1:]:
            for j in Vi[1:]:
                edge.append((i,j))

    chi = m.addVars(edge, vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)
    ############################################################################################################################
    for Vi in V:
        for i in Vi: 
            m.addConstr( 0 >= chi[i,i] )
            for j in Vi:
                if metricflag == 1: m.addConstr( chi[i,j] >= chi[j,i] ) ### metric case
                for k in Vi:
                    m.addConstr( chi[i,j]+chi[j,k] >= chi[i,k] )
                    
    m.addConstr( quicksum(quicksum(chi[i, j] for i,j in set(ET[k])) for k in range(1, Npat+1)) <= 1 )
    ############################################################################################################################
    m.addConstr( xi0[N] - xi0[N//2] >= quicksum( chi[0,i] + chi[i,0] for i in vertex if mapvertex[i] == N//2 ) )

    for k in range(len(nums)-1):
        m.addConstr( xi0[nums[k+1]] - xi0[nums[k]] >= quicksum( chi[0,i] + chi[i,0] for i in vertex if mapvertex[i] == nums[k]) )
    ############################################################################################################################
    phi1 = m.addVars(nums, vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)

    phi2 = m.addVars(nums, vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)
    phi2i = m.addVars(list(itertools.product(nums, range(1, Npat+1))), vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)
    m.addConstrs( phi2[j] >= quicksum( phi2i[j, k] for k in range(1, Npat+1)) for j in nums )
    
    phi3 = m.addVars(nums, vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)
    phi3i = m.addVars(list(itertools.product(nums, range(1, Npat+1))), vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)

    m.addConstrs( phi3[j] >= quicksum( phi3i[j, k] for k in range(1, Npat+1)) for j in nums )

    if randomflag == 1:
        phi4 = m.addVars(nums, vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)
        phi4i = m.addVars(list(itertools.product(nums, range(1, Npat+1))), vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)

        m.addConstrs( phi4[j] >= quicksum( phi4i[j, k] for k in range(1, Npat+1)) for j in nums )

    for atleastdemand in nums:
        m.addConstr( phi1[atleastdemand] <= quicksum(chi[0, i] + chi[i, 0] for i in vertex if mapvertex[i] >= atleastdemand) )

        for k in range(1, Npat+1):
            bigcustomer = []
            for i in V[k][1:]:
                if mapvertex[i] >= atleastdemand:
                    bigcustomer.append(i)

            if len(bigcustomer) == 2:
                m.addConstr( phi2i[atleastdemand, k] >= chi[bigcustomer[0],0]+chi[0,bigcustomer[1]] - chi[bigcustomer[0],bigcustomer[1]] )
            elif len(bigcustomer) == 3:
                m.addConstr( phi3i[atleastdemand, k] >= chi[bigcustomer[0],0]+chi[0,bigcustomer[1]] - chi[bigcustomer[0],bigcustomer[1]] )
                m.addConstr( phi3i[atleastdemand, k] >= chi[bigcustomer[0],0]+chi[0,bigcustomer[2]] - chi[bigcustomer[0],bigcustomer[2]] )
                m.addConstr( phi3i[atleastdemand, k] >= chi[bigcustomer[1],0]+chi[0,bigcustomer[2]] - chi[bigcustomer[1],bigcustomer[2]] )

                if randomflag == 1:
                    m.addConstr( phi4i[atleastdemand, k] >= 2*chi[bigcustomer[0],0]+2*chi[0,bigcustomer[2]]
                                                                    + chi[bigcustomer[1],0]+chi[0,bigcustomer[1]]
                                                                    - chi[bigcustomer[0],bigcustomer[1]]-chi[bigcustomer[0],bigcustomer[2]]-chi[bigcustomer[1],bigcustomer[2]] )
    ############################################################################################################################
    if generalflag == 0:
        for iM in range(M+1):
            gamma = iM/M
            p = math.e**(-gamma)
            for iN in set(nums):
                delta = iN / N
                if delta <= 1/3:
                    for iiN in (x for x in set(nums) if x >= iN):
                        if randomflag == 0:
                            m.addConstr( y <= gamma + p*((1/(1-delta))*xi1[iN]
                                                            + 2/(1-delta)*(xi1[iiN]-xi1[iN]) - delta/(1-delta)*(xi0[iiN]-xi0[iN]))
                                                            + p*phi1[iiN] - max(0, 2*p-1)*phi2[iiN] - max(0, 2*p-1)*phi3[iiN] )                            
                        else:
                            if p <= 1/2:
                                m.addConstr( y <= gamma + p*((1/(1-delta))*xi1[iN]
                                                                + 2/(1-delta)*(xi1[iiN]-xi1[iN]) - delta/(1-delta)*(xi0[iiN]-xi0[iN]))
                                                                + p*phi1[iiN] )
                            elif p <= 2/3:
                                    m.addConstr( y <= gamma + p*((1/(1-delta))*xi1[iN]
                                                                    + 2/(1-delta)*(xi1[iiN]-xi1[iN]) - delta/(1-delta)*(xi0[iiN]-xi0[iN]))
                                                                    + p*phi1[iiN] - (2*p-1)*phi2[iiN] - (2*p-1)*phi3[iiN] )

                                    m.addConstr( y <= gamma + p*((1/(1-delta))*xi1[iN]
                                                                    + 2/(1-delta)*(xi1[iiN]-xi1[iN]) - delta/(1-delta)*(xi0[iiN]-xi0[iN]))
                                                                    + p*phi1[iiN] - (2*p-1)*phi2[iiN] - (2*p-1)/2*phi4[iiN] )
                            else:
                                    m.addConstr( y <= gamma + p*((1/(1-delta))*xi1[iN]
                                                                    + 2/(1-delta)*(xi1[iiN]-xi1[iN]) - delta/(1-delta)*(xi0[iiN]-xi0[iN]))
                                                                    + p*phi1[iiN] - (2*p-1)*phi2[iiN] - (2*p-1)*phi3[iiN] )
                                    
                                    m.addConstr( y <= gamma + p*((1/(1-delta))*xi1[iN]
                                                                    + 2/(1-delta)*(xi1[iiN]-xi1[iN]) - delta/(1-delta)*(xi0[iiN]-xi0[iN]))
                                                                    + p*phi1[iiN] - (2*p-1)*phi2[iiN] - (3*p-2)*phi3[iiN] - (1-p)/2*phi4[iiN] )
                                    
                                    m.addConstr( y <= gamma + p*((1/(1-delta))*xi1[iN]
                                                                    + 2/(1-delta)*(xi1[iiN]-xi1[iN]) - delta/(1-delta)*(xi0[iiN]-xi0[iN]))
                                                                    + p*phi1[iiN] - (2*p-1)*phi2[iiN] - (3*p-1)/6*phi4[iiN] )
    ############################################################################################################################
    else:
        xi0epsilonN = m.addVar(vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)
        xi1epsilonN = m.addVar(vtype=GRB.CONTINUOUS, lb=0.0, ub=GRB.INFINITY)
        m.addConstr( xi0[1] - xi0epsilonN >= 0 )
        m.addConstr( xi1[1] - xi1epsilonN >= 0 )
        m.addConstr( xi1[1] - xi1epsilonN <= 1/N*(xi0[1] - xi0epsilonN) )
        
        for iM in range(M+1):
            gamma = iM/M
            p = math.e**(-gamma)
            for iN in set(nums) | {-1}:
                delta = iN / N
                if delta <= 1/3:
                    for iiN in (x for x in set(nums) if x >= iN):
                        if iN > 0:
                            if randomflag == 0:
                                m.addConstr( y <= gamma + (1/(1-delta))*xi1epsilonN + p*((1/(1-delta))*(xi1[iN]-xi1epsilonN)
                                                                + 2/(1-delta)*(xi1[iiN]-xi1[iN]) - delta/(1-delta)*(xi0[iiN]-xi0[iN]))
                                                                + p*phi1[iiN] - max(0, 2*p-1)*phi2[iiN] - max(0, 2*p-1)*phi3[iiN] )
                            else:
                                if p <= 1/2:
                                    m.addConstr( y <= gamma + (1/(1-delta))*xi1epsilonN + p*((1/(1-delta))*(xi1[iN]-xi1epsilonN)
                                                                    + 2/(1-delta)*(xi1[iiN]-xi1[iN]) - delta/(1-delta)*(xi0[iiN]-xi0[iN]))
                                                                    + p*phi1[iiN] )
                                elif p <= 2/3:
                                        m.addConstr( y <= gamma + (1/(1-delta))*xi1epsilonN + p*((1/(1-delta))*(xi1[iN]-xi1epsilonN)
                                                                        + 2/(1-delta)*(xi1[iiN]-xi1[iN]) - delta/(1-delta)*(xi0[iiN]-xi0[iN]))
                                                                        + p*phi1[iiN] - (2*p-1)*phi2[iiN] - (2*p-1)*phi3[iiN] )
                                        
                                        m.addConstr( y <= gamma + (1/(1-delta))*xi1epsilonN + p*((1/(1-delta))*(xi1[iN]-xi1epsilonN)
                                                                        + 2/(1-delta)*(xi1[iiN]-xi1[iN]) - delta/(1-delta)*(xi0[iiN]-xi0[iN]))
                                                                        + p*phi1[iiN] - (2*p-1)*phi2[iiN] - (2*p-1)/2*phi4[iiN] )
                                else:
                                        m.addConstr( y <= gamma + (1/(1-delta))*xi1epsilonN + p*((1/(1-delta))*(xi1[iN]-xi1epsilonN)
                                                                        + 2/(1-delta)*(xi1[iiN]-xi1[iN]) - delta/(1-delta)*(xi0[iiN]-xi0[iN]))
                                                                        + p*phi1[iiN] - (2*p-1)*phi2[iiN] - (2*p-1)*phi3[iiN] )
                                        
                                        m.addConstr( y <= gamma + (1/(1-delta))*xi1epsilonN + p*((1/(1-delta))*(xi1[iN]-xi1epsilonN)
                                                                        + 2/(1-delta)*(xi1[iiN]-xi1[iN]) - delta/(1-delta)*(xi0[iiN]-xi0[iN]))
                                                                        + p*phi1[iiN] - (2*p-1)*phi2[iiN] - (3*p-2)*phi3[iiN] - (1-p)/2*phi4[iiN] )
                                        
                                        m.addConstr( y <= gamma + (1/(1-delta))*xi1epsilonN + p*((1/(1-delta))*(xi1[iN]-xi1epsilonN)
                                                                        + 2/(1-delta)*(xi1[iiN]-xi1[iN]) - delta/(1-delta)*(xi0[iiN]-xi0[iN]))
                                                                        + p*phi1[iiN] - (2*p-1)*phi2[iiN] - (3*p-1)/6*phi4[iiN] )
                        else:
                            if randomflag == 0:
                                m.addConstr( y <= gamma + xi1epsilonN + p*(2*(xi1[iiN]-xi1epsilonN))
                                                                + p*phi1[iiN] - max(0, 2*p-1)*phi2[iiN] - max(0, 2*p-1)*phi3[iiN] )
                            else:
                                if p <= 1/2:
                                    m.addConstr( y <= gamma + xi1epsilonN + p*(2*(xi1[iiN]-xi1epsilonN))
                                                                    + p*phi1[iiN] )
                                elif p <= 2/3:
                                    m.addConstr( y <= gamma + xi1epsilonN + p*(2*(xi1[iiN]-xi1epsilonN))
                                                                    + p*phi1[iiN] - (2*p-1)*phi2[iiN] - (2*p-1)*phi3[iiN] )

                                    m.addConstr( y <= gamma + xi1epsilonN + p*(2*(xi1[iiN]-xi1epsilonN))
                                                                    + p*phi1[iiN] - (2*p-1)*phi2[iiN] - (2*p-1)/2*phi4[iiN] )
                                else:
                                    m.addConstr( y <= gamma + xi1epsilonN + p*(2*(xi1[iiN]-xi1epsilonN))
                                                                    + p*phi1[iiN] - (2*p-1)*phi2[iiN] - (2*p-1)*phi3[iiN] )
                                    
                                    m.addConstr( y <= gamma + xi1epsilonN + p*(2*(xi1[iiN]-xi1epsilonN))
                                                                    + p*phi1[iiN] - (2*p-1)*phi2[iiN] - (3*p-2)*phi3[iiN] - (1-p)/2*phi4[iiN] )
                                    
                                    m.addConstr( y <= gamma + xi1epsilonN + p*(2*(xi1[iiN]-xi1epsilonN))
                                                                    + p*phi1[iiN] - (2*p-1)*phi2[iiN] - (3*p-1)/6*phi4[iiN] )
    ############################################################################################################################
    m.setObjective(y, GRB.MAXIMIZE); m.optimize()
    if m.status == GRB.OPTIMAL:
        return m.objVal
    else: return -1

if __name__ == "__main__":
    generalflag = int(sys.argv[1])
    N = int(sys.argv[2])
    randomflag = int(sys.argv[3])

    start = time.time()

    res = lp(generalflag, N, randomflag)

    end = time.time()
    runtime = end - start

    with open("results.txt", "a") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        f.write(f"{generalflag} {N} {randomflag} {res} {runtime}\n")
        f.flush()
        fcntl.flock(f, fcntl.LOCK_UN)
    