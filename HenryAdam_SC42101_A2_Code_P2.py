import numpy as np 
import cvxpy as cvx
from scipy.signal import StateSpace, cont2discrete 
import matplotlib.pyplot as plt

#since I can't isntall cvx on matlab, I'm just gonna write this one with python
#rewrite all the setup from MATLAB
a = 6; #elements of student ID
b = 4; 
c = 2; 

#Continuous system setup 
A = np.array([[0, 0.5-c], [0.2 + a - b, -1]])
B = np.array([[1], [0]])
C = np.eye(2)
D = np.zeros([2,1])

K1 = np.array([3.0000, -0.5909])
K2 = np.array([0.4000, 0.6000])

# Statistical properties: 
q=0.3 
p= 0.2

#Method from LMI Lecture
epsilon = 1e-6

#Range of hs to iterate through
h_range = np.linspace(0.01, 10, 1000)
is_feasible = np.ones(np.shape(h_range))
for i in range(len(h_range)):
    #This will all go inside the 
    h1 = h_range[i]
    h2 = 2 * h1 
    # Discretize using the state‐space matrices directly
    sys_d1 = cont2discrete((A, B, C, D), h1)
    F1, G1, _, _, _ = sys_d1
    sys_d2 = cont2discrete((A, B, C, D), h2)
    F2, G2, _, _, _ = sys_d2

    A1 = (F1 - G1 * K1)
    A2 = (F1 - G1 * K2)
    A3 = (F2 - G2 * K2)
    A4 = F2

    #Declare the P variables
    P1 = cvx.Variable((2,2), symmetric=True)
    P2 = cvx.Variable((2,2), symmetric=True)
    P3 = cvx.Variable((2,2), symmetric=True)
    P4 = cvx.Variable((2,2), symmetric=True)

    # Define the LMIs constituting the constraints of the optimization problem (including enforcing strict positive definite-ness on P matrices) 
    cons = [P1 - epsilon * np.eye(2) >> 0,
            P2 - epsilon * np.eye(2) >> 0,
            P3 - epsilon * np.eye(2) >> 0,  
            P4 - epsilon * np.eye(2) >> 0,      
            P1 - q * A3.T @ P3 @ A3 - (1-q) * A4.T @ P4 @ A4 >= epsilon * np.eye(2), 
            P2 - q * A3.T @ P3 @ A3 - (1-q) * A4.T @ P4 @ A4 >= epsilon * np.eye(2), 
            P3 - p * A1.T @ P1 @ A1 - (1-p) * A2.T @ P2 @ A2 >= epsilon * np.eye(2), 
            P4 - p * A1.T @ P1 @ A1 - (1-p) * A2.T @ P2 @ A2 >= epsilon * np.eye(2)]
    
    # Solve the system of LMIs
    prob = cvx.Problem(cvx.Minimize(0), cons)
    prob.solve(solver=cvx.SCS)
    #print(prob.status)
    #check if the solution is actually feasible(aka solution can be found, aka system is stable)
    if prob.status == cvx.INFEASIBLE or prob.status.lower() == "infeasible":
        #for infeasibility, mark as zero 
        is_feasible[i] = 0


plt.figure()
plt.plot(h_range, is_feasible)
plt.title("Feasibility of LMI System vs Sampling Time")
plt.xlabel("h [s]")
plt.xticks(np.arange(0, 11))
plt.yticks([0, 1], ["Unstable", "Stable"])
plt.ylim(-0.1, 1.1)
plt.savefig("Figures/q2_stability.pdf")
plt.show()


