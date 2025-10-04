import numpy as np 
import cvxpy as cvx
from scipy.signal import StateSpace, cont2discrete 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

"""
#Re-Do the function I did in MATLAB: 
def get_delayed_dtsys(h, tau):
    Fx = np.array([
        [np.exp(2.3 * h), 0],
        [5*np.exp(2.5*h) - 5 * np.exp(2.3 * h), np.exp(2.5 * h)]
                   ])
    
    Fu = np.array([
        [(1/2.3)*(np.exp(2.3 * h) - np.exp(2.3 * (h-tau)))],
        [2 * (np.exp(2.5 * h) - np.exp(2.5 * (h - tau))) + (5/2.3)*(np.exp(2.3 * (h - tau)) - np.exp(2.3 * h))]
        ])

    G1 = np.array([
        [0.435 * np.exp(2.3*(h -tau)) - 0.4348],
        [2 * np.exp(2.5 * (h - tau)) - 2.174 * np.exp(2.3 * (h - tau)) + 0.174]
        ])

    #Extend the state to include the previous input: 
    # Stack Fx (2×2) and Fu (2×1) side by side, then append a 1×3 zero‐row
    F = np.vstack([
        np.hstack([Fx, Fu]),
        np.zeros((1, Fx.shape[1] + Fu.shape[1]))
    ])
    
    G = np.vstack((G1, 1))
    
    return F, G 

F, G = get_delayed_dtsys(0.5, 0.5)
print(f"F = {F}\n\nG={G}")
"""


def get_polytopic_vertex(h, delta1, delta2):
    F0 = np.array([[np.exp(2.3 * h), 0, (1/2.3) * np.exp(2.3 * h)],
                  [5* np.exp(2.5 * h) -  5* np.exp(2.3 * h), np.exp(2.5*h), (2 * np.exp(2.5 *h) - (5/2.3) * np.exp(2.3 *h))],
                  [0, 0, 0]])
    
    F1 = np.array([[0,0,(-1/2.3)],
                   [0,0,(5/2.3)],
                   [0,0,0]])
    
    F2 = np.array([[0,0,0],
                   [0,0,-2],
                   [0,0,0]])
    
    G0 = np.array([[-0.4348],
                   [0.174],
                   [1]])
    
    G1 = np.array([[0.435],
                   [-2.174],
                   [0]])
    
    G2 = np.array([[0],
                   [2],
                   [0]])
    
    HF = F0 + delta1 * F1 + delta2 * F2
    HG = G0 + delta1 * G1 + delta2 * G2

    return HF, HG

def get_alpha1(h, tau_k):
    return np.exp(2.3 * (h - tau_k))

def get_alpha2(h, tau_k):
    return np.exp(2.5 * (h - tau_k))

#Method from LMI Lecture
epsilon = 2e-6

K = np.array([[8.8, 21.25, 0]])

#Range of hs to iterate through
h_range = np.linspace(0.01, 0.5, 100)
is_feasible = np.zeros((len(h_range), len(h_range)))
for i in range(len(h_range)):
    #get the h that we're testing at
    h = h_range[i]

    #get the min and max for tau_k
    tau_min = 0
    tau_max = h
    
    #get the max and min alphas
    alpha1_underbar = get_alpha1(h, tau_max) #smallest alpha1 
    alpha1_bar = get_alpha1(h, tau_min) #biggest alpha1 
    alpha2_underbar = get_alpha2(h, tau_max) #smallest alpha2 
    alpha2_bar = get_alpha2(h, tau_min) #biggest alpha2

    #get all the combinations of Hf and Hg: 
    Hf1, Hg1 = get_polytopic_vertex(h, delta1=alpha1_underbar, delta2=alpha2_underbar)
    Hf2, Hg2 = get_polytopic_vertex(h, delta1=alpha1_underbar, delta2=alpha2_bar)
    Hf3, Hg3 = get_polytopic_vertex(h, delta1=alpha1_bar, delta2=alpha2_underbar)
    Hf4, Hg4 = get_polytopic_vertex(h, delta1=alpha1_bar, delta2=alpha2_bar)

    #Declare the P variables
    P = cvx.Variable((3,3), PSD=True)
    a = cvx.Parameter(nonneg=True, value=1e-3)  # fixed parameter with assigned value

    # Define the LMIs constituting the constraints of the optimization problem (including enforcing strict positive definite-ness on P matrices) 
    cons = [  # a is now a Parameter, so its bound is enforced externally
            P - epsilon * np.eye(3) >> 0,
            a * np.eye(2) << np.eye(2),
            (Hf1 - Hg1 @ K).T @ P @ (Hf1 - Hg1 * K) - P <= -a*P,
            (Hf2 - Hg2 @ K).T @ P @ (Hf2 - Hg2 * K) - P <= -a*P,
            (Hf3 - Hg3 @ K).T @ P @ (Hf3 - Hg3 * K) - P <= -a*P,
            (Hf4 - Hg4 @ K).T @ P @ (Hf4 - Hg4 * K) - P <= -a*P ]
    
    # Solve the system of LMIs
    prob = cvx.Problem(cvx.Minimize(0), cons)
    prob.solve(solver=cvx.SCS)
    #print(prob.status)
    #check if the solution is actually feasible(aka solution can be found, aka system is stable)
    if prob.status == cvx.OPTIMAL or prob.status.lower() == "optimal":
        #for infeasibility, mark as zero 
        is_feasible[1:i,i] = 1
    else: 
        is_feasible[1:i,i] = 0.5

# Create grid for heatmap
X, Y = np.meshgrid(h_range, h_range)

fig, ax = plt.subplots(figsize=(8,6))
pcm = ax.pcolormesh(X, Y, is_feasible, shading="auto", cmap='inferno')

# Add a colorbar with custom ticks and labels
cbar = fig.colorbar(pcm, ax=ax, ticks=[0, 0.5, 1])
cbar.ax.set_yticklabels(['tau>h', 'Unstable', 'Stable'])
cbar.set_label('Feasibility')

ax.set_title("Feasibility of LMI System vs Sampling Time", fontsize=16)
ax.set_xlabel("h [s]", fontsize=14)
ax.set_ylabel("Tau [s]", fontsize=14)
ax.set_xlim((0.01, 0.5))
ax.set_ylim((0.01, 0.5))
ax.tick_params(axis='both', labelsize=12)

# enlarge colorbar ticks and label
cbar.ax.tick_params(labelsize=12)

plt.savefig("Figures/q3_jordan.pdf")
plt.show()

