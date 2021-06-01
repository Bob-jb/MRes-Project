import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

def gen_ising(adj,local_t=False):
    #Create a hardware native spin glass of the form 
    # \sum_{i,j} J_{i,j} Z_i Z_j+h_i Z_i with J_{i,j} being 1 or -1 
    # from a numpy array or list of lists. Returns a numpy array

    adj=np.array(adj)
    compj=np.random.choice([-1,1],size=(len(adj),len(adj)))
    j=np.multiply(compj,adj)
    j=np.transpose(np.triu(j))+np.triu(j)

    if local_t:
        j=j+np.diag([np.random.choice([-1,1]) for value in np.diag(j)])
        

    return j

def plot_ising(adj):
    # Create plot of Ising model
    
    G=nx.from_numpy_array(adj)
    pos = nx.shell_layout(G)
    nx.draw(G, pos, with_labels=True, font_weight='bold')
    nx.draw_networkx_edge_labels(G, pos)
    plt.show()

def solve_ising(adj):
    #Solve Ising by exhaustion
    n=len(adj)
    h=np.diag(adj)
    adj=np.array(adj/2)-np.diag(h)/2
    print(adj)
    costs=[]

    for k in range(2**n):
        bits=list(format(k, '0'+str(n)+'b'))
        zeig=np.array([1-2*int(x) for x in bits])
        value=np.dot(zeig,np.dot(adj,zeig))+np.dot(np.array(h),zeig)
        print(np.dot(zeig,np.dot(adj,zeig)))
        costs.append(value)
        
    min_value=min(costs)
    place=[format(i, '0'+str(n)+'b') for i,j in enumerate(costs) if j == min_value]

    return place






#mat=[[0,1,0,0,0,0,0],
#     [1,0,1,1,0,0,0],
#     [0,1,0,0,0,0,0],
#     [0,1,0,0,0,1,0],
#     [0,0,0,0,0,1,0],
#     [0,0,0,1,1,0,1],
#     [0,0,0,0,0,1,0]]
mat=[[0,1],[1,0]]

y=gen_ising(mat,True)

y=np.array([[1,1],[1,-1]])
print(y)
print(solve_ising(y))