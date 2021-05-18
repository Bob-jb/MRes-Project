import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

def gen_ising(adj):
    #Create a hardware native spin glass of the form 
    # \sum_{i,j} J_{i,j} Z_i Z_j with J_{i,j} being 1 or -1 
    # from a numpy array or list of lists. Returns a numpy array

    adj=np.array(adj)
    compj=np.random.choice([-1,1],size=(len(adj),len(adj)))
    j=np.multiply(compj,adj)
    j=np.transpose(np.triu(j))+np.triu(j)

    return j

def plot_ising(adj):
    # Create plot of Ising model
    
    G=nx.from_numpy_array(adj)
    pos = nx.spring_layout(G, seed=63)
    nx.draw(G, pos, with_labels=True, font_weight='bold')
    nx.draw_networkx_edge_labels(G, pos)
    plt.show()

def solve_ising(adj):
    #Solve Ising 


mat=[[0,1],[1,0]]
y=gen(mat)
print(y)
plot_ising(y)
