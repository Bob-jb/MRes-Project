import numpy as np
import networkx as nx
import math
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
        costs.append(value)
        
    min_value=min(costs)
    place=[format(i, '0'+str(n)+'b') for i,j in enumerate(costs) if j == min_value]

    return place



def cost_function(beta,gamma):
    return (beta-2)**2+(gamma-5)**2

def exhaustive_search(resolution=30):
    step=list(range(resolution))
    dbeta=[b/resolution*math.pi for b in step]
    dgamma=[2*g/resolution*math.pi for g in step]
    cost=[]
    for gvalue in dgamma:
        cost.append([cost_function(bvalue,gvalue) for bvalue in dbeta])
    
    cost=np.array(cost)
    [m,n] = np.shape(cost)
    plt.figure()
    plt.imshow(cost, alpha=0.8)
    plt.xticks(np.arange(n))
    plt.yticks(np.arange(m))
    plt.xlabel('Numbers')
    plt.ylabel('Value')
    plt.title('Color Maps')
    plt.show()

    min=np.min(cost)
    place=[(dbeta[m],dgamma[i]) for i,j in enumerate(cost) for m,n in enumerate(j) if n == min]
    print(place)
    return place



def gradient(bi=0,gi=0,threshold=0.00001,max_iterations=500,learning_rate=0.05,momentum=0):
    i=0
    diff=1.0e10
    b=bi
    g=gi
    cbgn= cost_function(b,g)

    while i<=max_iterations and diff>threshold:

       cbg=cbgn
       d = [cost_function(b+learning_rate,g) - cbg,
       cost_function(b,g+learning_rate) - cbg]
       delta=[-1*dvalue+momentum*dvalue for dvalue in d]
       b=b+delta[0]
       g=g+delta[1]
       cbgn= cost_function(b,g)
       i=i+1
       diff=np.abs(cbg-cbgn)


    print(b,g)
    return((b,g))



mat=[[0,1,0,0,0,0,0],
     [1,0,1,1,0,0,0],
     [0,1,0,0,0,0,0],
     [0,1,0,0,0,1,0],
     [0,0,0,0,0,1,0],
     [0,0,0,1,1,0,1],
     [0,0,0,0,0,1,0]]
#mat=[[0,1],[1,0]]

#y=gen_ising(mat,True)
#print(y)
#print(solve_ising(y))

gradient()