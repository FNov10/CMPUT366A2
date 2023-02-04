from search.algorithms import State
from search.map import Map
import getopt
import sys
import heapq
def main():
    """
    Function for testing your implementation. Run it with a -help option to see the options available. 
    """
    optlist, _ = getopt.getopt(sys.argv[1:], 'h:m:r:', ['testinstances', 'plots', 'help'])

    plots = False
    for o, _ in optlist:
        if o in ("-help"):
            print("Examples of Usage:")
            print("Solve set of test instances: main.py --testinstances")
            print("Solve set of test instances and generate plots: main.py --testinstances --plots")
            exit()
        elif o in ("--plots"):
            plots = True
    test_instances = "test-instances/testinstances.txt"
    gridded_map = Map("dao-map/brc000d.map")
    
    nodes_expanded_biastar = []   
    nodes_expanded_astar = []   
    nodes_expanded_mm = []
    
    start_states = []
    goal_states = []
    solution_costs = []
       
    file = open(test_instances, "r")
    for instance_string in file:
        list_instance = instance_string.split(",")
        start_states.append(State(int(list_instance[0]), int(list_instance[1])))
        goal_states.append(State(int(list_instance[2]), int(list_instance[3])))
        
        solution_costs.append(float(list_instance[4]))
    file.close()
        
    for i in range(0, len(start_states)):   

        start = start_states[i]
        goal = goal_states[i]
    
        cost, expanded_astar = dj(start,goal,gridded_map) # Replace None, None with a call to your implementation of A*
        nodes_expanded_astar.append(expanded_astar)

        if cost != solution_costs[i]:
            print("There is a mismatch in the solution cost found by A* and what was expected for the problem:")
            print("Start state: ", start)
            print("Goal state: ", goal)
            print("Solution cost encountered: ", cost)
            print("Solution cost expected: ", solution_costs[i])
            print()

        cost, expanded_mm = None, None # Replace None, None with a call to your implementation of MM
        nodes_expanded_mm.append(expanded_mm)
        
        if cost != solution_costs[i]:
            print("There is a mismatch in the solution cost found by MM and what was expected for the problem:")
            print("Start state: ", start)
            print("Goal state: ", goal)
            print("Solution cost encountered: ", cost)
            print("Solution cost expected: ", solution_costs[i])
            print()

        cost, expanded_biastar = bibs(start,goal,gridded_map)# Replace None, None with a call to your implementation of Bi-A*
        nodes_expanded_biastar.append(expanded_biastar)
        
        if cost != solution_costs[i]:
            print("There is a mismatch in the solution cost found by Bi-A* and what was expected for the problem:")
            print("Start state: ", start)
            print("Goal state: ", goal)
            print("Solution cost encountered: ", cost)
            print("Solution cost expected: ", solution_costs[i])
            print()
    
    print('Finished running all tests. The implementation of an algorithm is likely correct if you do not see mismatch messages for it.')

    if plots:
        from search.plot_results import PlotResults
        plotter = PlotResults()
        plotter.plot_results(nodes_expanded_mm, nodes_expanded_astar, "Nodes Expanded (MM)", "Nodes Expanded (A*)", "nodes_expanded_mm_astar")
        plotter.plot_results(nodes_expanded_mm, nodes_expanded_biastar, "Nodes Expanded (MM)", "Nodes Expanded (Bi-A*)", "nodes_expanded_mm_biastar")


def dj(s,g,gridded_map):
    #gridded_map = Map("dao-map/brc000d.map")
    OPEN = []
    CLOSED = {}
    heapq.heappush(OPEN,s)
    NodesExpanded=0
    #print(s._g)
    dx = abs(s._x-g._x)
    dy = abs(s._y-g._y)
    hs = 1.5*(min(dx,dy)) + abs(dx-dy)
    s._cost = hs
    CLOSED[s.state_hash()] = s
    while (len(OPEN) != 0):
        n = OPEN.pop(0)

        NodesExpanded+=1
        if n == g:

            return n._cost,NodesExpanded
        children = gridded_map.successors(n)
        nx,ny = n._x, n._y

        for x in children:
            hash = x.state_hash()
            cx,cy = x._x,x._y
            dx = abs(cx-g._x)
            dy = abs(cy-g._y)
            hs = 1.5*(min(dx,dy)) + abs(dx-dy) 
            x._cost +=x._g 
                
            if x.state_hash() not in CLOSED:
                heapq.heappush(OPEN,x)
                CLOSED[hash] = x
                
            if hash in CLOSED and x._g< CLOSED[hash]._g:
                heapq.heappush(OPEN,x)
                CLOSED[hash] = x
   
    
    #gridded_map.plot_map(CLOSED,s,g,'ponisDL')
    return float(-1),None    

def bibs(s,g,gridded_map):
    #gridded_map = Map("dao-map/brc000d.map")
    ###########################FORWARD
    OPENf = []
    CLOSEDf = {}
    heapq.heappush(OPENf,s)
    dx = abs(s._x-g._x)
    dy = abs(s._y-g._y)
    hs = 1.5*(min(dx,dy)) + abs(dx-dy)
    s._cost = hs

    dx = abs(g._x-g._x)
    dy = abs(g._y-g._y)
    hs = 1.5*(min(dx,dy)) + abs(dx-dy)
    g._cost = hs
    CLOSEDf[s.state_hash()] = s
    ##########################BACKWARD
    OPENb = []
    CLOSEDb = {}
    heapq.heappush(OPENb,g)
    CLOSEDb[g.state_hash()] = g
    NodesExpanded=0
    u=float('inf')
    while (len(OPENf) != 0) and (len(OPENb) !=0):
        #stopping condition
        if u<= (OPENf[0]._cost + OPENb[0]._cost):
            #print(NodesExpanded)
            #gridded_map.plot_map(CLOSEDb|CLOSEDf,s,g,'ponisbibs')
            return u,NodesExpanded
        
        #expanding forward search
        if OPENf[0]._g<OPENb[0]._g:

            n =OPENf.pop(0)
            NodesExpanded+=1
            children = gridded_map.successors(n)
            nx,ny = n._x, n._y
            for x in children:
                
                hash = x.state_hash()
                cx,cy = x._x,x._y
                dx = abs(cx-g._x)
                dy = abs(cy-g._y)
                hs = 1.5*(min(dx,dy)) + abs(dx-dy) 
                x._cost +=x._g 
                #Found a solution path going through x
                if hash in  CLOSEDb:
                    u = min(u, x._cost +CLOSEDb[hash]._cost)
                if hash not in CLOSEDf:
                    heapq.heappush(OPENf,x)
                    CLOSEDf[hash] = x

                #If it has found a better path
                if hash in CLOSEDf and x._g< CLOSEDf[hash]._g:
                    heapq.heappush(OPENf,x)
                    CLOSEDf[hash]._g = x._g
                    heapq.heapify(OPENf)
        else:
            #Same as abpove but with OPENb
            n =OPENb.pop(0)
            NodesExpanded+=1
            children = gridded_map.successors(n)
            nx,ny = n._x, n._y
            for x in children:
                hash = x.state_hash()
                cx,cy = x._x,x._y
                dx = abs(cx-g._x)
                dy = abs(cy-g._y)
                hs = 1.5*(min(dx,dy)) + abs(dx-dy) 
                x._cost +=x._g 
                if hash in  CLOSEDf:
                    u = min(u, x._cost +CLOSEDf[hash]._cost)
                if hash not in CLOSEDb:
                    heapq.heappush(OPENb,x)
                    CLOSEDb[hash] = x
                if hash in CLOSEDb and x._g< CLOSEDb[hash]._g:
                    heapq.heappush(OPENb,x)
                    CLOSEDb[hash]._g = x._g
                    heapq.heapify(OPENb)
    #gridded_map.plot_map(CLOSEDb|CLOSEDf,s,g,'ponisbibsL')
    return -1,None
        

if __name__ == "__main__":
    main()