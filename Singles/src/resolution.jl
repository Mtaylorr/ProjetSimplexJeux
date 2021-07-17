# This file contains methods to solve an instance (heuristically or with CPLEX)
using CPLEX
using JuMP
using Random

include("generation.jl")

TOL = 0.00001

"""
Solve an instance with CPLEX
"""
function cplexSolve(grid::Array{Int64,2})

    sz = size(grid,1)
    possibilites=[(1,1),(1,2)]
    start = time()
    #println(possibilites)
    for (posi , posj) in possibilites

        #println(posi," ",posj)
        # Create the model
        m = Model(with_optimizer(CPLEX.Optimizer))
        @variable(m,x[1:sz,1:sz],Bin)
        @variable(m,flowH[1:sz,1:(sz-1),1:2]>=0,Int)
        @variable(m,flowV[1:(sz-1),1:sz,1:2]>=0,Int)
        @variable(m,flowSink[1:(sz),1:(sz)], Bin)
        @variable(m,phi>=0,Int)
        for i in 1:sz
            for j in 1:sz-1
                @constraint(m,x[i,j]+x[i,j+1] >= 1)
            end
        end
        for i in 1:sz
            for j in 1:sz-1
                @constraint(m,x[j,i]+x[j+1,i] >= 1)
            end
        end
        for i in 1:sz
            for j in 1:sz
                @constraint(m,sum(x[i,k] for k in 1:sz if grid[i,k]==j)<=1)
                @constraint(m,sum(x[k,i] for k in 1:sz if grid[k,i]==j)<=1)
            end
        end
        for i in 1:(sz)
            for j in 1:(sz)
                if(i==posi && j==posj)
                    continue    
                end
                if(i==1 && j==1)
                    @constraint(m,flowH[i,j,2]+flowV[i,j,2]-flowH[i,j,1]-flowV[i,j,1]-flowSink[i,j]==0)
                elseif (i==1 && j==sz)
                    @constraint(m,flowH[1,sz-1,1]+flowV[1,sz,2]-flowH[1,sz-1,2]-flowV[1,sz,1]-flowSink[i,j]==0)
                elseif(i==1)
                    @constraint(m,flowH[1,j,2]+flowH[1,j-1,1]+flowV[1,j,2]-flowH[1,j,1]-flowH[1,j-1,2]-flowV[1,j,1]-flowSink[i,j]==0)
                elseif(i==sz && j==1)
                    @constraint(m,flowH[i,j,2]+flowV[sz-1,j,1]-flowH[i,j,1]-flowV[sz-1,j,2]-flowSink[i,j]==0)
                elseif (i==sz && j==sz)
                    @constraint(m,flowH[sz,sz-1,1]+flowV[sz-1,sz,1]-flowH[sz,sz-1,2]-flowV[sz-1,sz,2]-flowSink[i,j]==0)
                elseif(i==sz)
                    @constraint(m,flowH[i,j,2]+flowH[i,j-1,1]+flowV[i-1,j,1]-(flowH[i,j,1]+flowH[i,j-1,2]+flowV[i-1,j,2])-flowSink[i,j]==0)
                elseif(j==1)
                    @constraint(m,flowV[i-1,j,1]+flowV[i,j,2]+flowH[i,j,2]-(flowV[i-1,j,2]+flowV[i,j,1]+flowH[i,j,1])-flowSink[i,j]==0)
                elseif(j==sz)
                    @constraint(m,flowV[i-1,j,1]+flowV[i,j,2]+flowH[i,j-1,1]-(flowV[i-1,j,2]+flowV[i,j,1]+flowH[i,j-1,2])-flowSink[i,j]==0)
                else
                    @constraint(m,flowV[i-1,j,1]+flowV[i,j,2]+flowH[i,j-1,1]+flowH[i,j,2]-(flowV[i-1,j,2]+flowV[i,j,1]+flowH[i,j-1,2]+flowH[i,j,1])-flowSink[i,j]==0)
                end
            end
        end 
        for i in 1:(sz)
            for j in 1:(sz-1)
                for k in 1:2
                    @constraint(m,flowH[i,j,k] <= 10000*x[i,j+1])
                    @constraint(m,flowH[i,j,k] <= 10000*x[i,j])
                    @constraint(m,flowV[j,i,k] <= 10000*x[j,i])
                    @constraint(m,flowV[j,i,k] <= 10000*x[j+1,i])
                end
            end
        end

        i=posi
        j=posj
        if(i==1 && j==1)
            @constraint(m,phi-flowH[1,1,1]-flowV[1,1,1]-flowSink[i,j]==0)
        elseif (i==1 && j==sz)
            @constraint(m,phi-flowH[1,sz-1,2]-flowV[1,sz,1]-flowSink[i,j]==0)
        elseif(i==1)
            @constraint(m,phi-flowH[1,j,1]-flowH[1,j-1,2]-flowV[1,j,1]-flowSink[i,j]==0)
        elseif(i==sz && j==1)
            @constraint(m,phi-flowH[i,j,1]-flowV[sz-1,j,2]-flowSink[i,j]==0)
        elseif (i==sz && j==sz)
            @constraint(m,phi-flowH[sz,sz-1,2]-flowV[sz-1,sz,2]-flowSink[i,j]==0)
        elseif(i==sz)
            @constraint(m,phi-(flowH[i,j,1]+flowH[i,j-1,2]+flowV[i-1,j,2])-flowSink[i,j]==0)
        elseif(j==1)
            @constraint(m,phi-(flowV[i-1,j,2]+flowV[i,j,1]+flowH[i,j,1])-flowSink[i,j]==0)
        elseif(j==sz)
            @constraint(m,phi-(flowV[i-1,j,2]+flowV[i,j,1]+flowH[i,j-1,2])-flowSink[i,j]==0)
        else
            @constraint(m,phi-(flowV[i-1,j,2]+flowV[i,j,1]+flowH[i,j-1,2]+flowH[i,j,1])-flowSink[i,j]==0)
        end

        @constraint(m,phi-sum(x[i,j] for i in 1:(sz) for j in 1:(sz))==0)
        @constraint(m,phi-sum(flowSink[i,j] for i in 1:sz for j in 1:sz) == 0)

        for i in 1:sz
            for j in 1:sz
                @constraint(m,flowSink[i,j] == x[i,j])
            end
        end
        @objective(m,Max,phi)
        # Start a chronometer
        

        # Solve the model
        optimize!(m)

        # Return:
        # 1 - true if an optimum is found
        # 2 - the resolution time
        if JuMP.primal_status(m) == JuMP.MathOptInterface.FEASIBLE_POINT
            return JuMP.primal_status(m) == JuMP.MathOptInterface.FEASIBLE_POINT, time() - start,JuMP.value.(x)
        end
    end 
    return false, time() - start, []
end

"""
Heuristically solve an instance
"""

function notValid(mask::Array{Int64,2},grid::Array{Int64,2},x,y)
    sz = size(mask,1)
    cntCols = 0
    cntLines = 0
    if mask[x,y] == 0
        for k in 1:y
            if grid[x,y]==grid[x,k] && mask[x,k]==0
                cntCols +=1
            end
        end

        for k in 1:x
            if grid[x,y] == grid[k,y] && mask[k,y]==0
                cntLines +=1
            end
        end
        #println(cntCols," ",cntLines)
        if cntCols > 1 || cntLines > 1
            return true
        end
    end
    for i in 1:y-1
        if mask[x,i]==1 && mask[x,i+1]==1
            return true
        end
    end
    for i in 1:x-1
        if mask[i,y]==1 && mask[i+1,y]==1
            return true
        end
    end
    return false
end

function merge(id::Array{Tuple{Int64,Int64},2},x,y,i,j)
    id1 = getId(x,y,id)
    id[id1[1],id1[2]] = getId(i,j,id)
end

function getId(x,y,id::Array{Tuple{Int64,Int64},2})
    if id[x,y] == (x,y)
        return (x,y)
    end
    id[x,y] = getId(id[x,y][1],id[x,y][2],id)
    return id[x,y]
end

function nbOfComponents(mask::Array{Int64,2},grid::Array{Int64,2})
    sz = size(mask,1)
    id = Array{Tuple{Int64,Int64},2}(undef,sz,sz)
    for i in 1:sz
        for j in 1:sz
            id[i,j] = (i,j)
        end
    end
    for i in 1:sz
        for j in 1:sz
            if i>1 && mask[i,j]==0 && mask[i-1,j]==0
                merge(id,i,j,i-1,j)
            end
            if j>1 && mask[i,j] == 0 && mask[i,j-1]==0
                merge(id,i,j,i,j-1)
            end
        end
    end

    nbCompt = 0
    for i in 1:sz
        for j in 1:sz
            if getId(i,j,id) == (i,j) && mask[i,j]==0
                nbCompt+=1
            end
        end
    end
    #println(nbCompt)
    return nbCompt
end

function valid(mask::Array{Int64,2},grid::Array{Int64,2})
    nb = nbOfComponents(mask,grid)
    #println(nb)
    return nb==1
end

function solve(x::Int64,y::Int64,mask::Array{Int64,2},hidden::Array{Int64,2},grid::Array{Int64,2},start)
    sz = size(mask,1)
    if time()-start > 100
        return false
    end
    if nbOfComponents(mask,grid)>1
        return false
    end
    if x==sz && y==sz
        if valid(mask,grid)
            return true
        end
        return false
    end
    if hidden[x,y] == 0
        if y < sz
            return solve(x,y+1,mask,hidden,grid,start)
        else
            return solve(x+1,1,mask,hidden,grid,start)
        end
    else
        if y < sz
            if !notValid(mask,grid,x,y) && solve(x,y+1,mask,hidden,grid,start)
                return true
            end
        else
            if !notValid(mask,grid,x,y) && solve(x+1,1,mask,hidden,grid,start)
                return true
            end
        end
        mask[x,y] = 1
        if y < sz
            if !notValid(mask,grid,x,y) && solve(x,y+1,mask,hidden,grid,start)
                return true
            end
        else
            if !notValid(mask,grid,x,y) && solve(x+1,1,mask,hidden,grid,start)
                return true
            end
        end
        mask[x,y] = 0
    end
    return false
end

function heuristicSolve(grid::Array{Int64,2})
    sz = size(grid,1)
    mask = zeros(Int64,sz,sz)
    hidden = zeros(Int64,sz,sz)

    for i in 1:sz
        for j in 1:sz
            cntCols = 0
            cntLines = 0
            for k in 1:sz
                if grid[i,j]==grid[i,k]
                    cntCols +=1
                end
                if grid[i,j] == grid[k,j]
                    cntLines +=1
                end
            end
            if cntCols > 1 || cntLines > 1
                hidden[i,j] = 1
            end
        end
    end
    start = time()
    if solve(1,1,mask,hidden,grid,start)
        return true, time() - start , mask
    end
    return false, time() - start, []
end
"""
Solve all the instances contained in "../data" through CPLEX and heuristics

The results are written in "../res/cplex" and "../res/heuristic"

Remark: If an instance has previously been solved (either by cplex or the heuristic) it will not be solved again
"""
function solveDataSet()

    dataFolder = "../data/"
    resFolder = "../res/"

    # Array which contains the name of the resolution methods
    resolutionMethod = ["cplex"]
    #resolutionMethod = ["heuristique"]
    #resolutionMethod = ["cplex", "heuristique"]

    # Array which contains the result folder of each resolution method
    resolutionFolder = resFolder .* resolutionMethod

    # Create each result folder if it does not exist
    for folder in resolutionFolder
        if !isdir(folder)
            mkdir(folder)
        end
    end
            
    global isOptimal = false
    global solveTime = -1

    # For each instance
    # (for each file in folder dataFolder which ends by ".txt")
    for file in filter(x->occursin(".txt", x), readdir(dataFolder))  

        println("-- Resolution of ", file)
        grid = readInputFile(dataFolder * file)

        # TODO
        println("In file resolution.jl, in method solveDataSet(), TODO: read value returned by readInputFile()")
        
        # For each resolution method
        for methodId in 1:size(resolutionMethod, 1)
            
            outputFile = resolutionFolder[methodId] * "/" * file

            # If the instance has not already been solved by this method
            if !isfile(outputFile)
                
                fout = open(outputFile, "w")  

                resolutionTime = -1
                isOptimal = false
                
                # If the method is cplex
                if resolutionMethod[methodId] == "cplex"
                    
                    # TODO 
                    #println("In file resolution.jl, in method solveDataSet(), TODO: fix cplexSolve() arguments and returned values")
                    
                    # Solve it and get the results
                    isOptimal, resolutionTime, X = cplexSolve(grid)
                    #println(X)
                    # If a solution is found, write it
                    if isOptimal
                        solveTime = resolutionTime
                        for i in 1:size(X,1)
                            for j in 1:size(X,2)
                                if(X[i,j]>=0.5)
                                    print(fout,grid[i,j])
                                else
                                    print(fout,0)
                                end
                                if(j!=size(X,2))
                                    print(fout,",")
                                end
                            end
                            println(fout,"")
                        end
                    else
                        rm(dataFolder*file)
                    end
                # If the method is one of the heuristics
                else
                        isOptimal,solveTime, mask = heuristicSolve(grid)
                    if isOptimal
                        sz = size(mask,1)
                        for i in 1:sz
                            for j in 1:sz
                                if(mask[i,j]==0)
                                    print(fout,grid[i,j])
                                else
                                    print(fout,0)
                                end
                                if(j!=sz)
                                    print(fout,",")
                                end
                            end
                            println(fout,"")
                        end
                    end
                end

                println(fout, "solveTime = ", solveTime)
                println(fout, "isOptimal = ", isOptimal)
                close(fout)
            end


            # Display the results obtained with the method on the current instance
            include(outputFile)
            displaySolution(outputFile)
            println(resolutionMethod[methodId], " optimal: ", isOptimal)
            println(resolutionMethod[methodId], " time: " * string(round(solveTime, sigdigits=2)) * "s\n")
        end         
    end 
end
