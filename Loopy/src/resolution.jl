# This file contains methods to solve an instance (heuristically or with CPLEX)
using CPLEX
using JuMP
include("generation.jl")

TOL = 0.00001

"""
Solve an instance with CPLEX
"""
function cplexSolve(grid::Array{Int64,2})
    sz = size(grid,1)
    posi=0
    posj=0
    posval=0
    for i in 1:size(grid,1)
        for j in 1:size(grid,2)
            if grid[i,j] == 3
                posi=i
                posj=j
                posval=3
            elseif posval!=3 && grid[i,j]!=5
                posi=i
                posj=j
                posval=grid[i,j]
            end
        end
    end
    possibilites = []
    if(posval==3)
        possibilites=[(posi,posj)]
    else
        possibilites=[(posi,posj),(posi+1,posj),(posi,posj+1),(posi+1,posj+1)]
    end
    start = time()
    for (posi , posj) in possibilites 
        # Create the model
        m = Model(with_optimizer(CPLEX.Optimizer))
        @variable(m,degree[1:(sz+1),1:(sz+1)]>=0,Int)
        @variable(m,degreeCond[1:(sz+1),1:(sz+1)],Bin)
        @variable(m,edgesH[1:(sz+1),1:sz],Bin)
        @variable(m,edgesV[1:(sz),1:(sz+1)],Bin)
        @variable(m,flowH[1:(sz+1),1:sz,1:2]>=0,Int)
        @variable(m,flowV[1:(sz),1:(sz+1),1:2]>=0,Int)
        @variable(m,capH[1:(sz+1),1:sz]>=0,Int)
        @variable(m,capV[1:(sz),1:(sz+1)]>=0,Int)
        @variable(m,flowSink[1:(sz+1),1:(sz+1)], Bin)
        @variable(m,phi>=0,Int)
        for i in 1:sz
            for j in 1:sz
                if(grid[i,j]!=5)
                    @constraint(m,edgesH[i,j]+edgesH[i+1,j]+edgesV[i,j]+edgesV[i,j+1]==grid[i,j])
                else
                    @constraint(m,edgesH[i,j]+edgesH[i+1,j]+edgesV[i,j]+edgesV[i,j+1]<=3)  
                end
            end
        end
        for i in 1:(sz+1)
            for j in 1:(sz+1)
                @constraint(m,degree[i,j]-2*degreeCond[i,j]==0)
                if(i==posi && j==posj)
                    continue    
                end
                if(i==1 && j==1) 
                    @constraint(m,edgesH[i,j]+edgesV[i,j]-degree[i,j]==0)
                    @constraint(m,flowH[i,j,2]+flowV[i,j,2]-flowH[i,j,1]-flowV[i,j,1]-flowSink[i,j]==0)
                elseif (i==1 && j==sz+1)
                    @constraint(m,edgesH[i,j-1]+edgesV[i,j]-degree[i,j]==0)
                    @constraint(m,flowH[i,j-1,1]+flowV[i,j,2]-flowH[i,j-1,2]-flowV[i,j,1]-flowSink[i,j]==0)
                elseif(i==1)
                    @constraint(m,edgesH[i,j]+edgesH[i,j-1]+edgesV[i,j]-degree[i,j]==0)
                    @constraint(m,flowH[i,j,2]+flowH[i,j-1,1]+flowV[i,j,2]-flowH[i,j,1]-flowH[i,j-1,2]-flowV[i,j,1]-flowSink[i,j]==0)
                elseif(i==sz+1 && j==1) 
                    @constraint(m,edgesH[i,j]+edgesV[i-1,j]-degree[i,j]==0)
                    @constraint(m,flowH[i,j,2]+flowV[i-1,j,1]-flowH[i,j,1]-flowV[i-1,j,2]-flowSink[i,j]==0)
                elseif (i==sz+1 && j==sz+1)
                    @constraint(m,edgesH[i,j-1]+edgesV[i-1,j]-degree[i,j]==0)
                    @constraint(m,flowH[i,j-1,1]+flowV[i-1,j,1]-flowH[i,j-1,2]-flowV[i-1,j,2]-flowSink[i,j]==0)
                elseif(i==sz+1)
                    @constraint(m,edgesH[i,j]+edgesH[i,j-1]+edgesV[i-1,j]-degree[i,j]==0)
                    @constraint(m,flowH[i,j,2]+flowH[i,j-1,1]+flowV[i-1,j,1]-(flowH[i,j,1]+flowH[i,j-1,2]+flowV[i-1,j,2])-flowSink[i,j]==0)
                elseif(j==1)
                    @constraint(m,edgesV[i-1,j]+edgesV[i,j]+edgesH[i,j]-degree[i,j]==0)
                    @constraint(m,flowV[i-1,j,1]+flowV[i,j,2]+flowH[i,j,2]-(flowV[i-1,j,2]+flowV[i,j,1]+flowH[i,j,1])-flowSink[i,j]==0)
                elseif(j==sz+1)
                    @constraint(m,edgesV[i-1,j]+edgesV[i,j]+edgesH[i,j-1]-degree[i,j]==0)
                    @constraint(m,flowV[i-1,j,1]+flowV[i,j,2]+flowH[i,j-1,1]-(flowV[i-1,j,2]+flowV[i,j,1]+flowH[i,j-1,2])-flowSink[i,j]==0)
                else
                    @constraint(m,edgesV[i-1,j]+edgesV[i,j]+edgesH[i,j-1]+edgesH[i,j]-degree[i,j]==0)
                    @constraint(m,flowV[i-1,j,1]+flowV[i,j,2]+flowH[i,j-1,1]+flowH[i,j,2]-(flowV[i-1,j,2]+flowV[i,j,1]+flowH[i,j-1,2]+flowH[i,j,1])-flowSink[i,j]==0)
                end
            end
        end 
        i=posi
        j=posj
        if(i==1 && j==1) 
            @constraint(m,edgesH[1,1]+edgesV[1,1]-degree[i,j]==0)
            @constraint(m,phi-flowH[1,1,1]-flowV[1,1,1]-flowSink[i,j]==0)
        elseif (i==1 && j==sz+1)
            @constraint(m,edgesH[1,sz]+edgesV[1,sz+1]-degree[i,j]==0)
            @constraint(m,phi-flowH[1,sz,2]-flowV[1,sz,1]-flowSink[i,j]==0)
        elseif(i==1)
            @constraint(m,edgesH[1,j]+edgesH[1,j-1]+edgesV[1,j]-degree[i,j]==0)
            @constraint(m,phi-flowH[1,j,1]-flowH[1,j-1,2]-flowV[1,j,1]-flowSink[i,j]==0)
        elseif(i==sz+1 && j==1) 
            @constraint(m,edgesH[i,j]+edgesV[sz,j]-degree[i,j]==0)
            @constraint(m,phi-flowH[i,j,1]-flowV[sz,j,2]-flowSink[i,j]==0)
        elseif (i==sz+1 && j==sz+1)
            @constraint(m,edgesH[sz+1,sz]+edgesV[sz,sz+1]-degree[i,j]==0)
            @constraint(m,phi-flowH[sz+1,sz,2]-flowV[sz,sz+1,2]-flowSink[i,j]==0)
        elseif(i==sz+1)
            @constraint(m,edgesH[i,j]+edgesH[i,j-1]+edgesV[i-1,j]-degree[i,j]==0)
            @constraint(m,phi-(flowH[i,j,1]+flowH[i,j-1,2]+flowV[i-1,j,2])-flowSink[i,j]==0)
        elseif(j==1)
            @constraint(m,edgesV[i-1,j]+edgesV[i,j]+edgesH[i,j]-degree[i,j]==0)
            @constraint(m,phi-(flowV[i-1,j,2]+flowV[i,j,1]+flowH[i,j,1])-flowSink[i,j]==0)
        elseif(j==sz+1)
            @constraint(m,edgesV[i-1,j]+edgesV[i,j]+edgesH[i,j-1]-degree[i,j]==0)
            @constraint(m,phi-(flowV[i-1,j,2]+flowV[i,j,1]+flowH[i,j-1,2])-flowSink[i,j]==0)
        else
            @constraint(m,edgesV[i-1,j]+edgesV[i,j]+edgesH[i,j-1]+edgesH[i,j]-degree[i,j]==0)
            @constraint(m,phi-(flowV[i-1,j,2]+flowV[i,j,1]+flowH[i,j-1,2]+flowH[i,j,1])-flowSink[i,j]==0)
        end

        @constraint(m,phi-sum(flowSink[i,j] for i in 1:(sz+1) for j in 1:(sz+1))==0)
        @constraint(m,phi-sum(degreeCond[i,j] for i in 1:(sz+1) for j in 1:(sz+1))==0)
        
        for i in 1:sz+1
            for j in 1:sz 
                for k in 1:2
                    @constraint(m,flowH[i,j,k]-capH[i,j]<=0)
                end
                @constraint(m,capH[i,j]-10000*edgesH[i,j]<=0)
            end
        end
        for i in 1:sz
            for j in 1:sz+1  
                for k in 1:2
                    @constraint(m,flowV[i,j,k]-capV[i,j]<=0)
                end
                @constraint(m,capV[i,j]-10000*edgesV[i,j]<=0)
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
            return JuMP.primal_status(m) == JuMP.MathOptInterface.FEASIBLE_POINT, time() - start, JuMP.value.(edgesH),JuMP.value.(edgesV)
        end
    end 
    return false, time() - start,[],[]
end


function getid(x::Int64, y::Int64, id::Array{Tuple{Int64,Int64},2})
    if(id[x,y]==(x,y))
        return (x,y)
    else 
        id[x,y]=getid(id[x,y][1], id[x,y][2], id)
        return id[x,y]
    end
end

function uni(cell1::Tuple{Int64,Int64}, cell2::Tuple{Int64,Int64}, id::Array{Tuple{Int64,Int64},2})
    id1 = getid(cell1[1],cell1[2],id)
    id2 = getid(cell2[1],cell2[2],id)
    if(id1!=id2)
        id[id1[1],id1[2]]=id2
    end
end


function nbComponents(mask::Array{Int64,2}, grid::Array{Int64,2})
    sz = size(mask,1)
    id = Array{Tuple{Int64,Int64},2}(undef,sz+1,sz+1)
    for i in 1:sz
        for j in 1:sz
            id[i,j]=(i,j)
        end
    end
    id[sz+1,sz+1]=(sz+1,sz+1)
    nb=0
    for i in 1:sz
        for j in 1:sz
            if(i>1 && mask[i,j]==mask[i-1,j])
                uni((i,j),(i-1,j),id)
            end
            if(i<sz && mask[i,j]==mask[i+1,j])
                uni((i,j),(i+1,j),id)
            end
            if(j>1 && mask[i,j]==mask[i,j-1])
                uni((i,j),(i,j-1),id)
            end
            if(j<sz && mask[i,j]==mask[i,j+1])
                uni((i,j),(i,j+1),id)
            end
        end
    end
    for i in 1:sz
        if(mask[1,i]==1)
            uni((1,i),(sz+1,sz+1),id)
        end
        if(mask[i,1]==1)
            uni((i,1),(sz+1,sz+1),id)
        end
        if(mask[sz,i]==1)
            uni((sz,i),(sz+1,sz+1),id)
        end
        if(mask[i,sz]==1)
            uni((i,sz),(sz+1,sz+1),id)
        end
    end
    for i in 1:sz
        for j in 1:sz
            if(getid(i,j,id)==(i,j))
                nb+=1
            end
        end
    end
    if(getid(sz+1,sz+1,id)==(sz+1,sz+1))
        nb+=1
    end
    return nb
end



function valid(mask ::Array{Int64,2},grid ::Array{Int64,2})
    sz = size(mask,1)
    neighbors = Array{Int64}(zeros(sz,sz))
    for i in 1:sz
        for j in 1:sz
            if(i>1 && mask[i,j]==mask[i-1,j])
                neighbors[i,j]+=1
            end
            if(i<sz && mask[i,j]==mask[i+1,j])
                neighbors[i,j]+=1
            end
            if(j>1 && mask[i,j]==mask[i,j-1])
                neighbors[i,j]+=1
            end
            if(j<sz && mask[i,j]==mask[i,j+1])
                neighbors[i,j]+=1
            end
        end
    end
    for i in 1:sz
        if(mask[1,i]==1)
            neighbors[1,i]+=1
        end
        if(mask[i,1]==1)
            neighbors[i,1]+=1
        end
        if(mask[sz,i]==1)
            neighbors[sz,i]+=1
        end
        if(mask[i,sz]==1)
            neighbors[i,sz]+=1
        end
    end
    for i in 1:sz
        for j in 1:sz
            if(grid[i,j]<5 && neighbors[i,j]!=4-grid[i,j])
                return false
            end
        end
    end
    return nbComponents(mask,grid)==2
end


function solve(x::Int64, y::Int64, nbCells::Int64, mask ::Array{Int64,2},grid ::Array{Int64,2}, start)
    if(time()-start>60)
        return false
    end
    if(nbCells<0)
        return false
    end
    if(nbComponents(mask,grid)>2)
        return false
    end
    sz = size(mask,1)
    if( x==sz && y==sz)
        if(valid(mask,grid))
            return true
        end
        return false
    end

    if(y<sz)
        if(solve(x,y+1, nbCells, mask,grid,start))
            return true
        end
        mask[x,y]=1
        if(solve(x,y+1, nbCells-1,mask,grid,start))
            return true
        end
        mask[x,y]=0
    else
        if(solve(x+1,1, nbCells, mask,grid,start))
            return true
        end
        mask[x,y]=1
        if(solve(x+1,1, nbCells-1,mask,grid,start))
            return true
        end
        mask[x,y]=0
    end
    return false
end

function buildSolution(mask ::Array{Int64,2},edgesH ::Array{Int64,2}, edgesV)
    sz = size(mask,1)
    for i in 1:sz
        for j in 1:sz
            if(i>1 && mask[i,j]!=mask[i-1,j])
                edgesH[i,j]=1
            end
            if(i<sz && mask[i,j]!=mask[i+1,j])
                edgesH[i+1,j]=1
            end
            if(j>1 && mask[i,j]!=mask[i,j-1])
                edgesV[i,j]=1
            end
            if(j<sz && mask[i,j]!=mask[i,j+1])
                edgesV[i,j+1]=1
            end
        end
    end
    for i in 1:sz
        if(mask[1,i]==0)
            edgesH[1,i]=1
        end
        if(mask[i,1]==0)
            edgesV[i,1]=1
        end
        if(mask[sz,i]==0)
             edgesH[sz+1,i]=1
        end
        if(mask[i,sz]==0)
            edgesV[i,sz+1]=1
        end
    end

end

"""
Heuristically solve an instance
"""
function heuristicSolve(grid::Array{Int64,2})
    sz = size(grid,1)
    nbCells = sz*sz
    if(sz>=10)
        return false, 60,[],[]
    end
    if(sz>=6)
        nbCells = floor(Int64, nbCells*1/3)
    end
    start = time()
    #println(nbCells)
    #println(sz)
    mask = Array{Int64}(zeros(sz,sz))
    if solve(1,1,nbCells,mask, grid, start)
        edgesH = Array{Int64}(zeros(sz+1,sz))
        edgesV = Array{Int64}(zeros(sz,sz+1))
        buildSolution(mask,edgesH, edgesV)
        return true, time()-start, edgesH, edgesV
    end
    return false, 60,[],[]
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
    #resolutionMethod = ["cplex"]
    resolutionMethod = ["cplex", "heuristique"]

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
        sz = size(grid,1)

        # For each resolution method
        for methodId in 1:size(resolutionMethod, 1)
            
            outputFile = resolutionFolder[methodId] * "/" * file
            solveTime =-1 
            isOptimal = false
            # If the instance has not already been solved by this method
            if !isfile(outputFile)
                
                fout = open(outputFile, "w")  

                
                
                # If the method is cplex
                if resolutionMethod[methodId] == "cplex"
                    
                    # Solve it and get the results
                    isOptimal, solveTime, edgesH, edgesV = cplexSolve(grid)
                    
                    # If a solution is found, write it
                    if isOptimal
                        isOptimal = true
                        for i in 1:2*sz+1
                            if i%2==0
                                for j in 1:2*sz+1
                                    if j%2==1
                                        print(fout,round(Int64,edgesV[floor(Int,i/2),floor(Int,j/2+1)]))
                                    else
                                        print(fout,round(Int64,grid[floor(Int,i/2),floor(Int,j/2)]))
                                    end
                                    if(j!=2*sz+1)
                                        print(fout,",")
                                    end
                                end
                            else
                                for j in 1:sz
                                    print(fout,round(Int64,edgesH[floor(Int,i/2+1),j]))
                                    if(j!=sz)
                                        print(fout,",")
                                    end
                                end
                            end
                            print(fout,"\n")

                        end
                       
                    end

                # If the method is one of the heuristics
                else
                    
                    isSolved = false

                    # Solve it and get the results
                    isOptimal, solveTime, edgesH, edgesV = heuristicSolve(grid)
                    # Write the solution (if any)
                    if isOptimal
                        isOptimal = true
                        for i in 1:2*sz+1
                            if i%2==0
                                for j in 1:2*sz+1
                                    if j%2==1
                                        print(fout,round(Int64,edgesV[floor(Int,i/2),floor(Int,j/2+1)]))
                                    else
                                        print(fout,round(Int64,grid[floor(Int,i/2),floor(Int,j/2)]))
                                    end
                                    if(j!=2*sz+1)
                                        print(fout,",")
                                    end
                                end
                            else
                                for j in 1:sz
                                    print(fout,round(Int64,edgesH[floor(Int,i/2+1),j]))
                                    if(j!=sz)
                                        print(fout,",")
                                    end
                                end
                            end
                            print(fout,"\n")

                        end
                        
                    end 
                end

                println(fout, "solveTime = ", solveTime) 
                println(fout, "isOptimal = ", isOptimal)
               close(fout)
            end


            # Display the results obtained with the method on the current instance
            #include(outputFile)
            println(resolutionMethod[methodId], " optimal: ", isOptimal)
            println(resolutionMethod[methodId], " time: " * string(round(solveTime, sigdigits=2)) * "s\n")
        end         
    end 
end
