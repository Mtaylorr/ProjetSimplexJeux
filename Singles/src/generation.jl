# This file contains methods to generate a data set of instances (i.e., sudoku grids)
include("io.jl")

"""
Generate an n*n grid with a given density

Argument
- n: size of the grid
- density: percentage in [0, 1] of initial values in the grid
"""
lignes = []
colonnes = []

function getK(l::Int64,c::Int64,n::Int64)
    global lignes
    global colonnes
    k = 0
    for i in 1:n
        if lignes[l,i]==0 && colonnes[c,i]==0
            k = i
            lignes[l,i] += 1
            colonnes[c,i] += 1
            break
        end
    end
    return k
end

function generateInstance(n::Int64,density::Float64)
    startX = rand(1:n)
    startY = rand(1:n)
    nodeAdded = 0
    nodeNeeded = Int64(n*n - floor(n*n*density))
    nextNode = [(startX,startY)]
    usedNodes = [(startX,startY)]
    global lignes
    global colonnes
    lignes = zeros(Int64,n,n)
    colonnes = zeros(Int64,n,n)

    instance = zeros(Int64, n, n)

    while nodeAdded < nodeNeeded && length(nextNode) != 0
        currentPos = nextNode[rand(1:length(nextNode))]
        i = currentPos[1]
        j = currentPos[2]
        k = getK(i,j,n)

        instance[i,j] = k

        nodeAdded = nodeAdded+1
        push!(usedNodes,currentPos)
        nextNode = filter!(e->e!=currentPos,nextNode)
        possibleNodes = []
        if i == 1 && j == 1
            push!(possibleNodes,(1,2))
            push!(possibleNodes,(2,1))
        elseif i== 1 && j==n
            push!(possibleNodes,(1,j-1))
            push!(possibleNodes,(i+1,j))
        elseif i==1
            push!(possibleNodes,(i,j-1))
            push!(possibleNodes,(i,j+1))
            push!(possibleNodes,(i+1,j))
        elseif i==n && j == 1
            push!(possibleNodes,(i-1,j))
            push!(possibleNodes,(i,j+1))
        elseif i==n && j==n
            push!(possibleNodes,(i-1,j))
            push!(possibleNodes,(i,j-1))
        elseif i==n
            push!(possibleNodes,(i-1,j))
            push!(possibleNodes,(i,j+1))
            push!(possibleNodes,(i,j-1))
        elseif j == 1
            push!(possibleNodes,(i-1,j))
            push!(possibleNodes,(i+1,j))
            push!(possibleNodes,(i,j+1))
        elseif j==n
            push!(possibleNodes,(i,j-1))
            push!(possibleNodes,(i-1,j))
            push!(possibleNodes,(i+1,j))
        else
            push!(possibleNodes,(i,j+1))
            push!(possibleNodes,(i,j-1))
            push!(possibleNodes,(i-1,j))
            push!(possibleNodes,(i+1,j))
        end
        if issubset(possibleNodes,usedNodes)
            continue
        else

            nx = possibleNodes[rand(1:length(possibleNodes))]
            while nx in usedNodes
                nx = possibleNodes[rand(1:length(possibleNodes))]
            end
            push!(nextNode,nx)
            possibleNodes = filter!(e->e!=nx,possibleNodes)
            for k in 1:length(possibleNodes)
                if !(possibleNodes[k] in usedNodes) && rand()>=0.5
                    push!(usedNodes,possibleNodes[k])
                    push!(nextNode,possibleNodes[k])
                end
            end
        end
    end
    for i in 1:n
        for j in 1:n
            if instance[i,j] == 0
                instance[i,j] = rand(1:n)
            end
        end
    end
    return instance
end 

"""
Generate all the instances

Remark: a grid is generated only if the corresponding output file does not already exist
"""
function generateDataSet(n::Int64)
    dataFolder = "../data/"
    fileName = dataFolder*"instanceTest"
    for i in 1:n

            k=1
            while isfile(fileName*string(k)*".txt")
                k=k+1
            end
            outputFile = fileName *string(k)*".txt"
            fout = open(outputFile, "w")
            n = rand(5:15)
            density = 0.2+rand()*0.2
            X = generateInstance(n,density)
            for i in 1:n
                for j in 1:n
                    print(fout,X[i,j])
                    if(j!=n)
                        print(fout,",")
                    end
                end
                println(fout,"")
            end
            close(fout)

            #println("n ",n, " density ",density)
            #displayGrid(fileName*string(k)*".txt")
    end
end



