# This file contains methods to generate a data set of instances (i.e., sudoku grids)
include("io.jl")

"""
Generate an n*n grid with a given density

Argument
- n: size of the grid
- density: percentage in [0, 1] of initial values in the grid
"""
function generateInstance(n::Int64, density::Float64)
	
	grid = Array{Int64}(ones(n,n))
	grid = grid .* 5

	for i in 1:n
		for j in 1:n 
			rn = rand()
			if rn < density 
				grid[i,j] = round(Int,rand()*2)+1
			end
		end
	end

	datafolder = "../data/"
	filename = "instanceTest"
	i=1
	while filename * string(i) * ".txt" in readdir(datafolder) 
		i = i+1
	end
	filename = datafolder*filename * string(i) * ".txt"
	fout = open(filename, "w")
	for i in 1:n
		for j in 1:n
			print(fout,grid[i,j])
			if j!=n 
				print(fout,",")
			end
		end
		print(fout,"\n")
	end
	close(fout)
end 

"""
Generate all the instances

Remark: a grid is generated only if the corresponding output file does not already exist
"""
function generateDataSet()

	for i in 5:7
   		for j in 1:5
   			generateInstance(i,max(min(rand(),0.2),0.1))
   		end
   	end 
    
end



