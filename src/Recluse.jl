module Recluse

using RecipesBase

export Hammock, savegwl

"""
```julia
Hammock(length, width, z, numthreads; [rotation, center])
```

Cover a rectangle with size `length` in dimension 1 and size `width` in dimension 2
with `numthreads` line segments along dimension 2. These line segments will be written
on a plane with constant coordinate `z` in dimension 3.

# Optional keyword arguments
- `rotation`: rotate the rectangle about its center, default is zero

- `center`: translate the center of the rectangle, default is `[0,0]`
"""
struct Hammock
    points::Matrix{Float64} #GWL scripts only do floats

    function Hammock(length::Number,width::Number,z,numthreads::Number;rotation=0,center=[0,0])

        #first generate points with no rotation centered on (0,0)
        ypoints=(width/2) * vcat(map(1:numthreads) do i
                                     if isodd(i)
                                         return [-1,1]
                                     else
                                         return [1,-1]
                                     end
                                 end...)
        spacing = length / (numthreads-1)
        xpoints = (-length/2) .+ (spacing * vcat((repeat([j],2) for j in 0:(numthreads-1))...))
        pointsvec=map(zip(xpoints,ypoints)) do (x,y) #get a vector of [x,y] coordinate vectors
            [x,y]
        end
        #build a rotation matrix
        rotmatrix=[cos(rotation) -sin(rotation)
                   sin(rotation) cos(rotation)]
        #do the rotation and translation, then add on our z coordinate
        transformedpointsvec=map(pointsvec) do p        
            vcat(center + (rotmatrix * p),z)
        end
        #collapse our points vector into a nx2 matrix
        pointsmat=vcat((permutedims(p) for p in transformedpointsvec)...)
        return new(pointsmat)
    end
end

#add a plot recipe
@recipe function plothammock(h::Hammock)
    aspect_ratio --> 1
    (h.points[:,1], h.points[:,2])
end

function Base.print(io::IO,h::Hammock)
    for i in 1:size(h.points)[1]
        println(io,join(h.points[i,:],"\t"))
    end
    println(io,"write")
end

function savegwl(filename::String,h::Hammock)
    open(filename,"w") do io
        print(io,h)
    end
end

end # module Recluse
