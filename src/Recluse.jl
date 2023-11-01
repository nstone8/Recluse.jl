module Recluse

using RecipesBase

export Hammock, savegwl, Window

abstract type GWLObject end #supertype for things representable by gwl code

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
struct Hammock <: GWLObject
    points::Matrix{Float64} #GWL scripts only do floats

    function Hammock(length::Number,width::Number,z,numthreads::Number,rotation::Number;center=[0,0])

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

function Hammock(length::Number,width::Number,z,numthreads::Number;rotation=0,center=[0,0])
    Hammock(length,width,z,numthreads,rotation;center)
end
#add a plot recipe
@recipe function plothammock(h::Hammock)
    aspect_ratio --> 1
    (h.points[:,1], h.points[:,2])
end

function Base.print(io::IO,h::Hammock;bufsize=200)
    @assert iseven(bufsize)
    for i in 1:size(h.points)[1]
        println(io,join(h.points[i,:],"\t"))
        if i%bufsize == 0
            println(io,"write")
        end
    end
    if size(h.points)[1]%bufsize != 0
        println(io,"write")
    end
end

"""
```julia
savegwl(filename,element)
```
Save a .gwl file
"""
function savegwl(filename::String,h::GWLObject)
    open(filename,"w") do io
        print(io,h)
    end
end

function savegwl(filename::String,h::Hammock;bufsize)
    open(filename,"w") do io
        print(io,h;bufsize)
    end
end

"""
```julia
Window(width, height, z0, numthreads; [rotation, center, taper])
```

Cover a rectangle with size `width` in dimension 1 and size `height` in dimension 3
with `numthreads` line segments along dimension 2. These line segments will be written
on the xz plane such that the bottom edge of the rectangle has coordinate z0 in dimension 3.

If the `taper` argument is used, `width` refers to the width of the window along its bottom
edge.

# Optional keyword arguments
- `rotation`: if passed a scalar, rotate the resulting rectangle about the z axis. If provided
  a (length 2) vector, perform a rotation about the center of the bottom edge, first around the y
  axis, then around z.

- `center`: translate the center of the rectangle's bottom edge parallel to the xy plane,
  default is `[0,0]`

- `taper`: Create a trapezoid rather than a rectangle
"""
struct Window <: GWLObject
    points::Matrix{Float64} #GWL scripts only do floats

    function Window(width::Number,height::Number,z0,numthreads::Number,rotation::Vector{<:Number},taper::Number;center=[0,0])
        @assert length(rotation) == 2 "rotation vector must have length 2"
        #first generate points with no rotation centered on (0,0)
        xpoints=(width/2) * vcat(map(1:numthreads) do i
                                     if isodd(i)
                                         return [-1,1]
                                     else
                                         return [1,-1]
                                     end
                                 end...)
        spacing = height / (numthreads-1)
        zpoints = (-height/2) .+ (spacing * vcat((repeat([j],2) for j in 0:(numthreads-1))...))
        pointsvec=map(zip(xpoints,zpoints)) do (x,z) #get a vector of [x,y,z] coordinate vectors, setting z to zero for now
            [x,0,z]
        end

        #do the tapering
        pointsvec=map(pointsvec) do p
            heightfrombottom=p[3] + height/2
            return [p[1] + sign(p[1])*heightfrombottom*tan(taper), p[2], p[3]]
        end

        #build a rotation matrix around z axis
        zrotmatrix=[cos(rotation[1]) -sin(rotation[1]) 0
                    sin(rotation[1]) cos(rotation[1])  0
                    0             0              1]
        #rotation around x
        xrotmatrix=[1 0                0
                    0 cos(rotation[2]) -sin(rotation[2])
                    0 sin(rotation[2]) cos(rotation[2])]
        
        #=rotation around y
        yrotmatrix=[cos(rotation[2])  0 sin(rotation[2])
                    0                 1 0
                    -sin(rotation[2]) 0 cos(rotation[2])]
        =#
        #do the rotation and translation
        transformedpointsvec=map(pointsvec) do p
            #first translate the center of the bottom edge to the origin
            p += [0,0,height/2]
            #now do the rotation (about x first so lines are always parallel to the xy plane) then add on z0
            vcat(center,z0) + (zrotmatrix * xrotmatrix * p)
        end
        #collapse our points vector into a nx3 matrix
        pointsmat=vcat((permutedims(p) for p in transformedpointsvec)...)
        return new(pointsmat)
    end
end

function Window(width::Number,height::Number,z0,numthreads::Number,rotation::Number,taper;kwargs...)
    Window(width,height,z0,numthreads,[rotation,0],taper;kwargs...)
end

function Window(width::Number,height::Number,z0,numthreads::Number;rotation=0,taper=0,kwargs...)
    Window(width,height,z0,numthreads,rotation,taper;kwargs...)
end

#add a plot recipe
@recipe function plotwindow(w::Window)
    aspect_ratio --> 1
    #need to insert 'missing' between each pair of lines
    points=permutedims(w.points)
    points=reshape(points,3,2,:)
    m=Array{Missing}(missing,size(points)...)
    points=cat(points,m,dims=4)
    points=permutedims(points,[1,2,4,3])
    points=reshape(points,3,:)
    #all that had the effect of inserting two 'missing' points between each line
    (points[1,:],points[2,:],points[3,:])
end

function Base.print(io::IO,w::Window)
    r=reshape(permutedims(w.points),3,2,:)
    for li in 1:size(r)[3]
        println(io,join(r[:,1,li],"\t"))
        println(io,join(r[:,2,li],"\t"))
        println(io,"write")
    end
end
end # module Recluse
