module Recluse

using RecipesBase

export Hammock, savegwl, Window, Box, Cantilever

abstract type GWLObject end #supertype for things representable by gwl code

#helper function to rotate vectors around the z axis
function zrotate(v,rotation)
    rotmatrix=[cos(rotation) -sin(rotation)
               sin(rotation) cos(rotation)]
    return rotmatrix*v
end

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
        #do the rotation and translation, then add on our z coordinate
        transformedpointsvec=map(pointsvec) do p        
            vcat(center + zrotate(p,rotation),z)
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
    legend --> false
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
savegwl(filename,objects...)
```
Save a .gwl file
"""
function savegwl(filename::String,h::GWLObject...)
    open(filename,"w") do io
        for hi in h
            print(io,hi)
        end
    end
end

function savegwl(filename::String,h::Hammock;bufsize=200)
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
    legend --> false
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

"""
```julia
Box(length, width, height, dslice, dhatch; [rotation,center,chamfer])
```
Build a 3D box with the provided dimensions. `dslice` and `dhatch`
are the maximum allowable slicing and hatching distances (the number
of required layers and lines will be rounded up). When `rotation` is
zero, the dimension corresponding to `length` will be along the x axis.
`chamfer` should be a matrix with size (2,2) specifying the angle with
which the edges of the box should be tapered. chamfer[1,:] specifies
tapering normal to the 'length' dimension, chamfer[2,:] specifies
tapering normal to `width`. If `chamfer` is provided the box has a
crossection of ``length × width`` at the z coordinate corresponding to
the center of the object.
"""
struct Box <: GWLObject
    hammocks::Vector{Hammock} #a box is just a bunch of hammocks
    function Box(length::Number, width::Number, height::Number,
                 dslice::Number, dhatch::Number;
                 rotation = 0, center=[0,0,0],
                 chamfer=zeros(Float64,2,2))
        @assert size(center) == (3,)
        @assert size(chamfer) == (2,2)
        #get the z position of all of our layers
        numz = ceil(Int,height/dslice) #minimum value for dslice
        zpos=range(center[3] - height/2,
                   center[3] + height/2, length = numz) |> collect
        #build all of the hammocks
        hammocks = map(1:Base.length(zpos)) do i
            #apply the chamfer
            zoffset=zpos[i] - center[3]
            #amount of material added to each edge
            #negative sign so positive chamfers correspond to shapes
            #without overhang
            added = -1*(zoffset * tan.(chamfer))
            #the change in the center of the crosssection
            deltacenter = (added[:,2] - added[:,1]) / 2
            #we need to apply the eventual rotation to deltacenter
            deltacenter = zrotate(deltacenter,rotation)
            thiscenter=center[1:2] + deltacenter
            (lprime,wprime) = [length,width] + sum(added,dims=2)
            #we have to swap the length and width if i is odd
            rot = isodd(i) ? rotation : (rotation + pi/2)
            l = isodd(i) ? lprime : wprime
            w = isodd(i) ? wprime : lprime
            #calculate the number of lines
            nt=ceil(Int,l/dhatch)
            Hammock(l,w,zpos[i],nt,rotation=rot,center=thiscenter)
        end
        new(hammocks)
    end
end

function Base.print(io::IO,b::Box)
    print.(io,b.hammocks)
    return nothing
end

@recipe function plotbox(b::Box)
    legend --> false
    aspect_ratio --> 1
    for h in b.hammocks
        @series begin
            (h.points[:,1], h.points[:,2], h.points[:,3])
        end
    end
end


"""
```julia
Cantilever(length,width,height,
           segmentlength,dslice,dhatch,
           overlap, overlapangle;
           [rotation, center, chamfer])
```
Write a segmented cantilever beam which is assumed to be supported on the face
normal to the `length` direction with the lowest coordinate value. All arguments
shared with `Box` have the same interpretation. `segmentlength` is the maximum
length of a single segment, `overlap` is the amount of overlap between neighboring
segments, `overlapangle` is the angle that the segments are split with.
"""
struct Cantilever <: GWLObject
    boxes::Vector{Box}
    function Cantilever(length::Number, width::Number, height::Number,
                        segmentlength::Number,dslice::Number,
                        dhatch::Number,overlap::Number,overlapangle::Number;
                        rotation = 0, center=[0,0,0],
                        chamfer=zeros(Float64,2,2))

        nsegments = ceil(Int,length/(segmentlength-(overlap/2)))
        #if nsegments is one, this could just be a box
        if nsegments < 2
            @warn "creating a Cantilever with one segment, consider replacing with a Box"
            return new([Box(length,width,height,dslice,dhatch;
                            rotation, center, chamfer)])
        end
        #this is the chamfer for the middle segments
        internalchamfer = hcat([-overlapangle, overlapangle],chamfer[2,:]) |> permutedims
        #calculate the length of one segment
        lseg = (length + (nsegments - 1)*overlap) / nsegments
        #get all of our center x coordinates assuming we are centered
        #at the origin
        centersx = range((lseg-length)/2,
                         (length-lseg)/2,length = nsegments) |> collect
        #rotate these points
        rotcenters = zrotate.(collect([cx, 0] for cx in centersx),
                              rotation)
        #translate these points
        centers = map(rotcenters) do rc
            vcat(rc,[0]) + center
        end
        boxes = map(1:Base.length(centers)) do i
            if i == 1
                thischamfer=hcat(chamfer[:,1],internalchamfer[:,2])
            elseif i == nsegments
                thischamfer=hcat(internalchamfer[:,1],chamfer[:,2])
            else
                thischamfer = internalchamfer
            end
            Box(lseg,width,height,dslice,dhatch;
                rotation,center=centers[i],chamfer=thischamfer)
        end
        return new(boxes)
    end
end

function Base.print(io::IO,c::Cantilever)
    print.(io,c.boxes)
    return nothing
end

@recipe function plotbox(c::Cantilever)
    aspect_ratio --> 1
    legend --> false
    for b in c.boxes
        @series begin
            b
        end
    end
end

end # module Recluse
