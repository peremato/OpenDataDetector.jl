using Accessors


"""
    Vector3D with doubles
"""
struct Vector3d
    x::Float64
    y::Float64
    z::Float64
    Vector3d(x=0,y=0,z=0) = new(x,y,z)
end
"""
    Vector3D with floats
"""
struct Vector3f
    x::Float32
    y::Float32
    z::Float32
    Vector3f(x=0,y=0,z=0) = new(x,y,z)
end
"""
    Vector2D with Int32
"""
struct Vector2i
    a::Int32
    b::Int32
    Vector2i(a=0,b=0) = new(a,b)
end

"""
    Description: "The Monte Carlo particle - based on the lcio::MCParticle."
    Author: "F.Gaede, DESY"
"""

#---MCParticleIdx--------------------------------------------------------------
struct MCParticleIdx <: Integer
    idx::Int32
end
Base.iszero(x::MCParticleIdx) = x.idx == 0
Base.show(io::IO, x::MCParticleIdx) = print(io, "MCParticles$(x.idx)")

#---For implementation of OneToManyRelation
struct MCParticles
    idxs::Vector{MCParticleIdx}
    MCParticles() = new(MCParticleIdx[])
end
Base.show(io::IO, c::MCParticles) = print(io, "MCParticle#$([Int64(p.idx) for p in c.idxs])")

struct MCParticle
    idx::MCParticleIdx
    #  Members
    PDG::Int32                         # PDG code of the particle
    generatorStatus::Int32             # status of the particle as defined by the generator
    simulatorStatus::Int32             # status of the particle from the simulation program - use BIT constants below
    charge::Float32                    # particle charge
    time::Float32                      # creation time of the particle in [ns] wrt. the event, e.g. for preassigned decays or decays in flight from the simulator.
    mass::Float64                      # mass of the particle in [GeV]
    vertex::Vector3d                   # production vertex of the particle in [mm].
    endpoint::Vector3d                 # endpoint of the particle in [mm]
    momentum::Vector3f                 # particle 3-momentum at the production vertex in [GeV]
    momentumAtEndpoint::Vector3f       # particle 3-momentum at the endpoint in [GeV]
    spin::Vector3f                     # spin (helicity) vector of the particle.
    colorFlow::Vector2i                # color flow as defined by the generator
    # OneToManyRelations
    parents::MCParticles
    daughters::MCParticles
    #edm4hep::MCParticle parents #  The parents of this particle.
    #edm4hep::MCParticle daughters #  The daughters this particle.
end

function MCParticle(;pdg=0, generatorStatus=0, simulatorStatus=0, charge=0, time=0, mass=0,
                    vertex=Vector3d(), endpoint=Vector3d(), momentum=Vector3f(), momentumAtEndpoint=Vector3f(),
                    spin=Vector3f(), colorFlow=Vector2i(), parents=MCParticles(), daughters=MCParticles())
    MCParticle(0, pdg,generatorStatus, simulatorStatus, charge, time, mass, vertex, endpoint, momentum, momentumAtEndpoint, spin, colorFlow, 
               parents, daughters)
end

Base.convert(::Type{MCParticle}, p::MCParticleIdx) = iszero(p.idx) ? nothing : @inbounds mcparticles[p.idx]
Base.convert(::Type{MCParticleIdx}, p::MCParticle) = iszero(p.idx) ? register(p) : return p.idx

Base.push!(c::MCParticles, p::MCParticle) = push!(c.idxs, p)
Base.iterate(c::MCParticles, i=1) = i > length(c.idxs) ? nothing : (convert(MCParticle, c.idxs[i]), i + 1)
Base.getindex(c::MCParticles, i::Int64) = convert(MCParticle, c.idxs[i])
Base.size(c::MCParticles) = size(c.idxs)
Base.length(c::MCParticles) = length(c.idxs)
Base.eltype(::Type{MCParticles}) = MCParticle


const mcparticles = MCParticle[]
const mcparticles_idxs = MCParticleIdx[]

function register(p::MCParticle)
    !iszero(p.idx) && error("Registering an already registered MCParticle $p")
    len = length(mcparticles)
    p = @set p.idx = MCParticleIdx(len + 1)
    push!(mcparticles, p)
    return p.idx
end


function Base.getproperty(obj::MCParticle, sym::Symbol)
    if sym === :parent
        return iszero(obj.parent_idx) ? nothing : convert(MCParticle, obj.parent_idx)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end



p1 = MCParticle(pdg=10)
p2 = MCParticle(pdg=11)
push!(p2.parents, p1)
push!(p1.daughters, p2)
push!(p2.daughters, MCParticle(pdg=99))
push!(p2.daughters, MCParticle(pdg=100))