using JLD2: jldopen

export save_results, load_results, reparse_u_v1_as_u_v2

# in the migration to support non-uniform tubes, the data structures are 
# going to change significantly across the next several commits
# I'm adding more explicit versioning on file formats to keep track of the changes
# and migrate the gold standard data for regression testing to match the new formats

const RESULTS_FORMAT_VERSION = 2
const RESULTS_HISTORY_FIELDS = (:tube_hist_t, :tube_hist_u, :tube_hist_θwall, :plate_T_hist)

# Convert array wrappers to ordinary arrays and reject executable/object payloads.
_results_numeric(x::Real) = x
_results_numeric(x::AbstractArray{<:Real}) = Array(x)
_results_numeric(x::AbstractArray) = map(_results_numeric, x)
_results_numeric(x) = throw(ArgumentError("Numeric histories cannot contain $(typeof(x))"))

function _results_history(result)
    NamedTuple{RESULTS_HISTORY_FIELDS}(Tuple(
        _results_numeric(getproperty(result, field)) for field in RESULTS_HISTORY_FIELDS))
end

"""
    save_results(path, result; kind=:history)

Save a versioned JLD2 result. `:history` stores only numeric timestamps, tube
states, wall temperatures, and plate temperatures as ordinary arrays; it can
also save a history previously returned by `load_results`. `:checkpoint` stores
the complete `SimulationResult`, including integrators and their functions.
Checkpoints remain dependent on compatible Julia/package definitions; the
file version does not make serialized functions portable.

The writer emits only the current schema. Future schema changes should add a
version-specific reader rather than relabeling older payloads.
"""
function save_results(path::AbstractString, result; kind::Symbol=:history)
    kind in (:history, :checkpoint) || throw(ArgumentError("Unknown result kind: $kind"))
    kind == :checkpoint && !(result isa SimulationResult) &&
        throw(ArgumentError("A checkpoint requires a SimulationResult"))
    payload = kind == :history ? _results_history(result) : result
    jldopen(path, "w") do file
        file["format"] = "OscillatingHeatPipe.Results"
        file["format_version"] = RESULTS_FORMAT_VERSION
        file["kind"] = String(kind)
        if kind == :history
            for field in RESULTS_HISTORY_FIELDS
                file["history/$(field)"] = getproperty(payload, field)
            end
        else
            file["checkpoint"] = payload
        end
    end
    return path
end

"""
    load_results(path; kind=:auto)

Load a versioned result as a numeric-history named tuple or a complete
`SimulationResult`. `kind=:history` can read either file kind, extracting the
histories from a checkpoint. `kind=:checkpoint` requires a checkpoint file.
Unknown versions and unversioned legacy files are rejected before reading the
payload. Regenerate legacy regression fixtures with `save_results`.
"""
function load_results(path::AbstractString; kind::Symbol=:auto)
    kind in (:auto, :history, :checkpoint) || throw(ArgumentError("Unknown result kind: $kind"))
    jldopen(path, "r") do file
        haskey(file, "format") && file["format"] == "OscillatingHeatPipe.Results" ||
            throw(ArgumentError("Not a versioned OHP results file; regenerate or explicitly migrate legacy data"))
        haskey(file, "format_version") || throw(ArgumentError("Missing results format version"))
        version = file["format_version"]
        if (version == 1)
            return _load_results_v1_v2(file, kind)
        elseif (version == 2)
            return _load_results_v1_v2(file, kind)
        else
            throw(ArgumentError("Unsupported results format version: $version"))
        end
    end
end

# there are two 'kinds' of result within the file:
# 'history' contains numerical values only -> used for regression testing
# 'checkpoint' contains numerical values and integrators -> used for postprocessing
# Version 1 uses the SimulationResult struct from 205ff6bab4d94a44f9c0c082e9f7e04b0580975d
# Version 2 uses the SimulationResult struct from 0484a7d65f8b94f9d86cae09ccd94396f858a6ad
# we can reuse this function for both version because they contain the same RESULTS_HISTORY_FIELDS
# only the number/meaning of values in the state vectors are different

function _load_results_v1_v2(file, kind)
    stored_kind = file["kind"]
    stored_kind in ("history", "checkpoint") || throw(ArgumentError("Unknown stored result kind: $stored_kind"))
    kind == :checkpoint && stored_kind != "checkpoint" &&
        throw(ArgumentError("Numeric histories do not contain a checkpoint"))
    if stored_kind == "checkpoint"
        result = file["checkpoint"]
        return kind == :history ? _results_history(result) : result
    end
    NamedTuple{RESULTS_HISTORY_FIELDS}(Tuple(
        _results_numeric(file["history/$(field)"]) for field in RESULTS_HISTORY_FIELDS))
end

"""
    u_v1_to_XMδL(u) -> Vector{Tuple}, Vector{Tuple}, Vector, Vector, Vector, Vector, Vector, Vector, Vector

Disassemble the state vector `u` into vectors `Xp`, `dXdt` of the coordinates and velocities of the slug interfaces, and
mass `M` of the vapor regions, `δstart` and `δend` of film thicknesses, and `Lfilm_start`
and `Lfilm_end` of film lengths.
"""
function u_v1_to_XMδL(u::Array{Float64,1})

    maxindex = div(length(u),9)
    maxliquidindex = 0
    if mod(length(u),9) == 4
        maxliquidindex = maxindex + 1
    elseif mod(length(u),9) == 0
        maxliquidindex = maxindex
    else
        return error("The state vector given has ", length(u), " elements. We expect groupings of either 9n-5, or 9n, where n is number of liquid slugs.")
    end

    Xp =   map(tuple, zeros(maxliquidindex), zeros(maxliquidindex))
    dXdt = map(tuple, zeros(maxliquidindex), zeros(maxliquidindex))
    for i = 1:maxliquidindex
        # input Xp
        Xp[i] = (u[2*i-1],u[2*i])
        # input dXdt
        dXdt[i] = (u[2*maxliquidindex + 2*i-1],u[2*maxliquidindex + 2*i])
    end

    M =             u[4*maxliquidindex+1:           4*maxliquidindex+maxindex]
    δstart =        u[4*maxliquidindex+maxindex+1:  4*maxliquidindex+2maxindex]
    δend =          u[4*maxliquidindex+2maxindex+1: 4*maxliquidindex+3maxindex]
    Lfilm_start =   u[4*maxliquidindex+3maxindex+1: 4*maxliquidindex+4maxindex]
    Lfilm_end =     u[4*maxliquidindex+4maxindex+1: 4*maxliquidindex+5maxindex]
    return Xp,dXdt,M,δstart,δend,Lfilm_start,Lfilm_end
end

"""
Take a complete state vector in the V1 specification and return a flat V2 state vector.
Convert only the mechanical prefix; preserve all temperature separators and profiles unchanged.
The block lengths below exclude the separator-delimited temperature suffix.

V1 State Vector
- u =  [Xp; dXdt; M; δstart; δend; Lfilm_start; Lfilm_end]
- Xp is tuple, first value left interface, second value is right interface
- dXdt is tuple, duplicated information in both values
- total length is 4*NL + 5*NV

V2 State Vector
- u = [ζ_L; ζ_R; Q_B; M; Af_left; Af_right; Lfilm_left; Lfilm_right]
- total length is 3*NL + 5*NV
"""
function reparse_u_v1_as_u_v2(u_v1::AbstractVector, Ac::Real, Dh::Real)
    separator = findfirst(==(SEPERATION_VAR), u_v1)
    isnothing(separator) && throw(ArgumentError("V1 state must contain a temperature separator"))
    mechanical = Vector{Float64}(u_v1[1:separator-1])
    temperatures = u_v1[separator:end]
    Xp,dXdt,M,δstart,δend,Lfilm_start,Lfilm_end = u_v1_to_XMδL(mechanical)
    ζ_L = Ac * [first(ξp) for ξp in Xp]
    ζ_R = Ac * [last(ξp) for ξp in Xp]
    Q_B = Ac * [first(Vp) for Vp in dXdt]
    Af_L = getδarea(Ac,Dh,δstart)
    Af_R = getδarea(Ac,Dh,δend)
    u_v2 = [ζ_L; ζ_R; Q_B; M; Af_L; Af_R; Lfilm_start; Lfilm_end; temperatures]
    return u_v2
end