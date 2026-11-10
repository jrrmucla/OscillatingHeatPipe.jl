using JLD2: jldopen

export save_results, load_results

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
        version == RESULTS_FORMAT_VERSION || throw(ArgumentError("Unsupported results format version: $version"))
        # Dispatch future supported versions here, before decoding their payloads.
        version == 1 && _load_results_v1(file, kind)
    end
end

# Version 1 of the save format
# there are two 'kinds' of result within the file:
# 'history' contains numerical values only -> used for regression testing
# 'checkpoint' contains numerical values and integrators -> used for postprocessing
# Version 1 uses the SimulationResult struct from 205ff6bab4d94a44f9c0c082e9f7e04b0580975d
# Version 2 uses the SimulationResult struct from 0484a7d65f8b94f9d86cae09ccd94396f858a6ad

function _load_results_v1(file, kind)
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
