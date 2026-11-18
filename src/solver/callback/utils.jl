# All of the callbacks expect interface values in arc coordinates
# and paired together in the old [(L,R)] format
# these adapters oblige, so we can keep changes as small as possible at first

_ξp(sys::UniformPhpSystem) =
    collect(zip(sys.liquid.ζ_L ./ sys.tube.Ac, sys.liquid.ζ_R ./ sys.tube.Ac))

function _set_ξp!(sys::UniformPhpSystem, interfaces)
    sys.liquid.ζ_L = sys.tube.Ac .* first.(interfaces)
    sys.liquid.ζ_R = sys.tube.Ac .* last.(interfaces)
    return sys
end

_reject_nonuniform(sys::UniformPhpSystem) = sys