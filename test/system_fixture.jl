using OscillatingHeatPipe

"""
    get_adiabatic_system_fixture(Dh, type; tubeshape="square", kwargs...)

Return a UniformPhpSystem for `:UNIFORM` or a VariablePhpSystem for `:VARIABLE`,
with identical constant channel geometry and an adiabatic plate. Reset the random
seed before each call when comparing initial states.
"""
function get_adiabatic_system_fixture(Dh, type; tubeshape="square",
        g=OscillatingHeatPipe.DEFAULT_G,
        closedornot=OscillatingHeatPipe.DEFAULT_CLOSED, kwargs...)
    type in (:UNIFORM, :VARIABLE) ||
        throw(ArgumentError("Select either :UNIFORM or :VARIABLE for the fixture type"))
    plate = get_adiabatic_plate_fixture()
    if type == :UNIFORM
        return initialize_ohpsys(plate, p_fluid, 0.0; d=Dh, tubeshape, g, closedornot, kwargs...)
    else
        shape = OscillatingHeatPipe._get_ohp_from_forcing_list(plate).shape
        L = arclength(shape)
        Pc, Ac = OscillatingHeatPipe.Pc_Ac(Dh, tubeshape)
        tube = OscillatingHeatPipe.VariableTube([0.0, L], [Dh, Dh],
            [Pc, Pc], [Ac, Ac], g, closedornot, numpts(shape))
        return initialize_ohpsys(plate, p_fluid, 0.0; tube, kwargs...)
    end
end
