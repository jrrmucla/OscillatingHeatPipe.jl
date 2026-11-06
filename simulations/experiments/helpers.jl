export InsulatedStaticPolygon

using OscillatingHeatPipe

function get_qbplus(t,x,base_cache,phys_params,motions)
    nrm = normals(base_cache)
    qbplus = zeros_surface(base_cache)
    return qbplus
end

function get_qbminus(t,x,base_cache,phys_params,motions)
    nrm = normals(base_cache)
    qbminus = zeros_surface(base_cache)
    return qbminus
end

function InsulatedStaticPolygon(xbound, ybound, Δs)
    body = Polygon(xbound,ybound,Δs)
    
    # stationary rigid-body motion object for the object
    X = MotionTransform([0,0],0)
    joint = Joint(X)
    m = RigidBodyMotion(joint,body)
    x = zero_motion_state(body,m)
    update_body!(body,x,m)

    # no heat flux at the boundaries
    bcdict = Dict("exterior" => get_qbplus,"interior" => get_qbminus)
    return body, m, bcdict
end