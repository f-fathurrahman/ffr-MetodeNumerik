using LinearAlgebra: I

mutable struct mpmMaterialPoint_2D_Classic
    mass::Float64
    initial_volume::Float64
    volume::Float64
    elastic_modulus::Float64
    poisson_ratio::Float64
    yield_stress::Float64
    centroid::Vector{Float64}
    velocity::Vector{Float64}
    momentum::Vector{Float64}
    ext_force::Vector{Float64}
    # 0.0=no restraint, 1.0=fully restrained
    restraint::Array{Float64}
    corner::Matrix{Float64} # corner position
    deform_grad::Matrix{Float64}
    deform_grad_incr::Matrix{Float64}
    strain::Array{Float64} # xx, yy, zz, xy, yz, zx
    plastic_strain::Array{Float64} # xx, yy, zz, xy, yz, zx
    equiv_plastic_strain::Float64 # equivalent plastic strain
    stress::Array{Float64}
end

# constructor
function mpmMaterialPoint_2D_Classic()

    mass = 1.0
    initial_volume = 1.0
    volume = 1.0
    
    elastic_modulus = 1.0
    poisson_ratio = 0.3
    yield_stress = 1e24

    NDIM = 2

    centroid = zeros(Float64, NDIM)
    velocity = zeros(Float64, NDIM)
    momentum = zeros(Float64, NDIM)
    ext_force = zeros(Float64, NDIM)
    restraint = zeros(Float64, NDIM) # FIXME: change to Bool?

    corner = zeros(Float64, NDIM, 4)
    deform_grad = Matrix{Float64}(I(2))
    deform_grad_incr = Matrix{Float64}(I(2))
    strain = zeros(Float64, 3)
    plastic_strain = zeros(Float64, 3)
    equiv_plastic_strain = 0.0
    stress = zeros(Float64, 3)

    return mpmMaterialPoint_2D_Classic(
        mass,
        initial_volume,
        volume,
        elastic_modulus,
        poisson_ratio,
        yield_stress,
        centroid,
        velocity,
        momentum,
        ext_force,
        restraint,
        corner,
        deform_grad,
        deform_grad_incr,
        strain,
        plastic_strain,
        equiv_plastic_strain,
        stress
    )

end


function createMaterialDomain_Circle(
    fCenter::Array{Float64},
    fRadius::Float64,
    fOffset::Float64
)
    # just in case radius is not a multiple of offset
    fRadius = floor(fRadius/fOffset) * fOffset
    Npoints = 0
    for fy = -fRadius-0.5*fOffset:fOffset:+fRadius+0.5*fOffset
        for fx = -fRadius-0.5*fOffset:fOffset:+fRadius+0.5*fOffset
            if(fx^2 + fy^2 <= fRadius^2)
                Npoints +=1
            end
        end
    end
    println("Npoints circle = ", Npoints)

    thisMaterialDomain = Vector{mpmMaterialPoint_2D_Classic}(undef,Npoints)
    ip = 0
    for fy = -fRadius-0.5*fOffset:fOffset:+fRadius+0.5*fOffset
        for fx = -fRadius-0.5*fOffset:fOffset:+fRadius+0.5*fOffset
            if(fx^2 + fy^2 <= fRadius^2)
                ip += 1
                thisMaterialDomain[ip] = mpmMaterialPoint_2D_Classic()
                # Modify centroid
                thisMaterialDomain[ip].centroid[:] .= [fCenter[1] + fx; fCenter[2] + fy]
            end
        end
    end
    return thisMaterialDomain
end



function createMaterialDomain_Rectangle(
    v2Center::Array{Float64},
    fWidth::Real,
    fHeight::Real,
    fOffset::Real
)

    fWidth = floor(fWidth/fOffset) * fOffset    #just in case width is not a multiple of offset
    fHeight = floor(fHeight/fOffset) * fOffset    #just in case height is not a multiple of offset
    Npoints = 0
    for _ = -0.5*fHeight+0.5*fOffset:fOffset:+0.5*fHeight-0.5*fOffset
        for _ = -0.5*fWidth+0.5*fOffset:fOffset:+0.5*fWidth-0.5*fOffset
            Npoints += 1
        end
    end
    println("Npoints rectangle = ", Npoints)

    thisMaterialDomain = Array{mpmMaterialPoint_2D_Classic}(undef,Npoints)
    ip = 0
    for fy in -0.5*fHeight+0.5*fOffset:fOffset:+0.5*fHeight-0.5*fOffset
        for fx in -0.5*fWidth+0.5*fOffset:fOffset:+0.5*fWidth-0.5*fOffset
            ip += 1
            thisMaterialDomain[ip] = mpmMaterialPoint_2D_Classic()
            thisMaterialDomain[ip].centroid[:] .= v2Center .+ [fx; fy]
        end
    end
    return thisMaterialDomain
end


function getStressIncrement_Elastic(fE::Float64, fNu::Float64, v3StrainIncrement::Array{Float64})
    v3Result = zeros(3)

    fConstant = fE/(1.0 + fNu)/(1.0 - 2.0*fNu)

    v3Result[1] = fConstant * ((1.0-fNu)*v3StrainIncrement[1] + fNu*v3StrainIncrement[2])
    v3Result[2] = fConstant * ((1.0-fNu)*v3StrainIncrement[2] + fNu*v3StrainIncrement[1])
    v3Result[3] = fConstant * ((0.5-fNu)*v3StrainIncrement[3])

    return v3Result
end



function getIncrement_Plastic(
    fE::Float64,
    fNu::Float64,
    fK0::Float64,
    fAlphaCurrent::Float64,
    v3StressCurrent::Array{Float64},
    v3StrainCurrent::Array{Float64},
    v3PlasticStrainCurrent::Array{Float64},
    v3StrainIncrement::Array{Float64}
)
    v32Result = zeros(3,2)

    v3StressIncrement = zeros(3)
    v3PlasticStrainIncrement = zeros(3)
    fAlphaIncrement = 0.0

    mu      = fE/2.0/(1.0+fNu);    # shear modulus
    lambda  = fE*fNu/((1.0+fNu)*(1.0-2.0*fNu));
    kappa   = lambda + mu;

    eye2   = [1.0; 1.0; 0.0];
    eye2x2 = [1.0 1.0 0.0; 1.0 1.0 0.0; 0.0 0.0 0.0];
    I_dev  = Matrix(I(3)) - 0.5*eye2x2;
    I33      = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 0.5]; #0.5 to make engineering strain to physical one
    Iinv   = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 2.0];

    # Compute trial stress
    epsilon_dev  = I_dev * v3StrainCurrent;
    s_trial      = 2.0 * mu * I33 * (epsilon_dev - v3PlasticStrainCurrent);
    norm_s_trial = sqrt(s_trial[1]^2 + s_trial[2]^2 + 2*s_trial[3]^2);
    sigma_trial  = kappa*sum(v3StrainCurrent[1] + v3StrainCurrent[2])*eye2 + s_trial;

    # Check yield condition
    # alpha0 = 0.0 # sina, perfect plasticity for now
    k1 = 0.0 #sina, perfect plasticity for now
    f_trial = norm_s_trial - (k1*fAlphaCurrent + fK0);

    if f_trial <= 0.0 # elastic update
        fAlphaIncrement = 0.0
        v3PlasticStrainIncrement = zeros(3)
        v3StressIncrement = sigma_trial - v3StressCurrent
    else # plastic step
        normal = s_trial/norm_s_trial;
        lambda = (norm_s_trial - k1*fAlphaCurrent - fK0)/(2.0*mu + k1);
        fAlphaIncrement = lambda
        # Update plastic strain and stress
        v3PlasticStrainIncrement = lambda*Iinv*normal;
        v3StressIncrement = kappa*sum(v3StrainCurrent[1] + v3StrainCurrent[2])*eye2 + s_trial - 2.0*mu*lambda*normal;
        v3StressIncrement -= v3StressCurrent
    end

    # why use hcat?
    v32Result = hcat(v3StressIncrement, v3PlasticStrainIncrement)
    v32Result = hcat(v32Result, [fAlphaIncrement; fAlphaIncrement; fAlphaIncrement])
    return v32Result
end

