"""
    one_day!(du,u,p,t)

Describes the change in the state variables during one day.
"""
function one_day!(du, u, p, t)

    adj = 365
    t_adj = (t - 1) - ((t - 1) ÷ adj * adj) + 1

    # --------------------- biomass dynamics
    if t_adj < 80 || t_adj > 300
        du.Biomass .= 0
    else
        for pa in 1:p.npatches
            # -------------- potential growth
            pot_growth = Growth.potential_growth(
                p, 
                u.Biomass[pa, :];
                PAR = PAR[t_adj - 79]
            )

            # -------------- senescence
            sen = Growth.senescence(; 
                biomass=u.Biomass[pa, :], 
                t=t_adj,
                p)

            for i in 1:p.nspecies
                du.Biomass[pa, i] = pot_growth[i] - sen[i]
            end
        end
    end

    # --------------------- water dynamics
    for pa in 1:p.npatches
        du.Water[pa] = - 0.01 .* u.Water[pa]
    end

    return nothing
end

function initial_conditions(p)
    return ca.ComponentVector(
        Biomass=fill(
            10.0  + 1000*rand(), 
            p.npatches, 
            p.nspecies
        ),
        Water=fill(5.0, p.npatches)
    )
end

function initialize_parameters()

    slas = [2.213,3,4,5,6]
    ldmcs =  [0.4, 0.2, 0.3, 0.4, 0.2]
    μs = (0.023 .+ 0.382 .* slas .- 0.0000629 .* ldmcs) .^ (-1)

    p = (
        npatches=2,
        nspecies=5,
        species = (
            sla = slas,
            ldmc = ldmcs,
            μ = μs
        ))

    return p
end

"""
    discrete_prob()

Initialise the difference equation problem.
"""
function discrete_prob(; tmax=100, input_file)
    global PAR = load(input_file, "PAR")
    
    p = initialize_parameters()

    return prob = DiscreteProblem(
        one_day!,
        initial_conditions(p),
        (0, tmax),
        p;
    )
end

"""
    run_model(prob)

Solve the difference equation and return a solution.

See also [`transform_solution`](@ref) for further details about the solution type.
"""
function run_model(prob)
    sol = solve(prob, FunctionMap(scale_by_time=true))
    d = transform_solution(sol, prob)

    return d
end


"""
    transform_solution(sol, prob)

Transfroms `ODESolution` into a `NamedTuple` that can be easily used for further analysis. 

The `NamedTuple` contains the biomass of each species in each patch (`Biomass[patch, species, time]`), the soil water of each patch (`Water[patch, time]`), the `time` as an array in days and a `NamedTuple` of all parameters (`p`). 

See also [`run_model`](@ref).
"""
function transform_solution(sol, prob)
    b = hcat([u.Biomass for u in sol.u]...)
    b = reshape(b, prob.p.npatches, prob.p.nspecies, :);

    w = hcat([u.Water for u in sol.u]...)

    return (
        Biomass=b, 
        Water=w,
        time=sol.t,
        p=prob.p
    )
end