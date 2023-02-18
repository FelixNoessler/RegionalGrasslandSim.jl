"""
    one_day!(du,u,p,t)

Change in the state variables during one day.
"""
function one_day!(du,u,p,t)

    for pa in 1:p.npatches

        # -------------- potential growth
        pot_growth = Growth.potential_growth(
            p, 
            u.Biomass[pa, :]
        )
        du.Biomass[pa, :] += pot_growth

        for i in 1:p.nspecies
            du.Biomass[pa, i] = pot_growth[i] #-0.02 * u.Biomass[pa,i]
        end

        du.Water[pa] = -rand()*2
    end

    return nothing
end

function initial_conditions(p)
    return ca.ComponentVector(
        Biomass=fill(
            100.0 + 2*rand(), 
            p.npatches, 
            p.nspecies
        ),
        Water=fill(5.0, p.npatches)
    )
end

function initialize_parameters()
    return (
        npatches=2,
        nspecies=5,
        species = (
            sla = [2.213,3,4,5,6],
            ldmc = [0.4, 0.2, 0.3, 0.4, 0.2]
        )
    )
end

"""
    discrete_prob()

Initialise the difference equation problem.
"""
function discrete_prob()
    p = initialize_parameters()

    return DiscreteProblem(
        one_day!,
        initial_conditions(p),
        (0, 100),
        p;
    )
end


function transform_solution(sol, prob)
    b = hcat([u.Biomass for u in sol.u]...)
    b = reshape(b, prob.p.npatches, prob.p.nspecies, :);

    w = hcat([u.Water for u in sol.u]...)

    return ca.ComponentVector(
        Biomass=b, 
        Water=w
    )
end