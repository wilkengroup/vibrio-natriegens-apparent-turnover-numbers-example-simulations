
function get_membrane_genes(model)
    is_membrane_reaction(model, rid) = begin
        cs = unique([
            model.metabolites[m].compartment for
            m in keys(model.reactions[rid].stoichiometry)
        ])
        length(cs) != 1
    end

    rids = [rid for rid in keys(model.reactions) if is_membrane_reaction(model, rid)]
    gids = String[]
    for rid in rids
        isos = model.reactions[rid].gene_association_dnf
        isnothing(isos) && continue
        for gs in isos
            append!(gids, gs)
        end
    end
    unique(gids)
end

function find_rxns_with_gene(model, gid)
    rxns = String[]

    for (rid, rxn) in model.reactions
        grrs = rxn.gene_association_dnf
        isnothing(grrs) && continue
        for grr in grrs
            for g in grr
                if gid == g
                    push!(rxns, rid)
                end
            end
        end
    end

    unique!(rxns)
    rxns
end


function unidirectional_enzyme_constrained_flux_balance_constraints(
    model::AbstractFBCModels.AbstractFBCModel;
    reaction_isozymes::Dict{String,Dict{String,COBREXA.IsozymeT{R}}},
    gene_product_molar_masses::Dict{String,Float64},
    capacity,
) where {R<:Real}
    # prepare some accessor functions for the later stuff
    # TODO: might be nicer to somehow parametrize the fwd/rev directions out.
    # Also there is a lot of conversion between symbols and strings, might be
    # nicer to have that sorted out in some better way.
    function isozyme_forward_ids(rid)
        haskey(reaction_isozymes, String(rid)) || return nothing
        return [
            Symbol(k) for
            (k, i) in reaction_isozymes[String(rid)] if !isnothing(i.kcat_forward)
        ]
    end
    function isozyme_reverse_ids(rid)
        haskey(reaction_isozymes, String(rid)) || return nothing
        return [
            Symbol(k) for
            (k, i) in reaction_isozymes[String(rid)] if !isnothing(i.kcat_reverse)
        ]
    end
    kcat_forward(rid, iso_id) = reaction_isozymes[String(rid)][String(iso_id)].kcat_forward
    kcat_reverse(rid, iso_id) = reaction_isozymes[String(rid)][String(iso_id)].kcat_reverse
    isozyme_gene_product_stoichiometry(rid, iso_id) = Dict(
        Symbol(k) => v for (k, v) in
        reaction_isozymes[String(rid)][String(iso_id)].gene_product_stoichiometry
    )
    gene_ids = Symbol.(keys(gene_product_molar_masses))
    gene_product_molar_mass(gid) = get(gene_product_molar_masses, String(gid), 0.0)

    isforward(c) = c.bound.lower >= 0

    # allocate all variables and build the system
    constraints = COBREXA.flux_balance_constraints(model)

    constraints *=
        :fluxes_forward^ConstraintTrees.ConstraintTree(
            k => c for (k, c) in constraints.fluxes if k ∉ [:BIOMASS, :ATPM] && isforward(c)
        )
    constraints *=
        :fluxes_reverse^ConstraintTrees.ConstraintTree(
            k => ConstraintTrees.Constraint(-c.value, ConstraintTrees.Between(0, 1000)) for
            (k, c) in constraints.fluxes if k ∉ [:BIOMASS, :ATPM] && !isforward(c)
        )

    constraints += enzyme_variables(;
        fluxes_forward = constraints.fluxes_forward,
        fluxes_reverse = constraints.fluxes_reverse,
        isozyme_forward_ids,
        isozyme_reverse_ids,
    )

    return constraints * enzyme_constraints(;
        fluxes_forward = constraints.fluxes_forward,
        fluxes_reverse = constraints.fluxes_reverse,
        isozyme_forward_amounts = constraints.isozyme_forward_amounts,
        isozyme_reverse_amounts = constraints.isozyme_reverse_amounts,
        kcat_forward,
        kcat_reverse,
        isozyme_gene_product_stoichiometry,
        gene_product_molar_mass,
        capacity_limits = expand_enzyme_capacity(capacity, gene_ids),
    )
end

function fg!(F, G, pvs, enz_kkt, enz, opt, parameters)
    parameter_values = Dict(k => exp(v) for (k, v) in zip(parameters, pvs)) # update

    _sol = DifferentiableMetabolism.optimized_values(
        enz,
        parameter_values;
        objective = enz.loss_total.value,
        optimizer = opt,
        sense = COBREXA.Minimal,
    )

    if !isnothing(G)
        dx_dk = DifferentiableMetabolism.differentiate_solution(
            enz_kkt,
            _sol.primal_values,
            _sol.equality_dual_values,
            _sol.inequality_dual_values,
            parameter_values,
        ) # derivative of variables wrt parameters

        dL_dx = DifferentiableMetabolism.differentiate_function(
            enz_kkt,
            :loss_total,
            _sol.primal_values,
            _sol.equality_dual_values,
            _sol.inequality_dual_values,
            parameter_values,
        ) # derivative of loss function wrt variables

        G .= ((dx_dk' * dL_dx) .* exp.(pvs)) # gradient

    end

    if !isnothing(F)
        return _sol.tree.loss_total
    end
end

function cb(
    st,
    enz,
    parameters,
    opt,
    iterstate,
    loss_state,
    tstart,
    flux_res,
    prot_res,
    loss_res,
    kapp_res;
    save_every = 10,
)

    # prevent getting stuck
    if now() - tstart > Millisecond(1000 * 60 * 20) # 10 min max per problem
        return true
    end

    # update iteration counter
    iterstate[1] += 1
    iter = first(iterstate)

    if !(iter <= 10 || iter % save_every == 0)
        return false
    end

    # solve problem
    parameter_values = Dict(k => exp(v) for (k, v) in zip(parameters, st.x)) # update

    _sol = DifferentiableMetabolism.optimized_values(
        enz,
        parameter_values;
        objective = enz.loss_total.value,
        optimizer = opt,
        sense = COBREXA.Minimal,
    )

    obj = _sol.tree.loss_total

    # write results
    loss_res_dict = Dict(
        id => _sol.tree[id] for
        id in filter(startswith("loss_"), string.(collect(keys(_sol.tree))))
    )
    for (k, v) in loss_res_dict
        push!(loss_res, (iter, k, v))
    end

    flux_res_dict = Dict(string(k) => v for (k, v) in _sol.tree.fluxes)
    for (k, v) in flux_res_dict
        push!(flux_res, (iter, k, v))
    end

    prot_res_dict = Dict(string(k) => v for (k, v) in _sol.tree.gene_product_amounts)
    for (k, v) in prot_res_dict
        push!(prot_res, (iter, k, v))
    end


    for (i, p) in enumerate(parameters)
        push!(kapp_res, (iter, string(p), parameter_values[p] / unit_correction))
    end

    @info """Progress:
    Iterations= $iter
    Objective = $obj"""

    # termination
    push!(loss_state, obj)
    if length(loss_state) > 2
        l1 = loss_state[end]
        l2 = loss_state[end-1]
        abs(l1 - l2) <= 1e-8 && return true
    end

    false
end
