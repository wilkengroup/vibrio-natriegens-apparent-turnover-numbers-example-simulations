using Serialization, FastDifferentiation, Statistics
using VibrioNatriegens, COBREXA, JSONFBCModels, AbstractFBCModels
using Gurobi, CairoMakie, AlgebraOfGraphics
using DifferentiableMetabolism, ConstraintTrees
using LinearAlgebra, DataFrames, DataFramesMeta
using Optim, LineSearches, Dates, NLSolversBase

const Ex = FastDifferentiation.Node
const unit_correction = (3600 * 1e-6)
const num_threads = 1
BLAS.set_num_threads(num_threads)

include("utils.jl")

# choose from Glucose, Succinate, Alanine, Ribose, Glutamate, Acetate
condition = "Glucose"

prot_mass_fraction = 0.45 # g/gDW
membrane_fraction = 0.12 # of total proteome
atpm_measured = 97.5 # from 13C data

# import molar masses
mw = deserialize("tmp/molar-masses.js")

# load model and convert to a easy-to-use format
model =
    convert(AbstractFBCModels.CanonicalModel.Model, deserialize("tmp/$condition-model.js"))

# load proteomics data
massfd = deserialize("tmp/$condition-mass-fractions.js")
molfd = deserialize("tmp/$condition-mole-fractions.js")

# determine total amount of protein mass accounted for in the model
total_proteome_in_model =
    sum(massfd[gid] for gid in keys(model.genes) if haskey(massfd, gid))
total_capacity = prot_mass_fraction * total_proteome_in_model * 10^6 # ug/gDW [this is to rescale the units]

# load condition specific data
cond_fluxes_dict = deserialize("tmp/$condition-fluxes.js")
mu = cond_fluxes_dict[:BIOMASS]
delete!(cond_fluxes_dict, :BIOMASS)

# flux observations
Nfluxes =
    count(true for (vid, vhat) in cond_fluxes_dict if haskey(model.reactions, string(vid)))

# get initial parameter values
reaction_isozymes = deserialize("tmp/$condition-reaction-isozymes.js")

# membrane bound
membrane_gids = get_membrane_genes(model)

# protein observations
measured_prots = collect(gid for (gid, molf) in molfd if haskey(model.genes, string(gid)))

# find unmeasured and not involved in a complex that is recorded
unmeasured_prots = setdiff(Symbol.(collect(keys(model.genes))), measured_prots)

get_stoich(iso) = Symbol.(collect(keys(iso.gene_product_stoichiometry)))

protein_in_active_isozymes = unique(
    vcat(
        [
            get_stoich(iso) for isos in values(reaction_isozymes) for
            iso in values(isos) if any(in.(get_stoich(iso), Ref(measured_prots)))
        ]...,
    ),
)

# at least one protein from an isozyme must be present, ideally the largest.
# If a subunit is used in other enzymes, then use the next biggest one to make sure each isozyme is represented
measured_biggest = Symbol[]
for isos in values(reaction_isozymes)
    for iso in values(isos)
        gids = get_stoich(iso)
        mfs = [get(massfd, string(gid), 0) for gid in gids]
        idxs = sortperm(mfs, rev = true)
        for idx in idxs
            mfs[idx] == 0 && break
            gid = gids[idx]
            if !(gid in measured_biggest)
                push!(measured_biggest, gid)
                break
            end
        end
    end
end

# to this measured_biggest group, also add the proteins in the upper quartile
pthresh = quantile(collect(values(massfd)), 0.75)
append!(measured_biggest, [Symbol(k) for (k, v) in massfd if v >= pthresh])

measured_prots = intersect(measured_prots, measured_biggest)

excluded_prots = setdiff(unmeasured_prots, protein_in_active_isozymes) # these should not be descended on
rxns_no_descent =
    vcat([find_rxns_with_gene(model, string(gid)) for gid in excluded_prots]...)
for rid in rxns_no_descent # these are freebies for the model, but it causes problems when they are not excluded
    delete!(reaction_isozymes, rid)
end

Nprots = length(measured_prots)

# make parameterized version
parameter_symbols = Dict(Symbol(rid) => Ex(Symbol(rid)) for rid in keys(reaction_isozymes))

parameter_values = Dict(
    Symbol(rid) => begin
        kf = first(values(iso)).kcat_forward
        kr = first(values(iso)).kcat_reverse
        isnothing(kf) ? kr : kf
    end for (rid, iso) in reaction_isozymes
)

parameter_reaction_isozymes = Dict(
    rid => Dict(
        iso_id => COBREXA.IsozymeT{Ex}(
            iso.gene_product_stoichiometry,
            isnothing(iso.kcat_forward) ? nothing : parameter_symbols[Symbol(rid)],
            isnothing(iso.kcat_reverse) ? nothing : parameter_symbols[Symbol(rid)],
        ) for (iso_id, iso) in isos
    ) for (rid, isos) in reaction_isozymes
)

# simplify constraint tree model to get rid of bidirectional hangovers
enz = unidirectional_enzyme_constrained_flux_balance_constraints(
    model;
    reaction_isozymes = parameter_reaction_isozymes,
    gene_product_molar_masses = mw,
    capacity = Dict(
        :total => (
            Symbol.(AbstractFBCModels.genes(model)),
            ConstraintTrees.Between(0.0, total_capacity),
        ),
        :membrane => (
            Symbol.(membrane_gids),
            ConstraintTrees.Between(0.0, membrane_fraction * total_capacity),
        ),
    ),
)
enz.fluxes.ATPM.bound = ConstraintTrees.Between(0, 1000) # this makes the optimization smoother when in the objective instead of as a hard constraint

# add loss function in pieces

# add a function to account for 13c mfa special cases in case they are still in model
special_fluxes(enz, rid) = begin
    if rid == Symbol("15841") # ZWF
        if haskey(enz.fluxes, Symbol("38215"))
            ConstraintTrees.value(enz.fluxes[rid]) +
            ConstraintTrees.value(enz.fluxes[Symbol("38215")])
        else
            ConstraintTrees.value(enz.fluxes[rid])
        end
    elseif rid == Symbol("10116") # GND
        if haskey(enz.fluxes, Symbol("33023"))
            ConstraintTrees.value(enz.fluxes[rid]) +
            ConstraintTrees.value(enz.fluxes[Symbol("33023")])
        else
            ConstraintTrees.value(enz.fluxes[rid])
        end
    elseif rid == Symbol("18253") # ME
        if haskey(enz.fluxes, Symbol("12653"))
            ConstraintTrees.value(enz.fluxes[rid]) +
            ConstraintTrees.value(enz.fluxes[Symbol("12653")])
        else
            ConstraintTrees.value(enz.fluxes[rid])
        end
    else
        ConstraintTrees.value(enz.fluxes[rid])
    end
end

enz *=
    :loss_flux^ConstraintTrees.Constraint(; # squared relative error mean
        value = 1 / Nfluxes * sum(
            ConstraintTrees.squared(1 - special_fluxes(enz, vid) / vhat) for
            (vid, vhat) in cond_fluxes_dict if haskey(enz.fluxes, vid)
        ),
        bound = nothing,
    )

enz *=
    :loss_growth^ConstraintTrees.Constraint(; # squared relative error mean
        value = ConstraintTrees.squared(1 - ConstraintTrees.value(enz.fluxes.BIOMASS) / mu),
        bound = nothing,
    )

enz *=
    :loss_measured_proteome^ConstraintTrees.Constraint(; # squared relative error mean
        value = 1 / Nprots * sum(
            ConstraintTrees.squared(
                1 - ConstraintTrees.value(enz.gene_product_amounts[gid]) / molfd[gid],
            ) for gid in measured_prots if haskey(enz.gene_product_amounts, gid)
        ),
        bound = nothing,
    )

enz *=
    :loss_atpm^ConstraintTrees.Constraint(;
        value = ConstraintTrees.squared(1 - enz.fluxes.ATPM.value / atpm_measured),
        bound = nothing,
    )

enz *=
    :loss_total^ConstraintTrees.Constraint(;
        value = enz.loss_flux.value +
                enz.loss_measured_proteome.value +
                enz.loss_growth.value +
                enz.loss_atpm.value,
        bound = nothing,
    )

# test if everything works
sol = DifferentiableMetabolism.optimized_values(
    enz,
    parameter_values;
    objective = enz.loss_total.value,#
    optimizer = Gurobi.Optimizer,
    sense = COBREXA.Minimal,
)
isnothing(sol) && @warn "Descent infeasibility $variant"

parameters = collect(keys(parameter_symbols))

enz_kkt, _ = DifferentiableMetabolism.differentiate_prepare_kkt(
    enz,
    enz.loss_total.value,
    parameters;
    additional_derivatives = Dict(:loss_total => enz.loss_total.value), # outer problem
);

(f_A, f_B, fs, parameters, _A, _B, _xs) = enz_kkt;


kapp_lb = 0.1 # 1/s
kapp_ub = 1e7 # 1/s

GRB_ENV = Gurobi.Env(
    Dict{String,Any}("NumericFocus" => 3, "OutputFlag" => 0, "ThreadLimit" => num_threads),
)
opt = () -> Gurobi.Optimizer(GRB_ENV)

save_dir = "tmp/$condition/"
isdir(save_dir) || mkdir(save_dir)

# init
pvs = [log(parameter_values[p]) for p in parameters]
iterstate = [0]
loss_state = Float64[]
tstart = now()

# add box constraints to make optimization easier
lower = fill(log(kapp_lb * unit_correction), length(pvs))
upper = fill(log(kapp_ub * unit_correction), length(pvs))

# create storate containers
flux_res = DataFrame(Iteration = Int[], Reaction = String[], Flux = Float64[])
prot_res = DataFrame(Iteration = Int[], Protein = String[], Concentration = Float64[])
loss_res = DataFrame(Iteration = Int[], LossType = String[], LossValue = Float64[])
kapp_res = DataFrame(Iteration = Int[], Reaction = String[], Kapp = Float64[])

# create optimization program
inner_opt = LBFGS(; m = 5)
opts = Optim.Options(
    iterations = 2_000,
    outer_iterations = 1_000,
    f_abstol = 1e-8,
    g_tol = 1e-8,
    outer_g_abstol = 1e-8,
    time_limit = 60 * 20, # timeout after 20 minutes
    callback = x -> cb(
        x,
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
    ),
)

# finally optimize!
res = optimize(
    NLSolversBase.only_fg!((F, G, x) -> fg!(F, G, x, enz_kkt, enz, opt, parameters)),
    lower,
    upper,
    pvs,
    Fminbox(inner_opt),
    opts,
)

# inspect results (the flux, proteome, and kapp results can also be plotted like shown below)
(data(loss_res) * mapping(:Iteration, :LossValue, row = :LossType) * visual(Lines)) |> draw(
    axis = (xlabel = "Iteration", ylabel = "Loss", yscale = log10),
    facet = (linkyaxes = :none,),
    figure = (size = (800, 800),),
)
