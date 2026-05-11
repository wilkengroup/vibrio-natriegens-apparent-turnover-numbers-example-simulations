using DataFrames,
    DataFramesMeta, CSV, Statistics, CairoMakie, AlgebraOfGraphics, Serialization
using VibrioNatriegens, Gurobi, COBREXA, AbstractFBCModels, JSONFBCModels, ConstraintTrees
const scale_kcat = 3600 * 1e-6
include("utils.jl")

df = DataFrame(
    Condition = String[],
    Kinetics = String[],
    Growth = Float64[],
    MembraneCapacity = Float64[],
    TotalCapacity = Float64[],
)

# import molar masses
mw = deserialize("tmp/molar-masses.js")

for condition in ["Glucose", "Succinate", "Alanine", "Ribose", "Glutamate", "Acetate"]

    # load kinetics
    kapps = CSV.File("data/best_fit_kapps.csv") |> DataFrame
    @rsubset!(kapps, :Condition == condition,)
    kapps = @by kapps :RID begin
        :Kapp = mean(:Kapp)
    end
    kapp_dict = Dict(r.RID => r.Kapp for r in eachrow(kapps))

    kmaxs = CSV.File("data/best_fit_kapps.csv") |> DataFrame
    kmaxs = @by kmaxs :RID begin
        :Kmax = maximum(:Kapp)
    end
    kmax_dict = Dict(r.RID => r.Kmax for r in eachrow(kmaxs))

    # load model and add constraints
    model = convert(
        AbstractFBCModels.CanonicalModel.Model,
        deserialize("tmp/$condition-model.js"),
    )

    cond_fluxes_dict = deserialize("tmp/$condition-fluxes.js")
    for (k, v) in cond_fluxes_dict
        !startswith(string(k), "EX_") && continue
        model.reactions[string(k)].lower_bound = v
        model.reactions[string(k)].upper_bound = v
    end

    membrane_gids = get_membrane_genes(model)
    total_capacity = 160_000.0 # ug/gDW 
    membrane_fraction = 0.12

    # update kinetics to best fit values
    reaction_isozymes_kmax = deserialize("tmp/$condition-reaction-isozymes.js")
    for (k, v) in kmax_dict
        if haskey(reaction_isozymes_kmax, k)
            for iso in values(reaction_isozymes_kmax[k])
                iso.kcat_forward = v
                iso.kcat_reverse = v
            end
        end
    end

    reaction_isozymes_kapp = deserialize("tmp/$condition-reaction-isozymes.js")
    for (k, v) in kapp_dict
        if haskey(reaction_isozymes_kapp, k)
            for iso in values(reaction_isozymes_kapp[k])
                iso.kcat_forward = v
                iso.kcat_reverse = v
            end
        end
    end

    # run simulations
    for (kinetics_id, kinetics) in
        [("kmax", reaction_isozymes_kmax), ("kapp", reaction_isozymes_kapp)]

        sol = enzyme_constrained_flux_balance_analysis(
            model,
            optimizer = Gurobi.Optimizer,
            reaction_isozymes = kinetics,
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

        push!(
            df,
            (
                condition,
                kinetics_id,
                isnothing(sol) ? missing : sol.objective,
                isnothing(sol) ? missing :
                sol.gene_product_capacity.membrane / (total_capacity*membrane_fraction),
                isnothing(sol) ? missing : sol.gene_product_capacity.total / total_capacity,
            ),
        )
    end
end

# plot the results
data(df) *
mapping(:Condition, :Growth, color = :Kinetics, dodge = :Kinetics) *
visual(BarPlot) |> draw
