using CSV, DataFrames, DataFramesMeta, Serialization
using VibrioNatriegens, COBREXA, JSONFBCModels, AbstractFBCModels

isdir("tmp") || mkdir("tmp") # create dir if necessary

# save molar mass data
mw = Dict(
    CSV.File(
        "data/vnat_molar_masses.csv",
        drop = ["Length"],
        types = [String, Float64, Float64],
    ),
)
serialize("tmp/molar-masses.js", mw)

# prepare condition specific data
for condition in [
    "Acetate"
    "Alanine"
    "Glutamate"
    "Succinate"
    "Ribose"
    "Glucose"
]

    # serialize the condition specific models
    model = AbstractFBCModels.load(
        JSONFBCModels.JSONFBCModel,
        "data/"*condition*"-vibrio-natriegens-model.json",
    )
    serialize("tmp/"*condition*"-model.js", model)

    # serialize a specific flux sample [units: mmol/gDW/h]
    fluxes = Dict(
        Symbol(r.Reaction) => r.Flux for r in CSV.File("data/"*condition*"-fluxes.csv")
    )
    serialize("tmp/"*condition*"-fluxes.js", fluxes)

    # serialize the proteomics mass fraction data
    mass = Dict(
        r.Protein => r.MassFraction for
        r in CSV.File("data/"*condition*"-proteomics-mass-fraction.csv")
    )
    serialize("tmp/"*condition*"-mass-fractions.js", mass)

    # serialize the proteomics absolute mole fraction data [units: ng/gDW]
    mole = Dict(
        Symbol(r.Protein) => r.AbsoluteMoleFraction for
        r in CSV.File("data/"*condition*"-proteomics-absolute-mole-fraction.csv")
    )
    serialize("tmp/"*condition*"-mole-fractions.js", mole)

    # prepare the reaction isozyme data which includes the gene reaction rules and kinetics [kinetics unit: M/h]
    reaction_isozymes = Dict{String,Dict{String,COBREXA.Isozyme}}()
    for r in CSV.File("data/"*condition*"-reaction-isozymes.csv")
        if haskey(reaction_isozymes, r.Reaction)
            if haskey(reaction_isozymes[r.Reaction], r.Isozyme)
                iso = reaction_isozymes[r.Reaction][r.Isozyme]
                iso.gene_product_stoichiometry[r.Protein] = r.Stoichiometry
            else
                iso = COBREXA.Isozyme(
                    Dict(r.Protein => r.Stoichiometry),
                    r.KcatForward,
                    r.KcatReverse,
                )
                reaction_isozymes[r.Reaction] = Dict(r.Isozyme => iso)
            end
        else
            isos = Dict{String,COBREXA.Isozyme}()
            isos[r.Isozyme] = COBREXA.Isozyme(
                Dict(r.Protein => r.Stoichiometry),
                r.KcatForward,
                r.KcatReverse,
            )
            reaction_isozymes[r.Reaction] = isos
        end
    end
    reaction_isozymes
    serialize("tmp/"*condition*"-reaction-isozymes.js", reaction_isozymes)
end
