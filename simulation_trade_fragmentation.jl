using CSV
using DataFrames
using LinearAlgebra
using Statistics


"""
This is the main function called in the background for the simulations.
It uses the Alvarez and Lucas (2007) method to find a solution.
"""
function eqm_AL(τ_hat, c_hat, ξ_hat, ϕ, ψ, X; numeraire=0)
    Y = sum(X, dims=2)
    E = transpose(sum(X, dims=1))
    N = size(X, 1)
    p = ones(N, 1) / N
    P = ones(N, 1) / N
    p_next = Array{Float64}(undef, N, 1)
    P_next = Array{Float64}(undef, N, 1)

    maxit = 10000  
    for iter in 1:maxit
        for i in 1:N
            rhs = 0.0
            for j in 1:N
                rhs += X[i, j] / Y[i] * τ_hat[i, j]^(-ϕ) * P[j]^(ϕ) * ξ_hat[j] * c_hat[j] * p[j] * (p[j]/P[j])^ψ
            end
            temp = P[i]^ψ * rhs
            p_next[i] = temp^(1/(1+ψ+ϕ))
        end
        p_next .= p_next ./ sum(p_next)
        prueba_p = norm(p_next .- p)
        if prueba_p < 1e-15
            break
        else
            p .= p_next
        end

        for i in 1:N
            rhs = 0.0
            for j in 1:N
                rhs += X[j, i] / E[i] * τ_hat[j, i]^(-ϕ) * p[j]^(-ϕ)
            end
            temp = rhs
            P_next[i] = temp^(-1/ϕ)
        end
        prueba_P = norm(P_next .- P)
        if prueba_P < 1e-15
            break
        else
            P .= P_next
        end
        if iter == maxit
            println("Reached maximum iterations")
        end
    end

    # Price normalization
    if numeraire == 0 # change in mean gross output is normalized to 1
        t = 1 / mean(p .* c_hat .* (p./P).^ψ)
        p .= t * p
        P .= t * P
    else # specific price is chosen
        p .= p ./ p[numeraire]
        for i in 1:N
            rhs = 0.0
            for j in 1:N
                rhs += X[j, i] / E[i] * τ_hat[j, i]^(-ϕ) * p[j]^(-ϕ)
            end
            temp = rhs
            P[i] = temp^(-1/ϕ)
        end
    end

    # Income change
    Y_hat = p .* c_hat .* (p./P).^ψ

    # Calculate trade flows that are correct up to scale
    E_hat = ξ_hat .* Y_hat  # initial guess for E_hat to figure out Ξ_hat
    X_hat = Array{Float64}(undef, N, N)
    for i in 1:N
        for j in 1:N
            X_hat[i, j] = τ_hat[i, j]^(-ϕ) * p[i]^(-ϕ) * P[j]^ϕ * E_hat[j]
        end
    end
    X_c = X_hat .* X  # correct up to scale
    E_c = sum(X_c, dims=2)  # correct up to scale
    Y_c = Y_hat .* Y  # true counterfactual

    # Figure out scale for E and trade flows
    excess_E = E_c ./ Y_c
    Ξ_hat = 1.0 / excess_E[1]

    # Re-scale trade flows and expenditure
    X_hat .= Ξ_hat .* X_hat
    E_hat .= Ξ_hat .* E_hat

    # Relative price
    W = p ./ P

    return p, P, W, X_hat, Y_hat, E_hat, Ξ_hat
end

"""
    counterfactual_AL(X, τ_hat, ϕ, ψ; numeraire=0)

Compute the change of all endogenous variables in response to a change in trade frictions.

### Input
* `X`: matrix of baseline trade flows
* `τ_hat`: change in trade costs (`τ_hat = τ_counterfactual/τ_baseline`)
* `ϕ`: demand elasticity
* `ψ`: supply elasticity
* the optional argument `numeraire` chooses which country is used as the numeraire. Choosing `numeraire=0` normalizes the change in gross world output to one.
### Outputs
All outputs are "hat" variables (the ratio of counterfactual value to baseline value).
* `X_hat`: trade flows
* `Y_hat`: real output
* `p_hat`: prices
* `P_hat`: price indexes
* `λ_hat`: scalar to balance
* `W_hat`: welfare
"""
function counterfactual_AL(X, τ_hat, ϕ, ψ; numeraire=0)
    N = size(X, 1)
    c_hat = ones(N, 1)
    ξ_hat = ones(N, 1)
    p_hat, P_hat, W_hat, X_hat, Y_hat, E_hat, Ξ_hat = eqm_AL(τ_hat, c_hat, ξ_hat, ϕ, ψ, X; numeraire=numeraire)
    return X_hat, Y_hat, p_hat, P_hat, Ξ_hat, W_hat
end

"""
Load the list of country ISOs and trade matrix and replace missing values with zeros.
"""
function load_X(filename)
    filepath = joinpath(@__DIR__, "source_data" , filename)
    df = DataFrame(CSV.File(filepath))
    df_zeros = coalesce.(df, 0.0)  # replace missing with zeros
    df_wide = unstack(df_zeros, :iso_o, :iso_d, :trade)
    iso, trade = extract_iso_var(df_wide)
    return iso, trade
end

"""
Extract the first column of df as ISOs and the remaining columns as a matrix of numbers.
"""
function extract_iso_var(df)
    iso = Array(df[!, names(df)[1]])
    df = df[!, names(df)[2:end]]
    @assert iso == names(df)  # Check that rows and cols are sorted alike
    return iso, Matrix{Float64}(df)
end

"""
Load a specific scenario.
"""
function load_τ_hat(scenario_number)
    filepath = joinpath(@__DIR__, "source_data" , "tau_hat.csv")
    df = DataFrame(CSV.File(filepath))
    desired_col = string("tau_hat_", scenario_number)
    df_scenario = select(df, :iso_o, :iso_d, desired_col)
    df_wide = unstack(df_scenario, :iso_o, :iso_d, desired_col)
    iso, τ_hat = extract_iso_var(df_wide)
    return iso, τ_hat
end

"""
Load bloc definitions and return them as a dictionary of the form {bloc: array of ISOs}.
"""
function load_bloc_definitions()
    filepath = joinpath(@__DIR__, "source_data" , "def_blocs.csv")
    df = DataFrame(CSV.File(filepath))
    df_blocs = select(df, :iso, :bloc)
    categories = ["Western", "Eastern", "Neutral"]
    blocs = Dict()
    for c in categories
        a = Array(df_blocs[df_blocs[!, :bloc] .== c, :iso])
        blocs[Symbol(c)] = a
    end
    return blocs
end

"""
Produce a masking array for a2 indicating whether the elements in a2 are in a1.
"""
function mask_array(a1, a2)
    N = length(a2)
    res = Array{Bool}(undef, N)
    for i in 1:N
        is_present = (a2[i] in a1)
        res[i] = is_present
    end
    return res
end

"""
Simulate a scenario.

Scenario definitions:
0: Russia in autarky
1: Increase of MATR with full data
2: Increase of MATR with PROD data
3: Increase of MATR + East exits WTO + TA break with full data
4: Increase of MATR + East exits WTO + TA break with PROD data
"""
function simulate(scenario, ϕ, ψ; numeraire=0)
    # Load appropriate data
    if scenario in [0, 1, 3]
        filename = "trade_2019.csv"
    elseif scenario in [2, 4]
        filename = "trade_prod_2018.csv"
    end
    iso, X = load_X(filename)

    if scenario == 0  # Russia in autarky
        N = length(iso)
        τ_hat = ones(Float64, N, N)
        mask_RUS = (iso .== "RUS")
        mask_not_RUS = (iso .!= "RUS")
        τ_hat[mask_RUS, mask_not_RUS] .= 1e9  # effectively infinity
        τ_hat[mask_not_RUS, mask_RUS] .= 1e9  # effectively infinity
        iso_τ = iso
    elseif scenario in [1, 3]  # scenarios with full data
        iso_τ, τ_hat = load_τ_hat(scenario)
    elseif scenario in [2, 4]  # scenarios with PROD data
        iso_τ, τ_hat = load_τ_hat(scenario-1)
        # Align trade costs with the smaller trade dataset
        mask = mask_array(iso, iso_τ)
        iso_τ = iso_τ[mask]
        τ_hat = τ_hat[mask, mask]
        mask_reverse = mask_array(iso_τ, iso)
        iso = iso[mask_reverse]
        X = X[mask_reverse, mask_reverse]
    end

    @assert iso == iso_τ  # Check that country order is the same

    X_hat, Y_hat, p_hat, P_hat, λ_hat, W_hat = counterfactual_AL(X, τ_hat, ϕ, ψ; numeraire=numeraire)
    @assert p_hat .* (p_hat ./ P_hat).^(ψ) ≈ Y_hat  # price x output = income

    X_1 = X_hat .* X  # counterfactual trade flows
    return iso, X, X_1, p_hat, P_hat
end

"""
Calculate aggregate trade flows between blocs. Assumes that the initial trade matrix is square.
"""
function bloc_trade_matrix(X, exporting_blocs, importing_blocs; partition=true, exclude_domestic_trade=true)
    N = size(X, 1)  # number of countries
    B_1 = length(exporting_blocs)  # Number of exporting blocs
    B_2 = length(importing_blocs)  # Number of importing blocs

    # Check that the trade matrix is square
    @assert size(X) == (N, N)

    if partition  # default option
        # Check that bloc_partition is a complete partition or rows/columns in the trade matrix
        all_idx_1 = sort(reduce(vcat, [b for b in exporting_blocs]))
        all_idx_2 = sort(reduce(vcat, [b for b in importing_blocs]))
        @assert all_idx_1 == Array(1:N)
        @assert all_idx_2 == Array(1:N)
    end

    if exclude_domestic_trade  # default option
        # Keep only international trade by setting the diagonal to zero
        X = (ones(N, N) - I) .* X  # does not mutate X in the function argument
    end

    # Calculate aggregate trade flows between blocs
    X_blocs = Array{Float64}(undef, B_1, B_2)
    for (i, bi) in enumerate(exporting_blocs)
        for (j, bj) in enumerate(importing_blocs)
            X_ij = X[bi, bj] 
            X_blocs[i, j] = sum(X_ij)
        end
    end
    if partition
        # Verify adding up
        @assert sum(X) ≈ sum(X_blocs)
    end
    return X_blocs
end

"""
Summarize main (non-bilateral) variables for all countries
"""
function summarize(iso, X, X_1, p_hat, P_hat, ψ)
    N = size(X, 1)

    # International trade
    trade_0 = (ones(N, N) - I) .* X
    trade_1 = (ones(N, N) - I) .* X_1

    # Domestic trade
    domestic_0 = diag(X)
    domestic_1 = diag(X_1)

    # Total exports
    tot_exports_0 = sum(trade_0, dims=2)
    tot_exports_1 = sum(trade_1, dims=2)

    # Total imports
    tot_imports_0 = transpose(sum(trade_0, dims=1))
    tot_imports_1 = transpose(sum(trade_1, dims=1))

    # Total international trade
    tot_trade_0 = tot_exports_0 .+ tot_imports_0
    tot_trade_1 = tot_exports_1 .+ tot_imports_1

    # Percent changes
    g_exports = 100 .* tot_exports_1[:] ./ tot_exports_0[:] .- 100
    g_imports = 100 .* tot_imports_1[:] ./ tot_imports_0[:] .- 100
    g_trade = 100 .* tot_trade_1[:] ./ tot_trade_0[:] .- 100
    g_domestic = 100 .* domestic_1 ./ domestic_0 .- 100
    welfare = 100 .* (p_hat[:] ./ P_hat[:]).^(1+ψ) .- 100
    output = 100 .* (p_hat[:] ./ P_hat[:]).^(ψ) .- 100

    df = DataFrame(iso=iso, Trade=g_trade, Output=output, Welfare=welfare, Exports=g_exports, Imports=g_imports, Domestic=g_domestic)
    return df
end

"""
Convert dict {bloc: array of ISOs} to dict {bloc: indices of countries in an array of ISOs}
"""
function iso2idx(blocs, iso)
    blocs_idx = Dict()
    for (k, v) in blocs
        idx = findfirst.(isequal.(v), (iso,))
        blocs_idx[k] = idx
    end
    return blocs_idx
end

"""
Summarize trade flows between countries and blocs
"""
function summarize_blocs(iso, X, X_1)
    # Blocs
    blocs = load_bloc_definitions()
    blocs_idx = iso2idx(blocs, iso)
    bloc_partition = [v for v in values(blocs_idx)]

    # International trade
    N = size(X, 1)
    trade_0 = (ones(N, N) - I) .* X
    trade_1 = (ones(N, N) - I) .* X_1

    # Aggregate blocs
    # X_blocs_0 = bloc_trade_matrix(trade_0, bloc_partition, bloc_partition)
    # X_blocs_1 = bloc_trade_matrix(trade_1, bloc_partition, bloc_partition)
    # g_blocs = 100 .* X_blocs_1 ./ X_blocs_0 .- 100

    # Blocs with detailed countries as exporters
    country_partition = [[c] for c in 1:N]
    X_blocs_0_x = bloc_trade_matrix(trade_0, country_partition, bloc_partition)
    X_blocs_1_x = bloc_trade_matrix(trade_1, country_partition, bloc_partition)

    # Blocs with detailed countries as importers
    X_blocs_0_m = bloc_trade_matrix(trade_0, bloc_partition, country_partition)
    X_blocs_1_m = bloc_trade_matrix(trade_1, bloc_partition, country_partition)

    # Average exports and imports
    X_blocs_0 = (X_blocs_0_x .+ transpose(X_blocs_0_m)) ./ 2
    X_blocs_1 = (X_blocs_1_x .+ transpose(X_blocs_1_m)) ./ 2
    g_blocs = 100 .* X_blocs_1 ./ X_blocs_0 .- 100

    # Construct dataframe
    rows = iso
    cols = [string(k) for k in keys(blocs)]
    df = DataFrame(g_blocs, cols)
    insertcols!(df, 1, "iso" => rows)

    return df
end

"""
Summarize trade flows between blocs
"""
function summarize_aggregate_blocs(iso, X, X_1)
    # Blocs
    blocs = load_bloc_definitions()
    blocs_idx = iso2idx(blocs, iso)
    bloc_partition = [v for v in values(blocs_idx)]

    # International trade
    N = size(X, 1)
    trade_0 = (ones(N, N) - I) .* X
    trade_1 = (ones(N, N) - I) .* X_1

    # Aggregate blocs
    X_blocs_0 = bloc_trade_matrix(trade_0, bloc_partition, bloc_partition)
    X_blocs_1 = bloc_trade_matrix(trade_1, bloc_partition, bloc_partition)
    g_blocs = 100 .* X_blocs_1 ./ X_blocs_0 .- 100

    # Construct dataframe
    rows = [string(k) for k in keys(blocs)]
    cols = [string(k) for k in keys(blocs)]
    df = DataFrame(g_blocs, cols)
    insertcols!(df, 1, "Bloc" => rows)

    return df
end

# Choose parameters

# There are different options for the supply elasticity ψ, e.g.,
# ψ = 3.76  (Easton and Kortum, 2002)
# ψ = 1.0  (Alvarez and Lucas, 2007)

# We calculate the supply elasticity as the average of the elasticities implied by
# the 10-90 percentiles of intermediate shares in KLEMS reported in Huo, Levchenko, Pandalai-Nayar (2021)
η_10 = 0.31
η_90 = 0.67
ψ_10 = η_10/(1-η_10)
ψ_90 = η_90/(1-η_90)
ψ = ψ_10/2 + ψ_90/2  # 1.24

ϕ = 3.8  # Trade elasticity taken from the meta analysis for the Armington elasticity

# Load bloc information
filepath = joinpath(@__DIR__, "source_data" , "def_blocs.csv")
df_blocs = DataFrame(CSV.File(filepath))
df_blocs = select(df_blocs, :iso, :bloc)  # Note: df[!, :iso] is not necessarily ordered in the same way as the iso array produced by load_X() and simulate(...)
rename!(df_blocs, :bloc => :Bloc)

# Scenario: 0 (Russia in autarky)

"""
Calculate domestic expenditure shares
"""
function domestic_expenditure_share(iso)
    # Benchmark data
    iso_a, X_a = load_X("trade_2019.csv")
    idx_a = findfirst(isequal(iso), iso_a)
    # TiVA data
    iso_b, X_b = load_X("trade_prod_2018.csv")
    idx_b = findfirst(isequal(iso), iso_b)

    # Calculate with benchmark data
    λ_a = X_a[:, idx_a] ./ sum(X_a[:, idx_a])
    # Calculate with TiVA data
    λ_b = X_b[:, idx_b] ./ sum(X_b[:, idx_b])

    return λ_a[idx_a], λ_b[idx_b]
end

"""
Autarky welfare loss with ACR formula
"""
function welfare_autarky(λ, ϕ)
    W_hat = λ^(1/ϕ)
    return 100 * (W_hat - 1)
end

"""
Autarky welfare loss with positive supply
"""
function welfare_autarky(λ, ϕ, ψ)
    W_hat = λ^((1+ψ)/ϕ)
    return 100 * (W_hat - 1)
end



# Obtain welfare for RUS
λ_RUS = domestic_expenditure_share("RUS")
W_ACR = welfare_autarky.(λ_RUS, ϕ)
W_model = welfare_autarky.(λ_RUS, ϕ, ψ)

# Now for the rest of the world
iso, X, X_1, p_hat, P_hat = simulate(0, ϕ, ψ)
df_0 = summarize(iso, X, X_1, p_hat, P_hat, ψ)  # all but domestic trade for RUS is correct

idx_RUS = findfirst(isequal("RUS"), iso)
abs(df_0[idx_RUS, :Trade] - (-100.0))
@assert df_0[idx_RUS, :Welfare] ≈ W_model[1]  # Check that welfare for RUS is correct

# Domestic trade for Russia is incorrect because the country is now disconnected from world prices
# To obtain it, we redo the exercise by setting Russia's price as the numeraire
# Note, however, that results for international/domestic trade for all other countries will be incorrect in this exercise because of the same disconnect in prices
_, _, X_1_RUS, p_hat_RUS, P_hat_RUS = simulate(0, ϕ, ψ; numeraire=idx_RUS)
df_0_RUS = summarize(iso, X, X_1_RUS, p_hat_RUS, P_hat_RUS, ψ)  # only domestic trade for RUS is correct
# Check that relative prices agree
@assert p_hat ./ P_hat ≈ p_hat_RUS ./ P_hat_RUS
# If relative prices agree, then output and welfare must also agree for all countries
@assert Array(df_0[:, [:Output, :Welfare]]) ≈ Array(df_0_RUS[:, [:Output, :Welfare]])

# Now replace domestic trade for RUS in the counterfactual trade matrix...
X_1[idx_RUS, idx_RUS] = X_1_RUS[idx_RUS, idx_RUS]
# ...to obtain correct summary data
df_0 = summarize(iso, X, X_1, p_hat, P_hat, ψ)


df_0 = select(df_0, [:iso, :Trade, :Output, :Welfare])

filepath = joinpath(@__DIR__, "results" , "fragmentation_0.csv")
CSV.write(filepath, df_0)

# Scenario: 1 (MATR increases between opposing blocs)

iso, X, X_1, p_hat, P_hat = simulate(1, ϕ, ψ)

df_1 = summarize(iso, X, X_1, p_hat, P_hat, ψ)
df_1b = summarize_blocs(iso, X, X_1)
df_1 = innerjoin(df_1b, df_1, on=:iso)
df_1 = innerjoin(df_blocs, df_1, on=:iso)
df_1 = select(df_1, [:iso, :Bloc, :Eastern, :Western, :Neutral, :Trade, :Output, :Welfare])

filepath = joinpath(@__DIR__, "results" , "fragmentation_1.csv")
CSV.write(filepath, df_1)

df_blocs_1 = summarize_aggregate_blocs(iso, X, X_1)

filepath = joinpath(@__DIR__, "results" , "fragmentation_blocs_1.csv")
CSV.write(filepath, df_blocs_1)

# Scenario: 3 (Scenario 1 + WTO exit of the Eastern bloc + no TA between opposing blocs)
iso, X, X_1, p_hat, P_hat = simulate(3, ϕ, ψ)

df_3 = summarize(iso, X, X_1, p_hat, P_hat, ψ)
df_3b = summarize_blocs(iso, X, X_1)
df_3 = innerjoin(df_3b, df_3, on=:iso)
df_3 = innerjoin(df_blocs, df_3, on=:iso)
df_3 = select(df_3, [:iso, :Bloc, :Eastern, :Western, :Neutral, :Trade, :Output, :Welfare])
filepath = joinpath(@__DIR__, "results" , "fragmentation_3.csv")
CSV.write(filepath, df_3)

df_blocs_3 = summarize_aggregate_blocs(iso, X, X_1)
filepath = joinpath(@__DIR__, "results" , "fragmentation_blocs_3.csv")
CSV.write(filepath, df_blocs_3)

# Robustness checks

# Scenario 1 with PROD data (scenario=2)

iso, X, X_1, p_hat, P_hat = simulate(2, ϕ, ψ)

df_2 = summarize(iso, X, X_1, p_hat, P_hat, ψ)
df_2 = select(df_2, [:iso, :Trade, :Output, :Welfare])

rename!(df_2, :Trade => :Trade_prod, :Output => :Output_prod, :Welfare => :Welfare_prod)
df_1_comp = innerjoin(df_1, df_2, on=:iso)

filepath = joinpath(@__DIR__, "results" , "fragmentation_1_prod.csv")
CSV.write(filepath, df_1_comp)

# Scenario 3 with PROD data (scenario=4)

iso, X, X_1, p_hat, P_hat = simulate(4, ϕ, ψ)

df_4 = summarize(iso, X, X_1, p_hat, P_hat, ψ)
df_4 = select(df_4, [:iso, :Trade, :Output, :Welfare])

rename!(df_4, :Trade => :Trade_prod, :Output => :Output_prod, :Welfare => :Welfare_prod)
df_3_comp = innerjoin(df_3, df_4, on=:iso)

filepath = joinpath(@__DIR__, "results" , "fragmentation_3_prod.csv")
CSV.write(filepath, df_3_comp)
