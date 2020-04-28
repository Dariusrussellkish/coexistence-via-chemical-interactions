## Well-mixed model for growth of interacting species
# UIC: Uniform initial condition
# Ex: explicitly including the mediators
# MT: multi-target mediators
# ExMT2: corrected the error in ExMT, now rIntMat only includes links in R
using Random;
import Future;
using SharedArrays
using Distributed
using Serialization
using JLD
using JSON
using UUIDs
@everywhere using Logging
@everywhere using Distributions
@everywhere include("DistInteractionStrengthMT_PB.jl")
@everywhere include("NetworkConfig_Binomial.jl")
@everywhere include("WellmixedInteraction_DpMM_ExMT4.jl")

# Main harness function
function main(
    nSample,
    nGen,
    nInitialCell,
    dilTh,
    nCellType,
    nMediator,
    kSatLevel,
    extTh,
    ri0,
    posIntRatio,
    tauf,
    dtau,
    at,
    bt,
    qp,
    qc,
    r,
)
    GenPerRound = log(dilTh / nInitialCell) / log(2)
    nRound = Integer(round(nGen / GenPerRound)) # number of rounds of propagation

    tau0 = 0 # in hours
    dirichlet = Dirichlet(ones(nCellType))

    r0 = 0.08 .+ 0.04 .* rand(r[myid()], nCellType, 1) # population reproduction rates, per hour
    kSatVector = kSatLevel * (0.5 .+ rand(r[myid()], nMediator, 1)) # population levels for influence saturation

    ## Parameters
    # Network configuration
    # NetworkConfig_Balanced(Nc,Nm,q): link between Nc and Nm present with
    # a probability q
    R = NetworkConfig_Binomial(nCellType, nMediator, qc, r)
    P = NetworkConfig_Binomial(nCellType, nMediator, qp, r)

    # interaction matrix
    alpha = at .* (0.5 .+ rand(r[myid()], nCellType, nMediator)) # consumption rates 0.5at to 1.5at
    beta = bt .* (0.5 .+ rand(r[myid()], nCellType, nMediator)) # production rates 05.bt to 1.5bt
    A = (R .* alpha)' # nc x nm function was defined as nm x nc
    B = (P .* beta)' # ^

    rIntMatA =
        R .* DistInteractionStrengthMT_PB(nCellType, nMediator, ri0, posIntRatio, r) # matrix of interaction coefficients, 50/50

    CMPs = SharedArray{Float64, 2}((nCellType, nSample), init=zeros(nCellType, nSample))
    @sync @distributed for ns = 1:nSample

        cellRatioArray = reshape(rand(r[myid()], dirichlet), (1, :)) # cell distribution population ratios

        NeAD, CmpAD = WellmixedInteraction_DpMM_ExMT4(
            nRound,
            r0,
            deepcopy(cellRatioArray),
            rIntMatA,
            nInitialCell,
            kSatVector,
            A,
            B,
            kSatLevel,
            extTh,
            dilTh,
            tauf,
            dtau,
        )

        Cmp0AD = zeros(1, nCellType)
        Cmp0AD[NeAD] = CmpAD
        CMPs[:, ns] = Cmp0AD

        @info "pid: $(myid()) ns: $ns, Comp: $Cmp0AD"
    end
    return CMPs
end

# Start of "script" portion
@info "Starting Simulation with $(nworkers()) worker$(nworkers() > 1 ? "s" : "")"

# get a starting seed
rndseed0 = convert(Int32, trunc(time()))

# generate sufficient independent RNGs
r = let m = MersenneTwister(rndseed0)
    [m; accumulate(Future.randjump, fill(big(10)^20, nworkers() + 1), init = m)]
end

# read in config file
open(ARGS[1], "r") do io
    params = JSON.parse(io)

    nSample = params["nSample"]
    nGen = params["nGen"]
    nInitialCell = params["nInitialCell"] # total initial cells
    dilTh = params["dilTh"] # coculture dilution threshold
    nCellType = params["nCellType"] # # of cell types in the initial pool
    nMediator = params["nMediator"] # # of mediators
    kSatLevel = params["kSatLevel"] # interaction strength saturation level of each population
    extTh = params["extTh"] # population extinction threshold
    ri0 = params["ri0"] # maximum interaction strength, 1/hr
    posIntRatio = params["posIntRatio"] # fraction of interactions that are positive
    tau0 = params["tau0"] # in hours
    tauf = params["tauf"] # in hours
    dtau = params["dtau"] # in hours, cell growth update and uptake timescale
    at = params["at"] # avg. consumption values (fmole per cell) alpha_ij: population i, resource j
    bt = params["bt"] # avg. production rates (fmole per cell per hour) beta_ij: population i, resource j
    qp = params["qp"] # probability of production link per population
    qc = params["qc"] # probability of influence link per population

    # filename generation
    stripChar = (s, r) -> replace(s, Regex("[$r]") => "")
    filename_uuid = stripChar(string(uuid4()), "-")
    basename = "nGen_$(nGen)_nCellType_$(nCellType)_nMediator_$(nMediator)_ri0_$(ri0)_posIntRatio_$(posIntRatio)_at_$(at)_bt_$(bt)_qp_$(qp)_qc_$(qc)_seed_$(rndseed0)_$(filename_uuid)"

    # harness call
    CMPs = main(
        nSample,
        nGen,
        nInitialCell,
        dilTh,
        nCellType,
        nMediator,
        kSatLevel,
        extTh,
        ri0,
        posIntRatio,
        tauf,
        dtau,
        at,
        bt,
        qp,
        qc,
        r,
    )

    # serialize
    save(
        "$(basename).jld",
        "CMPs",
        CMPs,
        "params",
        params,
    )
end
