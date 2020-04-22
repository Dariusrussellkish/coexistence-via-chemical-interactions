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
    nRound = round(nGen / GenPerRound) # number of rounds of propagation

    tau0 = 0 # in hours

    # T – record keeping over replicates
    r0T = SharedArray{Float64,2}(
        (nCellType, nSample),
        init = zeros(nCellType, nSample),
    )
    #r0T = zeros(nCellType, nSample) #matrix of zeros: 4 rows because 4 species and 1000 columns because 1000 samples being screened
    SiT = SharedArray{Float64,2}(
        (nMediator, nSample),
        init = zeros(nMediator, nSample),
    )
    #SiT = zeros(nMediator, nSample)
    initCellRatioArray = SharedArray{Float64,2}(
        (nCellType, nSample),
        init = zeros(nCellType, nSample),
    )
    #initCellRatioArray = zeros(nCellType, nSample)

    AT = SharedArray{Float64,3}(
        (nMediator, nCellType, nSample),
        init = zeros(nMediator, nCellType, nSample),
    )
    #AT = zeros(nMediator, nCellType, nSample)
    BT = SharedArray{Float64,3}(
        (nMediator, nCellType, nSample),
        init = zeros(nMediator, nCellType, nSample),
    )
    #BT = zeros(nMediator, nCellType, nSample)
    DT = SharedArray{Float64,3}(
        (nMediator, nCellType, nSample),
        init = zeros(nMediator, nCellType, nSample),
    )
    #DT = zeros(nMediator, nSample)

    CmpADT = SharedArray{Float64,2}(
        (nCellType, nSample),
        init = zeros(nCellType, nSample),
    )
    #CmpADT = zeros(nCellType, nSample)
    CmpBDT = SharedArray{Float64,2}(
        (nCellType, nSample),
        init = zeros(nCellType, nSample),
    )
    #CmpBDT = zeros(nCellType, nSample)
    CmpEDT = SharedArray{Float64,2}(
        (nCellType, nSample),
        init = zeros(nCellType, nSample),
    )
    #CmpEDT = zeros(nCellType, nSample)

    rintAT = SharedArray{Float64,3}(
        (nMediator, nCellType, nSample),
        init = zeros(nMediator, nCellType, nSample),
    )
    rintBT = SharedArray{Float64,3}(
        (nMediator, nCellType, nSample),
        init = zeros(nMediator, nCellType, nSample),
    )
    rintET = SharedArray{Float64,3}(
        (nMediator, nCellType, nSample),
        init = zeros(nMediator, nCellType, nSample),
    )
    # rintAT = zeros(nMediator, nCellType, nSample)
    # rintBT = zeros(nMediator, nCellType, nSample)
    # rintET = zeros(nMediator, nCellType, nSample)

    V0DT = SharedArray{Float64,3}(
        (3, nCellType, nSample),
        init = zeros(3, nCellType, nSample),
    )
    VDT = SharedArray{Float64,3}(
        (3, nCellType, nSample),
        init = zeros(3, nCellType, nSample),
    )
    # V0DT = zeros(3, nCellType, nSample)
    # VDT = zeros(3, nCellType, nSample)

    # primary variable of interest, tracks number of existing species at
    #end of simulation
    NE0D = SharedArray{Float64,2}((3, nSample), init = zeros(3, nSample))

    @sync @distributed for ns = 1:nSample

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
            R .* DistInteractionStrengthMT_PB(nCellType, nMediator, ri0, 0.5, r) # matrix of interaction coefficients, 50/50
        rIntMatB =
            R .* DistInteractionStrengthMT_PB(
                nCellType,
                nMediator,
                ri0,
                posIntRatio,
                r,
            ) # matrix of interaction coefficients, more negative
        rIntMatE =
            R .* DistInteractionStrengthMT_PB(
                nCellType,
                nMediator,
                ri0,
                1 - posIntRatio,
                r,
            ) # matrix of interaction coefficients, more positive

        # TODO: Change ones distribution to random distribution
        cellRatioArray = 1 / nCellType .* ones(1, nCellType) # cell distribution population ratios


        ## Simulating dynamics, Dp, depletable
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
        NeBD, CmpBD = WellmixedInteraction_DpMM_ExMT4(
            nRound,
            r0,
            deepcopy(cellRatioArray),
            rIntMatB,
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
        NeED, CmpED = WellmixedInteraction_DpMM_ExMT4(
            nRound,
            r0,
            deepcopy(cellRatioArray),
            rIntMatE,
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

        # fine grained tracking of surviving species
        V0AD = zeros(1, nCellType)
        V0BD = zeros(1, nCellType)
        V0ED = zeros(1, nCellType)
        V0AD[NeAD] .= 1
        V0BD[NeBD] .= 1
        V0ED[NeED] .= 1

        @info "pid: $(myid()) ns: $ns, NeED: $NeED, NeAD: $NeAD, NeBD: $NeBD"

        NE0D[:, ns] = sum([V0AD; V0BD; V0ED], dims = 2)

        Cmp0AD = zeros(1, nCellType)
        Cmp0BD = zeros(1, nCellType)
        Cmp0ED = zeros(1, nCellType)
        Cmp0AD[NeAD] = CmpAD
        Cmp0BD[NeBD] = CmpBD
        Cmp0ED[NeED] = CmpED

        CmpADT[:, ns] = Cmp0AD
        CmpBDT[:, ns] = Cmp0BD
        CmpEDT[:, ns] = Cmp0ED

        r0T[:, ns] = r0
        SiT[:, ns] = kSatVector
        AT[:, :, ns] = A
        BT[:, :, ns] = B
        rintAT[:, :, ns] = rIntMatA'
        rintBT[:, :, ns] = rIntMatB'
        rintET[:, :, ns] = rIntMatE'
        V0DT[:, :, ns] = [V0AD; V0BD; V0ED]
        initCellRatioArray[:, ns] = cellRatioArray
    end
    return NE0D,
    CmpADT,
    CmpBDT,
    CmpEDT,
    r0T,
    SiT,
    AT,
    BT,
    rintAT,
    rintBT,
    rintET,
    V0DT,
    initCellRatioArray
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
    NE0D,
    CmpADT,
    CmpBDT,
    CmpEDT,
    r0T,
    SiT,
    AT,
    BT,
    rintAT,
    rintBT,
    rintET,
    V0DT,
    initCellRatioArray = main(
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
        "NE0D",
        NE0D,
        "CmpADT",
        CmpADT,
        "CmpBDT",
        CmpBDT,
        "CmpEDT",
        CmpEDT,
        "r0T",
        r0T,
        "SiT",
        SiT,
        "AT",
        AT,
        "BT",
        BT,
        "rintAT",
        rintAT,
        "rintBT",
        rintBT,
        "rintET",
        rintET,
        "V0DT",
        V0DT,
        "params",
        params,
    )

    println(maximum(NE0D))
end
