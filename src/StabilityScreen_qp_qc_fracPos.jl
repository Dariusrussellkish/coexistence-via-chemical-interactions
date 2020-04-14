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
@everywhere using Logging
@everywhere include("DistInteractionStrengthMT_PB.jl")
@everywhere include("NetworkConfig_Binomial.jl")
@everywhere include("WellmixedInteraction_DpMM_ExMT4.jl")

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
    coarseness,
    minFracPos,
    maxFracPos,
    deltaFracPos,
)
    GenPerRound = log(dilTh / nInitialCell) / log(2)
    nRound = round(nGen / GenPerRound) # number of rounds of propagation

    tau0 = 0 # in hours

    fracPosArray = collect(minFracPos:deltaFracPos:maxFracPos)
    nFracPos = size(fracPosArray)[1]

    # T – record keeping over replicates
    r0T = SharedArray{Float64,5}(
        (nFracPos, coarseness, coarseness, nCellType, nSample),
        init = zeros(nFracPos, coarseness, coarseness, nCellType, nSample),
    )
    SiT = SharedArray{Float64,5}(
        (nFracPos, coarseness, coarseness, nMediator, nSample),
        init = zeros(nFracPos, coarseness, coarseness, nMediator, nSample),
    )
    initCellRatioArray = SharedArray{Float64,5}(
        (nFracPos, coarseness, coarseness, nCellType, nSample),
        init = zeros(nFracPos, coarseness, coarseness, nCellType, nSample),
    )
    AT = SharedArray{Float64,6}(
        (nFracPos, coarseness, coarseness, nMediator, nCellType, nSample),
        init = zeros(
            nFracPos,
            coarseness,
            coarseness,
            nMediator,
            nCellType,
            nSample,
        ),
    )
    BT = SharedArray{Float64,6}(
        (nFracPos, coarseness, coarseness, nMediator, nCellType, nSample),
        init = zeros(
            nFracPos,
            coarseness,
            coarseness,
            nMediator,
            nCellType,
            nSample,
        ),
    )
    CmpT = SharedArray{Float64,5}(
        (nFracPos, coarseness, coarseness, nCellType, nSample),
        init = zeros(nFracPos, coarseness, coarseness, nCellType, nSample),
    )
    rintT = SharedArray{Float64,6}(
        (nFracPos, coarseness, coarseness, nMediator, nCellType, nSample),
        init = zeros(
            nFracPos,
            coarseness,
            coarseness,
            nMediator,
            nCellType,
            nSample,
        ),
    )
    V0DT = SharedArray{Float64,6}(
        (nFracPos, coarseness, coarseness, 1, nCellType, nSample),
        init = zeros(nFracPos, coarseness, coarseness, 1, nCellType, nSample),
    )
    VDT = SharedArray{Float64,6}(
        (nFracPos, coarseness, coarseness, 1, nCellType, nSample),
        init = zeros(nFracPos, coarseness, coarseness, 1, nCellType, nSample),
    )
    # V0DT = zeros(3, nCellType, nSample)
    # VDT = zeros(3, nCellType, nSample)

    NE0D = SharedArray{Float64,5}(
        (nFracPos, coarseness, coarseness, 1, nSample),
        init = zeros(nFracPos, coarseness, coarseness, 1, nSample),
    )
    #NE0D = zeros(3, nSample)
    for qp = 1:coarseness, qc = 1:coarseness, fpi = 1:nFracPos
        qpv = (qp - 1) / (coarseness - 1)
        qcv = (qc - 1) / (coarseness - 1)
        @sync @distributed for ns = 1:nSample
            # @time begin
            #println("Sample: $ns")
            # tic
            posIntRatio = fracPosArray[fpi]
            # rand('twister',rndseed(ns)) # prevents RNG loop?
            r0 = 0.08 .+ 0.04 .* rand(r[myid()], nCellType, 1) # population reproduction rates, per hour
            kSatVector = kSatLevel * (0.5 .+ rand(r[myid()], nMediator, 1)) # population levels for influence saturation

            ## Parameters
            # Network configuration
            # NetworkConfig_Balanced(Nc,Nm,q): link between Nc and Nm present with
            # a probability q
            R = NetworkConfig_Binomial(nCellType, nMediator, qcv, r)
            P = NetworkConfig_Binomial(nCellType, nMediator, qpv, r)

            # interaction matrix
            alpha = at .* (0.5 .+ rand(r[myid()], nCellType, nMediator)) # consumption rates 0.5at to 1.5at
            beta = bt .* (0.5 .+ rand(r[myid()], nCellType, nMediator)) # production rates 05.bt to 1.5bt
            A = (R .* alpha)' # nc x nm function was defined as nm x nc
            B = (P .* beta)' # ^

            rIntMatA =
                R .* DistInteractionStrengthMT_PB(
                    nCellType,
                    nMediator,
                    ri0,
                    posIntRatio,
                    r,
                ) # matrix of interaction coefficients

            cellRatioArray = 1 / nCellType .* ones(1, nCellType) # cell distribution population ratios

            # @info "pid: $(myid()) cellRatioArray: $cellRatioArray"

            ## Simulating dynamics, Dp, depletable
            # @info "pid: $(myid()) ns: $ns"
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
            # @info "pid: $(myid()) ns: $ns, NeED: $NeED, cellRatioArray: $cellRatioArray"

            V0AD = zeros(1, nCellType)
            V0AD[NeAD] .= 1

            @info "pid: $(myid()) qp: $qpv qc: $qcv fp: $posIntRatio ns: $ns, NeAD: $NeAD"

            #println(NeED)
            NE0D[fpi, qp, qc, :, ns] .= sum(V0AD)

            Cmp0AD = zeros(1, nCellType)
            Cmp0AD[NeAD] = CmpAD

            CmpT[fpi, qp, qc, :, ns] = Cmp0AD

            r0T[fpi, qp, qc, :, ns] = r0
            SiT[fpi, qp, qc, :, ns] = kSatVector
            AT[fpi, qp, qc, :, :, ns] = A
            BT[fpi, qp, qc, :, :, ns] = B
            rintT[fpi, qp, qc, :, :, ns] = rIntMatA'
            V0DT[fpi, qp, qc, :, :, ns] .= V0AD
            initCellRatioArray[fpi, qp, qc, :, ns] = cellRatioArray
        end
    end
    return NE0D,
    CmpT,
    r0T,
    SiT,
    AT,
    BT,
    rintT,
    V0DT,
    initCellRatioArray
end

@info "Starting Simulation with $(nworkers()) worker$(nworkers() > 1 ? "s" : "")"

rndseed0 = convert(Int32, trunc(time()))

r = let m = MersenneTwister(rndseed0)
    [m; accumulate(Future.randjump, fill(big(10)^20, nworkers() + 1), init = m)]
end

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
    coarseness = params["coarseness"]
    minFracPos = params["minFracPos"]
    maxFracPos = params["maxFracPos"]
    deltaFracPos = params["deltaFracPos"]

    basename = "nGen_$(nGen)_nCellType_$(nCellType)_nMediator_$(nMediator)_ri0_$(ri0)_posIntRatio_$(posIntRatio)_at_$(at)_bt_$(bt)_CS3390_HEATMAP_seed_$(rndseed0)"

    NE0D,
    CmpT,
    r0T,
    SiT,
    AT,
    BT,
    rintT,
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
        coarseness + 1,
        minFracPos,
        maxFracPos,
        deltaFracPos,
    )

    save(
        "$(basename).jld",
        "NE0D",
        NE0D,
        "CmpT",
        CmpT,
        "r0T",
        r0T,
        "SiT",
        SiT,
        "AT",
        AT,
        "BT",
        BT,
        "rintT",
        rintT,
        "V0DT",
        V0DT,
        "params",
        params,
    )

    println(maximum(NE0D))
end
# serialize("$basename.save",
# (
# NE0D = NE0D,
# CmpADT = CmpADT,
# CmpBDT = CmpBDT,
# CmpEDT = CmpEDT,
# r0T = r0T,
# SiT = SiT,
# AT = AT,
# BT = BT,
# rintAT = rintAT,
# rintBT = rintBT,
# rintET = rintET,
# V0DT = V0DT,
# ))
