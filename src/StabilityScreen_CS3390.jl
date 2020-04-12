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
)
    GenPerRound = log(dilTh / nInitialCell) / log(2)
    nRound = round(nGen / GenPerRound) # number of rounds of propagation

    tau0 = 0 # in hours

    # T – record keeping over replicates
    r0T = SharedArray{Float64,4}(
        (coarseness, coarseness, nCellType, nSample),
        init = zeros(coarseness, coarseness, nCellType, nSample),
    )
    #r0T = zeros(nCellType, nSample) #matrix of zeros: 4 rows because 4 species and 1000 columns because 1000 samples being screened
    SiT = SharedArray{Float64,4}(
        (coarseness, coarseness, nMediator, nSample),
        init = zeros(coarseness, coarseness, nMediator, nSample),
    )
    #SiT = zeros(nMediator, nSample)
    initCellRatioArray = SharedArray{Float64,4}(
        (coarseness, coarseness, nCellType, nSample),
        init = zeros(coarseness, coarseness, nCellType, nSample),
    )
    #initCellRatioArray = zeros(nCellType, nSample)

    AT = SharedArray{Float64,5}(
        (coarseness, coarseness, nMediator, nCellType, nSample),
        init = zeros(coarseness, coarseness, nMediator, nCellType, nSample),
    )
    #AT = zeros(nMediator, nCellType, nSample)
    BT = SharedArray{Float64,5}(
        (coarseness, coarseness, nMediator, nCellType, nSample),
        init = zeros(coarseness, coarseness, nMediator, nCellType, nSample),
    )
    #BT = zeros(nMediator, nCellType, nSample)
    DT = SharedArray{Float64,5}(
        (coarseness, coarseness, nMediator, nCellType, nSample),
        init = zeros(coarseness, coarseness, nMediator, nCellType, nSample),
    )
    #DT = zeros(nMediator, nSample)

    CmpADT = SharedArray{Float64,4}(
        (coarseness, coarseness, nCellType, nSample),
        init = zeros(coarseness, coarseness, nCellType, nSample),
    )
    #CmpADT = zeros(nCellType, nSample)
    CmpBDT = SharedArray{Float64,4}(
        (coarseness, coarseness, nCellType, nSample),
        init = zeros(coarseness, coarseness, nCellType, nSample),
    )
    #CmpBDT = zeros(nCellType, nSample)
    CmpEDT = SharedArray{Float64,4}(
        (coarseness, coarseness, nCellType, nSample),
        init = zeros(coarseness, coarseness, nCellType, nSample),
    )
    #CmpEDT = zeros(nCellType, nSample)

    rintAT = SharedArray{Float64,5}(
        (coarseness, coarseness, nMediator, nCellType, nSample),
        init = zeros(coarseness, coarseness, nMediator, nCellType, nSample),
    )
    rintBT = SharedArray{Float64,5}(
        (coarseness, coarseness, nMediator, nCellType, nSample),
        init = zeros(coarseness, coarseness, nMediator, nCellType, nSample),
    )
    rintET = SharedArray{Float64,5}(
        (coarseness, coarseness, nMediator, nCellType, nSample),
        init = zeros(coarseness, coarseness, nMediator, nCellType, nSample),
    )
    # rintAT = zeros(nMediator, nCellType, nSample)
    # rintBT = zeros(nMediator, nCellType, nSample)
    # rintET = zeros(nMediator, nCellType, nSample)

    V0DT = SharedArray{Float64,5}(
        (coarseness, coarseness, 3, nCellType, nSample),
        init = zeros(coarseness, coarseness, 3, nCellType, nSample),
    )
    VDT = SharedArray{Float64,5}(
        (coarseness, coarseness, 3, nCellType, nSample),
        init = zeros(coarseness, coarseness, 3, nCellType, nSample),
    )
    # V0DT = zeros(3, nCellType, nSample)
    # VDT = zeros(3, nCellType, nSample)

    NE0D = SharedArray{Float64,4}((coarseness, coarseness, 3, nSample), init = zeros(coarseness, coarseness, 3, nSample))
    #NE0D = zeros(3, nSample)
    for qp=1:coarseness, qc=1:coarseness
        qpv = qp/coarseness
        qcv = qc/coarseness
        @sync @distributed for ns = 1:nSample
            # @time begin
            #println("Sample: $ns")
            # tic

            # rand('twister',rndseed(ns)) # prevents RNG loop?
            r0 = 0.08 .+ 0.04 .* rand(r[myid()], nCellType, 1) # population reproduction rates, per hour
            kSatVector = kSatLevel * (0.5 .+ rand(r[myid()], nMediator, 1)) # population levels for influence saturation

            ## Parameters
            # Network configuration
            # NetworkConfig_Balanced(Nc,Nm,q): link between Nc and Nm present with
            # a probability q
            R = NetworkConfig_Binomial(nCellType, nMediator, qcv, r)
            P = NetworkConfig_Binomial(nCellType, nMediator, qpv, r)

            # display(R)
            # println("\n")
            # display(P)
            # println("\n")

            # interaction matrix
            alpha = at .* (0.5 .+ rand(r[myid()], nCellType, nMediator)) # consumption rates 0.5at to 1.5at
            beta = bt .* (0.5 .+ rand(r[myid()], nCellType, nMediator)) # production rates 05.bt to 1.5bt
            A = (R .* alpha)' # nc x nm function was defined as nm x nc
            B = (P .* beta)' # ^
            # display(A)
            # println("\n")
            # display(B)
            # println("\n")
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
            # display(map(x -> round(x, digits=2), rIntMatA))
            # println("\n")
            # display(map(x -> round(x, digits=2), rIntMatB))
            # println("\n")
            # display(map(x -> round(x, digits=2), rIntMatE))
            # println("\n")
            # TODO: Change ones distribution to random distribution
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
            # @info "pid: $(myid()) ns: $ns, NeAD: $NeAD"
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
            # @info "pid: $(myid()) ns: $ns, NeBD: $NeBD"
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
            # @info "pid: $(myid()) ns: $ns, NeED: $NeED, cellRatioArray: $cellRatioArray"

            V0AD = zeros(1, nCellType)
            V0BD = zeros(1, nCellType)
            V0ED = zeros(1, nCellType)
            V0AD[NeAD] .= 1
            V0BD[NeBD] .= 1
            V0ED[NeED] .= 1

            @info "pid: $(myid()) qp: $qpv qc: $qcv ns: $ns, NeED: $NeED, NeAD: $NeAD, NeBD: $NeBD"

            #println(NeED)
            NE0D[qp, qc, :, ns] = sum([V0AD; V0BD; V0ED], dims = 2)

            Cmp0AD = zeros(1, nCellType)
            Cmp0BD = zeros(1, nCellType)
            Cmp0ED = zeros(1, nCellType)
            Cmp0AD[NeAD] = CmpAD
            Cmp0BD[NeBD] = CmpBD
            Cmp0ED[NeED] = CmpED

            CmpADT[qp, qc, :, ns] = Cmp0AD
            CmpBDT[qp, qc, :, ns] = Cmp0BD
            CmpEDT[qp, qc, :, ns] = Cmp0ED

            r0T[qp, qc, :, ns] = r0
            SiT[qp, qc, :, ns] = kSatVector
            AT[qp, qc, :, :, ns] = A
            BT[qp, qc, :, :, ns] = B
            rintAT[qp, qc, :, :, ns] = rIntMatA'
            rintBT[qp, qc, :, :, ns] = rIntMatB'
            rintET[qp, qc, :, :, ns] = rIntMatE'
            V0DT[qp, qc, :, :, ns] = [V0AD; V0BD; V0ED]
            initCellRatioArray[qp, qc, :, ns] = cellRatioArray
        end
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

@info "Starting Simulation with $(nworkers()) worker$(nworkers() > 1 ? "s" : "")"

rndseed0 = convert(Int32, trunc(time()))

r = let m = MersenneTwister(rndseed0)
    [m; accumulate(Future.randjump, fill(big(10)^20, nworkers() + 1), init = m)]
end

nGen = 200
nInitialCell = 1e4 # total initial cells
dilTh = 1e10 # coculture dilution threshold
nCellType = 20 # # of cell types in the initial pool
nMediator = 15 # # of mediators
kSatLevel = 1e4 # interaction strength saturation level of each population
extTh = 0.1 # population extinction threshold
ri0 = 0.2 # maximum interaction strength, 1/hr
posIntRatio = 0.1 # fraction of interactions that are positive
tau0 = 0 # in hours
tauf = 250 # in hours
dtau = 0.01 # in hours, cell growth update and uptake timescale
at = 1 # avg. consumption values (fmole per cell) alpha_ij: population i, resource j
bt = 0.1 # avg. production rates (fmole per cell per hour) beta_ij: population i, resource j
qp = 0.7 # probability of production link per population
qc = 0.7 # probability of influence link per population

basename = "nGen_$(nGen)_nCellType_$(nCellType)_nMediator_$(nMediator)_ri0_$(ri0)_posIntRatio_$(posIntRatio)_at_$(at)_bt_$(bt)_CS3390_HEATMAP_seed_$(rndseed0)"

@time NE0D,
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
    parse(Int64, ARGS[1]),
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
    parse(Int64, ARGS[2]),
)

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
)

println(maximum(NE0D))

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
