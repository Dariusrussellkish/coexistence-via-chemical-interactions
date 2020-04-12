
function WellmixedInteraction_DpMM_ExMT4(
    nRound,
    r0Vector,
    cellRatioArray,
    intMat,
    nInitialCell,
    kSatVector,
    A,
    B,
    kSatLevel,
    ExtTh,
    DilTh,
    tauf,
    dtau,
)
    tau0 = 0

    nCellType, nMediator = size(intMat)

    cMedVector = zeros(nMediator, 1)

    tauScaleArray = collect(tau0:dtau:tauf)
    nTauScale = length(tauScaleArray)

    nCellVector = nInitialCell * cellRatioArray'
    nCellVector0 = similar(nCellVector)

    count = 1

    onesNCellType = ones(1, nCellType)
    rIntPerCellVector = zeros(nCellType, 1)
    AMM = zeros(nMediator, nCellType)
    Ce = zeros(nMediator, nCellType)

    intMatLTZero = intMat .< 0
    intMatGTEZero = intMat .>= 0
    intMatPt1 = (intMatLTZero) .* intMat
    intMatPt2 = (intMatGTEZero) .* intMat

    rIntCMed1 = cMedVector ./ kSatVector
    rIntCMed21 = cMedVector + kSatVector
    rIntCMed2 = cMedVector ./ rIntCMed21

    rInt1 = (intMatPt1) * (rIntCMed1)
    rInt2 = (intMatPt2) * (rIntCMed2)
    rInt3 = rInt1 + rInt2
    rInt4 = r0Vector + rInt3

    cMedInner1 = B * nCellVector
    cMedInner2 = AMM * nCellVector
    cMedInner3 = cMedInner1 - cMedInner2
    cMedInner4 = dtau * cMedInner3
    cMedPart1 = cMedVector + cMedInner4
    cMedMask = cMedVector .< 0

    AMMPart1 = A .* Ce
    AMMPart2 = (Ce .+ kSatLevel)

    nCellVectorPart1 = rIntPerCellVector .* nCellVector
    nCellVectorPart2 = dtau * nCellVectorPart1
    nCellVectorPart3 = nCellVector + nCellVectorPart2
    nCellVectorMask = nCellVector .< ExtTh

    invNCellVecSum = 1.0 / sum(nCellVector)

    for iRound = 1:nRound
        # @info "Start Round: $iRound, nCellVector: $nCellVector"
        cMedVector = nInitialCell / sum(nCellVector) * cMedVector
        nCellVector = nInitialCell * cellRatioArray'

        tau0 = 0
        tau = tau0

        count = 1

        while (tau <= tauf - dtau) && (sum(nCellVector) < DilTh)
            tau = tauScaleArray[count]

            rIntCMed1 .= cMedVector ./ kSatVector
            rIntCMed21 .= cMedVector + kSatVector
            rIntCMed2 .= cMedVector ./ rIntCMed21

            rInt1 .= (intMatPt1) * (rIntCMed1)
            rInt2 .= (intMatPt2) * (rIntCMed2)
            rInt3 .= rInt1 .+ rInt2
            rInt4 .= r0Vector .+ rInt3
            rIntPerCellVector .= rInt4
            # println("rInt: $(size(rIntPerCellVector))")

            nCellVectorPart1 .= rIntPerCellVector .* nCellVector
            nCellVectorPart2 .= dtau .* nCellVectorPart1
            nCellVectorPart3 .= nCellVector .+ nCellVectorPart2
            nCellVector .= nCellVectorPart3

            nCellVectorMask .= nCellVector .< ExtTh
            nCellVector[nCellVectorMask] .= 0

            if count == 1
                nCellVector0 = deepcopy(nCellVector)
            end

            Ce .= cMedVector * onesNCellType # nM * nC matrix, each column is cMed

            AMMPart1 .= A .* Ce
            AMMPart2 .= (Ce .+ kSatLevel)

            AMM .= AMMPart1 ./ AMMPart2

            cMedInner1 .= B * nCellVector
            cMedInner2 .= AMM * nCellVector
            cMedInner3 .= cMedInner1 .- cMedInner2
            cMedInner4 .= dtau .* cMedInner3
            cMedPart1 .= cMedVector .+ cMedInner4

            cMedVector .= cMedPart1

            cMedMask .= (<).(cMedVector, 0)
            cMedVector[cMedMask] .= 0

            count += 1
        end

        invNCellVecSum = 1.0 / sum(nCellVector)
        cellRatioArray .= (*).(invNCellVecSum, nCellVector')
    end
    r = nCellVector ./ nCellVector0
    # @info "r: $r"
    maxr = maximum(x -> isnan(x) ? -Inf : x, r)[1]
    stp = r' .> (0.9 * maxr)
    indx = reshape(collect(1:nCellType), (1, nCellType))
    Ne = indx[stp] # vector of species ids that survive
    Cmp = cellRatioArray[Ne] # relative composition

    # get Cmp as percentage each cell type contributes to the total community
    if sum(Cmp) > 0
        Cmp_sum = zeros(1, size(Cmp, 2))
        Cmp_sum[1, :] .= sum(Cmp)

        Cmp = Cmp ./ Cmp_sum
    end
    return Ne, Cmp
end
