using Base

function DistInteractionStrengthMT_PB(Nc,Nm,ri0,fp, r)

    # Interaction matrix based on strength probability distribution B
    # Nc: number of interacting species
    # ri0: maximum interaction strength
    # fp: fraction for positive interactions
    return ri0*rand(r[myid()], Nc,Nm) .* sign.(rand(r[myid()], Nc,Nm) .- (1-fp))
end
