function NetworkConfig_Binomial(Nc, Nm, p, r)
    C = zeros(Nc,Nm)
    rndc = rand(r[myid()], Nc,Nm)
    C[map(x -> x <= p, rndc)] .= 1

    return C
end
